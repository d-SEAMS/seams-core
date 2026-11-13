/* Copyright (c) d-SEAMS developers
 * SPDX-License-Identifier: MIT
 *
 * C++ census for the Rust binding. The calls are the engine functions
 * behind seams_chill, seams_steinhardt, seams_cages, seams_signature_cages,
 * seams_neighbours_knn, seams_rdf, seams_running_cn, seams_coordination,
 * seams_hbonds, seams_density_axis, seams_fingerprint, seams_site_table,
 * seams_contact_pairs, seams_domains, seams_prisms, seams_monolayer, and
 * seams_rings_vertices.
 */

#include <bop.hpp>
#include <bond.hpp>
#include <cage_affiliation.hpp>
#include <cage_enum.hpp>
#include <chill_plus_scheme.hpp>
#include <cluster.hpp>
#include <density.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <rdf.hpp>
#include <seams_input.hpp>
#include <site.hpp>
#include <topo_bulk.hpp>
#include <topo_fingerprint.hpp>
#include <topo_one_dim.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

double meanOf(const std::vector<double> &values) {
  double total = 0.0;
  for (double value : values) {
    total += value;
  }
  return values.empty() ? 0.0 : total / static_cast<double>(values.size());
}

std::vector<std::vector<int>> sixRings(const std::vector<std::vector<int>> &rings) {
  std::vector<std::vector<int>> six;
  for (const auto &ring : rings) {
    if (ring.size() == 6) {
      six.push_back(ring);
    }
  }
  return six;
}

int hbondCount(const Cloud &cloud) {
  Cloud heavy;
  Cloud hydro;
  heavy.box = cloud.box;
  heavy.boxLow = cloud.boxLow;
  hydro.box = cloud.box;
  hydro.boxLow = cloud.boxLow;
  const int hydrogen_type = 2;
  for (int i = 0; i < cloud.nop; i++) {
    const auto &point = cloud.pts[static_cast<std::size_t>(i)];
    if (point.type == hydrogen_type) {
      hydro.pts.push_back(point);
    } else {
      heavy.pts.push_back(point);
    }
  }
  heavy.nop = static_cast<int>(heavy.pts.size());
  hydro.nop = static_cast<int>(hydro.pts.size());
  for (int i = 0; i < heavy.nop; i++) {
    heavy.idIndexMap[heavy.pts[static_cast<std::size_t>(i)].atomID] = i;
  }
  for (int i = 0; i < hydro.nop; i++) {
    hydro.idIndexMap[hydro.pts[static_cast<std::size_t>(i)].atomID] = i;
  }
  if (heavy.nop == 0 || hydro.nop == 0) {
    return 0;
  }
  const int typ = heavy.pts[0].type;
  auto rows = nneigh::neighListO(3.5, heavy, typ);
  auto net = bond::populateHbondsWithInputClouds(heavy, hydro, rows, 2.42, 30.0);
  int edges = 0;
  for (const auto &row : net) {
    if (row.size() > 1) {
      edges += static_cast<int>(row.size()) - 1;
    }
  }
  return edges / 2;
}

double coordinationOf(const Cloud &cloud, double rmax, int bins) {
  const auto gr = rdf::partialRdf(cloud, 1, 1, rmax, bins);
  const double rho = gr.volume > 0.0 ? static_cast<double>(gr.nJ) / gr.volume : 0.0;
  return rdf::coordinationNumber(gr, rmax, rho);
}

void printDesign(const char *name, Cloud cloud) {
  auto rows = nneigh::neighListO(3.5, cloud, 1);
  chill::getCorrel(cloud, rows, false);
  chill::getIceTypeNoPrint(cloud, rows, false);
  int chill[9] = {};
  for (int i = 0; i < cloud.nop; i++) {
    int label = static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType);
    if (label < 0 || label > 8) {
      label = 6;
    }
    chill[label] += 1;
  }
  std::cout << name << " chill0";
  for (int value : chill) {
    std::cout << " " << value;
  }
  std::cout << "\n";

  const auto q3 = chill::steinhardtQlCutoff(cloud, 3.5, 3, 1, false);
  const auto q4 = chill::steinhardtQlCutoff(cloud, 3.5, 4, 1, false);
  const auto q6 = chill::steinhardtQlCutoff(cloud, 3.5, 6, 1, false);
  const auto q6bar = chill::steinhardtQlCutoff(cloud, 3.5, 6, 1, true);
  const auto q8 = chill::steinhardtQlCutoff(cloud, 3.5, 8, 1, false);
  const auto q12 = chill::steinhardtQlCutoff(cloud, 3.5, 12, 1, false);
  std::cout.setf(std::ios::scientific);
  std::cout.precision(12);
  std::cout << name << " q3 " << meanOf(q3.ql) << "\n";
  std::cout << name << " q4 " << meanOf(q4.ql) << "\n";
  std::cout << name << " q6 " << meanOf(q6.ql) << "\n";
  std::cout << name << " q6bar " << meanOf(q6bar.qlBar) << "\n";
  std::cout << name << " q8 " << meanOf(q8.ql) << "\n";
  std::cout << name << " q12 " << meanOf(q12.ql) << "\n";
  std::cout.unsetf(std::ios::scientific);

  std::vector<int> csr_off;
  std::vector<int> csr_col;
  if (!nneigh::orthoFullCsr(cloud, nullptr, 3.5, csr_off, csr_col, false)) {
    throw std::runtime_error("cutoff list failed");
  }
  const auto rings = primitive::ringNetwork(csr_off, csr_col, 6);
  auto idx = nneigh::neighbourListByIndex(cloud, rows);
  const int nvert = rings.empty() ? 0 : static_cast<int>(rings[0].size());
  const int first = rings.empty() || rings[0].empty() ? -1 : rings[0][0];
  std::cout << name << " vert " << nvert << " " << first << "\n";

  auto graphs = nneigh::kNearestNeighbourPair(cloud, 4, 5.0, 1);
  auto idx_s = nneigh::neighbourListByIndex(cloud, graphs.first);
  auto idx_u = nneigh::neighbourListByIndex(cloud, graphs.second);
  auto six_s = sixRings(ring::contractibleRings(primitive::ringNetwork(idx_s, 6), cloud));
  auto six_u = sixRings(ring::contractibleRings(primitive::ringNetwork(idx_u, 6), cloud));
  ring::SeededAtomLabels aff;
  if (graphs.first == graphs.second) {
    const auto once = ring::cageAffiliation(six_s, idx_s);
    aff = ring::seededCageAffiliation(six_s, once, six_s, once, idx_s, false);
  } else {
    aff = ring::seededCageAffiliation(six_s, idx_s, six_u, idx_u, false);
  }
  int n_hex = 0;
  int n_cubic = 0;
  int n_water = 0;
  int n_both = 0;
  const int m = std::min(cloud.nop, static_cast<int>(aff.hc.size()));
  for (int i = 0; i < cloud.nop; i++) {
    const bool hexagonal = i < m && aff.hc[static_cast<std::size_t>(i)] != 0;
    const bool cubic = i < m && aff.ddc[static_cast<std::size_t>(i)] != 0;
    if (hexagonal && cubic) {
      n_both += 1;
    } else if (hexagonal) {
      n_hex += 1;
    } else if (cubic) {
      n_cubic += 1;
    } else {
      n_water += 1;
    }
  }
  std::cout << name << " cages " << n_hex << " " << n_cubic << " " << n_water << " " << n_both
            << "\n";

  auto knn = nneigh::kNearestNeighbourList(cloud, 4, 5.0, 1, true);
  auto knn_idx = nneigh::neighbourListByIndex(cloud, knn);
  int knn_nnz = 0;
  for (int i = 0; i < cloud.nop && i < static_cast<int>(knn_idx.size()); i++) {
    for (int id : knn_idx[static_cast<std::size_t>(i)]) {
      if (id != i) {
        knn_nnz += 1;
      }
    }
  }
  std::cout << name << " knn " << knn_nnz << "\n";

  const auto gr = rdf::partialRdf(cloud, 1, 1, 6.0, 8);
  const double g0 = gr.g.empty() ? 0.0 : gr.g[0];
  std::cout.setf(std::ios::scientific);
  std::cout.precision(12);
  std::cout << name << " rdf " << gr.r.size() << " " << gr.volume << " " << g0 << "\n";
  const double rho = gr.volume > 0.0 ? static_cast<double>(gr.nJ) / gr.volume : 0.0;
  const auto series = rdf::runningCN(gr, rho);
  std::cout << name << " rcn " << (series.empty() ? 0.0 : series.back()) << "\n";
  std::cout << name << " cn " << coordinationOf(cloud, 3.5, 8) << "\n";
  std::cout.unsetf(std::ios::scientific);
  std::cout << name << " hbond " << hbondCount(cloud) << "\n";

  const auto density = site::densityZ(cloud, 0, 4, 2);
  const double rho0 = density.rho.empty() ? 0.0 : density.rho[0];
  std::cout.setf(std::ios::scientific);
  std::cout.precision(12);
  std::cout << name << " rho " << density.rho.size() << " " << rho0 << "\n";
  std::cout.unsetf(std::ios::scientific);

  const auto fp = topo::fingerprint(idx, 2, 7, {});
  int ring_sum = 0;
  for (int count : fp.ringCensus) {
    ring_sum += count;
  }
  std::cout << name << " fp " << fp.key.size() << " " << ring_sum << "\n";
  const auto table = site::parseSiteSpec("1=waterO");
  const int kind0 = cloud.nop > 0 ? static_cast<int>(table.of(cloud.pts[0])) : -1;
  std::cout << name << " site " << kind0 << "\n";
  const auto ions = site::ionCloud(cloud, site::parseSiteSpec("1=cationHead"));
  const auto contacts = nneigh::mutualNearestUnlike(ions, 1, 2);
  int cations = 0;
  int anions = 0;
  for (const auto &point : ions.pts) {
    cations += point.type == 1 ? 1 : 0;
    anions += point.type == 2 ? 1 : 0;
  }
  std::cout << name << " pairs " << contacts.size() << " " << cations << " " << anions << "\n";
  std::vector<bool> mask(static_cast<std::size_t>(cloud.nop), false);
  for (const int index : site::indicesOf(cloud, site::parseSiteSpec("1=waterO"), site::Kind::apolar)) {
    if (index >= 0 && index < cloud.nop) {
      mask[static_cast<std::size_t>(index)] = true;
    }
  }
  auto by_index = nneigh::getNewNeighbourListByIndex(cloud, 3.5);
  std::vector<std::vector<int>> by_id(by_index.size());
  for (std::size_t i = 0; i < by_index.size(); ++i) {
    for (const int index : by_index[i]) {
      if (index >= 0 && index < cloud.nop) {
        by_id[i].push_back(cloud.pts[static_cast<std::size_t>(index)].atomID);
      }
    }
  }
  const auto domain = clump::largestDomain(cloud, by_id, mask);
  std::cout.setf(std::ios::scientific);
  std::cout.precision(12);
  std::cout << name << " domain " << domain.subset << " " << domain.largest << " "
            << domain.percolation << "\n";
  std::cout.unsetf(std::ios::scientific);

  auto knn_rows = nneigh::kNearestNeighbourList(cloud, 4, 5.0, 1, true);
  auto prism_idx = nneigh::neighbourListByIndex(cloud, knn_rows);
  const auto prism_rings = primitive::ringNetwork(prism_idx, 6);
  int perfect6 = 0;
  {
    std::vector<std::vector<int>> one;
    for (const auto &ring : prism_rings) {
      if (ring.size() == 6) {
        one.push_back(ring);
      }
    }
    std::vector<ring::strucType> kind(one.size());
    std::vector<double> rmsd(static_cast<std::size_t>(std::max(cloud.nop, 0)), -1.0);
    int n_perfect = 0;
    int n_imperfect = 0;
    if (!one.empty()) {
      ring::findPrisms(one, kind, n_perfect, n_imperfect, prism_idx, cloud, rmsd, false);
    }
    perfect6 = n_perfect;
  }
  std::cout << name << " prism " << perfect6 << "\n";
  int sheet6 = 0;
  for (const auto &ring : primitive::ringNetwork(idx, 6)) {
    if (ring.size() == 6) {
      sheet6 += 1;
    }
  }
  std::cout << name << " sheet " << sheet6 << "\n";
}

void printSignature(const char *name, Cloud cloud) {
  auto rows = nneigh::neighListO(3.5, cloud, 1);
  auto idx = nneigh::neighbourListByIndex(cloud, rows);
  const auto signature = cage::Signature::parse("512");
  const int depth = std::max(signature.maxRingSize(), 3);
  const auto rings = ring::cageRings(primitive::ringNetwork(idx, depth), cloud, 3.5);
  const auto window = cage::parseEdgeWindow("");
  std::vector<std::array<double, 3>> xyz(static_cast<std::size_t>(std::max(cloud.nop, 0)));
  for (int i = 0; i < cloud.nop; i++) {
    const auto &point = cloud.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(i)] = {point.x, point.y, point.z};
  }
  const auto found_of = cage::findBySignatures(rings, idx, std::vector<cage::Signature>{signature},
                                               {}, &xyz, &cloud.box, window);
  const auto &found = found_of.at(0);
  std::set<int> atoms;
  for (const auto &cage_i : found) {
    atoms.insert(cage_i.vertices.begin(), cage_i.vertices.end());
  }
  std::cout << name << " sig " << found.size() << " " << atoms.size() << "\n";
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: binding_oracle source-root\n";
    return 2;
  }
#ifdef _OPENMP
  omp_set_num_threads(1);
#endif
  const std::string root = argv[1];
  try {
    Cloud ih = sinp::openFrame(root + "/input/traj/ice1h.lammpstrj", 1, -1);
    std::cout << "ih n " << ih.nop << "\n";
    printDesign("ih", ih);
    Cloud guest = sinp::openFrame(root + "/input/traj/genice_sI.lammpstrj", 1, 1);
    std::cout << "si n " << guest.nop << "\n";
    printSignature("si", guest);
    Cloud xyz = sinp::openFrame(root + "/input/xyz/configuration.xyz", 1, -1);
    std::cout.setf(std::ios::scientific);
    std::cout.precision(12);
    std::cout << "xyz n " << xyz.nop << "\n";
    std::cout << "xyz cn " << coordinationOf(xyz, 6.0, 8) << "\n";
    Cloud con = sinp::openFrame(root + "/input/con/tiny_multi_cuh2.con", 1, -1);
    std::cout << "con n " << con.nop << "\n";
    std::cout.unsetf(std::ios::scientific);
    std::cout << "con hbond " << hbondCount(con) << "\n";
    Cloud tri = sinp::openFrame(root + "/input/traj/genice_sH.lammpstrj", 1, 1);
    const auto q6 = chill::steinhardtQlCutoff(tri, 3.5, 6, 1, false);
    std::cout.setf(std::ios::scientific);
    std::cout.precision(12);
    std::cout << "tri n " << tri.nop << "\n";
    std::cout << "tri q6 " << meanOf(q6.ql) << "\n";
    std::cout.unsetf(std::ios::scientific);
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << "\n";
    return 1;
  }
  return 0;
}
