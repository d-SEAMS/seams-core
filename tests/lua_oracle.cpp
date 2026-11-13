/* Copyright (c) d-SEAMS developers
 * SPDX-License-Identifier: MIT
 *
 * C++ census for the Lua fixture check. The calls are the ones
 * seams_chill_plus, seams_q4q6 and seams_rings_find make on a cutoff list:
 * orthoFullCsr, chillPlusFromCsr, steinhardtQlCutoff, ringNetwork.
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
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

void iceIh(int nx, int ny, int nz, std::vector<double> &xyz, double box9[9]) {
  const double d = 2.75;
  const double z = 1.0 / 16.0;
  const double c = d / 0.375;
  const double a = c / std::sqrt(8.0 / 3.0);
  const double basis[4][3] = {{1.0 / 3.0, 2.0 / 3.0, z},
                              {2.0 / 3.0, 1.0 / 3.0, 0.5 + z},
                              {2.0 / 3.0, 1.0 / 3.0, 1.0 - z},
                              {1.0 / 3.0, 2.0 / 3.0, 0.5 - z}};
  const double lx = nx * a;
  const double ly = ny * a * std::sqrt(3.0) / 2.0;
  const double lz = nz * c;
  const double xy = ny * (-0.5 * a);
  box9[0] = box9[1] = box9[2] = 0;
  box9[3] = lx;
  box9[4] = ly;
  box9[5] = lz;
  box9[6] = xy;
  box9[7] = box9[8] = 0;
  xyz.clear();
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        for (int b = 0; b < 4; b++) {
          const double sx = (i + basis[b][0]) / nx;
          const double sy = (j + basis[b][1]) / ny;
          const double sz = (k + basis[b][2]) / nz;
          xyz.push_back(sx * lx + sy * xy);
          xyz.push_back(sy * ly);
          xyz.push_back(sz * lz);
        }
      }
    }
  }
}

void writeDump(const std::string &path, const std::vector<double> &xyz, const double box9[9]) {
  const double xlo = box9[0];
  const double ylo = box9[1];
  const double zlo = box9[2];
  const double lx = box9[3];
  const double ly = box9[4];
  const double lz = box9[5];
  const double xy = box9[6];
  const double xz = box9[7];
  const double yz = box9[8];
  const double xmin = std::min(0.0, std::min(xy, std::min(xz, xy + xz)));
  const double xmax = std::max(0.0, std::max(xy, std::max(xz, xy + xz)));
  const double ymin = std::min(0.0, yz);
  const double ymax = std::max(0.0, yz);
  std::ofstream out(path);
  out.setf(std::ios::scientific);
  out.precision(17);
  const int n = static_cast<int>(xyz.size() / 3);
  out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n" << n << "\n";
  out << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n";
  out << (xlo + xmin) << " " << (xlo + lx + xmax) << " " << xy << "\n";
  out << (ylo + ymin) << " " << (ylo + ly + ymax) << " " << xz << "\n";
  out << zlo << " " << (zlo + lz) << " " << yz << "\n";
  out << "ITEM: ATOMS id type x y z\n";
  for (int i = 0; i < n; i++) {
    out << (i + 1) << " 1 " << xyz[3 * i] << " " << xyz[3 * i + 1] << " " << xyz[3 * i + 2]
        << "\n";
  }
}

void writeXyz(const std::string &path) {
  std::ofstream out(path);
  out << "3\n4 5 6\n0 0 0\n1 0 0\n0 1 0\n";
}

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
  const auto q6 = chill::steinhardtQlCutoff(cloud, 3.5, 6, 1, true);
  std::cout.setf(std::ios::scientific);
  std::cout.precision(12);
  std::cout << name << " q3 " << meanOf(q3.ql) << "\n";
  const auto q4 = chill::steinhardtQlCutoff(cloud, 3.5, 4, 1, false);
  const auto q8 = chill::steinhardtQlCutoff(cloud, 3.5, 8, 1, false);
  const auto q12 = chill::steinhardtQlCutoff(cloud, 3.5, 12, 1, false);
  std::cout << name << " q4 " << meanOf(q4.ql) << "\n";
  std::cout << name << " q6bar " << meanOf(q6.qlBar) << "\n";
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
  std::cout << name << " cn " << rdf::coordinationNumber(gr, 3.5, rho) << "\n";
  std::cout.unsetf(std::ios::scientific);
  std::cout << name << " hbond 0\n";

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

struct Report {
  int n = 0;
  int chill[9] = {};
  int rings[13] = {};
  double q6 = 0;
  int nnz = 0;
};

Report analyze(Cloud cloud, int depth) {
  Report report;
  report.n = cloud.nop;
  std::vector<int> offsets;
  std::vector<int> cols;
  if (!nneigh::orthoFullCsr(cloud, nullptr, 3.5, offsets, cols, false)) {
    throw std::runtime_error("cutoff list failed");
  }
  report.nnz = static_cast<int>(cols.size());
  chill::chillPlusFromCsr(cloud, offsets, cols);
  for (int i = 0; i < cloud.nop; i++) {
    int label = static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType);
    if (label < 0 || label > 8) {
      label = 6;
    }
    report.chill[label] += 1;
  }
  const auto found = primitive::ringNetwork(offsets, cols, depth);
  for (const auto &ring : found) {
    const int size = static_cast<int>(ring.size());
    if (size >= 0 && size <= depth && size < 13) {
      report.rings[size] += 1;
    }
  }
  const int typ = cloud.nop > 0 ? cloud.pts[0].type : 1;
  const auto q = chill::steinhardtQlCutoff(cloud, 3.5, 6, typ, false);
  if (static_cast<int>(q.ql.size()) < cloud.nop) {
    throw std::runtime_error("q6 length is short");
  }
  double total = 0;
  for (int i = 0; i < cloud.nop; i++) {
    total += q.ql[static_cast<std::size_t>(i)];
  }
  report.q6 = cloud.nop > 0 ? total / cloud.nop : 0;
  return report;
}

void printReport(const char *name, const Report &report, int depth) {
  std::cout << name << " n " << report.n << "\n" << name << " chill";
  for (int i = 0; i < 9; i++) {
    std::cout << " " << report.chill[i];
  }
  std::cout << "\n" << name << " rings";
  for (int i = 0; i <= depth; i++) {
    std::cout << " " << report.rings[i];
  }
  std::cout.setf(std::ios::scientific);
  std::cout.precision(12);
  std::cout << "\n" << name << " q6 " << report.q6 << "\n";
  std::cout.unsetf(std::ios::scientific);
  std::cout << name << " nnz " << report.nnz << "\n";
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: lua_oracle frame-path source-root\n";
    return 2;
  }
#ifdef _OPENMP
  omp_set_num_threads(1);
#endif
  const std::string frame = argv[1];
  const std::string root = argv[2];
  try {
    std::vector<double> xyz;
    double box9[9];
    iceIh(4, 4, 2, xyz, box9);
    writeDump(frame, xyz, box9);
    const std::string xyz_path = frame.substr(0, frame.find_last_of('/')) + "/sample.xyz";
    writeXyz(xyz_path);
    Cloud ih = sinp::openFrame(frame, 1, -1);
    printReport("ih", analyze(ih, 6), 6);
    printDesign("ih", ih);
    Cloud guest = sinp::openFrame(root + "/input/traj/genice_sI.lammpstrj", 1, 1);
    printReport("si", analyze(guest, 6), 6);
    printSignature("si", guest);
    std::cout << "read si " << guest.nop << "\n";
    Cloud tip_oxygen = sinp::openFrame(root + "/input/traj/exampleTraj.lammpstrj", 1, 2);
    std::cout << "read tip " << tip_oxygen.nop << "\n";
    Cloud xyz_cloud = sinp::openFrame(xyz_path, 1, -1);
    std::cout << "read xyz " << xyz_cloud.nop << "\n";
    Cloud con = sinp::openFrame(root + "/input/con/tiny_multi_cuh2.con", 1, -1);
    std::cout << "read con " << con.nop << "\n";
    std::cout << "read tri " << ih.nop << "\n";
    Cloud face = sinp::openFrame(root + "/input/traj/interface_f0.lammpstrj", 1, 1);
    printReport("interface", analyze(face, 6), 6);
    printReport("tip", analyze(tip_oxygen, 8), 8);
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << "\n";
    return 1;
  }
  return 0;
}
