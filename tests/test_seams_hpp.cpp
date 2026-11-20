/* Copyright (c) d-SEAMS developers
 * The seams.hpp wrappers against the C++ engine census.
 */
#include "seams.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::map<std::string, std::vector<double>> load_oracle(const std::string &oracle, const std::string &root) {
  const std::string command = "\"" + oracle + "\" \"" + root + "\"";
#if defined(_WIN32)
  _putenv_s("OMP_NUM_THREADS", "1");
  FILE *pipe = _popen(command.c_str(), "r");
#else
  setenv("OMP_NUM_THREADS", "1", 1);
  FILE *pipe = popen(command.c_str(), "r");
#endif
  if (pipe == nullptr) {
    throw std::runtime_error("oracle failed to start");
  }
  std::string text;
  char buf[4096];
  while (fgets(buf, sizeof(buf), pipe) != nullptr) {
    text += buf;
  }
#if defined(_WIN32)
  const int status = _pclose(pipe);
#else
  const int status = pclose(pipe);
#endif
  if (status != 0) {
    throw std::runtime_error("oracle failed");
  }
  std::map<std::string, std::vector<double>> rows;
  std::istringstream lines(text);
  std::string line;
  while (std::getline(lines, line)) {
    std::istringstream parts(line);
    std::string name;
    std::string field;
    parts >> name >> field;
    std::vector<double> values;
    double value = 0.0;
    while (parts >> value) {
      values.push_back(value);
    }
    rows[name + " " + field] = values;
  }
  return rows;
}

double mean(const std::vector<double> &values) {
  double total = 0.0;
  for (double value : values) {
    total += value;
  }
  return values.empty() ? 0.0 : total / static_cast<double>(values.size());
}

void close(double got, double want) {
  const double scale = std::max(std::abs(got), std::max(std::abs(want), 1.0));
  if (std::abs(got - want) > 1e-8 * scale) {
    throw std::runtime_error("value mismatch");
  }
}

void equal_ints(const std::vector<int> &got, const std::vector<double> &want) {
  if (got.size() != want.size()) {
    throw std::runtime_error("length mismatch");
  }
  for (std::size_t i = 0; i < got.size(); ++i) {
    if (static_cast<double>(got[i]) != want[i]) {
      throw std::runtime_error("integer mismatch");
    }
  }
}

seams::System read_frame(const std::string &path, int type_filter) {
  seams::System sys;
  sys.read(path, 1, type_filter);
  return sys;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: test_seams_hpp oracle source-root\n";
    return 2;
  }
  try {
    if (seams_set_num_threads(1) != SEAMS_OK) {
      return 1;
    }
    const std::string root = argv[2];
    const auto want = load_oracle(argv[1], root);
    auto ih = read_frame(root + "/input/traj/ice1h.lammpstrj", -1);
    if (static_cast<double>(ih.count()) != want.at("ih n")[0]) {
      return 1;
    }
    const auto chill = ih.chill(seams::ChillOptions{3.5, 0, seams::Graph::Cutoff, seams::ChillRule::Published});
    equal_ints(std::vector<int>(chill.counts.begin(), chill.counts.end()), want.at("ih chill0"));
    const auto q3 = ih.steinhardt(seams::SteinhardtOptions{3.5, 0, 3, false});
    const auto q4 = ih.steinhardt(seams::SteinhardtOptions{3.5, 0, 4, false});
    const auto q6 = ih.steinhardt(seams::SteinhardtOptions{3.5, 0, 6, true});
    const auto q8 = ih.steinhardt(seams::SteinhardtOptions{3.5, 0, 8, false});
    const auto q12 = ih.steinhardt(seams::SteinhardtOptions{3.5, 0, 12, false});
    close(mean(q3.local), want.at("ih q3")[0]);
    close(mean(q4.local), want.at("ih q4")[0]);
    close(mean(q6.local), want.at("ih q6")[0]);
    close(mean(q6.average), want.at("ih q6bar")[0]);
    close(mean(q8.local), want.at("ih q8")[0]);
    close(mean(q12.local), want.at("ih q12")[0]);
    const auto first = ih.rings(seams::RingOptions{3.5, 0, seams::Graph::Cutoff, 6}).vertices(0);
    if (static_cast<double>(first.size()) != want.at("ih vert")[0] ||
        static_cast<double>(first.empty() ? -1 : first[0]) != want.at("ih vert")[1]) {
      return 1;
    }
    const auto cage_counts = ih.cages(seams::CageOptions{3.5, 0, seams::Graph::Seeded, 4, false}).counts();
    equal_ints(std::vector<int>{cage_counts.hexagonal, cage_counts.cubic, cage_counts.water, cage_counts.both},
               want.at("ih cages"));
    if (static_cast<double>(ih.knn(seams::KnnOptions{4, 3.5, 0, true}).nnz()) != want.at("ih knn")[0]) {
      return 1;
    }
    const auto rdf = ih.rdf(seams::DistributionOptions{6.0, 8, 0, 0});
    if (static_cast<double>(rdf.radius.size()) != want.at("ih rdf")[0]) {
      return 1;
    }
    close(rdf.volume, want.at("ih rdf")[1]);
    close(rdf.g.empty() ? 0.0 : rdf.g[0], want.at("ih rdf")[2]);
    const auto running = ih.running_cn(seams::DistributionOptions{6.0, 8, 0, 0});
    close(running.coordination.empty() ? 0.0 : running.coordination.back(), want.at("ih rcn")[0]);
    close(ih.coordination(seams::DistributionOptions{3.5, 8, 0, 0}), want.at("ih cn")[0]);
    if (static_cast<double>(ih.hbonds(seams::BondOptions{3.5, 1, 2, 2.42, 30.0})) != want.at("ih hbond")[0]) {
      return 1;
    }
    const auto density = ih.density(seams::DensityOptions{0, 4, seams::Axis::Z});
    if (static_cast<double>(density.density.size()) != want.at("ih rho")[0]) {
      return 1;
    }
    close(density.density.empty() ? 0.0 : density.density[0], want.at("ih rho")[1]);
    auto fp = ih.fingerprint(seams::FingerprintOptions{3.5, 0, 4, seams::Graph::Cutoff, 2, false});
    int ring_sum = 0;
    for (int count : fp.rings()) {
      ring_sum += count;
    }
    if (static_cast<double>(fp.key().size()) != want.at("ih fp")[0] ||
        static_cast<double>(ring_sum) != want.at("ih fp")[1]) {
      return 1;
    }
    if (static_cast<double>(ih.site_table("1=waterO")[0]) != want.at("ih site")[0]) {
      return 1;
    }
    const auto pairs = ih.contact_pairs("1=cationHead");
    if (static_cast<double>(pairs.pairs) != want.at("ih pairs")[0] ||
        static_cast<double>(pairs.cations) != want.at("ih pairs")[1] ||
        static_cast<double>(pairs.anions) != want.at("ih pairs")[2]) {
      return 1;
    }
    const auto domain = ih.domains(seams::DomainOptions{3.5, "1=waterO", seams::DomainSubset::Apolar});
    if (static_cast<double>(domain.subset) != want.at("ih domain")[0] ||
        static_cast<double>(domain.largest) != want.at("ih domain")[1]) {
      return 1;
    }
    close(domain.percolation, want.at("ih domain")[2]);
    if (static_cast<double>(ih.prisms(seams::PrismOptions{3.5, 0}).perfect[6]) != want.at("ih prism")[0]) {
      return 1;
    }
    if (static_cast<double>(ih.monolayer(seams::SheetOptions{3.5, 0, 6, 100.0})[6]) != want.at("ih sheet")[0]) {
      return 1;
    }
    auto guest = read_frame(root + "/input/traj/genice_sI.lammpstrj", 1);
    if (static_cast<double>(guest.count()) != want.at("si n")[0]) {
      return 1;
    }
    const auto sig = guest.signature_cages(seams::SignatureOptions{3.5, 0, seams::Graph::Cutoff, "512"});
    if (static_cast<double>(sig.cages) != want.at("si sig")[0] ||
        static_cast<double>(sig.atoms) != want.at("si sig")[1]) {
      return 1;
    }
    auto xyz = read_frame(root + "/input/xyz/configuration.xyz", -1);
    if (static_cast<double>(xyz.count()) != want.at("xyz n")[0]) {
      return 1;
    }
    close(xyz.coordination(seams::DistributionOptions{6.0, 8, 0, 0}), want.at("xyz cn")[0]);
    auto con = read_frame(root + "/input/eon/tiny_multi_cuh2.con", -1);
    if (static_cast<double>(con.count()) != want.at("con n")[0] ||
        static_cast<double>(con.hbonds(seams::BondOptions{3.5, 0, 2, 2.42, 30.0})) != want.at("con hbond")[0]) {
      return 1;
    }
    auto tri = read_frame(root + "/input/traj/genice_sH.lammpstrj", 1);
    if (static_cast<double>(tri.count()) != want.at("tri n")[0]) {
      return 1;
    }
    close(mean(tri.steinhardt(seams::SteinhardtOptions{3.5, 0, 6, false}).local), want.at("tri q6")[0]);
    // Chordless six-cycles through the wrapper match the C call, and every
    // cycle row is counted once on each of its atoms.
    auto water = read_frame(root + "/input/traj/exampleTraj.lammpstrj", 2);
    const auto cycles = water.rings_chordless(seams::NeighbourOptions{3.5, 0}, 6);
    {
      // The same frame through the C calls the wrapper sits on.
      seams_system *raw = seams_system_create();
      seams_neighbours *nb = nullptr;
      std::size_t n_c = 0;
      std::vector<int> rows;
      bool ok = raw != nullptr &&
                seams_system_read(raw, (root + "/input/traj/exampleTraj.lammpstrj").c_str(), 1,
                                  2) == SEAMS_OK &&
                seams_neighbours_build(raw, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) == SEAMS_OK &&
                seams_rings_chordless(nb, 6, nullptr, 0, &n_c, nullptr, 0) == SEAMS_OK;
      if (ok) {
        rows.resize(n_c * 6);
        ok = seams_rings_chordless(nb, 6, rows.data(), rows.size(), &n_c, nullptr, 0) == SEAMS_OK;
      }
      seams_neighbours_destroy(nb);
      seams_system_destroy(raw);
      if (!ok || n_c == 0 || n_c != cycles.count() || rows != cycles.atoms) {
        return 1;
      }
    }
    std::vector<int> tally(cycles.per_atom.size(), 0);
    for (std::size_t c = 0; c < cycles.count(); ++c) {
      const auto row = cycles.cycle(c);
      if (row.size() != 6) {
        return 1;
      }
      for (const int atom : row) {
        tally[static_cast<std::size_t>(atom)] += 1;
      }
    }
    if (tally != cycles.per_atom) {
      return 1;
    }

    // Trajectory unwrapping, dipoles, order parameters, and TCF RAII wrappers.
    {
      const std::vector<double> w_xyz = {0.0, 0.0, 0.0, 0.75, 0.0, 0.58, -0.75, 0.0, 0.58};
      const std::vector<int> w_types = {1, 2, 2};
      const std::vector<int> w_mols = {7, 7, 7};
      const std::array<double, 3> w_origin = {0.0, 0.0, 0.0};
      const std::array<double, 3> w_lengths = {10.0, 10.0, 10.0};
      const std::array<double, 3> z_axis = {0.0, 0.0, 1.0};
      seams::System wsys;
      wsys.set_ortho(w_xyz, w_origin, w_lengths, w_types);
      wsys.set_mols(w_mols);
      if (wsys.mols() != w_mols) {
        return 1;
      }
      const auto mc = wsys.molecule_centres();
      const auto wv = wsys.water_vectors(1, 2);
      const std::vector<double> q = {-0.8476, 0.4238, 0.4238};
      const auto dip = wsys.molecule_dipoles({}, q);
      const auto ord = seams::orientation_order(wv.bisector, z_axis);
      if (mc.mol_ids.size() != 1 || mc.mol_ids[0] != 7 || wv.mol_ids.size() != 1 ||
          dip.dipoles.size() != 3 || std::abs(ord.mean_p1 - 1.0) > 1e-12 ||
          std::abs(ord.mean_p2 - 1.0) > 1e-12) {
        return 1;
      }

      const std::array<double, 9> box9 = {0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
      seams::Msd msd(1, 4, 1, 1.0, z_axis);
      for (int f = 0; f < 5; ++f) {
        std::vector<double> p = {std::fmod(1.0 + 2.0 * f, 10.0), 1.0,
                                 std::fmod(1.0 + 3.0 * f, 10.0)};
        msd.push(p, box9);
      }
      const auto mc_fit = msd.fit(1.0, 4.0);
      if (std::abs(mc_fit.alpha - 2.0) > 1e-10) {
        return 1;
      }

      seams::HbondTcf hb(4, 1, 1.0, false);
      const std::array<int, 1> d0 = {0};
      const std::array<int, 1> a1 = {1};
      hb.push(d0, a1);
      hb.push({}, {});
      hb.push(d0, a1);
      const auto hbc = hb.curve();
      if (hbc.lag_time.size() != 3 || std::abs(hbc.c_intermittent[2] - 1.0) > 1e-12 ||
          std::abs(hbc.s_continuous[2] - 0.0) > 1e-12) {
        return 1;
      }

      seams::Reorientation reo(1, 3, 1, 1.0);
      reo.push(std::vector<double>{0.0, 0.0, 1.0});
      reo.push(std::vector<double>{1.0, 0.0, 0.0});
      const auto rc = reo.curve();
      if (rc.lag_time.size() != 2 || std::abs(rc.c1[1]) > 1e-12 ||
          std::abs(rc.c2[1] + 0.5) > 1e-12) {
        return 1;
      }
    }
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << "\n";
    return 1;
  }
  std::cout << "seams.hpp ok\n";
  return 0;
}
