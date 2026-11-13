/* Copyright (c) d-SEAMS developers
 * The seams.hpp wrappers against the C++ engine census.
 */
#include "seams.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::map<std::string, std::vector<double>> load_oracle(const std::string &oracle, const std::string &root) {
  const std::string command = "OMP_NUM_THREADS=1 \"" + oracle + "\" \"" + root + "\"";
  FILE *pipe = popen(command.c_str(), "r");
  if (pipe == nullptr) {
    throw std::runtime_error("oracle failed to start");
  }
  std::string text;
  char buf[4096];
  while (fgets(buf, sizeof(buf), pipe) != nullptr) {
    text += buf;
  }
  const int status = pclose(pipe);
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
    const auto chill = ih.chill(3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED);
    equal_ints(chill.second, want.at("ih chill0"));
    const auto q3 = ih.steinhardt(3.5, 3, false);
    const auto q4 = ih.steinhardt(3.5, 4, false);
    const auto q6 = ih.steinhardt(3.5, 6, true);
    const auto q8 = ih.steinhardt(3.5, 8, false);
    const auto q12 = ih.steinhardt(3.5, 12, false);
    close(mean(q3.first), want.at("ih q3")[0]);
    close(mean(q4.first), want.at("ih q4")[0]);
    close(mean(q6.first), want.at("ih q6")[0]);
    close(mean(q6.second), want.at("ih q6bar")[0]);
    close(mean(q8.first), want.at("ih q8")[0]);
    close(mean(q12.first), want.at("ih q12")[0]);
    const auto first = ih.rings(3.5, 6).vertices(0);
    if (static_cast<double>(first.size()) != want.at("ih vert")[0] ||
        static_cast<double>(first.empty() ? -1 : first[0]) != want.at("ih vert")[1]) {
      return 1;
    }
    int n_hex = 0;
    int n_cubic = 0;
    int n_water = 0;
    int n_both = 0;
    ih.cages(3.5, SEAMS_GRAPH_SEEDED, 4, false).counts(n_hex, n_cubic, n_water, n_both);
    equal_ints(std::vector<int>{n_hex, n_cubic, n_water, n_both}, want.at("ih cages"));
    if (static_cast<double>(ih.knn(4, 3.5, true).nnz()) != want.at("ih knn")[0]) {
      return 1;
    }
    const auto rdf = ih.rdf(6.0, 8);
    if (static_cast<double>(std::get<0>(rdf).size()) != want.at("ih rdf")[0]) {
      return 1;
    }
    close(std::get<2>(rdf), want.at("ih rdf")[1]);
    close(std::get<1>(rdf).empty() ? 0.0 : std::get<1>(rdf)[0], want.at("ih rdf")[2]);
    const auto running = ih.running_cn(6.0, 8);
    close(running.empty() ? 0.0 : running.back(), want.at("ih rcn")[0]);
    close(ih.coordination(3.5, 8), want.at("ih cn")[0]);
    if (static_cast<double>(ih.hbonds(3.5, 1, 2, 2.42, 30.0)) != want.at("ih hbond")[0]) {
      return 1;
    }
    const auto density = ih.density_axis(0, 4, 2);
    if (static_cast<double>(density.second.size()) != want.at("ih rho")[0]) {
      return 1;
    }
    close(density.second.empty() ? 0.0 : density.second[0], want.at("ih rho")[1]);
    auto fp = ih.fingerprint(3.5, 2);
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
    if (static_cast<double>(std::get<0>(pairs)) != want.at("ih pairs")[0] ||
        static_cast<double>(std::get<1>(pairs)) != want.at("ih pairs")[1] ||
        static_cast<double>(std::get<2>(pairs)) != want.at("ih pairs")[2]) {
      return 1;
    }
    const auto domain = ih.domains(3.5, "1=waterO", 0);
    if (static_cast<double>(std::get<0>(domain)) != want.at("ih domain")[0] ||
        static_cast<double>(std::get<1>(domain)) != want.at("ih domain")[1]) {
      return 1;
    }
    close(std::get<2>(domain), want.at("ih domain")[2]);
    if (static_cast<double>(ih.prisms(3.5)[6]) != want.at("ih prism")[0]) {
      return 1;
    }
    if (static_cast<double>(ih.monolayer(3.5, 6, 100.0)[6]) != want.at("ih sheet")[0]) {
      return 1;
    }
    auto guest = read_frame(root + "/input/traj/genice_sI.lammpstrj", 1);
    if (static_cast<double>(guest.count()) != want.at("si n")[0]) {
      return 1;
    }
    const auto sig = guest.signature_cages(3.5, SEAMS_GRAPH_CUTOFF, "512");
    if (static_cast<double>(sig.first) != want.at("si sig")[0] ||
        static_cast<double>(sig.second) != want.at("si sig")[1]) {
      return 1;
    }
    auto xyz = read_frame(root + "/input/xyz/configuration.xyz", -1);
    if (static_cast<double>(xyz.count()) != want.at("xyz n")[0]) {
      return 1;
    }
    close(xyz.coordination(6.0, 8), want.at("xyz cn")[0]);
    auto con = read_frame(root + "/input/con/tiny_multi_cuh2.con", -1);
    if (static_cast<double>(con.count()) != want.at("con n")[0] ||
        static_cast<double>(con.hbonds(3.5, 0, 2, 2.42, 30.0)) != want.at("con hbond")[0]) {
      return 1;
    }
    auto tri = read_frame(root + "/input/traj/genice_sH.lammpstrj", 1);
    if (static_cast<double>(tri.count()) != want.at("tri n")[0]) {
      return 1;
    }
    close(mean(tri.steinhardt(3.5, 6, false).first), want.at("tri q6")[0]);
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << "\n";
    return 1;
  }
  std::cout << "seams.hpp ok\n";
  return 0;
}
