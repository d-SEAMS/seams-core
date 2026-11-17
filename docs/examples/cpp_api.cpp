/* Copyright (c) d-SEAMS developers
 * Idiomatic seams.hpp on the eight-atom diamond and on ice Ih.
 */
#include "seams.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

std::vector<double> diamond8(double &edge) {
  const double a = 2.75 * 4.0 / std::sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  std::vector<double> xyz;
  xyz.reserve(24);
  for (const auto &shift : fcc) {
    for (const auto &site : basis) {
      xyz.push_back((shift[0] + site[0]) * a);
      xyz.push_back((shift[1] + site[1]) * a);
      xyz.push_back((shift[2] + site[2]) * a);
    }
  }
  edge = a;
  return xyz;
}

}  // namespace

int main(int argc, char **argv) {
  try {
    seams::set_num_threads(1);
    if (seams::num_threads() != 1 || seams::version().empty()) {
      return 1;
    }
    double edge = 0;
    const auto xyz = diamond8(edge);
    seams::System cell;
    cell.set_ortho(xyz, edge, edge, edge);
    const auto cubic = cell.chill_plus(seams::ChillOptions{});
    const auto orders = cell.q4q6(3.5);
    double mean = 0;
    for (double value : orders.q6) {
      mean += value;
    }
    mean /= static_cast<double>(orders.q6.size());
    if (cubic.counts[0] != 8 || !(mean > 0.3)) {
      return 1;
    }
    std::printf("cpp api cubic %d q6 %.3f %s\n", cubic.counts[0], mean, seams::version().c_str());
    if (argc < 2) {
      return 0;
    }
    seams::System ice;
    ice.read(std::string(argv[1]) + "/input/traj/ice1h.lammpstrj", 1, -1);
    const auto labels = ice.chill_plus(seams::ChillOptions{});
    const auto q6 = ice.q6bar(seams::Q6BarOptions{});
    const auto shell = ice.symmetry(0.0);
    int n_ih = 0;
    for (int label : shell) {
      n_ih += label == 0 ? 1 : 0;
    }
    std::vector<int> flags(static_cast<std::size_t>(ice.count()), 1);
    const auto front = ice.brine(flags, seams::BrineOptions{});
    const auto pockets = front.pockets();
    std::printf("ice Ih atoms %d hexagonal %d q6bar %.12f shell_ih %d pockets %d\n", ice.count(),
                labels.counts[static_cast<std::size_t>(static_cast<int>(seams::Label::Hexagonal))],
                q6.mean, n_ih, pockets.pockets);
    if (ice.count() != 128 || labels.counts[1] != 128) {
      return 1;
    }
    if (std::abs(q6.mean - 0.562448843702) > 1e-9) {
      return 1;
    }
  } catch (const seams::Error &err) {
    std::fprintf(stderr, "seams::Error %d %s\n", err.code(), err.what());
    return 1;
  }
  return 0;
}
