/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Scaling benchmarks for the neighbour-list and bond-order-parameter hot
** paths.  Synthetic point clouds are generated at a fixed number density so
** that the mean coordination number stays constant across system sizes and
** the measured curve reflects the algorithm rather than the neighbourhood
** size.
**
** Build target: bench_scaling (not in the regular test suite).
** Run: bcpp/tests/bench_scaling [nAtoms ...]
*/

#include <bop.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

/// Number density of the mW water systems shipped with d-SEAMS, in atoms per
/// cubic Angstrom.  Holding this fixed keeps the coordination number, and so
/// the neighbour-list output size, proportional to the atom count.
constexpr double kNumberDensity = 0.0332;

/// Builds a cubic point cloud of @a nAtoms particles on a jittered lattice.
Cloud makeCloud(int nAtoms, int typeI) {
  Cloud cloud;
  const double boxLength = std::cbrt(nAtoms / kNumberDensity);
  const int perSide = static_cast<int>(std::ceil(std::cbrt(nAtoms)));
  const double spacing = boxLength / perSide;

  cloud.nop = nAtoms;
  cloud.currentFrame = 1;
  cloud.box = {boxLength, boxLength, boxLength};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.pts.reserve(nAtoms);

  // Deterministic jitter, so that repeated runs compare like with like
  unsigned long long state = 88172645463325252ULL;
  auto jitter = [&state, spacing]() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    const double unit = static_cast<double>(state >> 11) / 9007199254740992.0;
    return (unit - 0.5) * 0.3 * spacing;
  };

  for (int i = 0; i < nAtoms; i++) {
    molSys::Point<double> point;
    point.type = typeI;
    point.atomID = i + 1;
    point.molID = i + 1;
    point.x = ((i % perSide) + 0.5) * spacing + jitter();
    point.y = (((i / perSide) % perSide) + 0.5) * spacing + jitter();
    point.z = (((i / (perSide * perSide)) % perSide) + 0.5) * spacing +
              jitter();
    cloud.pts.push_back(point);
    cloud.idIndexMap[point.atomID] = i;
  }

  return cloud;
}

/// Runs @a fn @a reps times and returns the best wall time in milliseconds.
/// The minimum is reported rather than the mean because a loaded machine adds
/// only positive noise.
template <typename Fn> double bestMillis(Fn &&fn, int reps) {
  double best = std::numeric_limits<double>::max();
  for (int rep = 0; rep < reps; rep++) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    const double ms =
        std::chrono::duration<double, std::milli>(end - start).count();
    best = std::min(best, ms);
  }
  return best;
}

} // namespace

int main(int argc, char **argv) {
  std::vector<int> sizes;
  for (int arg = 1; arg < argc; arg++) {
    sizes.push_back(std::atoi(argv[arg]));
  }
  if (sizes.empty()) {
    sizes = {1000, 2000, 4000, 8000, 16000, 32000};
  }

  constexpr int typeI = 1;
  constexpr double cutoff = 3.5;

  std::cout << std::left << std::setw(10) << "nAtoms" << std::setw(16)
            << "neighListO/ms" << std::setw(16) << "byIndex/ms"
            << std::setw(16) << "getCorrelPlus/ms" << std::setw(16)
            << "ringNetwork/ms" << "\n";
  std::cout << std::string(74, '-') << "\n";

  for (int nAtoms : sizes) {
    Cloud cloud = makeCloud(nAtoms, typeI);
    const int reps = nAtoms <= 4000 ? 5 : 3;

    std::vector<std::vector<int>> nList;
    const double tNeigh = bestMillis(
        [&]() { nList = nneigh::neighListO(cutoff, cloud, typeI); }, reps);

    std::vector<std::vector<int>> byIndex;
    const double tIndex = bestMillis(
        [&]() { byIndex = nneigh::getNewNeighbourListByIndex(cloud, cutoff); },
        reps);

    const double tCorrel = bestMillis(
        [&]() { chill::getCorrelPlus(cloud, nList, false); }, reps);

    const double tRings = bestMillis(
        [&]() {
          volatile auto rings = primitive::ringNetwork(byIndex, 6);
          (void)rings;
        },
        nAtoms <= 8000 ? reps : 1);

    std::cout << std::left << std::setw(10) << nAtoms << std::setw(16)
              << std::fixed << std::setprecision(3) << tNeigh << std::setw(16)
              << tIndex << std::setw(16) << tCorrel << std::setw(16) << tRings
              << "\n";
  }

  return 0;
}
