/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Measures what d-SEAMS adds on top of the vesin cell list.  The cell list
** itself is a third-party kernel; the question this answers is how much of
** the neighbour-search wall time is spent in the kernel against the wrapper
** around it (type filtering, coordinate marshalling, identifier mapping and
** building the row-ordered output).
**
** Build target: bench_overhead.  Run: bcpp/tests/bench_overhead [nAtoms ...]
*/

#include <mol_sys.hpp>
#include <neighbours.hpp>

#ifdef SEAMS_HAS_VESIN
#include <vesin.h>
#endif

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

constexpr double kNumberDensity = 0.0332;

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
    point.z = (((i / (perSide * perSide)) % perSide) + 0.5) * spacing + jitter();
    cloud.pts.push_back(point);
    cloud.idIndexMap[point.atomID] = i;
  }
  return cloud;
}

template <typename Fn> double bestMillis(Fn &&fn, int reps) {
  double best = std::numeric_limits<double>::max();
  for (int rep = 0; rep < reps; rep++) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    best = std::min(
        best, std::chrono::duration<double, std::milli>(end - start).count());
  }
  return best;
}

} // namespace

int main(int argc, char **argv) {
#ifndef SEAMS_HAS_VESIN
  std::cout << "built without vesin; nothing to compare\n";
  return 0;
#else
  std::vector<int> sizes;
  for (int arg = 1; arg < argc; arg++) {
    sizes.push_back(std::atoi(argv[arg]));
  }
  if (sizes.empty()) {
    sizes = {4000, 16000, 64000, 256000};
  }

  constexpr int typeI = 1;
  constexpr double cutoff = 3.5;

  std::cout << std::left << std::setw(10) << "nAtoms" << std::setw(14)
            << "vesin/ms" << std::setw(14) << "neighListO/ms" << std::setw(12)
            << "overhead" << std::setw(16) << "Matoms/s" << std::setw(14)
            << "pairs/atom" << "\n";
  std::cout << std::string(80, '-') << "\n";

  for (int nAtoms : sizes) {
    Cloud cloud = makeCloud(nAtoms, typeI);
    const int reps = nAtoms <= 16000 ? 5 : 3;

    // Coordinates in the layout vesin wants, marshalled once and excluded
    // from the timing, so this measures the kernel alone
    std::vector<std::array<double, 3>> positions(nAtoms);
    for (int i = 0; i < nAtoms; i++) {
      positions[i] = {cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z};
    }
    double box[3][3] = {{cloud.box[0], 0.0, 0.0},
                        {0.0, cloud.box[1], 0.0},
                        {0.0, 0.0, cloud.box[2]}};
    bool periodic[3] = {true, true, true};

    VesinOptions options;
    options.cutoff = cutoff;
    options.full = true;
    options.sorted = false;
    options.algorithm = VesinAutoAlgorithm;
    options.return_shifts = false;
    options.return_distances = false;
    options.return_vectors = false;
    VesinDevice device = {VesinCPU, 0};

    size_t nPairs = 0;
    const double tVesin = bestMillis(
        [&]() {
          VesinNeighborList neighbors;
          const char *err = nullptr;
          vesin_neighbors(
              reinterpret_cast<const double (*)[3]>(positions.data()), nAtoms,
              box, periodic, device, options, &neighbors, &err);
          nPairs = neighbors.length;
          vesin_free(&neighbors);
        },
        reps);

    const double tFull = bestMillis(
        [&]() {
          volatile auto nList = nneigh::neighListO(cutoff, cloud, typeI);
          (void)nList;
        },
        reps);

    const double throughput = nAtoms / (tFull * 1000.0); // atoms per microsecond
    std::cout << std::left << std::setw(10) << nAtoms << std::setw(14)
              << std::fixed << std::setprecision(3) << tVesin << std::setw(14)
              << tFull << std::setw(12) << std::setprecision(2)
              << (tFull / tVesin) << std::setw(16) << std::setprecision(2)
              << throughput << std::setw(14) << std::setprecision(1)
              << (static_cast<double>(nPairs) / nAtoms) << "\n";
  }

  return 0;
#endif
}
