//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// One session, release build. Ideal ice Ic at 10648 and 110592 molecules,
// rings of size at most 6, the hop probe against the bounded-ball search,
// one thread and eight. Coordinates are the oxygen sites of cubic ice:
// lattice 4 * 2.75 / sqrt(3) angstrom, eight molecules per conventional
// cell. A 3.5 angstrom cutoff keeps the four hydrogen-bond neighbours and
// drops the second shell. Nothing is written except this log.

#include <franzblau.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

constexpr double kBond = 2.75;
constexpr double kCutoff = 3.5;
constexpr int kMaxDepth = 6;

struct Cloud {
  std::vector<std::array<double, 3>> xyz;
  double side = 0;
};

Cloud iceIc(int repetitions) {
  const double lattice = 4.0 * kBond / std::sqrt(3.0);
  const double basis[8][3] = {
      {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5},
      {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25}, {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  Cloud cloud;
  cloud.side = repetitions * lattice;
  cloud.xyz.reserve(static_cast<std::size_t>(8) * repetitions * repetitions * repetitions);
  for (int ix = 0; ix < repetitions; ix++) {
    for (int iy = 0; iy < repetitions; iy++) {
      for (int iz = 0; iz < repetitions; iz++) {
        for (const auto &b : basis) {
          cloud.xyz.push_back({(ix + b[0]) * lattice, (iy + b[1]) * lattice,
                               (iz + b[2]) * lattice});
        }
      }
    }
  }
  return cloud;
}

std::vector<std::vector<int>> bondsOf(const Cloud &cloud) {
  const int n = static_cast<int>(cloud.xyz.size());
  const double side = cloud.side;
  const int nCell = std::max(1, static_cast<int>(side / kCutoff));
  const double cell = side / nCell;
  std::vector<std::vector<int>> cells(static_cast<std::size_t>(nCell) * nCell * nCell);
  auto cellOf = [&](double x) {
    int c = static_cast<int>(x / cell);
    if (c < 0) {
      c = 0;
    }
    if (c >= nCell) {
      c = nCell - 1;
    }
    return c;
  };
  for (int i = 0; i < n; i++) {
    const int cx = cellOf(cloud.xyz[static_cast<std::size_t>(i)][0]);
    const int cy = cellOf(cloud.xyz[static_cast<std::size_t>(i)][1]);
    const int cz = cellOf(cloud.xyz[static_cast<std::size_t>(i)][2]);
    cells[(cx * nCell + cy) * nCell + cz].push_back(i);
  }
  const double cut2 = kCutoff * kCutoff;
  std::vector<std::vector<int>> nList(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    const auto &pi = cloud.xyz[static_cast<std::size_t>(i)];
    const int cx = cellOf(pi[0]);
    const int cy = cellOf(pi[1]);
    const int cz = cellOf(pi[2]);
    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        for (int dz = -1; dz <= 1; dz++) {
          const int nx = (cx + dx + nCell) % nCell;
          const int ny = (cy + dy + nCell) % nCell;
          const int nz = (cz + dz + nCell) % nCell;
          for (const int j : cells[(nx * nCell + ny) * nCell + nz]) {
            if (j <= i) {
              continue;
            }
            const auto &pj = cloud.xyz[static_cast<std::size_t>(j)];
            double r2 = 0;
            for (int a = 0; a < 3; a++) {
              double d = pi[a] - pj[a];
              if (d > side / 2) {
                d -= side;
              } else if (d < -side / 2) {
                d += side;
              }
              r2 += d * d;
            }
            if (r2 <= cut2 && r2 > 1e-12) {
              nList[static_cast<std::size_t>(i)].push_back(j);
              nList[static_cast<std::size_t>(j)].push_back(i);
            }
          }
        }
      }
    }
  }
  // The ring search treats the first entry of each row as the vertex
  // itself. Neighbours are appended before that prefix exists, so the
  // vertex is inserted once every pair has been recorded.
  for (int i = 0; i < n; i++) {
    nList[static_cast<std::size_t>(i)].insert(nList[static_cast<std::size_t>(i)].begin(), i);
  }
  return nList;
}

void setThreads(int n) {
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(n);
#else
  (void)n;
#endif
}

int threadCount() {
#ifdef SEAMS_HAS_OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

double millis(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

} // namespace

int main(int argc, char **argv) {
  std::vector<int> reps;
  int maxDepth = kMaxDepth;
  for (int i = 1; i < argc; i++) {
    const std::string arg = argv[i];
    if (arg == "--depth" && i + 1 < argc) {
      maxDepth = std::atoi(argv[++i]);
    } else {
      reps.push_back(std::atoi(argv[i]));
    }
  }
  if (reps.empty()) {
    reps = {11, 24};
  }
  const int threads[] = {1, 8};
  constexpr int kWarm = 1;
  constexpr int kRepeat = 6;
  std::cout << std::setprecision(6) << std::fixed;
  bool allEqual = true;
  for (const int rep : reps) {
    const Cloud cloud = iceIc(rep);
    const auto nList = bondsOf(cloud);
    const int n = static_cast<int>(cloud.xyz.size());
    int degree = 0;
    int minDegree = n;
    long edges = 0;
    for (const auto &row : nList) {
      const int d = static_cast<int>(row.size()) - 1;
      degree = std::max(degree, d);
      minDegree = std::min(minDegree, d);
      edges += d;
    }
    std::cout << "frame molecules " << n << " depth " << maxDepth << " side " << cloud.side
              << " min_degree "
              << minDegree << " max_degree " << degree << " edges " << edges / 2 << "\n";
    for (const int thr : threads) {
      setThreads(thr);
      std::vector<double> legacyMs;
      std::vector<double> hopMs;
      std::vector<std::vector<int>> legacyRings;
      std::vector<std::vector<int>> hopRings;
      for (int sample = 0; sample < kWarm + kRepeat; sample++) {
        const auto t0 = std::chrono::steady_clock::now();
        legacyRings = primitive::ringNetworkBalls(nList, maxDepth);
        const double legacy = millis(t0);
        const auto t1 = std::chrono::steady_clock::now();
        hopRings = primitive::ringNetwork(nList, maxDepth);
        const double hop = millis(t1);
        if (sample >= kWarm) {
          legacyMs.push_back(legacy);
          hopMs.push_back(hop);
          std::cout << "sample threads " << threadCount() << " molecules " << n
                    << " legacy_ms " << legacy << " hop_ms " << hop << " legacy_rings "
                    << legacyRings.size() << " hop_rings " << hopRings.size() << "\n";
        }
      }
      const bool equal = legacyRings == hopRings;
      allEqual = allEqual && equal;
      auto median = [](std::vector<double> v) {
        std::sort(v.begin(), v.end());
        const std::size_t m = v.size() / 2;
        if (v.size() % 2 == 0) {
          return 0.5 * (v[m - 1] + v[m]);
        }
        return v[m];
      };
      std::cout << "summary threads " << threadCount() << " molecules " << n
                << " legacy_median_ms " << median(legacyMs) << " hop_median_ms "
                << median(hopMs) << " rings " << hopRings.size() << " equal "
                << (equal ? 1 : 0) << " used_balls " << primitive::lastRingSearchUsedBalls()
                << "\n";
    }
  }
  return allEqual ? 0 : 1;
}
