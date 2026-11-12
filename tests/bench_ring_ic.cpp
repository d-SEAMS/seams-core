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
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <sstream>
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

molSys::PointCloud<molSys::Point<double>, double> asCloud(const Cloud &ice) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {ice.side, ice.side, ice.side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = static_cast<int>(ice.xyz.size());
  cloud.pts.resize(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; i++) {
    auto &p = cloud.pts[static_cast<std::size_t>(i)];
    p.x = ice.xyz[static_cast<std::size_t>(i)][0];
    p.y = ice.xyz[static_cast<std::size_t>(i)][1];
    p.z = ice.xyz[static_cast<std::size_t>(i)][2];
    p.type = 1;
    p.atomID = i;
    p.molID = i;
    cloud.idIndexMap[i] = i;
  }
  return cloud;
}

int main(int argc, char **argv) {
  std::vector<int> reps;
  int maxDepth = kMaxDepth;
  bool e2e = false;
  int onlyThreads = 0;
  std::string dumpPath;
  std::vector<int> frames;
  bool listRings = false;
  bool listEdges = false;
  for (int i = 1; i < argc; i++) {
    const std::string arg = argv[i];
    if (arg == "--depth" && i + 1 < argc) {
      maxDepth = std::atoi(argv[++i]);
    } else if (arg == "--e2e") {
      e2e = true;
    } else if (arg == "--threads" && i + 1 < argc) {
      onlyThreads = std::atoi(argv[++i]);
    } else if (arg == "--dump" && i + 1 < argc) {
      dumpPath = argv[++i];
    } else if (arg == "--frame" && i + 1 < argc) {
      frames.push_back(std::atoi(argv[++i]));
    } else if (arg == "--list") {
      listRings = true;
    } else if (arg == "--edges") {
      listEdges = true;
    } else {
      reps.push_back(std::atoi(argv[i]));
    }
  }
  if (reps.empty()) {
    reps = {11, 24};
  }
  const int threads[] = {1, 8};
  constexpr int kWarm = 1;
  constexpr int kRepeat = 8;
  std::cout << std::setprecision(6) << std::fixed;
  bool allEqual = true;
  if (!dumpPath.empty()) {
    std::ifstream in(dumpPath);
    if (!in) {
      std::cerr << "cannot read " << dumpPath << "\n";
      return 1;
    }
    std::string line;
    int frame = 0;
    while (std::getline(in, line)) {
      if (line.find("NUMBER OF ATOMS") == std::string::npos) {
        continue;
      }
      int n = 0;
      in >> n;
      std::getline(in, line);
      frame++;
      const bool wanted = frames.empty() ||
                          std::find(frames.begin(), frames.end(), frame) != frames.end();
      double lo[3] = {}, hi[3] = {};
      for (int axis = 0; axis < 3; axis++) {
        std::getline(in, line);
        if (line.find("ITEM:") != std::string::npos) {
          axis--;
          continue;
        }
        std::istringstream box(line);
        box >> lo[axis] >> hi[axis];
      }
      std::getline(in, line); // ITEM: ATOMS ...
      std::vector<std::string> columns;
      {
        std::istringstream hdr(line);
        std::string word;
        hdr >> word;
        hdr >> word;
        while (hdr >> word) {
          columns.push_back(word);
        }
      }
      const auto column = [&](const char *name) {
        for (int c = 0; c < static_cast<int>(columns.size()); c++) {
          if (columns[static_cast<std::size_t>(c)] == name) {
            return c;
          }
        }
        return -1;
      };
      const int cx = column("x");
      const int cy = column("y");
      const int cz = column("z");
      const int cid = column("id");
      molSys::PointCloud<molSys::Point<double>, double> cloud;
      cloud.box = {hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
      cloud.boxLow = {lo[0], lo[1], lo[2]};
      cloud.nop = n;
      cloud.pts.resize(static_cast<std::size_t>(n));
      for (int i = 0; i < n; i++) {
        std::getline(in, line);
        if (!wanted) {
          continue;
        }
        std::istringstream row(line);
        std::vector<std::string> tok;
        std::string word;
        while (row >> word) {
          tok.push_back(word);
        }
        int id = i + 1;
        double x = 0, y = 0, z = 0;
        if (cx >= 0 && cy >= 0 && cz >= 0 &&
            cz < static_cast<int>(tok.size())) {
          x = std::atof(tok[static_cast<std::size_t>(cx)].c_str());
          y = std::atof(tok[static_cast<std::size_t>(cy)].c_str());
          z = std::atof(tok[static_cast<std::size_t>(cz)].c_str());
          if (cid >= 0) {
            id = std::atoi(tok[static_cast<std::size_t>(cid)].c_str());
          }
        } else if (tok.size() >= 4) {
          id = std::atoi(tok[0].c_str());
          x = std::atof(tok[tok.size() - 3].c_str());
          y = std::atof(tok[tok.size() - 2].c_str());
          z = std::atof(tok.back().c_str());
        } else if (tok.size() == 3) {
          x = std::atof(tok[0].c_str());
          y = std::atof(tok[1].c_str());
          z = std::atof(tok[2].c_str());
        }
        auto &p = cloud.pts[static_cast<std::size_t>(i)];
        p.atomID = id;
        p.molID = id;
        p.type = 1;
        p.x = x;
        p.y = y;
        p.z = z;
        cloud.idIndexMap[id] = i;
      }
      if (!wanted) {
        continue;
      }
      setThreads(onlyThreads > 0 ? onlyThreads : 1);
      const auto t0 = std::chrono::steady_clock::now();
      std::vector<int> offsets, cols;
      const bool ok = nneigh::orthoFullCsr(cloud, nullptr, kCutoff, offsets, cols);
      const auto rings = primitive::ringNetwork(offsets, cols, maxDepth);
      const double total = millis(t0);
      std::vector<int> bySize(static_cast<std::size_t>(maxDepth) + 1, 0);
      for (const auto &ring : rings) {
        if (ring.size() < bySize.size()) {
          bySize[ring.size()]++;
        }
      }
      std::cout << "dump frame " << frame << " molecules " << n << " depth " << maxDepth
                << " csr " << (ok ? 1 : 0) << " total_ms " << total << " rings "
                << rings.size();
      for (int s = 3; s <= maxDepth; s++) {
        std::cout << " n" << s << " " << bySize[static_cast<std::size_t>(s)];
      }
      long edges = 0;
      if (offsets.size() >= 2) {
        edges = offsets.back() / 2;
      }
      std::cout << " edges " << edges << "\n";
      if (listEdges) {
        const int n = static_cast<int>(offsets.size()) - 1;
        for (int i = 0; i < n; i++) {
          for (int k = offsets[static_cast<std::size_t>(i)];
               k < offsets[static_cast<std::size_t>(i) + 1]; k++) {
            const int j = cols[static_cast<std::size_t>(k)];
            if (i < j) {
              std::cout << i << " " << j << "\n";
            }
          }
        }
      }
      if (listRings) {
        for (const auto &ring : rings) {
          for (std::size_t k = 0; k < ring.size(); k++) {
            if (k) {
              std::cout << " ";
            }
            std::cout << ring[k];
          }
          std::cout << "\n";
        }
      }
    }
    return 0;
  }
  if (e2e) {
    for (const int rep : reps) {
      const Cloud ice = iceIc(rep);
      auto cloud = asCloud(ice);
      const auto nList = bondsOf(ice);
      const auto reference = primitive::ringNetwork(nList, maxDepth);
      const int n = cloud.nop;
      const int threads[] = {1, 8};
      for (const int thr : threads) {
        if (onlyThreads != 0 && thr != onlyThreads) {
          continue;
        }
        setThreads(thr);
        std::vector<double> totalMs;
        std::vector<std::vector<int>> rings;
        for (int sample = 0; sample < kWarm + kRepeat; sample++) {
          const auto t0 = std::chrono::steady_clock::now();
          std::vector<int> offsets;
          std::vector<int> cols;
          const bool ok = nneigh::orthoFullCsr(cloud, nullptr, kCutoff, offsets, cols);
          const auto t1 = std::chrono::steady_clock::now();
          rings = ok ? primitive::ringNetwork(offsets, cols, maxDepth)
                     : primitive::ringNetwork(nList, maxDepth);
          const auto t2 = std::chrono::steady_clock::now();
          if (sample >= kWarm) {
            const double neigh = std::chrono::duration<double, std::milli>(t1 - t0).count();
            const double ring = std::chrono::duration<double, std::milli>(t2 - t1).count();
            totalMs.push_back(neigh + ring);
            const auto split = primitive::lastRingSearchSplit();
            std::cout << "sample e2e threads " << threadCount() << " molecules " << n
                      << " depth " << maxDepth << " neigh_ms " << neigh << " ring_ms " << ring
                      << " total_ms " << neigh + ring << " rings " << rings.size()
                      << " csr " << (ok ? 1 : 0) << " balls_ms " << split.balls_ms
                      << " enum_ms " << split.enumerate_ms << " merge_ms " << split.merge_ms
                      << "\n";
          }
        }
        const bool equal = rings == reference;
        allEqual = allEqual && equal;
        std::sort(totalMs.begin(), totalMs.end());
        const double med = totalMs[totalMs.size() / 2];
        std::cout << "summary e2e threads " << threadCount() << " molecules " << n
                  << " depth " << maxDepth << " total_median_ms " << med << " rings "
                  << rings.size() << " equal " << (equal ? 1 : 0) << " used_balls "
                  << primitive::lastRingSearchUsedBalls() << "\n";
      }
    }
    return allEqual ? 0 : 1;
  }
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
