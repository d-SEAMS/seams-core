/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Shift-and-permute fuzz of the cage count. For every frame of every
** LAMMPS dump in a directory, the hexagonal and double-diamond cage counts
** are compared with the unshifted frame under
**
**   - every half-box translation (each axis 0 or 1/2),
**   - N random translations, rewrapped into the box,
**   - N random permutations of the atoms.
**
** A translation keeps atom indices, so the cage set must match, not only
** the count. A permutation changes the atom-index key that selects one
** maximum packing among several, so the count must match and the
** particular cages may differ.
**
**   fuzz_cage_invariance TRAJ_DIR [trials]
**
** Exits 0 when every frame agrees, 1 on a mismatch, 2 when TRAJ_DIR
** cannot be read.
*/

#include <cage.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
using Rings = std::vector<std::vector<int>>;

struct Result {
  int cages = 0;
  int hc = 0;
  int ddc = 0;
  int six = 0;
  std::vector<std::pair<int, std::vector<std::vector<int>>>> keys;
};

int majorityType(const std::string &path, int frame) {
  std::ifstream in(path);
  std::string line;
  int seen = 0;
  while (std::getline(in, line)) {
    if (line.rfind("ITEM: TIMESTEP", 0) != 0) {
      continue;
    }
    seen++;
    std::getline(in, line);
    std::getline(in, line);
    std::getline(in, line);
    const int n = std::stoi(line);
    if (seen != frame) {
      for (int i = 0; i < n + 5; i++) {
        std::getline(in, line);
      }
      continue;
    }
    for (int i = 0; i < 4; i++) {
      std::getline(in, line);
    }
    std::getline(in, line);
    std::istringstream hdr(line);
    std::string tok;
    hdr >> tok >> tok;
    int typeCol = -1;
    int col = 0;
    while (hdr >> tok) {
      if (tok == "type") {
        typeCol = col;
      }
      col++;
    }
    std::map<int, int> count;
    for (int i = 0; i < n; i++) {
      std::getline(in, line);
      if (typeCol < 0) {
        continue;
      }
      std::istringstream row(line);
      std::string field;
      for (int c = 0; c <= typeCol; c++) {
        row >> field;
      }
      count[std::stoi(field)]++;
    }
    int bestType = 1;
    int best = -1;
    for (const auto &kv : count) {
      if (kv.second > best) {
        best = kv.second;
        bestType = kv.first;
      }
    }
    return bestType;
  }
  return 1;
}

int frameCount(const std::string &path) {
  std::ifstream in(path);
  std::string line;
  int n = 0;
  while (std::getline(in, line)) {
    if (line.rfind("ITEM: TIMESTEP", 0) == 0) {
      n++;
    }
  }
  return n;
}

struct GraphSize {
  long edges = 0;
  int maxDegree = 0;
};

GraphSize graphSize(Cloud &cloud, double cutoff, int atomType) {
  GraphSize g;
  if (cloud.nop <= 0) {
    return g;
  }
  const auto nl = nneigh::neighListO(cutoff, cloud, atomType);
  const auto idx = nneigh::neighbourListByIndex(cloud, nl);
  for (const auto &row : idx) {
    const int d = static_cast<int>(row.size()) > 0 ? static_cast<int>(row.size()) - 1 : 0;
    if (d > g.maxDegree) {
      g.maxDegree = d;
    }
    if (d > 0) {
      g.edges += d;
    }
  }
  g.edges /= 2;
  return g;
}

Result evaluate(Cloud &cloud, double cutoff, int atomType) {
  Result out;
  if (cloud.nop <= 0) {
    return out;
  }
  const auto nl = nneigh::neighListO(cutoff, cloud, atomType);
  const auto idx = nneigh::neighbourListByIndex(cloud, nl);
  Rings six;
  for (const auto &r : primitive::ringNetwork(idx, 6)) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  out.six = static_cast<int>(six.size());
  if (six.empty()) {
    return out;
  }
  const auto index = ring::buildRingSearchIndex(six, cloud);
  std::vector<ring::strucType> types(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cages;
  const auto listHC = ring::findHC(six, types, idx, cages, index);
  (void)ring::findDDC(six, types, listHC, cages, index);
  out.cages = static_cast<int>(cages.size());
  for (const auto &cg : cages) {
    const bool hc = cg.type == cage::cageType::HexC;
    (hc ? out.hc : out.ddc)++;
    std::vector<std::vector<int>> faces;
    for (const int ri : cg.rings) {
      auto face = six[static_cast<size_t>(ri)];
      std::sort(face.begin(), face.end());
      faces.push_back(std::move(face));
    }
    std::sort(faces.begin(), faces.end());
    out.keys.emplace_back(hc ? 0 : 1, std::move(faces));
  }
  std::sort(out.keys.begin(), out.keys.end());
  return out;
}

void shiftCloud(Cloud &cloud, double fx, double fy, double fz) {
  // Lattice translation, then wrap in fractional coordinates. An
  // orthorhombic box reduces to an independent wrap of each axis. A
  // triclinic box does not: the dump spans are bounding edges, and a
  // Cartesian wrap of those spans is not a symmetry of the cell.
  const gen::FracBox b = gen::makeFracBox(cloud);
  if (!b.ok) {
    return;
  }
  const double add[3] = {fx, fy, fz};
  for (int i = 0; i < cloud.nop; i++) {
    auto &p = cloud.pts[static_cast<size_t>(i)];
    double sz = (p.z - b.oz) / b.lz;
    double sy = (p.y - b.oy - b.yz * sz) / b.ly;
    double sx = (p.x - b.ox - b.xy * sy - b.xz * sz) / b.lx;
    const double s[3] = {sx + add[0], sy + add[1], sz + add[2]};
    sx = s[0] - std::floor(s[0]);
    sy = s[1] - std::floor(s[1]);
    sz = s[2] - std::floor(s[2]);
    p.x = b.ox + b.lx * sx + b.xy * sy + b.xz * sz;
    p.y = b.oy + b.ly * sy + b.yz * sz;
    p.z = b.oz + b.lz * sz;
  }
}

Cloud permuteCloud(const Cloud &cloud, const std::vector<int> &perm) {
  Cloud out = cloud;
  const int n = cloud.nop;
  for (int a = 0; a < n; a++) {
    out.pts[static_cast<size_t>(perm[static_cast<size_t>(a)])] = cloud.pts[static_cast<size_t>(a)];
  }
  out.idIndexMap.clear();
  for (int i = 0; i < n; i++) {
    out.pts[static_cast<size_t>(i)].atomID = i + 1;
    out.idIndexMap[i + 1] = i;
  }
  return out;
}

bool sameCount(const Result &a, const Result &b) {
  return a.cages == b.cages && a.hc == b.hc && a.ddc == b.ddc && a.six == b.six;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "fuzz_cage_invariance TRAJ_DIR [trials]\n");
    return 2;
  }
  const fs::path dir(argv[1]);
  const int trials = argc > 2 ? std::atoi(argv[2]) : 1000;
  const double cutoff = 3.5;
  if (!fs::is_directory(dir)) {
    std::fprintf(stderr, "fuzz_cage_invariance: %s is not a directory\n", dir.c_str());
    return 2;
  }
  std::vector<fs::path> dumps;
  for (const auto &ent : fs::directory_iterator(dir)) {
    if (ent.path().extension() == ".lammpstrj") {
      dumps.push_back(ent.path());
    }
  }
  std::sort(dumps.begin(), dumps.end());
  if (dumps.empty()) {
    std::fprintf(stderr, "fuzz_cage_invariance: no dumps in %s\n", dir.c_str());
    return 2;
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(1);
#endif
  const unsigned nThreads = std::max(1u, std::thread::hardware_concurrency());

  int failures = 0;
  int frames = 0;
  std::uint64_t shiftTrials = 0;
  std::uint64_t permTrials = 0;
  for (const auto &path : dumps) {
    const int nFrames = frameCount(path.string());
    for (int frame = 1; frame <= nFrames; frame++) {
      frames++;
      const int atomType = majorityType(path.string(), frame);
      Cloud cloud;
      cloud = sinp::readLammpsTrjO(path.string(), frame, cloud, atomType);
      if (cloud.nop <= 0) {
        std::printf("%s frame %d type %d EMPTY\n", path.filename().c_str(), frame, atomType);
        continue;
      }
      const GraphSize g0 = graphSize(cloud, cutoff, atomType);
      // A cutoff graph this dense is not an ice cage frame: ring enumeration
      // either does not return within half a minute or emits on the order of
      // a million six-rings. 1000 enumerations are not a cage-count test.
      if (g0.maxDegree > 40) {
        std::printf("%s frame %d type %d nop %d maxdeg %d edges %ld SKIP dense cutoff graph\n",
                    path.filename().c_str(), frame, atomType, cloud.nop, g0.maxDegree, g0.edges);
        std::fflush(stdout);
        continue;
      }
      const Result ref = evaluate(cloud, cutoff, atomType);
      bool ok = true;
      for (int ix = 0; ix < 2 && ok; ix++) {
        for (int iy = 0; iy < 2 && ok; iy++) {
          for (int iz = 0; iz < 2 && ok; iz++) {
            Cloud moved = cloud;
            shiftCloud(moved, 0.5 * ix, 0.5 * iy, 0.5 * iz);
            const Result got = evaluate(moved, cutoff, atomType);
            shiftTrials++;
            if (!sameCount(ref, got) || got.keys != ref.keys) {
              std::printf("FAIL %s frame %d half %d%d%d cages %d/%d (hc %d/%d ddc %d/%d) six %d/%d set %s\n",
                          path.filename().c_str(), frame, ix, iy, iz, ref.cages, got.cages, ref.hc,
                          got.hc, ref.ddc, got.ddc, ref.six, got.six,
                          got.keys == ref.keys ? "same" : "DIFFERS");
              ok = false;
            }
          }
        }
      }
      std::atomic<int> bad{0};
      auto runRange = [&](int begin, int end) {
        for (int t = begin; t < end && bad.load() == 0; t++) {
          std::mt19937 rng(20261101u + static_cast<unsigned>(frames) * 9973u +
                           static_cast<unsigned>(t) * 10007u);
          std::uniform_real_distribution<double> unit(0.0, 1.0);
          Cloud moved = cloud;
          shiftCloud(moved, unit(rng), unit(rng), unit(rng));
          const Result got = evaluate(moved, cutoff, atomType);
          if (!sameCount(ref, got) || got.keys != ref.keys) {
            std::printf("FAIL %s frame %d shift trial %d cages %d/%d (hc %d/%d ddc %d/%d)\n",
                        path.filename().c_str(), frame, t, ref.cages, got.cages, ref.hc, got.hc,
                        ref.ddc, got.ddc);
            bad.store(1);
            return;
          }
          const int n = cloud.nop;
          std::vector<int> perm(static_cast<size_t>(n));
          std::iota(perm.begin(), perm.end(), 0);
          std::shuffle(perm.begin(), perm.end(), rng);
          Cloud permuted = permuteCloud(cloud, perm);
          const Result pgot = evaluate(permuted, cutoff, atomType);
          if (!sameCount(ref, pgot)) {
            std::printf("FAIL %s frame %d perm trial %d cages %d/%d (hc %d/%d ddc %d/%d) six %d/%d\n",
                        path.filename().c_str(), frame, t, ref.cages, pgot.cages, ref.hc, pgot.hc,
                        ref.ddc, pgot.ddc, ref.six, pgot.six);
            bad.store(1);
            return;
          }
        }
      };
      std::vector<std::thread> pool;
      const int chunk = std::max(1, (trials + static_cast<int>(nThreads) - 1) / static_cast<int>(nThreads));
      for (unsigned w = 0; w < nThreads; w++) {
        const int begin = static_cast<int>(w) * chunk;
        const int end = std::min(trials, begin + chunk);
        if (begin >= end) {
          break;
        }
        pool.emplace_back(runRange, begin, end);
      }
      for (auto &th : pool) {
        th.join();
      }
      shiftTrials += static_cast<std::uint64_t>(trials);
      permTrials += static_cast<std::uint64_t>(trials);
      if (bad.load() != 0) {
        ok = false;
      }
      if (!ok) {
        failures++;
      }
      std::printf("%s frame %d type %d nop %d six %d cages %d (hc %d ddc %d) %s\n",
                  path.filename().c_str(), frame, atomType, cloud.nop, ref.six, ref.cages, ref.hc,
                  ref.ddc, ok ? "OK" : "FAIL");
      std::fflush(stdout);
    }
  }
  std::printf("frames %d failures %d shiftTrials %llu permTrials %llu\n", frames, failures,
              static_cast<unsigned long long>(shiftTrials),
              static_cast<unsigned long long>(permTrials));
  return failures == 0 ? 0 : 1;
}
