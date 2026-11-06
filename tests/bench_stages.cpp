/*
** Stage times for one frame: I/O, neighbours, rings, affiliation.
**   bench_stages TRAJ [frame] [atomType] [reps]
*/

#include <cage_affiliation.hpp>
#include <edge_set.hpp>
#include <franzblau.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
using seams::EdgeSet;
using seams::edgeSetOf;
using seams::provenEqual;

std::vector<std::vector<int>> onlySix(const std::vector<std::vector<int>> &rings) {
  std::vector<std::vector<int>> six;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  return six;
}

template <typename Fn> double bestMs(Fn &&fn, int reps) {
  double best = std::numeric_limits<double>::max();
  for (int i = 0; i < reps; i++) {
    const auto t0 = Clock::now();
    fn();
    const auto t1 = Clock::now();
    const double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    if (ms < best) {
      best = ms;
    }
  }
  return best;
}
} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr,
                 "usage: bench_stages TRAJ [frame] [atomType] [reps]\n");
    return 2;
  }
  const std::string traj = argv[1];
  const int frame = argc > 2 ? std::atoi(argv[2]) : 1;
  const int typeI = argc > 3 ? std::atoi(argv[3]) : 1;
  const int reps = argc > 4 ? std::atoi(argv[4]) : 3;

  molSys::PointCloud<molSys::Point<double>, double> cloud;
  if (std::getenv("SEAMS_TIME_READ") != nullptr) {
    const auto t0 = Clock::now();
    cloud = sinp::readLammpsTrjO(traj, frame, cloud, typeI);
    const auto t1 = Clock::now();
    if (cloud.nop == 0) {
      std::fprintf(stderr, "empty frame %d of %s\n", frame, traj.c_str());
      return 1;
    }
    const double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::printf("# traj %s frame %d nop %d rings 0 six 0 ddc 0\n", traj.c_str(),
                frame, cloud.nop);
    std::printf("io_ms %.3f\nneigh_ms 0\nknn_ms 0\nknn_pair_ms 0\n"
                "rings_ms 0\naffil_ms 0\nseeded_ms 0\n",
                ms);
    return 0;
  }

  const double tIo = bestMs(
      [&]() { cloud = sinp::readLammpsTrjO(traj, frame, cloud, typeI); },
      reps);
  if (cloud.nop == 0) {
    std::fprintf(stderr, "empty frame %d of %s\n", frame, traj.c_str());
    return 1;
  }

  std::vector<std::vector<int>> nList, idx, rings, six;
  const double tNeigh = bestMs(
      [&]() {
        nList = nneigh::neighListO(3.5, cloud, typeI);
        idx = nneigh::neighbourListByIndex(cloud, nList);
      },
      reps);
  const double tKnn = bestMs(
      [&]() {
        auto knn = nneigh::kNearestNeighbourList(cloud, 4, 5.5, typeI, true);
        (void)knn;
      },
      reps);
  const double tKnnPair = bestMs(
      [&]() {
        auto both = nneigh::kNearestNeighbourPair(cloud, 4, 5.5, typeI);
        (void)both;
      },
      reps);
  const double tRings = bestMs(
      [&]() { rings = primitive::ringNetwork(idx, 6); }, reps);
  six.clear();
  for (const auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  ring::CageAffiliation aff;
  const double tAff = bestMs(
      [&]() { aff = ring::cageAffiliation(six, idx); }, reps);
  auto knnPair = nneigh::kNearestNeighbourPair(cloud, 4, 5.5, typeI);
  const auto &knn = knnPair.first;
  const auto &uni = knnPair.second;
  auto idxS = nneigh::neighbourListByIndex(cloud, knn);
  auto idxU = nneigh::neighbourListByIndex(cloud, uni);
  int edgesSu = 0;
  int edgesCutFlag = 0;
  const double tSeeded = bestMs(
      [&]() {
        const EdgeSet edgesCut = edgeSetOf(idx);
        const EdgeSet edgesS = edgeSetOf(idxS);
        const EdgeSet edgesU = edgeSetOf(idxU);
        // Hash first. A match is checked again by comparing every edge.
        const bool strictIsPermissive = provenEqual(edgesS, edgesU);
        const bool strictIsCutoff = provenEqual(edgesS, edgesCut);
        edgesSu = strictIsPermissive ? 1 : 0;
        edgesCutFlag = strictIsCutoff ? 1 : 0;
        if (strictIsPermissive && strictIsCutoff) {
          (void)ring::seededCageAffiliation(six, aff, six, aff, idx);
          return;
        }
        if (strictIsPermissive) {
          const auto s6 = onlySix(primitive::ringNetwork(idxS, 6));
          const auto once = ring::cageAffiliation(s6, idxS);
          (void)ring::seededCageAffiliation(s6, once, s6, once, idxS);
          return;
        }
        const auto s6 = onlySix(primitive::ringNetwork(idxS, 6));
        const auto u6 = onlySix(primitive::ringNetwork(idxU, 6));
        (void)ring::seededCageAffiliation(s6, idxS, u6, idxU);
      },
      reps);
  std::printf("edges_su %d\nedges_cut %d\n", edgesSu, edgesCutFlag);

  if (const char *dumpPath = std::getenv("SEAMS_DUMP")) {
    const auto seeded = [&]() {
      const EdgeSet cut = edgeSetOf(idx);
      const EdgeSet es = edgeSetOf(idxS);
      const EdgeSet eu = edgeSetOf(idxU);
      if (provenEqual(es, eu) && provenEqual(es, cut)) {
        return ring::seededCageAffiliation(six, aff, six, aff, idx);
      }
      if (provenEqual(es, eu)) {
        const auto s6 = onlySix(primitive::ringNetwork(idxS, 6));
        const auto once = ring::cageAffiliation(s6, idxS);
        return ring::seededCageAffiliation(s6, once, s6, once, idxS);
      }
      const auto s6 = onlySix(primitive::ringNetwork(idxS, 6));
      const auto u6 = onlySix(primitive::ringNetwork(idxU, 6));
      return ring::seededCageAffiliation(s6, idxS, u6, idxU);
    }();
    std::vector<char> molHc(static_cast<std::size_t>(cloud.nop), 0);
    std::vector<char> molDdc(static_cast<std::size_t>(cloud.nop), 0);
    for (std::size_t r = 0; r < six.size(); ++r) {
      for (const int a : six[r]) {
        if (a < 0 || a >= cloud.nop) {
          continue;
        }
        if (r < aff.hc.size() && aff.hc[r]) {
          molHc[static_cast<std::size_t>(a)] = 1;
        }
        if (r < aff.ddc.size() && aff.ddc[r]) {
          molDdc[static_cast<std::size_t>(a)] = 1;
        }
      }
    }
    auto canon = [](std::vector<int> ring) {
      if (ring.empty()) {
        return ring;
      }
      const auto least = std::min_element(ring.begin(), ring.end());
      std::rotate(ring.begin(), least, ring.end());
      if (ring.size() > 1 && ring.back() < ring[1]) {
        std::reverse(ring.begin() + 1, ring.end());
      }
      return ring;
    };
    std::vector<std::vector<int>> canonRings;
    canonRings.reserve(rings.size());
    for (const auto &ring : rings) {
      canonRings.push_back(canon(ring));
    }
    std::sort(canonRings.begin(), canonRings.end());
    std::FILE *fh = std::fopen(dumpPath, "w");
    if (fh == nullptr) {
      std::perror(dumpPath);
      return 1;
    }
    std::fprintf(fh, "rings %zu\n", canonRings.size());
    for (const auto &ring : canonRings) {
      for (std::size_t i = 0; i < ring.size(); ++i) {
        std::fprintf(fh, "%s%d", i ? " " : "", ring[i]);
      }
      std::fprintf(fh, "\n");
    }
    std::fprintf(fh, "molecules %d\n", cloud.nop);
    for (int a = 0; a < cloud.nop; ++a) {
      const int sh = a < static_cast<int>(seeded.hc.size()) && seeded.hc[static_cast<std::size_t>(a)];
      const int sd = a < static_cast<int>(seeded.ddc.size()) && seeded.ddc[static_cast<std::size_t>(a)];
      std::fprintf(fh, "%d %d %d %d %d\n", a, molHc[static_cast<std::size_t>(a)],
                   molDdc[static_cast<std::size_t>(a)], sh, sd);
    }
    std::fclose(fh);
  }

  int nDdc = 0;
  for (const bool f : aff.ddc) {
    nDdc += f ? 1 : 0;
  }
  std::printf("# traj %s frame %d nop %d rings %zu six %zu ddc %d\n",
              traj.c_str(), frame, cloud.nop, rings.size(), six.size(), nDdc);
  std::printf("io_ms %.3f\nneigh_ms %.3f\nknn_ms %.3f\nknn_pair_ms %.3f\n"
              "rings_ms %.3f\naffil_ms %.3f\nseeded_ms %.3f\n",
              tIo, tNeigh, tKnn, tKnnPair, tRings, tAff, tSeeded);
  return 0;
}
