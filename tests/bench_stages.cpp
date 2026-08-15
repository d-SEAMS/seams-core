/*
** Stage times for one frame: I/O, neighbours, rings, affiliation.
**   bench_stages TRAJ [frame] [atomType] [reps]
*/

#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;

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
  auto knn = nneigh::kNearestNeighbourList(cloud, 4, 5.5, typeI, true);
  auto uni = nneigh::kNearestNeighbourList(cloud, 4, 5.5, typeI, false);
  auto idxS = nneigh::neighbourListByIndex(cloud, knn);
  auto idxU = nneigh::neighbourListByIndex(cloud, uni);
  const double tSeeded = bestMs(
      [&]() {
        auto sR = primitive::ringNetwork(idxS, 6);
        auto uR = primitive::ringNetwork(idxU, 6);
        std::vector<std::vector<int>> s6, u6;
        for (const auto &r : sR) {
          if (r.size() == 6) {
            s6.push_back(r);
          }
        }
        for (const auto &r : uR) {
          if (r.size() == 6) {
            u6.push_back(r);
          }
        }
        (void)ring::seededCageAffiliation(s6, idxS, u6, idxU);
      },
      reps);

  int nDdc = 0;
  for (const bool f : aff.ddc) {
    nDdc += f ? 1 : 0;
  }
  std::printf("# traj %s frame %d nop %d rings %zu six %zu ddc %d\n",
              traj.c_str(), frame, cloud.nop, rings.size(), six.size(), nDdc);
  std::printf("io_ms %.3f\nneigh_ms %.3f\nknn_ms %.3f\nrings_ms %.3f\n"
              "affil_ms %.3f\nseeded_ms %.3f\n",
              tIo, tNeigh, tKnn, tRings, tAff, tSeeded);
  return 0;
}
