/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Time a chunked frame-parallel walk: read + neighbours + ringNetwork
** + cageAffiliation per frame. Incremental updaters are not used;
** each chunk is independent. Usage:
**   bench_lammps_walk TRAJ [frames] [atomType] [threads]
*/

#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;

double ms(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

int classifyFrame(molSys::PointCloud<molSys::Point<double>, double> &cloud,
                  int atomType) {
  auto nList = nneigh::neighListO(3.5, cloud, atomType);
  auto idx = nneigh::neighbourListByIndex(cloud, nList);
  auto rings = primitive::ringNetwork(idx, 6);
  std::vector<std::vector<int>> six;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  const auto aff = ring::cageAffiliation(six, idx);
  int ddc = 0;
  for (const bool flag : aff.ddc) {
    ddc += flag ? 1 : 0;
  }
  return ddc;
}
} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr,
                 "usage: bench_lammps_walk TRAJ [frames] [atomType] [threads]\n");
    return 2;
  }
  const std::string traj = argv[1];
  const int want = argc > 2 ? std::atoi(argv[2]) : 0;
  const int atomType = argc > 3 ? std::atoi(argv[3]) : 1;
  const int threads = argc > 4 ? std::atoi(argv[4]) : 1;

  const int nframes = sinp::nLammpsFrames(traj);
  if (nframes <= 0) {
    std::fprintf(stderr, "no frames in %s\n", traj.c_str());
    return 1;
  }
  const int frames = (want > 0 && want < nframes) ? want : nframes;
  std::atomic<int> lastDdc{0};
  std::atomic<int> empty{0};

  const auto t0 = Clock::now();
  sinp::forEachLammpsFrame(
      traj, 1, frames, atomType,
      [&](int, molSys::PointCloud<molSys::Point<double>, double> &cloud) {
        if (cloud.nop == 0) {
          empty.fetch_add(1);
          return;
        }
        lastDdc.store(classifyFrame(cloud, atomType));
      },
      threads);
  const auto t1 = Clock::now();
  if (empty.load() > 0) {
    std::fprintf(stderr, "empty frames in %s\n", traj.c_str());
    return 1;
  }

  std::printf("# frames nop threads wall_ms ms_per_frame last_ddc_rings\n");
  std::printf("%d %d %d %.3f %.3f %d\n", frames, atomType == 0 ? 0 : 1,
              threads, ms(t0, t1), ms(t0, t1) / frames, lastDdc.load());
  return 0;
}
