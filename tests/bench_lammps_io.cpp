/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Time LAMMPS dump I/O only: nLammpsFrames, a sequential load_frame
** walk, and a last-frame seek. Usage:
**   bench_lammps_io TRAJ [frames] [atomType] [threads]
** frames <= 0 means every ITEM: TIMESTEP. atomType <= 0 reads every atom.
** threads > 1 walks with forEachLammpsFrame (one handle per worker).
*/

#include <seams_input.hpp>

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace {
using Clock = std::chrono::steady_clock;

double ms(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}
} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr,
                 "usage: bench_lammps_io TRAJ [frames] [atomType] [threads]\n");
    return 2;
  }
  const std::string traj = argv[1];
  const int want = argc > 2 ? std::atoi(argv[2]) : 0;
  const int atomType = argc > 3 ? std::atoi(argv[3]) : 0;
  const int threads = argc > 4 ? std::atoi(argv[4]) : 1;

  sinp::dropLammpsDumpIndex(traj);

  const auto t0 = Clock::now();
  const int nframes = sinp::nLammpsFrames(traj);
  const auto t1 = Clock::now();
  if (nframes <= 0) {
    std::fprintf(stderr, "no frames in %s\n", traj.c_str());
    return 1;
  }
  const int frames = (want > 0 && want < nframes) ? want : nframes;

  molSys::PointCloud<molSys::Point<double>, double> cloud;
  std::atomic<int> nop{0};
  std::atomic<int> empty{0};
  const auto t2 = Clock::now();
  if (threads > 1) {
    sinp::forEachLammpsFrame(
        traj, 1, frames, atomType,
        [&](int, molSys::PointCloud<molSys::Point<double>, double> &c) {
          if (c.nop == 0) {
            empty.fetch_add(1);
          } else {
            nop.store(c.nop);
          }
        },
        threads);
  } else {
    for (int frame = 1; frame <= frames; frame++) {
      if (atomType > 0) {
        cloud = sinp::readLammpsTrjO(traj, frame, cloud, atomType);
      } else {
        cloud = sinp::readLammpsTrj(traj, frame, cloud);
      }
      if (cloud.nop == 0) {
        std::fprintf(stderr, "empty frame %d of %s\n", frame, traj.c_str());
        return 1;
      }
      nop.store(cloud.nop);
    }
  }
  const auto t3 = Clock::now();
  if (empty.load() > 0) {
    std::fprintf(stderr, "empty frames in %s\n", traj.c_str());
    return 1;
  }

  if (atomType > 0) {
    cloud = sinp::readLammpsTrjO(traj, frames, cloud, atomType);
  } else {
    cloud = sinp::readLammpsTrj(traj, frames, cloud);
  }
  const auto t4 = Clock::now();

  std::printf("# nframes counted walk_frames nop threads index_ms walk_ms "
              "walk_ms_per_frame last_seek_ms\n");
  std::printf("%d %d %d %d %.3f %.3f %.3f %.3f\n", nframes, frames, nop.load(),
              threads, ms(t0, t1), ms(t2, t3), ms(t2, t3) / frames,
              ms(t3, t4));
  return 0;
}
