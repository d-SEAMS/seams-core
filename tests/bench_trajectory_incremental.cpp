/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Frame-by-frame validation and timing of the incremental topological path
** on a real trajectory: RingUpdater against a full ringNetwork by canonical
** ring-set equality, and AffiliationUpdater against batch cageAffiliation by
** flag equality. Emits one whitespace-separated record per frame plus a
** summary, machine-parsable for the reproducibility workflow.
**
** Run from the input/ directory:
**   bench_trajectory_incremental [trajectory] [frames] [atomType] [maxDepth]
*/

#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

double ms(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

std::set<std::vector<int>> canon(const std::vector<std::vector<int>> &rings) {
  std::set<std::vector<int>> out;
  for (auto r : rings) {
    auto it = std::min_element(r.begin(), r.end());
    std::rotate(r.begin(), it, r.end());
    if (r.size() > 2 && r[1] > r.back()) {
      std::reverse(r.begin() + 1, r.end());
    }
    out.insert(r);
  }
  return out;
}

} // namespace

int main(int argc, char **argv) {
  const std::string traj =
      argc > 1 ? argv[1] : "traj/mW_cubic.lammpstrj";
  const int frames = argc > 2 ? std::atoi(argv[2]) : 11;
  const int atomType = argc > 3 ? std::atoi(argv[3]) : 1;
  const int maxDepth = argc > 4 ? std::atoi(argv[4]) : 6;

  primitive::RingUpdater ringUpd(maxDepth);
  ring::AffiliationUpdater affilUpd;
  nneigh::SkinNeighborList skin(3.5, 2.0, atomType);

  double ringFullSum = 0.0, ringIncSum = 0.0;
  double affilBatchSum = 0.0, affilIncSum = 0.0;
  int steadyFrames = 0;
  bool allEqual = true;

  std::printf("# frame nAtoms nRings ringFull_ms ringInc_ms ringRecomp "
              "affilBatch_ms affilInc_ms affilReclass ringsEqual "
              "affilEqual\n");

  for (int frame = 1; frame <= frames; frame++) {
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    cloud = sinp::readLammpsTrjO(traj, frame, cloud, atomType);
    if (cloud.nop == 0) {
      std::fprintf(stderr, "could not read frame %d of %s\n", frame,
                   traj.c_str());
      return 1;
    }
    auto nList = skin.update(cloud);
    auto idx = nneigh::neighbourListByIndex(cloud, nList);

    const auto t0 = Clock::now();
    auto full = primitive::ringNetwork(idx, maxDepth);
    const auto t1 = Clock::now();
    const auto &inc = ringUpd.update(idx);
    const auto t2 = Clock::now();
    const bool ringsEqual = canon(full) == canon(inc);

    // The cage stage runs on the six-membered rings of the incremental set
    std::vector<std::vector<int>> six;
    for (const auto &r : inc) {
      if (r.size() == 6) {
        six.push_back(r);
      }
    }
    const auto t3 = Clock::now();
    const auto batchAffil = ring::cageAffiliation(six, idx);
    const auto t4 = Clock::now();
    const auto &incAffil = affilUpd.update(six, idx);
    const auto t5 = Clock::now();
    const bool affilEqual =
        batchAffil.hc == incAffil.hc && batchAffil.ddc == incAffil.ddc;

    allEqual = allEqual && ringsEqual && affilEqual;
    std::printf("%d %d %zu %.3f %.3f %d %.3f %.3f %d %s %s\n", frame,
                cloud.nop, full.size(), ms(t0, t1), ms(t1, t2),
                ringUpd.lastRecomputedSources(), ms(t3, t4), ms(t4, t5),
                affilUpd.lastReclassified(), ringsEqual ? "yes" : "NO",
                affilEqual ? "yes" : "NO");

    if (frame > 1) {
      ringFullSum += ms(t0, t1);
      ringIncSum += ms(t1, t2);
      affilBatchSum += ms(t3, t4);
      affilIncSum += ms(t4, t5);
      steadyFrames++;
    }
  }

  if (steadyFrames > 0) {
    std::printf("# steady ringFull_ms %.3f ringInc_ms %.3f affilBatch_ms "
                "%.3f affilInc_ms %.3f allEqual %s\n",
                ringFullSum / steadyFrames, ringIncSum / steadyFrames,
                affilBatchSum / steadyFrames, affilIncSum / steadyFrames,
                allEqual ? "yes" : "NO");
  }
  return allEqual ? 0 : 1;
}
