//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

/*
** Per-frame cost of RingUpdater and AffiliationUpdater against the batch
** passes as a function of how much of a frame changes. Frame A of a
** trajectory is the reference; in each test frame the molecules within a
** radius of the box centre take their frame-B positions and the rest keep
** their frame-A ones. The updaters alternate between A and the test frame,
** so every timed update sees the same change; each time is the median of the
** repetitions. Every incremental result is checked against the batch one.
**
** Usage: bench_churn_sweep TRAJ [frameA] [frameB] [cutoff] [atomType] [reps]
**                          [ringFullPassFraction] [affilBatchFraction]
*/

#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using Rings = std::vector<std::vector<int>>;
using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

double ms(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}

Rings sixOf(const Rings &all) {
  Rings out;
  for (const auto &r : all) {
    if (r.size() == 6) {
      out.push_back(r);
    }
  }
  return out;
}

Rings canon(const Rings &rings) {
  Rings out = rings;
  for (auto &r : out) {
    std::sort(r.begin(), r.end());
  }
  std::sort(out.begin(), out.end());
  return out;
}

Rings graphOf(const Cloud &cloud, double cutoff, int atomType) {
  return nneigh::neighbourListByIndex(
      cloud, nneigh::neighListO(cutoff, cloud, atomType));
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: bench_churn_sweep TRAJ [frameA] [frameB] "
                         "[cutoff] [atomType] [reps] [ringFraction] "
                         "[affilFraction]\n");
    return 2;
  }
  const std::string traj = argv[1];
  const int frameA = argc > 2 ? std::atoi(argv[2]) : 1;
  const int frameB = argc > 3 ? std::atoi(argv[3]) : 2;
  const double cutoff = argc > 4 ? std::atof(argv[4]) : 3.5;
  const int atomType = argc > 5 ? std::atoi(argv[5]) : 1;
  const int reps = argc > 6 ? std::atoi(argv[6]) : 7;

  Cloud a;
  Cloud b;
  a = sinp::readLammpsTrjO(traj, frameA, a, atomType);
  b = sinp::readLammpsTrjO(traj, frameB, b, atomType);
  std::unordered_map<int, int> inB;
  for (int i = 0; i < b.nop; i++) {
    inB[b.pts[static_cast<std::size_t>(i)].atomID] = i;
  }
  const Rings idxA = graphOf(a, cutoff, atomType);
  const Rings sixA = sixOf(primitive::ringNetwork(idxA, 6));
  std::vector<std::vector<int>> rowsA(idxA.size());
  for (std::size_t v = 0; v < idxA.size(); v++) {
    rowsA[v].assign(idxA[v].begin() + 1, idxA[v].end());
    std::sort(rowsA[v].begin(), rowsA[v].end());
  }

  std::printf("# %s frames %d->%d cutoff %.2f nop %d reps %d\n", traj.c_str(),
              frameA, frameB, cutoff, a.nop, reps);
  std::printf("# radius moved changedAtoms bondChurn nSix ringRecomp "
              "affilReclass ringBatch_ms ringInc_ms affilBatch_ms "
              "affilInc_ms ringsEqual affilEqual\n");

  const double radii[] = {0, 3, 5, 7, 9, 11, 13, 15, 18, 21, 25, 1e9};
  bool allEqual = true;
  for (const double radius : radii) {
    Cloud c = a;
    int moved = 0;
    for (std::size_t i = 0; i < static_cast<std::size_t>(c.nop); i++) {
      double d2 = 0.0;
      const double p[3] = {c.pts[i].x, c.pts[i].y, c.pts[i].z};
      for (std::size_t k = 0; k < 3; k++) {
        double d = p[k] - (c.boxLow[k] + 0.5 * c.box[k]);
        d -= c.box[k] * std::round(d / c.box[k]);
        d2 += d * d;
      }
      const auto it = inB.find(c.pts[i].atomID);
      if (d2 < radius * radius && it != inB.end()) {
        const auto j = static_cast<std::size_t>(it->second);
        c.pts[i].x = b.pts[j].x;
        c.pts[i].y = b.pts[j].y;
        c.pts[i].z = b.pts[j].z;
        ++moved;
      }
    }
    const Rings idxC = graphOf(c, cutoff, atomType);
    int changedAtoms = 0;
    int churn = 0;
    std::vector<int> row;
    for (std::size_t v = 0; v < idxC.size(); v++) {
      row.assign(idxC[v].begin() + 1, idxC[v].end());
      std::sort(row.begin(), row.end());
      if (row != rowsA[v]) {
        ++changedAtoms;
        std::vector<int> diff;
        std::set_symmetric_difference(row.begin(), row.end(),
                                      rowsA[v].begin(), rowsA[v].end(),
                                      std::back_inserter(diff));
        churn += static_cast<int>(diff.size());
      }
    }
    churn /= 2;

    std::vector<double> tRingBatch, tRingInc, tAffBatch, tAffInc;
    Rings batchRings;
    for (int r = 0; r < reps; r++) {
      // Construction only. The previous result is released after the
      // timer, matching bench_trajectory_incremental.
      const auto t0 = Clock::now();
      auto fresh = primitive::ringNetwork(idxC, 6);
      const auto t1 = Clock::now();
      tRingBatch.push_back(ms(t0, t1));
      batchRings.swap(fresh);
    }
    const Rings sixC = sixOf(batchRings);
    ring::CageAffiliation batchAff;
    for (int r = 0; r < reps; r++) {
      const auto t0 = Clock::now();
      auto fresh = ring::cageAffiliation(sixC, idxC);
      const auto t1 = Clock::now();
      tAffBatch.push_back(ms(t0, t1));
      batchAff = std::move(fresh);
    }

    primitive::RingUpdater ringUpd(6);
    ring::AffiliationUpdater affUpd;
    if (argc > 7) {
      ringUpd.setFullPassFraction(std::atof(argv[7]));
    }
    if (argc > 8) {
      affUpd.setBatchFraction(std::atof(argv[8]));
    }
    (void)ringUpd.update(idxA);
    (void)affUpd.update(sixA, idxA);
    bool ringsEqual = true;
    bool affEqual = true;
    int ringRecomp = 0;
    int affReclass = 0;
    for (int r = 0; r < reps; r++) {
      auto t0 = Clock::now();
      const Rings &inc = ringUpd.update(idxC);
      tRingInc.push_back(ms(t0, Clock::now()));
      ringRecomp = ringUpd.lastRecomputedSources();
      if (r == 0) {
        ringsEqual = canon(inc) == canon(batchRings);
      }
      t0 = Clock::now();
      const auto &aff = affUpd.update(sixC, idxC);
      tAffInc.push_back(ms(t0, Clock::now()));
      affReclass = affUpd.lastReclassified();
      affEqual = affEqual && aff.hc == batchAff.hc && aff.ddc == batchAff.ddc;
      (void)ringUpd.update(idxA);
      (void)affUpd.update(sixA, idxA);
    }
    allEqual = allEqual && ringsEqual && affEqual;
    std::printf("%g %d %d %d %zu %d %d %.3f %.3f %.3f %.3f %s %s\n",
                radius >= 1e9 ? -1.0 : radius, moved, changedAtoms, churn,
                sixC.size(), ringRecomp, affReclass, median(tRingBatch),
                median(tRingInc), median(tAffBatch), median(tAffInc),
                ringsEqual ? "yes" : "NO", affEqual ? "yes" : "NO");
  }
  return allEqual ? 0 : 1;
}
