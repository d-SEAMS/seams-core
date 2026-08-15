/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Report ice scores for every published bond graph on one frame.
**   report_ice_graphs TRAJ [frame] [atomType]
*/

#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

std::vector<std::vector<int>>
sixOf(const std::vector<std::vector<int>> &rings) {
  std::vector<std::vector<int>> six;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  return six;
}

void tallyRings(const char *name, const ring::CageAffiliation &aff, int nSix) {
  int hc = 0;
  int ddc = 0;
  for (std::size_t i = 0; i < aff.hc.size(); i++) {
    hc += aff.hc[i] ? 1 : 0;
    ddc += aff.ddc[i] ? 1 : 0;
  }
  std::printf("rings  %-10s nSix %d hc %d ddc %d\n", name, nSix, hc, ddc);
}

void tallyAtoms(const char *name, const std::vector<bool> &hc,
                const std::vector<bool> &ddc) {
  int nHC = 0;
  int nDDC = 0;
  int nBoth = 0;
  int nNone = 0;
  const int n = static_cast<int>(hc.size());
  for (int i = 0; i < n; i++) {
    const bool h = hc[static_cast<std::size_t>(i)];
    const bool d = ddc[static_cast<std::size_t>(i)];
    if (h && d) {
      ++nBoth;
    } else if (h) {
      ++nHC;
    } else if (d) {
      ++nDDC;
    } else {
      ++nNone;
    }
  }
  std::printf("atoms  %-10s nop %d ih %d ic %d both %d water %d\n", name, n,
              nHC, nDDC, nBoth, nNone);
}

void atomsFromRings(const std::vector<std::vector<int>> &six,
                    const ring::CageAffiliation &aff, int nop,
                    std::vector<bool> &hc, std::vector<bool> &ddc) {
  hc.assign(static_cast<std::size_t>(nop), false);
  ddc.assign(static_cast<std::size_t>(nop), false);
  for (std::size_t r = 0; r < six.size(); r++) {
    for (const int a : six[r]) {
      if (a >= 0 && a < nop) {
        hc[static_cast<std::size_t>(a)] =
            hc[static_cast<std::size_t>(a)] || aff.hc[r];
        ddc[static_cast<std::size_t>(a)] =
            ddc[static_cast<std::size_t>(a)] || aff.ddc[r];
      }
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  const std::string traj =
      argc > 1 ? argv[1] : "traj/mW_cubic.lammpstrj";
  const int frame = argc > 2 ? std::atoi(argv[2]) : 1;
  const int typeI = argc > 3 ? std::atoi(argv[3]) : 1;

  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO(traj, frame, cloud, typeI);
  if (cloud.nop == 0) {
    std::fprintf(stderr, "empty frame %d of %s\n", frame, traj.c_str());
    return 1;
  }

  const double cutoff = 3.5;
  const double cand = cutoff + 2.0;
  const int k = 4;
  const int nop = cloud.nop;

  auto cutoffRows = nneigh::neighListO(cutoff, cloud, typeI);
  auto knnRows = nneigh::kNearestNeighbourList(cloud, k, cand, typeI, true);
  auto unionRows = nneigh::kNearestNeighbourList(cloud, k, cand, typeI, false);
  auto idxC = nneigh::neighbourListByIndex(cloud, cutoffRows);
  auto idxK = nneigh::neighbourListByIndex(cloud, knnRows);
  auto idxU = nneigh::neighbourListByIndex(cloud, unionRows);
  auto sixC = sixOf(primitive::ringNetwork(idxC, 6));
  auto sixK = sixOf(primitive::ringNetwork(idxK, 6));
  auto sixU = sixOf(primitive::ringNetwork(idxU, 6));
  const auto affC = ring::cageAffiliation(sixC, idxC);
  const auto affK = ring::cageAffiliation(sixK, idxK);
  const auto affU = ring::cageAffiliation(sixU, idxU);
  const auto seeded = ring::seededCageAffiliation(sixK, idxK, sixU, idxU);

  std::printf("# %s frame %d nop %d\n", traj.c_str(), frame, nop);
  tallyRings("cutoff", affC, static_cast<int>(sixC.size()));
  tallyRings("knn", affK, static_cast<int>(sixK.size()));
  tallyRings("knn-union", affU, static_cast<int>(sixU.size()));

  std::vector<bool> hc, ddc;
  atomsFromRings(sixC, affC, nop, hc, ddc);
  tallyAtoms("cutoff", hc, ddc);
  atomsFromRings(sixK, affK, nop, hc, ddc);
  tallyAtoms("knn", hc, ddc);
  atomsFromRings(sixU, affU, nop, hc, ddc);
  tallyAtoms("knn-union", hc, ddc);
  tallyAtoms("seeded", seeded.hc, seeded.ddc);
  return 0;
}
