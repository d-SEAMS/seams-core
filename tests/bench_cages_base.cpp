/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Cage-stage bench against the pre-index findHC/findDDC signatures, so the
** baseline tree of a two-tree comparison can be timed with the same driver
** shape as bench_cages. Compiled directly against the baseline library by
** scripts/elja_paper_benches.sh rather than through the meson test setup,
** since the baseline tree predates the bench targets.
**
** Run from the input/ directory.
*/

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

template <typename Fn> double bestMillis(Fn &&fn, int reps) {
  double best = std::numeric_limits<double>::max();
  for (int rep = 0; rep < reps; rep++) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    best = std::min(
        best, std::chrono::duration<double, std::milli>(end - start).count());
  }
  return best;
}

} // namespace

int main(int argc, char **argv) {
  const std::string traj = argc > 1 ? argv[1] : "traj/mW_cubic.lammpstrj";
  const int frame = argc > 2 ? std::atoi(argv[2]) : 1;
  const int reps = argc > 3 ? std::atoi(argv[3]) : 3;
  const int atomType = argc > 4 ? std::atoi(argv[4]) : 1;

  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO(traj, frame, yCloud, atomType);
  if (yCloud.nop == 0) {
    std::cerr << "could not read " << traj << "\n";
    return 1;
  }

  auto nList = nneigh::neighListO(3.5, yCloud, atomType);
  auto hbondIdx = nneigh::neighbourListByIndex(yCloud, nList);

  std::vector<std::vector<int>> rings = primitive::ringNetwork(hbondIdx, 7);
  const double tRings =
      bestMillis([&]() { rings = primitive::ringNetwork(hbondIdx, 7); }, reps);

  std::vector<std::vector<int>> sixRings;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      sixRings.push_back(r);
    }
  }

  std::cout << "atoms      " << yCloud.nop << "\n"
            << "rings      " << rings.size() << "\n"
            << "six-rings  " << sixRings.size() << "\n\n";

  std::vector<int> listHC, listDDC;
  std::vector<cage::Cage> cageList;
  std::vector<ring::strucType> ringType;

  const double tHC = bestMillis(
      [&]() {
        ringType.assign(sixRings.size(), ring::strucType::unclassified);
        cageList.clear();
        listHC = ring::findHC(sixRings, ringType, hbondIdx, cageList);
      },
      reps);

  const double tDDC = bestMillis(
      [&]() {
        std::vector<ring::strucType> rt(sixRings.size(),
                                        ring::strucType::unclassified);
        std::vector<cage::Cage> cl;
        auto hc = ring::findHC(sixRings, rt, hbondIdx, cl);
        listDDC = ring::findDDC(sixRings, rt, hc, cl);
      },
      reps);

  ringType.assign(sixRings.size(), ring::strucType::unclassified);
  cageList.clear();
  listHC = ring::findHC(sixRings, ringType, hbondIdx, cageList);
  listDDC = ring::findDDC(sixRings, ringType, listHC, cageList);

  std::cout << "ringNetwork/ms " << tRings << "\n"
            << "findHC/ms      " << tHC << "\n"
            << "findDDC/ms     " << tDDC << "\n"
            << "nHCrings       " << listHC.size() << "\n"
            << "nDDCrings      " << listDDC.size() << "\n"
            << "nCages         " << cageList.size() << "\n";
  return 0;
}
