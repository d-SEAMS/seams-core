/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Cage-stage bench against the pre-index findHC/findDDC signatures, so the
** baseline tree of a two-tree comparison can be timed with the same driver
** shape as bench_cages. Compiled directly against the baseline library by
** the paper bench script rather than through the meson test setup,
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

// setup() runs before every rep, outside the timer.
template <typename Setup, typename Fn>
double bestMillis(Setup &&setup, Fn &&fn, int reps) {
  double best = std::numeric_limits<double>::max();
  for (int rep = 0; rep < reps; rep++) {
    setup();
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    best = std::min(
        best, std::chrono::duration<double, std::milli>(end - start).count());
  }
  return best;
}

template <typename Fn> double bestMillis(Fn &&fn, int reps) {
  return bestMillis([]() {}, fn, reps);
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

  // findDDC reads the ring types and cages findHC left; that state is
  // computed once here and restored outside the timer before every rep.
  ringType.assign(sixRings.size(), ring::strucType::unclassified);
  cageList.clear();
  listHC = ring::findHC(sixRings, ringType, hbondIdx, cageList);
  const std::vector<ring::strucType> ringTypeAfterHC = ringType;
  const std::vector<cage::Cage> cageListAfterHC = cageList;

  const double tDDC = bestMillis(
      [&]() {
        ringType = ringTypeAfterHC;
        cageList = cageListAfterHC;
      },
      [&]() { listDDC = ring::findDDC(sixRings, ringType, listHC, cageList); },
      reps);

  std::cout << "ringNetwork/ms " << tRings << "\n"
            << "findHC/ms      " << tHC << "\n"
            << "findDDC/ms     " << tDDC << "\n"
            << "pipeline/ms    " << (tRings + tHC + tDDC) << "\n"
            << "nHCrings       " << listHC.size() << "\n"
            << "nDDCrings      " << listDDC.size() << "\n"
            << "nCages         " << cageList.size() << "\n";
  return 0;
}
