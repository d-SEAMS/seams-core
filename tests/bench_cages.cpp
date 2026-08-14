/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Cost of the cage-detection stage: given a primitive ring network, how long
** does it take to classify hexagonal rings into double-diamond and hexagonal
** cages.  Uses the real mW trajectory rather than a synthetic cell, because
** the cage search is sensitive to the actual ring topology in a way a jittered
** lattice does not reproduce.
**
** Build target: bench_cages.  Run from the input/ directory.
*/

#include <bop.hpp>
#include <bulkTUM.hpp>
#include <cage.hpp>
#include <cage_canon.hpp>
#include <franzblau.hpp>
#include <ira_sofi.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <sphericart_ylm.hpp>
#include <topo_bulk.hpp>

#include <chrono>
#include <iomanip>
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
  const std::string traj =
      argc > 1 ? argv[1] : "traj/mW_cubic.lammpstrj";
  const int frame = argc > 2 ? std::atoi(argv[2]) : 1;
  const int reps = argc > 3 ? std::atoi(argv[3]) : 3;

  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO(traj, frame, yCloud, 1);
  if (yCloud.nop == 0) {
    std::cerr << "could not read " << traj << "\n";
    return 1;
  }

  // mW is monatomic, so the ring network is built from the oxygen neighbour
  // list directly rather than from a hydrogen-bond network
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto hbondIdx = nneigh::neighbourListByIndex(yCloud, nList);

  std::vector<std::vector<int>> rings = primitive::ringNetwork(hbondIdx, 7);
  const double tRings =
      bestMillis([&]() { rings = primitive::ringNetwork(hbondIdx, 7); }, reps);

  // Only the six-membered rings take part in DDC and HC detection
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

  // Steinhardt parameters, which the OpenMP build spreads over threads
  auto qRef = chill::steinhardtQl(yCloud, nList, 6);
  const double tSteinhardt = bestMillis(
      [&]() {
        volatile auto q = chill::steinhardtQl(yCloud, nList, 6);
        (void)q;
      },
      reps);

  double tNauty = -1.0;
  if (cage::nautyAvailable()) {
    tNauty = bestMillis(
        [&]() {
          volatile auto cert = cage::canonicalCertificate(sixRings);
          (void)cert;
        },
        reps);
  }

  double tMatchHC = -1.0;
  int nMatched = 0;
  const char *hcTemplate = argc > 4 ? argv[4] : "../templates/hc.xyz";
  Eigen::MatrixXd refHC = tum3::buildRefHC(hcTemplate);
  if (refHC.rows() > 0) {
    tMatchHC = bestMillis(
        [&]() {
          nMatched = 0;
          for (const auto &c : cageList) {
            if (c.type != cage::cageType::HexC) {
              continue;
            }
            std::vector<double> quat;
            double rmsd = 0.0;
            (void)tum3::shapeMatchHC(yCloud, refHC, c, sixRings, hbondIdx, quat,
                                     rmsd);
            nMatched++;
          }
        },
        reps);
  }

  std::cout << "backends   vesin=on"
            << " sphericart=" << (seams::sphericart_ylm::available() ? "on" : "off")
            << " nauty=" << (cage::nautyAvailable() ? "on" : "off")
            << " ira=" << (ira::available() ? "on" : "off") << "\n\n";

  std::cout << std::left << std::setw(28) << "steinhardtQl l=6/ms"
            << std::fixed << std::setprecision(3) << tSteinhardt << "\n";
  std::cout << std::left << std::setw(28) << "ringNetwork/ms" << std::fixed
            << std::setprecision(3) << tRings << "\n"
            << std::setw(28) << "findHC/ms" << tHC << "\n"
            << std::setw(28) << "findHC+findDDC/ms" << tDDC << "\n"
            << std::setw(28) << "findDDC alone/ms" << (tDDC - tHC) << "\n";
  if (tNauty >= 0.0) {
    std::cout << std::setw(28) << "nauty cert six-rings/ms" << tNauty << "\n";
  }
  if (tMatchHC >= 0.0) {
    std::cout << std::setw(28) << "shapeMatchHC/ms" << tMatchHC << "\n"
              << "matched HCs " << nMatched << "\n";
  }
  std::cout << "\nHC rings   " << listHC.size() << "\n"
            << "DDC rings  " << listDDC.size() << "\n";

  return 0;
}
