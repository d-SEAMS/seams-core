/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
** Copyright (c) d-SEAMS developers
**
** Time findDDC alone. Neighbours, rings and findHC run before the timer.
** One warm-up and the timed repeats are the arguments after the trajectory.
**
**   bench_ddc_ladder TRAJ WARMUP REPEATS
*/

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "usage: bench_ddc_ladder TRAJ WARMUP REPEATS\n";
    return 2;
  }
  const std::string traj = argv[1];
  const int warmup = std::atoi(argv[2]);
  const int reps = std::atoi(argv[3]);
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO(traj, 1, cloud, 1);
  if (cloud.nop == 0) {
    std::cerr << "could not read " << traj << "\n";
    return 1;
  }
  auto neighbours = nneigh::neighListO(3.5, cloud, 1);
  auto indexed = nneigh::neighbourListByIndex(cloud, neighbours);
  auto rings = primitive::ringNetwork(indexed, 7);
  std::vector<std::vector<int>> six;
  for (const auto &ring : rings) {
    if (ring.size() == 6) {
      six.push_back(ring);
    }
  }
  const auto index = ring::buildRingSearchIndex(six, cloud);
  std::vector<ring::strucType> ringType(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cages;
  auto listHC = ring::findHC(six, ringType, indexed, cages, index);
  const auto ringTypeAfterHC = ringType;
  const auto cagesAfterHC = cages;
  auto once = [&]() {
    ringType = ringTypeAfterHC;
    cages = cagesAfterHC;
    return ring::findDDC(six, ringType, listHC, cages, index);
  };
  for (int i = 0; i < warmup; i++) {
    once();
  }
  int nDDC = 0;
  int nCages = 0;
  std::cout << "atoms " << cloud.nop << "\n";
  std::cout << "six " << six.size() << "\n";
  for (int i = 0; i < reps; i++) {
    ringType = ringTypeAfterHC;
    cages = cagesAfterHC;
    const auto start = std::chrono::steady_clock::now();
    const auto listDDC = ring::findDDC(six, ringType, listHC, cages, index);
    const auto stop = std::chrono::steady_clock::now();
    nDDC = static_cast<int>(listDDC.size());
    nCages = 0;
    for (const auto &cage : cages) {
      if (cage.type == cage::cageType::DoubleDiaC) {
        nCages++;
      }
    }
    const double ms =
        std::chrono::duration<double, std::milli>(stop - start).count();
    std::cout << ms << "\n";
  }
  std::cout << "ddc_rings " << nDDC << "\n";
  std::cout << "ddc_cages " << nCages << "\n";
  return 0;
}
