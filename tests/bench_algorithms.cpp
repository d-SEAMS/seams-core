/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Performance benchmarks using Catch2's BENCHMARK macro.
** Build target: bench_algorithms (not in the regular test suite).
**
** Run: builddir/tests/bench_algorithms
*/

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <simd_distance.hpp>
#include <bond.hpp>

#include <cmath>
#include <random>
#include <string>
#include <vector>

namespace {

/// Helper: load the example trajectory (oxygen atoms, frame 1).
molSys::PointCloud<molSys::Point<double>, double> loadCloud() {
  std::string traj = "traj/exampleTraj.lammpstrj";
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  return sinp::readLammpsTrjreduced(traj, 1, yCloud, 2, false,
                                     {0, 0, 0}, {0, 0, 0});
}

}  // namespace

TEST_CASE("Benchmark neighListO", "[!benchmark]") {
  auto cloud = loadCloud();

  BENCHMARK("neighListO rcutoff=3.5") {
    return nneigh::neighListO(3.5, cloud, 2);
  };

  BENCHMARK("neighListO rcutoff=4.0") {
    return nneigh::neighListO(4.0, cloud, 2);
  };
}

TEST_CASE("Benchmark ringNetwork", "[!benchmark]") {
  auto cloud = loadCloud();
  auto nList = nneigh::neighListO(3.5, cloud, 2);
  std::string traj = "traj/exampleTraj.lammpstrj";
  auto hbonds = bond::populateHbonds(traj, cloud, nList, 1, 1);
  auto hbondsIdx = nneigh::neighbourListByIndex(cloud, hbonds);

  BENCHMARK("ringNetwork maxDepth=6") {
    return primitive::ringNetwork(hbondsIdx, 6);
  };

  BENCHMARK("ringNetwork maxDepth=4") {
    return primitive::ringNetwork(hbondsIdx, 4);
  };
}

TEST_CASE("Benchmark BatchPeriodicDistSq", "[!benchmark]") {
  constexpr size_t N = 10000;
  std::mt19937_64 rng(42);
  std::uniform_real_distribution<double> dist(-25.0, 25.0);

  std::vector<double> dx(N), dy(N), dz(N), out(N);
  for (size_t i = 0; i < N; i++) {
    dx[i] = dist(rng);
    dy[i] = dist(rng);
    dz[i] = dist(rng);
  }
  double bx = 50.0, by = 50.0, bz = 50.0;

  BENCHMARK("BatchPeriodicDistSq N=10000") {
    seams::BatchPeriodicDistSq(dx.data(), dy.data(), dz.data(),
                               bx, by, bz, out.data(), N);
    return out[0];
  };

  // Scalar baseline for comparison
  BENCHMARK("Scalar periodic dist N=10000") {
    for (size_t i = 0; i < N; i++) {
      double ddx = std::fabs(dx[i]);
      double ddy = std::fabs(dy[i]);
      double ddz = std::fabs(dz[i]);
      ddx -= bx * std::round(ddx / bx);
      ddy -= by * std::round(ddy / by);
      ddz -= bz * std::round(ddz / bz);
      out[i] = ddx * ddx + ddy * ddy + ddz * ddz;
    }
    return out[0];
  };
}
