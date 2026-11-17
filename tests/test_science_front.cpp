//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//-----------------------------------------------------------------------------------

#include <science_front.hpp>
#include <seams_input.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

using Catch::Approx;

TEST_CASE("a block of ice has thickness four and two fronts", "[science]") {
  const std::vector<double> mid = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5};
  const std::vector<double> phi = {0, 0, 0, 1, 1, 1, 1, 0, 0, 0};
  const auto fit = science::fitSolidFraction(phi, mid, 10.0);
  REQUIRE(fit.defined);
  REQUIRE(fit.thickness == Approx(4.0));
  REQUIRE(fit.fronts.size() == 2);
  REQUIRE(fit.fronts[0].position == Approx(3.0));
  REQUIRE(fit.fronts[1].position == Approx(7.0));
  REQUIRE(fit.fronts[0].width == Approx(0.8));
  REQUIRE(fit.fronts[1].width == Approx(0.8));
}

TEST_CASE("a flat profile has no front", "[science]") {
  const std::vector<double> mid = {0.5, 1.5, 2.5};
  const std::vector<double> phi = {0.2, 0.2, 0.2};
  const auto fit = science::fitSolidFraction(phi, mid, 3.0);
  REQUIRE_FALSE(fit.defined);
  REQUIRE(fit.fronts.empty());
}

TEST_CASE("points on a line keep the slope under the bootstrap", "[science]") {
  const std::vector<double> x = {0, 1, 2, 3, 4, 5};
  const std::vector<double> y = {1, 3, 5, 7, 9, 11};
  const auto slope = science::bootstrapSlope(x, y, 400, 11);
  REQUIRE(slope.slope == Approx(2.0));
  REQUIRE(slope.lo == Approx(2.0));
  REQUIRE(slope.hi == Approx(2.0));
  REQUIRE(science::bootstrapSlope({}, {}, 10, 1).n == 0);
}

TEST_CASE("the critical size is the mean of the band around the minimum", "[science]") {
  const std::vector<double> series = {8, 9, 8, 10, 9, 40, 70, 100};
  const auto size = science::criticalSize(series, 500, 23);
  REQUIRE(size.defined);
  REQUIRE(size.upperBound);
  REQUIRE(size.onset == 0);
  REQUIRE(size.window == 5);
  REQUIRE(size.nStar == Approx(8.8));
  REQUIRE(size.sigma > 0.0);
  REQUIRE(size.lo <= size.nStar);
  REQUIRE(size.hi >= size.nStar);
  REQUIRE(size.nStar < 20.0);

  const auto bound = science::criticalSize({10, 20, 40, 80}, 100, 5);
  REQUIRE(bound.upperBound);
  REQUIRE(bound.onset == 0);
  REQUIRE(bound.nStar == Approx(10.0));
  REQUIRE(bound.lo == Approx(10.0));
  REQUIRE(bound.hi == Approx(10.0));
  REQUIRE_FALSE(science::criticalSize({}, 10, 1).defined);

  const auto flat = science::criticalSize({4, 4, 4, 4}, 50, 3);
  REQUIRE(flat.defined);
  REQUIRE_FALSE(flat.upperBound);
  REQUIRE(flat.nStar == Approx(4.0));
}

TEST_CASE("ideal ice Ih is hexagonal under CHILL+ and the shell label", "[science]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO("traj/ice1h.lammpstrj", 1, cloud, 1);
  const auto layers = science::layersOf(cloud, -1, 1.0);
  REQUIRE(layers.n == 128);
  int hex = 0;
  int shellHex = 0;
  double q6 = 0;
  int bins = 0;
  for (const auto &row : layers.rows) {
    hex += row.hexagonal;
    shellHex += row.shellHex;
    if (row.n > 0) {
      q6 += row.q6bar;
      ++bins;
    }
  }
  REQUIRE(hex == 128);
  REQUIRE(shellHex == 128);
  REQUIRE(bins > 0);
  REQUIRE(q6 / static_cast<double>(bins) > 0.4);

  auto labeled = science::nucleusOf(cloud, 1);
  REQUIRE(labeled.chillHex == 128);
  REQUIRE(labeled.shellHex == 128);
  REQUIRE(labeled.polytypeBoth == 128);
  REQUIRE(labeled.polytypeAgree == 128);
  REQUIRE(labeled.extended.size == 128);
}

TEST_CASE("disordering one half lowers the solid fraction there", "[science]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO("traj/ice1h.lammpstrj", 1, cloud, 1);
  REQUIRE(cloud.box.size() >= 3);
  const double mid = 0.5 * cloud.box[2];
  const double lx = cloud.box[0];
  const double ly = cloud.box[1];
  const double lz = cloud.box[2];
  for (auto &p : cloud.pts) {
    if (p.z < mid) {
      continue;
    }
    const double s = static_cast<double>(p.atomID);
    p.x = std::fmod(p.x * 1.7 + 0.37 * s, lx);
    p.y = std::fmod(p.y * 1.3 + 0.19 * s, ly);
    p.z = mid + std::fmod(0.41 * s, lz - mid);
    if (p.x < 0) {
      p.x += lx;
    }
    if (p.y < 0) {
      p.y += ly;
    }
  }
  const auto layers = science::layersOf(cloud, 2, 2.0);
  double lowSolid = 0;
  double highSolid = 0;
  int lowN = 0;
  int highN = 0;
  for (const auto &row : layers.rows) {
    const int solid = row.cubic + row.hexagonal + row.interfacial;
    if (row.mid < mid) {
      lowSolid += solid;
      lowN += row.n;
    } else {
      highSolid += solid;
      highN += row.n;
    }
  }
  REQUIRE(lowN > 0);
  REQUIRE(highN > 0);
  REQUIRE(lowSolid / static_cast<double>(lowN) > highSolid / static_cast<double>(highN));
}
