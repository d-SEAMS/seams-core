//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <brine_front.hpp>

#include <cmath>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
using Catch::Matchers::WithinAbs;

Cloud emptyCloud(double lx, double ly, double lz) {
  Cloud cloud;
  cloud.box = {lx, ly, lz};
  cloud.boxLow = {0, 0, 0};
  return cloud;
}

int addPoint(Cloud &cloud, double x, double y, double z, int type) {
  molSys::Point<double> point;
  point.x = x;
  point.y = y;
  point.z = z;
  point.type = type;
  point.atomID = static_cast<int>(cloud.pts.size()) + 1;
  point.molID = point.atomID;
  cloud.pts.push_back(point);
  cloud.nop = static_cast<int>(cloud.pts.size());
  return cloud.nop - 1;
}

} // namespace

TEST_CASE("front crossings and the rejection coefficient on a labelled slab", "[brine]") {
  auto cloud = emptyCloud(10, 10, 10);
  std::vector<bool> ice;
  std::vector<int> waters;
  std::vector<int> ions;
  for (int i = 0; i < 10; ++i) {
    const int index = addPoint(cloud, 5, 5, i + 0.5, 1);
    waters.push_back(index);
    ice.push_back(i >= 2 && i <= 5);
  }
  ions.push_back(addPoint(cloud, 5, 5, 8.5, 3));
  ions.push_back(addPoint(cloud, 5, 5, 4.5, 3));
  ice.push_back(false);
  ice.push_back(false);

  const auto profile = site::concentrationProfile(cloud, ice, waters, ions, 10, 2, 2);
  REQUIRE(profile.bins.size() == 10);
  REQUIRE(profile.bins[2].iceFraction == 1);
  REQUIRE(profile.bins[8].iceFraction == 0);
  REQUIRE(profile.bins[4].nIon == 1);
  REQUIRE(profile.bins[8].nIon == 1);

  const auto fronts = site::locateFronts(profile, 0.5);
  REQUIRE(fronts.found);
  REQUIRE_FALSE(fronts.wrapped);
  REQUIRE_THAT(fronts.zEnter, WithinAbs(2.0, 1e-9));
  REQUIRE_THAT(fronts.zExit, WithinAbs(6.0, 1e-9));
  REQUIRE_THAT(fronts.thickness, WithinAbs(4.0, 1e-9));

  const auto rejection = site::rejectionCoefficient(profile);
  REQUIRE(rejection.defined);
  REQUIRE(rejection.nIceIon == 1);
  REQUIRE(rejection.nLiquidIon == 1);
  REQUIRE(rejection.nIceWater == 4);
  REQUIRE(rejection.nLiquidWater == 6);
  REQUIRE_THAT(rejection.k, WithinAbs(1.5, 1e-9));
}

TEST_CASE("empty bins between crystal planes stay inside one slab", "[brine]") {
  auto cloud = emptyCloud(10, 10, 10);
  std::vector<bool> ice;
  std::vector<int> waters;
  for (int bin : {0, 2, 4, 7, 8}) {
    waters.push_back(addPoint(cloud, 1, 1, bin + 0.5, 1));
    ice.push_back(bin <= 4);
  }
  const auto profile = site::concentrationProfile(cloud, ice, waters, {}, 10, 2, 2);
  const auto fronts = site::locateFronts(profile, 0.5);
  REQUIRE(fronts.found);
  REQUIRE_FALSE(fronts.wrapped);
  REQUIRE_THAT(fronts.zEnter, WithinAbs(9.5, 1e-9));
  REQUIRE_THAT(fronts.zExit, WithinAbs(6.0, 1e-9));
  REQUIRE_THAT(fronts.thickness, WithinAbs(6.5, 1e-9));
}

TEST_CASE("an ice slab that crosses the periodic boundary is one slab", "[brine]") {
  auto cloud = emptyCloud(10, 10, 10);
  std::vector<bool> ice;
  std::vector<int> waters;
  for (int i = 0; i < 10; ++i) {
    waters.push_back(addPoint(cloud, 1, 1, i + 0.5, 1));
    ice.push_back(i <= 1 || i >= 8);
  }
  const auto profile = site::concentrationProfile(cloud, ice, waters, {}, 10, 2, 2);
  const auto fronts = site::locateFronts(profile, 0.5);
  REQUIRE(fronts.found);
  REQUIRE(fronts.wrapped);
  REQUIRE_THAT(fronts.zEnter, WithinAbs(8.0, 1e-9));
  REQUIRE_THAT(fronts.zExit, WithinAbs(2.0, 1e-9));
  REQUIRE_THAT(fronts.thickness, WithinAbs(4.0, 1e-9));
}

TEST_CASE("ordinary least squares reports the slope and its standard error", "[brine]") {
  const auto exact = site::linearFit({0, 1, 2, 3}, {0, 2, 4, 6});
  REQUIRE(exact.n == 4);
  REQUIRE_THAT(exact.slope, WithinAbs(2.0, 1e-12));
  REQUIRE_THAT(exact.intercept, WithinAbs(0.0, 1e-12));
  REQUIRE_THAT(exact.slopeStdErr, WithinAbs(0.0, 1e-12));

  const auto noisy = site::linearFit({0, 1, 2}, {0, 0, 1});
  REQUIRE_THAT(noisy.slope, WithinAbs(0.5, 1e-12));
  REQUIRE_THAT(noisy.slopeStdErr, WithinAbs(1.0 / std::sqrt(12.0), 1e-12));
}

TEST_CASE("a liquid cluster inside the ice is a brine pocket and holds its ion", "[brine]") {
  auto cloud = emptyCloud(20, 10, 20);
  std::vector<bool> ice;
  std::vector<int> waters;
  for (int i = 0; i < 5; ++i) {
    waters.push_back(addPoint(cloud, i, 0, 0, 1));
    ice.push_back(false);
  }
  for (int i = 0; i < 3; ++i) {
    waters.push_back(addPoint(cloud, i, 0, 8, 1));
    ice.push_back(false);
  }
  const int ion = addPoint(cloud, 0.2, 0, 8, 3);
  ice.push_back(false);
  const auto census = site::brinePockets(cloud, ice, waters, {ion}, 1.5);
  REQUIRE(census.reservoirWaters == 5);
  REQUIRE(census.nPockets == 1);
  REQUIRE(census.sizes == std::vector<int>{3});
  REQUIRE(census.pockets[0].nIons == 1);
}

TEST_CASE("the front layer keeps its hexagon and counts the missing neighbours", "[brine]") {
  auto cloud = emptyCloud(20, 20, 20);
  std::vector<int> waters;
  constexpr double kPi = 3.14159265358979323846;
  for (int i = 0; i < 6; ++i) {
    const double angle = i * kPi / 3.0;
    waters.push_back(addPoint(cloud, 10 + std::cos(angle), 10 + std::sin(angle), 5, 1));
  }
  site::FrontPair fronts;
  fronts.found = true;
  fronts.zEnter = 5;
  fronts.zExit = 15;
  const auto layer = site::frontLayer(cloud, waters, fronts, 1.0, 1.2, 6, 2);
  REQUIRE(layer.waters.size() == 6);
  REQUIRE(layer.nDefect == 6);
  REQUIRE(layer.nFour == 0);
  REQUIRE(layer.ringCensus[6] == 1);
}

TEST_CASE("moving into ice is a trapping event", "[brine]") {
  site::IonEnvironment before;
  site::IonEnvironment after;
  before.state = {site::IonState::liquid, site::IonState::front, site::IonState::ice};
  after.state = {site::IonState::ice, site::IonState::liquid, site::IonState::ice};
  const auto events = site::trappingTransitions(before, after);
  REQUIRE(events.size() == 2);
  REQUIRE(events[0].ion == 0);
  REQUIRE(events[0].from == site::IonState::liquid);
  REQUIRE(events[0].to == site::IonState::ice);
  REQUIRE(events[1].to == site::IonState::liquid);
}
