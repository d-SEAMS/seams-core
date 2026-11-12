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

#include <nucleation_pathway.hpp>

#include <cmath>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

Cloud diamond(int nrep) {
  const double oo = 2.75;
  const double a = oo * 4.0 / std::sqrt(3.0);
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  Cloud cloud;
  cloud.box = {nrep * a, nrep * a, nrep * a};
  cloud.boxLow = {0, 0, 0};
  for (int i = 0; i < nrep; ++i) {
    for (int j = 0; j < nrep; ++j) {
      for (int k = 0; k < nrep; ++k) {
        for (const auto &shift : fcc) {
          for (const auto &site : basis) {
            molSys::Point<double> p;
            p.x = (i + shift[0] + site[0]) * a;
            p.y = (j + shift[1] + site[1]) * a;
            p.z = (k + shift[2] + site[2]) * a;
            p.type = 1;
            p.atomID = cloud.nop + 1;
            p.molID = p.atomID;
            cloud.idIndexMap[p.atomID] = cloud.nop;
            cloud.pts.push_back(p);
            cloud.nop += 1;
          }
        }
      }
    }
  }
  return cloud;
}

Cloud dilute(int n) {
  Cloud cloud;
  cloud.box = {80, 80, 80};
  cloud.boxLow = {0, 0, 0};
  for (int i = 0; i < n; ++i) {
    molSys::Point<double> p;
    p.x = 3.0 * i;
    p.y = 1.0 * ((i * 5) % 17);
    p.z = 2.0 * ((i * 3) % 19);
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    cloud.idIndexMap[p.atomID] = cloud.nop;
    cloud.pts.push_back(p);
    cloud.nop += 1;
  }
  return cloud;
}

} // namespace

TEST_CASE("stack fraction is cubic over cubic plus hexagonal", "[nucleation]") {
  REQUIRE_THAT(nucleation::stackFraction(3, 1),
               Catch::Matchers::WithinAbs(0.75, 1e-12));
  REQUIRE(std::isnan(nucleation::stackFraction(0, 0)));
}

TEST_CASE("final ascent stops at the last point before the strict climb",
          "[nucleation]") {
  const nucleation::Ascent climb = nucleation::finalAscent({10, 10, 10, 50});
  REQUIRE(climb.climbed);
  REQUIRE(climb.index == 2);
  REQUIRE(climb.maxIndex == 3);
  REQUIRE_THAT(climb.baseline, Catch::Matchers::WithinAbs(10, 1e-12));
  REQUIRE_THAT(climb.maximum, Catch::Matchers::WithinAbs(50, 1e-12));

  const nucleation::Ascent dipped = nucleation::finalAscent({5, 8, 3, 4, 9});
  REQUIRE(dipped.index == 2);
  REQUIRE(dipped.maxIndex == 4);
  REQUIRE_THAT(dipped.baseline, Catch::Matchers::WithinAbs(3, 1e-12));

  const nucleation::Ascent flat = nucleation::finalAscent({4, 4, 4});
  REQUIRE_FALSE(flat.climbed);
  REQUIRE(flat.index == 0);
}

TEST_CASE("ring fraction leads ice by the constructed lag", "[nucleation]") {
  // The ring differences are the ice differences shifted three samples earlier.
  std::vector<double> rings = {0.20, 0.20, 0.25, 0.45, 0.55, 0.57, 0.57,
                               0.57, 0.57, 0.57, 0.57, 0.57};
  std::vector<double> ice = {0.10, 0.10, 0.10, 0.10, 0.10, 0.15, 0.35,
                             0.45, 0.47, 0.47, 0.47, 0.47};
  const nucleation::Lead lead = nucleation::ringLead(rings, ice);
  REQUIRE(lead.ok);
  REQUIRE(lead.ascentFrames == 3);
  REQUIRE(lead.halfFrames == 3);
  REQUIRE(lead.corrLag == 3);
  REQUIRE(lead.corr > 0.9);
}

TEST_CASE("half rise is the first sample past half the gain", "[nucleation]") {
  const nucleation::HalfRise rise = nucleation::halfRise({0, 0, 1, 2, 3, 4});
  REQUIRE(rise.ok);
  REQUIRE(rise.index == 3);
  REQUIRE_THAT(rise.level, Catch::Matchers::WithinAbs(2.0, 1e-12));
  const nucleation::HalfRise flat = nucleation::halfRise({3, 3, 2});
  REQUIRE_FALSE(flat.ok);
}

TEST_CASE("definition shift is the difference of onset size and frame",
          "[nucleation]") {
  const nucleation::Shift shift = nucleation::definitionShift(40, 12, 55, 14);
  REQUIRE(shift.dSize == 15);
  REQUIRE(shift.dFrame == 2);
}

TEST_CASE("a periodic diamond is cubic ice under both CHILL+ rules",
          "[nucleation]") {
  auto crystal = diamond(2);
  const auto ice =
      nucleation::analyzeFrame(crystal, 3.5, 6, 1, nucleation::kQ6BarIceCut, 1);
  REQUIRE(ice.nop == 64);
  REQUIRE(ice.chillCubic == 64);
  REQUIRE(ice.ovitoCubic == 64);
  REQUIRE(ice.chillHex == 0);
  REQUIRE(ice.ringCount[6] > 0);
  REQUIRE(ice.symCubic == 64);
  REQUIRE(ice.nChillCore == 64);
  REQUIRE(ice.q6barMean > nucleation::kQ6BarIceCut);
  REQUIRE(ice.q6barIce == 64);
  REQUIRE_THAT(nucleation::stackFraction(ice.cChillCore, ice.hChillCore),
               Catch::Matchers::WithinAbs(1.0, 1e-12));

  auto gas = dilute(64);
  const auto bare =
      nucleation::analyzeFrame(gas, 3.5, 6, 1, nucleation::kQ6BarIceCut, 1);
  REQUIRE(bare.q6barMean < ice.q6barMean);
  REQUIRE(bare.chillCubic < ice.chillCubic);
  REQUIRE(bare.symCubic < ice.symCubic);
}
