#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <cmath>
#include <complex>
#include <numeric>

TEST_CASE("radialCoord converts Cartesian to spherical angles", "[bop]") {
  // Point along the z-axis: (0, 0, 1)
  // Polar angle (theta) should be 0, azimuth (phi) should be 0
  std::array<double, 3> zAxis = {0.0, 0.0, 1.0};
  auto angles = sph::radialCoord(zAxis);

  // angles[1] is the polar angle (theta from z-axis)
  REQUIRE_THAT(angles[1], Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("radialCoord for point along x-axis", "[bop]") {
  std::array<double, 3> xAxis = {1.0, 0.0, 0.0};
  auto angles = sph::radialCoord(xAxis);

  // Polar angle should be pi/2 (perpendicular to z)
  REQUIRE_THAT(angles[1], Catch::Matchers::WithinAbs(M_PI / 2.0, 1e-10));
}

TEST_CASE("spheriHarmo returns correct vector length for l=3", "[bop]") {
  std::array<double, 2> angles = {0.5, 1.0}; // phi, theta
  auto result = sph::spheriHarmo(3, angles);

  // 2*3+1 = 7 components
  REQUIRE(result.size() == 7);
}

TEST_CASE("spheriHarmo returns correct vector length for l=6", "[bop]") {
  std::array<double, 2> angles = {0.5, 1.0};
  auto result = sph::spheriHarmo(6, angles);

  // 2*6+1 = 13 components
  REQUIRE(result.size() == 13);
}

TEST_CASE("spheriHarmo Y_l0 at theta=0 is real and matches known value",
          "[bop]") {
  // At theta=0 (along z-axis), Y_lm for m!=0 should vanish
  // Y_30(0,phi) = sqrt(7/(4*pi)) for theta=0
  std::array<double, 2> angles = {0.0, 0.0}; // phi=0, theta=0
  auto result = sph::spheriHarmo(3, angles);

  // m ranges from -3 to 3; index 3 corresponds to m=0
  // For m!=0, values should be zero at theta=0
  for (int k = 0; k < 7; k++) {
    if (k == 3)
      continue; // skip m=0
    REQUIRE_THAT(std::abs(result[k]), Catch::Matchers::WithinAbs(0.0, 1e-10));
  }

  // m=0 component should be nonzero
  double y30_expected = 0.5 * std::sqrt(7.0 / M_PI);
  REQUIRE_THAT(result[3].real(),
               Catch::Matchers::WithinAbs(y30_expected, 1e-10));
  REQUIRE_THAT(result[3].imag(), Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("lookupTableQ3Vec matches spheriHarmo for l=3", "[bop]") {
  std::array<double, 2> angles = {0.7, 1.2}; // arbitrary phi, theta
  auto boostResult = sph::spheriHarmo(3, angles);
  auto lookupResult = sph::lookupTableQ3Vec(angles);

  REQUIRE(boostResult.size() == lookupResult.size());

  for (size_t i = 0; i < boostResult.size(); i++) {
    REQUIRE_THAT(boostResult[i].real(),
                 Catch::Matchers::WithinAbs(lookupResult[i].real(), 1e-8));
    REQUIRE_THAT(boostResult[i].imag(),
                 Catch::Matchers::WithinAbs(lookupResult[i].imag(), 1e-8));
  }
}

TEST_CASE("lookupTableQ6Vec matches spheriHarmo for l=6", "[bop]") {
  std::array<double, 2> angles = {0.3, 0.8};
  auto boostResult = sph::spheriHarmo(6, angles);
  auto lookupResult = sph::lookupTableQ6Vec(angles);

  REQUIRE(boostResult.size() == lookupResult.size());

  for (size_t i = 0; i < boostResult.size(); i++) {
    REQUIRE_THAT(boostResult[i].real(),
                 Catch::Matchers::WithinAbs(lookupResult[i].real(), 1e-8));
    REQUIRE_THAT(boostResult[i].imag(),
                 Catch::Matchers::WithinAbs(lookupResult[i].imag(), 1e-8));
  }
}

TEST_CASE("getCorrel populates c_ij for a simple tetrahedral system", "[bop]") {
  // Build a 5-atom system: central atom + 4 tetrahedral neighbours
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  // Central atom at origin
  double coords[5][3] = {
      {10.0, 10.0, 10.0},                       // center
      {11.0, 11.0, 11.0},                        // vertex 1
      {11.0, 9.0, 9.0},                          // vertex 2
      {9.0, 11.0, 9.0},                          // vertex 3
      {9.0, 9.0, 11.0}                           // vertex 4
  };

  for (int i = 0; i < 5; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.molID = i;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 5;

  // Build neighbour list
  auto nList = nneigh::neighListO(3.0, &cloud, 1);

  // Get correlations
  cloud = chill::getCorrel(&cloud, nList, false);

  // The central atom (index 0) should now have c_ij entries
  REQUIRE_FALSE(cloud.pts[0].c_ij.empty());

  // Each c_ij entry should have a real c_value in [-1, 1]
  for (const auto &cij : cloud.pts[0].c_ij) {
    REQUIRE(cij.c_value >= -1.0);
    REQUIRE(cij.c_value <= 1.0);
  }
}
