#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <order_parameter.hpp>

#include <cmath>
#include <vector>

// Helper to build a PointCloud from a list of (x,y,z) coordinates
static molSys::PointCloud<molSys::Point<double>, double>
makeCloud(const std::vector<std::array<double, 3>> &coords,
          double boxLen = 100.0) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {boxLen, boxLen, boxLen};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  for (int i = 0; i < static_cast<int>(coords.size()); i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i;
    p.molID = 0;
    p.x = coords[i][0];
    p.y = coords[i][1];
    p.z = coords[i][2];
    cloud.pts.push_back(p);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = static_cast<int>(coords.size());
  return cloud;
}

// Regression test for projAreaSingleRing return order.
// The function should return {areaXY, areaXZ, areaYZ} (index 0=XY, 1=XZ, 2=YZ).
// The caller calcCoverageArea reads [0]=XY, [1]=XZ, [2]=YZ.
TEST_CASE("projAreaSingleRing returns areas in XY, XZ, YZ order",
          "[order_parameter]") {
  // A 2x3 rectangle in the XZ plane (y=5 for all points).
  // Vertices: (0,5,0), (2,5,0), (2,5,3), (0,5,3)
  // XY projected area: all y are the same, so area = 0
  // XZ projected area: 2 * 3 = 6
  // YZ projected area: all y are the same, so area = 0
  auto cloud = makeCloud({{0, 5, 0}, {2, 5, 0}, {2, 5, 3}, {0, 5, 3}});
  std::vector<int> ring = {0, 1, 2, 3};

  auto areas = topoparam::projAreaSingleRing(cloud, ring);

  REQUIRE(areas.size() == 3);
  // [0] = XY area = 0 (flat in XY since all y identical)
  REQUIRE_THAT(areas[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
  // [1] = XZ area = 6.0
  REQUIRE_THAT(areas[1], Catch::Matchers::WithinAbs(6.0, 1e-10));
  // [2] = YZ area = 0 (flat in YZ since all y identical)
  REQUIRE_THAT(areas[2], Catch::Matchers::WithinAbs(0.0, 1e-10));
}

// Additional test: rectangle in XY plane
TEST_CASE("projAreaSingleRing XY plane rectangle", "[order_parameter]") {
  // A 4x3 rectangle in the XY plane (z=0 for all).
  // Vertices: (0,0,0), (4,0,0), (4,3,0), (0,3,0)
  // XY area: 4*3 = 12
  // XZ area: 0 (all z identical)
  // YZ area: 0 (all z identical)
  auto cloud = makeCloud({{0, 0, 0}, {4, 0, 0}, {4, 3, 0}, {0, 3, 0}});
  std::vector<int> ring = {0, 1, 2, 3};

  auto areas = topoparam::projAreaSingleRing(cloud, ring);

  REQUIRE_THAT(areas[0], Catch::Matchers::WithinAbs(12.0, 1e-10));
  REQUIRE_THAT(areas[1], Catch::Matchers::WithinAbs(0.0, 1e-10));
  REQUIRE_THAT(areas[2], Catch::Matchers::WithinAbs(0.0, 1e-10));
}
