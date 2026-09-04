#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <order_parameter.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
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

TEST_CASE("rodgerF4 is cos 3 phi on a known H-O-O-H dihedral",
          "[order_parameter]") {
  // Two waters. Outer hydrogens give a 90 degree H-O-O-H dihedral
  // (phi = pi/2, cos 3phi = 0) or a planar 0 degree one (cos 3phi = 1).
  auto twoWaters = [](double h2y, double h2z) {
    // O1 (0,0,0) mol 1; H_outer (0,1,0); H_inner (1,0,0)
    // O2 (3,0,0) mol 2; H_outer (3,h2y,h2z); H_inner (2,0,0)
    auto cloud = makeCloud({{0, 0, 0}, {0, 1, 0}, {1, 0, 0},
                            {3, 0, 0}, {3, h2y, h2z}, {2, 0, 0}});
    cloud.pts[0].type = 1;
    cloud.pts[0].molID = 1;
    cloud.pts[0].atomID = 1;
    cloud.pts[1].type = 2;
    cloud.pts[1].molID = 1;
    cloud.pts[1].atomID = 2;
    cloud.pts[2].type = 2;
    cloud.pts[2].molID = 1;
    cloud.pts[2].atomID = 3;
    cloud.pts[3].type = 1;
    cloud.pts[3].molID = 2;
    cloud.pts[3].atomID = 4;
    cloud.pts[4].type = 2;
    cloud.pts[4].molID = 2;
    cloud.pts[4].atomID = 5;
    cloud.pts[5].type = 2;
    cloud.pts[5].molID = 2;
    cloud.pts[5].atomID = 6;
    cloud.idIndexMap.clear();
    for (int i = 0; i < 6; i++) {
      cloud.idIndexMap[cloud.pts[static_cast<std::size_t>(i)].atomID] = i;
    }
    return cloud;
  };
  // nList by atom ID, leading self
  const std::vector<std::vector<int>> nList = {
      {1, 4}, {2}, {3}, {4, 1}, {5}, {6}};
  auto planar = twoWaters(1.0, 0.0);
  const auto f0 = topoparam::rodgerF4(planar, nList, 1, 2);
  REQUIRE(std::isfinite(f0[0]));
  REQUIRE_THAT(f0[0], Catch::Matchers::WithinAbs(1.0, 1e-9));
  REQUIRE_THAT(f0[3], Catch::Matchers::WithinAbs(1.0, 1e-9));
  REQUIRE_THAT(topoparam::meanFinite(f0), Catch::Matchers::WithinAbs(1.0, 1e-9));

  auto right = twoWaters(0.0, 1.0);
  const auto f90 = topoparam::rodgerF4(right, nList, 1, 2);
  REQUIRE_THAT(f90[0], Catch::Matchers::WithinAbs(0.0, 1e-9));
  REQUIRE_THAT(f90[3], Catch::Matchers::WithinAbs(0.0, 1e-9));
}

TEST_CASE("rodgerF4 is NaN on mW with no hydrogens", "[order_parameter]") {
  auto cloud = makeCloud({{0, 0, 0}, {3, 0, 0}});
  cloud.pts[0].atomID = 1;
  cloud.pts[1].atomID = 2;
  cloud.idIndexMap = {{1, 0}, {2, 1}};
  const std::vector<std::vector<int>> nList = {{1, 2}, {2, 1}};
  const auto f4 = topoparam::rodgerF4(cloud, nList, 1, 2);
  REQUIRE_FALSE(std::isfinite(f4[0]));
  REQUIRE_FALSE(std::isfinite(f4[1]));
  REQUIRE_FALSE(std::isfinite(topoparam::meanFinite(f4)));
}

TEST_CASE("CHILL+ layerCubicity reproduces the literature I_sd string",
          "[order_parameter]") {
  // Four layers along z at 0, 3.7, 7.4, 11.1 in a 14.8 box: C H C H
  std::vector<std::array<double, 3>> coords;
  for (int k = 0; k < 4; k++) {
    coords.push_back({1.0, 1.0, k * 3.7});
    coords.push_back({2.0, 2.0, k * 3.7});
  }
  auto cloud = makeCloud(coords, 14.8);
  for (int i = 0; i < 8; i++) {
    cloud.pts[static_cast<std::size_t>(i)].iceType =
        (i / 2) % 2 == 0 ? molSys::atom_state_type::cubic
                         : molSys::atom_state_type::hexagonal;
  }
  const auto st = topoparam::layerCubicity(cloud, 2, 3.7);
  REQUIRE(st.sequence == "CHCH");
  REQUIRE_THAT(st.phiC, Catch::Matchers::WithinAbs(0.5, 1e-12));
  REQUIRE(st.phiC > 0.0);
  REQUIRE(st.phiC < 1.0);
}

TEST_CASE("TUM stacking uses ring planes and stays empty without basal rings",
          "[order_parameter]") {
  // Two disjoint six-rings stacked along z, marked basal. Their centroids
  // fall in two H layers. A five-ring (clathrate face) does not vote.
  auto cloud = makeCloud({{0, 0, 0}, {2, 0, 0}, {3, 1.5, 0},
                          {2, 3, 0}, {0, 3, 0}, {-1, 1.5, 0},
                          {0, 0, 7.4}, {2, 0, 7.4}, {3, 1.5, 7.4},
                          {2, 3, 7.4}, {0, 3, 7.4}, {-1, 1.5, 7.4},
                          {1, 1, 3.7}, {2, 1, 3.7}, {2.5, 2, 3.7},
                          {1.5, 2.8, 3.7}, {0.5, 2, 3.7}},
                         14.8);
  const std::vector<std::vector<int>> rings = {
      {0, 1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11}, {12, 13, 14, 15, 16}};
  std::vector<bool> basal = {true, true, false};
  std::vector<bool> equatorial = {false, false, false};
  const auto tum = topoparam::tumLayerStack(cloud, rings, basal, equatorial, 2, 3.7);
  REQUIRE(tum.sequence.size() == 4);
  REQUIRE(tum.sequence[0] == 'H');
  REQUIRE(tum.sequence[2] == 'H');
  REQUIRE(tum.phiC == 0.0);
  // the five-ring is not a plane: no C letter
  REQUIRE(tum.sequence.find('C') == std::string::npos);

  // A clathrate 5-ring is not a stacking plane even if a caller flags it.
  std::vector<bool> flaggedFive = {true, true, true};
  const auto five = topoparam::tumLayerStack(cloud, rings, flaggedFive,
                                            equatorial, 2, 3.7);
  REQUIRE(five.sequence == tum.sequence);
  REQUIRE(five.phiC == 0.0);

  std::vector<bool> eq = {false, true, false};
  const auto mixed = topoparam::tumLayerStack(cloud, rings, basal, eq, 2, 3.7);
  REQUIRE(mixed.phiC > 0.0);
  REQUIRE(mixed.phiC < 1.0);
  REQUIRE(mixed.sequence.find('C') != std::string::npos);
  REQUIRE(mixed.sequence.find('H') != std::string::npos);
}

TEST_CASE("ions stay off the type-filtered oxygen neighbour graph",
          "[order_parameter]") {
  auto cloud = makeCloud({{0, 0, 0}, {2.8, 0, 0}, {1.4, 2.4, 0},
                          {0, 0, 2.8}, {10, 10, 10}},
                         20.0);
  cloud.pts[4].type = 3;
  cloud.pts[4].atomID = 99;
  cloud.idIndexMap.clear();
  for (int i = 0; i < cloud.nop; i++) {
    cloud.idIndexMap[cloud.pts[static_cast<std::size_t>(i)].atomID] = i;
  }
  const auto nList = nneigh::kNearestNeighbourList(cloud, 4, 5.5, 1, true);
  REQUIRE(nList.size() == static_cast<std::size_t>(cloud.nop));
  REQUIRE(nList[4].size() == 1);
  REQUIRE(nList[4][0] == 99);
  for (int i = 0; i < 4; i++) {
    REQUIRE(std::find(nList[static_cast<std::size_t>(i)].begin(),
                      nList[static_cast<std::size_t>(i)].end(),
                      99) == nList[static_cast<std::size_t>(i)].end());
  }
}

TEST_CASE("normHeightPercent uses recovered lz not tilt",
          "[order_parameter]") {
  auto cloud = makeCloud({{0.0, 0.0, 0.0}});
  cloud.box = {61.0, 12.0, 50.0, 60.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  const double h = topoparam::normHeightPercent(cloud, 10, 2.5);
  REQUIRE_THAT(h, Catch::Matchers::WithinAbs(50.0, 1e-10));
}
