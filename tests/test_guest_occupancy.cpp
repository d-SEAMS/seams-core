#include <catch2/catch_test_macros.hpp>

#include <mol_sys.hpp>
#include <site.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

// Eight "water" vertices on a cube of edge 4 whose centre sits at (2, 2, 2)
// and a second cube shifted along x by the box length minus its edge, so
// that its vertices straddle the periodic boundary; guests are appended.
Cloud twoCubes(const std::vector<std::array<double, 3>> &guests) {
  Cloud cloud;
  const double L = 20.0;
  int id = 0;
  auto add = [&](double x, double y, double z, int type) {
    molSys::Point<double> p;
    p.x = std::fmod(x + L, L);
    p.y = std::fmod(y + L, L);
    p.z = std::fmod(z + L, L);
    p.type = type;
    p.atomID = id + 1;
    p.molID = id + 1;
    cloud.pts.push_back(p);
    cloud.idIndexMap[id + 1] = id;
    ++id;
  };
  for (int c = 0; c < 2; c++) {
    const double ox = c == 0 ? 0.0 : L - 2.0;  // second cube spans x in [18, 22) -> wraps
    for (int i = 0; i < 8; i++) {
      add(ox + 4.0 * (i & 1), 4.0 * ((i >> 1) & 1), 4.0 * ((i >> 2) & 1), 1);
    }
  }
  for (const auto &g : guests) {
    add(g[0], g[1], g[2], 2);
  }
  cloud.nop = id;
  cloud.box = {L, L, L};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

std::vector<int> range(int lo, int hi) {
  std::vector<int> v;
  for (int i = lo; i < hi; i++) {
    v.push_back(i);
  }
  return v;
}

} // namespace

TEST_CASE("the periodic centroid of a cage across the boundary is inside it", "[guests]") {
  const auto cloud = twoCubes({});
  const auto c0 = site::periodicCentroid(cloud, range(0, 8));
  REQUIRE(std::fabs(c0[0] - 2.0) < 1e-9);
  REQUIRE(std::fabs(c0[1] - 2.0) < 1e-9);
  REQUIRE(std::fabs(c0[2] - 2.0) < 1e-9);
  const auto c1 = site::periodicCentroid(cloud, range(8, 16));
  // the wrapped cube spans x in [18, 22): centre at 20 == 0 modulo the box
  const double wrapped = std::fmod(c1[0] + 20.0, 20.0);
  REQUIRE((std::fabs(wrapped) < 1e-9 || std::fabs(wrapped - 20.0) < 1e-9));
  REQUIRE(std::fabs(c1[1] - 2.0) < 1e-9);
}

TEST_CASE("guests are assigned to the cage whose centre they are nearest", "[guests]") {
  // one guest at the first centre, one near the wrapped centre seen from
  // the far side of the box, one in open space
  const auto cloud = twoCubes({{2.0, 2.0, 2.0}, {19.5, 2.0, 2.2}, {10.0, 10.0, 10.0}});
  const std::vector<std::vector<int>> cages = {range(0, 8), range(8, 16)};
  const auto occ = site::guestOccupancy(cloud, cages, {16, 17, 18}, 3.0);
  REQUIRE(occ.cageOfGuest == std::vector<int>{0, 1, -1});
  REQUIRE(occ.guestsPerCage == std::vector<int>{1, 1});
  REQUIRE(occ.occupied == 2);
  REQUIRE(occ.multiply == 0);
  REQUIRE(occ.free == 1);
  REQUIRE(std::fabs(occ.centreDistance[0]) < 1e-9);
  REQUIRE(std::fabs(occ.centreDistance[1] - std::sqrt(0.25 + 0.04)) < 1e-9);
  REQUIRE(occ.centreDistance[2] == -1.0);
}

TEST_CASE("ray-parity occupancy keeps a centre guest and rejects one outside",
          "[guests]") {
  const auto cloud = twoCubes({{2.0, 2.0, 2.0}, {10.0, 10.0, 10.0}, {2.0, 2.0, -0.5}});
  // Cube 0: vertices 0-7. Faces as 6 rings of 4.
  const std::vector<std::vector<int>> rings = {
      {0, 1, 3, 2}, {4, 5, 7, 6}, {0, 1, 5, 4},
      {2, 3, 7, 6}, {0, 2, 6, 4}, {1, 3, 7, 5},
      {8, 9, 11, 10}, {12, 13, 15, 14}, {8, 9, 13, 12},
      {10, 11, 15, 14}, {8, 10, 14, 12}, {9, 11, 15, 13}};
  const std::vector<std::vector<int>> faces = {
      {0, 1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11}};
  const auto occ = site::guestOccupancyInside(cloud, rings, faces, {16, 17, 18});
  REQUIRE(occ.cageOfGuest == std::vector<int>{0, -1, -1});
  REQUIRE(occ.guestsPerCage == std::vector<int>{1, 0});
  REQUIRE(occ.occupied == 1);
  REQUIRE(occ.free == 2);
}

TEST_CASE("two guests in one cage count as multiple occupancy", "[guests]") {
  const auto cloud = twoCubes({{1.5, 2.0, 2.0}, {2.5, 2.0, 2.0}});
  const std::vector<std::vector<int>> cages = {range(0, 8), range(8, 16)};
  const auto occ = site::guestOccupancy(cloud, cages, {16, 17}, 3.0);
  REQUIRE(occ.guestsPerCage == std::vector<int>{2, 0});
  REQUIRE(occ.occupied == 1);
  REQUIRE(occ.multiply == 1);
  REQUIRE(occ.free == 0);
  REQUIRE(occ.occupancyHistogram.size() >= 3);
  REQUIRE(occ.occupancyHistogram[0] == 1);
  REQUIRE(occ.occupancyHistogram[2] == 1);
  // a radius too small leaves every guest free
  const auto none = site::guestOccupancy(cloud, cages, {16, 17}, 0.1);
  REQUIRE(none.free == 2);
  REQUIRE(none.occupied == 0);
}
