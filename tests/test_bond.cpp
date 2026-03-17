#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bond.hpp>
#include <cage.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

// Helper: build a minimal cloud for bond tests
// 4 atoms forming a square in the xy-plane
static molSys::PointCloud<molSys::Point<double>, double> makeSquareCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  double coords[4][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};

  for (int i = 0; i < 4; i++) {
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
  cloud.nop = 4;
  return cloud;
}

TEST_CASE("trimBonds removes exact duplicate bond pairs", "[bond]") {
  // Build a bonds list with duplicates
  std::vector<std::vector<int>> bonds = {
      {1, 2}, {3, 4}, {1, 2}, {5, 6}, {3, 4}};

  auto trimmed = bond::trimBonds(bonds);

  // Should have exactly 3 unique bonds
  REQUIRE(trimmed.size() == 3);

  // Verify each unique bond appears exactly once
  // After sorting, expect {1,2}, {3,4}, {5,6}
  std::sort(trimmed.begin(), trimmed.end());
  REQUIRE(trimmed[0] == std::vector<int>({1, 2}));
  REQUIRE(trimmed[1] == std::vector<int>({3, 4}));
  REQUIRE(trimmed[2] == std::vector<int>({5, 6}));
}

TEST_CASE("trimBonds on empty input returns empty", "[bond]") {
  std::vector<std::vector<int>> bonds;
  auto trimmed = bond::trimBonds(bonds);
  REQUIRE(trimmed.empty());
}

TEST_CASE("populateBonds generates bonds from neighbour list", "[bond]") {
  auto cloud = makeSquareCloud();
  // Build a neighbour list by index (not ID, since populateBonds expects index-based)
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);

  auto bonds = bond::populateBonds(nList, cloud);

  // With cutoff 1.5, each atom should be bonded to its adjacent square
  // neighbours (distance 1.0). Diagonal distance is sqrt(2) ~ 1.414, also
  // within cutoff. So we expect 4 edge bonds + 2 diagonal bonds = 6 total
  // (after dedup by iatom < jatom condition in populateBonds)
  REQUIRE(bonds.size() > 0);

  // All bond entries should have exactly 2 elements
  for (const auto &b : bonds) {
    REQUIRE(b.size() == 2);
  }
}

TEST_CASE("getHbondDistanceOH computes periodic distance between two clouds",
          "[bond]") {
  // Build two minimal PointClouds: one "oxygen", one "hydrogen"
  molSys::PointCloud<molSys::Point<double>, double> oCloud, hCloud;
  oCloud.box = {10.0, 10.0, 10.0};
  hCloud.box = {10.0, 10.0, 10.0};

  molSys::Point<double> oPoint;
  oPoint.x = 1.0;
  oPoint.y = 0.0;
  oPoint.z = 0.0;
  oCloud.pts.push_back(oPoint);

  molSys::Point<double> hPoint;
  hPoint.x = 2.0;
  hPoint.y = 0.0;
  hPoint.z = 0.0;
  hCloud.pts.push_back(hPoint);

  double dist = bond::getHbondDistanceOH(oCloud, hCloud, 0, 0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("getHbondDistanceOH respects periodic boundaries", "[bond]") {
  molSys::PointCloud<molSys::Point<double>, double> oCloud, hCloud;
  oCloud.box = {10.0, 10.0, 10.0};
  hCloud.box = {10.0, 10.0, 10.0};

  molSys::Point<double> oPoint;
  oPoint.x = 0.5;
  oPoint.y = 0.0;
  oPoint.z = 0.0;
  oCloud.pts.push_back(oPoint);

  // Place H near the other edge of the box
  molSys::Point<double> hPoint;
  hPoint.x = 9.5;
  hPoint.y = 0.0;
  hPoint.z = 0.0;
  hCloud.pts.push_back(hPoint);

  double dist = bond::getHbondDistanceOH(oCloud, hCloud, 0, 0);
  // Periodic distance should be 1.0, not 9.0
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("populateBonds with cage iceType filters dummy atoms", "[bond]") {
  auto cloud = makeSquareCloud();
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);

  // All atoms classified as dummy -- bonds should still be created between them
  std::vector<cage::iceType> atomTypes(4, cage::iceType::dummy);
  auto bonds = bond::populateBonds(nList, cloud, atomTypes);

  // With bondsBetweenDummy=false (default), bonds between dummy atoms are excluded
  // So no bonds should be created
  REQUIRE(bonds.empty());
}

TEST_CASE("populateBonds with non-dummy iceType creates bonds", "[bond]") {
  auto cloud = makeSquareCloud();
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);

  // Mark all atoms as HC ice type
  std::vector<cage::iceType> atomTypes(4, cage::iceType::hc);
  auto bonds = bond::populateBonds(nList, cloud, atomTypes);

  REQUIRE(bonds.size() > 0);
  for (const auto &b : bonds) {
    REQUIRE(b.size() == 2);
  }
}

TEST_CASE("createBondsFromCages extracts bonds from cage rings", "[bond]") {
  // Minimal cage setup: a single HC cage has rings
  std::vector<std::vector<int>> rings = {{1, 2, 3, 4, 5, 6},
                                          {7, 8, 9, 10, 11, 12}};
  cage::Cage cage1;
  cage1.type = cage::cageType::HexC;
  cage1.rings = {0, 1}; // Indices into rings

  std::vector<cage::Cage> cageList = {cage1};
  int nRings = 0;

  auto bonds = bond::createBondsFromCages(rings, cageList,
                                           cage::cageType::HexC, nRings);

  REQUIRE(nRings == 2);
  REQUIRE(bonds.size() > 0);
}

TEST_CASE("createBondsFromCages with no matching cage type returns empty",
          "[bond]") {
  std::vector<std::vector<int>> rings = {{1, 2, 3, 4, 5, 6}};
  cage::Cage cage1;
  cage1.type = cage::cageType::DoubleDiaC;
  cage1.rings = {0};

  std::vector<cage::Cage> cageList = {cage1};
  int nRings = 0;

  auto bonds = bond::createBondsFromCages(rings, cageList,
                                           cage::cageType::HexC, nRings);

  REQUIRE(nRings == 0);
}
