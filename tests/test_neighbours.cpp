#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <neighbours.hpp>

// Helper: build a 4-atom system in a 10x10x10 box
// Atoms at (0,0,0), (1,0,0), (0,1,0), (5,5,5)
// With cutoff 1.5, atoms 0,1,2 are mutual neighbours; atom 3 is isolated
static molSys::PointCloud<molSys::Point<double>, double> makeFourAtomCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  double coords[4][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {5.0, 5.0, 5.0}};

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

TEST_CASE("neighListO produces full neighbour list by atom ID", "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  double cutoff = 1.5;

  auto nList = nneigh::neighListO(cutoff, &cloud, 1);

  // Should have one row per atom
  REQUIRE(nList.size() == 4);

  // Atom 0 (ID 0): neighbours should include IDs 1 and 2
  // First element is the atom's own ID
  REQUIRE(nList[0][0] == 0);
  REQUIRE(nList[0].size() >= 3); // self + 2 neighbours
  // Check that 1 and 2 appear in the neighbour list of atom 0
  bool found1 =
      std::find(nList[0].begin() + 1, nList[0].end(), 1) != nList[0].end();
  bool found2 =
      std::find(nList[0].begin() + 1, nList[0].end(), 2) != nList[0].end();
  REQUIRE(found1);
  REQUIRE(found2);

  // Atom 3 (ID 3): should have no neighbours (distance > cutoff from all)
  REQUIRE(nList[3][0] == 3);
  REQUIRE(nList[3].size() == 1); // only self
}

TEST_CASE("halfNeighList produces half neighbour list", "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  double cutoff = 1.5;

  auto nList = nneigh::halfNeighList(cutoff, &cloud, 1);

  REQUIRE(nList.size() == 4);

  // In a half list, if 0 lists 1 as neighbour, 1 should NOT list 0
  // Atom 0 should have 1 and 2 as neighbours
  bool zeroHas1 =
      std::find(nList[0].begin() + 1, nList[0].end(), 1) != nList[0].end();
  REQUIRE(zeroHas1);

  // Atom 1 should NOT have 0 as neighbour (half list)
  bool oneHas0 =
      std::find(nList[1].begin() + 1, nList[1].end(), 0) != nList[1].end();
  REQUIRE_FALSE(oneHas0);
}

TEST_CASE("neighbourListByIndex converts ID list to index list",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  double cutoff = 1.5;

  // Get a neighbour list by ID
  auto nListByID = nneigh::neighListO(cutoff, &cloud, 1);
  // Convert to index-based
  auto nListByIdx = nneigh::neighbourListByIndex(&cloud, nListByID);

  REQUIRE(nListByIdx.size() == 4);

  // Since atomID == index in this test cloud, the lists should be identical
  for (int i = 0; i < 4; i++) {
    REQUIRE(nListByIdx[i].size() == nListByID[i].size());
  }
}

TEST_CASE("getNewNeighbourListByIndex builds index-based list from scratch",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  double cutoff = 1.5;

  auto nList = nneigh::getNewNeighbourListByIndex(&cloud, cutoff);

  REQUIRE(nList.size() == 4);
  // Atom 0 (index 0) should have indices 1 and 2 as neighbours
  REQUIRE(nList[0][0] == 0);
  bool found1 =
      std::find(nList[0].begin() + 1, nList[0].end(), 1) != nList[0].end();
  bool found2 =
      std::find(nList[0].begin() + 1, nList[0].end(), 2) != nList[0].end();
  REQUIRE(found1);
  REQUIRE(found2);
}

TEST_CASE("clearNeighbourList empties the list", "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  auto nList = nneigh::neighListO(1.5, &cloud, 1);
  REQUIRE_FALSE(nList.empty());

  nneigh::clearNeighbourList(nList);
  REQUIRE(nList.empty());
}
