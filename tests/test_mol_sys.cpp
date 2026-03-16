#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>

// Helper to build a small PointCloud with known data
static molSys::PointCloud<molSys::Point<double>, double> makeSmallCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  for (int i = 0; i < 4; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = 100 + i;
    pt.molID = 10 + i;
    pt.x = 1.0 * i;
    pt.y = 2.0 * i;
    pt.z = 3.0 * i;
    cloud.pts.push_back(pt);
    cloud.idIndexMap[pt.atomID] = i;
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  return cloud;
}

TEST_CASE("Point default construction", "[mol_sys]") {
  molSys::Point<double> pt;
  // Default iceType should be unclassified
  REQUIRE(pt.iceType == molSys::atom_state_type::unclassified);
  // Default inSlice should be true
  REQUIRE(pt.inSlice == true);
  // c_ij should be empty
  REQUIRE(pt.c_ij.empty());
}

TEST_CASE("PointCloud construction and field access", "[mol_sys]") {
  auto cloud = makeSmallCloud();

  REQUIRE(cloud.nop == 4);
  REQUIRE(cloud.currentFrame == 1);
  REQUIRE(cloud.box.size() == 3);
  REQUIRE(cloud.boxLow.size() == 3);

  // Check that the points were stored correctly
  REQUIRE(cloud.pts[0].atomID == 100);
  REQUIRE(cloud.pts[3].atomID == 103);
  REQUIRE_THAT(cloud.pts[2].x, Catch::Matchers::WithinAbs(2.0, 1e-12));
  REQUIRE_THAT(cloud.pts[2].y, Catch::Matchers::WithinAbs(4.0, 1e-12));

  // Check idIndexMap
  REQUIRE(cloud.idIndexMap.size() == 4);
  REQUIRE(cloud.idIndexMap.at(101) == 1);
}

TEST_CASE("clearPointCloud empties all vectors", "[mol_sys]") {
  auto cloud = makeSmallCloud();
  REQUIRE(cloud.nop == 4);
  REQUIRE_FALSE(cloud.pts.empty());

  cloud = molSys::clearPointCloud(&cloud);

  REQUIRE(cloud.pts.empty());
  REQUIRE(cloud.box.empty());
  REQUIRE(cloud.boxLow.empty());
  REQUIRE(cloud.idIndexMap.empty());
}

TEST_CASE("createIDMolIDmap maps atomID to molID", "[mol_sys]") {
  auto cloud = makeSmallCloud();
  auto idMolMap = molSys::createIDMolIDmap(&cloud);

  REQUIRE(idMolMap.size() == 4);
  // atomID 100 -> molID 10
  REQUIRE(idMolMap[100] == 10);
  // atomID 103 -> molID 13
  REQUIRE(idMolMap[103] == 13);
}

TEST_CASE("createMolIDAtomIDMultiMap maps molID to atomID", "[mol_sys]") {
  auto cloud = makeSmallCloud();
  auto multiMap = molSys::createMolIDAtomIDMultiMap(&cloud);

  REQUIRE(multiMap.size() == 4);
  // Each molID should map to exactly one atomID in this simple case
  auto range = multiMap.equal_range(11);
  int count = 0;
  for (auto it = range.first; it != range.second; ++it) {
    REQUIRE(it->second == 101);
    count++;
  }
  REQUIRE(count == 1);
}

TEST_CASE("searchMolList finds correct index", "[mol_sys]") {
  // Build a mock molList: each row starts with a molID
  std::vector<std::vector<int>> molList = {
      {10, 0, 1}, {11, 2, 3}, {12, 4, 5}};

  REQUIRE(molSys::searchMolList(molList, 10) == 0);
  REQUIRE(molSys::searchMolList(molList, 12) == 2);
  // Not found returns -1
  REQUIRE(molSys::searchMolList(molList, 99) == -1);
}
