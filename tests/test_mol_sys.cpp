#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>

#include <algorithm>

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

  cloud = molSys::clearPointCloud(cloud);

  REQUIRE(cloud.pts.empty());
  REQUIRE(cloud.box.empty());
  REQUIRE(cloud.boxLow.empty());
  REQUIRE(cloud.idIndexMap.empty());
  REQUIRE(cloud.nop == 0);
}

TEST_CASE("createIDMolIDmap maps atomID to molID", "[mol_sys]") {
  auto cloud = makeSmallCloud();
  auto idMolMap = molSys::createIDMolIDmap(cloud);

  REQUIRE(idMolMap.size() == 4);
  // atomID 100 -> molID 10
  REQUIRE(idMolMap[100] == 10);
  // atomID 103 -> molID 13
  REQUIRE(idMolMap[103] == 13);
}

TEST_CASE("createMolIDAtomIDMultiMap maps molID to atomID", "[mol_sys]") {
  auto cloud = makeSmallCloud();
  auto multiMap = molSys::createMolIDAtomIDMultiMap(cloud);

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

// The list keeps the first two hydrogens in hydrogen-cloud order, including
// when that order is not the atom-id order and a molecule has 0, 1 or 3.
std::vector<std::vector<int>>
scanHydrogens(molSys::PointCloud<molSys::Point<double>, double> &hydrogen,
              molSys::PointCloud<molSys::Point<double>, double> &oxygen) {
  std::vector<std::vector<int>> rows;
  for (int iatom = 0; iatom < oxygen.nop; iatom++) {
    const int molID = oxygen.pts[static_cast<std::size_t>(iatom)].molID;
    std::vector<int> row;
    row.push_back(molID);
    int nHatoms = 0;
    for (int jatom = 0; jatom < hydrogen.nop; jatom++) {
      if (hydrogen.pts[static_cast<std::size_t>(jatom)].molID == molID) {
        row.push_back(jatom);
        nHatoms++;
        if (nHatoms == 2) {
          break;
        }
      }
    }
    rows.push_back(std::move(row));
  }
  return rows;
}

molSys::Point<double> placed(int atomID, int molID, double x) {
  molSys::Point<double> pt;
  pt.type = 1;
  pt.atomID = atomID;
  pt.molID = molID;
  pt.x = x;
  pt.y = 0.0;
  pt.z = 0.0;
  return pt;
}

TEST_CASE("hAtomMolList matches a scan of a shuffled hydrogen cloud",
          "[mol_sys]") {
  molSys::PointCloud<molSys::Point<double>, double> oxygen;
  molSys::PointCloud<molSys::Point<double>, double> hydrogen;
  oxygen.box = hydrogen.box = {10.0, 10.0, 10.0};
  oxygen.boxLow = hydrogen.boxLow = {0.0, 0.0, 0.0};
  // mol 4 has no hydrogen, mol 5 has one, mol 6 has three, mol 8 has two.
  // Two oxygens share mol 6. Hydrogen order is not molecule order: the
  // third hydrogen of mol 6 sits between the first and the second.
  oxygen.pts.push_back(placed(1, 6, 0.0));
  oxygen.pts.push_back(placed(2, 4, 1.0));
  oxygen.pts.push_back(placed(3, 5, 2.0));
  oxygen.pts.push_back(placed(4, 8, 3.0));
  oxygen.pts.push_back(placed(5, 6, 4.0));
  oxygen.nop = static_cast<int>(oxygen.pts.size());
  hydrogen.pts.push_back(placed(11, 8, 0.2));
  hydrogen.pts.push_back(placed(12, 6, 0.4));
  hydrogen.pts.push_back(placed(13, 6, 0.6));
  hydrogen.pts.push_back(placed(14, 5, 0.8));
  hydrogen.pts.push_back(placed(15, 6, 1.2));
  hydrogen.pts.push_back(placed(16, 8, 1.4));
  hydrogen.nop = static_cast<int>(hydrogen.pts.size());

  const auto scanned = scanHydrogens(hydrogen, oxygen);
  REQUIRE(molSys::hAtomMolList(hydrogen, oxygen) == scanned);
  REQUIRE(scanned[0] == std::vector<int>({6, 1, 2}));
  REQUIRE(scanned[1] == std::vector<int>({4}));
  REQUIRE(scanned[2] == std::vector<int>({5, 3}));
  REQUIRE(scanned[3] == std::vector<int>({8, 0, 5}));
  REQUIRE(scanned[4] == std::vector<int>({6, 1, 2}));

  std::reverse(hydrogen.pts.begin(), hydrogen.pts.end());
  REQUIRE(molSys::hAtomMolList(hydrogen, oxygen) ==
          scanHydrogens(hydrogen, oxygen));
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
