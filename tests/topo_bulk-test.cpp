// Internal
#include <cage.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>
#include <topo_bulk.hpp>

// Standard
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>

#include <catch2/catch_test_macros.hpp>
#include <rang.hpp>

namespace fs = std::filesystem;

// Helper: build the 12-atom HC geometry (reusable across tests)
static void buildHCCloud(molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  yCloud.box = {53.9690018, 54.5289994, 51.257};
  yCloud.boxLow = {0.0, 0.0, 0.0};
  yCloud.currentFrame = 1;
  molSys::Point<double> iPoint;
  iPoint.type = 1;
  double coords[12][3] = {
      {8.995, 10.3859997, 15.0939999}, {6.7459998, 14.2810001, 15.0939999},
      {4.4970002, 10.3859997, 15.0939999}, {8.995, 12.9829998, 14.1949997},
      {6.7459998, 9.0880003, 14.1949997}, {4.4970002, 12.9829998, 14.1949997},
      {8.995, 12.9829998, 11.4329996}, {6.7459998, 9.0880003, 11.4329996},
      {4.4970002, 12.9829998, 11.4329996}, {8.995, 10.3859997, 10.5340004},
      {6.7459998, 14.2810001, 10.5340004}, {4.4970002, 10.3859997, 10.5340004}};
  for (int i = 0; i < 12; i++) {
    iPoint.atomID = i;
    iPoint.x = coords[i][0]; iPoint.y = coords[i][1]; iPoint.z = coords[i][2];
    yCloud.pts.push_back(iPoint);
    yCloud.idIndexMap[i] = i;
  }
  yCloud.nop = 12;
}

SCENARIO("Test the HC algorithm for a single hexagonal cage.", "[topo]") {
  GIVEN("A pointCloud") {
    // Hard-coded example of a single tetragonal prism
    // Here, we manually input a tetragonal prism's coordinates
    molSys::PointCloud<molSys::Point<double>, double> yCloud;  // pointCloud
    molSys::Point<double> iPoint;                              // A single point
    std::vector<double> box = {53.9690018, 54.5289994, 51.257};  // Box lengths
    yCloud.box = box;  // Update the box lengths. z dim is axial
    std::vector<std::vector<int>> nList;  // Neighbour list
    std::vector<std::vector<int>> rings;  // Rings
    std::vector<ring::strucType>
        ringType;  // This vector will have a value for each ring inside
    std::vector<int> listHC;  // Contains atom indices of atoms making up HCs
    // Make a list of all the DDCs and HCs
    std::vector<cage::Cage> cageList;
    // --------------------
    // Building the hexagonal cage
    //
    iPoint.type = 1;  // Same for all the points here
    // Element {0}
    iPoint.atomID = 0;      // iatom
    iPoint.x = 8.995;       // x
    iPoint.y = 10.3859997;  // y
    iPoint.z = 15.0939999;  // z
    yCloud.pts.push_back(iPoint);
    // Element {1}
    iPoint.atomID = 1;      // iatom
    iPoint.x = 6.7459998;   // x
    iPoint.y = 14.2810001;  // y
    iPoint.z = 15.0939999;  // z
    yCloud.pts.push_back(iPoint);
    // Element {2}
    iPoint.atomID = 2;      // iatom
    iPoint.x = 4.4970002;   // x
    iPoint.y = 10.3859997;  // y
    iPoint.z = 15.0939999;  // z
    yCloud.pts.push_back(iPoint);
    // Element {3}
    iPoint.atomID = 3;      // iatom
    iPoint.x = 8.995;       // x
    iPoint.y = 12.9829998;  // y
    iPoint.z = 14.1949997;  // z
    yCloud.pts.push_back(iPoint);
    // Element {4}
    iPoint.atomID = 4;      // iatom
    iPoint.x = 6.7459998;   // x
    iPoint.y = 9.0880003;   // y
    iPoint.z = 14.1949997;  // z
    yCloud.pts.push_back(iPoint);
    // Element {5}
    iPoint.atomID = 5;      // iatom
    iPoint.x = 4.4970002;   // x
    iPoint.y = 12.9829998;  // y
    iPoint.z = 14.1949997;  // z
    yCloud.pts.push_back(iPoint);
    // Element {6}
    iPoint.atomID = 6;      // iatom
    iPoint.x = 8.995;       // x
    iPoint.y = 12.9829998;  // y
    iPoint.z = 11.4329996;  // z
    yCloud.pts.push_back(iPoint);
    // Element {7}
    iPoint.atomID = 7;      // iatom
    iPoint.x = 6.7459998;   // x
    iPoint.y = 9.0880003;   // y
    iPoint.z = 11.4329996;  // z
    yCloud.pts.push_back(iPoint);
    // Element {8}
    iPoint.atomID = 8;      // iatom
    iPoint.x = 4.4970002;   // x
    iPoint.y = 12.9829998;  // y
    iPoint.z = 11.4329996;  // z
    yCloud.pts.push_back(iPoint);
    // Element {9}
    iPoint.atomID = 9;      // iatom
    iPoint.x = 8.995;       // x
    iPoint.y = 10.3859997;  // y
    iPoint.z = 10.5340004;  // z
    yCloud.pts.push_back(iPoint);
    // Element {10}
    iPoint.atomID = 10;     // iatom
    iPoint.x = 6.7459998;   // x
    iPoint.y = 14.2810001;  // y
    iPoint.z = 10.5340004;  // z
    yCloud.pts.push_back(iPoint);
    // Element {11}
    iPoint.atomID = 11;     // iatom
    iPoint.x = 4.4970002;   // x
    iPoint.y = 10.3859997;  // y
    iPoint.z = 10.5340004;  // z
    yCloud.pts.push_back(iPoint);
    // Update nop
    yCloud.nop = static_cast<int>(yCloud.pts.size());
    // Update the unordered map
    for (int iatom = 0; iatom < yCloud.nop; iatom++) {
      yCloud.idIndexMap[iatom] = iatom;
    }  // end of filling the map
    // --------------------
    WHEN("Given a pointCloud, and a neighbour list") {
      // Calculate a neighbour list
      nList = nneigh::neighListO(3.5, yCloud, 1);
      // Neighbour list by index
      nList = nneigh::neighbourListByIndex(yCloud, nList);
      // Find the vector of vector of rings
      rings = primitive::ringNetwork(nList, 7);
      THEN("There should be exactly one hexagonal cage.") {
        ringType.resize(
            rings.size());  // Has a value for each ring. init to zero.
        // Find the number of hexagonal cages
        listHC = ring::findHC(rings, ringType, nList, cageList);

        // Assert the number of cages
        REQUIRE(cageList.size() ==
                1);  // Evaluate condition for a single tetragonal prism
      }
    }  // End of getting the neighbour list
  }    // End of given
}  // End of scenario

// DDC determination
SCENARIO("Test the DDC algorithm for a single double-diamond cage.", "[topo]") {
  GIVEN("A pointCloud") {
    // Hard-coded example of a single tetragonal prism
    // Here, we manually input a tetragonal prism's coordinates
    molSys::PointCloud<molSys::Point<double>, double> yCloud;  // pointCloud
    molSys::Point<double> iPoint;                              // A single point
    std::vector<double> box = {53.9690018, 54.5289994, 51.257};  // Box lengths
    yCloud.box = box;  // Update the box lengths. z dim is axial
    std::vector<std::vector<int>> nList;  // Neighbour list
    std::vector<std::vector<int>> rings;  // Rings
    std::vector<ring::strucType>
        ringType;  // This vector will have a value for each ring inside
    std::vector<int> listDDC;  // Contains atom indices of atoms making up HCs
    // Make a list of all the DDCs and HCs
    std::vector<cage::Cage> cageList;
    std::vector<int> listHC;  // Empty
    // --------------------
    // Building the double-diamond cage
    //
    iPoint.type = 1;  // Same for all the points here
    // Element {0}
    iPoint.atomID = 0;      // iatom
    iPoint.x = 25.4799996;  // x
    iPoint.y = 25.4799996;  // y
    iPoint.z = 19.1100006;  // z
    yCloud.pts.push_back(iPoint);
    // Element {1}
    iPoint.atomID = 1;      // iatom
    iPoint.x = 25.4799996;  // x
    iPoint.y = 25.4799996;  // y
    iPoint.z = 25.4799996;  // z
    yCloud.pts.push_back(iPoint);
    // Element {2}
    iPoint.atomID = 2;      // iatom
    iPoint.x = 25.4799996;  // x
    iPoint.y = 22.2950001;  // y
    iPoint.z = 22.2950001;  // z
    yCloud.pts.push_back(iPoint);
    // Element {3}
    iPoint.atomID = 3;      // iatom
    iPoint.x = 25.4799996;  // x
    iPoint.y = 28.6650009;  // y
    iPoint.z = 22.2950001;  // z
    yCloud.pts.push_back(iPoint);
    // Element {4}
    iPoint.atomID = 4;      // iatom
    iPoint.x = 22.2950001;  // x
    iPoint.y = 22.2950001;  // y
    iPoint.z = 25.4799996;  // z
    yCloud.pts.push_back(iPoint);
    // Element {5}
    iPoint.atomID = 5;      // iatom
    iPoint.x = 22.2950001;  // x
    iPoint.y = 25.4799996;  // y
    iPoint.z = 22.2950001;  // z
    yCloud.pts.push_back(iPoint);
    // Element {6}
    iPoint.atomID = 6;      // iatom
    iPoint.x = 28.6650009;  // x
    iPoint.y = 25.4799996;  // y
    iPoint.z = 22.2950001;  // z
    yCloud.pts.push_back(iPoint);
    // Element {7}
    iPoint.atomID = 7;      // iatom
    iPoint.x = 23.8880006;  // x
    iPoint.y = 20.7029991;  // y
    iPoint.z = 23.8880006;  // z
    yCloud.pts.push_back(iPoint);
    // Element {8}
    iPoint.atomID = 8;      // iatom
    iPoint.x = 23.8880006;  // x
    iPoint.y = 27.073;      // y
    iPoint.z = 23.8880006;  // z
    yCloud.pts.push_back(iPoint);
    // Element {9}
    iPoint.atomID = 9;      // iatom
    iPoint.x = 27.073;      // x
    iPoint.y = 27.073;      // y
    iPoint.z = 20.7029991;  // z
    yCloud.pts.push_back(iPoint);
    // Element {10}
    iPoint.atomID = 10;     // iatom
    iPoint.x = 20.7029991;  // x
    iPoint.y = 23.8880006;  // y
    iPoint.z = 23.8880006;  // z
    yCloud.pts.push_back(iPoint);
    // Element {11}
    iPoint.atomID = 11;     // iatom
    iPoint.x = 27.073;      // x
    iPoint.y = 23.8880006;  // y
    iPoint.z = 23.8880006;  // z
    yCloud.pts.push_back(iPoint);
    // Element {12}
    iPoint.atomID = 12;     // iatom
    iPoint.x = 23.8880006;  // x
    iPoint.y = 23.8880006;  // y
    iPoint.z = 20.7029991;  // z
    yCloud.pts.push_back(iPoint);
    // Element {12}
    iPoint.atomID = 13;     // iatom
    iPoint.x = 23.8880006;  // x
    iPoint.y = 23.8880006;  // y
    iPoint.z = 27.073;      // z
    yCloud.pts.push_back(iPoint);
    // Update nop
    yCloud.nop = static_cast<int>(yCloud.pts.size());
    // Update the unordered map
    for (int iatom = 0; iatom < yCloud.nop; iatom++) {
      yCloud.idIndexMap[iatom] = iatom;
    }  // end of filling the map
    // --------------------
    WHEN("Given a pointCloud") {
      // Calculate a neighbour list
      nList = nneigh::neighListO(3.5, yCloud, 1);
      // Neighbour list by index
      nList = nneigh::neighbourListByIndex(yCloud, nList);
      // Find the vector of vector of rings
      rings = primitive::ringNetwork(nList, 7);
      THEN("There should be exactly one double-diamond cage.") {
        ringType.resize(
            rings.size());  // Has a value for each ring. init to zero.
        // Find the number of hexagonal cages
        listDDC = ring::findDDC(rings, ringType, listHC, cageList);
        // Assert the number of cages
        REQUIRE(cageList.size() ==
                1);  // Evaluate condition for a single tetragonal prism
      }
    }  // End of getting the neighbour list
  }    // End of given
} // End of scenario

// -- Additional tests for wrapper functions --

TEST_CASE("topoBulkAnalysis integration with HC geometry", "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_topobulk_analysis/").string();

  int ret = ring::topoBulkAnalysis(tmpPath, rings, nList, yCloud, 1, true);
  REQUIRE(ret == 0);

  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}


TEST_CASE("getStrucNumbers counts HC/DDC structures correctly", "[topo_bulk]") {
  // Set up known ring types: 5 rings, 2 basal HC, 3 prismatic
  std::vector<ring::strucType> ringType = {
      ring::strucType::HCbasal, ring::strucType::HCbasal,
      ring::strucType::HCprismatic, ring::strucType::HCprismatic,
      ring::strucType::HCprismatic};

  // One HC cage referencing rings 0-4
  cage::Cage hcCage;
  hcCage.type = cage::cageType::HexC;
  hcCage.rings = {0, 1, 2, 3, 4};
  std::vector<cage::Cage> cageList = {hcCage};

  int numHC = 0, numDDC = 0, mixedRings = 0, prismaticRings = 0,
      basalRings = 0;
  int ret = ring::getStrucNumbers(ringType, cageList, numHC, numDDC,
                                   mixedRings, prismaticRings, basalRings);

  REQUIRE(ret == 0);
  REQUIRE(numHC == 1);
  REQUIRE(numDDC == 0);
  REQUIRE(basalRings == 2);
  REQUIRE(prismaticRings == 3);
}

TEST_CASE("getAtomTypesTopoBulk assigns atom types from ring types",
          "[topo_bulk]") {
  // 6-membered rings with known atoms
  std::vector<std::vector<int>> rings = {{0, 1, 2, 3, 4, 5},
                                          {6, 7, 8, 9, 10, 11}};
  std::vector<ring::strucType> ringType = {ring::strucType::HCbasal,
                                            ring::strucType::DDC};
  std::vector<cage::iceType> atomTypes(12, cage::iceType::dummy);

  int ret = ring::getAtomTypesTopoBulk(rings, ringType, atomTypes);
  REQUIRE(ret == 0);

  // Atoms in ring 0 (HC basal) should be hc type
  REQUIRE(atomTypes[0] == cage::iceType::hc);
  // Atoms in ring 1 (DDC) should be ddc type
  REQUIRE(atomTypes[6] == cage::iceType::ddc);
}

TEST_CASE("findMixedRings identifies shared DDC/HC rings", "[topo_bulk]") {
  // Build the HC geometry, find rings, then test
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<ring::strucType> ringType(rings.size());
  std::vector<cage::Cage> cageList;
  auto listHC = ring::findHC(rings, ringType, nList, cageList);
  auto listDDC = ring::findDDC(rings, ringType, listHC, cageList);

  // Find mixed rings (with a single HC, there should be none)
  auto mixedList = ring::findMixedRings(rings, ringType, listDDC, listHC);
  REQUIRE(mixedList.empty());
}

TEST_CASE("buildRingSearchIndex lists every ring that contains an atom",
          "[topo_bulk]") {
  const std::vector<std::vector<int>> rings = {
      {0, 1, 2}, {2, 3, 4}, {0, 4, 5}, {9}};
  const auto index = ring::buildRingSearchIndex(rings, 6);

  REQUIRE(index.ringsContainingAtom.size() == 6);
  REQUIRE(index.ringsContainingAtom[0] == std::vector<int>{0, 2});
  REQUIRE(index.ringsContainingAtom[1] == std::vector<int>{0});
  REQUIRE(index.ringsContainingAtom[2] == std::vector<int>{0, 1});
  REQUIRE(index.ringsContainingAtom[3] == std::vector<int>{1});
  REQUIRE(index.ringsContainingAtom[4] == std::vector<int>{1, 2});
  REQUIRE(index.ringsContainingAtom[5] == std::vector<int>{2});
}

TEST_CASE("buildRingSearchIndex ignores atoms outside numAtoms",
          "[topo_bulk]") {
  const std::vector<std::vector<int>> rings = {{0, 1, 99}, {-1, 2}};
  const auto index = ring::buildRingSearchIndex(rings, 3);
  REQUIRE(index.ringsContainingAtom.size() == 3);
  REQUIRE(index.ringsContainingAtom[0] == std::vector<int>{0});
  REQUIRE(index.ringsContainingAtom[1] == std::vector<int>{0});
  REQUIRE(index.ringsContainingAtom[2] == std::vector<int>{1});
}

TEST_CASE("conditionOneDDC rejects an atom outside the ring index",
          "[topo_bulk]") {
  const std::vector<std::vector<int>> rings = {{0, 1, 2, 3, 4, 5}};
  ring::RingSearchIndex index;
  index.ringsContainingAtom.resize(3);
  std::vector<int> peripheral;
  REQUIRE_FALSE(ring::conditionOneDDC(rings, peripheral, 0, index));
}

TEST_CASE("findHC on the twelve-atom HC reports one HexC of five rings",
          "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<ring::strucType> ringType(rings.size());
  std::vector<cage::Cage> cageList;
  auto listHC = ring::findHC(rings, ringType, nList, cageList);

  REQUIRE(cageList.size() == 1);
  REQUIRE(cageList[0].type == cage::cageType::HexC);
  REQUIRE(cageList[0].rings.size() == 5);
  REQUIRE_FALSE(listHC.empty());
}

TEST_CASE("indexed findHC and findDDC match the convenience overloads",
          "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<std::vector<int>> sixRings;
  int maxAtom = 0;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      sixRings.push_back(r);
      for (int atom : r) {
        maxAtom = std::max(maxAtom, atom);
      }
    }
  }

  const auto index =
      ring::buildRingSearchIndex(sixRings, maxAtom + 1);

  std::vector<ring::strucType> typeWrap(sixRings.size());
  std::vector<ring::strucType> typeIdx(sixRings.size());
  std::vector<cage::Cage> cagesWrap;
  std::vector<cage::Cage> cagesIdx;

  auto hcWrap = ring::findHC(sixRings, typeWrap, nList, cagesWrap);
  auto hcIdx = ring::findHC(sixRings, typeIdx, nList, cagesIdx, index);
  REQUIRE(hcWrap == hcIdx);
  REQUIRE(typeWrap == typeIdx);

  auto ddcWrap = ring::findDDC(sixRings, typeWrap, hcWrap, cagesWrap);
  auto ddcIdx = ring::findDDC(sixRings, typeIdx, hcIdx, cagesIdx, index);
  REQUIRE(ddcWrap == ddcIdx);
  REQUIRE(typeWrap == typeIdx);
}

TEST_CASE("topoBulkAnalysis on mW cubic trajectory", "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_topobulk_mw/").string();
  int ret = ring::topoBulkAnalysis(tmpPath, rings, nList, yCloud, 1, true);
  REQUIRE(ret == 0);
  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}

TEST_CASE("bulkPolygonRingAnalysis on mW cubic trajectory", "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_bulkpoly_mw/").string();
  int ret = ring::bulkPolygonRingAnalysis(tmpPath, rings, nList, yCloud, 7, 1);
  REQUIRE(ret == 0);
  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}

// -- prism3 namespace tests --

TEST_CASE("basalPrismConditions returns true for connected basal rings",
          "[topo_bulk][prism3]") {
  // Build a prism: 8 atoms, two 4-membered basal rings connected
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10, 10, 50};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  double coords[8][3] = {{0, 0, 0},       {2.75, 0, 0},    {2.75, 2.75, 0},
                          {0, 2.75, 0},    {0, 0, 2.75},    {2.75, 0, 2.75},
                          {2.75, 2.75, 2.75}, {0, 2.75, 2.75}};
  for (int i = 0; i < 8; i++) {
    pt.atomID = i;
    pt.x = coords[i][0]; pt.y = coords[i][1]; pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;

  auto nList = nneigh::neighListO(3.5, cloud, 1);
  nList = nneigh::neighbourListByIndex(cloud, nList);

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  bool result = prism3::basalPrismConditions(nList, basal1, basal2);
  // Each atom in basal1 should have a neighbour in basal2
  REQUIRE(result == true);
}

TEST_CASE("basalPrismConditions returns false for disconnected rings",
          "[topo_bulk][prism3]") {
  // Two rings far apart with no bonds between them
  // Build nList manually: atoms 0-3 are in one plane, 4-7 in another far away
  // nList where atoms only know their own basal ring neighbours
  std::vector<std::vector<int>> nList = {
      {0, 1, 3}, {1, 0, 2}, {2, 1, 3}, {3, 0, 2},
      {4, 5, 7}, {5, 4, 6}, {6, 5, 7}, {7, 4, 6}};

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  bool result = prism3::basalPrismConditions(nList, basal1, basal2);
  REQUIRE(result == false);
}

TEST_CASE("relaxedPrismConditions returns true when at least one bond exists",
          "[topo_bulk][prism3]") {
  // nList where only atom 0 has atom 4 as neighbour (one bond between rings)
  std::vector<std::vector<int>> nList = {
      {0, 1, 3, 4}, {1, 0, 2}, {2, 1, 3}, {3, 0, 2},
      {4, 5, 7, 0}, {5, 4, 6}, {6, 5, 7}, {7, 4, 6}};

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  bool result = prism3::relaxedPrismConditions(nList, basal1, basal2);
  REQUIRE(result == true);
}

TEST_CASE("relaxedPrismConditions returns false when no bonds exist",
          "[topo_bulk][prism3]") {
  std::vector<std::vector<int>> nList = {
      {0, 1, 3}, {1, 0, 2}, {2, 1, 3}, {3, 0, 2},
      {4, 5, 7}, {5, 4, 6}, {6, 5, 7}, {7, 4, 6}};

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  bool result = prism3::relaxedPrismConditions(nList, basal1, basal2);
  REQUIRE(result == false);
}

TEST_CASE("basalRingsSeparation returns true for close rings",
          "[topo_bulk][prism3]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10, 10, 50};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  // Two parallel 4-rings separated by 2.75 in z
  double coords[8][3] = {{0, 0, 0},       {2.75, 0, 0},    {2.75, 2.75, 0},
                          {0, 2.75, 0},    {0, 0, 2.75},    {2.75, 0, 2.75},
                          {2.75, 2.75, 2.75}, {0, 2.75, 2.75}};
  for (int i = 0; i < 8; i++) {
    pt.atomID = i;
    pt.x = coords[i][0]; pt.y = coords[i][1]; pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  // Height cutoff 8.0 (default), separation is 2.75
  bool result = prism3::basalRingsSeparation(cloud, basal1, basal2, 8.0);
  REQUIRE(result == true);
}

TEST_CASE("basalRingsSeparation returns false for far apart rings",
          "[topo_bulk][prism3]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {100, 100, 100};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  // Two rings separated by 20 in z
  double coords[8][3] = {{0, 0, 0},    {2, 0, 0},    {2, 2, 0},    {0, 2, 0},
                          {0, 0, 20},   {2, 0, 20},   {2, 2, 20},   {0, 2, 20}};
  for (int i = 0; i < 8; i++) {
    pt.atomID = i;
    pt.x = coords[i][0]; pt.y = coords[i][1]; pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  bool result = prism3::basalRingsSeparation(cloud, basal1, basal2, 8.0);
  REQUIRE(result == false);
}

TEST_CASE("findBulkPrisms on ideal prism geometry", "[topo_bulk][prism3]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10, 10, 50};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  double coords[8][3] = {{0, 0, 0},       {2.75, 0, 0},    {2.75, 2.75, 0},
                          {0, 2.75, 0},    {0, 0, 2.75},    {2.75, 0, 2.75},
                          {2.75, 2.75, 2.75}, {0, 2.75, 2.75}};
  for (int i = 0; i < 8; i++) {
    pt.atomID = i;
    pt.x = coords[i][0]; pt.y = coords[i][1]; pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;

  auto nList = nneigh::neighListO(3.5, cloud, 1);
  nList = nneigh::neighbourListByIndex(cloud, nList);
  auto rings = primitive::ringNetwork(nList, 5);

  // Get only 4-membered rings
  auto fourRings = ring::getSingleRingSize(rings, 4);
  if (fourRings.size() >= 2) {
    std::vector<ring::strucType> ringType(fourRings.size(),
                                           ring::strucType::unclassified);
    std::vector<double> rmsdPerAtom(static_cast<size_t>(cloud.nop), -1.0);

    int ret = prism3::findBulkPrisms(fourRings, ringType, nList, cloud,
                                      rmsdPerAtom, 8.0);
    REQUIRE(ret == 0);
  }
}

TEST_CASE("findPrismatic identifies prismatic rings between HC basals",
          "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<ring::strucType> ringType(rings.size());
  std::vector<cage::Cage> cageList;
  auto listHC = ring::findHC(rings, ringType, nList, cageList);

  // The HC should have basal rings; find prismatic between them
  REQUIRE(cageList.size() == 1);
  if (cageList[0].rings.size() >= 2) {
    std::vector<int> prismaticRings;
    int ret = ring::findPrismatic(rings, listHC, ringType,
                                   cageList[0].rings[0], cageList[0].rings[1],
                                   prismaticRings);
    REQUIRE(ret == 0);
    // HC has 3 prismatic rings between its 2 basal rings
    REQUIRE(prismaticRings.size() == 3);
  }
}

TEST_CASE("DDC cage assembly is independent of ring enumeration order",
          "[topo_bulk]") {
  // The greedy assembly claims peripheral rings as it accepts cages, so
  // without a canonical visiting order the cage-object count depends on the
  // order the ring network was enumerated in. The per-ring classification
  // and the cage partition must both survive any permutation of the input.
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(idx, 7);
  std::vector<std::vector<int>> six;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }

  size_t referenceCages = 0;
  size_t referenceDDCRings = 0;
  for (int variant = 0; variant < 3; variant++) {
    auto shuffled = six;
    if (variant == 1) {
      std::sort(shuffled.begin(), shuffled.end());
    } else if (variant == 2) {
      std::sort(shuffled.rbegin(), shuffled.rend());
    }
    std::vector<ring::strucType> rt(shuffled.size(),
                                    ring::strucType::unclassified);
    std::vector<cage::Cage> cl;
    const auto index = ring::buildRingSearchIndex(shuffled, yCloud);
    auto hc = ring::findHC(shuffled, rt, idx, cl, index);
    auto ddc = ring::findDDC(shuffled, rt, hc, cl, index);
    if (variant == 0) {
      referenceCages = cl.size();
      referenceDDCRings = ddc.size();
      REQUIRE(referenceDDCRings == 8192);
    } else {
      REQUIRE(cl.size() == referenceCages);
      REQUIRE(ddc.size() == referenceDDCRings);
    }
  }
}

TEST_CASE("Cage counts are invariant under atom relabelling",
          "[topo_bulk]") {
  // The frame is relabelled at random and its neighbour list and rings are
  // rebuilt. The cage count is the size of a maximum non-conflicting set,
  // so it does not depend on the labels. The particular packing is the one
  // selected by the atom-index key.
  using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
  Cloud yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);
  const int n = yCloud.nop;

  using Cages = std::vector<std::pair<int, std::vector<std::vector<int>>>>;
  auto cagesOf = [](Cloud &cloud, const std::vector<int> &label) {
    auto nList = nneigh::neighListO(3.5, cloud, 1);
    auto idx = nneigh::neighbourListByIndex(cloud, nList);
    std::vector<std::vector<int>> six;
    for (const auto &r : primitive::ringNetwork(idx, 6)) {
      if (r.size() == 6) {
        six.push_back(r);
      }
    }
    const auto index = ring::buildRingSearchIndex(six, cloud);
    std::vector<ring::strucType> rt(six.size(), ring::strucType::unclassified);
    std::vector<cage::Cage> cl;
    const auto hc = ring::findHC(six, rt, idx, cl, index);
    (void)ring::findDDC(six, rt, hc, cl, index);
    Cages out;
    for (const auto &c : cl) {
      std::vector<std::vector<int>> faces;
      for (const int ri : c.rings) {
        std::vector<int> face;
        for (const int a : six[static_cast<size_t>(ri)]) {
          face.push_back(label[static_cast<size_t>(a)]);
        }
        std::sort(face.begin(), face.end());
        faces.push_back(face);
      }
      std::sort(faces.begin(), faces.end());
      out.emplace_back(c.type == cage::cageType::HexC ? 0 : 1, faces);
    }
    std::sort(out.begin(), out.end());
    return out;
  };

  std::vector<int> identity(static_cast<size_t>(n));
  std::iota(identity.begin(), identity.end(), 0);
  const Cages reference = cagesOf(yCloud, identity);
  REQUIRE(!reference.empty());
  int refHc = 0;
  int refDdc = 0;
  for (const auto &c : reference) {
    (c.first == 0 ? refHc : refDdc)++;
  }

  for (const unsigned seed : {1U, 2U, 3U}) {
    std::vector<int> perm(identity);
    std::mt19937 rng(seed);
    std::shuffle(perm.begin(), perm.end(), rng);
    Cloud relabelled = yCloud;
    std::vector<int> back(static_cast<size_t>(n));
    for (int a = 0; a < n; a++) {
      const auto to = static_cast<size_t>(perm[static_cast<size_t>(a)]);
      relabelled.pts[to] = yCloud.pts[static_cast<size_t>(a)];
      back[to] = a;
    }
    relabelled.idIndexMap.clear();
    for (int i = 0; i < n; i++) {
      relabelled.pts[static_cast<size_t>(i)].atomID = i + 1;
      relabelled.idIndexMap[i + 1] = i;
    }
    const Cages got = cagesOf(relabelled, back);
    int hc = 0;
    int ddc = 0;
    for (const auto &c : got) {
      (c.first == 0 ? hc : ddc)++;
    }
    // Cardinality is the maximum packing. The particular cages follow the
    // atom-index key, so they are not the same set after a relabeling.
    REQUIRE(hc == refHc);
    REQUIRE(ddc == refDdc);
    REQUIRE(got.size() == reference.size());
  }
}

TEST_CASE("Cage objects are invariant under a periodic translation",
          "[topo_bulk]") {
  // Shifting the same frame by a fraction of the box and rewrapping does
  // not change the bonds. The cage set is a maximum packing keyed by atom
  // indices and the unwrapped shift, so it does not move with the origin.
  using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
  auto cagesOf = [](Cloud &cloud) {
    auto nList = nneigh::neighListO(3.5, cloud, 1);
    auto idx = nneigh::neighbourListByIndex(cloud, nList);
    std::vector<std::vector<int>> six;
    for (const auto &r : primitive::ringNetwork(idx, 6)) {
      if (r.size() == 6) {
        six.push_back(r);
      }
    }
    const auto index = ring::buildRingSearchIndex(six, cloud);
    std::vector<ring::strucType> rt(six.size(), ring::strucType::unclassified);
    std::vector<cage::Cage> cl;
    const auto hc = ring::findHC(six, rt, idx, cl, index);
    (void)ring::findDDC(six, rt, hc, cl, index);
    std::vector<std::pair<int, std::vector<std::vector<int>>>> out;
    for (const auto &c : cl) {
      std::vector<std::vector<int>> faces;
      for (const int ri : c.rings) {
        auto face = six[static_cast<size_t>(ri)];
        std::sort(face.begin(), face.end());
        faces.push_back(std::move(face));
      }
      std::sort(faces.begin(), faces.end());
      out.emplace_back(c.type == cage::cageType::HexC ? 0 : 1, faces);
    }
    std::sort(out.begin(), out.end());
    return out;
  };
  Cloud cloud;
  cloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, cloud, 1);
  REQUIRE(cloud.nop == 4096);
  const auto reference = cagesOf(cloud);
  int ddc = 0;
  for (const auto &c : reference) {
    ddc += c.first == 1 ? 1 : 0;
  }
  // One equator for every two molecules: the conflict graph is the line
  // graph of a 4-regular bipartite graph on the molecules' rings, so the
  // maximum packing has n/2 cages.
  REQUIRE(ddc == 2048);
  const double frac[3] = {0.37, 0.21, 0.53};
  for (int k = 0; k < 3; k++) {
    const double lo = k < static_cast<int>(cloud.boxLow.size()) ? cloud.boxLow[static_cast<size_t>(k)] : 0.0;
    const double len = cloud.box[static_cast<size_t>(k)];
    for (int i = 0; i < cloud.nop; i++) {
      double *coord = k == 0 ? &cloud.pts[static_cast<size_t>(i)].x
                              : k == 1 ? &cloud.pts[static_cast<size_t>(i)].y
                                       : &cloud.pts[static_cast<size_t>(i)].z;
      double v = *coord - lo + frac[k] * len;
      v -= len * std::floor(v / len);
      *coord = lo + v;
    }
  }
  REQUIRE(cagesOf(cloud) == reference);
}

TEST_CASE("positionRank rejects two atoms at one wrapped position",
          "[topo_bulk]") {
  // An atom and the periodic image of another wrap onto one point, which no
  // position order separates, so the frame is rejected
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.pts.resize(3);
  const double xs[3] = {2.0, 1.0, 3.0};
  for (size_t i = 0; i < 3; i++) {
    cloud.pts[i].x = xs[i];
    cloud.pts[i].y = 1.0;
    cloud.pts[i].z = 1.0;
  }
  cloud.nop = 3;
  REQUIRE(ring::positionRank(cloud) == std::vector<int>{1, 0, 2});
  cloud.pts[2].x = 12.0;
  REQUIRE_THROWS_AS(ring::positionRank(cloud), std::invalid_argument);
}

namespace {

std::vector<ring::EquatorChoice>
conflictChoices(int n, const std::vector<std::pair<int, int>> &edges) {
  std::vector<std::vector<int>> adj(static_cast<size_t>(n));
  for (const auto &e : edges) {
    adj[static_cast<size_t>(e.first)].push_back(e.second);
    adj[static_cast<size_t>(e.second)].push_back(e.first);
  }
  std::vector<ring::EquatorChoice> periph(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    auto &row = periph[static_cast<size_t>(i)];
    if (adj[static_cast<size_t>(i)].empty()) {
      row.n = 1;
      row.ring[0] = -1;
      continue;
    }
    row.n = static_cast<int>(adj[static_cast<size_t>(i)].size());
    for (int k = 0; k < row.n; k++) {
      row.ring[k] = adj[static_cast<size_t>(i)][static_cast<size_t>(k)];
    }
  }
  return periph;
}

int bruteIndependent(int n, const std::vector<std::pair<int, int>> &edges) {
  int best = 0;
  const int limit = 1 << n;
  for (int mask = 0; mask < limit; mask++) {
    bool ok = true;
    for (const auto &e : edges) {
      if ((mask & (1 << e.first)) && (mask & (1 << e.second))) {
        ok = false;
        break;
      }
    }
    if (ok) {
      best = std::max(best, __builtin_popcount(static_cast<unsigned>(mask)));
    }
  }
  return best;
}

} // namespace

TEST_CASE("equatorial search is exact on K4 minus an edge", "[topo_bulk]") {
  for (int a = 0; a < 4; a++) {
    for (int b = a + 1; b < 4; b++) {
      std::vector<std::pair<int, int>> edges;
      for (int u = 0; u < 4; u++) {
        for (int v = u + 1; v < 4; v++) {
          if (u == a && v == b) {
            continue;
          }
          edges.push_back({u, v});
        }
      }
      const std::vector<int> order{0, 1, 2, 3};
      ring::PackingCertificate cert;
      const auto kept = ring::maximumEquators(conflictChoices(4, edges), order, &cert);
      REQUIRE(static_cast<int>(kept.size()) == 2);
      REQUIRE(cert.method == 2);
    }
  }
}

TEST_CASE("equatorial search matches brute force on random graphs", "[topo_bulk]") {
  std::mt19937 rng(20261102);
  std::uniform_real_distribution<double> coin(0.0, 1.0);
  int checked = 0;
  for (int n = 4; n <= 8; n++) {
    for (int sample = 0; sample < 20; sample++) {
      std::vector<std::pair<int, int>> edges;
      std::vector<int> deg(static_cast<size_t>(n), 0);
      for (int u = 0; u < n; u++) {
        for (int v = u + 1; v < n; v++) {
          if (coin(rng) < 0.45 && deg[static_cast<size_t>(u)] < 6 &&
              deg[static_cast<size_t>(v)] < 6) {
            edges.push_back({u, v});
            deg[static_cast<size_t>(u)]++;
            deg[static_cast<size_t>(v)]++;
          }
        }
      }
      std::vector<int> order(static_cast<size_t>(n));
      std::iota(order.begin(), order.end(), 0);
      const auto kept = ring::maximumEquators(conflictChoices(n, edges), order, nullptr);
      REQUIRE(static_cast<int>(kept.size()) == bruteIndependent(n, edges));
      checked++;
    }
  }
  REQUIRE(checked == 100);
}

namespace {

void mixCage(std::uint64_t &h, std::uint64_t v) {
  h ^= v;
  h *= 0x100000001b3ULL;
}

void mixCageInt(std::uint64_t &h, int v) {
  mixCage(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(v)));
}

// Same 8x8x8 ice Ic cell the bench driver times.
molSys::PointCloud<molSys::Point<double>, double> diamondIce4096() {
  constexpr double a = 6.37;
  constexpr double basis[8][3] = {
      {0.00, 0.00, 0.00}, {0.50, 0.50, 0.00}, {0.50, 0.00, 0.50},
      {0.00, 0.50, 0.50}, {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
      {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75},
  };
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  constexpr int cells = 8;
  const double L = a * static_cast<double>(cells);
  cloud.box = {L, L, L};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  int id = 0;
  for (int ix = 0; ix < cells; ++ix) {
    for (int iy = 0; iy < cells; ++iy) {
      for (int iz = 0; iz < cells; ++iz) {
        for (const auto &b : basis) {
          molSys::Point<double> pt;
          pt.type = 1;
          pt.molID = id;
          pt.atomID = id;
          pt.x = (static_cast<double>(ix) + b[0]) * a;
          pt.y = (static_cast<double>(iy) + b[1]) * a;
          pt.z = (static_cast<double>(iz) + b[2]) * a;
          cloud.pts.push_back(pt);
          cloud.idIndexMap[id] = id;
          ++id;
        }
      }
    }
  }
  cloud.nop = id;
  return cloud;
}

} // namespace

TEST_CASE("cubic ice cage census matches the pre-sort fingerprint",
          "[topo_bulk]") {
  auto cloud = diamondIce4096();
  REQUIRE(cloud.nop == 4096);
  const auto byId = nneigh::neighListO(3.5, cloud, 1);
  const auto byIndex = nneigh::neighbourListByIndex(cloud, byId);
  const auto rings = primitive::ringNetwork(byIndex, 6);
  const auto six = ring::getSingleRingSize(rings, 6);
  REQUIRE(six.size() == 8192);

  std::vector<ring::strucType> ringType(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cages;
  const auto index = ring::buildRingSearchIndex(six, cloud);
  const auto hc = ring::findHC(six, ringType, byIndex, cages, index);
  const auto ddc = ring::findDDC(six, ringType, hc, cages, index);

  std::uint64_t h = 0xcbf29ce484222325ULL;
  mixCageInt(h, static_cast<int>(hc.size()));
  for (const int id : hc) {
    mixCageInt(h, id);
  }
  mixCageInt(h, static_cast<int>(ddc.size()));
  for (const int id : ddc) {
    mixCageInt(h, id);
  }
  mixCageInt(h, static_cast<int>(cages.size()));
  for (const auto &cageUnit : cages) {
    mixCageInt(h, static_cast<int>(cageUnit.type));
    for (const int ringId : cageUnit.rings) {
      mixCageInt(h, ringId);
    }
  }
  for (const auto kind : ringType) {
    mixCageInt(h, static_cast<int>(kind));
  }
  // Recorded from the cage search before the unique-count and cycle-key
  // rewrite. A changed assignment fails this, a faster path does not.
  REQUIRE(h == 0xd88fddd29c1c1fddULL);
}
