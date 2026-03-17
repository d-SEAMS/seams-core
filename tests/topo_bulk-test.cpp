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
#include <filesystem>
#include <iostream>

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
    yCloud.nop = yCloud.pts.size();
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
    yCloud.nop = yCloud.pts.size();
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

  std::string tmpPath = "/tmp/dseams_test_topobulk_analysis/";

  int ret = ring::topoBulkAnalysis(tmpPath, rings, nList, yCloud, 1, true);
  REQUIRE(ret == 0);

  fs::remove_all(tmpPath);
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

TEST_CASE("topoBulkAnalysis on mW cubic trajectory", "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = "/tmp/dseams_test_topobulk_mw/";
  int ret = ring::topoBulkAnalysis(tmpPath, rings, nList, yCloud, 1, true);
  REQUIRE(ret == 0);
  fs::remove_all(tmpPath);
}

TEST_CASE("bulkPolygonRingAnalysis on mW cubic trajectory", "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = "/tmp/dseams_test_bulkpoly_mw/";
  int ret = ring::bulkPolygonRingAnalysis(tmpPath, rings, nList, yCloud, 7, 1);
  REQUIRE(ret == 0);
  fs::remove_all(tmpPath);
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
