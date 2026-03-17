// Internal
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>
#include <topo_one_dim.hpp>

// Standard
#include <filesystem>
#include <iostream>

#include <catch2/catch_test_macros.hpp>
#include <rang.hpp>

namespace fs = std::filesystem;

// Helper: build the 8-atom prism geometry (reusable)
static void buildPrismCloud(molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  yCloud.box = {10, 10, 50};
  yCloud.boxLow = {0.0, 0.0, 0.0};
  yCloud.currentFrame = 1;
  molSys::Point<double> iPoint;
  iPoint.type = 1;
  double coords[8][3] = {
      {0, 0, 0}, {2.75, 0, 0}, {2.75, 2.75, 0}, {0, 2.75, 0},
      {0, 0, 2.75}, {2.75, 0, 2.75}, {2.75, 2.75, 2.75}, {0, 2.75, 2.75}};
  for (int i = 0; i < 8; i++) {
    iPoint.atomID = i;
    iPoint.x = coords[i][0]; iPoint.y = coords[i][1]; iPoint.z = coords[i][2];
    yCloud.pts.push_back(iPoint);
    yCloud.idIndexMap[i] = i;
  }
  yCloud.nop = 8;
}

SCENARIO("Test the prism identification scheme for a single tetragonal prism.",
         "[topo]") {
  GIVEN("A pointCloud") {
    // Hard-coded example of a single tetragonal prism
    // Here, we manually input a tetragonal prism's coordinates
    molSys::PointCloud<molSys::Point<double>, double> yCloud; // pointCloud
    molSys::Point<double> iPoint;                             // A single point
    std::vector<double> box = {10, 10, 50};                   // Box lengths
    yCloud.box = box; // Update the box lengths. z dim is axial
    std::vector<std::vector<int>> nList; // Neighbour list
    std::vector<std::vector<int>> rings; // Rings
    std::vector<ring::strucType>
        ringType; // This vector will have a value for each ring inside
    std::vector<int> listPrism; // Contains indices of rings in each prism
    // --------------------
    // Building the prism block (Assuming an ideal tetragonal prism)
    // Each side length is 2.75
    //
    iPoint.type = 1; // Same for all the points here
    // Element {0}
    iPoint.atomID = 0; // iatom
    iPoint.x = 0;      // x
    iPoint.y = 0;      // y
    iPoint.z = 0;      // z
    yCloud.pts.push_back(iPoint);
    // Element {1}
    iPoint.atomID = 1; // iatom
    iPoint.x = 2.75;   // x
    iPoint.y = 0;      // y
    iPoint.z = 0;      // z
    yCloud.pts.push_back(iPoint);
    // Element {2}
    iPoint.atomID = 2; // iatom
    iPoint.x = 2.75;   // x
    iPoint.y = 2.75;   // y
    iPoint.z = 0;      // z
    yCloud.pts.push_back(iPoint);
    // Element {3}
    iPoint.atomID = 3; // iatom
    iPoint.x = 0;      // x
    iPoint.y = 2.75;   // y
    iPoint.z = 0;      // z
    yCloud.pts.push_back(iPoint);
    // Element {4}
    iPoint.atomID = 4; // iatom
    iPoint.x = 0;      // x
    iPoint.y = 0;      // y
    iPoint.z = 2.75;   // z
    yCloud.pts.push_back(iPoint);
    // Element {5}
    iPoint.atomID = 5; // iatom
    iPoint.x = 2.75;   // x
    iPoint.y = 0;      // y
    iPoint.z = 2.75;   // z
    yCloud.pts.push_back(iPoint);
    // Element {6}
    iPoint.atomID = 6; // iatom
    iPoint.x = 2.75;   // x
    iPoint.y = 2.75;   // y
    iPoint.z = 2.75;   // z
    yCloud.pts.push_back(iPoint);
    // Element {7}
    iPoint.atomID = 7; // iatom
    iPoint.x = 0;      // x
    iPoint.y = 2.75;   // y
    iPoint.z = 2.75;   // z
    yCloud.pts.push_back(iPoint);
    // Update nop
    yCloud.nop = yCloud.pts.size();
    // Update the unordered map
    for (int iatom = 0; iatom < yCloud.nop; iatom++) {
      yCloud.idIndexMap[iatom] = iatom;
    } // end of filling the map
    // --------------------
    WHEN("Given a pointCloud, and a neighbour list") {
      // Calculate a neighbour list
      nList = nneigh::neighListO(3.5, yCloud, 1);
      // Find the vector of vector of rings
      rings = primitive::ringNetwork(nList, 5);
      THEN("There should be exactly one tetragonal prism.") {
        int nPrisms = 0; // The number of perfect prisms
        int nDeformed = 0;
        // Qualifier for the type of atom it is:
        std::vector<double> rmsdPerAtom; // Not used here
        ringType.resize(
            rings.size()); // Has a value for each ring. init to zero.
        // Find the number of tetragonal prisms
        listPrism = ring::findPrisms(rings, ringType, nPrisms, nDeformed,
                                     nList, yCloud, rmsdPerAtom, false);
        // Assert the number of prism blocks
        REQUIRE(nPrisms ==
                1); // Evaluate condition for a single tetragonal prism
      }
    } // End of getting the neighbour list
  }   // End of given
} // End of scenario

// -- Additional tests for wrapper functions --

TEST_CASE("prismAnalysis integration with tetragonal prism", "[topo_one_dim]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildPrismCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto rings = primitive::ringNetwork(nList, 5);

  std::string tmpPath = "/tmp/dseams_test_prismanalysis/";
  int atomID = 0;

  int ret = ring::prismAnalysis(tmpPath, rings, nList, yCloud, 5, atomID, 1, 1,
                                 false);
  REQUIRE(ret == 0);

  fs::remove_all(tmpPath);
}


TEST_CASE("relaxedPrismConditions detects prism with relaxed criteria",
          "[topo_one_dim]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildPrismCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 5);

  // Find two 4-membered rings (basal rings of the tetragonal prism)
  // The prism should have basal rings
  std::vector<std::vector<int>> fourRings;
  for (const auto &r : rings) {
    if (r.size() == 4) fourRings.push_back(r);
  }

  if (fourRings.size() >= 2) {
    bool relaxed =
        ring::relaxedPrismConditions(nList, fourRings[0], fourRings[1]);
    // At least one bond between basal rings should exist
    REQUIRE(relaxed == true);
  }
}

TEST_CASE("assignPrismType assigns atom types for known prism list",
          "[topo_one_dim]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildPrismCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto rings = primitive::ringNetwork(nList, 5);

  std::vector<ring::strucType> ringType(rings.size());
  int nPrisms = 0, nDeformed = 0;
  std::vector<double> rmsdPerAtom;
  auto listPrism = ring::findPrisms(rings, ringType, nPrisms, nDeformed,
                                     nList, yCloud, rmsdPerAtom, false);

  // Now assign prism types
  std::vector<int> atomTypes(yCloud.nop, 1); // 1 = dummy
  std::vector<ring::strucType> atomState(yCloud.nop);

  int ret = ring::assignPrismType(rings, listPrism, 4, ringType, atomTypes,
                                   atomState);
  REQUIRE(ret == 0);

  // At least some atoms should be classified as prism type (not dummy)
  bool hasNonDummy = false;
  for (int i = 0; i < yCloud.nop; i++) {
    if (atomTypes[i] != 1) hasNonDummy = true;
  }
  REQUIRE(hasNonDummy);
}
