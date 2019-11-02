// Internal
#include <cage.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <topo_bulk.hpp>

// Standard
#include <iostream>

// Conan
#include <catch2/catch.hpp>
#include <rang.hpp>

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
      nList = nneigh::neighListO(3.5, &yCloud, 1);
      // Neighbour list by index
      nList = nneigh::neighbourListByIndex(&yCloud, nList);
      // Find the vector of vector of rings
      rings = primitive::ringNetwork(nList, 7);
      THEN("There should be exactly one hexagonal cage.") {
        ringType.resize(
            rings.size());  // Has a value for each ring. init to zero.
        // Find the number of hexagonal cages
        listHC = ring::findHC(rings, &ringType, nList, &cageList);
        // Assert the number of cages
        REQUIRE(cageList.size() ==
                1);  // Evaluate condition for a single tetragonal prism
      }
    }  // End of getting the neighbour list
  }    // End of given
}  // End of scenario
