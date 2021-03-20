// Internal
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <topo_one_dim.hpp>

// Standard
#include <iostream>

#include <catch2/catch.hpp>
#include <rang.hpp>

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
      nList = nneigh::neighListO(3.5, &yCloud, 1);
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
        listPrism = ring::findPrisms(rings, &ringType, &nPrisms, &nDeformed,
                                     nList, &yCloud, &rmsdPerAtom, false);
        // Assert the number of prism blocks
        REQUIRE(nPrisms ==
                1); // Evaluate condition for a single tetragonal prism
      }
    } // End of getting the neighbour list
  }   // End of given
} // End of scenario
