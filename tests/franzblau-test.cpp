// Internal
#include <franzblau.hpp>

// Standard
#include <iostream>

#include <catch2/catch.hpp>
#include <rang.hpp>

SCENARIO(
    "Test the number of rings formed when there is one 4-membered ring and one "
    "3-membered ring.",
    "[ring]") {
  GIVEN("An adjacency or neighbour list") {
    // Hard-coded example of a known system
    // There are 2 franzblau rings
    // For simplicity, we write everything in terms of indices
    std::vector<std::vector<int>>
        nList;                // Vector of vectors of the neighbour list
    std::vector<int> iNeigh;  // Neighbours of iatom
    primitive::Graph
        fullGraph;  // Graph object that will contain the information of nList
    int maxDepth = 7;  // Maximum depth upto which the search will be conducted
    int nRings;        // Total number of rings
    // --------------------
    // Building the neighbour list
    // {element} #neighbors [neighbor list]
    // {0} 2 [1 2 ]
    // {1} 2 [0 2 ]
    // {2} 2 [0 1 ]
    // {3} 2 [4 6 ]
    // {4} 2 [3 5 ]
    // {5} 2 [4 6 ]
    // {6} 2 [3 5 ]
    //
    // Element {0}
    iNeigh.push_back(0);      // iatom
    iNeigh.push_back(1);      // neighbour
    iNeigh.push_back(2);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {1}
    iNeigh.clear();
    iNeigh.push_back(1);      // iatom
    iNeigh.push_back(0);      // neighbour
    iNeigh.push_back(2);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {2}
    iNeigh.clear();
    iNeigh.push_back(2);      // iatom
    iNeigh.push_back(0);      // neighbour
    iNeigh.push_back(1);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {3}
    iNeigh.clear();
    iNeigh.push_back(3);      // iatom
    iNeigh.push_back(4);      // neighbour
    iNeigh.push_back(6);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {4}
    iNeigh.clear();
    iNeigh.push_back(4);      // iatom
    iNeigh.push_back(3);      // neighbour
    iNeigh.push_back(5);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {5}
    iNeigh.clear();
    iNeigh.push_back(5);      // iatom
    iNeigh.push_back(4);      // neighbour
    iNeigh.push_back(6);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {6}
    iNeigh.clear();
    iNeigh.push_back(6);      // iatom
    iNeigh.push_back(3);      // neighbour
    iNeigh.push_back(5);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // --------------------
    WHEN("Calculate every possible ring using backtracking") {
      // Calculate the number of rings
      fullGraph = primitive::countAllRingsFromIndex(nList, maxDepth);
      THEN("The number of rings should be 2 for this case.") {
        // Find the number of rings inside the graph object
        nRings = fullGraph.rings.size();  // Number of rings
        REQUIRE(nRings == 2);             // Evaluate condition
      }
    }  // End of getting all the rings (non-primitive included)
  }    // End of given
}  // End of scenario

SCENARIO(
    "Test the number of rings formed when there are 8 Franzblau rings and 5 "
    "primitive rings.",
    "[ring]") {
  GIVEN(
      "A known adjacency or neighbour list with the number of all rings "
      "obtained from backtracking not equal to the number of primitive rings") {
    // Hard-coded example of a known system
    // There are 2 franzblau rings
    // For simplicity, we write everything in terms of indices
    std::vector<std::vector<int>>
        nList;                // Vector of vectors of the neighbour list
    std::vector<int> iNeigh;  // Neighbours of iatom
    primitive::Graph
        fullGraph;  // Graph object that will contain the information of nList
    int maxDepth = 7;  // Maximum depth upto which the search will be conducted
    int nAllRings;     // Total number of rings from backtracking
    int nPrimitiveRings;  // Number of primitive rings
    // --------------------
    // Building the neighbour list
    // {element} #neighbors [neighbor list]
    // {0} 2 [1 2 ]
    // {1} 2 [0 2 ]
    // {2} 2 [0 1 ]
    // {3} 3 [4 5 6 ]
    // {4} 3 [3 5 6 ]
    // {5} 3 [3 4 6 ]
    // {6} 3 [3 4 5 ]
    //
    // Element {0}
    iNeigh.push_back(0);      // iatom
    iNeigh.push_back(1);      // neighbour
    iNeigh.push_back(2);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {1}
    iNeigh.clear();
    iNeigh.push_back(1);      // iatom
    iNeigh.push_back(0);      // neighbour
    iNeigh.push_back(2);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {2}
    iNeigh.clear();
    iNeigh.push_back(2);      // iatom
    iNeigh.push_back(0);      // neighbour
    iNeigh.push_back(1);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {3}
    iNeigh.clear();
    iNeigh.push_back(3);      // iatom
    iNeigh.push_back(4);      // neighbour
    iNeigh.push_back(5);      // neighbour
    iNeigh.push_back(6);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {4}
    iNeigh.clear();
    iNeigh.push_back(4);      // iatom
    iNeigh.push_back(3);      // neighbour
    iNeigh.push_back(5);      // neighbour
    iNeigh.push_back(6);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {5}
    iNeigh.clear();
    iNeigh.push_back(5);      // iatom
    iNeigh.push_back(3);      // neighbour
    iNeigh.push_back(4);      // neighbour
    iNeigh.push_back(6);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // Element {6}
    iNeigh.clear();
    iNeigh.push_back(6);      // iatom
    iNeigh.push_back(3);      // neighbour
    iNeigh.push_back(4);      // neighbour
    iNeigh.push_back(5);      // neighbour
    nList.push_back(iNeigh);  // Add to neighbour list
    // --------------------
    WHEN(
        "Calculating primitive rings according to Phys. Rev. B, "
        "44(10):4925-4930 (1991).") {
      // Calculate the number of rings
      fullGraph = primitive::countAllRingsFromIndex(nList, maxDepth);
      THEN("The number of primitive rings should be 5 for this case.") {
        // Find the number of rings inside the graph object
        nAllRings =
            fullGraph.rings.size();  // Number of rings after backtracking
        // The number of all rings from backtracking should be 8
        REQUIRE(nAllRings == 8);
        // Get rid of all non-SP rings
        fullGraph = primitive::removeNonSPrings(&fullGraph);
        // Get the number of primitive rings
        nPrimitiveRings = fullGraph.rings.size();
        // For this case, the number of primitive rings should be 5
        REQUIRE(nPrimitiveRings == 5);  // Evaluate condition
      }
    }  // End of getting all the rings (non-primitive included)
  }    // End of given
} // End of scenario
