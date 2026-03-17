#include <catch2/catch_test_macros.hpp>

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>

#include <vector>

// Helper: build a small cloud with known atom IDs for populateGraphFromNListID
static molSys::PointCloud<molSys::Point<double>, double> makeSmallCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  // Triangle: 3 atoms mutually within cutoff
  double coords[3][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
  for (int i = 0; i < 3; i++) {
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
  cloud.nop = 3;
  return cloud;
}

// -- populateGraphFromNListID tests --

TEST_CASE("populateGraphFromNListID builds graph from ID-based nList",
          "[franzblau]") {
  auto cloud = makeSmallCloud();
  // Build nList by atom ID (neighListO format: first elem is self ID)
  auto nList = nneigh::neighListO(1.5, cloud, 1);
  REQUIRE(nList.size() == 3);

  auto graph = primitive::populateGraphFromNListID(cloud, nList);

  // Graph should have 3 vertices
  REQUIRE(graph.pts.size() == 3);
  // Each vertex should have the correct atomIndex
  for (int i = 0; i < 3; i++) {
    REQUIRE(graph.pts[i].atomIndex == i);
  }
  // Atom 0 should list atoms 1 and 2 as neighbours (by index)
  REQUIRE(graph.pts[0].neighListIndex.size() == 2);
}

// -- clearGraph tests --

TEST_CASE("clearGraph empties a graph", "[franzblau]") {
  // Build a graph with some data
  std::vector<std::vector<int>> nList = {{0, 1, 2}, {1, 0, 2}, {2, 0, 1}};
  auto graph = primitive::populateGraphFromIndices(nList);
  REQUIRE(graph.pts.size() == 3);

  // Add a synthetic ring
  graph.rings.push_back({0, 1, 2});
  REQUIRE(graph.rings.size() == 1);

  auto cleared = primitive::clearGraph(graph);

  // After clearing, the original graph should be empty (swapped)
  REQUIRE(graph.pts.empty());
  REQUIRE(graph.rings.empty());
}

// -- ring::compareRings tests --

TEST_CASE("compareRings returns true for same elements in different order",
          "[ring]") {
  std::vector<int> ring1 = {3, 1, 2, 0};
  std::vector<int> ring2 = {0, 1, 2, 3};
  REQUIRE(ring::compareRings(ring1, ring2) == true);
}

TEST_CASE("compareRings returns false for different elements", "[ring]") {
  std::vector<int> ring1 = {0, 1, 2};
  std::vector<int> ring2 = {0, 1, 3};
  REQUIRE(ring::compareRings(ring1, ring2) == false);
}

TEST_CASE("compareRings returns false for different sizes", "[ring]") {
  std::vector<int> ring1 = {0, 1, 2};
  std::vector<int> ring2 = {0, 1, 2, 3};
  REQUIRE(ring::compareRings(ring1, ring2) == false);
}

TEST_CASE("compareRings with empty rings returns true", "[ring]") {
  std::vector<int> ring1 = {};
  std::vector<int> ring2 = {};
  REQUIRE(ring::compareRings(ring1, ring2) == true);
}
