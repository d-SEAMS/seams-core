#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cluster.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

// Helper: build a PointCloud with some "ice" and some "water" particles
// to test the clustering machinery at the linked-list level.
// We test singleClusterLinkedList since largestIceCluster and clusterAnalysis
// require file I/O (sout::writeClusterStats), making them harder to unit test.
static molSys::PointCloud<molSys::Point<double>, double>
makeClusterCloud(int nAtoms, double spacing) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {100.0, 100.0, 100.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  // Place atoms in a line along x, separated by spacing
  for (int i = 0; i < nAtoms; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.molID = i;
    pt.x = spacing * i;
    pt.y = 0.0;
    pt.z = 0.0;
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = nAtoms;
  return cloud;
}

TEST_CASE("singleClusterLinkedList groups all connected atoms", "[cluster]") {
  // 4 atoms in a line with spacing 1.0; cutoff 1.5 connects consecutive pairs
  auto cloud = makeClusterCloud(4, 1.0);
  auto nList = nneigh::getNewNeighbourListByIndex(&cloud, 1.5);

  std::vector<int> linkedList;
  int result = clump::singleClusterLinkedList(&cloud, nList, &linkedList);

  REQUIRE(result == 0);
  REQUIRE(linkedList.size() == 4);

  // In a linked list for a single cluster, following the links from any
  // start should visit all 4 atoms before returning to start.
  int start = 0;
  int current = start;
  int count = 0;
  do {
    current = linkedList[current];
    count++;
  } while (current != start && count <= 4);

  // All 4 atoms should be in one cluster
  REQUIRE(count == 4);
  REQUIRE(current == start);
}

TEST_CASE("singleClusterLinkedList with disconnected atom", "[cluster]") {
  // 3 atoms close together + 1 far away
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {100.0, 100.0, 100.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  double coords[4][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {50.0, 50.0, 50.0}};

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

  auto nList = nneigh::getNewNeighbourListByIndex(&cloud, 1.5);

  std::vector<int> linkedList;
  clump::singleClusterLinkedList(&cloud, nList, &linkedList);

  REQUIRE(linkedList.size() == 4);

  // The isolated atom (index 3) should point to itself
  // (no neighbour to swap with)
  REQUIRE(linkedList[3] == 3);
}

TEST_CASE("recenterClusterCloud shifts centroid toward box center",
          "[cluster]") {
  // Build a small cluster offset from the box center
  auto cloud = makeClusterCloud(3, 1.0);
  // Atoms at x=0,1,2; box center is at 50

  auto nList = nneigh::getNewNeighbourListByIndex(&cloud, 1.5);
  int result = clump::recenterClusterCloud(&cloud, nList);

  REQUIRE(result == 0);

  // After recentering, the centroid should be at the box center (50,50,50)
  double cx = 0.0, cy = 0.0, cz = 0.0;
  for (int i = 0; i < cloud.nop; i++) {
    cx += cloud.pts[i].x;
    cy += cloud.pts[i].y;
    cz += cloud.pts[i].z;
  }
  cx /= cloud.nop;
  cy /= cloud.nop;
  cz /= cloud.nop;

  REQUIRE_THAT(cx, Catch::Matchers::WithinAbs(50.0, 1e-8));
  REQUIRE_THAT(cy, Catch::Matchers::WithinAbs(50.0, 1e-8));
  REQUIRE_THAT(cz, Catch::Matchers::WithinAbs(50.0, 1e-8));
}
