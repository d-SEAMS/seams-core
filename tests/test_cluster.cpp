#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
#include <cluster.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <filesystem>
#include <unordered_map>

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
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);

  std::vector<int> linkedList;
  int result = clump::singleClusterLinkedList(cloud, nList, linkedList);

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

  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);

  std::vector<int> linkedList;
  clump::singleClusterLinkedList(cloud, nList, linkedList);

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

  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);
  int result = clump::recenterClusterCloud(cloud, nList);

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

TEST_CASE("largestIceCluster identifies clusters and writes stats",
          "[cluster]") {
  // Build a system with a few "ice" atoms close together and one "water" atom
  // far away
  auto cloud = makeClusterCloud(6, 1.0);
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);

  // Mark first 4 atoms as ice, last 2 as water
  std::vector<bool> isIce = {true, true, true, true, false, false};

  // Build an "ice" PointCloud from ice atoms
  molSys::PointCloud<molSys::Point<double>, double> iceCloud;
  iceCloud.box = cloud.box;
  iceCloud.boxLow = cloud.boxLow;
  iceCloud.currentFrame = cloud.currentFrame;
  for (int i = 0; i < cloud.nop; i++) {
    if (isIce[i]) {
      iceCloud.pts.push_back(cloud.pts[i]);
      iceCloud.idIndexMap[cloud.pts[i].atomID] =
          static_cast<int>(iceCloud.pts.size()) - 1;
    }
  }
  iceCloud.nop = static_cast<int>(iceCloud.pts.size());

  // Build an ID-based neighbour list for the ice cloud
  auto iceNList = nneigh::neighListO(1.5, iceCloud, 1);

  std::vector<int> linkedList;
  std::vector<int> nClusters;
  std::unordered_map<int, int> indexNumber;

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_largestcluster/").string();

  int ret = clump::largestIceCluster(tmpPath, cloud, iceCloud, iceNList, isIce,
                                      linkedList, nClusters, indexNumber, 1);

  REQUIRE(ret == 0);

  std::error_code _ec_; std::filesystem::remove_all(tmpPath, _ec_);
}

TEST_CASE("singleClusterLinkedList with single atom", "[cluster]") {
  auto cloud = makeClusterCloud(1, 1.0);
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);

  std::vector<int> linkedList;
  int result = clump::singleClusterLinkedList(cloud, nList, linkedList);

  REQUIRE(result == 0);
  REQUIRE(linkedList.size() == 1);
  REQUIRE(linkedList[0] == 0); // Points to itself
}

TEST_CASE("clusterAnalysis with Q6 on mW cubic trajectory", "[cluster]") {
  // Read O atoms from mW cubic ice
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);

  // Build neighbour list
  auto nList = nneigh::neighListO(3.5, yCloud, 1);

  // Run CHILL correlation + ice type classification
  yCloud = chill::getCorrel(yCloud, nList, false);
  yCloud = chill::getIceTypeNoPrint(yCloud, nList, false);

  // clusterAnalysis expects the full yCloud and returns the largest ice cluster
  molSys::PointCloud<molSys::Point<double>, double> iceCloud;
  std::vector<std::vector<int>> iceNeighbourList;
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_clusteranalysis/").string();

  int ret = clump::clusterAnalysis(tmpPath, iceCloud, yCloud, nList,
                                    iceNeighbourList, 3.5, 1, "q6");
  REQUIRE(ret == 0);
  // mW cubic ice should produce a non-trivial ice cluster
  REQUIRE(iceCloud.nop > 0);

  std::error_code _ec_; std::filesystem::remove_all(tmpPath, _ec_);
}
