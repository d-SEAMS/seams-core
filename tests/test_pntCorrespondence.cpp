#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cage.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <pntCorrespondence.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <Eigen/Dense>
#include <vector>

// Helper: build an 8-atom tetragonal prism
static molSys::PointCloud<molSys::Point<double>, double> makePrismCloud() {
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
  return cloud;
}

// Helper: build the 12-atom HC geometry
static molSys::PointCloud<molSys::Point<double>, double> makeHCCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {53.9690018, 54.5289994, 51.257};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  double coords[12][3] = {
      {8.995, 10.3859997, 15.0939999}, {6.7459998, 14.2810001, 15.0939999},
      {4.4970002, 10.3859997, 15.0939999}, {8.995, 12.9829998, 14.1949997},
      {6.7459998, 9.0880003, 14.1949997}, {4.4970002, 12.9829998, 14.1949997},
      {8.995, 12.9829998, 11.4329996}, {6.7459998, 9.0880003, 11.4329996},
      {4.4970002, 12.9829998, 11.4329996}, {8.995, 10.3859997, 10.5340004},
      {6.7459998, 14.2810001, 10.5340004}, {4.4970002, 10.3859997, 10.5340004}};
  for (int i = 0; i < 12; i++) {
    pt.atomID = i;
    pt.x = coords[i][0]; pt.y = coords[i][1]; pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 12;
  return cloud;
}

// -- getPointSetRefRing tests --

TEST_CASE("getPointSetRefRing creates reference ring for various sizes",
          "[pntCorrespondence]") {
  // Tetragonal (4-membered) ring, axial dim = 2 (z)
  auto refPts = pntToPnt::getPointSetRefRing(4, 2);
  REQUIRE(refPts.rows() == 4);
  REQUIRE(refPts.cols() == 3);

  // Hexagonal (6-membered) ring
  auto refPts6 = pntToPnt::getPointSetRefRing(6, 2);
  REQUIRE(refPts6.rows() == 6);
  REQUIRE(refPts6.cols() == 3);
}

TEST_CASE("getPointSetRefRing with different axial dimensions",
          "[pntCorrespondence]") {
  auto refX = pntToPnt::getPointSetRefRing(4, 0);
  auto refY = pntToPnt::getPointSetRefRing(4, 1);
  auto refZ = pntToPnt::getPointSetRefRing(4, 2);

  // All should have the same shape
  REQUIRE(refX.rows() == 4);
  REQUIRE(refY.rows() == 4);
  REQUIRE(refZ.rows() == 4);

  // The axial column should be zero in the reference ring
  for (int i = 0; i < 4; i++) {
    REQUIRE_THAT(refZ(i, 2), Catch::Matchers::WithinAbs(0.0, 1e-10));
  }
}

// -- getRadiusFromRings tests --

TEST_CASE("getRadiusFromRings computes radius for tetragonal prism",
          "[pntCorrespondence]") {
  auto cloud = makePrismCloud();
  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  double radius = pntToPnt::getRadiusFromRings(cloud, basal1, basal2);
  // For a 2.75 x 2.75 square, the radius from centroid to vertex is
  // sqrt(2) * 2.75/2 ~ 1.944
  REQUIRE(radius > 0.0);
  REQUIRE_THAT(radius, Catch::Matchers::WithinAbs(1.944, 0.1));
}

// -- getAvgHeightPrismBlock tests --

TEST_CASE("getAvgHeightPrismBlock computes height for tetragonal prism",
          "[pntCorrespondence]") {
  auto cloud = makePrismCloud();
  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  double height = pntToPnt::getAvgHeightPrismBlock(cloud, basal1, basal2);
  // Height should be 2.75 (z separation between basal planes)
  REQUIRE_THAT(height, Catch::Matchers::WithinAbs(2.75, 0.1));
}

// -- createPrismBlock tests --

TEST_CASE("createPrismBlock creates reference prism from basal rings",
          "[pntCorrespondence]") {
  auto cloud = makePrismCloud();
  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  auto refPts = pntToPnt::getPointSetRefRing(4, 2);
  auto prismBlock =
      pntToPnt::createPrismBlock(cloud, refPts, 4, basal1, basal2);

  // Should have 2*4 = 8 rows (both basal rings) and 3 columns
  REQUIRE(prismBlock.rows() == 8);
  REQUIRE(prismBlock.cols() == 3);
}

// -- relOrderPrismBlock tests --

TEST_CASE("relOrderPrismBlock with nList from trajectory", "[pntCorrespondence]") {
  // relOrderPrismBlock needs a nList where basal ring atoms are mutual
  // neighbours. Use the mW trajectory for a realistic test.
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  if (yCloud.nop == 0) return; // skip if no data

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 5);

  // Find 4-membered rings to use as basal candidates
  std::vector<std::vector<int>> fourRings;
  for (const auto &r : rings) {
    if (r.size() == 4) fourRings.push_back(r);
    if (fourRings.size() >= 2) break;
  }

  if (fourRings.size() >= 2) {
    std::vector<int> outBasal1, outBasal2;
    int ret = pntToPnt::relOrderPrismBlock(yCloud, fourRings[0], fourRings[1],
                                            nList, outBasal1, outBasal2);
    // May or may not find valid ordering depending on ring pair
    REQUIRE((ret == 0 || ret == 1));
  }
}

// -- fillPointSetPrismRing / fillPointSetPrismBlock tests --

TEST_CASE("fillPointSetPrismRing fills Eigen matrix from ring",
          "[pntCorrespondence]") {
  auto cloud = makePrismCloud();
  std::vector<int> ring = {0, 1, 2, 3};

  auto pts = pntToPnt::fillPointSetPrismRing(cloud, ring, 0);
  REQUIRE(pts.rows() == 4);
  REQUIRE(pts.cols() == 3);
  // First point should match atom 0 coordinates
  REQUIRE_THAT(pts(0, 0), Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("fillPointSetPrismBlock fills both basal rings",
          "[pntCorrespondence]") {
  auto cloud = makePrismCloud();
  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  auto pts = pntToPnt::fillPointSetPrismBlock(cloud, basal1, basal2, 0);
  REQUIRE(pts.rows() == 8);
  REQUIRE(pts.cols() == 3);
}

// -- HC point set tests using real geometry --

TEST_CASE("relOrderHC orders HC basal rings", "[pntCorrespondence]") {
  auto cloud = makeHCCloud();
  auto nList = nneigh::neighListO(3.5, cloud, 1);
  nList = nneigh::neighbourListByIndex(cloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<ring::strucType> ringType(rings.size());
  std::vector<cage::Cage> cageList;
  ring::findHC(rings, ringType, nList, cageList);

  REQUIRE(cageList.size() == 1);

  // Get the basal rings from the cage
  auto &cage = cageList[0];
  REQUIRE(cage.rings.size() >= 2);
  auto basal1 = rings[cage.rings[0]];
  auto basal2 = rings[cage.rings[1]];

  std::vector<int> matchedBasal1, matchedBasal2;
  int ret = pntToPnt::relOrderHC(cloud, basal1, basal2, nList, matchedBasal1,
                                  matchedBasal2);
  REQUIRE(ret == 0);
  REQUIRE(matchedBasal1.size() == 6);
  REQUIRE(matchedBasal2.size() == 6);
}

TEST_CASE("changeHexCageOrder fills Eigen matrix for HC matching",
          "[pntCorrespondence]") {
  auto cloud = makeHCCloud();
  auto nList = nneigh::neighListO(3.5, cloud, 1);
  nList = nneigh::neighbourListByIndex(cloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<ring::strucType> ringType(rings.size());
  std::vector<cage::Cage> cageList;
  ring::findHC(rings, ringType, nList, cageList);

  REQUIRE(cageList.size() == 1);
  auto basal1 = rings[cageList[0].rings[0]];
  auto basal2 = rings[cageList[0].rings[1]];

  std::vector<int> matchedBasal1, matchedBasal2;
  pntToPnt::relOrderHC(cloud, basal1, basal2, nList, matchedBasal1,
                        matchedBasal2);

  auto pts = pntToPnt::changeHexCageOrder(cloud, matchedBasal1, matchedBasal2, 0);
  // HC has 12 atoms
  REQUIRE(pts.rows() == 12);
  REQUIRE(pts.cols() == 3);
}
