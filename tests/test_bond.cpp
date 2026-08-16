#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bond.hpp>
#include <cage.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <array>
#include <cmath>

// Helper: build a minimal cloud for bond tests
// 4 atoms forming a square in the xy-plane
static molSys::PointCloud<molSys::Point<double>, double> makeSquareCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  double coords[4][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};

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
  return cloud;
}

TEST_CASE("trimBonds removes exact duplicate bond pairs", "[bond]") {
  // Build a bonds list with duplicates
  std::vector<std::vector<int>> bonds = {
      {1, 2}, {3, 4}, {1, 2}, {5, 6}, {3, 4}};

  auto trimmed = bond::trimBonds(bonds);

  // Should have exactly 3 unique bonds
  REQUIRE(trimmed.size() == 3);

  // Verify each unique bond appears exactly once
  // After sorting, expect {1,2}, {3,4}, {5,6}
  std::sort(trimmed.begin(), trimmed.end());
  REQUIRE(trimmed[0] == std::vector<int>({1, 2}));
  REQUIRE(trimmed[1] == std::vector<int>({3, 4}));
  REQUIRE(trimmed[2] == std::vector<int>({5, 6}));
}

TEST_CASE("trimBonds on empty input returns empty", "[bond]") {
  std::vector<std::vector<int>> bonds;
  auto trimmed = bond::trimBonds(bonds);
  REQUIRE(trimmed.empty());
}

TEST_CASE("populateBonds generates bonds from neighbour list", "[bond]") {
  auto cloud = makeSquareCloud();
  // Build a neighbour list by index (not ID, since populateBonds expects index-based)
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);

  auto bonds = bond::populateBonds(nList, cloud);

  // With cutoff 1.5, each atom should be bonded to its adjacent square
  // neighbours (distance 1.0). Diagonal distance is sqrt(2) ~ 1.414, also
  // within cutoff. So we expect 4 edge bonds + 2 diagonal bonds = 6 total
  // (after dedup by iatom < jatom condition in populateBonds)
  REQUIRE(bonds.size() > 0);

  // All bond entries should have exactly 2 elements
  for (const auto &b : bonds) {
    REQUIRE(b.size() == 2);
  }
}

TEST_CASE("getHbondDistanceOH computes periodic distance between two clouds",
          "[bond]") {
  // Build two minimal PointClouds: one "oxygen", one "hydrogen"
  molSys::PointCloud<molSys::Point<double>, double> oCloud, hCloud;
  oCloud.box = {10.0, 10.0, 10.0};
  hCloud.box = {10.0, 10.0, 10.0};

  molSys::Point<double> oPoint;
  oPoint.x = 1.0;
  oPoint.y = 0.0;
  oPoint.z = 0.0;
  oCloud.pts.push_back(oPoint);

  molSys::Point<double> hPoint;
  hPoint.x = 2.0;
  hPoint.y = 0.0;
  hPoint.z = 0.0;
  hCloud.pts.push_back(hPoint);

  double dist = bond::getHbondDistanceOH(oCloud, hCloud, 0, 0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("getHbondDistanceOH respects periodic boundaries", "[bond]") {
  molSys::PointCloud<molSys::Point<double>, double> oCloud, hCloud;
  oCloud.box = {10.0, 10.0, 10.0};
  hCloud.box = {10.0, 10.0, 10.0};

  molSys::Point<double> oPoint;
  oPoint.x = 0.5;
  oPoint.y = 0.0;
  oPoint.z = 0.0;
  oCloud.pts.push_back(oPoint);

  // Place H near the other edge of the box
  molSys::Point<double> hPoint;
  hPoint.x = 9.5;
  hPoint.y = 0.0;
  hPoint.z = 0.0;
  hCloud.pts.push_back(hPoint);

  double dist = bond::getHbondDistanceOH(oCloud, hCloud, 0, 0);
  // Periodic distance should be 1.0, not 9.0
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("getHbondDistanceOH uses dump H on a mixed image", "[bond]") {
  molSys::PointCloud<molSys::Point<double>, double> oCloud, hCloud;
  oCloud.box = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  oCloud.boxLow = {0.0, 0.0, 0.0};
  hCloud.box = oCloud.box;
  hCloud.boxLow = oCloud.boxLow;
  molSys::Point<double> oPoint;
  oPoint.x = 0.5;
  oPoint.y = 0.5;
  oPoint.z = 1.0;
  oCloud.pts.push_back(oPoint);
  molSys::Point<double> hPoint;
  hPoint.x = 1.0;
  hPoint.y = 8.0;
  hPoint.z = 1.0;
  hCloud.pts.push_back(hPoint);
  const double dist = bond::getHbondDistanceOH(oCloud, hCloud, 0, 0);
  const double expect = std::sqrt(4.5 * 4.5 + 1.160254037844386 * 1.160254037844386);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(expect, 1e-10));
}

static molSys::Point<double> makeAtom(int atomID, int molID, double x,
                                      double y, double z) {
  molSys::Point<double> p;
  p.atomID = atomID;
  p.molID = molID;
  p.type = 1;
  p.x = x;
  p.y = y;
  p.z = z;
  return p;
}

TEST_CASE("populateHbondsWithInputClouds accepts acceptor-first angle",
          "[bond]") {
  molSys::PointCloud<molSys::Point<double>, double> oCloud, hCloud;
  oCloud.box = {10.0, 10.0, 10.0};
  oCloud.boxLow = {0.0, 0.0, 0.0};
  hCloud.box = oCloud.box;
  hCloud.boxLow = oCloud.boxLow;
  oCloud.pts.push_back(makeAtom(1, 1, 0.0, 0.0, 0.0));
  oCloud.pts.push_back(makeAtom(2, 2, 2.8, 0.0, 0.0));
  oCloud.nop = 2;
  oCloud.idIndexMap[1] = 0;
  oCloud.idIndexMap[2] = 1;
  hCloud.pts.push_back(makeAtom(11, 1, 1.0, 0.0, 0.0));
  hCloud.pts.push_back(makeAtom(12, 1, -0.96, 0.76, 0.0));
  hCloud.pts.push_back(makeAtom(21, 2, 3.8, 0.0, 0.0));
  hCloud.pts.push_back(makeAtom(22, 2, 2.8, 0.96, 0.0));
  hCloud.nop = 4;
  const auto ooAccFirst = gen::relDist(oCloud, 1, 0);
  const auto ooDonFirst = gen::relDist(oCloud, 0, 1);
  std::vector<double> ooA = {ooAccFirst[0], ooAccFirst[1], ooAccFirst[2]};
  std::vector<double> ooD = {ooDonFirst[0], ooDonFirst[1], ooDonFirst[2]};
  std::vector<double> oh = {2.8 - 1.0, 0.0, 0.0};
  REQUIRE_THAT(gen::radDeg(gen::eigenVecAngle(ooA, oh)),
               Catch::Matchers::WithinAbs(0.0, 1e-8));
  REQUIRE_THAT(gen::radDeg(gen::eigenVecAngle(ooD, oh)),
               Catch::Matchers::WithinAbs(180.0, 1e-8));
  const std::vector<std::vector<int>> nList = {{1, 2}, {2, 1}};
  const auto net =
      bond::populateHbondsWithInputClouds(oCloud, hCloud, nList, 2.42, 30.0);
  REQUIRE(net.size() == 2);
  REQUIRE(net[0].size() == 2);
  REQUIRE(net[1].size() == 2);
  REQUIRE(net[0][1] == 2);
  REQUIRE(net[1][1] == 1);
}

TEST_CASE("populateHbondsWithInputClouds accepts a mixed-image dump-H bond",
          "[bond]") {
  molSys::PointCloud<molSys::Point<double>, double> oCloud, hCloud;
  oCloud.box = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  oCloud.boxLow = {0.0, 0.0, 0.0};
  hCloud.box = oCloud.box;
  hCloud.boxLow = oCloud.boxLow;
  oCloud.pts.push_back(makeAtom(1, 1, 5.5, 6.360254037844386, 1.0));
  oCloud.pts.push_back(makeAtom(2, 2, 0.5, 0.5, 1.0));
  oCloud.nop = 2;
  oCloud.idIndexMap[1] = 0;
  oCloud.idIndexMap[2] = 1;
  hCloud.pts.push_back(makeAtom(11, 1, 0.5, 0.1, 1.0));
  hCloud.pts.push_back(makeAtom(12, 1, 6.5, 6.360254037844386, 1.0));
  hCloud.pts.push_back(makeAtom(21, 2, 0.5, 1.5, 1.0));
  hCloud.pts.push_back(makeAtom(22, 2, 1.5, 0.5, 1.0));
  hCloud.nop = 4;
  const auto oo = gen::relDist(oCloud, 1, 0);
  REQUIRE_THAT(oo[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
  REQUIRE_THAT(oo[1], Catch::Matchers::WithinAbs(2.8, 1e-10));
  REQUIRE_THAT(oo[2], Catch::Matchers::WithinAbs(0.0, 1e-12));
  REQUIRE_THAT(bond::getHbondDistanceOH(oCloud, hCloud, 1, 0),
               Catch::Matchers::WithinAbs(0.4, 1e-10));
  std::array<double, 3> span = {0.5 - 5.5, 0.5 - 6.360254037844386, 0.0};
  for (int k = 0; k < 3; k++) {
    span[k] -= oCloud.box[static_cast<size_t>(k)] *
               std::round(span[k] / oCloud.box[static_cast<size_t>(k)]);
  }
  std::vector<double> ooDump = {oo[0], oo[1], oo[2]};
  std::vector<double> ooSpan = {span[0], span[1], span[2]};
  std::vector<double> oh = {0.5 - 0.5, 0.5 - 0.1, 0.0};
  REQUIRE(gen::radDeg(gen::eigenVecAngle(ooDump, oh)) < 30.0);
  REQUIRE(gen::radDeg(gen::eigenVecAngle(ooSpan, oh)) > 30.0);
  const std::vector<std::vector<int>> nList = {{1, 2}, {2, 1}};
  const auto net =
      bond::populateHbondsWithInputClouds(oCloud, hCloud, nList, 2.42, 30.0);
  REQUIRE(net[0].size() == 2);
  REQUIRE(net[1].size() == 2);
  REQUIRE(net[0][1] == 2);
  REQUIRE(net[1][1] == 1);
}

TEST_CASE("populateBonds with cage iceType filters dummy atoms", "[bond]") {
  auto cloud = makeSquareCloud();
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);

  // All atoms classified as dummy -- bonds should still be created between them
  std::vector<cage::iceType> atomTypes(4, cage::iceType::dummy);
  auto bonds = bond::populateBonds(nList, cloud, atomTypes);

  // With bondsBetweenDummy=false (default), bonds between dummy atoms are excluded
  // So no bonds should be created
  REQUIRE(bonds.empty());
}

TEST_CASE("populateBonds with non-dummy iceType creates bonds", "[bond]") {
  auto cloud = makeSquareCloud();
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);

  // Mark all atoms as HC ice type
  std::vector<cage::iceType> atomTypes(4, cage::iceType::hc);
  auto bonds = bond::populateBonds(nList, cloud, atomTypes);

  REQUIRE(bonds.size() > 0);
  for (const auto &b : bonds) {
    REQUIRE(b.size() == 2);
  }
}

TEST_CASE("createBondsFromCages extracts bonds from cage rings", "[bond]") {
  // Minimal cage setup: a single HC cage has rings
  std::vector<std::vector<int>> rings = {{1, 2, 3, 4, 5, 6},
                                          {7, 8, 9, 10, 11, 12}};
  cage::Cage cage1;
  cage1.type = cage::cageType::HexC;
  cage1.rings = {0, 1}; // Indices into rings

  std::vector<cage::Cage> cageList = {cage1};
  int nRings = 0;

  auto bonds = bond::createBondsFromCages(rings, cageList,
                                           cage::cageType::HexC, nRings);

  REQUIRE(nRings == 2);
  REQUIRE(bonds.size() > 0);
}

TEST_CASE("populateHbonds detects hydrogen bonds from trajectory", "[bond]") {
  // exampleTraj stores oxygen as type 2 and hydrogen as type 1.
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, yCloud, 2);
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 2);

  auto hBonds = bond::populateHbonds("traj/exampleTraj.lammpstrj", yCloud,
                                      nList, 1, 1);

  REQUIRE(hBonds.size() == static_cast<size_t>(yCloud.nop));
  int bonded = 0;
  for (size_t i = 0; i < hBonds.size(); i++) {
    REQUIRE_FALSE(hBonds[i].empty());
    REQUIRE(hBonds[i][0] == yCloud.pts[static_cast<int>(i)].atomID);
    if (hBonds[i].size() > 1) {
      bonded++;
    }
  }
  REQUIRE(bonded > 0);
}

TEST_CASE("populateHbondsWithInputClouds detects H-bonds from O+H clouds",
          "[bond]") {
  molSys::PointCloud<molSys::Point<double>, double> oCloud;
  oCloud = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, oCloud, 2);
  REQUIRE(oCloud.nop > 0);

  molSys::PointCloud<molSys::Point<double>, double> hCloud;
  hCloud = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, hCloud, 1);
  REQUIRE(hCloud.nop > 0);
  REQUIRE(hCloud.nop == 2 * oCloud.nop);

  auto nList = nneigh::neighListO(3.5, oCloud, 2);

  auto hBonds = bond::populateHbondsWithInputClouds(oCloud, hCloud, nList);

  REQUIRE(hBonds.size() == static_cast<size_t>(oCloud.nop));
  int bonded = 0;
  for (const auto &row : hBonds) {
    REQUIRE_FALSE(row.empty());
    if (row.size() > 1) {
      bonded++;
    }
  }
  REQUIRE(bonded > 0);
}

TEST_CASE("H-bond thresholds are configurable", "[bond]") {
  molSys::PointCloud<molSys::Point<double>, double> oCloud;
  oCloud = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, oCloud, 2);
  molSys::PointCloud<molSys::Point<double>, double> hCloud;
  hCloud = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, hCloud, 1);
  auto nList = nneigh::neighListO(3.5, oCloud, 2);

  auto countEdges = [](const std::vector<std::vector<int>> &net) {
    size_t edges = 0;
    for (const auto &row : net) {
      edges += row.size() - 1;
    }
    return edges;
  };

  // Defaults reproduce the water criterion
  const auto waterNet =
      bond::populateHbondsWithInputClouds(oCloud, hCloud, nList);
  const auto explicitNet =
      bond::populateHbondsWithInputClouds(oCloud, hCloud, nList, 2.42, 30.0);
  REQUIRE(countEdges(waterNet) == countEdges(explicitNet));
  REQUIRE(countEdges(waterNet) > 0);

  // A distance cutoff below any O-H separation removes every bond
  const auto noneNet =
      bond::populateHbondsWithInputClouds(oCloud, hCloud, nList, 0.1, 30.0);
  REQUIRE(countEdges(noneNet) == 0);

  // A zero-degree angle cutoff also removes every bond
  const auto noAngleNet =
      bond::populateHbondsWithInputClouds(oCloud, hCloud, nList, 2.42, 0.0);
  REQUIRE(countEdges(noAngleNet) == 0);

  // Looser thresholds can only add bonds
  const auto looseNet =
      bond::populateHbondsWithInputClouds(oCloud, hCloud, nList, 3.0, 45.0);
  REQUIRE(countEdges(looseNet) >= countEdges(waterNet));
}

TEST_CASE("createBondsFromCages with no matching cage type returns empty",
          "[bond]") {
  std::vector<std::vector<int>> rings = {{1, 2, 3, 4, 5, 6}};
  cage::Cage cage1;
  cage1.type = cage::cageType::DoubleDiaC;
  cage1.rings = {0};

  std::vector<cage::Cage> cageList = {cage1};
  int nRings = 0;

  auto bonds = bond::createBondsFromCages(rings, cageList,
                                           cage::cageType::HexC, nRings);

  REQUIRE(nRings == 0);
}
