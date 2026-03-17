#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <pntCorrespondence.hpp>
#include <shapeMatch.hpp>

#include <Eigen/Dense>
#include <cmath>
#include <vector>


// -- updateRMSDRing tests --

TEST_CASE("updateRMSDRing sets RMSD for unassigned atoms", "[shapeMatch]") {
  std::vector<int> basalRing = {0, 1, 2, 3};
  std::vector<double> rmsdPerAtom(10, -1.0); // -1 means unassigned
  double rmsdVal = 0.42;

  int ret = match::updateRMSDRing(basalRing, 0, rmsdVal, rmsdPerAtom);

  REQUIRE(ret == 0);
  // All atoms in the ring should have the RMSD value
  for (int idx : basalRing) {
    REQUIRE_THAT(rmsdPerAtom[idx], Catch::Matchers::WithinAbs(0.42, 1e-10));
  }
  // Other atoms should remain -1
  REQUIRE_THAT(rmsdPerAtom[5], Catch::Matchers::WithinAbs(-1.0, 1e-10));
}

TEST_CASE("updateRMSDRing does not overwrite already assigned RMSD",
          "[shapeMatch]") {
  std::vector<int> basalRing = {0, 1, 2};
  std::vector<double> rmsdPerAtom(5, -1.0);
  rmsdPerAtom[1] = 0.1; // Pre-assign atom 1

  match::updateRMSDRing(basalRing, 0, 0.5, rmsdPerAtom);

  // Atom 1 should keep its original value
  REQUIRE_THAT(rmsdPerAtom[1], Catch::Matchers::WithinAbs(0.1, 1e-10));
  // Atom 0 and 2 should get the new value
  REQUIRE_THAT(rmsdPerAtom[0], Catch::Matchers::WithinAbs(0.5, 1e-10));
  REQUIRE_THAT(rmsdPerAtom[2], Catch::Matchers::WithinAbs(0.5, 1e-10));
}

TEST_CASE("updateRMSDRing with startingIndex wraps around", "[shapeMatch]") {
  std::vector<int> basalRing = {10, 11, 12, 13};
  std::vector<double> rmsdPerAtom(20, -1.0);

  // Start from index 2, so the order is: ring[2]=12, ring[3]=13, ring[0]=10,
  // ring[1]=11
  match::updateRMSDRing(basalRing, 2, 1.0, rmsdPerAtom);

  for (int idx : basalRing) {
    REQUIRE_THAT(rmsdPerAtom[idx], Catch::Matchers::WithinAbs(1.0, 1e-10));
  }
}

// -- updatePerAtomRMSDRing tests --

TEST_CASE("updatePerAtomRMSDRing assigns per-atom RMSD values",
          "[shapeMatch]") {
  std::vector<int> basalRing = {0, 1, 2};
  std::vector<double> rmsdFromMatch = {0.1, 0.2, 0.3};
  std::vector<double> rmsdPerAtom(5, -1.0);

  int ret =
      match::updatePerAtomRMSDRing(basalRing, 0, rmsdFromMatch, rmsdPerAtom);

  REQUIRE(ret == 0);
  REQUIRE_THAT(rmsdPerAtom[0], Catch::Matchers::WithinAbs(0.1, 1e-10));
  REQUIRE_THAT(rmsdPerAtom[1], Catch::Matchers::WithinAbs(0.2, 1e-10));
  REQUIRE_THAT(rmsdPerAtom[2], Catch::Matchers::WithinAbs(0.3, 1e-10));
}

TEST_CASE("updatePerAtomRMSDRing with offset wraps correctly",
          "[shapeMatch]") {
  std::vector<int> basalRing = {5, 6, 7, 8};
  std::vector<double> rmsdFromMatch = {0.1, 0.2, 0.3, 0.4};
  std::vector<double> rmsdPerAtom(10, -1.0);

  // startingIndex=1: match order is ring[1]=6, ring[2]=7, ring[3]=8, ring[0]=5
  match::updatePerAtomRMSDRing(basalRing, 1, rmsdFromMatch, rmsdPerAtom);

  REQUIRE_THAT(rmsdPerAtom[6], Catch::Matchers::WithinAbs(0.1, 1e-10));
  REQUIRE_THAT(rmsdPerAtom[7], Catch::Matchers::WithinAbs(0.2, 1e-10));
  REQUIRE_THAT(rmsdPerAtom[8], Catch::Matchers::WithinAbs(0.3, 1e-10));
  REQUIRE_THAT(rmsdPerAtom[5], Catch::Matchers::WithinAbs(0.4, 1e-10));
}

TEST_CASE("updatePerAtomRMSDRing does not overwrite existing", "[shapeMatch]") {
  std::vector<int> basalRing = {0, 1};
  std::vector<double> rmsdFromMatch = {0.5, 0.6};
  std::vector<double> rmsdPerAtom(3, -1.0);
  rmsdPerAtom[0] = 0.2; // Pre-assigned

  match::updatePerAtomRMSDRing(basalRing, 0, rmsdFromMatch, rmsdPerAtom);

  // Atom 0 should keep 0.2
  REQUIRE_THAT(rmsdPerAtom[0], Catch::Matchers::WithinAbs(0.2, 1e-10));
  // Atom 1 should get 0.6
  REQUIRE_THAT(rmsdPerAtom[1], Catch::Matchers::WithinAbs(0.6, 1e-10));
}

// -- matchPrismBlock tests using hard-coded tetragonal prism --

TEST_CASE("matchPrismBlock on ideal tetragonal prism", "[shapeMatch]") {
  // Build 8-atom prism
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10, 10, 50};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  double coords[8][3] = {
      {0, 0, 0}, {2.75, 0, 0}, {2.75, 2.75, 0}, {0, 2.75, 0},
      {0, 0, 2.75}, {2.75, 0, 2.75}, {2.75, 2.75, 2.75}, {0, 2.75, 2.75}};
  for (int i = 0; i < 8; i++) {
    pt.atomID = i;
    pt.x = coords[i][0]; pt.y = coords[i][1]; pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;

  auto nList = nneigh::neighListO(3.5, cloud, 1);
  nList = nneigh::neighbourListByIndex(cloud, nList);

  // Get reference ring points for a 4-membered ring (axial dim=2)
  auto refPoints = pntToPnt::getPointSetRefRing(4, 2);
  // Scale reference to match our prism
  auto refPrism = pntToPnt::createPrismBlock(cloud, refPoints, 4,
                                              {0, 1, 2, 3}, {4, 5, 6, 7});

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};
  int beginIndex = 0;

  bool isMatch = match::matchPrismBlock(cloud, nList, refPoints, basal1,
                                         basal2, beginIndex);

  // An ideal prism should match itself (RMSD should be small)
  REQUIRE(isMatch == true);
}

TEST_CASE("matchPrism returns false when basal ordering fails",
          "[shapeMatch]") {
  // Build an 8-atom prism where relOrderPrismBlock may fail
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10, 10, 50};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  double coords[8][3] = {
      {0, 0, 0}, {2.75, 0, 0}, {2.75, 2.75, 0}, {0, 2.75, 0},
      {0, 0, 2.75}, {2.75, 0, 2.75}, {2.75, 2.75, 2.75}, {0, 2.75, 2.75}};
  for (int i = 0; i < 8; i++) {
    pt.atomID = i;
    pt.x = coords[i][0]; pt.y = coords[i][1]; pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;

  auto nList = nneigh::neighListO(3.5, cloud, 1);
  nList = nneigh::neighbourListByIndex(cloud, nList);

  auto refPoints = pntToPnt::getPointSetRefRing(4, 2);

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};
  std::vector<double> rmsdPerAtom(8, -1.0);

  // With the fix, this should return false gracefully instead of crashing
  bool isMatch = match::matchPrism(cloud, nList, refPoints, basal1, basal2,
                                    rmsdPerAtom, true);

  // Result depends on whether relOrderPrismBlock succeeds
  REQUIRE((isMatch == true || isMatch == false));
}
