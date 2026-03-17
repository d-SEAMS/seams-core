#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <neighbours.hpp>
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
