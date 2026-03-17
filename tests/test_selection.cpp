#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <selection.hpp>

#include <array>
#include <vector>

// Helper: build a cloud with atoms at known positions and molecule IDs
static molSys::PointCloud<molSys::Point<double>, double>
makeSelectionCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  // 6 atoms: types 1 and 2, grouped into 2 molecules
  // Molecule 0: atoms 0,1,2 (type 1, 2, 1)
  // Molecule 1: atoms 3,4,5 (type 2, 1, 2)
  int types[] = {1, 2, 1, 2, 1, 2};
  int molIDs[] = {0, 0, 0, 1, 1, 1};
  double coords[6][3] = {
      {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0},
      {7.0, 7.0, 7.0}, {8.0, 7.0, 7.0}, {9.0, 7.0, 7.0}};

  for (int i = 0; i < 6; i++) {
    molSys::Point<double> pt;
    pt.type = types[i];
    pt.atomID = i;
    pt.molID = molIDs[i];
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    pt.inSlice = true;
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 6;
  return cloud;
}

// -- getPointCloudOneAtomType tests --

TEST_CASE("getPointCloudOneAtomType filters by type", "[selection]") {
  auto cloud = makeSelectionCloud();
  molSys::PointCloud<molSys::Point<double>, double> outCloud;

  auto result = gen::getPointCloudOneAtomType(
      cloud, outCloud, 1, false, {0, 0, 0}, {0, 0, 0});

  // Type 1 atoms: indices 0, 2, 4
  REQUIRE(result.nop == 3);
  for (int i = 0; i < result.nop; i++) {
    REQUIRE(result.pts[i].type == 1);
  }
}

TEST_CASE("getPointCloudOneAtomType with slice filtering", "[selection]") {
  auto cloud = makeSelectionCloud();
  molSys::PointCloud<molSys::Point<double>, double> outCloud;

  // Slice: x in [0,5], y in [0,5], z in [0,5]
  std::array<double, 3> lo = {0.0, 0.0, 0.0};
  std::array<double, 3> hi = {5.0, 5.0, 5.0};

  auto result =
      gen::getPointCloudOneAtomType(cloud, outCloud, 1, true, lo, hi);

  // Only type 1 atoms in the slice: atoms 0 (1,1,1) and 2 (3,1,1)
  REQUIRE(result.nop == 2);
}

TEST_CASE("getPointCloudOneAtomType preserves box dimensions", "[selection]") {
  auto cloud = makeSelectionCloud();
  molSys::PointCloud<molSys::Point<double>, double> outCloud;

  auto result = gen::getPointCloudOneAtomType(
      cloud, outCloud, 2, false, {0, 0, 0}, {0, 0, 0});

  REQUIRE(result.box.size() == 3);
  REQUIRE_THAT(result.box[0], Catch::Matchers::WithinAbs(10.0, 1e-10));
}

// -- atomsInSingleSlice tests --

TEST_CASE("atomsInSingleSlice marks atoms inside slice", "[selection]") {
  auto cloud = makeSelectionCloud();

  // Slice: x in [0,5], y in [0,5], z in [0,5]
  std::array<double, 3> lo = {0.0, 0.0, 0.0};
  std::array<double, 3> hi = {5.0, 5.0, 5.0};

  gen::atomsInSingleSlice(cloud, true, lo, hi);

  // Atoms 0,1,2 are at (1,1,1), (2,1,1), (3,1,1) -- all in slice
  REQUIRE(cloud.pts[0].inSlice == true);
  REQUIRE(cloud.pts[1].inSlice == true);
  REQUIRE(cloud.pts[2].inSlice == true);

  // Atoms 3,4,5 are at (7,7,7), (8,7,7), (9,7,7) -- all outside slice
  REQUIRE(cloud.pts[3].inSlice == false);
  REQUIRE(cloud.pts[4].inSlice == false);
  REQUIRE(cloud.pts[5].inSlice == false);
}

TEST_CASE("atomsInSingleSlice clearPreviousSliceSelection=false preserves",
          "[selection]") {
  auto cloud = makeSelectionCloud();

  // First, mark only the upper region
  std::array<double, 3> lo = {6.0, 6.0, 6.0};
  std::array<double, 3> hi = {10.0, 10.0, 10.0};
  gen::atomsInSingleSlice(cloud, true, lo, hi);

  REQUIRE(cloud.pts[0].inSlice == false);
  REQUIRE(cloud.pts[3].inSlice == true);

  // Now add lower region without clearing
  std::array<double, 3> lo2 = {0.0, 0.0, 0.0};
  std::array<double, 3> hi2 = {5.0, 5.0, 5.0};
  gen::atomsInSingleSlice(cloud, false, lo2, hi2);

  // Atom 0 should now be in slice (from second call)
  REQUIRE(cloud.pts[0].inSlice == true);
}

// -- moleculesInSingleSlice tests --

TEST_CASE("moleculesInSingleSlice includes whole molecule", "[selection]") {
  auto cloud = makeSelectionCloud();

  // Slice that includes atom 0 (1,1,1) but not atom 1 (2,1,1):
  // use a tight slice
  std::array<double, 3> lo = {0.5, 0.5, 0.5};
  std::array<double, 3> hi = {1.5, 1.5, 1.5};

  gen::moleculesInSingleSlice(cloud, true, lo, hi);

  // Atom 0 is in slice, so all atoms of molecule 0 (atoms 0,1,2) should be
  // inSlice
  REQUIRE(cloud.pts[0].inSlice == true);
  REQUIRE(cloud.pts[1].inSlice == true);
  REQUIRE(cloud.pts[2].inSlice == true);

  // Molecule 1 (atoms 3,4,5) should not be in slice
  REQUIRE(cloud.pts[3].inSlice == false);
  REQUIRE(cloud.pts[4].inSlice == false);
  REQUIRE(cloud.pts[5].inSlice == false);
}

// -- setAtomsWithSameMolID tests --

TEST_CASE("setAtomsWithSameMolID sets all atoms of a molecule", "[selection]") {
  auto cloud = makeSelectionCloud();

  // Set all atoms to not in slice first
  for (auto &pt : cloud.pts) {
    pt.inSlice = false;
  }

  auto molMap = molSys::createMolIDAtomIDMultiMap(cloud);
  gen::setAtomsWithSameMolID(cloud, molMap, 0, true);

  // All atoms of molecule 0 (atoms 0,1,2) should be true
  REQUIRE(cloud.pts[0].inSlice == true);
  REQUIRE(cloud.pts[1].inSlice == true);
  REQUIRE(cloud.pts[2].inSlice == true);

  // Molecule 1 should be unaffected
  REQUIRE(cloud.pts[3].inSlice == false);
}
