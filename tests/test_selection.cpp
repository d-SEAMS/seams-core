#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <selection.hpp>
#include <seams_output.hpp>

#include <array>
#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

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

// Helper: build a cloud with rings for edge molecule tests
// 6 atoms at known positions, all same type and molecule ID
static molSys::PointCloud<molSys::Point<double>, double>
makeRingSelectionCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  // 6 atoms forming a ring, some inside slice and some outside
  // Atoms 0-3 at low coords (in slice), atoms 4-5 at high coords (out of slice)
  double coords[6][3] = {
      {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0},
      {4.0, 1.0, 1.0}, {8.0, 8.0, 8.0}, {9.0, 8.0, 8.0}};

  for (int i = 0; i < 6; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.molID = i; // each atom in its own molecule for simplicity
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    pt.inSlice = false;
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 6;
  return cloud;
}

// -- getEdgeMoleculesInRings tests --

TEST_CASE("getEdgeMoleculesInRings marks edge atoms in rings", "[selection]") {
  auto oCloud = makeRingSelectionCloud();
  auto yCloud = makeRingSelectionCloud();

  // A ring containing atoms 0, 1, 4 (atom 0 and 1 in slice, 4 outside)
  std::vector<std::vector<int>> rings = {{0, 1, 4}};

  std::array<double, 3> lo = {0.0, 0.0, 0.0};
  std::array<double, 3> hi = {5.0, 5.0, 5.0};

  // identicalCloud=true means oCloud and yCloud are the same
  ring::getEdgeMoleculesInRings(rings, oCloud, yCloud, lo, hi, true);

  // Atom 0 and 1 are in the slice, so the ring is "in slice"
  // Therefore atom 4 (outside slice but in ring) should be marked inSlice
  REQUIRE(oCloud.pts[0].inSlice == true);
  REQUIRE(oCloud.pts[1].inSlice == true);
  // Atom 4 is part of a ring that has atoms in the slice
  REQUIRE(oCloud.pts[4].inSlice == true);
}

// -- printSliceGetEdgeMoleculesInRings tests --

TEST_CASE("printSliceGetEdgeMoleculesInRings runs without crash",
          "[selection]") {
  auto oCloud = makeRingSelectionCloud();
  auto yCloud = makeRingSelectionCloud();

  std::vector<std::vector<int>> rings = {{0, 1, 2}};

  std::array<double, 3> lo = {0.0, 0.0, 0.0};
  std::array<double, 3> hi = {5.0, 5.0, 5.0};

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_printslice/").string();
  fs::create_directories(tmpPath);

  // Should not crash; exercises the full pipeline
  ring::printSliceGetEdgeMoleculesInRings(tmpPath, rings, oCloud, yCloud, lo,
                                           hi, true);

  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}
