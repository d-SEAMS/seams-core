#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bulkTUM.hpp>
#include <cage.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

// Read the mW cubic trajectory (single atom type, ideal for bulk ice analysis)
static molSys::PointCloud<molSys::Point<double>, double>
readMwCubicFrame1() {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  return sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
}

TEST_CASE("buildRefHC reads HC template from XYZ file", "[bulkTUM]") {
  auto refPnts = tum3::buildRefHC("../templates/hc.xyz");
  // HC has 12 atoms, 3 dimensions
  REQUIRE(refPnts.rows() == 12);
  REQUIRE(refPnts.cols() == 3);
}

TEST_CASE("buildRefDDC reads DDC template from XYZ file", "[bulkTUM]") {
  auto refPnts = tum3::buildRefDDC("../templates/ddc.xyz");
  // DDC has 14 atoms, 3 dimensions
  REQUIRE(refPnts.rows() == 14);
  REQUIRE(refPnts.cols() == 3);
}

TEST_CASE("topoBulkCriteria finds cages in mW cubic ice", "[bulkTUM]") {
  auto yCloud = readMwCubicFrame1();
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  int numHC = 0, numDDC = 0;
  std::vector<ring::strucType> ringType;

  std::string tmpPath = "/tmp/dseams_test_topobulkcrit/";

  auto cageList = tum3::topoBulkCriteria(tmpPath, rings, nList, yCloud, 1,
                                          numHC, numDDC, ringType);

  // mW cubic ice should have DDC cages
  REQUIRE((numHC + numDDC) > 0);

  fs::remove_all(tmpPath);
}

TEST_CASE("topoUnitMatchingBulk runs full TUM pipeline on mW cubic",
          "[bulkTUM]") {
  auto yCloud = readMwCubicFrame1();
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = "/tmp/dseams_test_tum_full/";

  int ret = tum3::topoUnitMatchingBulk(tmpPath, rings, nList, yCloud, 1, false,
                                         true, "../templates");
  REQUIRE(ret == 0);

  fs::remove_all(tmpPath);
}

TEST_CASE("topoUnitMatchingBulk with clustering enabled on mW cubic",
          "[bulkTUM]") {
  auto yCloud = readMwCubicFrame1();
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = "/tmp/dseams_test_tum_cluster/";

  // printClusters=true exercises clusterCages
  int ret = tum3::topoUnitMatchingBulk(tmpPath, rings, nList, yCloud, 1, true,
                                         true, "../templates");
  REQUIRE(ret == 0);

  fs::remove_all(tmpPath);
}

TEST_CASE("topoUnitMatchingBulk with onlyTetrahedral=false on mW cubic",
          "[bulkTUM]") {
  auto yCloud = readMwCubicFrame1();
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = "/tmp/dseams_test_tum_nonttet/";

  // onlyTetrahedral=false exercises the 5-membered ring (PNC) path
  int ret = tum3::topoUnitMatchingBulk(tmpPath, rings, nList, yCloud, 1, false,
                                         false, "../templates");
  REQUIRE(ret == 0);

  fs::remove_all(tmpPath);
}

TEST_CASE("atomsFromCages extracts unique atom indices from cage list",
          "[bulkTUM]") {
  // Synthetic cage data
  std::vector<std::vector<int>> rings = {{0, 1, 2, 3, 4, 5},
                                          {6, 7, 8, 9, 10, 11},
                                          {0, 3, 6, 9, 1, 7}};
  cage::Cage cage1;
  cage1.type = cage::cageType::HexC;
  cage1.rings = {0, 1, 2};
  std::vector<cage::Cage> cageList = {cage1};
  std::vector<int> clusterCageIndices = {0};

  auto atoms = tum3::atomsFromCages(rings, cageList, clusterCageIndices);

  // Should have unique atoms from all 3 rings
  REQUIRE(atoms.size() == 12); // 0-11, all unique
}
