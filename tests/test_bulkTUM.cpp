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

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_topobulkcrit/").string();

  auto cageList = tum3::topoBulkCriteria(tmpPath, rings, nList, yCloud, 1,
                                          numHC, numDDC, ringType);

  // mW cubic ice should have DDC cages
  REQUIRE((numHC + numDDC) > 0);

  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}

TEST_CASE("topoUnitMatchingBulk runs full TUM pipeline on mW cubic",
          "[bulkTUM]") {
  auto yCloud = readMwCubicFrame1();
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_tum_full/").string();

  int ret = tum3::topoUnitMatchingBulk(tmpPath, rings, nList, yCloud, 1, false,
                                         true, "../templates");
  REQUIRE(ret == 0);

  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}

TEST_CASE("topoUnitMatchingBulk with clustering enabled on mW cubic",
          "[bulkTUM]") {
  auto yCloud = readMwCubicFrame1();
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_tum_cluster/").string();

  // printClusters=true exercises clusterCages
  int ret = tum3::topoUnitMatchingBulk(tmpPath, rings, nList, yCloud, 1, true,
                                         true, "../templates");
  REQUIRE(ret == 0);

  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}

TEST_CASE("topoUnitMatchingBulk with onlyTetrahedral=false on mW cubic",
          "[bulkTUM]") {
  auto yCloud = readMwCubicFrame1();
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_tum_nonttet/").string();

  // onlyTetrahedral=false exercises the 5-membered ring (PNC) path
  int ret = tum3::topoUnitMatchingBulk(tmpPath, rings, nList, yCloud, 1, false,
                                         false, "../templates");
  REQUIRE(ret == 0);

  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}

TEST_CASE("shapeMatchHC computes RMSD for a known HC cage", "[bulkTUM]") {
  // Build the 12-atom HC geometry
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud.box = {53.9690018, 54.5289994, 51.257};
  yCloud.boxLow = {0.0, 0.0, 0.0};
  yCloud.currentFrame = 1;
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
    yCloud.pts.push_back(pt);
    yCloud.idIndexMap[i] = i;
  }
  yCloud.nop = 12;

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  // Find the HC cage
  std::vector<ring::strucType> ringType(rings.size());
  std::vector<cage::Cage> cageList;
  ring::findHC(rings, ringType, nList, cageList);
  REQUIRE(cageList.size() == 1);

  // Build reference HC
  auto refPnts = tum3::buildRefHC("../templates/hc.xyz");
  REQUIRE(refPnts.rows() == 12);

  // Run shapeMatchHC
  std::vector<double> quat;
  double rmsd = -1.0;
  int ret = tum3::shapeMatchHC(yCloud, refPnts, cageList[0], rings, nList,
                                quat, rmsd);
  REQUIRE(ret == 0);
  // RMSD should be non-negative
  REQUIRE(rmsd >= 0.0);
  // Quaternion should have 4 elements
  REQUIRE(quat.size() == 4);
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
