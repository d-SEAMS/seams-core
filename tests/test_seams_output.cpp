#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bond.hpp>
#include <cage.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_output.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// Reusable test cloud builder
static molSys::PointCloud<molSys::Point<double>, double>
makeTestCloud(int nop, double boxLen = 10.0) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {boxLen, boxLen, boxLen};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  for (int i = 0; i < nop; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i + 1; // 1-based IDs
    pt.molID = i + 1;
    pt.x = 1.0 * (i % 4);
    pt.y = 1.0 * ((i / 4) % 4);
    pt.z = 1.0 * (i / 16);
    pt.inSlice = (i < nop / 2); // half in slice
    cloud.pts.push_back(pt);
    cloud.idIndexMap[pt.atomID] = i;
  }
  cloud.nop = nop;
  return cloud;
}

// -- makePath tests --

TEST_CASE("makePath creates directories", "[seams_output]") {
  std::string tmpDir = "/tmp/dseams_test_makepath/nested/dir";
  int ret = sout::makePath(tmpDir);
  REQUIRE(ret == 0);
  REQUIRE(fs::is_directory(tmpDir));

  // Calling again on existing dir should also return 0
  ret = sout::makePath(tmpDir);
  REQUIRE(ret == 0);

  fs::remove_all("/tmp/dseams_test_makepath");
}

// -- printRDF tests --

TEST_CASE("printRDF writes RDF data to file", "[seams_output]") {
  std::string tmpFile = "/tmp/dseams_test_rdf.dat";
  std::vector<double> rdfValues = {0.0, 0.5, 1.0, 0.8, 0.3};
  double binwidth = 0.5;
  int nbin = 5;

  int ret = sout::printRDF(tmpFile, rdfValues, binwidth, nbin);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpFile));
  REQUIRE(fs::file_size(tmpFile) > 0);

  fs::remove(tmpFile);
}

// -- writeRings tests --

TEST_CASE("writeRings writes ring data", "[seams_output]") {
  // This writes to ../output/ relative to CWD
  std::string outDir = "../output";
  fs::create_directories(outDir);

  std::vector<std::vector<int>> rings = {{1, 2, 3}, {4, 5, 6}};
  int ret = sout::writeRings(rings);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(outDir + "/rings.dat"));

  fs::remove_all(outDir);
}

// -- writeClusterStats tests --

TEST_CASE("writeClusterStats writes cluster statistics", "[seams_output]") {
  std::string tmpPath = "/tmp/dseams_test_cluster_stats/";

  int ret = sout::writeClusterStats(tmpPath, 1, 100, 5, 10, 25.0, 1);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "clusterStats.dat"));
  REQUIRE(fs::file_size(tmpPath + "clusterStats.dat") > 0);

  fs::remove_all(tmpPath);
}

TEST_CASE("writeClusterStats appends for subsequent frames", "[seams_output]") {
  std::string tmpPath = "/tmp/dseams_test_cluster_append/";

  sout::writeClusterStats(tmpPath, 1, 100, 5, 10, 25.0, 1);
  sout::writeClusterStats(tmpPath, 2, 120, 4, 15, 30.0, 1);

  // File should have header + 2 data lines
  std::ifstream f(tmpPath + "clusterStats.dat");
  int lineCount = 0;
  std::string line;
  while (std::getline(f, line)) lineCount++;
  REQUIRE(lineCount == 3); // header + 2 frames

  fs::remove_all(tmpPath);
}

// -- writeCluster tests --

TEST_CASE("writeCluster writes ice cluster info", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  std::string tmpFile = "/tmp/dseams_test_writeCluster.txt";

  int ret = sout::writeCluster(cloud, tmpFile, false, 42);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpFile));

  fs::remove(tmpFile);
}

// -- writeDump tests --

TEST_CASE("writeDump writes LAMMPS dump format", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  std::string tmpPath = "/tmp/dseams_test_dump/";

  int ret = sout::writeDump(cloud, tmpPath, "test.lammpstrj");
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "test.lammpstrj"));
  REQUIRE(fs::file_size(tmpPath + "test.lammpstrj") > 0);

  fs::remove_all(tmpPath);
}

// -- writeTopoBulkData tests --

TEST_CASE("writeTopoBulkData writes bulk topology data", "[seams_output]") {
  std::string tmpPath = "/tmp/dseams_test_topobulk/";

  int ret = sout::writeTopoBulkData(tmpPath, 1, 10, 5, 3, 8, 12, 1);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "bulkTopo/cageData.dat"));

  fs::remove_all(tmpPath);
}

// -- writePrismNum tests --

TEST_CASE("writePrismNum writes prism statistics", "[seams_output]") {
  std::string tmpPath = "/tmp/dseams_test_prismnum/";
  std::vector<int> nPrisms = {5, 3, 1};
  std::vector<int> nDefPrisms = {2, 1, 0};
  std::vector<double> heightPercent = {50.0, 30.0, 10.0};

  int ret = sout::writePrismNum(tmpPath, nPrisms, nDefPrisms, heightPercent, 5,
                                 1, 1);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "topoINT/nPrisms.dat"));

  fs::remove_all(tmpPath);
}

// -- writeRingNum tests --

TEST_CASE("writeRingNum writes ring coverage data for monolayer",
          "[seams_output]") {
  std::string tmpPath = "/tmp/dseams_test_ringnum/";
  std::vector<int> nRings = {10, 5, 2};
  std::vector<double> covXY = {0.5, 0.3, 0.1};
  std::vector<double> covXZ = {0.4, 0.2, 0.1};
  std::vector<double> covYZ = {0.3, 0.1, 0.05};

  int ret = sout::writeRingNum(tmpPath, 1, nRings, covXY, covXZ, covYZ, 5, 1);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "topoMonolayer/coverageAreaXY.dat"));
  REQUIRE(fs::exists(tmpPath + "topoMonolayer/coverageAreaXZ.dat"));
  REQUIRE(fs::exists(tmpPath + "topoMonolayer/coverageAreaYZ.dat"));

  fs::remove_all(tmpPath);
}

// -- writeRingNumBulk tests --

TEST_CASE("writeRingNumBulk writes bulk ring data", "[seams_output]") {
  std::string tmpPath = "/tmp/dseams_test_ringnumbulk/";
  std::vector<int> nRings = {10, 5, 2};

  int ret = sout::writeRingNumBulk(tmpPath, 1, nRings, 5, 1);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "bulkTopo/num_rings.dat"));

  fs::remove_all(tmpPath);
}

// -- writeMoleculeIDsInSlice tests --

TEST_CASE("writeMoleculeIDsInSlice writes molecule selection", "[seams_output]") {
  auto cloud = makeTestCloud(6);
  std::string tmpPath = "/tmp/dseams_test_molids/";

  int ret = sout::writeMoleculeIDsInSlice(tmpPath, cloud);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "selection/IDtextFiles/molID-1.dat"));

  fs::remove_all(tmpPath);
}

// -- writeMoleculeIDsExpressionSelectOVITO tests --

TEST_CASE("writeMoleculeIDsExpressionSelectOVITO writes OVITO selection",
          "[seams_output]") {
  auto cloud = makeTestCloud(6);
  std::string tmpPath = "/tmp/dseams_test_ovito/";

  int ret = sout::writeMoleculeIDsExpressionSelectOVITO(tmpPath, cloud);
  REQUIRE(ret == 0);
  REQUIRE(
      fs::exists(tmpPath + "selection/IDovitoFiles/ovito-molIDSelect-1.dat"));

  fs::remove_all(tmpPath);
}

// -- writeLAMMPSdumpSlice tests --

TEST_CASE("writeLAMMPSdumpSlice writes slice dump", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  std::string tmpPath = "/tmp/dseams_test_dumpslice/";

  int ret = sout::writeLAMMPSdumpSlice(cloud, tmpPath);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "selection/dumpFiles/dump-1.lammpstrj"));

  fs::remove_all(tmpPath);
}

// -- writeLAMMPSdumpINT tests --

TEST_CASE("writeLAMMPSdumpINT writes prism dump with RMSD", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  std::vector<double> rmsdPerAtom = {0.1, 0.2, 0.3, 0.4};
  std::vector<int> atomTypes = {1, 2, 1, 3};
  std::string tmpPath = "/tmp/dseams_test_dumpint/";

  int ret = sout::writeLAMMPSdumpINT(cloud, rmsdPerAtom, atomTypes, 6, tmpPath);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "topoINT/dumpFiles/dump-1.lammpstrj"));

  fs::remove_all(tmpPath);
}

// -- writeLAMMPSdumpCages tests --

TEST_CASE("writeLAMMPSdumpCages writes cage dump with RMSD", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  std::vector<double> rmsdPerAtom = {0.1, 0.2, 0.3, 0.4};
  std::vector<int> atomTypes = {1, 2, 3, 4};
  std::string tmpPath = "/tmp/dseams_test_dumpcages/";

  int ret =
      sout::writeLAMMPSdumpCages(cloud, rmsdPerAtom, atomTypes, tmpPath, 1);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "bulkTopo/dumpFiles/dump-1.lammpstrj"));
  // First frame also creates typeInfo.dat
  REQUIRE(fs::exists(tmpPath + "bulkTopo/typeInfo.dat"));

  fs::remove_all(tmpPath);
}

// -- writeLAMMPSdataAllPrisms tests --

TEST_CASE("writeLAMMPSdataAllPrisms writes data file", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);
  std::vector<int> atomTypes = {1, 1, 1, 1};
  std::string tmpPath = "/tmp/dseams_test_dataprisms/";

  int ret =
      sout::writeLAMMPSdataAllPrisms(cloud, nList, atomTypes, 6, tmpPath);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "topoINT/dataFiles/system-prisms-1.data"));

  fs::remove_all(tmpPath);
}

// -- writeLAMMPSdataAllRings tests --

TEST_CASE("writeLAMMPSdataAllRings writes ring data file", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);
  std::vector<int> atomTypes = {1, 1, 1, 1};
  std::string tmpPath = "/tmp/dseams_test_datarings/";

  int ret =
      sout::writeLAMMPSdataAllRings(cloud, nList, atomTypes, 6, tmpPath);
  REQUIRE(ret == 0);
  REQUIRE(
      fs::exists(tmpPath + "topoMonolayer/dataFiles/system-rings-1.data"));

  fs::remove_all(tmpPath);
}

TEST_CASE("writeLAMMPSdataAllRings bulk mode writes to bulkTopo",
          "[seams_output]") {
  auto cloud = makeTestCloud(4);
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);
  std::vector<int> atomTypes = {1, 1, 1, 1};
  std::string tmpPath = "/tmp/dseams_test_datarings_bulk/";

  int ret = sout::writeLAMMPSdataAllRings(cloud, nList, atomTypes, 6, tmpPath,
                                           false);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "bulkTopo/dataFiles/system-rings-1.data"));

  fs::remove_all(tmpPath);
}

// -- writeLAMMPSdataTopoBulk tests --

TEST_CASE("writeLAMMPSdataTopoBulk writes bulk topo data file",
          "[seams_output]") {
  auto cloud = makeTestCloud(4);
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);
  std::vector<cage::iceType> atomTypes(4, cage::iceType::dummy);
  std::string tmpPath = "/tmp/dseams_test_datatopobulk/";

  int ret = sout::writeLAMMPSdataTopoBulk(cloud, nList, atomTypes, tmpPath);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "bulkTopo/dataFiles/system-1.data"));

  fs::remove_all(tmpPath);
}

// -- writeLAMMPSdata tests --

TEST_CASE("writeLAMMPSdata returns 1 for empty rings", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  std::vector<std::vector<int>> rings;
  std::vector<std::vector<int>> bonds;

  int ret = sout::writeLAMMPSdata(cloud, rings, bonds);
  REQUIRE(ret == 1);
}

TEST_CASE("writeLAMMPSdata writes data for valid rings", "[seams_output]") {
  // Build a cloud where all atoms are used in rings (no padding needed)
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  for (int i = 0; i < 3; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = 1.0 * i;
    pt.y = 0.0;
    pt.z = 0.0;
    cloud.pts.push_back(pt);
    cloud.idIndexMap[pt.atomID] = i;
  }
  cloud.nop = 3;

  std::string outDir = "../output";
  fs::create_directories(outDir);

  std::vector<std::vector<int>> rings = {{1, 2, 3}};
  std::vector<std::vector<int>> bonds = {{1, 2}, {2, 3}};

  int ret = sout::writeLAMMPSdata(cloud, rings, bonds);
  REQUIRE(ret == 0);

  fs::remove_all(outDir);
}

// -- writeXYZcluster tests --

TEST_CASE("writeXYZcluster writes XYZ cluster file", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  std::vector<int> atoms = {0, 1, 2};
  std::string tmpPath = "/tmp/dseams_test_xyzcluster/";

  int ret =
      sout::writeXYZcluster(tmpPath, cloud, atoms, 1, cage::cageType::HexC);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "bulkTopo/clusterXYZ/frame-1/cluster-1.xyz"));

  fs::remove_all(tmpPath);
}

// -- writeHisto tests --

TEST_CASE("writeHisto writes cij q3 q6 files", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);
  // Need c_ij populated
  std::vector<double> avgQ6 = {0.1, 0.2, 0.3, 0.4};

  // Populate c_ij for each atom (needed by writeHisto)
  for (int i = 0; i < cloud.nop; i++) {
    int nNeigh = nList[i].size() - 1;
    for (int j = 0; j < nNeigh; j++) {
      molSys::Result res;
      res.c_value = 0.5;
      res.classifier = molSys::bond_type::staggered;
      cloud.pts[i].c_ij.push_back(res);
    }
  }

  // writeHisto writes to CWD, so we check it doesn't crash
  int ret = sout::writeHisto(cloud, nList, avgQ6);
  REQUIRE(ret == 0);

  // Cleanup files written to CWD
  fs::remove("cij.txt");
  fs::remove("q3.txt");
  fs::remove("q6.txt");
}

// -- writeAllCages tests --

TEST_CASE("writeAllCages returns 1 for empty cage list", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);
  std::vector<cage::Cage> cageList;
  std::vector<std::vector<int>> rings = {{1, 2, 3}};
  std::string tmpPath = "/tmp/dseams_test_allcages/";

  int ret = sout::writeAllCages(tmpPath, cageList, rings, nList, cloud, 1);
  REQUIRE(ret == 1);

  fs::remove_all(tmpPath);
}
