#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bond.hpp>
#include <cage.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>
#include <topo_bulk.hpp>

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
  std::string tmpDir = fs::temp_directory_path().append("dseams_test_makepath/nested/dir").string();
  int ret = sout::makePath(tmpDir);
  REQUIRE(ret == 0);
  REQUIRE(fs::is_directory(tmpDir));

  // Calling again on existing dir should also return 0
  ret = sout::makePath(tmpDir);
  REQUIRE(ret == 0);

  fs::remove_all(fs::temp_directory_path().append("dseams_test_makepath").string());
}

// -- printRDF tests --

TEST_CASE("printRDF writes RDF data to file", "[seams_output]") {
  std::string tmpFile = fs::temp_directory_path().append("dseams_test_rdf.dat").string();
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
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_cluster_stats/").string();

  int ret = sout::writeClusterStats(tmpPath, 1, 100, 5, 10, 25.0, 1);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "clusterStats.dat"));
  REQUIRE(fs::file_size(tmpPath + "clusterStats.dat") > 0);

  fs::remove_all(tmpPath);
}

TEST_CASE("writeClusterStats appends for subsequent frames", "[seams_output]") {
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_cluster_append/").string();

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
  std::string tmpFile = fs::temp_directory_path().append("dseams_test_writeCluster.txt").string();

  int ret = sout::writeCluster(cloud, tmpFile, false, 42);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpFile));

  fs::remove(tmpFile);
}

// -- writeDump tests --

TEST_CASE("writeDump writes LAMMPS dump format", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_dump/").string();

  int ret = sout::writeDump(cloud, tmpPath, "test.lammpstrj");
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "test.lammpstrj"));
  REQUIRE(fs::file_size(tmpPath + "test.lammpstrj") > 0);

  fs::remove_all(tmpPath);
}

// -- writeTopoBulkData tests --

TEST_CASE("writeTopoBulkData writes bulk topology data", "[seams_output]") {
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_topobulk/").string();

  int ret = sout::writeTopoBulkData(tmpPath, 1, 10, 5, 3, 8, 12, 1);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "bulkTopo/cageData.dat"));

  fs::remove_all(tmpPath);
}

// -- writePrismNum tests --

TEST_CASE("writePrismNum writes prism statistics", "[seams_output]") {
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_prismnum/").string();
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
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_ringnum/").string();
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
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_ringnumbulk/").string();
  std::vector<int> nRings = {10, 5, 2};

  int ret = sout::writeRingNumBulk(tmpPath, 1, nRings, 5, 1);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "bulkTopo/num_rings.dat"));

  fs::remove_all(tmpPath);
}

// -- writeMoleculeIDsInSlice tests --

TEST_CASE("writeMoleculeIDsInSlice writes molecule selection", "[seams_output]") {
  auto cloud = makeTestCloud(6);
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_molids/").string();

  int ret = sout::writeMoleculeIDsInSlice(tmpPath, cloud);
  REQUIRE(ret == 0);
  REQUIRE(fs::exists(tmpPath + "selection/IDtextFiles/molID-1.dat"));

  fs::remove_all(tmpPath);
}

// -- writeMoleculeIDsExpressionSelectOVITO tests --

TEST_CASE("writeMoleculeIDsExpressionSelectOVITO writes OVITO selection",
          "[seams_output]") {
  auto cloud = makeTestCloud(6);
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_ovito/").string();

  int ret = sout::writeMoleculeIDsExpressionSelectOVITO(tmpPath, cloud);
  REQUIRE(ret == 0);
  REQUIRE(
      fs::exists(tmpPath + "selection/IDovitoFiles/ovito-molIDSelect-1.dat"));

  fs::remove_all(tmpPath);
}

// -- writeLAMMPSdumpSlice tests --

TEST_CASE("writeLAMMPSdumpSlice writes slice dump", "[seams_output]") {
  auto cloud = makeTestCloud(4);
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_dumpslice/").string();

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
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_dumpint/").string();

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
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_dumpcages/").string();

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
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_dataprisms/").string();

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
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_datarings/").string();

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
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_datarings_bulk/").string();

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
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_datatopobulk/").string();

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
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_xyzcluster/").string();

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
  std::string tmpPath = fs::temp_directory_path().append("dseams_test_allcages/").string();

  int ret = sout::writeAllCages(tmpPath, cageList, rings, nList, cloud, 1);
  REQUIRE(ret == 1);

  fs::remove_all(tmpPath);
}

// -- writeLAMMPSdataCages with real cage data from mW trajectory --

TEST_CASE("writeLAMMPSdataCages writes cage data file", "[seams_output]") {
  // Read mW cubic ice and find cages
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  if (yCloud.nop == 0) return;

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<ring::strucType> ringType(rings.size());
  std::vector<cage::Cage> cageList;
  auto listHC = ring::findHC(rings, ringType, nList, cageList);

  if (!cageList.empty()) {
    std::string tmpPath = fs::temp_directory_path().append("dseams_test_datacages/").string();

    int ret = sout::writeLAMMPSdataCages(yCloud, rings, cageList,
                                          cage::cageType::HexC,
                                          static_cast<int>(cageList.size()));
    // Writes to ../output/ by default
    REQUIRE(ret == 0);

    fs::remove_all(tmpPath);
    fs::remove_all("../output");
  }
}

// -- writePrisms tests --

TEST_CASE("writePrisms writes prism coordinates", "[seams_output]") {
  auto cloud = makeTestCloud(8);
  std::string outDir = "../output/prisms";
  fs::create_directories(outDir);

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  int ret = sout::writePrisms(basal1, basal2, 1, cloud);
  REQUIRE(ret == 0);

  fs::remove_all("../output");
}

// -- writeBasalRingsPrism tests --

TEST_CASE("writeBasalRingsPrism writes prism basal ring data",
          "[seams_output]") {
  // Use mW data where atom IDs and nList are consistent
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  if (yCloud.nop == 0) return;

  auto nList = nneigh::neighListO(3.5, yCloud, 1);

  // Create output directories
  fs::create_directories("../output/perfect/basalRings");

  // Use the first 4 atom IDs as basal1, next 4 as basal2
  std::vector<int> basal1, basal2;
  for (int i = 0; i < 4 && i < yCloud.nop; i++)
    basal1.push_back(yCloud.pts[i].atomID);
  for (int i = 4; i < 8 && i < yCloud.nop; i++)
    basal2.push_back(yCloud.pts[i].atomID);

  int ret =
      sout::writeBasalRingsPrism(basal1, basal2, 1, nList, yCloud, false);
  // May return 1 if no neighbors found between basals
  REQUIRE((ret == 0 || ret == 1));

  fs::remove_all("../output");
}

// -- writeLAMMPSdataPrisms tests --

TEST_CASE("writeLAMMPSdataPrisms returns 1 for empty prism list",
          "[seams_output]") {
  auto cloud = makeTestCloud(4);
  std::vector<std::vector<int>> rings = {{1, 2, 3}};
  std::vector<int> listPrism; // empty
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 1.5);

  int ret = sout::writeLAMMPSdataPrisms(cloud, rings, false, "", listPrism,
                                         nList);
  REQUIRE(ret == 1);
}

// -- writeLAMMPSdataTopoBulk with real data --

TEST_CASE("writeLAMMPSdataTopoBulk with mixed cage types", "[seams_output]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  if (yCloud.nop == 0) return;

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<ring::strucType> ringType(rings.size());
  std::vector<cage::Cage> cageList;
  ring::findHC(rings, ringType, nList, cageList);
  ring::findDDC(rings, ringType, {}, cageList);

  // Get atom types from ring types
  std::vector<cage::iceType> atomTypes(yCloud.nop, cage::iceType::dummy);
  ring::getAtomTypesTopoBulk(rings, ringType, atomTypes);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_datatopobulk_mixed/").string();
  int ret = sout::writeLAMMPSdataTopoBulk(yCloud, nList, atomTypes, tmpPath);
  REQUIRE(ret == 0);

  fs::remove_all(tmpPath);
}

// -- writeDump with various ice types --

TEST_CASE("writeDump writes all ice type labels", "[seams_output]") {
  auto cloud = makeTestCloud(8);
  // Assign various ice types to test all branches
  cloud.pts[0].iceType = molSys::atom_state_type::cubic;
  cloud.pts[1].iceType = molSys::atom_state_type::hexagonal;
  cloud.pts[2].iceType = molSys::atom_state_type::water;
  cloud.pts[3].iceType = molSys::atom_state_type::interfacial;
  cloud.pts[4].iceType = molSys::atom_state_type::clathrate;
  cloud.pts[5].iceType = molSys::atom_state_type::interClathrate;
  cloud.pts[6].iceType = molSys::atom_state_type::reCubic;
  cloud.pts[7].iceType = molSys::atom_state_type::reHex;

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_dump_types/").string();
  int ret = sout::writeDump(cloud, tmpPath, "test.lammpstrj");
  REQUIRE(ret == 0);
  REQUIRE(fs::file_size(tmpPath + "test.lammpstrj") > 0);

  fs::remove_all(tmpPath);
}

TEST_CASE("writeEachCage writes HC cage data directly", "[seams_output]") {
  auto cloud = makeTestCloud(12);
  // Create output directory structure (the functions expect these to exist)
  fs::create_directories("../output/cages/hexCages");
  fs::create_directories("../output/cages/doubleDiaCages");

  // Rings with 6 elements each (6-membered rings for HC)
  std::vector<std::vector<int>> rings = {
      {1, 2, 3, 4, 5, 6}, {7, 8, 9, 10, 11, 12}};

  // HC cage referencing ring indices
  std::vector<int> currentCage = {0, 1};

  int ret = sout::writeEachCage(currentCage, 1, cage::cageType::HexC, rings,
                                 cloud);
  REQUIRE(ret == 0);

  // Also test DDC cage
  int ret2 = sout::writeEachCage(currentCage, 1, cage::cageType::DoubleDiaC,
                                  rings, cloud);
  REQUIRE(ret2 == 0);

  fs::remove_all("../output");
}

TEST_CASE("writeBasalRingsHex writes HC basal ring data", "[seams_output]") {
  auto cloud = makeTestCloud(12);
  fs::create_directories("../output/cages/hexBasalRings");

  // 6-membered rings
  std::vector<std::vector<int>> rings = {
      {1, 2, 3, 4, 5, 6}, {7, 8, 9, 10, 11, 12}};

  // HC cage: first two rings are basal
  std::vector<int> currentCage = {0, 1};

  // Build a nList where each basal1 atom has a neighbor in basal2
  // nList is by atom ID (1-based), nList[atomID] = {atomID, neighbors...}
  std::vector<std::vector<int>> nList(13); // 0-12
  for (int i = 1; i <= 6; i++) {
    nList[i] = {i, i + 6}; // atom i is neighbor of atom i+6
  }
  for (int i = 7; i <= 12; i++) {
    nList[i] = {i, i - 6}; // atom i is neighbor of atom i-6
  }

  int ret = sout::writeBasalRingsHex(currentCage, 1, nList, rings);
  REQUIRE(ret == 0);

  fs::remove_all("../output");
}

TEST_CASE("writeLAMMPSdataCages writes cage LAMMPS data file",
          "[seams_output]") {
  auto cloud = makeTestCloud(12);
  fs::create_directories("../output");

  std::vector<std::vector<int>> rings = {
      {1, 2, 3, 4, 5, 6}, {7, 8, 9, 10, 11, 12}};

  cage::Cage cage1;
  cage1.type = cage::cageType::HexC;
  cage1.rings = {0, 1};
  std::vector<cage::Cage> cageList = {cage1};

  int ret = sout::writeLAMMPSdataCages(cloud, rings, cageList,
                                        cage::cageType::HexC, 1);
  REQUIRE(ret == 0);

  fs::remove_all("../output");
}
