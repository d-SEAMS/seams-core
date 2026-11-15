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
#include <tum_offload.hpp>

#ifdef SEAMS_HAS_OFFLOAD
#include <omp.h>
#endif
#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <array>
#include <cstdlib>
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

TEST_CASE("SEAMS_OFFLOAD TUM cage counts match host on mW cubic",
          "[bulkTUM][offload]") {
  auto yCloud = readMwCubicFrame1();
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);

  setenv("SEAMS_OFFLOAD", "0", 1);
  const auto host = tum::hostCageCounts(nList);
  REQUIRE(host.nSix > 0);
  REQUIRE((host.nHcRings + host.nDdcRings) > 0);
  REQUIRE((host.nHcAtoms + host.nDdcAtoms) > 0);

  setenv("SEAMS_OFFLOAD", "1", 1);
  const auto off = tum::specializedCageCounts(nList);
  REQUIRE(off.retries == 0);
  REQUIRE_FALSE(off.hostFallback);
  REQUIRE(off.nSix == host.nSix);
  REQUIRE(off.nHcRings == host.nHcRings);
  REQUIRE(off.nDdcRings == host.nDdcRings);
  REQUIRE(off.nHcAtoms == host.nHcAtoms);
  REQUIRE(off.nDdcAtoms == host.nDdcAtoms);
  REQUIRE(off.atomHc.size() == host.atomHc.size());
  REQUIRE(off.atomDdc.size() == host.atomDdc.size());
  REQUIRE(off.atomHc == host.atomHc);
  REQUIRE(off.atomDdc == host.atomDdc);

#ifdef SEAMS_HAS_OFFLOAD
  if (omp_get_num_devices() > 0) {
    REQUIRE(off.usedDevice);
  }
#endif
}

namespace {

// Point-line incidence graph of PG(2, q), q prime, with leading self
// entries: girth 6, q^3 (q+1) / 2 six-rings through every vertex
std::vector<std::vector<int>> leviGraph(int q) {
  std::vector<std::array<int, 3>> v;
  for (int a = 0; a < q; ++a) {
    for (int b = 0; b < q; ++b) {
      for (int c = 0; c < q; ++c) {
        const int lead = a != 0 ? a : (b != 0 ? b : c);
        if (lead == 1) {
          v.push_back({a, b, c});
        }
      }
    }
  }
  const int m = static_cast<int>(v.size());
  std::vector<std::vector<int>> nList(static_cast<std::size_t>(2 * m));
  for (int i = 0; i < 2 * m; ++i) {
    nList[static_cast<std::size_t>(i)].push_back(i);
  }
  for (int p = 0; p < m; ++p) {
    const auto &vp = v[static_cast<std::size_t>(p)];
    for (int l = 0; l < m; ++l) {
      const auto &vl = v[static_cast<std::size_t>(l)];
      if ((vp[0] * vl[0] + vp[1] * vl[1] + vp[2] * vl[2]) % q == 0) {
        nList[static_cast<std::size_t>(p)].push_back(m + l);
        nList[static_cast<std::size_t>(m + l)].push_back(p);
      }
    }
  }
  return nList;
}

void requireSameLabels(const tum::CageCounts &a, const tum::CageCounts &b) {
  REQUIRE(a.nSix == b.nSix);
  REQUIRE(a.nHcRings == b.nHcRings);
  REQUIRE(a.nDdcRings == b.nDdcRings);
  REQUIRE(a.atomSix == b.atomSix);
  REQUIRE(a.atomHc == b.atomHc);
  REQUIRE(a.atomDdc == b.atomDdc);
}

tum::SpecializedCaps oneEntryCaps() {
  tum::SpecializedCaps caps;
  caps.rowCap = 1;
  caps.ringCap = 1;
  caps.pairCap = 1;
  return caps;
}

} // namespace

TEST_CASE("hop-bound path grows overflowed buffers to the host labels",
          "[bulkTUM][offload]") {
  auto yCloud = readMwCubicFrame1();
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  const auto host = tum::hostCageCounts(nList);
  for (const bool target : {false, true}) {
    auto caps = oneEntryCaps();
    caps.forceTarget = target;
    const auto grown = tum::specializedCageCounts(nList, caps);
    INFO("target " << target);
    REQUIRE(grown.retries > 0);
    REQUIRE_FALSE(grown.hostFallback);
    requireSameLabels(grown, host);
  }
}

TEST_CASE("hop-bound path on high-incidence Levi graphs matches the host",
          "[bulkTUM][offload]") {
  // q = 2: 12 rings per vertex, inside the row cap. q = 3 and 5: 54 and
  // 375, past it, and for q = 5 past 16n rings in total.
  for (const int q : {2, 3, 5}) {
    INFO("PG(2," << q << ")");
    const auto nList = leviGraph(q);
    const int n = static_cast<int>(nList.size());
    const int perVertex = q * q * q * (q + 1) / 2;
    const auto host = tum::hostCageCounts(nList);
    REQUIRE(host.nSix == perVertex * n / 6);
    REQUIRE(std::all_of(host.atomSix.begin(), host.atomSix.end(),
                        [&](int v) { return v == perVertex; }));
    for (const auto &caps : {tum::SpecializedCaps{}, oneEntryCaps()}) {
      const auto dev = tum::specializedCageCounts(nList, caps);
      REQUIRE(dev.hostFallback == (perVertex > 32));
      requireSameLabels(dev, host);
    }
  }
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

TEST_CASE("clusterCages accepts zero cages", "[bulkTUM]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud.box = {10.0, 10.0, 10.0};
  yCloud.boxLow = {0.0, 0.0, 0.0};
  yCloud.nop = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  pt.atomID = 1;
  yCloud.pts.push_back(pt);
  std::string tmpPath =
      fs::temp_directory_path().append("dseams_test_tum_noring/").string();
  std::vector<std::vector<int>> rings;
  std::vector<cage::Cage> cageList;
  REQUIRE(tum3::clusterCages(yCloud, tmpPath, rings, cageList, 0, 0) == 0);
  std::error_code ec;
  fs::remove_all(tmpPath, ec);
}

TEST_CASE("updateRMSDatom ignores a failed Horn sentinel", "[bulkTUM]") {
  std::vector<std::vector<int>> rings = {{0, 1, 2, 3, 4, 5}};
  cage::Cage cageUnit;
  cageUnit.type = cage::cageType::HexC;
  cageUnit.rings = {0};
  std::vector<double> rmsdPerAtom(6, 0.5);
  std::vector<int> noOfCommonAtoms(6, 1);
  std::vector<cage::iceType> atomTypes(6, cage::iceType::hc);
  REQUIRE(tum3::updateRMSDatom(rings, cageUnit, -1.0, rmsdPerAtom,
                               noOfCommonAtoms, atomTypes) == 0);
  for (int i = 0; i < 6; i++) {
    REQUIRE_THAT(rmsdPerAtom[static_cast<std::size_t>(i)],
                 Catch::Matchers::WithinAbs(0.5, 1e-15));
    REQUIRE(noOfCommonAtoms[static_cast<std::size_t>(i)] == 1);
  }
}
