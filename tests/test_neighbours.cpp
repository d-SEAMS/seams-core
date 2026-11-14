#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#ifdef SEAMS_HAS_VESIN
#include <vesin.h>
#endif

#include <bop.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <array>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

// Helper: build a 4-atom system in a 10x10x10 box
// Atoms at (0,0,0), (1,0,0), (0,1,0), (5,5,5)
// With cutoff 1.5, atoms 0,1,2 are mutual neighbours; atom 3 is isolated
static molSys::PointCloud<molSys::Point<double>, double> makeFourAtomCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  double coords[4][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {5.0, 5.0, 5.0}};

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

TEST_CASE("neighListO produces full neighbour list by atom ID", "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  double cutoff = 1.5;

  auto nList = nneigh::neighListO(cutoff, cloud, 1);

  // Should have one row per atom
  REQUIRE(nList.size() == 4);

  // Atom 0 (ID 0): neighbours should include IDs 1 and 2
  // First element is the atom's own ID
  REQUIRE(nList[0][0] == 0);
  REQUIRE(nList[0].size() >= 3); // self + 2 neighbours
  // Check that 1 and 2 appear in the neighbour list of atom 0
  bool found1 =
      std::find(nList[0].begin() + 1, nList[0].end(), 1) != nList[0].end();
  bool found2 =
      std::find(nList[0].begin() + 1, nList[0].end(), 2) != nList[0].end();
  REQUIRE(found1);
  REQUIRE(found2);

  // Atom 3 (ID 3): should have no neighbours (distance > cutoff from all)
  REQUIRE(nList[3][0] == 3);
  REQUIRE(nList[3].size() == 1); // only self
}

TEST_CASE("halfNeighList produces half neighbour list", "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  double cutoff = 1.5;

  auto nList = nneigh::halfNeighList(cutoff, cloud, 1);

  REQUIRE(nList.size() == 4);

  // In a half list, if 0 lists 1 as neighbour, 1 should NOT list 0
  // Atom 0 should have 1 and 2 as neighbours
  bool zeroHas1 =
      std::find(nList[0].begin() + 1, nList[0].end(), 1) != nList[0].end();
  REQUIRE(zeroHas1);

  // Atom 1 should NOT have 0 as neighbour (half list)
  bool oneHas0 =
      std::find(nList[1].begin() + 1, nList[1].end(), 0) != nList[1].end();
  REQUIRE_FALSE(oneHas0);
}

TEST_CASE("neighbourListByIndex converts ID list to index list",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  double cutoff = 1.5;

  // Get a neighbour list by ID
  auto nListByID = nneigh::neighListO(cutoff, cloud, 1);
  // Convert to index-based
  auto nListByIdx = nneigh::neighbourListByIndex(cloud, nListByID);

  REQUIRE(nListByIdx.size() == 4);

  // Since atomID == index in this test cloud, the lists should be identical
  for (std::size_t i = 0; i < 4; i++) {
    REQUIRE(nListByIdx[i].size() == nListByID[i].size());
  }
}

TEST_CASE("getNewNeighbourListByIndex builds index-based list from scratch",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  double cutoff = 1.5;

  auto nList = nneigh::getNewNeighbourListByIndex(cloud, cutoff);

  REQUIRE(nList.size() == 4);
  // Atom 0 (index 0) should have indices 1 and 2 as neighbours
  REQUIRE(nList[0][0] == 0);
  bool found1 =
      std::find(nList[0].begin() + 1, nList[0].end(), 1) != nList[0].end();
  bool found2 =
      std::find(nList[0].begin() + 1, nList[0].end(), 2) != nList[0].end();
  REQUIRE(found1);
  REQUIRE(found2);
}

// Helper: build a cloud with two atom types
static molSys::PointCloud<molSys::Point<double>, double>
makeTwoTypeCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  // 4 atoms: types 1 and 2, close together
  int types[] = {1, 2, 1, 2};
  double coords[4][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {5.0, 5.0, 5.0}};

  for (int i = 0; i < 4; i++) {
    molSys::Point<double> pt;
    pt.type = types[i];
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

TEST_CASE("neighList produces neighbour list for two atom types",
          "[neighbours]") {
  auto cloud = makeTwoTypeCloud();
  double cutoff = 1.5;

  // Find neighbours between type 1 and type 2
  auto nList = nneigh::neighList(cutoff, cloud, 1, 2);

  // Should have one row per atom (all atom types)
  REQUIRE(nList.size() == 4);

  // Atom 0 (type 1) should have atom 1 (type 2) as neighbour (distance 1.0)
  bool found1 =
      std::find(nList[0].begin() + 1, nList[0].end(), 1) != nList[0].end();
  REQUIRE(found1);

  // Same-type pair 0-2 is inside the cutoff but must not appear on an I-J list
  bool foundSameType =
      std::find(nList[0].begin() + 1, nList[0].end(), 2) != nList[0].end();
  REQUIRE_FALSE(foundSameType);

  // Atom 3 (type 2, far away) should have no type-1 neighbours within cutoff
  REQUIRE(nList[3].size() == 1); // only self
}

TEST_CASE("clearNeighbourList empties the list", "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  auto nList = nneigh::neighListO(1.5, cloud, 1);
  REQUIRE_FALSE(nList.empty());

  nneigh::clearNeighbourList(nList);
  REQUIRE(nList.empty());
}

// ---------------------------------------------------------------------------
// A cell list enumerates periodic images; the minimum image convention does
// not. Once the box stops being larger than twice the cutoff the two disagree
// unless the cell-list output is reduced: a neighbour reachable through more
// than one image would otherwise be listed more than once, and below the
// cutoff a particle becomes its own neighbour. Every other fixture in this
// suite uses a box comfortably larger than 2*cutoff, so nothing else covers
// this.
// ---------------------------------------------------------------------------

namespace {

molSys::PointCloud<molSys::Point<double>, double>
smallBoxCloud(double boxLength) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const std::vector<std::array<double, 3>> pos = {
      {0.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {1.0, 1.0, 1.0}};

  cloud.nop = static_cast<int>(pos.size());
  cloud.currentFrame = 1;
  cloud.box = {boxLength, boxLength, boxLength};
  cloud.boxLow = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < pos.size(); i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = static_cast<int>(i) + 1;
    p.molID = p.atomID;
    p.x = pos[i][0];
    p.y = pos[i][1];
    p.z = pos[i][2];
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = static_cast<int>(i);
  }
  return cloud;
}

//! Neighbours of iatom under the minimum image convention, sorted, no self
std::vector<int> minimumImageNeighbours(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud, int iatom,
    double cutoff) {
  std::vector<int> ref;
  for (int j = 0; j < cloud.nop; j++) {
    if (j == iatom) {
      continue;
    }
    if (gen::periodicDist(cloud, iatom, j) <= cutoff) {
      ref.push_back(j);
    }
  }
  std::sort(ref.begin(), ref.end());
  return ref;
}

} // namespace

TEST_CASE("neighbour list matches minimum image for a tight box",
          "[neighbours]") {
  const double cutoff = 3.5;

  // 30.0 is the comfortable case; 6.0 puts two images of a pair inside the
  // sphere; 3.0 puts a particle inside its own image
  for (double boxLength : {30.0, 6.0, 3.0}) {
    auto cloud = smallBoxCloud(boxLength);
    auto nList = nneigh::getNewNeighbourListByIndex(cloud, cutoff);

    REQUIRE(nList.size() == static_cast<size_t>(cloud.nop));

    for (int iatom = 0; iatom < cloud.nop; iatom++) {
      INFO("box = " << boxLength << ", atom " << iatom);

      // Row header is the particle itself
      const auto &row = nList[static_cast<std::size_t>(iatom)];
      REQUIRE(row.size() >= 1);
      REQUIRE(row[0] == iatom);

      std::vector<int> got(row.begin() + 1, row.end());
      std::sort(got.begin(), got.end());

      // No duplicates, and never the particle itself
      REQUIRE(std::adjacent_find(got.begin(), got.end()) == got.end());
      REQUIRE(std::find(got.begin(), got.end(), iatom) == got.end());

      REQUIRE(got == minimumImageNeighbours(cloud, iatom, cutoff));
    }
  }
}

TEST_CASE("neighbour builders return empty when the cloud has no box",
          "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = 2;
  cloud.currentFrame = 1;
  for (int i = 0; i < 2; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.x = static_cast<double>(i);
    pt.y = 0.0;
    pt.z = 0.0;
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }

  REQUIRE(cloud.box.empty());
  REQUIRE(nneigh::neighListO(3.5, cloud, 1).empty());
  REQUIRE(nneigh::neighList(3.5, cloud, 1, 1).empty());
  REQUIRE(nneigh::neighListPair(3.5, cloud, 1, 2).empty());
  REQUIRE(nneigh::halfNeighList(3.5, cloud, 1).empty());
  REQUIRE(nneigh::getNewNeighbourListByIndex(cloud, 3.5).empty());
}

TEST_CASE("neighListPair finds unlike pair in an orthogonal cell",
          "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double lx = 10.0;
  cloud.box = {lx, lx, lx};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 2;
  cloud.currentFrame = 1;
  const double coords[2][3] = {{0.0, 0.0, 0.0}, {0.5 * lx, 0.0, 0.0}};
  const int types[2] = {1, 2};
  for (int i = 0; i < 2; i++) {
    molSys::Point<double> pt;
    pt.type = types[i];
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }

  auto unlike = nneigh::neighListPair(6.0, cloud, 1, 2);
  REQUIRE(unlike.size() == 2);
  REQUIRE(unlike[0][0] == 1);
  REQUIRE(unlike[1][0] == 2);
  const bool oneHasTwo =
      std::find(unlike[0].begin() + 1, unlike[0].end(), 2) != unlike[0].end();
  const bool twoHasOne =
      std::find(unlike[1].begin() + 1, unlike[1].end(), 1) != unlike[1].end();
  REQUIRE(oneHasTwo);
  REQUIRE(twoHasOne);

  auto like = nneigh::neighListPair(6.0, cloud, 1, 1);
  REQUIRE(like.size() == 2);
  REQUIRE(like[0][0] == 1);
  REQUIRE(like[0].size() == 1);
  const bool likeHasType2 =
      std::find(like[0].begin() + 1, like[0].end(), 2) != like[0].end();
  REQUIRE_FALSE(likeHasType2);
}

TEST_CASE("neighListPair sees the sheared a-image unlike pair",
          "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 2;
  const double coords[2][3] = {{0.2, 0.1, 1.0}, {9.7, 0.1, 1.0}};
  const int types[2] = {1, 2};
  for (int i = 0; i < 2; i++) {
    molSys::Point<double> pt;
    pt.type = types[i];
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }
  REQUIRE_THAT(gen::periodicDistSq(cloud, 0, 1),
               Catch::Matchers::WithinAbs(0.25, 1e-9));
  auto nList = nneigh::neighListPair(1.0, cloud, 1, 2);
  REQUIRE(nList.size() == 2);
  const bool oneHasTwo =
      std::find(nList[0].begin() + 1, nList[0].end(), 2) != nList[0].end();
  const bool twoHasOne =
      std::find(nList[1].begin() + 1, nList[1].end(), 1) != nList[1].end();
  REQUIRE(oneHasTwo);
  REQUIRE(twoHasOne);
}

TEST_CASE("neighList same type does not self-include or double-count",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  auto nList = nneigh::neighList(1.5, cloud, 1, 1);
  REQUIRE(nList.size() == 4);
  for (int i = 0; i < 4; i++) {
    REQUIRE(nList[static_cast<std::size_t>(i)][0] == i);
    auto neigh = nList[static_cast<std::size_t>(i)];
    std::sort(neigh.begin() + 1, neigh.end());
    REQUIRE(std::adjacent_find(neigh.begin() + 1, neigh.end()) == neigh.end());
    REQUIRE(std::find(neigh.begin() + 1, neigh.end(), i) == neigh.end());
  }
  REQUIRE(nList[0].size() == 3);
  REQUIRE(nList[3].size() == 1);
}

TEST_CASE("neighbourListByIndex skips empty rows and unknown IDs",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  cloud.idIndexMap.erase(2);
  std::vector<std::vector<int>> nList = {{0, 1, 2}, {}, {2, 0}, {3, 99}};
  auto idx = nneigh::neighbourListByIndex(cloud, nList);
  REQUIRE(idx.size() == 4);
  REQUIRE(idx[0][0] == 0);
  REQUIRE(std::find(idx[0].begin() + 1, idx[0].end(), 1) != idx[0].end());
  REQUIRE(std::find(idx[0].begin() + 1, idx[0].end(), 2) == idx[0].end());
  REQUIRE(idx[1].empty());
  REQUIRE(idx[2].empty());
  REQUIRE(idx[3].size() == 1);
  REQUIRE(idx[3][0] == 3);
}

TEST_CASE("neighListO leaves an unmapped atom as an empty row",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  cloud.idIndexMap.erase(2);

  auto nList = nneigh::neighListO(1.5, cloud, 1);
  REQUIRE(nList.size() == 4);
  REQUIRE_FALSE(nList[0].empty());
  REQUIRE(nList[0][0] == 0);
  REQUIRE(nList[2].empty());
}

TEST_CASE("SkinNeighborList matches neighListO on a static cloud",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  nneigh::SkinNeighborList skin(1.5, 0.5, 1, nneigh::BondGraph::Cutoff);
  const auto &bonds = skin.update(cloud);
  auto hard = nneigh::neighListO(1.5, cloud, 1);
  REQUIRE(skin.lastRebuilt());
  REQUIRE(bonds.size() == hard.size());
  for (std::size_t i = 0; i < bonds.size(); i++) {
    auto a = bonds[i];
    auto b = hard[i];
    std::sort(a.begin() + 1, a.end());
    std::sort(b.begin() + 1, b.end());
    REQUIRE(a == b);
  }
}

TEST_CASE("SkinNeighborList rebuilds only after the Verlet trigger",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  nneigh::SkinNeighborList skin(1.5, 1.0, 1, nneigh::BondGraph::Cutoff);
  (void)skin.update(cloud);
  REQUIRE(skin.lastRebuilt());

  cloud.pts[1].x += 0.2; // skin/2 = 0.5
  (void)skin.update(cloud);
  REQUIRE_FALSE(skin.lastRebuilt());

  cloud.pts[1].x += 0.5;
  (void)skin.update(cloud);
  REQUIRE(skin.lastRebuilt());
}

TEST_CASE("SkinNeighborList drops a bond as soon as it leaves the cutoff",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  nneigh::SkinNeighborList skin(1.5, 2.0, 1, nneigh::BondGraph::Cutoff);
  (void)skin.update(cloud);
  REQUIRE(std::find(skin.bonds()[0].begin() + 1, skin.bonds()[0].end(), 1) !=
          skin.bonds()[0].end());

  cloud.pts[1].x = 1.8; // 1.8 > 1.5, still inside the skin halo
  (void)skin.update(cloud);
  REQUIRE_FALSE(skin.lastRebuilt());
  REQUIRE(std::find(skin.bonds()[0].begin() + 1, skin.bonds()[0].end(), 1) ==
          skin.bonds()[0].end());
}

TEST_CASE("neighListO returns empty when the cloud has no atoms",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  cloud.nop = 0;
  cloud.pts.clear();
  cloud.idIndexMap.clear();

  auto nList = nneigh::neighListO(3.5, cloud, 1);
  REQUIRE(nList.empty());
  auto knn = nneigh::kNearestNeighbourList(cloud, 4, 5.5, 1, true);
  REQUIRE(knn.empty());
  auto both = nneigh::kNearestNeighbourPair(cloud, 4, 5.5, 1);
  REQUIRE(both.first.empty());
  REQUIRE(both.second.empty());
}

TEST_CASE("k-nearest graph reduces exactly to the cutoff graph when the "
          "shell-separation certificate holds",
          "[neighbours]") {
  // Theorem: when max_i d_k(i) <= rcutoff <= min_i d_{k+1}(i), the cutoff
  // neighbourhood of every particle is exactly its k nearest, so the
  // union-symmetrized k-nearest graph equals the cutoff graph edge for edge
  // and every downstream graph predicate is identical. Checked here on the
  // thermal mW frame at k = 4, rcutoff = 3.5.
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);

  const auto [maxKth, minNext] = nneigh::shellSeparation(yCloud, 4, 1);
  REQUIRE(maxKth > 0.0);
  REQUIRE(maxKth <= 3.5);
  REQUIRE(minNext >= 3.5);

  auto cutoffRows = nneigh::neighListO(3.5, yCloud, 1);
  auto knnRows = nneigh::kNearestNeighbourList(yCloud, 4, 3.5, 1);
  REQUIRE(cutoffRows.size() == knnRows.size());
  for (size_t i = 0; i < cutoffRows.size(); i++) {
    REQUIRE_FALSE(cutoffRows[i].empty());
    REQUIRE(cutoffRows[i][0] == knnRows[i][0]);
    std::vector<int> a(cutoffRows[i].begin() + 1, cutoffRows[i].end());
    std::vector<int> b(knnRows[i].begin() + 1, knnRows[i].end());
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    REQUIRE(a == b);
  }
}

TEST_CASE("k-nearest pair matches two separate k-nearest builds",
          "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  auto mutual = nneigh::kNearestNeighbourList(yCloud, 4, 5.5, 1, true);
  auto uni = nneigh::kNearestNeighbourList(yCloud, 4, 5.5, 1, false);
  auto both = nneigh::kNearestNeighbourPair(yCloud, 4, 5.5, 1);
  REQUIRE(both.first == mutual);
  REQUIRE(both.second == uni);
}

TEST_CASE("k-nearest graph is exact under an undersized candidate cutoff",
          "[neighbours]") {
  // The brute-force fallback recomputes any atom whose k-th neighbour lies
  // beyond the candidate cutoff, so the nominations are exact regardless
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);

  auto generous = nneigh::kNearestNeighbourList(yCloud, 4, 3.5, 1);
  auto starved = nneigh::kNearestNeighbourList(yCloud, 4, 1.0, 1);
  REQUIRE(generous.size() == starved.size());
  for (size_t i = 0; i < generous.size(); i++) {
    std::vector<int> a(generous[i].begin() + 1, generous[i].end());
    std::vector<int> b(starved[i].begin() + 1, starved[i].end());
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    REQUIRE(a == b);
  }
}

TEST_CASE("mutual and union k-nearest graphs coincide on the crystal",
          "[neighbours]") {
  // In a crystal the first shell is mutual, so the intersection and union
  // symmetrizations produce the same graph, and both reduce to the cutoff
  // graph under the shell-separation certificate
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  auto unionRows = nneigh::kNearestNeighbourList(yCloud, 4, 3.5, 1, false);
  auto mutualRows = nneigh::kNearestNeighbourList(yCloud, 4, 3.5, 1, true);
  REQUIRE(unionRows.size() == mutualRows.size());
  for (size_t i = 0; i < unionRows.size(); i++) {
    std::vector<int> a(unionRows[i].begin() + 1, unionRows[i].end());
    std::vector<int> b(mutualRows[i].begin() + 1, mutualRows[i].end());
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    REQUIRE(a == b);
  }
}

TEST_CASE("bondGraphFromName accepts the published graph names",
          "[neighbours]") {
  REQUIRE(nneigh::bondGraphFromName("cutoff") == nneigh::BondGraph::Cutoff);
  REQUIRE(nneigh::bondGraphFromName("knn") == nneigh::BondGraph::KnnMutual);
  REQUIRE(nneigh::bondGraphFromName("knn-union") ==
          nneigh::BondGraph::KnnUnion);
  REQUIRE(std::string(nneigh::bondGraphName(nneigh::BondGraph::KnnMutual)) ==
          "knn");
  REQUIRE_THROWS_AS(nneigh::bondGraphFromName("not-a-graph"),
                    std::invalid_argument);
}

TEST_CASE("SkinNeighborList default is the mutual four-nearest graph",
          "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  nneigh::SkinNeighborList skin(3.5, 2.0, 1);
  const auto &bonds = skin.update(yCloud);
  auto knn = nneigh::kNearestNeighbourList(yCloud, 4, 5.5, 1, true);
  REQUIRE(bonds.size() == knn.size());
  for (std::size_t i = 0; i < bonds.size(); i++) {
    std::vector<int> a(bonds[i].begin() + 1, bonds[i].end());
    std::vector<int> b(knn[i].begin() + 1, knn[i].end());
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    REQUIRE(a == b);
  }
}

TEST_CASE("dumpAxialDim uses recovered lx ly lz not tilt", "[neighbours]") {
  double lengths[3];
  nneigh::dumpCellLengths({61.0, 12.0, 50.0, 60.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
                          lengths);
  REQUIRE_THAT(lengths[0], Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(lengths[1], Catch::Matchers::WithinAbs(12.0, 1e-12));
  REQUIRE_THAT(lengths[2], Catch::Matchers::WithinAbs(50.0, 1e-12));
  REQUIRE(nneigh::dumpAxialDim({61.0, 12.0, 50.0, 60.0, 0.0, 0.0},
                               {0.0, 0.0, 0.0}) == 2);
  REQUIRE(nneigh::dumpAxialDim({10.0, 10.0, 50.0}, {0.0, 0.0, 0.0}) == 2);
}

TEST_CASE("empty cloud k-nearest does not throw", "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.nop = 0;
  auto rows = nneigh::kNearestNeighbourList(cloud, 4, 5.5, 1, true);
  REQUIRE(rows.empty());
  auto both = nneigh::kNearestNeighbourPair(cloud, 4, 5.5, 1);
  REQUIRE(both.first.empty());
  REQUIRE(both.second.empty());
}

#ifdef SEAMS_HAS_LINKCELL
TEST_CASE("LAMMPS dump bounds convert to lattice H",
          "[neighbours][linkcell]") {
  // Dump stores bound spans, not edge lengths. These bounds recover
  // lx=10, ly=11, lz=12, xy=1.5, xz=0.4, yz=-0.2, origin (1,2,3).
  const lc_cell box = nneigh::lammpsBoxToLcCell(
      {11.9, 11.2, 12.0, 1.5, 0.4, -0.2}, {1.0, 1.8, 3.0});
  REQUIRE_THAT(box.ax, Catch::Matchers::WithinAbs(10.0, 1e-12));
  REQUIRE(box.ay == 0.0);
  REQUIRE(box.az == 0.0);
  REQUIRE(box.bx == 1.5);
  REQUIRE_THAT(box.by, Catch::Matchers::WithinAbs(11.0, 1e-12));
  REQUIRE(box.bz == 0.0);
  REQUIRE(box.cx == 0.4);
  REQUIRE(box.cy == -0.2);
  REQUIRE(box.cz == 12.0);
  REQUIRE_THAT(box.ox, Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(box.oy, Catch::Matchers::WithinAbs(2.0, 1e-12));
  REQUIRE(box.oz == 3.0);
}

TEST_CASE("sheared box k-nearest uses the periodic image",
          "[neighbours][linkcell]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  // Bound span 15 with xy=5 recovers lx=10, b=(5, 8.66, 0).
  cloud.box = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 2;
  const double coords[2][3] = {{0.2, 0.1, 1.0}, {9.7, 0.1, 1.0}};
  for (int i = 0; i < 2; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }
  auto knn = nneigh::kNearestNeighbourList(cloud, 1, 5.5, 1, true);
  REQUIRE(knn.size() == 2);
  REQUIRE(knn[0].size() >= 2);
  REQUIRE(knn[1].size() >= 2);
  REQUIRE(knn[0][0] == 1);
  REQUIRE(knn[1][0] == 2);
  REQUIRE(knn[0][1] == 2);
  REQUIRE(knn[1][1] == 1);
  // Bound span 15 is not lx. The a-image distance is 0.5.
  REQUIRE_THAT(gen::periodicDistSq(cloud, 0, 1),
               Catch::Matchers::WithinAbs(0.25, 1e-9));
}
#endif

TEST_CASE("tilt dump cutoff list finds the a-image pair", "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 2;
  const double coords[2][3] = {{0.2, 0.1, 1.0}, {9.7, 0.1, 1.0}};
  for (int i = 0; i < 2; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }
  REQUIRE_THAT(gen::periodicDistSq(cloud, 0, 1),
               Catch::Matchers::WithinAbs(0.25, 1e-9));
  auto nList = nneigh::neighListO(1.0, cloud, 1);
  REQUIRE(nList.size() == 2);
  REQUIRE(nList[0][0] == 1);
  REQUIRE(nList[1][0] == 2);
  const bool zeroHasOne =
      std::find(nList[0].begin() + 1, nList[0].end(), 2) != nList[0].end();
  const bool oneHasZero =
      std::find(nList[1].begin() + 1, nList[1].end(), 1) != nList[1].end();
  REQUIRE(zeroHasOne);
  REQUIRE(oneHasZero);

  auto byIndex = nneigh::getNewNeighbourListByIndex(cloud, 1.0);
  REQUIRE(byIndex[0][0] == 0);
  const bool idx01 =
      std::find(byIndex[0].begin() + 1, byIndex[0].end(), 1) !=
      byIndex[0].end();
  REQUIRE(idx01);

  nneigh::SkinNeighborList skin(1.0, 0.5, 1, nneigh::BondGraph::Cutoff);
  const auto &skinList = skin.update(cloud);
  REQUIRE(skin.lastRebuilt());
  const bool skin01 =
      std::find(skinList[0].begin() + 1, skinList[0].end(), 2) !=
      skinList[0].end();
  REQUIRE(skin01);
  cloud.box[3] = 4.0;
  skin.update(cloud);
  REQUIRE(skin.lastRebuilt());
}

// Type 1 / type 2 ion vertices (cation, anion). Not a cutoff shell.
static molSys::PointCloud<molSys::Point<double>, double>
makeFourIonCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 4;
  cloud.currentFrame = 1;
  // 0 --1.0-- 1 --2.0-- 2 --3.0-- 3
  // type 1    2         1         2
  // Mutual nearest unlike: (0, 1). Ion 2's nearest unlike is 1,
  // whose nearest type-1 is 0, so (2, 1) is not mutual.
  const int types[4] = {1, 2, 1, 2};
  const double x[4] = {0.0, 1.0, 3.0, 6.0};
  const int atomIDs[4] = {10, 20, 30, 40};
  for (int i = 0; i < 4; i++) {
    molSys::Point<double> pt;
    pt.type = types[i];
    pt.atomID = atomIDs[i];
    pt.molID = atomIDs[i];
    pt.x = x[i];
    pt.y = 0.0;
    pt.z = 0.0;
    cloud.pts.push_back(pt);
    cloud.idIndexMap[pt.atomID] = i;
  }
  return cloud;
}

TEST_CASE("nearestUnlike on 2+2 ions has one mutual pair",
          "[neighbours]") {
  auto cloud = makeFourIonCloud();
  const auto nearest = nneigh::nearestUnlike(cloud, 1, 2);
  REQUIRE(nearest.size() == 2);
  REQUIRE(std::get<0>(nearest[0]) == 0);
  REQUIRE(std::get<1>(nearest[0]) == 1);
  REQUIRE_THAT(std::get<2>(nearest[0]),
               Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE(std::get<0>(nearest[1]) == 2);
  REQUIRE(std::get<1>(nearest[1]) == 1);
  REQUIRE_THAT(std::get<2>(nearest[1]),
               Catch::Matchers::WithinAbs(2.0, 1e-12));

  const auto pairs = nneigh::mutualNearestUnlike(cloud, 1, 2);
  REQUIRE(pairs.size() == 1);
  REQUIRE(pairs[0].first == 0);
  REQUIRE(pairs[0].second == 1);
}

TEST_CASE("nearestUnlike picks the sheared a-image unlike ion",
          "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 3;
  // Type-1 at the origin side. Type-2 at 9.7 is the a-image partner
  // (MIC 0.5). The decoy type-2 at 1.0 is closer in the unwrapped
  // cell (0.8) and must not win.
  const double coords[3][3] = {
      {0.2, 0.1, 1.0}, {9.7, 0.1, 1.0}, {1.0, 0.1, 1.0}};
  const int types[3] = {1, 2, 2};
  for (int i = 0; i < 3; i++) {
    molSys::Point<double> pt;
    pt.type = types[i];
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }
  REQUIRE_THAT(gen::periodicDistSq(cloud, 0, 1),
               Catch::Matchers::WithinAbs(0.25, 1e-9));
  REQUIRE(gen::periodicDistSq(cloud, 0, 2) > 0.25);

  const auto nearest = nneigh::nearestUnlike(cloud, 1, 2);
  REQUIRE(nearest.size() == 1);
  REQUIRE(std::get<0>(nearest[0]) == 0);
  REQUIRE(std::get<1>(nearest[0]) == 1);
  REQUIRE_THAT(std::get<2>(nearest[0]),
               Catch::Matchers::WithinAbs(0.5, 1e-9));

  const auto pairs = nneigh::mutualNearestUnlike(cloud, 1, 2);
  REQUIRE(pairs.size() == 1);
  REQUIRE(pairs[0].first == 0);
  REQUIRE(pairs[0].second == 1);
}

TEST_CASE("nearestUnlike finds an image the tilted wrap overshoots",
          "[neighbours]") {
  // yz = ly/2. The partner at 4.686 wraps to 7.07; the decoy at 4.8 wraps
  // exactly and sits below half of every H diagonal, not half the width.
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 15.0, 10.0, 0.0, 0.0, 5.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 3;
  const double coords[3][3] = {{5.0, 5.0, 2.0}, {5.0, 1.4, 5.0}, {9.8, 5.0, 2.0}};
  const int types[3] = {1, 2, 2};
  for (int i = 0; i < 3; i++) {
    molSys::Point<double> pt;
    pt.type = types[i];
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }
  const auto nearest = nneigh::nearestUnlike(cloud, 1, 2);
  REQUIRE(nearest.size() == 1);
  REQUIRE(std::get<1>(nearest[0]) == 1);
  REQUIRE_THAT(std::get<2>(nearest[0]),
               Catch::Matchers::WithinAbs(std::sqrt(3.6 * 3.6 + 9.0), 1e-12));
}

TEST_CASE("nearestUnlike matches the brute-force minimum image at scale",
          "[neighbours]") {
  // 300 ions of each type, orthorhombic and tilted on every axis; in the
  // 2 A slab most nearest partners lie past the Smith ball.
  const std::array<double, 6> frames[] = {{20, 20, 20, 0, 0, 0},
                                          {20, 20, 20, 3, -2, 4},
                                          {60, 60, 2, 9, -6, 0.5}};
  for (const auto &f : frames) {
    const double xy = f[3], xz = f[4], yz = f[5];
    const double xmin = std::min({0.0, xy, xz, xy + xz});
    const double xmax = std::max({0.0, xy, xz, xy + xz});
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    cloud.box = {f[0] + xmax - xmin, f[1] + std::max(0.0, yz) - std::min(0.0, yz),
                 f[2], xy, xz, yz};
    cloud.boxLow = {xmin, std::min(0.0, yz), 0.0};
    if (xy == 0.0 && xz == 0.0 && yz == 0.0) {
      cloud.box.resize(3);
    }
    std::mt19937 rng(5);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (int i = 0; i < 600; i++) {
      const double sx = u(rng);
      const double sy = u(rng);
      const double sz = u(rng);
      molSys::Point<double> pt;
      pt.type = 1 + i % 2;
      pt.atomID = i + 1;
      pt.x = f[0] * sx + xy * sy + xz * sz;
      pt.y = f[1] * sy + yz * sz;
      pt.z = f[2] * sz;
      cloud.pts.push_back(pt);
      cloud.idIndexMap[i + 1] = i;
    }
    cloud.nop = 600;
    const gen::FracBox frame = gen::makeFracBox(cloud);
    const auto got = nneigh::nearestUnlike(cloud, 1, 2);
    REQUIRE(got.size() == 300);
    for (const auto &[i, j, d] : got) {
      int bestJ = -1;
      double bestD2 = std::numeric_limits<double>::infinity();
      for (int k = 1; k < cloud.nop; k += 2) {
        const double d2 = gen::periodicDistSq(frame, cloud, i, k);
        if (d2 < bestD2) {
          bestD2 = d2;
          bestJ = k;
        }
      }
      REQUIRE(j == bestJ);
      REQUIRE_THAT(d, Catch::Matchers::WithinAbs(std::sqrt(bestD2), 1e-9));
    }
  }
}

TEST_CASE("cell-list neighbours match the minimum image on a tilted frame",
          "[neighbours]") {
  // Water density and fewer atoms than the threaded rows take, so both
  // lists come from cutoffRows across many cells.
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {44.0, 40.0, 36.0, 5.0, -3.0, 4.0};
  cloud.boxLow = {-3.0, 0.0, 0.0};
  constexpr int n = 1500;
  std::mt19937 rng(17);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int i = 0; i < n; i++) {
    const double sx = u(rng);
    const double sy = u(rng);
    const double sz = u(rng);
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = 36.0 * sx + 5.0 * sy - 3.0 * sz;
    pt.y = 36.0 * sy + 4.0 * sz;
    pt.z = 36.0 * sz;
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }
  cloud.nop = n;
  const double cutoff = 3.0;
  const gen::FracBox frame = gen::makeFracBox(cloud);
  std::vector<std::vector<int>> wantFull(n), wantHalf(n);
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (gen::periodicDistSq(frame, cloud, i, j) <= cutoff * cutoff) {
        wantFull[static_cast<std::size_t>(i)].push_back(j + 1);
        wantFull[static_cast<std::size_t>(j)].push_back(i + 1);
        wantHalf[static_cast<std::size_t>(i)].push_back(j + 1);
      }
    }
  }
  const auto rows = [](const std::vector<std::vector<int>> &list) {
    std::vector<std::vector<int>> out;
    for (const auto &row : list) {
      std::vector<int> ids(row.begin() + 1, row.end());
      std::sort(ids.begin(), ids.end());
      out.push_back(ids);
    }
    return out;
  };
  for (auto &row : wantFull) {
    std::sort(row.begin(), row.end());
  }
  REQUIRE(rows(nneigh::neighListO(cutoff, cloud, 1)) == wantFull);
  REQUIRE(rows(nneigh::halfNeighList(cutoff, cloud, 1)) == wantHalf);
}

#ifdef SEAMS_HAS_LINKCELL
TEST_CASE("residentFrameCell uses dump H when nBox is 6", "[neighbours]") {
  const double box[6] = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  linkcell::Cell cell = linkcell::Cell::ortho(1.0, 1.0, 1.0);
  const double *frameLens = box;
  nneigh::residentFrameCell(box, nullptr, 6, cell, frameLens);
  REQUIRE(frameLens == nullptr);
  REQUIRE_THAT(cell.a[0], Catch::Matchers::WithinAbs(10.0, 1e-12));
  REQUIRE_THAT(cell.b[0], Catch::Matchers::WithinAbs(5.0, 1e-12));
}

TEST_CASE("residentFrameCell stays ortho for three lengths", "[neighbours]") {
  const double box[3] = {10.0, 11.0, 12.0};
  linkcell::Cell cell = linkcell::Cell::ortho(1.0, 1.0, 1.0);
  const double *frameLens = nullptr;
  nneigh::residentFrameCell(box, nullptr, 3, cell, frameLens);
  REQUIRE(frameLens == box);
  REQUIRE_THAT(cell.a[0], Catch::Matchers::WithinAbs(10.0, 1e-12));
  REQUIRE_THAT(cell.b[1], Catch::Matchers::WithinAbs(11.0, 1e-12));
  REQUIRE_THAT(cell.c[2], Catch::Matchers::WithinAbs(12.0, 1e-12));
}
#endif

namespace {

// A jittered simple-cubic lattice of `side`^3 atoms at spacing `a`, in an
// orthorhombic box or a sheared one (tilts xy, xz, yz in Angstrom).
molSys::PointCloud<molSys::Point<double>, double>
jitteredLattice(int side, double a, double xy, double xz, double yz) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double L = side * a;
  const bool tilted = xy != 0.0 || xz != 0.0 || yz != 0.0;
  if (tilted) {
    // LAMMPS bound spans: lo/hi bounds widen by the tilt extents
    const double xmin = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
    const double xmax = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
    cloud.box = {L + xmax - xmin, L + std::max(0.0, yz) - std::min(0.0, yz), L, xy, xz, yz};
    cloud.boxLow = {xmin, std::min(0.0, yz), 0.0};
  } else {
    cloud.box = {L, L, L};
    cloud.boxLow = {0.0, 0.0, 0.0};
  }
  unsigned long long state = 0x9E3779B97F4A7C15ULL;
  int id = 0;
  for (int i = 0; i < side; i++) {
    for (int j = 0; j < side; j++) {
      for (int k = 0; k < side; k++) {
        double f[3];
        for (double &fk : f) {
          state ^= state << 13;
          state ^= state >> 7;
          state ^= state << 17;
          fk = (static_cast<double>(state >> 11) / 9007199254740992.0 - 0.5) * 0.3 / side;
        }
        const double fa = (i + 0.5) / side + f[0];
        const double fb = (j + 0.5) / side + f[1];
        const double fc = (k + 0.5) / side + f[2];
        molSys::Point<double> p;
        p.type = 1;
        p.atomID = id + 1;
        p.molID = id + 1;
        // r = fa a + fb b + fc c with a = (L,0,0), b = (xy,L,0), c = (xz,yz,L)
        p.x = fa * L + fb * xy + fc * xz;
        p.y = fb * L + fc * yz;
        p.z = fc * L;
        cloud.pts.push_back(p);
        cloud.idIndexMap[p.atomID] = id;
        ++id;
      }
    }
  }
  cloud.nop = id;
  return cloud;
}

} // namespace

TEST_CASE("threaded cell-list rows are the minimum-image neighbours", "[neighbours]") {
  const double cutoff = 3.5;
#ifdef SEAMS_HAS_OPENMP
  // 4096 atoms and more than one thread: the rows run in parallel.
  const int threads = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  for (const auto tilt : {std::array<double, 3>{0.0, 0.0, 0.0}, std::array<double, 3>{4.0, 2.0, 3.0}}) {
    const auto cloud = jitteredLattice(16, 3.0, tilt[0], tilt[1], tilt[2]);
    REQUIRE(cloud.nop == 16 * 16 * 16);
    std::vector<int> all(static_cast<std::size_t>(cloud.nop));
    std::iota(all.begin(), all.end(), 0);
    std::vector<std::vector<int>> rows;
    REQUIRE(nneigh::cellListRowsThreaded(cloud, all, cutoff, rows));
    REQUIRE(rows.size() == all.size());
    for (int iatom = 0; iatom < cloud.nop; iatom += 37) {
      INFO("tilt " << tilt[0] << " atom " << iatom);
      REQUIRE(rows[static_cast<std::size_t>(iatom)] == minimumImageNeighbours(cloud, iatom, cutoff));
    }
    // the public builders agree with the brute-force reference row by row
    const auto byIndex = nneigh::getNewNeighbourListByIndex(cloud, cutoff);
    const auto byID = nneigh::neighbourListByIndex(cloud, nneigh::neighListO(cutoff, cloud, 1));
    REQUIRE(byIndex.size() == all.size());
    REQUIRE(byID.size() == all.size());
    for (int iatom = 0; iatom < cloud.nop; iatom += 53) {
      std::vector<int> a(byIndex[static_cast<std::size_t>(iatom)].begin() + 1,
                         byIndex[static_cast<std::size_t>(iatom)].end());
      std::vector<int> b(byID[static_cast<std::size_t>(iatom)].begin() + 1,
                         byID[static_cast<std::size_t>(iatom)].end());
      std::sort(a.begin(), a.end());
      std::sort(b.begin(), b.end());
      REQUIRE(a == minimumImageNeighbours(cloud, iatom, cutoff));
      REQUIRE(b == a);
    }
  }
  // a cell with fewer than three cutoff widths on an axis refuses, and the
  // callers fall back to the other builders
  const auto tight = jitteredLattice(3, 3.0, 0.0, 0.0, 0.0);
  std::vector<int> all(static_cast<std::size_t>(tight.nop));
  std::iota(all.begin(), all.end(), 0);
  std::vector<std::vector<int>> rows;
  REQUIRE_FALSE(nneigh::cellListRowsThreaded(tight, all, cutoff, rows));
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(threads);
#endif
}

TEST_CASE("index rows list self then the minimum image in ascending order",
          "[neighbours]") {
#ifdef SEAMS_HAS_OPENMP
  const int threads = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  const auto rowsMatch = [](const auto &cloud, double cutoff, int step) {
    const auto byIndex = nneigh::getNewNeighbourListByIndex(cloud, cutoff);
    const auto byID =
        nneigh::neighbourListByIndex(cloud, nneigh::neighListO(cutoff, cloud, 1));
    REQUIRE(byIndex.size() == static_cast<std::size_t>(cloud.nop));
    REQUIRE(byID.size() == static_cast<std::size_t>(cloud.nop));
    for (int iatom = 0; iatom < cloud.nop; iatom += step) {
      std::vector<int> want{iatom};
      const auto ref = minimumImageNeighbours(cloud, iatom, cutoff);
      want.insert(want.end(), ref.begin(), ref.end());
      INFO("cutoff " << cutoff << " atom " << iatom);
      REQUIRE(byIndex[static_cast<std::size_t>(iatom)] == want);
      REQUIRE(byID[static_cast<std::size_t>(iatom)] == want);
    }
  };
  for (const auto tilt : {std::array<double, 3>{0.0, 0.0, 0.0}, std::array<double, 3>{4.0, 2.0, 3.0}}) {
    rowsMatch(jitteredLattice(16, 3.0, tilt[0], tilt[1], tilt[2]), 3.5, 7);
  }
  // Just under and just over half the narrowest face separation, where a
  // pair can sit inside the cutoff through more than one image
  const auto small = jitteredLattice(4, 3.0, 3.0, 1.0, -2.0);
  const double half = gen::makeFracBox(small).halfMin;
  for (const double cutoff : {half - 0.05, half + 0.05}) {
    rowsMatch(small, cutoff, 1);
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(threads);
#endif
}

#ifdef SEAMS_HAS_OPENMP
TEST_CASE("neighListO on four threads leaves an unmapped atom as an empty row",
          "[neighbours]") {
  auto cloud = jitteredLattice(16, 3.0, 0.0, 0.0, 0.0);
  REQUIRE(cloud.nop == 4096);
  cloud.idIndexMap.erase(101);
  const int threads = omp_get_max_threads();
  omp_set_num_threads(4);
  const auto nList = nneigh::neighListO(3.5, cloud, 1);
  omp_set_num_threads(threads);
  REQUIRE(nList[100].empty());
  bool unmappedPartner = false;
  for (const auto &row : nList) {
    unmappedPartner |= std::find(row.begin(), row.end(), -1) != row.end();
  }
  REQUIRE_FALSE(unmappedPartner);
}

TEST_CASE("type-filtered cutoff rows agree on uneven thread counts",
          "[neighbours]") {
  // 4913 rows split unevenly into blocks on three and four threads, and
  // alternating types make every type subset differ from the cloud order
  const double cutoff = 3.5;
  const int threads = omp_get_max_threads();
  for (const auto tilt : {std::array<double, 3>{0.0, 0.0, 0.0}, std::array<double, 3>{4.0, 2.0, 3.0}}) {
    auto cloud = jitteredLattice(17, 3.0, tilt[0], tilt[1], tilt[2]);
    for (auto &p : cloud.pts) {
      p.type = 1 + p.atomID % 2;
    }
    std::vector<std::array<std::vector<std::vector<int>>, 4>> lists;
    for (const int width : {1, 3, 4}) {
      omp_set_num_threads(width);
      lists.push_back({nneigh::neighListO(cutoff, cloud, 2),
                       nneigh::halfNeighList(cutoff, cloud, 1),
                       nneigh::neighListPair(cutoff, cloud, 1, 2),
                       nneigh::getNewNeighbourListByIndex(cloud, cutoff)});
    }
    omp_set_num_threads(threads);
    INFO("tilt " << tilt[0]);
    REQUIRE((lists[1] == lists[0]));
    REQUIRE((lists[2] == lists[0]));
    // the type-2 rows on four threads are the type-2 minimum image, ascending
    const auto byIndex = nneigh::neighbourListByIndex(cloud, lists[2][0]);
    for (int iatom = 0; iatom < cloud.nop; iatom += 7) {
      if (cloud.pts[static_cast<std::size_t>(iatom)].type != 2) {
        continue;
      }
      std::vector<int> want{iatom};
      for (const int j : minimumImageNeighbours(cloud, iatom, cutoff)) {
        if (cloud.pts[static_cast<std::size_t>(j)].type == 2) {
          want.push_back(j);
        }
      }
      INFO("atom " << iatom);
      REQUIRE(byIndex[static_cast<std::size_t>(iatom)] == want);
    }
  }
}
#endif

TEST_CASE("water-type mask keeps Ag out of the 4-NN list", "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  // Five waters in a plus, Ag sitting closer to the centre than one water.
  const double coords[6][3] = {{0, 0, 0},
                               {3, 0, 0},
                               {-3, 0, 0},
                               {0, 3, 0},
                               {0, -3, 0},
                               {1.0, 0, 0}};
  const int types[6] = {1, 1, 1, 1, 1, 3};
  for (int i = 0; i < 6; i++) {
    molSys::Point<double> pt;
    pt.type = types[i];
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }
  cloud.nop = 6;
  const auto nList =
      nneigh::kNearestNeighbourList(cloud, 4, 6.0, std::vector<int>{1}, true);
  REQUIRE(nList.size() == static_cast<std::size_t>(cloud.nop));
  for (const auto &row : nList) {
    for (std::size_t m = 1; m < row.size(); m++) {
      const auto it = cloud.idIndexMap.find(row[m]);
      REQUIRE(it != cloud.idIndexMap.end());
      REQUIRE(cloud.pts[static_cast<std::size_t>(it->second)].type == 1);
    }
  }
  // Ag has no water-graph row bonds
  REQUIRE(nList[5].size() <= 1);
  const auto both =
      nneigh::kNearestNeighbourPair(cloud, 4, 6.0, std::vector<int>{1});
  REQUIRE(both.first[5].size() <= 1);
  REQUIRE(both.second[5].size() <= 1);
  for (const auto &row : both.first) {
    for (std::size_t m = 1; m < row.size(); m++) {
      const auto it = cloud.idIndexMap.find(row[m]);
      REQUIRE(it != cloud.idIndexMap.end());
      REQUIRE(cloud.pts[static_cast<std::size_t>(it->second)].type == 1);
    }
  }
}

TEST_CASE("water-type set bonds two water types and drops Ag", "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  // Type 1 and type 2 are both water; type 3 is Ag closer than one water.
  const double coords[6][3] = {{0, 0, 0},
                               {3, 0, 0},
                               {-3, 0, 0},
                               {0, 3, 0},
                               {0, -3, 0},
                               {0.4, 0, 0}};
  const int types[6] = {1, 2, 1, 2, 1, 3};
  for (int i = 0; i < 6; i++) {
    molSys::Point<double> pt;
    pt.type = types[i];
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }
  cloud.nop = 6;
  const std::vector<int> water{1, 2};
  const auto nList =
      nneigh::kNearestNeighbourList(cloud, 4, 6.0, water, true);
  REQUIRE(nList.size() == static_cast<std::size_t>(cloud.nop));
  bool sawType2 = false;
  for (const auto &row : nList) {
    for (std::size_t m = 1; m < row.size(); m++) {
      const auto it = cloud.idIndexMap.find(row[m]);
      REQUIRE(it != cloud.idIndexMap.end());
      const int t = cloud.pts[static_cast<std::size_t>(it->second)].type;
      REQUIRE((t == 1 || t == 2));
      REQUIRE(t != 3);
      if (t == 2) {
        sawType2 = true;
      }
    }
  }
  REQUIRE(sawType2);
  REQUIRE(nList[5].size() <= 1);
  const auto empty =
      nneigh::kNearestNeighbourList(cloud, 4, 6.0, std::vector<int>{}, true);
  REQUIRE(empty.empty());
  const auto pair = nneigh::kNearestNeighbourPair(cloud, 4, 6.0, water);
  REQUIRE(pair.first.size() == static_cast<std::size_t>(cloud.nop));
  for (const auto &row : pair.first) {
    for (std::size_t m = 1; m < row.size(); m++) {
      const auto it = cloud.idIndexMap.find(row[m]);
      REQUIRE(it != cloud.idIndexMap.end());
      REQUIRE(cloud.pts[static_cast<std::size_t>(it->second)].type != 3);
    }
  }
}

TEST_CASE("mixed AgI+water dump keeps Ag and I out of the water 4-NN list",
          "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj("traj/agi_water_tiny.lammpstrj", 1, cloud);
  REQUIRE(cloud.nop == 7);
  int nAg = 0;
  int nI = 0;
  int nOw = 0;
  for (const auto &pt : cloud.pts) {
    nOw += pt.type == 1 ? 1 : 0;
    nAg += pt.type == 3 ? 1 : 0;
    nI += pt.type == 4 ? 1 : 0;
  }
  REQUIRE(nOw == 5);
  REQUIRE(nAg == 1);
  REQUIRE(nI == 1);
  const auto nList =
      nneigh::kNearestNeighbourList(cloud, 4, 6.0, std::vector<int>{1}, true);
  REQUIRE(nList.size() == static_cast<std::size_t>(cloud.nop));
  for (const auto &row : nList) {
    for (std::size_t m = 1; m < row.size(); m++) {
      const auto it = cloud.idIndexMap.find(row[m]);
      REQUIRE(it != cloud.idIndexMap.end());
      const int t = cloud.pts[static_cast<std::size_t>(it->second)].type;
      REQUIRE(t == 1);
      REQUIRE(t != 3);
      REQUIRE(t != 4);
    }
  }
  int agIdx = -1;
  int iIdx = -1;
  for (int i = 0; i < cloud.nop; i++) {
    if (cloud.pts[static_cast<std::size_t>(i)].type == 3) {
      agIdx = i;
    }
    if (cloud.pts[static_cast<std::size_t>(i)].type == 4) {
      iIdx = i;
    }
  }
  REQUIRE(agIdx >= 0);
  REQUIRE(iIdx >= 0);
  REQUIRE(nList[static_cast<std::size_t>(agIdx)].size() <= 1);
  REQUIRE(nList[static_cast<std::size_t>(iIdx)].size() <= 1);
}

TEST_CASE("a non-positive box does not return a neighbour list", "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = 2;
  cloud.box = {0, 0, 0};
  cloud.boxLow = {0, 0, 0};
  cloud.pts.resize(2);
  cloud.pts[0].type = 1;
  cloud.pts[0].atomID = 1;
  cloud.pts[1].type = 1;
  cloud.pts[1].atomID = 2;
  cloud.pts[1].x = 1;
  cloud.idIndexMap[1] = 0;
  cloud.idIndexMap[2] = 1;
  const auto nList = nneigh::neighListO(3.5, cloud, 1);
  REQUIRE(nList.empty());
}

TEST_CASE("a pair exactly at the cutoff is not a neighbour", "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = 2;
  cloud.box = {20, 20, 20};
  cloud.boxLow = {0, 0, 0};
  cloud.pts.resize(2);
  cloud.pts[0].type = 1;
  cloud.pts[0].atomID = 1;
  cloud.pts[1].type = 1;
  cloud.pts[1].atomID = 2;
  cloud.pts[1].x = 4.0;
  cloud.idIndexMap[1] = 0;
  cloud.idIndexMap[2] = 1;
  const auto on = nneigh::neighListO(4.0, cloud, 1);
  REQUIRE(on[0].size() == 1);
  const auto inside = nneigh::neighListO(4.0 + 1e-6, cloud, 1);
  REQUIRE(inside[0].size() == 2);
}

namespace {

molSys::PointCloud<molSys::Point<double>, double>
cloudFrom(const std::vector<std::array<double, 3>> &pos, double box) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {box, box, box};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  cloud.nop = static_cast<int>(pos.size());
  for (int i = 0; i < cloud.nop; ++i) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i;
    p.molID = i;
    p.x = pos[static_cast<std::size_t>(i)][0];
    p.y = pos[static_cast<std::size_t>(i)][1];
    p.z = pos[static_cast<std::size_t>(i)][2];
    cloud.pts.push_back(p);
    cloud.idIndexMap[i] = i;
  }
  return cloud;
}

std::vector<std::vector<int>> brutePairs(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud, double cutoff) {
  const double rc2 = cutoff * cutoff;
  const gen::FracBox frame = gen::makeFracBox(cloud);
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    rows[static_cast<std::size_t>(i)].push_back(i);
    for (int j = 0; j < cloud.nop; ++j) {
      if (i == j) {
        continue;
      }
      if (gen::periodicDistSq(frame, cloud, i, j) < rc2) {
        rows[static_cast<std::size_t>(i)].push_back(j);
      }
    }
    std::sort(rows[static_cast<std::size_t>(i)].begin() + 1,
              rows[static_cast<std::size_t>(i)].end());
  }
  return rows;
}

std::vector<std::array<double, 3>> iceIc(int n, double a) {
  const double b[8][3] = {
      {0, 0, 0},         {0, 0.5, 0.5},     {0.5, 0, 0.5},     {0.5, 0.5, 0},
      {0.25, 0.25, 0.25},{0.25, 0.75, 0.75},{0.75, 0.25, 0.75},{0.75, 0.75, 0.25}};
  std::vector<std::array<double, 3>> pos;
  pos.reserve(static_cast<std::size_t>(n * n * n * 8));
  for (int ix = 0; ix < n; ++ix) {
    for (int iy = 0; iy < n; ++iy) {
      for (int iz = 0; iz < n; ++iz) {
        for (int k = 0; k < 8; ++k) {
          pos.push_back({(ix + b[k][0]) * a, (iy + b[k][1]) * a, (iz + b[k][2]) * a});
        }
      }
    }
  }
  return pos;
}

} // namespace

TEST_CASE("cutoff pairs match the minimum image at every thread count",
          "[neighbours]") {
  const double cutoff = 3.5;
  const auto ice = iceIc(4, 6.35);
  auto crystal = cloudFrom(ice, 4 * 6.35);
  const auto expect = brutePairs(crystal, cutoff);
  for (int i = 0; i < crystal.nop; ++i) {
    REQUIRE(expect[static_cast<std::size_t>(i)].size() == 5);
  }

  std::vector<std::array<double, 3>> tight = {
      {0.1, 0.1, 0.1}, {7.9, 0.2, 0.2}, {0.1, 4.0, 0.1}, {3.0, 3.0, 3.0}};
  auto edge = cloudFrom(tight, 8.0);
  const auto edgeExpect = brutePairs(edge, 3.9);

  std::vector<std::array<double, 3>> grid;
  const int side = 15;
  for (int ix = 0; ix < side; ++ix) {
    for (int iy = 0; iy < side; ++iy) {
      for (int iz = 0; iz < side; ++iz) {
        grid.push_back({10.0 * ix, 10.0 * iy, 10.0 * iz});
      }
    }
  }
  grid[1] = {3.5, 0.0, 0.0};
  auto boundary = cloudFrom(grid, 10.0 * side);
  const auto boundaryExpect = brutePairs(boundary, 3.5);
  REQUIRE(boundaryExpect[0].size() == 1);

#ifdef SEAMS_HAS_OPENMP
  const int saved = omp_get_max_threads();
  const int widths[] = {1, 2, 4, 8};
#else
  const int widths[] = {1};
#endif
  for (int width : widths) {
#ifdef SEAMS_HAS_OPENMP
    omp_set_num_threads(width);
#endif
    INFO("threads " << width);
    REQUIRE(nneigh::neighListO(cutoff, crystal, 1) == expect);
    REQUIRE(nneigh::getNewNeighbourListByIndex(edge, 3.9) == edgeExpect);
    REQUIRE(nneigh::neighListO(3.5, boundary, 1) == boundaryExpect);
    auto half = nneigh::halfNeighList(cutoff, crystal, 1);
    for (int i = 0; i < crystal.nop; ++i) {
      for (std::size_t k = 1; k < half[static_cast<std::size_t>(i)].size(); ++k) {
        REQUIRE(half[static_cast<std::size_t>(i)][k] > i);
      }
    }
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(saved);
#endif
}

#ifdef SEAMS_HAS_VESIN
TEST_CASE("cutoff pairs match vesin on ice Ic", "[neighbours]") {
  const auto ice = iceIc(3, 6.35);
  auto cloud = cloudFrom(ice, 3 * 6.35);
  const double cutoff = 3.5;
  std::vector<std::array<double, 3>> xyz(ice.begin(), ice.end());
  double box[3][3] = {};
  box[0][0] = box[1][1] = box[2][2] = 3 * 6.35;
  const bool periodic[3] = {true, true, true};
  VesinOptions options{};
  options.cutoff = cutoff;
  options.full = true;
  options.sorted = false;
  options.algorithm = VesinCellList;
  options.return_shifts = false;
  options.return_distances = false;
  options.return_vectors = false;
  VesinNeighborList neighbors{};
  const char *error = nullptr;
  const VesinDevice device{VesinCPU, 0};
  REQUIRE(vesin_neighbors(reinterpret_cast<const double (*)[3]>(xyz.data()), xyz.size(),
                          box, periodic, device, options, &neighbors, &error) == 0);
  std::vector<std::pair<int, int>> got;
  got.reserve(neighbors.length);
  for (std::size_t k = 0; k < neighbors.length; ++k) {
    got.emplace_back(static_cast<int>(neighbors.pairs[k][0]),
                     static_cast<int>(neighbors.pairs[k][1]));
  }
  vesin_free(&neighbors);
  std::sort(got.begin(), got.end());
  got.erase(std::unique(got.begin(), got.end()), got.end());

  const auto rows = nneigh::neighListO(cutoff, cloud, 1);
  std::vector<std::pair<int, int>> ours;
  for (int i = 0; i < cloud.nop; ++i) {
    for (std::size_t k = 1; k < rows[static_cast<std::size_t>(i)].size(); ++k) {
      ours.emplace_back(i, rows[static_cast<std::size_t>(i)][k]);
    }
  }
  std::sort(ours.begin(), ours.end());
  REQUIRE(ours == got);
  REQUIRE(ours.size() == static_cast<std::size_t>(4 * cloud.nop));
}
#endif

TEST_CASE("triclinic cutoff pairs match the minimum image", "[neighbours]") {
  const double L = 4 * 6.35;
  const double xy = 0.25 * L;
  const double xz = 0.05 * L;
  const double yz = 0.12 * L;
  const double xmin = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
  const double xmax = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
  const double ymin = std::min(0.0, yz);
  const double ymax = std::max(0.0, yz);
  auto frac = iceIc(4, 1.0);
  for (auto &s : frac) {
    s[0] /= 4.0;
    s[1] /= 4.0;
    s[2] /= 4.0;
  }
  std::vector<std::array<double, 3>> cart;
  cart.reserve(frac.size());
  for (const auto &s : frac) {
    cart.push_back({xmin + L * s[0] + xy * s[1] + xz * s[2],
                    ymin + L * s[1] + yz * s[2], s[2] * L});
  }
  auto cloud = cloudFrom(cart, 1.0);
  cloud.box = {L + xmax - xmin, L + ymax - ymin, L, xy, xz, yz};
  cloud.boxLow = {xmin, ymin, 0.0};
  const double cutoff = 3.5;
  const auto expect = brutePairs(cloud, cutoff);
#ifdef SEAMS_HAS_OPENMP
  const int saved = omp_get_max_threads();
  for (int width : {1, 2, 8}) {
    omp_set_num_threads(width);
    REQUIRE(nneigh::neighListO(cutoff, cloud, 1) == expect);
  }
  omp_set_num_threads(saved);
#else
  REQUIRE(nneigh::neighListO(cutoff, cloud, 1) == expect);
#endif
}

TEST_CASE("half cutoff list is the lower index of the full list", "[neighbours]") {
  auto cloud = cloudFrom(iceIc(3, 6.35), 3 * 6.35);
  std::vector<int> fullOff, fullCols, halfOff, halfCols;
  REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, 3.5, fullOff, fullCols, false));
  REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, 3.5, halfOff, halfCols, true));
  std::vector<std::pair<int, int>> full;
  std::vector<std::pair<int, int>> half;
  for (int i = 0; i < cloud.nop; ++i) {
    for (int k = fullOff[static_cast<std::size_t>(i)];
         k < fullOff[static_cast<std::size_t>(i) + 1]; ++k) {
      const int j = fullCols[static_cast<std::size_t>(k)];
      if (i < j) {
        full.emplace_back(i, j);
      }
    }
    for (int k = halfOff[static_cast<std::size_t>(i)];
         k < halfOff[static_cast<std::size_t>(i) + 1]; ++k) {
      half.emplace_back(i, halfCols[static_cast<std::size_t>(k)]);
    }
  }
  REQUIRE(half == full);
  REQUIRE(half.size() == static_cast<std::size_t>(2 * cloud.nop));
}

TEST_CASE("a skin rebuild matches a fresh cutoff list", "[neighbours]") {
  auto cloud = cloudFrom(iceIc(3, 6.35), 3 * 6.35);
  cloud.pts[0].x = 0.05;
  cloud.pts[1].x = 3.40;
  cloud.pts[1].y = 0.0;
  cloud.pts[1].z = 0.0;
  nneigh::CutoffSkin skin(3.5, 0.8);
  REQUIRE(skin.update(cloud));
  REQUIRE_FALSE(skin.lastReused());
  cloud.pts[1].x = 3.55;
  REQUIRE(skin.update(cloud));
  REQUIRE(skin.lastReused());
  std::vector<int> freshOff, freshCols;
  REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, 3.5, freshOff, freshCols, false));
  REQUIRE(skin.offsets() == freshOff);
  REQUIRE(skin.cols() == freshCols);
  cloud.pts[1].x += 1.0;
  REQUIRE(skin.update(cloud));
  REQUIRE_FALSE(skin.lastReused());
  REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, 3.5, freshOff, freshCols, false));
  REQUIRE(skin.offsets() == freshOff);
  REQUIRE(skin.cols() == freshCols);
}

namespace {

molSys::PointCloud<molSys::Point<double>, double> shearedIce(int cells) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double bond = 2.75;
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  const double side = cells * lattice;
  const double xy = 4.0;
  const double xz = 1.5;
  const double yz = 0.8;
  const double basis[8][3] = {
      {0, 0, 0},          {0.5, 0.5, 0},         {0.5, 0, 0.5},         {0, 0.5, 0.5},
      {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},    {0.75, 0.25, 0.75},    {0.25, 0.75, 0.75}};
  int id = 1;
  for (int i = 0; i < cells; ++i) {
    for (int j = 0; j < cells; ++j) {
      for (int k = 0; k < cells; ++k) {
        for (const auto &b : basis) {
          const double X = (i + b[0]) * lattice;
          const double Y = (j + b[1]) * lattice;
          const double Z = (k + b[2]) * lattice;
          molSys::Point<double> p;
          p.x = X + (xy / side) * Y + (xz / side) * Z;
          p.y = Y + (yz / side) * Z;
          p.z = Z;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          cloud.idIndexMap[id] = cloud.nop;
          cloud.pts.push_back(p);
          ++cloud.nop;
          ++id;
        }
      }
    }
  }
  const double xmin = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
  const double xmax = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
  const double ymin = std::min(0.0, yz);
  const double ymax = std::max(0.0, yz);
  cloud.box = {side + xmax - xmin, side + ymax - ymin, side, xy, xz, yz};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

std::uint64_t csrDigest(const std::vector<int> &offsets, const std::vector<int> &cols) {
  std::uint64_t h = 1469598103934665603ull;
  auto mix = [&](std::uint64_t v) {
    h ^= v;
    h *= 1099511628211ull;
  };
  mix(offsets.size());
  for (int v : offsets) {
    mix(static_cast<std::uint64_t>(v));
  }
  mix(cols.size());
  for (int v : cols) {
    mix(static_cast<std::uint64_t>(v));
  }
  return h;
}

} // namespace

TEST_CASE("sheared ice keeps the shared shift-table pairs", "[neighbours][triclinic]") {
  std::vector<int> offsets, cols;
  auto small = shearedIce(24);
  REQUIRE(nneigh::orthoFullCsr(small, nullptr, 3.5, offsets, cols));
  REQUIRE(small.nop == 110592);
  REQUIRE(cols.size() == 442368);
  REQUIRE(csrDigest(offsets, cols) == 0x647177f805d1e1c6ull);
  auto large = shearedIce(48);
  REQUIRE(nneigh::orthoFullCsr(large, nullptr, 3.5, offsets, cols));
  REQUIRE(large.nop == 884736);
  REQUIRE(cols.size() == 3538944);
  REQUIRE(csrDigest(offsets, cols) == 0xba88eb7da3372fc6ull);
}

TEST_CASE("a shell longer than the stash matches the brute-force cutoff", "[neighbours]") {
#ifdef SEAMS_HAS_OPENMP
  const int saved = omp_get_max_threads();
  omp_set_num_threads(saved < 8 ? saved : 8);
#endif
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {8.0, 8.0, 8.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  const int nSide = 4;
  const int n = nSide * nSide * nSide;
  for (int iz = 0; iz < nSide; ++iz) {
    for (int iy = 0; iy < nSide; ++iy) {
      for (int ix = 0; ix < nSide; ++ix) {
        molSys::Point<double> pt;
        pt.type = 1;
        pt.atomID = static_cast<int>(cloud.pts.size());
        pt.molID = pt.atomID;
        pt.x = ix * 2.0;
        pt.y = iy * 2.0;
        pt.z = iz * 2.0;
        cloud.idIndexMap[pt.atomID] = pt.atomID;
        cloud.pts.push_back(pt);
      }
    }
  }
  cloud.nop = n;
  const double cutoff = 6.0;
  const double rc2 = cutoff * cutoff;
  const double edge = 8.0;
  std::vector<int> offsets;
  std::vector<int> cols;
  REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, cutoff, offsets, cols));
  REQUIRE(offsets.size() == static_cast<std::size_t>(n + 1));
  for (int i = 0; i < n; ++i) {
    int brute = 0;
    for (int j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      const double dx = cloud.pts[static_cast<std::size_t>(j)].x - cloud.pts[static_cast<std::size_t>(i)].x;
      const double dy = cloud.pts[static_cast<std::size_t>(j)].y - cloud.pts[static_cast<std::size_t>(i)].y;
      const double dz = cloud.pts[static_cast<std::size_t>(j)].z - cloud.pts[static_cast<std::size_t>(i)].z;
      const double wx = dx - std::nearbyint(dx / edge) * edge;
      const double wy = dy - std::nearbyint(dy / edge) * edge;
      const double wz = dz - std::nearbyint(dz / edge) * edge;
      if (wx * wx + wy * wy + wz * wz < rc2) {
        ++brute;
      }
    }
    REQUIRE(brute > 48);
    REQUIRE(offsets[static_cast<std::size_t>(i + 1)] - offsets[static_cast<std::size_t>(i)] == brute);
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(saved);
#endif
}
