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

TEST_CASE("shellSeparation fails when a site has fewer than k others",
          "[neighbours]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  const double xs[4] = {1.0, 2.0, 3.0, 4.0};
  for (int i = 0; i < 4; i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i + 1;
    p.x = xs[i];
    p.y = 5.0;
    p.z = 5.0;
    cloud.pts.push_back(p);
  }
  cloud.nop = 4;
  // Three others per site: k = 4 has no k-th neighbour.
  const auto [maxKth, minNext] = nneigh::shellSeparation(cloud, 4, 1);
  REQUIRE(std::isinf(maxKth));
  REQUIRE(minNext == 0.0);
  // k = 3 is exactly the others: the k-th distance is 3 (from an end site),
  // and with no (k+1)-th distance there is no upper bound either.
  const auto [maxThird, nextAfter] = nneigh::shellSeparation(cloud, 3, 1);
  REQUIRE_THAT(maxThird, Catch::Matchers::WithinAbs(3.0, 1e-12));
  REQUIRE(nextAfter == 0.0);
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
  yCloud = sinp::readLammpsTrjO("traj/cubic-real-7.fixture.lammpstrj", 1, yCloud, 1);
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
  yCloud = sinp::readLammpsTrjO("traj/cubic-real-7.fixture.lammpstrj", 1, yCloud, 1);
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
  yCloud = sinp::readLammpsTrjO("traj/cubic-real-7.fixture.lammpstrj", 1, yCloud, 1);

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
  yCloud = sinp::readLammpsTrjO("traj/cubic-real-7.fixture.lammpstrj", 1, yCloud, 1);
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
  yCloud = sinp::readLammpsTrjO("traj/cubic-real-7.fixture.lammpstrj", 1, yCloud, 1);
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
  // Five waters in a plus, Ag sitting closer to the center than one water.
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

molSys::PointCloud<molSys::Point<double>, double> hexIce(int ncell) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double bond = 2.76;
  const double a = bond * std::sqrt(8.0 / 3.0);
  const double cLen = a * std::sqrt(8.0 / 3.0);
  const double b = a * std::sqrt(3.0);
  const double basis[4][3] = {
      {1.0 / 3.0, 2.0 / 3.0, 0.0},
      {2.0 / 3.0, 1.0 / 3.0, 0.5},
      {1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0},
      {2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0},
  };
  const double lx = ncell * a;
  const double ly = ncell * b;
  const double lz = ncell * cLen;
  int id = 1;
  for (int ix = 0; ix < ncell; ++ix) {
    for (int iy = 0; iy < ncell; ++iy) {
      for (int iz = 0; iz < ncell; ++iz) {
        for (double shift : {0.0, 0.5}) {
          for (const auto &site : basis) {
            molSys::Point<double> p;
            const double x = shift * a + site[0] * a - site[1] * a / 2.0 + ix * a;
            const double y = shift * b + site[1] * b / 2.0 + iy * b;
            p.x = x - lx * std::floor(x / lx);
            p.y = y - ly * std::floor(y / ly);
            p.z = (site[2] + iz) * cLen;
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
  }
  cloud.box = {lx, ly, lz};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

std::vector<std::pair<int, int>> bruteOrthoPairs(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud, double cutoff) {
  const double rc2 = cutoff * cutoff;
  const double lx = cloud.box[0];
  const double ly = cloud.box[1];
  const double lz = cloud.box[2];
  std::vector<std::pair<int, int>> pairs;
  for (int i = 0; i < cloud.nop; ++i) {
    const auto &a = cloud.pts[static_cast<std::size_t>(i)];
    for (int j = i + 1; j < cloud.nop; ++j) {
      const auto &b = cloud.pts[static_cast<std::size_t>(j)];
      double dx = a.x - b.x;
      double dy = a.y - b.y;
      double dz = a.z - b.z;
      dx -= lx * std::nearbyint(dx / lx);
      dy -= ly * std::nearbyint(dy / ly);
      dz -= lz * std::nearbyint(dz / lz);
      if (dx * dx + dy * dy + dz * dz < rc2) {
        pairs.emplace_back(i, j);
      }
    }
  }
  return pairs;
}

std::vector<std::pair<int, int>> csrPairs(const std::vector<int> &offsets,
                                         const std::vector<int> &cols, bool lowerOnly) {
  std::vector<std::pair<int, int>> pairs;
  const int n = static_cast<int>(offsets.size()) - 1;
  for (int i = 0; i < n; ++i) {
    for (int k = offsets[static_cast<std::size_t>(i)];
         k < offsets[static_cast<std::size_t>(i) + 1]; ++k) {
      const int j = cols[static_cast<std::size_t>(k)];
      if (!lowerOnly || i < j) {
        pairs.emplace_back(i, j);
      }
    }
  }
  if (lowerOnly) {
    std::sort(pairs.begin(), pairs.end());
  }
  return pairs;
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

TEST_CASE("the half-stencil stash gives the pair list at one and eight threads",
          "[neighbours]") {
  // Hexagonal ice has four neighbours at 3.5 A and stays in the eight
  // stash slots; at 6.0 A every shell overflows them and the list is
  // recounted. Both have to be the brute-force pair set, and eight
  // threads have to write the bytes one thread writes.
  const struct {
    int cells;
    double cutoff;
  } cases[] = {{8, 3.5}, {4, 6.0}, {2, 3.5}};
#ifdef SEAMS_HAS_OPENMP
  const int saved = omp_get_max_threads();
#endif
  for (const auto &row : cases) {
    INFO("cells " << row.cells << " cutoff " << row.cutoff);
    auto cloud = hexIce(row.cells);
    std::vector<int> serialOff, serialCols;
#ifdef SEAMS_HAS_OPENMP
    omp_set_num_threads(1);
#endif
    REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, row.cutoff, serialOff, serialCols, false));
    REQUIRE(csrPairs(serialOff, serialCols, true) == bruteOrthoPairs(cloud, row.cutoff));
#ifdef SEAMS_HAS_OPENMP
    omp_set_num_threads(8);
    std::vector<int> wideOff, wideCols;
    REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, row.cutoff, wideOff, wideCols, false));
    REQUIRE(wideOff == serialOff);
    REQUIRE(wideCols == serialCols);
#endif
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(saved);
#endif
}

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

namespace {

std::uint64_t rowsDigest(const std::vector<std::vector<int>> &rows) {
  std::uint64_t h = 1469598103934665603ull;
  auto mix = [&](std::uint64_t v) {
    h ^= v;
    h *= 1099511628211ull;
  };
  mix(rows.size());
  for (const auto &row : rows) {
    mix(row.size());
    for (int v : row) {
      mix(static_cast<std::uint64_t>(v));
    }
  }
  return h;
}

//! Squared distance of the nearest image of j to i, over every image
//! within two cells along the periodic axes. Valid past half the box.
double nearestImageSq(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
                      int i, int j, const std::array<bool, 3> &periodic) {
  double H[3][3] = {};
  double origin[3] = {};
  if (cloud.box.size() >= 3) {
    nneigh::dumpBoundsToH(cloud.box, cloud.boxLow, H, origin);
  }
  const auto &a = cloud.pts[static_cast<std::size_t>(i)];
  const auto &b = cloud.pts[static_cast<std::size_t>(j)];
  const int reach[3] = {periodic[0] ? 2 : 0, periodic[1] ? 2 : 0, periodic[2] ? 2 : 0};
  double best = std::numeric_limits<double>::infinity();
  for (int na = -reach[0]; na <= reach[0]; ++na) {
    for (int nb = -reach[1]; nb <= reach[1]; ++nb) {
      for (int nc = -reach[2]; nc <= reach[2]; ++nc) {
        const double dx = b.x - a.x + na * H[0][0] + nb * H[1][0] + nc * H[2][0];
        const double dy = b.y - a.y + nb * H[1][1] + nc * H[2][1];
        const double dz = b.z - a.z + nc * H[2][2];
        best = std::min(best, dx * dx + dy * dy + dz * dz);
      }
    }
  }
  return best;
}

//! Rows of cloud indices, ascending, without self, from nearestImageSq.
std::vector<std::vector<int>> bruteImageRows(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud, double cutoff,
    const std::array<bool, 3> &periodic, bool half) {
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    for (int j = half ? i + 1 : 0; j < cloud.nop; ++j) {
      if (j != i && nearestImageSq(cloud, i, j, periodic) < cutoff * cutoff) {
        rows[static_cast<std::size_t>(i)].push_back(j);
      }
    }
  }
  return rows;
}

std::vector<std::vector<int>> csrRows(const linkcell::PairCsr &csr) {
  std::vector<std::vector<int>> rows(csr.offsets.size() - 1);
  for (std::size_t i = 0; i + 1 < csr.offsets.size(); ++i) {
    rows[i].assign(csr.cols.begin() + csr.offsets[i], csr.cols.begin() + csr.offsets[i + 1]);
  }
  return rows;
}

molSys::PointCloud<molSys::Point<double>, double> randomCloud(int n, double lx, double ly,
                                                             double lz, unsigned seed) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int i = 0; i < n; ++i) {
    molSys::Point<double> p;
    p.type = 1 + (i % 3 == 0 ? 1 : 0);
    p.atomID = i + 1;
    p.molID = i + 1;
    p.x = lx * u(rng);
    p.y = ly * u(rng);
    p.z = lz * u(rng);
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = i;
  }
  cloud.nop = n;
  cloud.currentFrame = 1;
  return cloud;
}

std::vector<int> bruteKnnPacked(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
                                int k) {
  const int n = cloud.nop;
  std::vector<int> packed(static_cast<std::size_t>(n) * static_cast<std::size_t>(k), -1);
  for (int i = 0; i < n; ++i) {
    std::vector<std::pair<double, int>> dist;
    for (int j = 0; j < n; ++j) {
      if (i != j) {
        dist.emplace_back(nearestImageSq(cloud, i, j, {true, true, true}), j);
      }
    }
    std::sort(dist.begin(), dist.end());
    const int take = std::min(k, static_cast<int>(dist.size()));
    for (int t = 0; t < take; ++t) {
      packed[static_cast<std::size_t>(i) * static_cast<std::size_t>(k) +
             static_cast<std::size_t>(t)] = dist[static_cast<std::size_t>(t)].second;
    }
  }
  return packed;
}

std::vector<int> threadWidths() {
#ifdef SEAMS_HAS_OPENMP
  return {1, 2, 8};
#else
  return {1};
#endif
}

void setWidth(int width) {
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(width);
#else
  (void)width;
#endif
}

} // namespace

TEST_CASE("cutoff entry points keep the rows develop built", "[neighbours][linkcell]") {
  // Digests of the rows the in-tree builders gave on develop 13e26c97,
  // the same at one and eight threads. The linkcell rows have to be
  // those bytes: same partners, same order, same self entry.
  struct Want {
    std::uint64_t full;
    std::uint64_t half;
    std::uint64_t index;
  };
  const struct {
    const char *name;
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    Want want;
  } cases[] = {
      {"sheared 4096", shearedIce(8),
       {0xad66965d5892eef9ull, 0xa1e773c82de3faf9ull, 0x9b707739aec69379ull}},
      {"hexagonal 4096", hexIce(8),
       {0x882e581977793a71ull, 0x017f54531508a399ull, 0x2bb91f7456b13219ull}},
      {"sheared 110592", shearedIce(24),
       {0xc1db08d2407deb79ull, 0xffab258f688a4099ull, 0x1697434ab75036f9ull}},
      {"hexagonal 110592", hexIce(24),
       {0x3f2057f2b42a5af1ull, 0xc37fc364ab8c2ed9ull, 0x36f4ee9046abde19ull}},
  };
#ifdef SEAMS_HAS_OPENMP
  const int saved = omp_get_max_threads();
#endif
  for (const auto &c : cases) {
    for (int width : {1, 8}) {
      setWidth(width);
      INFO(c.name << " threads " << width);
      REQUIRE(rowsDigest(nneigh::neighListO(3.5, c.cloud, 1)) == c.want.full);
      REQUIRE(rowsDigest(nneigh::halfNeighList(3.5, c.cloud, 1)) == c.want.half);
      REQUIRE(rowsDigest(nneigh::getNewNeighbourListByIndex(c.cloud, 3.5)) == c.want.index);
      std::vector<std::vector<int>> refilled;
      REQUIRE(nneigh::refillNeighListO(3.5, c.cloud, 1, refilled));
      REQUIRE(rowsDigest(refilled) == c.want.full);
      REQUIRE(nneigh::refillIndexRows(c.cloud, 3.5, refilled));
      REQUIRE(rowsDigest(refilled) == c.want.index);
    }
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(saved);
#endif
}

TEST_CASE("cutoff builders reject an invalid cutoff and non-finite coordinates",
          "[neighbours][linkcell]") {
  auto cloud = makeFourAtomCloud();
  REQUIRE_THROWS_AS(nneigh::neighListO(0.0, cloud, 1), std::invalid_argument);
  REQUIRE_THROWS_AS(nneigh::halfNeighList(-1.0, cloud, 1), std::invalid_argument);
  REQUIRE_THROWS_AS(nneigh::getNewNeighbourListByIndex(
                        cloud, std::numeric_limits<double>::quiet_NaN()),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(
      nneigh::neighListO(std::numeric_limits<double>::infinity(), cloud, 1),
      std::invalid_argument);

  cloud.pts[1].x = std::numeric_limits<double>::infinity();
  REQUIRE_THROWS_AS(nneigh::neighListO(1.5, cloud, 1), linkcell::Error);
  REQUIRE_THROWS_AS(nneigh::halfNeighList(1.5, cloud, 1), linkcell::Error);
  REQUIRE_THROWS_AS(nneigh::getNewNeighbourListByIndex(cloud, 1.5), linkcell::Error);
  REQUIRE_THROWS_AS(nneigh::neighListPair(1.5, cloud, 1, 2), linkcell::Error);
}

TEST_CASE("SkinNeighborList rejects a cloud without a box", "[neighbours][linkcell]") {
  auto cloud = makeFourAtomCloud();
  cloud.box.clear();
  cloud.boxLow.clear();
  nneigh::SkinNeighborList skin(1.5, 0.5, 1, nneigh::BondGraph::Cutoff);
  REQUIRE_THROWS_AS(skin.update(cloud), std::invalid_argument);
}

TEST_CASE("a tiny cutoff in a large box stays under the linkcell mesh limit",
          "[neighbours][linkcell]") {
  // One cell per 0.01 A of a 20 A box is 8e9 cells, past linkcell's limit.
  auto cloud = randomCloud(64, 20.0, 20.0, 20.0, 31);
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.pts[1].x = cloud.pts[0].x + 0.005;
  cloud.pts[1].y = cloud.pts[0].y;
  cloud.pts[1].z = cloud.pts[0].z;
  auto loose = cloud;
  loose.box.clear();
  loose.boxLow.clear();
  const double cutoff = 0.01;
  linkcell::PairCsr csr;
  for (bool half : {false, true}) {
    nneigh::linkcellCutoffCsr(cloud, cutoff, half, {true, true, true}, {}, 1, csr);
    REQUIRE(csrRows(csr) == bruteImageRows(cloud, cutoff, {true, true, true}, half));
    nneigh::linkcellCutoffCsr(loose, cutoff, half, {false, false, false}, {}, 1, csr);
    REQUIRE(csrRows(csr) == bruteImageRows(loose, cutoff, {false, false, false}, half));
  }
  const auto rows = nneigh::neighListPair(cutoff, cloud, 1, 2);
  REQUIRE(rows.size() == static_cast<std::size_t>(cloud.nop));
}

TEST_CASE("a cutoff past half a small triclinic cell keeps every image pair",
          "[neighbours][linkcell][triclinic]") {
  // 3.5 is more than half of every width of this cell, so one pair can
  // meet through several images; the row lists the nearest one once.
  auto cloud = randomCloud(40, 6.2, 6.4, 6.6, 11);
  for (auto &p : cloud.pts) {
    p.type = 1;
  }
  REQUIRE(nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, 0.0, 0.0, 0.0, 6.2, 6.4,
                                        6.6, 1.5, 1.0, -1.1));
  const double cutoff = 3.5;
  const auto full = bruteImageRows(cloud, cutoff, {true, true, true}, false);
  const auto half = bruteImageRows(cloud, cutoff, {true, true, true}, true);
  std::vector<std::vector<int>> wantIndex(full.size());
  std::vector<std::vector<int>> wantFull(full.size());
  std::vector<std::vector<int>> wantHalf(full.size());
  for (std::size_t i = 0; i < full.size(); ++i) {
    wantIndex[i] = {static_cast<int>(i)};
    wantIndex[i].insert(wantIndex[i].end(), full[i].begin(), full[i].end());
    wantFull[i] = {static_cast<int>(i) + 1};
    for (int j : full[i]) {
      wantFull[i].push_back(j + 1);
    }
    wantHalf[i] = {static_cast<int>(i) + 1};
    for (int j : half[i]) {
      wantHalf[i].push_back(j + 1);
    }
  }
#ifdef SEAMS_HAS_OPENMP
  const int saved = omp_get_max_threads();
#endif
  for (int width : threadWidths()) {
    setWidth(width);
    INFO("threads " << width);
    REQUIRE(nneigh::getNewNeighbourListByIndex(cloud, cutoff) == wantIndex);
    REQUIRE(nneigh::neighListO(cutoff, cloud, 1) == wantFull);
    REQUIRE(nneigh::halfNeighList(cutoff, cloud, 1) == wantHalf);
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(saved);
#endif
}

TEST_CASE("linkcell CSR rows match brute force on open and mixed boundaries",
          "[neighbours][linkcell]") {
  const double cutoff = 2.4;
  auto boxed = randomCloud(300, 12.0, 11.0, 10.0, 5);
  boxed.box = {12.0, 11.0, 10.0};
  boxed.boxLow = {0.0, 0.0, 0.0};
  auto tilted = randomCloud(300, 12.0, 11.0, 10.0, 6);
  REQUIRE(nneigh::setBoxFromLammpsEdges(tilted.box, tilted.boxLow, 0.0, 0.0, 0.0, 12.0,
                                        11.0, 10.0, 2.0, -1.0, 1.5));
  auto loose = randomCloud(300, 12.0, 11.0, 10.0, 7);
  const std::array<std::array<bool, 3>, 4> layouts = {{{true, true, true},
                                                       {false, false, false},
                                                       {true, false, true},
                                                       {false, true, false}}};
  linkcell::PairCsr csr;
#ifdef SEAMS_HAS_OPENMP
  const int saved = omp_get_max_threads();
#endif
  for (int width : threadWidths()) {
    setWidth(width);
    for (const auto *cloud : {&boxed, &tilted}) {
      for (const auto &periodic : layouts) {
        INFO("threads " << width << " periodic " << periodic[0] << periodic[1]
                        << periodic[2] << " tilted " << (cloud == &tilted));
        for (bool half : {false, true}) {
          nneigh::linkcellCutoffCsr(*cloud, cutoff, half, periodic, {}, width, csr);
          REQUIRE(csrRows(csr) == bruteImageRows(*cloud, cutoff, periodic, half));
        }
      }
    }
    INFO("threads " << width << " no box");
    nneigh::linkcellCutoffCsr(loose, cutoff, false, {false, false, false}, {}, width, csr);
    REQUIRE(csrRows(csr) == bruteImageRows(loose, cutoff, {false, false, false}, false));
    REQUIRE_THROWS_AS(
        nneigh::linkcellCutoffCsr(loose, cutoff, false, {true, false, false}, {}, width, csr),
        std::invalid_argument);

    std::vector<int> mask(static_cast<std::size_t>(boxed.nop), 0);
    for (int i = 0; i < boxed.nop; ++i) {
      mask[static_cast<std::size_t>(i)] = boxed.pts[static_cast<std::size_t>(i)].type == 1;
    }
    nneigh::linkcellCutoffCsr(boxed, cutoff, false, {true, true, true}, mask, width, csr);
    auto want = bruteImageRows(boxed, cutoff, {true, true, true}, false);
    for (int i = 0; i < boxed.nop; ++i) {
      auto &row = want[static_cast<std::size_t>(i)];
      if (mask[static_cast<std::size_t>(i)] == 0) {
        row.clear();
      }
      row.erase(std::remove_if(row.begin(), row.end(),
                               [&](int j) { return mask[static_cast<std::size_t>(j)] == 0; }),
                row.end());
    }
    REQUIRE(csrRows(csr) == want);
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(saved);
#endif
}

TEST_CASE("linkcell k-nearest nominations match the brute-force minimum image",
          "[neighbours][linkcell][knn]") {
  auto ortho = randomCloud(200, 12.0, 11.0, 10.0, 21);
  ortho.box = {12.0, 11.0, 10.0};
  ortho.boxLow = {0.0, 0.0, 0.0};
  auto tilted = randomCloud(200, 12.0, 11.0, 10.0, 22);
  REQUIRE(nneigh::setBoxFromLammpsEdges(tilted.box, tilted.boxLow, -1.0, 0.5, 0.0, 12.0,
                                        11.0, 10.0, 3.4, -1.2, 2.1));
  for (auto *cloud : {&ortho, &tilted}) {
    for (auto &p : cloud->pts) {
      p.type = 1;
    }
  }
#ifdef SEAMS_HAS_OPENMP
  const int saved = omp_get_max_threads();
#endif
  for (int width : threadWidths()) {
    setWidth(width);
    for (const auto *cloud : {&ortho, &tilted}) {
      for (int k : {1, 4, 12}) {
        INFO("threads " << width << " k " << k << " tilted " << (cloud == &tilted));
        std::vector<double> xyz;
        std::vector<int> got;
        nneigh::kNearestPacked(*cloud, k, 3.5, 1, xyz, got);
        REQUIRE(got == bruteKnnPacked(*cloud, k));
      }
    }
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(saved);
#endif
}
