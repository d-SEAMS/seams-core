#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <array>
#include <vector>

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
  for (int i = 0; i < 4; i++) {
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
      REQUIRE(nList[iatom].size() >= 1);
      REQUIRE(nList[iatom][0] == iatom);

      std::vector<int> got(nList[iatom].begin() + 1, nList[iatom].end());
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
  REQUIRE(nneigh::halfNeighList(3.5, cloud, 1).empty());
  REQUIRE(nneigh::getNewNeighbourListByIndex(cloud, 3.5).empty());
}

TEST_CASE("neighList same type does not self-include or double-count",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  auto nList = nneigh::neighList(1.5, cloud, 1, 1);
  REQUIRE(nList.size() == 4);
  for (int i = 0; i < 4; i++) {
    REQUIRE(nList[i][0] == i);
    auto neigh = nList[i];
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
  nneigh::SkinNeighborList skin(1.5, 0.5, 1);
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
  nneigh::SkinNeighborList skin(1.5, 1.0, 1);
  (void)skin.update(cloud);
  REQUIRE(skin.lastRebuilt());

  cloud.pts[1].x += 0.2; // skin/2 = 0.5
  (void)skin.update(cloud);
  REQUIRE_FALSE(skin.lastRebuilt());

  cloud.pts[1].x += 0.5;
  (void)skin.update(cloud);
  REQUIRE(skin.lastRebuilt());
}

TEST_CASE("SkinNeighborList keeps a bond until cutoff plus skin",
          "[neighbours]") {
  auto cloud = makeFourAtomCloud();
  nneigh::SkinNeighborList skin(1.5, 0.6, 1);
  (void)skin.update(cloud);
  REQUIRE(std::find(skin.bonds()[0].begin() + 1, skin.bonds()[0].end(), 1) !=
          skin.bonds()[0].end());

  cloud.pts[1].x = 1.8; // 1.8 > 1.5, still < 2.1
  (void)skin.update(cloud);
  REQUIRE(std::find(skin.bonds()[0].begin() + 1, skin.bonds()[0].end(), 1) !=
          skin.bonds()[0].end());

  cloud.pts[1].x = 2.3; // 2.3 > 2.1
  (void)skin.update(cloud);
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
