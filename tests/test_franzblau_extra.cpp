#include <catch2/catch_test_macros.hpp>

#include <franzblau.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <set>
#include <utility>
#include <vector>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>

#include <vector>

// Helper: build a small cloud with known atom IDs for populateGraphFromNListID
static molSys::PointCloud<molSys::Point<double>, double> makeSmallCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  // Triangle: 3 atoms mutually within cutoff
  double coords[3][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
  for (int i = 0; i < 3; i++) {
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
  cloud.nop = 3;
  return cloud;
}

// -- populateGraphFromNListID tests --

TEST_CASE("populateGraphFromNListID builds graph from ID-based nList",
          "[franzblau]") {
  auto cloud = makeSmallCloud();
  // Build nList by atom ID (neighListO format: first elem is self ID)
  auto nList = nneigh::neighListO(1.5, cloud, 1);
  REQUIRE(nList.size() == 3);

  auto graph = primitive::populateGraphFromNListID(cloud, nList);

  // Graph should have 3 vertices
  REQUIRE(graph.pts.size() == 3);
  // Each vertex should have the correct atomIndex
  for (int i = 0; i < 3; i++) {
    REQUIRE(graph.pts[static_cast<std::size_t>(i)].atomIndex == i);
  }
  // Atom 0 should list atoms 1 and 2 as neighbours (by index)
  REQUIRE(graph.pts[0].neighListIndex.size() == 2);
}

// -- clearGraph tests --

TEST_CASE("clearGraph empties a graph", "[franzblau]") {
  // Build a graph with some data
  std::vector<std::vector<int>> nList = {{0, 1, 2}, {1, 0, 2}, {2, 0, 1}};
  auto graph = primitive::populateGraphFromIndices(nList);
  REQUIRE(graph.pts.size() == 3);

  // Add a synthetic ring
  graph.rings.push_back({0, 1, 2});
  REQUIRE(graph.rings.size() == 1);

  auto cleared = primitive::clearGraph(graph);

  // After clearing, the original graph should be empty (swapped)
  REQUIRE(graph.pts.empty());
  REQUIRE(graph.rings.empty());
}

// -- ring::compareRings tests --

TEST_CASE("compareRings returns true for same elements in different order",
          "[ring]") {
  std::vector<int> ring1 = {3, 1, 2, 0};
  std::vector<int> ring2 = {0, 1, 2, 3};
  REQUIRE(ring::compareRings(ring1, ring2) == true);
}

TEST_CASE("compareRings returns false for different elements", "[ring]") {
  std::vector<int> ring1 = {0, 1, 2};
  std::vector<int> ring2 = {0, 1, 3};
  REQUIRE(ring::compareRings(ring1, ring2) == false);
}

TEST_CASE("compareRings returns false for different sizes", "[ring]") {
  std::vector<int> ring1 = {0, 1, 2};
  std::vector<int> ring2 = {0, 1, 2, 3};
  REQUIRE(ring::compareRings(ring1, ring2) == false);
}

TEST_CASE("compareRings with empty rings returns true", "[ring]") {
  std::vector<int> ring1 = {};
  std::vector<int> ring2 = {};
  REQUIRE(ring::compareRings(ring1, ring2) == true);
}

// ---------------------------------------------------------------------------
// The shortest-path-pair construction in primitive::ringNetwork and the
// generate-then-filter route through primitive::countAllRingsFromIndex plus
// primitive::removeNonSPrings must describe the same ring set. The comparison
// is over canonical forms, since the two routes emit members in different
// rotations, and runs at several disorder levels because an ordered lattice
// cannot distinguish them: every candidate the backtracking generates there is
// already primitive.
// ---------------------------------------------------------------------------

namespace {

std::vector<int> canonicalRing(std::vector<int> r) {
  const int n = static_cast<int>(r.size());
  std::vector<int> best;
  for (int dir = 0; dir < 2; dir++) {
    for (int s = 0; s < n; s++) {
      std::vector<int> cand(static_cast<std::size_t>(n));
      for (int i = 0; i < n; i++) {
        cand[static_cast<std::size_t>(i)] =
            r[static_cast<std::size_t>((dir ? (s - i + 2 * n) : (s + i)) % n)];
      }
      if (best.empty() || cand < best) {
        best = cand;
      }
    }
  }
  return best;
}

// Jittered cubic network at the mW number density, matching the scaling
// benchmarks, with the neighbour list already by index
std::vector<std::vector<int>> jitteredNetwork(int nAtoms, double jitterFrac) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double density = 0.0332;
  const double boxLength = std::cbrt(nAtoms / density);
  const int perSide = static_cast<int>(std::ceil(std::cbrt(nAtoms)));
  const double spacing = boxLength / perSide;
  cloud.nop = nAtoms;
  cloud.currentFrame = 1;
  cloud.box = {boxLength, boxLength, boxLength};
  cloud.boxLow = {0.0, 0.0, 0.0};
  unsigned long long state = 88172645463325252ULL;
  auto jitter = [&state, spacing, jitterFrac]() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    const double unit = static_cast<double>(state >> 11) / 9007199254740992.0;
    return (unit - 0.5) * jitterFrac * spacing;
  };
  for (int i = 0; i < nAtoms; i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    p.x = ((i % perSide) + 0.5) * spacing + jitter();
    p.y = (((i / perSide) % perSide) + 0.5) * spacing + jitter();
    p.z = (((i / (perSide * perSide)) % perSide) + 0.5) * spacing + jitter();
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = i;
  }
  auto nList = nneigh::neighListO(3.5, cloud, 1);
  return nneigh::neighbourListByIndex(cloud, nList);
}

} // namespace

TEST_CASE("ringNetwork matches the generate-then-filter route", "[franzblau]") {
  for (const double jitterFrac : {0.05, 0.30, 0.60}) {
    for (const int maxDepth : {3, 4, 5, 6, 7, 8}) {
      INFO("jitter " << jitterFrac << ", maxDepth " << maxDepth);
      auto idx = jitteredNetwork(400, jitterFrac);

      auto fast = primitive::ringNetwork(idx, maxDepth);

      auto graph = primitive::countAllRingsFromIndex(idx, maxDepth);
      primitive::removeNonSPrings(graph);

      std::set<std::vector<int>> fastSet, refSet;
      for (auto &r : fast) {
        fastSet.insert(canonicalRing(r));
      }
      for (auto &r : graph.rings) {
        refSet.insert(canonicalRing(r));
      }

      // No duplicates from the directed construction
      REQUIRE(fastSet.size() == fast.size());
      REQUIRE(fastSet == refSet);
    }
  }
}

TEST_CASE("ringNetwork lists its rings in canonical sorted order",
          "[franzblau]") {
  // Callers and the incremental updater compare ring lists element by
  // element, so the order may depend on the graph only: each ring from its
  // lowest member toward the lesser neighbour, the list sorted, and the same
  // list for any order of the neighbour rows
  const auto idx = jitteredNetwork(2000, 0.30);
  auto shuffled = idx;
  unsigned long long state = 977ULL;
  for (auto &row : shuffled) {
    for (std::size_t k = row.size() - 1; k > 1; k--) {
      state ^= state << 13;
      state ^= state >> 7;
      state ^= state << 17;
      const auto j = 1 + static_cast<std::size_t>(state % k);
      std::swap(row[k], row[j]);
    }
  }
  for (int maxDepth = 3; maxDepth <= 8; maxDepth++) {
    INFO("maxDepth " << maxDepth);
    const auto rings = primitive::ringNetwork(idx, maxDepth);
    REQUIRE(std::is_sorted(rings.begin(), rings.end()));
    for (const auto &r : rings) {
      REQUIRE(r == canonicalRing(r));
    }
    REQUIRE(primitive::ringNetwork(shuffled, maxDepth) == rings);
  }
}

TEST_CASE("ringNetwork stays linear in the number of vertices", "[franzblau]") {
  // A pass per source over every lower index would take minutes here
  std::vector<std::vector<int>> isolated(1000000);
  for (int i = 0; i < static_cast<int>(isolated.size()); i++) {
    isolated[static_cast<std::size_t>(i)].push_back(i);
  }
  const auto t0 = std::chrono::steady_clock::now();
  const auto rings = primitive::ringNetwork(isolated, 6);
  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  REQUIRE(rings.empty());
  REQUIRE(seconds < 10.0);
}

// ---------------------------------------------------------------------------
// The incremental updater must be indistinguishable from a full recomputation
// on every frame, while re-enumerating only sources within the locality
// radius of a change.
// ---------------------------------------------------------------------------

namespace {

std::vector<std::vector<int>> networkFromCloud(
    molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  auto nList = nneigh::neighListO(3.5, cloud, 1);
  return nneigh::neighbourListByIndex(cloud, nList);
}

molSys::PointCloud<molSys::Point<double>, double>
jitteredCloud(int nAtoms, double jitterFrac, unsigned long long seed) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double density = 0.0332;
  const double boxLength = std::cbrt(nAtoms / density);
  const int perSide = static_cast<int>(std::ceil(std::cbrt(nAtoms)));
  const double spacing = boxLength / perSide;
  cloud.nop = nAtoms;
  cloud.currentFrame = 1;
  cloud.box = {boxLength, boxLength, boxLength};
  cloud.boxLow = {0.0, 0.0, 0.0};
  unsigned long long state = seed;
  auto jitter = [&state, spacing, jitterFrac]() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    const double unit = static_cast<double>(state >> 11) / 9007199254740992.0;
    return (unit - 0.5) * jitterFrac * spacing;
  };
  for (int i = 0; i < nAtoms; i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    p.x = ((i % perSide) + 0.5) * spacing + jitter();
    p.y = (((i / perSide) % perSide) + 0.5) * spacing + jitter();
    p.z = (((i / (perSide * perSide)) % perSide) + 0.5) * spacing + jitter();
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = i;
  }
  return cloud;
}

std::set<std::vector<int>> canonicalSet(const std::vector<std::vector<int>> &rings) {
  std::set<std::vector<int>> out;
  for (auto r : rings) {
    out.insert(canonicalRing(r));
  }
  return out;
}

} // namespace

TEST_CASE("RingUpdater equals full recomputation across perturbed frames",
          "[franzblau]") {
  // The locality radius is 2*maxLvl - 1 hops; the update only wins once the
  // system dwarfs that ball, so the test must too
  const int n = 2000;
  auto cloud = jitteredCloud(n, 0.30, 88172645463325252ULL);
  primitive::RingUpdater updater(6);

  // Frame one: a full pass
  auto idx = networkFromCloud(cloud);
  auto incr = updater.update(idx);
  REQUIRE(updater.lastRecomputedSources() == n);
  REQUIRE(canonicalSet(incr) == canonicalSet(primitive::ringNetwork(idx, 6)));

  // Same frame again: nothing to do
  incr = updater.update(idx);
  REQUIRE(updater.lastRecomputedSources() == 0);
  REQUIRE(canonicalSet(incr) == canonicalSet(primitive::ringNetwork(idx, 6)));

  // Three frames of local perturbation: displace a few atoms enough to move
  // bonds, and require exact agreement with a full recomputation each time,
  // with strictly fewer sources re-enumerated than a full pass
  unsigned long long state = 424242ULL;
  auto pick = [&state](int bound) {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return static_cast<int>((state >> 11) % static_cast<unsigned long long>(bound));
  };
  for (int frame = 0; frame < 3; frame++) {
    for (int k = 0; k < 3; k++) {
      auto &p = cloud.pts[static_cast<std::size_t>(pick(n))];
      p.x += 0.8 * ((frame + k) % 2 ? 1.0 : -1.0);
      p.y += 0.4;
    }
    idx = networkFromCloud(cloud);
    incr = updater.update(idx);
    INFO("frame " << frame << ", recomputed "
                  << updater.lastRecomputedSources() << " of " << n);
    const auto full = primitive::ringNetwork(idx, 6);
    REQUIRE(incr == full);
    REQUIRE(canonicalSet(incr) == canonicalSet(full));
    REQUIRE(updater.lastRecomputedSources() < n);
    REQUIRE(updater.lastBallsRefreshed() < n);
  }
}

TEST_CASE("RingUpdater matches ringNetwork element by element on every frame",
          "[franzblau]") {
  const int n = 1200;
  auto cloud = jitteredCloud(n, 0.20, 88172645463325252ULL);
  primitive::RingUpdater updater(6);
  unsigned long long state = 99ULL;
  auto rnd = [&state]() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
  };
  for (int frame = 0; frame < 6; frame++) {
    const int moves = frame < 2 ? 0 : (frame < 4 ? 4 : n / 3);
    for (int k = 0; k < moves; k++) {
      auto &p =
          cloud.pts[static_cast<std::size_t>(rnd() % static_cast<unsigned long long>(n))];
      p.x += 1.1 * ((frame + k) % 2 ? 1.0 : -1.0);
      p.y += 0.6;
      p.z -= 0.4;
    }
    const auto idx = networkFromCloud(cloud);
    const auto &inc = updater.update(idx);
    const auto full = primitive::ringNetwork(idx, 6);
    INFO("frame " << frame << " bonds/molecule "
                  << updater.lastChangedBondFraction() << " recomputed "
                  << updater.lastRecomputedSources());
    REQUIRE(inc == full);
    if (frame < 2) {
      REQUIRE(updater.lastChangedBonds() == 0);
    }
    if (frame >= 4) {
      REQUIRE(updater.lastChangedBonds() > 0);
    }
  }
}

TEST_CASE("ringNetwork keeps one copy of a repeated edge", "[franzblau]") {
  const std::vector<std::vector<int>> nList = {
      {0, 1, 1, 2},
      {1, 0, 0, 2},
      {2, 0, 1},
  };
  const auto rings = primitive::ringNetwork(nList, 6, {});
  REQUIRE(rings.size() == 1);
  REQUIRE(rings[0].size() == 3);
}

TEST_CASE("populateGraphFromNListID skips a missing id and keeps the other vertices",
          "[franzblau]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = 3;
  cloud.pts.resize(3);
  for (int i = 0; i < 3; i++) {
    cloud.pts[static_cast<std::size_t>(i)].atomID = i + 1;
    cloud.idIndexMap[i + 1] = i;
  }
  const std::vector<std::vector<int>> nList = {{1, 2}, {2, 99}, {3, 1}};
  const auto graph = primitive::populateGraphFromNListID(cloud, nList);
  REQUIRE(graph.pts.size() == 3);
  REQUIRE(graph.pts[1].neighListIndex.empty());
  REQUIRE(graph.pts[0].neighListIndex.size() == 1);
}
