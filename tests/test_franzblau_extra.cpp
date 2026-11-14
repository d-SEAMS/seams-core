#include <catch2/catch_test_macros.hpp>

#include <franzblau.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <map>
#include <queue>
#include <set>
#include <thread>
#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace primitive {
std::vector<std::vector<int>>
ringNetworkAtCap(const std::vector<std::vector<int>> &nList, int maxDepth, int cap);
std::vector<int> ringCensusCounts(const std::vector<std::vector<int>> &nList, int maxDepth);
}
#include <utility>
#include <vector>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>

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

TEST_CASE("ringNetwork on a dense box matches the filtered search",
          "[franzblau]") {
  // Three complete bipartite layers. A vertex in the last layer has
  // width*width shortest paths of length 3, past the pairing cap, so the
  // chordless walk has to list the same primitive rings.
  const int width = 8;
  const int n = 1 + 3 * width;
  std::vector<std::vector<int>> idx(static_cast<std::size_t>(n));
  auto link = [&](int a, int b) {
    idx[static_cast<std::size_t>(a)].push_back(b);
    idx[static_cast<std::size_t>(b)].push_back(a);
  };
  for (int i = 0; i < n; i++) {
    idx[static_cast<std::size_t>(i)].push_back(i);
  }
  for (int layer = 0; layer < 3; layer++) {
    const int left = layer == 0 ? 0 : 1 + (layer - 1) * width;
    const int leftN = layer == 0 ? 1 : width;
    const int right = 1 + layer * width;
    for (int a = 0; a < leftN; a++) {
      for (int b = 0; b < width; b++) {
        link(left + a, right + b);
      }
    }
  }
  for (const int maxDepth : {5, 6}) {
    INFO("maxDepth " << maxDepth);
    const auto fast = primitive::ringNetwork(idx, maxDepth);
    auto graph = primitive::countAllRingsFromIndex(idx, maxDepth);
    primitive::removeNonSPrings(graph);
    std::set<std::vector<int>> fastSet;
    std::set<std::vector<int>> refSet;
    for (const auto &r : fast) {
      fastSet.insert(canonicalRing(r));
      REQUIRE(r == canonicalRing(r));
    }
    for (auto &r : graph.rings) {
      refSet.insert(canonicalRing(r));
    }
    REQUIRE(fastSet.size() == fast.size());
    REQUIRE(std::is_sorted(fast.begin(), fast.end()));
    REQUIRE(fastSet == refSet);
  }
}

TEST_CASE("a failing size-7 arm matches the pairing walk", "[franzblau]") {
  // 0-1-3-5-6-4-2 is a primitive 7-ring. A second arm 0-7-8-5 gives vertex
  // 5 two shortest paths, so cap 1 takes the level walk. Arm 0-9-10-5
  // already meets 6, through 10-6 and through 9-12-6, and is dropped
  // before it is paired.
  std::vector<std::vector<int>> idx(13);
  auto link = [&](int a, int b) {
    idx[static_cast<std::size_t>(a)].push_back(b);
    idx[static_cast<std::size_t>(b)].push_back(a);
  };
  for (int i = 0; i < 13; i++) {
    idx[static_cast<std::size_t>(i)].push_back(i);
  }
  for (const auto edge : {std::pair{0, 1}, {1, 3}, {3, 5}, {5, 6}, {6, 4}, {4, 2}, {2, 0},
                          {0, 7}, {7, 8}, {8, 5}, {0, 9}, {9, 10}, {10, 5}, {10, 6}, {9, 12},
                          {12, 6}}) {
    link(edge.first, edge.second);
  }
  const auto paired = primitive::ringNetworkAtCap(idx, 7, 0);
  const auto levels = primitive::ringNetworkAtCap(idx, 7, 1);
  REQUIRE(paired == levels);
  std::set<std::vector<int>> got;
  for (const auto &ring : levels) {
    got.insert(canonicalRing(ring));
  }
  REQUIRE(got.count(canonicalRing({0, 1, 3, 5, 6, 4, 2})) == 1);
  REQUIRE(got.count(canonicalRing({0, 9, 10, 5, 6, 4, 2})) == 0);
}

TEST_CASE("two ring searches at once keep their own rings", "[franzblau]") {
  std::vector<std::vector<int>> sparse(120);
  for (int c = 0; c < 20; c++) {
    for (int k = 0; k < 6; k++) {
      const int v = c * 6 + k;
      sparse[static_cast<std::size_t>(v)] = {v, c * 6 + (k + 5) % 6,
                                             c * 6 + (k + 1) % 6};
    }
  }
  const int width = 8;
  const int n = 1 + 3 * width;
  std::vector<std::vector<int>> dense(static_cast<std::size_t>(n));
  auto link = [&](int a, int b) {
    dense[static_cast<std::size_t>(a)].push_back(b);
    dense[static_cast<std::size_t>(b)].push_back(a);
  };
  for (int i = 0; i < n; i++) {
    dense[static_cast<std::size_t>(i)].push_back(i);
  }
  for (int layer = 0; layer < 3; layer++) {
    const int left = layer == 0 ? 0 : 1 + (layer - 1) * width;
    const int leftN = layer == 0 ? 1 : width;
    const int right = 1 + layer * width;
    for (int a = 0; a < leftN; a++) {
      for (int b = 0; b < width; b++) {
        link(left + a, right + b);
      }
    }
  }

  const auto serialSparse = primitive::ringNetwork(sparse, 7);
  const auto serialDense = primitive::ringNetwork(dense, 6);
  const auto censusSparse = primitive::ringCensusCounts(sparse, 7);
  const auto censusDense = primitive::ringCensusCounts(dense, 6);
  auto tally = [](const std::vector<std::vector<int>> &rings, int maxDepth) {
    std::vector<int> hist(static_cast<std::size_t>(maxDepth) + 1, 0);
    for (const auto &ring : rings) {
      if (static_cast<int>(ring.size()) <= maxDepth) {
        hist[ring.size()]++;
      }
    }
    return hist;
  };
  REQUIRE(serialSparse.size() == 20);
  REQUIRE(!serialDense.empty());
  REQUIRE(censusSparse == tally(serialSparse, 7));
  REQUIRE(censusDense == tally(serialDense, 6));

  for (int round = 0; round < 8; round++) {
    std::vector<std::vector<int>> gotSparse, gotDense;
    std::vector<int> gotCensusSparse, gotCensusDense;
    std::thread sparseThread([&] {
      gotSparse = primitive::ringNetwork(sparse, 7);
      gotCensusSparse = primitive::ringCensusCounts(sparse, 7);
    });
    std::thread denseThread([&] {
      gotDense = primitive::ringNetwork(dense, 6);
      gotCensusDense = primitive::ringCensusCounts(dense, 6);
    });
    sparseThread.join();
    denseThread.join();
    REQUIRE(gotSparse == serialSparse);
    REQUIRE(gotDense == serialDense);
    REQUIRE(gotCensusSparse == censusSparse);
    REQUIRE(gotCensusDense == censusDense);
  }
}

TEST_CASE("cropped dump-l rings match the pairing walk and the filtered search",
          "[franzblau]") {
  // The first atoms of the dense frame, in the original box. Small enough
  // that the pairing walk returns, and dense enough that a cap of one takes
  // the level walk.
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO("traj/dump-l.lammpstrj", 1, cloud, 1);
  REQUIRE(cloud.nop > 64);
  // Every eighth atom. The piece stays on the original box and cutoff,
  // and its shortest-path counts stay small enough for the pairing walk.
  std::vector<molSys::Point<double>> kept;
  for (int i = 0; i < cloud.nop; i += 8) {
    kept.push_back(cloud.pts[static_cast<std::size_t>(i)]);
  }
  cloud.pts.swap(kept);
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.idIndexMap.clear();
  for (int i = 0; i < cloud.nop; i++) {
    cloud.idIndexMap[cloud.pts[static_cast<std::size_t>(i)].atomID] = i;
  }
  const auto idx =
      nneigh::neighbourListByIndex(cloud, nneigh::neighListO(3.5, cloud, 1));
  int maxDeg = 0;
  for (const auto &row : idx) {
    maxDeg = std::max(maxDeg, static_cast<int>(row.size()) - 1);
  }
  REQUIRE(maxDeg > 8);

  const auto paired = primitive::ringNetworkAtCap(idx, 7, 0);
  const auto levels = primitive::ringNetworkAtCap(idx, 7, 1);
  REQUIRE(paired == levels);

  auto graph = primitive::countAllRingsFromIndex(idx, 7);
  primitive::removeNonSPrings(graph);
  std::set<std::vector<int>> got, ref;
  for (const auto &r : levels) {
    got.insert(canonicalRing(r));
    REQUIRE(r == canonicalRing(r));
  }
  for (auto &r : graph.rings) {
    ref.insert(canonicalRing(r));
  }
  REQUIRE(got.size() == levels.size());
  REQUIRE(got == ref);
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

TEST_CASE("atom ring counts match the stored rings", "[franzblau]") {
  constexpr int n = 72;
  constexpr int half = 9;
  constexpr int depth = 5;
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    auto &row = rows[static_cast<std::size_t>(i)];
    row.push_back(i);
    for (int d = 1; d <= half; d++) {
      row.push_back((i + d) % n);
      row.push_back((i + n - d) % n);
    }
  }
  const auto rings = primitive::ringNetwork(rows, depth);
  const auto counts = primitive::atomRingSizeCounts(rows, depth);
  const int stride = depth + 1;
  std::vector<std::map<int, int>> tally(static_cast<std::size_t>(n));
  for (const auto &ring : rings) {
    const int sz = static_cast<int>(ring.size());
    for (const int a : ring) {
      tally[static_cast<std::size_t>(a)][sz] += 1;
    }
  }
  for (int a = 0; a < n; a++) {
    for (int sz = 0; sz <= depth; sz++) {
      const int got = counts[static_cast<std::size_t>(a) * static_cast<std::size_t>(stride) +
                             static_cast<std::size_t>(sz)];
      const auto it = tally[static_cast<std::size_t>(a)].find(sz);
      const int expect = it == tally[static_cast<std::size_t>(a)].end() ? 0 : it->second;
      REQUIRE(got == expect);
    }
  }
}

TEST_CASE("a dense listing on several threads matches one thread", "[franzblau]") {
  constexpr int n = 72;
  constexpr int half = 9;
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    auto &row = rows[static_cast<std::size_t>(i)];
    row.push_back(i);
    for (int d = 1; d <= half; d++) {
      row.push_back((i + d) % n);
      row.push_back((i + n - d) % n);
    }
  }
#ifdef SEAMS_HAS_OPENMP
  const int prev = omp_get_max_threads();
  omp_set_num_threads(1);
  const auto serial = primitive::ringNetwork(rows, 5);
  for (int threads : {2, 8}) {
    omp_set_num_threads(threads);
    const auto parallel = primitive::ringNetwork(rows, 5);
    REQUIRE(parallel == serial);
  }
  omp_set_num_threads(prev);
#else
  REQUIRE_FALSE(primitive::ringNetwork(rows, 5).empty());
#endif
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

namespace {

// Shorter arc longer than the breadth-first distance: the ring is not primitive.
int nonPrimitive(const std::vector<std::vector<int>> &adj,
                 const std::vector<std::vector<int>> &rings) {
  auto dist = [&](int src, int limit) {
    std::vector<int> seen(adj.size(), -1);
    std::queue<int> wait;
    seen[static_cast<std::size_t>(src)] = 0;
    wait.push(src);
    while (!wait.empty()) {
      const int u = wait.front();
      wait.pop();
      if (seen[static_cast<std::size_t>(u)] == limit) {
        continue;
      }
      for (int v : adj[static_cast<std::size_t>(u)]) {
        if (seen[static_cast<std::size_t>(v)] >= 0) {
          continue;
        }
        seen[static_cast<std::size_t>(v)] = seen[static_cast<std::size_t>(u)] + 1;
        wait.push(v);
      }
    }
    return seen;
  };
  int bad = 0;
  for (const auto &ring : rings) {
    const int length = static_cast<int>(ring.size());
    bool ok = true;
    for (int i = 0; i < length && ok; i++) {
      const auto reached = dist(ring[static_cast<std::size_t>(i)], length / 2);
      for (int j = i + 1; j < length; j++) {
        const int arc = std::min(j - i, length - (j - i));
        const int d = reached[static_cast<std::size_t>(ring[static_cast<std::size_t>(j)])];
        if (d < 0 || d < arc) {
          ok = false;
          break;
        }
      }
    }
    if (!ok) {
      bad += 1;
    }
  }
  return bad;
}

} // namespace

TEST_CASE("K2,3 has three primitive 4-rings", "[franzblau]") {
  // Parts {0, 3} and {1, 2, 4}. Six edges, no chord, three 4-cycles.
  // OVITO FindRings reports two of them. tests/test_ovito_ring_counts.py
  // is that comparison. This case only needs the primitive side.
  // Each row starts with its own index, the same layout neighListO uses.
  const std::vector<std::vector<int>> adj = {
      {0, 1, 2, 4}, {1, 0, 3}, {2, 0, 3}, {3, 1, 2, 4}, {4, 0, 3},
  };
  const auto rings = primitive::ringNetwork(adj, 4);
  REQUIRE(rings.size() == 3);
  for (const auto &ring : rings) {
    REQUIRE(ring.size() == 4);
  }
  REQUIRE(nonPrimitive(adj, rings) == 0);
}

TEST_CASE("the interface frame rings are primitive on the 3.5 A graph", "[franzblau]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj("traj/interface_f0.lammpstrj", 1, cloud);
  REQUIRE(cloud.nop == 864);
  const auto rows = nneigh::neighListO(3.5, cloud, 1);
  const auto indexed = nneigh::neighbourListByIndex(cloud, rows);
  std::vector<std::vector<int>> adj(indexed.size());
  std::size_t bonds = 0;
  for (std::size_t i = 0; i < indexed.size(); i++) {
    for (int j : indexed[i]) {
      if (j == static_cast<int>(i)) {
        continue;
      }
      adj[i].push_back(j);
      if (static_cast<int>(i) < j) {
        bonds += 1;
      }
    }
  }
  REQUIRE(bonds == 1929u);
  // ringNetwork reads the owning index from the front of each row.
  std::vector<std::vector<int>> nList(adj.size());
  for (std::size_t i = 0; i < adj.size(); i++) {
    nList[i].push_back(static_cast<int>(i));
    nList[i].insert(nList[i].end(), adj[i].begin(), adj[i].end());
  }
  const auto rings = primitive::ringNetwork(nList, 10);
  std::map<int, int> hist;
  for (const auto &ring : rings) {
    hist[static_cast<int>(ring.size())] += 1;
  }
  REQUIRE(hist[3] == 142);
  REQUIRE(hist[4] == 188);
  REQUIRE(hist[5] == 377);
  REQUIRE(hist[6] == 1143);
  REQUIRE(hist[7] == 257);
  REQUIRE(hist[8] == 54);
  REQUIRE(hist[9] == 52);
  REQUIRE(hist[10] == 449);
  REQUIRE(nonPrimitive(adj, rings) == 0);
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
