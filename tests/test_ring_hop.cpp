//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <franzblau.hpp>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstdint>
#include <queue>
#include <set>
#include <vector>

namespace {

std::vector<std::vector<int>> withSelf(const std::vector<std::vector<int>> &adj) {
  std::vector<std::vector<int>> nList(adj.size());
  for (int i = 0; i < static_cast<int>(adj.size()); i++) {
    nList[static_cast<std::size_t>(i)].push_back(i);
    for (const int j : adj[static_cast<std::size_t>(i)]) {
      nList[static_cast<std::size_t>(i)].push_back(j);
    }
  }
  return nList;
}

int bfsDist(const std::vector<std::vector<int>> &adj, int src, int dst) {
  if (src == dst) {
    return 0;
  }
  const int n = static_cast<int>(adj.size());
  std::vector<int> dist(static_cast<std::size_t>(n), -1);
  std::queue<int> q;
  dist[static_cast<std::size_t>(src)] = 0;
  q.push(src);
  while (!q.empty()) {
    const int u = q.front();
    q.pop();
    for (const int v : adj[static_cast<std::size_t>(u)]) {
      if (dist[static_cast<std::size_t>(v)] == -1) {
        dist[static_cast<std::size_t>(v)] = dist[static_cast<std::size_t>(u)] + 1;
        if (v == dst) {
          return dist[static_cast<std::size_t>(v)];
        }
        q.push(v);
      }
    }
  }
  return 1 << 20;
}

//! Breadth-first distance that refuses to enter a vertex below `floor`.
//! This is the walk that kept the 62 non-primitive hexagons.
int bfsDistAbove(const std::vector<std::vector<int>> &adj, int src, int dst, int floor) {
  if (src == dst) {
    return 0;
  }
  const int n = static_cast<int>(adj.size());
  std::vector<int> dist(static_cast<std::size_t>(n), -1);
  std::queue<int> q;
  dist[static_cast<std::size_t>(src)] = 0;
  q.push(src);
  while (!q.empty()) {
    const int u = q.front();
    q.pop();
    for (const int v : adj[static_cast<std::size_t>(u)]) {
      if (v < floor || dist[static_cast<std::size_t>(v)] != -1) {
        continue;
      }
      dist[static_cast<std::size_t>(v)] = dist[static_cast<std::size_t>(u)] + 1;
      if (v == dst) {
        return dist[static_cast<std::size_t>(v)];
      }
      q.push(v);
    }
  }
  return 1 << 20;
}

std::vector<int> canonical(std::vector<int> ring) {
  auto it = std::min_element(ring.begin(), ring.end());
  std::rotate(ring.begin(), it, ring.end());
  if (ring.size() > 2 && ring[1] > ring.back()) {
    std::reverse(ring.begin() + 1, ring.end());
  }
  return ring;
}

bool primitiveCycle(const std::vector<std::vector<int>> &adj, const std::vector<int> &ring) {
  const int k = static_cast<int>(ring.size());
  for (int i = 0; i < k; i++) {
    for (int j = i + 1; j < k; j++) {
      const int sep = std::min(j - i, k - (j - i));
      if (bfsDist(adj, ring[static_cast<std::size_t>(i)], ring[static_cast<std::size_t>(j)]) !=
          sep) {
        return false;
      }
    }
  }
  return true;
}

// Every simple cycle up to maxDepth, kept when an independent breadth-first
// distance agrees with the ring distance. This is the specification the
// hop probe and the ball table both have to meet.
std::set<std::vector<int>> brutePrimitive(const std::vector<std::vector<int>> &adj,
                                          int maxDepth) {
  const int n = static_cast<int>(adj.size());
  std::set<std::vector<int>> found;
  for (int start = 0; start < n; start++) {
    std::vector<int> path{start};
    std::vector<char> used(static_cast<std::size_t>(n), 0);
    used[static_cast<std::size_t>(start)] = 1;
    auto dfs = [&](auto &&self) -> void {
      if (static_cast<int>(path.size()) > maxDepth) {
        return;
      }
      const int u = path.back();
      for (const int v : adj[static_cast<std::size_t>(u)]) {
        if (v == start && static_cast<int>(path.size()) >= 3 && path[1] < path.back()) {
          if (primitiveCycle(adj, path)) {
            found.insert(canonical(path));
          }
          continue;
        }
        if (v <= start || used[static_cast<std::size_t>(v)]) {
          continue;
        }
        if (static_cast<int>(path.size()) == maxDepth) {
          continue;
        }
        used[static_cast<std::size_t>(v)] = 1;
        path.push_back(v);
        self(self);
        path.pop_back();
        used[static_cast<std::size_t>(v)] = 0;
      }
    };
    dfs(dfs);
  }
  return found;
}

// Chordless cycles, the face set SQQ grows cages from. A chord is an edge
// between non-consecutive members; a longer shortcut is not a chord.
std::set<std::vector<int>> bruteChordless(const std::vector<std::vector<int>> &adj,
                                          int maxDepth) {
  const int n = static_cast<int>(adj.size());
  std::set<std::vector<int>> found;
  for (int start = 0; start < n; start++) {
    std::vector<int> path{start};
    auto dfs = [&](auto &&self) -> void {
      if (static_cast<int>(path.size()) >= maxDepth) {
        return;
      }
      const int u = path.back();
      for (const int v : adj[static_cast<std::size_t>(u)]) {
        if (v <= start) {
          continue;
        }
        bool chord = false;
        bool closes = false;
        for (int i = 0; i + 1 < static_cast<int>(path.size()); i++) {
          const int earlier = path[static_cast<std::size_t>(i)];
          const auto &row = adj[static_cast<std::size_t>(v)];
          if (std::find(row.begin(), row.end(), earlier) == row.end()) {
            continue;
          }
          if (earlier == start) {
            closes = true;
          } else {
            chord = true;
          }
        }
        if (chord) {
          continue;
        }
        path.push_back(v);
        if (closes && static_cast<int>(path.size()) >= 3 && path[1] < path.back()) {
          found.insert(canonical(path));
        } else if (!closes) {
          self(self);
        }
        path.pop_back();
      }
    };
    dfs(dfs);
  }
  return found;
}

std::set<std::vector<int>> asSet(const std::vector<std::vector<int>> &rings) {
  return {rings.begin(), rings.end()};
}

} // namespace

TEST_CASE("length-2 shortcut is not a chord and is not a primitive ring", "[ring]") {
  // Hexagon 0-1-2-3-4-5 plus vertex 6 joined to the antipodes 0 and 3.
  // The hexagon has no chord. Graph distance 0-3 is 2, ring distance is 3.
  std::vector<std::vector<int>> adj(7);
  const int hex[6] = {0, 1, 2, 3, 4, 5};
  for (int i = 0; i < 6; i++) {
    adj[hex[i]].push_back(hex[(i + 1) % 6]);
    adj[hex[(i + 1) % 6]].push_back(hex[i]);
  }
  adj[0].push_back(6);
  adj[6].push_back(0);
  adj[3].push_back(6);
  adj[6].push_back(3);
  const auto nList = withSelf(adj);
  const std::vector<int> hexagon = {0, 1, 2, 3, 4, 5};

  const auto chordless = bruteChordless(adj, 6);
  REQUIRE(chordless.count(hexagon) == 1);

  const auto hop = primitive::ringNetwork(nList, 6);
  REQUIRE(primitive::lastRingSearchUsedBalls() == 0);
  const auto balls = primitive::ringNetworkBalls(nList, 6);
  REQUIRE(primitive::lastRingSearchUsedBalls() == 1);
  REQUIRE(hop == balls);
  REQUIRE(asSet(hop).count(hexagon) == 0);
  REQUIRE(asSet(hop) == brutePrimitive(adj, 6));
  // The two pentagons through the spur are primitive. A chord-only search
  // would have kept the hexagon as well, and that set is strictly larger.
  REQUIRE(chordless.size() == asSet(hop).size() + 1);

  std::vector<int> offsets(adj.size() + 1, 0);
  std::vector<int> cols;
  for (std::size_t i = 0; i < adj.size(); i++) {
    offsets[i + 1] = offsets[i] + static_cast<int>(adj[i].size());
    cols.insert(cols.end(), adj[i].begin(), adj[i].end());
  }
  const auto fromCsr = primitive::ringNetwork(offsets, cols, 6);
  REQUIRE(fromCsr == hop);
}

TEST_CASE("full-graph distance check, the 62-hexagon shape", "[ring]") {
  // Hexagon 1-2-3-4-5-6, and vertex 0 joined to the antipodes 1 and 4.
  // On the first liquid frame, 62 hexagons had this shape: primitive in the
  // subgraph of vertices at or above the source, and not primitive once the
  // lower vertex is visible. The search that starts at 1 must see vertex 0.
  std::vector<std::vector<int>> adj(7);
  const int hex[6] = {1, 2, 3, 4, 5, 6};
  for (int i = 0; i < 6; i++) {
    adj[hex[i]].push_back(hex[(i + 1) % 6]);
    adj[hex[(i + 1) % 6]].push_back(hex[i]);
  }
  adj[0].push_back(1);
  adj[1].push_back(0);
  adj[0].push_back(4);
  adj[4].push_back(0);
  const auto nList = withSelf(adj);
  const auto hop = primitive::ringNetwork(nList, 6);
  const std::vector<int> hexagon = {1, 2, 3, 4, 5, 6};
  // The walk that stays at indices >= 1 never sees vertex 0, so the
  // antipodes look three apart. The full graph has a path of length 2.
  REQUIRE(bfsDistAbove(adj, 1, 4, 1) == 3);
  REQUIRE(bfsDist(adj, 1, 4) == 2);
  REQUIRE(bruteChordless(adj, 6).count(hexagon) == 1);
  REQUIRE(asSet(hop).count(hexagon) == 0);
  REQUIRE(brutePrimitive(adj, 6).count(hexagon) == 0);
  REQUIRE(asSet(hop) == brutePrimitive(adj, 6));
  REQUIRE(asSet(hop) == asSet(primitive::ringNetworkBalls(nList, 6)));
}

TEST_CASE("a chorded square is rejected by both searches", "[ring]") {
  std::vector<std::vector<int>> adj(4);
  for (int i = 0; i < 4; i++) {
    adj[i].push_back((i + 1) % 4);
    adj[(i + 1) % 4].push_back(i);
  }
  adj[0].push_back(2);
  adj[2].push_back(0);
  const auto nList = withSelf(adj);
  const auto hop = primitive::ringNetwork(nList, 4);
  const auto balls = primitive::ringNetworkBalls(nList, 4);
  REQUIRE(hop == balls);
  // The square is not primitive. Each triangle is.
  const auto got = asSet(hop);
  REQUIRE(got.count(std::vector<int>{0, 1, 2, 3}) == 0);
  REQUIRE(got.count(std::vector<int>{0, 1, 2}) == 1);
  REQUIRE(got.count(std::vector<int>{0, 2, 3}) == 1);
  REQUIRE(got == brutePrimitive(adj, 4));
  REQUIRE(got == bruteChordless(adj, 4));
}

TEST_CASE("hop probe matches the ball table and an independent cycle search", "[ring]") {
  std::uint64_t state = 0x9e3779b97f4a7c15ULL;
  auto next = [&state]() {
    state = state * 6364136223846793005ULL + 1;
    return state;
  };
  for (int n = 4; n <= 9; n++) {
    std::vector<std::vector<int>> adj(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        if ((next() % 3) == 0) {
          adj[static_cast<std::size_t>(i)].push_back(j);
          adj[static_cast<std::size_t>(j)].push_back(i);
        }
      }
    }
    const auto nList = withSelf(adj);
    for (const int depth : {5, 6, 7, 8, 9, 10}) {
      if (depth >= 9 && n > 6) {
        continue;
      }
      const auto hop = primitive::ringNetwork(nList, depth);
      const auto balls = primitive::ringNetworkBalls(nList, depth);
      REQUIRE(hop == balls);
      REQUIRE(asSet(hop) == brutePrimitive(adj, depth));
      (void)primitive::ringNetwork(nList, depth);
      if (depth <= 7) {
        REQUIRE(primitive::lastRingSearchUsedBalls() == 0);
      } else if (depth <= 9) {
        // Radius 3 stays on the sorted pair table.
        REQUIRE(primitive::lastRingSearchUsedBalls() == 1);
      } else if (depth == 10) {
        // Radius-2 pair table, meet in the middle for distance 4.
        REQUIRE(primitive::lastRingSearchUsedBalls() == 3);
      }
      (void)primitive::ringNetworkBalls(nList, depth);
      REQUIRE(primitive::lastRingSearchUsedBalls() == 1);
    }
  }
}

TEST_CASE("a changed bond only rewrites rings through its endpoints", "[ring]") {
  std::vector<std::vector<int>> adj(7);
  for (int i = 0; i < 6; i++) {
    adj[i].push_back((i + 1) % 6);
    adj[(i + 1) % 6].push_back(i);
  }
  adj[0].push_back(6);
  adj[6].push_back(0);
  adj[3].push_back(6);
  adj[6].push_back(3);
  primitive::RingUpdater updater(6);
  const auto first = updater.update(withSelf(adj));
  REQUIRE(first == primitive::ringNetwork(withSelf(adj), 6));
  // Drop the spur. The hexagon becomes primitive; vertices 1, 2, 4 and 5
  // are not endpoints of the deleted bonds.
  adj[0].erase(std::find(adj[0].begin(), adj[0].end(), 6));
  adj[3].erase(std::find(adj[3].begin(), adj[3].end(), 6));
  adj[6].clear();
  const auto second = updater.update(withSelf(adj));
  const auto full = primitive::ringNetwork(withSelf(adj), 6);
  REQUIRE(second == full);
  REQUIRE(asSet(full).count(std::vector<int>{0, 1, 2, 3, 4, 5}) == 1);
}
