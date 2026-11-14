// Internal
#include <franzblau.hpp>

// Standard
#include <algorithm>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

namespace {

std::vector<std::vector<int>>
undirected(std::vector<std::vector<int>> adj) {
  const int n = static_cast<int>(adj.size());
  std::vector<std::vector<int>> out(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    for (int j : adj[static_cast<std::size_t>(i)]) {
      if (j < 0 || j >= n || j == i) {
        continue;
      }
      out[static_cast<std::size_t>(i)].push_back(j);
      out[static_cast<std::size_t>(j)].push_back(i);
    }
  }
  for (auto &row : out) {
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
  }
  return out;
}

std::vector<std::vector<int>>
asNeighbourRows(const std::vector<std::vector<int>> &adj) {
  std::vector<std::vector<int>> rows(adj.size());
  for (int i = 0; i < static_cast<int>(adj.size()); i++) {
    rows[static_cast<std::size_t>(i)].push_back(i);
    for (int j : adj[static_cast<std::size_t>(i)]) {
      rows[static_cast<std::size_t>(i)].push_back(j);
    }
  }
  return rows;
}

int graphDistance(const std::vector<std::vector<int>> &adj, int src, int dst,
                  int cap) {
  if (src == dst) {
    return 0;
  }
  std::vector<int> dist(adj.size(), -1);
  std::queue<int> wait;
  dist[static_cast<std::size_t>(src)] = 0;
  wait.push(src);
  while (!wait.empty()) {
    const int u = wait.front();
    wait.pop();
    if (dist[static_cast<std::size_t>(u)] == cap) {
      continue;
    }
    for (const int w : adj[static_cast<std::size_t>(u)]) {
      if (dist[static_cast<std::size_t>(w)] != -1) {
        continue;
      }
      dist[static_cast<std::size_t>(w)] = dist[static_cast<std::size_t>(u)] + 1;
      if (w == dst) {
        return dist[static_cast<std::size_t>(w)];
      }
      wait.push(w);
    }
  }
  return cap + 1;
}

bool isPrimitive(const std::vector<std::vector<int>> &adj,
                 const std::vector<int> &ring) {
  const int k = static_cast<int>(ring.size());
  for (int i = 0; i < k; i++) {
    for (int j = i + 1; j < k; j++) {
      const int sep = j - i;
      const int ringDist = std::min(sep, k - sep);
      if (ringDist <= 1) {
        continue;
      }
      if (graphDistance(adj, ring[static_cast<std::size_t>(i)],
                        ring[static_cast<std::size_t>(j)], ringDist) <
          ringDist) {
        return false;
      }
    }
  }
  return true;
}

// Every simple cycle of length 3..maxDepth, once, in the same canonical
// writing ringNetwork uses: lowest index first, then the lesser neighbour.
std::vector<std::vector<int>>
bruteCycles(const std::vector<std::vector<int>> &adj, int maxDepth) {
  const int n = static_cast<int>(adj.size());
  std::vector<std::vector<int>> found;
  std::vector<char> seen(static_cast<std::size_t>(n), 0);
  std::vector<int> path;
  for (int start = 0; start < n; start++) {
    path.assign(1, start);
    seen[static_cast<std::size_t>(start)] = 1;
    const auto walk = [&](auto &&self, int v) -> void {
      for (const int w : adj[static_cast<std::size_t>(v)]) {
        if (w == start && static_cast<int>(path.size()) >= 3 &&
            static_cast<int>(path.size()) <= maxDepth && path[1] < path.back()) {
          found.push_back(path);
          continue;
        }
        if (w > start && !seen[static_cast<std::size_t>(w)] &&
            static_cast<int>(path.size()) < maxDepth) {
          seen[static_cast<std::size_t>(w)] = 1;
          path.push_back(w);
          self(self, w);
          path.pop_back();
          seen[static_cast<std::size_t>(w)] = 0;
        }
      }
    };
    walk(walk, start);
    seen[static_cast<std::size_t>(start)] = 0;
  }
  std::vector<std::vector<int>> primitive;
  for (const auto &ring : found) {
    if (isPrimitive(adj, ring)) {
      primitive.push_back(ring);
    }
  }
  std::sort(primitive.begin(), primitive.end());
  return primitive;
}

std::map<int, int> census(const std::vector<std::vector<int>> &rings) {
  std::map<int, int> counts;
  for (const auto &ring : rings) {
    counts[static_cast<int>(ring.size())] += 1;
  }
  return counts;
}

void expectSameCensus(const std::vector<std::vector<int>> &adj, int maxDepth) {
  const auto rows = asNeighbourRows(adj);
  auto got = primitive::ringNetwork(rows, maxDepth);
  std::sort(got.begin(), got.end());
  const auto brute = bruteCycles(adj, maxDepth);
  INFO("maxDepth " << maxDepth << " ringNetwork " << got.size() << " brute "
                   << brute.size());
  REQUIRE(census(got) == census(brute));
  REQUIRE(got == brute);
}

} // namespace

TEST_CASE("primitive census matches brute-force cycles", "[ring][census]") {
  SECTION("triangle beside a square") {
    // 0-1-2 and 3-4-5-6, the example in franzblau-test.
    auto adj = undirected({{1, 2}, {0, 2}, {0, 1}, {4, 6}, {3, 5}, {4, 6}, {3, 5}});
    expectSameCensus(adj, 6);
    REQUIRE(census(bruteCycles(adj, 6)) == std::map<int, int>{{3, 1}, {4, 1}});
  }

  SECTION("chorded diamond keeps the two triangles") {
    // Two triangles share the edge 0-1. The 4-cycle has that chord.
    auto adj = undirected({{1, 2, 3}, {0, 2, 3}, {0, 1}, {0, 1}});
    expectSameCensus(adj, 6);
    REQUIRE(census(bruteCycles(adj, 6)) == std::map<int, int>{{3, 2}});
  }

  SECTION("hexagon and a pendant edge") {
    auto adj = undirected(
        {{1, 5}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0, 6}, {5}});
    expectSameCensus(adj, 6);
    REQUIRE(census(bruteCycles(adj, 7)) == std::map<int, int>{{6, 1}});
    expectSameCensus(adj, 5);
    REQUIRE(census(bruteCycles(adj, 5)).empty());
  }

  SECTION("cube faces and equators") {
    // Vertices 0..7 of the cube. Faces are primitive 4-rings. Each pair of
    // opposite vertices leaves a primitive 6-ring.
    std::vector<std::vector<int>> edges(8);
    const int cube[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                             {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    for (const auto &e : cube) {
      edges[static_cast<std::size_t>(e[0])].push_back(e[1]);
    }
    auto adj = undirected(edges);
    expectSameCensus(adj, 4);
    REQUIRE(census(bruteCycles(adj, 4)) == std::map<int, int>{{4, 6}});
    expectSameCensus(adj, 6);
    REQUIRE(census(bruteCycles(adj, 6)).at(4) == 6);
    REQUIRE(census(bruteCycles(adj, 6)).at(6) == 4);
  }

  SECTION("a small cubic patch with one chord") {
    // A 6-cycle 0-1-2-3-4-5 plus the chord 0-3, and a triangle on 4-5-6.
    auto adj = undirected({{1, 5, 3},
                           {0, 2},
                           {1, 3},
                           {2, 4, 0},
                           {3, 5, 6},
                           {4, 0, 6},
                           {4, 5}});
    expectSameCensus(adj, 6);
    expectSameCensus(adj, 7);
  }
}
