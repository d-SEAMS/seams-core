// Copyright (c) d-SEAMS developers
#pragma once

// Brute-force reference for the primitive ring search.
//
// Every simple cycle of 3 to maxDepth vertices is listed by a plain
// depth-first walk from its lowest vertex. A cycle is kept when no pair of
// its members is joined through the graph by fewer edges than the shorter
// way around the cycle (Franzblau's shortest-path rings), with distances
// from a breadth-first search of the whole graph. Nothing here
// shares code with src/franzblau.cpp or src/chordless.cpp.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <set>
#include <utility>
#include <vector>

namespace oracle {

using Adjacency = std::vector<std::vector<int>>;
using Rings = std::vector<std::vector<int>>;

// Simple undirected graph from arbitrary rows: self entries, repeats and
// indices out of range are dropped, and every edge is made symmetric.
inline Adjacency simpleGraph(const std::vector<std::vector<int>> &rows, int n) {
  std::vector<std::set<int>> sets(static_cast<std::size_t>(n));
  for (int i = 0; i < n && i < static_cast<int>(rows.size()); ++i) {
    for (int j : rows[static_cast<std::size_t>(i)]) {
      if (j < 0 || j >= n || j == i) {
        continue;
      }
      sets[static_cast<std::size_t>(i)].insert(j);
      sets[static_cast<std::size_t>(j)].insert(i);
    }
  }
  Adjacency adj(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    adj[static_cast<std::size_t>(i)].assign(sets[static_cast<std::size_t>(i)].begin(),
                                            sets[static_cast<std::size_t>(i)].end());
  }
  return adj;
}

// Starts at the lowest member and runs toward the lesser of its two cycle
// neighbours.
inline std::vector<int> canonical(std::vector<int> cyc) {
  const auto lo = std::min_element(cyc.begin(), cyc.end());
  std::rotate(cyc.begin(), lo, cyc.end());
  if (cyc.size() > 2 && cyc.back() < cyc[1]) {
    std::reverse(cyc.begin() + 1, cyc.end());
  }
  return cyc;
}

inline bool hasEdge(const Adjacency &adj, int a, int b) {
  const auto &row = adj[static_cast<std::size_t>(a)];
  return std::binary_search(row.begin(), row.end(), b);
}

// Every simple cycle of 3..maxDepth vertices, once each, canonical.
inline Rings simpleCycles(const Adjacency &adj, int maxDepth) {
  Rings out;
  const int n = static_cast<int>(adj.size());
  std::vector<char> on(static_cast<std::size_t>(n), 0);
  std::vector<int> path;
  for (int s = 0; s < n; ++s) {
    path.assign(1, s);
    on[static_cast<std::size_t>(s)] = 1;
    // Explicit stack of neighbour cursors.
    std::vector<std::size_t> cursor(1, 0);
    while (!cursor.empty()) {
      const int u = path.back();
      const auto &row = adj[static_cast<std::size_t>(u)];
      if (cursor.back() >= row.size()) {
        on[static_cast<std::size_t>(u)] = 0;
        path.pop_back();
        cursor.pop_back();
        if (path.empty()) {
          on[static_cast<std::size_t>(s)] = 0;
        }
        continue;
      }
      const int v = row[cursor.back()++];
      if (v == s && path.size() >= 3 && path[1] < path.back()) {
        out.push_back(path);
        continue;
      }
      if (v <= s || on[static_cast<std::size_t>(v)] ||
          static_cast<int>(path.size()) >= maxDepth) {
        continue;
      }
      on[static_cast<std::size_t>(v)] = 1;
      path.push_back(v);
      cursor.push_back(0);
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

// Graph distances out to a radius, per vertex, as sorted (vertex, hops).
// A shortcut between two members of a ring of k <= maxDepth vertices is
// shorter than k / 2 hops, so a radius of maxDepth / 2 decides every pair.
struct Balls {
  int radius = 0;
  std::vector<std::vector<std::pair<int, int>>> rows;

  int distance(int a, int b) const {
    const auto &row = rows[static_cast<std::size_t>(a)];
    const auto it = std::lower_bound(row.begin(), row.end(), std::make_pair(b, -1));
    return (it != row.end() && it->first == b) ? it->second : radius + 1;
  }
};

inline Balls balls(const Adjacency &adj, int radius) {
  Balls out;
  out.radius = radius;
  out.rows.resize(adj.size());
  std::vector<int> d(adj.size(), -1);
  std::vector<int> seen;
  for (std::size_t s = 0; s < adj.size(); ++s) {
    seen.assign(1, static_cast<int>(s));
    d[s] = 0;
    for (std::size_t head = 0; head < seen.size(); ++head) {
      const int u = seen[head];
      if (d[static_cast<std::size_t>(u)] == radius) {
        continue;
      }
      for (int v : adj[static_cast<std::size_t>(u)]) {
        if (d[static_cast<std::size_t>(v)] < 0) {
          d[static_cast<std::size_t>(v)] = d[static_cast<std::size_t>(u)] + 1;
          seen.push_back(v);
        }
      }
    }
    auto &row = out.rows[s];
    for (int v : seen) {
      row.emplace_back(v, d[static_cast<std::size_t>(v)]);
      d[static_cast<std::size_t>(v)] = -1;
    }
    std::sort(row.begin(), row.end());
  }
  return out;
}

inline bool isShortestPathRing(const std::vector<int> &ring, const Balls &dist) {
  const int k = static_cast<int>(ring.size());
  for (int a = 0; a < k; ++a) {
    for (int b = a + 1; b < k; ++b) {
      const int around = std::min(b - a, k - (b - a));
      if (dist.distance(ring[static_cast<std::size_t>(a)], ring[static_cast<std::size_t>(b)]) <
          around) {
        return false;
      }
    }
  }
  return true;
}

inline Rings primitiveRings(const Adjacency &adj, int maxDepth) {
  const Balls dist = balls(adj, std::max(maxDepth, 0) / 2);
  Rings out;
  for (auto &cyc : simpleCycles(adj, maxDepth)) {
    if (isShortestPathRing(cyc, dist)) {
      out.push_back(std::move(cyc));
    }
  }
  return out;
}

inline bool isChordless(const std::vector<int> &ring, const Adjacency &adj) {
  const int k = static_cast<int>(ring.size());
  for (int a = 0; a < k; ++a) {
    for (int b = a + 2; b < k; ++b) {
      if (a == 0 && b == k - 1) {
        continue;
      }
      if (hasEdge(adj, ring[static_cast<std::size_t>(a)], ring[static_cast<std::size_t>(b)])) {
        return false;
      }
    }
  }
  return true;
}

inline Rings chordlessCycles(const Adjacency &adj, int size) {
  Rings out;
  for (auto &cyc : simpleCycles(adj, size)) {
    if (static_cast<int>(cyc.size()) == size && isChordless(cyc, adj)) {
      out.push_back(std::move(cyc));
    }
  }
  return out;
}

inline Rings canonicalSorted(Rings rings) {
  for (auto &r : rings) {
    r = canonical(r);
  }
  std::sort(rings.begin(), rings.end());
  return rings;
}

// Lattice rows a, b, c of a LAMMPS box9 (ox, oy, oz, lx, ly, lz, xy, xz, yz).
inline std::array<std::array<double, 3>, 3> lattice(const double *box9) {
  return {{{box9[3], 0.0, 0.0}, {box9[6], box9[4], 0.0}, {box9[7], box9[8], box9[5]}}};
}

// Cutoff graph under periodic images, by trying every image within two
// cells: i and j bond when any image of j lies closer than cutoff.
inline Adjacency cutoffGraph(const std::vector<double> &xyz, const double *box9,
                             double cutoff, double *closestGap = nullptr) {
  const int n = static_cast<int>(xyz.size() / 3);
  const auto h = lattice(box9);
  Adjacency adj(static_cast<std::size_t>(n));
  double gap = 1e300;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double best = 1e300;
      for (int a = -2; a <= 2; ++a) {
        for (int b = -2; b <= 2; ++b) {
          for (int c = -2; c <= 2; ++c) {
            double d2 = 0;
            for (int k = 0; k < 3; ++k) {
              const std::size_t kk = static_cast<std::size_t>(k);
              const double dx = xyz[static_cast<std::size_t>(3 * j + k)] -
                                xyz[static_cast<std::size_t>(3 * i + k)] + a * h[0][kk] +
                                b * h[1][kk] + c * h[2][kk];
              d2 += dx * dx;
            }
            best = std::min(best, d2);
          }
        }
      }
      const double d = std::sqrt(best);
      gap = std::min(gap, std::abs(d - cutoff));
      if (d < cutoff) {
        adj[static_cast<std::size_t>(i)].push_back(j);
        adj[static_cast<std::size_t>(j)].push_back(i);
      }
    }
  }
  for (auto &row : adj) {
    std::sort(row.begin(), row.end());
  }
  if (closestGap != nullptr) {
    *closestGap = gap;
  }
  return adj;
}

} // namespace oracle
