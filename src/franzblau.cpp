//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <numeric>
#include <vector>

#include <franzblau.hpp>
#include <par_team.hpp>
#include <thread_team.hpp>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

struct RingSplitTimes {
  double balls_ms = 0;
  double enumerate_ms = 0;
  double merge_ms = 0;
  double adj_ms = 0;
  double ball_count_ms = 0;
  double ball_fill_ms = 0;
};

RingSplitTimes gRingSplit;
int gUsedBalls = 0;

double millisSince(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
      .count();
}

/**
 * @brief Graph rows packed back to back: row u is cols[first[u]] up to
 *  cols[first[u + 1]]. Two allocations for the whole graph, where a vector
 *  per row costs one each to build and to free on the calling thread.
 */
struct PackedRows {
  std::vector<std::size_t> first;
  std::vector<int> cols;

  struct Row {
    const int *b;
    const int *e;
    const int *begin() const { return b; }
    const int *end() const { return e; }
  };

  Row operator[](std::size_t u) const {
    return {cols.data() + first[u], cols.data() + first[u + 1]};
  }
  std::size_t size() const { return first.empty() ? 0 : first.size() - 1; }
};

/**
 * @brief Sorted (vertex, distance) pairs within a fixed radius of each vertex.
 * @details One bounded breadth-first sweep per vertex for the whole ring
 *  search, cleared by walking what each sweep touched. Every shortcut question
 *  the primitivity test asks is then a sorted lookup, in place of a fresh
 *  sweep per candidate ring.
 */
struct BoundedBalls {
  // One packed table for every ball. off[v] is the first pair of vertex v
  // and off[v + 1] is the end. A vector per vertex was one allocation each.
  std::vector<int> off;
  std::vector<std::pair<int, int>> cols;

  void clear() {
    off.clear();
    cols.clear();
  }

  template <typename Graph>
  void appendBall(int v, const Graph &adjacency, int radius,
                  std::vector<int> &dist, std::vector<int> &touched,
                  std::vector<int> &frontier, std::vector<int> &next,
                  std::vector<std::pair<int, int>> &out) const {
    const int nVertices = static_cast<int>(adjacency.size());
    for (const int w : touched) {
      dist[w] = -1;
    }
    touched.clear();
    dist[v] = 0;
    touched.push_back(v);
    frontier.assign(1, v);
    for (int depth = 1; depth <= radius && !frontier.empty(); depth++) {
      next.clear();
      for (const int u : frontier) {
        for (const int w : adjacency[u]) {
          if (w >= 0 && w < nVertices && dist[w] == -1) {
            dist[w] = depth;
            touched.push_back(w);
            next.push_back(w);
          }
        }
      }
      frontier.swap(next);
    }
    const auto begin = out.size();
    for (const int w : touched) {
      if (w != v) {
        out.emplace_back(w, dist[w]);
      }
    }
    std::sort(out.begin() + static_cast<std::ptrdiff_t>(begin), out.end());
  }

  template <typename Graph> void fillAll(const Graph &adjacency, int radius) {
    std::vector<int> all(adjacency.size());
    std::iota(all.begin(), all.end(), 0);
    fillList(adjacency, radius, all);
  }

  // Same sweep as fillList. clear() keeps the column capacity, so a
  // later frame of the same size does not allocate the table again.
  template <typename Graph> void refillAll(const Graph &adjacency, int radius) {
    fillAll(adjacency, radius);
  }

  template <typename Graph>
  void fillList(const Graph &adjacency, int radius,
                const std::vector<int> &which) {
    const int nVertices = static_cast<int>(adjacency.size());
    const int nWhich = static_cast<int>(which.size());
    off.assign(static_cast<std::size_t>(nVertices) + 1, 0);
    cols.clear();
    if (nWhich == 0) {
      return;
    }
#ifdef SEAMS_HAS_OPENMP
    // Below this size the team costs more than the ball. 1000 molecules
    // on eight threads was the cell that lost to a serial ring code.
    const bool parallel = static_cast<int>(adjacency.size()) >= 2048 &&
                          !omp_in_parallel() && omp_get_max_threads() > 1;
#else
    const bool parallel = false;
#endif
    // Cells are contiguous runs of the vertex list. The run length is the
    // widest neighbour row, which is the coordination of the link-cell list.
    const auto countStart = std::chrono::steady_clock::now();
    int cell = 1;
    for (int v = 0; v < nVertices; v++) {
      int degree = 0;
      for (const int neighbour : adjacency[v]) {
        (void)neighbour;
        ++degree;
      }
      if (degree > cell) {
        cell = degree;
      }
    }
    const int nCells = (nWhich + cell - 1) / cell;
    std::vector<int> width(static_cast<std::size_t>(nVertices), 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (parallel)
#endif
    {
      std::vector<int> dist(nVertices, -1);
      std::vector<int> touched;
      std::vector<int> frontier;
      std::vector<int> next;
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
      for (int c = 0; c < nCells; c++) {
        const int lo = c * cell;
        const int hi = std::min(nWhich, lo + cell);
        for (int k = lo; k < hi; k++) {
          const int v = which[static_cast<std::size_t>(k)];
          for (const int w : touched) {
            dist[w] = -1;
          }
          touched.clear();
          dist[v] = 0;
          touched.push_back(v);
          frontier.assign(1, v);
          for (int depth = 1; depth <= radius && !frontier.empty(); depth++) {
            next.clear();
            for (const int u : frontier) {
              for (const int w : adjacency[u]) {
                if (w >= 0 && w < nVertices && dist[w] == -1) {
                  dist[w] = depth;
                  touched.push_back(w);
                  next.push_back(w);
                }
              }
            }
            frontier.swap(next);
          }
          width[static_cast<std::size_t>(v)] = static_cast<int>(touched.size()) - 1;
        }
      }
    }
    gRingSplit.ball_count_ms = millisSince(countStart);
    const auto fillStart = std::chrono::steady_clock::now();
    std::vector<int> scan(static_cast<std::size_t>(nVertices), 0);
    for (int v = 0; v < nVertices; v++) {
      scan[static_cast<std::size_t>(v)] = width[static_cast<std::size_t>(v)];
    }
    seams::exclusiveScan(scan);
    int total = 0;
    for (int v = 0; v < nVertices; v++) {
      off[static_cast<std::size_t>(v)] = scan[static_cast<std::size_t>(v)];
      total += width[static_cast<std::size_t>(v)];
    }
    off[static_cast<std::size_t>(nVertices)] = total;
    cols.resize(static_cast<std::size_t>(total));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (parallel)
#endif
    {
      std::vector<int> dist(nVertices, -1);
      std::vector<int> touched;
      std::vector<int> frontier;
      std::vector<int> next;
      std::vector<std::pair<int, int>> local;
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
      for (int c = 0; c < nCells; c++) {
        const int lo = c * cell;
        const int hi = std::min(nWhich, lo + cell);
        for (int k = lo; k < hi; k++) {
          const int v = which[static_cast<std::size_t>(k)];
          local.clear();
          appendBall(v, adjacency, radius, dist, touched, frontier, next, local);
          std::copy(local.begin(), local.end(),
                    cols.begin() + static_cast<std::ptrdiff_t>(off[static_cast<std::size_t>(v)]));
        }
      }
    }
    gRingSplit.ball_fill_ms = millisSince(fillStart);
    gRingSplit.balls_ms = gRingSplit.ball_count_ms + gRingSplit.ball_fill_ms;
  }

  // Balls that can change are those of vertices within `radius` of a
  // changed endpoint, in the union of the two frames' graphs. Unchanged
  // rows are identical, so the union is walked as the new row everywhere
  // and the old row as well on a changed vertex.
  template <typename Graph>
  int refreshNear(const Graph &oldAdj, const Graph &newAdj,
                  const std::vector<int> &changed, int radius) {
    const int nVertices = static_cast<int>(newAdj.size());
    if (static_cast<int>(off.size()) != nVertices + 1) {
      refillAll(newAdj, radius);
      return nVertices;
    }
    std::vector<char> isChanged(static_cast<std::size_t>(nVertices), 0);
    std::vector<int> dist(nVertices, -1);
    std::vector<int> frontier, next, affected;
    for (const int v : changed) {
      if (v >= 0 && v < nVertices && dist[v] == -1) {
        dist[v] = 0;
        frontier.push_back(v);
        affected.push_back(v);
        isChanged[static_cast<std::size_t>(v)] = 1;
      }
    }
    for (int depth = 1; depth <= radius && !frontier.empty(); depth++) {
      next.clear();
      for (const int u : frontier) {
        if (isChanged[static_cast<std::size_t>(u)]) {
          for (const int w : oldAdj[u]) {
            if (w >= 0 && w < nVertices && dist[w] == -1) {
              dist[w] = depth;
              affected.push_back(w);
              next.push_back(w);
            }
          }
        }
        for (const int w : newAdj[u]) {
          if (w >= 0 && w < nVertices && dist[w] == -1) {
            dist[w] = depth;
            affected.push_back(w);
            next.push_back(w);
          }
        }
      }
      frontier.swap(next);
    }
    std::vector<int> scratchDist(nVertices, -1);
    std::vector<int> touched, fr, nxt;
    std::vector<std::pair<int, int>> repl;
    std::vector<int> replAt(static_cast<std::size_t>(nVertices), -1);
    std::vector<int> replN(static_cast<std::size_t>(nVertices), 0);
    for (const int v : affected) {
      replAt[static_cast<std::size_t>(v)] = static_cast<int>(repl.size());
      const auto before = repl.size();
      appendBall(v, newAdj, radius, scratchDist, touched, fr, nxt, repl);
      replN[static_cast<std::size_t>(v)] =
          static_cast<int>(repl.size() - before);
    }
    std::vector<std::pair<int, int>> neu;
    neu.reserve(cols.size());
    std::vector<int> newOff(static_cast<std::size_t>(nVertices) + 1, 0);
    for (int v = 0; v < nVertices; v++) {
      newOff[static_cast<std::size_t>(v)] = static_cast<int>(neu.size());
      if (replAt[static_cast<std::size_t>(v)] >= 0) {
        const int b = replAt[static_cast<std::size_t>(v)];
        const int n = replN[static_cast<std::size_t>(v)];
        neu.insert(neu.end(), repl.begin() + b, repl.begin() + b + n);
      } else {
        neu.insert(neu.end(),
                   cols.begin() + off[static_cast<std::size_t>(v)],
                   cols.begin() + off[static_cast<std::size_t>(v + 1)]);
      }
    }
    newOff[static_cast<std::size_t>(nVertices)] = static_cast<int>(neu.size());
    cols.swap(neu);
    off.swap(newOff);
    return static_cast<int>(affected.size());
  }

  //! Graph distance from v to w when within the radius, else a large sentinel
  int dist(int v, int w) const {
    if (off.empty() || v < 0 || v + 1 >= static_cast<int>(off.size())) {
      return 1 << 20;
    }
    const auto b = cols.begin() + off[static_cast<std::size_t>(v)];
    const auto e = cols.begin() + off[static_cast<std::size_t>(v + 1)];
    const auto it = std::lower_bound(b, e, std::make_pair(w, 0));
    if (it != e && it->first == w) {
      return it->second;
    }
    return 1 << 20;
  }
};

//! Bounded breadth-first levels from src, scratch reused across sources
template <typename Graph>
void levelsFrom(const Graph &adjacency, int src, int maxLvl, int minVertex,
                std::vector<int> &lvl, std::vector<int> &touched,
                std::vector<int> &frontier, std::vector<int> &next) {
  for (const int w : touched) {
    lvl[w] = -1;
  }
  touched.clear();
  lvl[src] = 0;
  touched.push_back(src);
  frontier.assign(1, src);
  for (int depth = 1; depth <= maxLvl && !frontier.empty(); depth++) {
    next.clear();
    for (const int u : frontier) {
      for (const int w : adjacency[u]) {
        // Distances are on the whole graph. A shortcut may pass through a
        // vertex below the source without that vertex sitting on the ring.
        // minVertex 0 keeps those shortcuts visible. A positive minVertex
        // is only for a caller that has already restricted the component.
        if (w >= minVertex && w < static_cast<int>(lvl.size()) && lvl[w] == -1) {
          lvl[w] = depth;
          touched.push_back(w);
          next.push_back(w);
        }
      }
    }
    frontier.swap(next);
  }
}

/**
 * @brief Open-addressed ball of radius at most 4, with an occupancy bitset.
 * @details A primitive test on a cycle of length at most 10 asks only
 *  whether the distance is strictly less than the ring distance, and that
 *  distance is at most floor(k/2)-1. Pairs closer than 4 are still answered
 *  from the adjacency rows. Farther pairs probe this table. The occupancy
 *  word is tested before the key, so a miss does not read the slot.
 */
struct BitBalls {
  int radius = 1;
  std::vector<uint32_t> off;
  std::vector<uint32_t> mask;
  std::vector<int> keys;
  std::vector<uint8_t> dists;
  std::vector<uint64_t> occ;
  // Four occupancy words per vertex. A miss is one bit test.
  std::vector<uint64_t> bloom;

  void clear() {
    off.clear();
    mask.clear();
    keys.clear();
    dists.clear();
    occ.clear();
    bloom.clear();
  }

  bool maybe(int u, int v) const {
    const uint32_t h1 = static_cast<uint32_t>(v) * 0x9E3779B1u;
    const uint32_t h2 = static_cast<uint32_t>(v) * 0x85EBCA77u;
    const uint64_t *w = bloom.data() + static_cast<std::size_t>(u) * 4u;
    const bool a = ((w[(h1 >> 6) & 3u] >> (h1 & 63u)) & 1ull) != 0ull;
    const bool b = ((w[(h2 >> 6) & 3u] >> (h2 & 63u)) & 1ull) != 0ull;
    return a && b;
  }

  int lookup(int u, int v) const {
    if (off.empty() || u < 0 || u + 1 >= static_cast<int>(off.size())) {
      return 1 << 20;
    }
    const uint32_t begin = off[static_cast<std::size_t>(u)];
    const uint32_t end = off[static_cast<std::size_t>(u) + 1];
    if (begin == end || !maybe(u, v)) {
      return radius + 1;
    }
    const uint32_t m = mask[static_cast<std::size_t>(u)];
    uint32_t slot = (static_cast<uint32_t>(v) * 0x9E3779B1u) & m;
    for (uint32_t nprobe = 0; nprobe <= m; ++nprobe) {
      const uint32_t at = begin + slot;
      const uint32_t bit = at & 63u;
      if (((occ[at >> 6] >> bit) & 1ull) == 0ull) {
        return radius + 1;
      }
      if (keys[at] == v) {
        return dists[at];
      }
      slot = (slot + 1u) & m;
    }
    return radius + 1;
  }

  template <typename Graph>
  void fill(const Graph &adjacency, int rad, const std::vector<int> &which) {
    radius = rad;
    const int nVertices = static_cast<int>(adjacency.size());
    off.assign(static_cast<std::size_t>(nVertices) + 1, 0);
    mask.assign(static_cast<std::size_t>(nVertices), 0);
    bloom.assign(static_cast<std::size_t>(nVertices) * 4u, 0ull);
    std::vector<uint32_t> slots(static_cast<std::size_t>(nVertices), 0);
    std::vector<char> wanted(static_cast<std::size_t>(nVertices), 0);
    for (const int v : which) {
      if (v >= 0 && v < nVertices) {
        wanted[static_cast<std::size_t>(v)] = 1;
      }
    }
#ifdef SEAMS_HAS_OPENMP
    const bool parallel = nVertices >= 2048 && !omp_in_parallel() &&
                          omp_get_max_threads() > 1;
#pragma omp parallel if (parallel)
#endif
    {
      std::vector<int> dist(static_cast<std::size_t>(nVertices), -1);
      std::vector<int> touched, frontier, next;
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
      for (int v = 0; v < nVertices; v++) {
        if (!wanted[static_cast<std::size_t>(v)]) {
          continue;
        }
        for (const int w : touched) {
          dist[static_cast<std::size_t>(w)] = -1;
        }
        touched.clear();
        dist[static_cast<std::size_t>(v)] = 0;
        touched.push_back(v);
        frontier.assign(1, v);
        for (int depth = 1; depth <= rad && !frontier.empty(); depth++) {
          next.clear();
          for (const int u : frontier) {
            for (const int w : adjacency[u]) {
              if (w >= 0 && w < nVertices && dist[static_cast<std::size_t>(w)] == -1) {
                dist[static_cast<std::size_t>(w)] = depth;
                touched.push_back(w);
                next.push_back(w);
              }
            }
          }
          frontier.swap(next);
        }
        const auto n = static_cast<uint32_t>(touched.size() - 1);
        uint32_t size = 0;
        if (n != 0) {
          size = 64;
          while (size < n * 2u) {
            size <<= 1u;
          }
        }
        slots[static_cast<std::size_t>(v)] = size;
      }
    }
    off[0] = 0;
    for (int v = 0; v < nVertices; v++) {
      mask[static_cast<std::size_t>(v)] =
          slots[static_cast<std::size_t>(v)] == 0 ? 0u : slots[static_cast<std::size_t>(v)] - 1u;
      off[static_cast<std::size_t>(v) + 1] = off[static_cast<std::size_t>(v)] + slots[static_cast<std::size_t>(v)];
    }
    const uint32_t total = off[static_cast<std::size_t>(nVertices)];
    keys.assign(total, -1);
    dists.assign(total, 0);
    occ.assign((static_cast<std::size_t>(total) + 63u) / 64u, 0ull);
    // Each table length is a multiple of 64 and the tables are packed, so
    // occupancy words do not overlap and the inserts run side by side.
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (parallel)
#endif
    {
      std::vector<int> dist(static_cast<std::size_t>(nVertices), -1);
      std::vector<int> touched, frontier, next;
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
      for (int v = 0; v < nVertices; v++) {
        if (!wanted[static_cast<std::size_t>(v)]) {
          continue;
        }
        const uint32_t begin = off[static_cast<std::size_t>(v)];
        const uint32_t m = mask[static_cast<std::size_t>(v)];
        if (off[static_cast<std::size_t>(v) + 1] == begin) {
          continue;
        }
        for (const int w : touched) {
          dist[static_cast<std::size_t>(w)] = -1;
        }
        touched.clear();
        dist[static_cast<std::size_t>(v)] = 0;
        touched.push_back(v);
        frontier.assign(1, v);
        for (int depth = 1; depth <= rad && !frontier.empty(); depth++) {
          next.clear();
          for (const int u : frontier) {
            for (const int w : adjacency[u]) {
              if (w >= 0 && w < nVertices && dist[static_cast<std::size_t>(w)] == -1) {
                dist[static_cast<std::size_t>(w)] = depth;
                touched.push_back(w);
                next.push_back(w);
              }
            }
          }
          frontier.swap(next);
        }
        uint64_t *bloomv = bloom.data() + static_cast<std::size_t>(v) * 4u;
        for (const int w : touched) {
          if (w == v) {
            continue;
          }
          const uint32_t h1 = static_cast<uint32_t>(w) * 0x9E3779B1u;
          const uint32_t h2 = static_cast<uint32_t>(w) * 0x85EBCA77u;
          bloomv[(h1 >> 6) & 3u] |= 1ull << (h1 & 63u);
          bloomv[(h2 >> 6) & 3u] |= 1ull << (h2 & 63u);
          uint32_t slot = h1 & m;
          for (;;) {
            const uint32_t at = begin + slot;
            if (keys[at] < 0) {
              keys[at] = w;
              dists[at] = static_cast<uint8_t>(dist[static_cast<std::size_t>(w)]);
              occ[at >> 6] |= 1ull << (at & 63u);
              break;
            }
            slot = (slot + 1u) & m;
          }
        }
      }
    }
  }
};

/**
 * @brief Distance questions the primitivity test actually asks.
 * @details The antipodal arcs are shortest, so a failing pair has ring
 *  distance at least 2 and at most floor(k/2). A strict shortcut is then a
 *  walk of length at most floor(k/2)-1. Up to length 2 the adjacency rows
 *  answer it. Up to length 4, which covers rings of size 8, 9 and 10, the
 *  occupancy bitset of the ball answers it. The sorted pair table remains
 *  for larger rings and for the baseline. Pairs incident to the source, and
 *  on an even ring the antipode, are not tested again.
 */
struct DistOracle {
  const PackedRows *adj = nullptr;
  const BoundedBalls *balls = nullptr;
  const BitBalls *bits = nullptr;
  bool hop2 = false;

  bool shorter(int u, int v, int ringDist) const {
    if (ringDist <= 1) {
      return false;
    }
    if (hop2 && ringDist <= 3) {
      // A ring distance of 2 is a chord, and a chord is an edge. A row
      // longer than 48 neighbours is cheaper to answer from the ball.
      const auto row = (*adj)[static_cast<std::size_t>(u)];
      int seen = 0;
      bool capped = false;
      for (const int w : row) {
        if (w == v) {
          return true;
        }
        if (++seen == 48) {
          capped = true;
          break;
        }
      }
      if (!capped) {
        if (ringDist <= 2) {
          return false;
        }
        for (const int w : row) {
          for (const int x : (*adj)[static_cast<std::size_t>(w)]) {
            if (x == v) {
              return true;
            }
          }
        }
        return false;
      }
    }
    if (bits != nullptr) {
      return bits->lookup(u, v) < ringDist;
    }
    return balls->dist(u, v) < ringDist;
  }
};

bool ringIsPrimitive(const DistOracle &oracle, const std::vector<int> &ring,
                     int skipA, int skipB) {
  const int k = static_cast<int>(ring.size());
  for (int i = 0; i < k; i++) {
    if (i == skipA || i == skipB) {
      continue;
    }
    for (int j = i + 2; j < k; j++) {
      if (j == skipA || j == skipB) {
        continue;
      }
      const int sep = j - i;
      const int ringDist = std::min(sep, k - sep);
      if (ringDist <= 1) {
        continue;
      }
      if (oracle.shorter(ring[i], ring[j], ringDist)) {
        return false;
      }
    }
  }
  return true;
}

//! Append v if absent; ring sizes are tiny, so a scan beats a set
inline bool pushUniqueMember(std::vector<int> &ring, int v) {
  for (const int x : ring) {
    if (x == v) {
      return false;
    }
  }
  ring.push_back(v);
  return true;
}


/**
 * @brief Reusable scratch for the per-source ring enumeration.
 * @details Vertex v holds pathCount[v] shortest paths from the source, each
 *  of lvl[v] vertices (source excluded, v last), stored back to back in
 *  pathBuf from pathFirst[v].
 */
struct alignas(64) RingScratch {
  std::vector<int> lvl;
  std::vector<int> touched, frontier, next;
  std::vector<int> pathFirst, pathCount, pathBuf;
  std::vector<int> ring;
};

//! Shortest paths from src to every touched vertex, built in breadth-first
//! order as each parent's paths extended by the vertex. Parents are taken in
//! adjacency order, which lists the paths in the order a depth-first walk
//! back along decreasing levels would.
template <typename Graph>
void shortestPathsFrom(const Graph &adjacency, int src, RingScratch &scr) {
  const int nVertices = static_cast<int>(adjacency.size());
  if (static_cast<int>(scr.pathFirst.size()) != nVertices) {
    scr.pathFirst.assign(nVertices, 0);
    scr.pathCount.assign(nVertices, 0);
  }
  scr.pathBuf.clear();
    scr.pathFirst[src] = 0;
  scr.pathCount[src] = 1;
  for (const int v : scr.touched) {
    const int lv = scr.lvl[v];
    if (lv < 1) {
      continue;
    }
    const int first = static_cast<int>(scr.pathBuf.size());
    int count = 0;
    for (const int w : adjacency[v]) {
      if (w < 0 || w >= nVertices || scr.lvl[w] != lv - 1) {
        continue;
      }
      const int plen = lv - 1;
      for (int k = 0; k < scr.pathCount[w]; k++) {
        const int from = scr.pathFirst[w] + k * plen;
        const int dest = static_cast<int>(scr.pathBuf.size());
        scr.pathBuf.resize(scr.pathBuf.size() + static_cast<std::size_t>(lv));
        if (plen > 0) {
          std::copy_n(scr.pathBuf.begin() + from, plen,
                      scr.pathBuf.begin() + dest);
        }
        scr.pathBuf[static_cast<std::size_t>(dest + plen)] = v;
        count++;
      }
    }
    scr.pathFirst[v] = first;
    scr.pathCount[v] = count;
  }
}

/**
 * @details Enumerates every primitive ring whose lowest-indexed member is
 *  @a src, appending them to @a out in canonical form (from the lowest
 *  member toward its lesser neighbour) and sorted. With @a through, every
 *  primitive ring through @a src whose lowest flagged member is @a src
 *  instead, appended unsorted: the incremental updater collects those for
 *  the flagged vertices around a change.
 */
//! Logical prefix of a persistent ring buffer. push_back assigns into storage
//! that already has capacity, so a later frame does not allocate. One pool
//! occupies its own cache line: adjacent sources are often in flight on
//! different threads, and a shared line makes those writes wait on each other.
struct alignas(64) RingPool {
  std::vector<std::vector<int>> storage;
  std::size_t n = 0;

  void reset() { n = 0; }
  void release() {
    std::vector<std::vector<int>>().swap(storage);
    n = 0;
  }
  std::size_t size() const { return n; }

  void push_back(const std::vector<int> &ring) {
    if (n == storage.size()) {
      storage.emplace_back();
    }
    auto &dst = storage[n];
    if (dst.capacity() < ring.size()) {
      dst.reserve(ring.size());
    }
    dst.assign(ring.begin(), ring.end());
    ++n;
  }

  using iterator = std::vector<std::vector<int>>::iterator;
  iterator begin() { return storage.begin(); }
  iterator end() { return storage.begin() + static_cast<std::ptrdiff_t>(n); }
};

template <typename Graph, typename Out>
void enumerateFromSource(const Graph &adjacency, const DistOracle &oracle,
                         int src, int maxDepth, int maxLvl, RingScratch &scr,
                         Out &out, const std::vector<char> *through = nullptr,
                         bool legacyPairs = false) {
  const int nVertices = static_cast<int>(adjacency.size());
  if (static_cast<int>(scr.lvl.size()) != nVertices) {
    scr.lvl.assign(nVertices, -1);
    scr.touched.clear();
  }
  // Distances use every vertex. A shortcut through an index below the
  // source does not put that vertex on the ring, and hiding it made the
  // antipode look farther than the graph distance. Rings are still emitted
  // only from their lowest member: after the distances are known, vertices
  // below the source lose their level so no path is extended through them.
  // The updater passes a through-mask and keeps those levels.
  levelsFrom(adjacency, src, maxLvl, 0, scr.lvl, scr.touched, scr.frontier,
             scr.next);
  if (through == nullptr) {
    for (const int v : scr.touched) {
      if (v < src) {
        scr.lvl[v] = -1;
      }
    }
  }
  shortestPathsFrom(adjacency, src, scr);
  const std::size_t first = out.size();
  auto keep = [&](std::vector<int> &r) {
    if (through != nullptr) {
      for (const int v : r) {
        if (v < src && (*through)[static_cast<std::size_t>(v)]) {
          return;
        }
      }
      std::rotate(r.begin(), std::min_element(r.begin(), r.end()), r.end());
    }
    if (r.size() > 2 && r[1] > r.back()) {
      std::reverse(r.begin() + 1, r.end());
    }
    out.push_back(r);
  };
  std::vector<int> &ring = scr.ring;
  for (const int p : scr.touched) {
    const int lp = scr.lvl[p];
    if (lp < 1 || lp > maxLvl) {
      continue;
    }
    const int *pathsP = scr.pathBuf.data() + scr.pathFirst[p];
    const int nP = scr.pathCount[p];
    // Even rings: two vertex-disjoint shortest paths to an antipodal vertex
    if (2 * lp >= 3 && 2 * lp <= maxDepth) {
      for (int a = 0; a < nP; a++) {
        for (int b = a + 1; b < nP; b++) {
          const int *pa = pathsP + a * lp;
          const int *pb = pathsP + b * lp;
          ring.assign(1, src);
          bool ok = true;
          for (int i = 0; i < lp && ok; i++) {
            ok = pushUniqueMember(ring, pa[i]);
          }
          for (int i = lp - 2; i >= 0 && ok; i--) {
            ok = pushUniqueMember(ring, pb[i]);
          }
          if (ok && ringIsPrimitive(oracle, ring, legacyPairs ? -1 : 0,
                                    legacyPairs ? -1 : lp)) {
            keep(ring);
          }
        }
      }
    }
    // Odd rings: an antipodal edge between two vertices at the same level
    if (2 * lp + 1 <= maxDepth) {
      for (const int q : adjacency[p]) {
        if (q <= p || q >= nVertices || scr.lvl[q] != lp) {
          continue;
        }
        const int *pathsQ = scr.pathBuf.data() + scr.pathFirst[q];
        for (int a = 0; a < nP; a++) {
          for (int b = 0; b < scr.pathCount[q]; b++) {
            const int *pa = pathsP + a * lp;
            const int *pb = pathsQ + b * lp;
            ring.assign(1, src);
            bool ok = true;
            for (int i = 0; i < lp && ok; i++) {
              ok = pushUniqueMember(ring, pa[i]);
            }
            for (int i = lp - 1; i >= 0 && ok; i--) {
              ok = pushUniqueMember(ring, pb[i]);
            }
            if (ok && ringIsPrimitive(oracle, ring, legacyPairs ? -1 : 0, -1)) {
              keep(ring);
            }
          }
        }
      }
    }
  }
  if (through == nullptr) {
    std::sort(out.begin() + static_cast<std::ptrdiff_t>(first), out.end());
  }
}

} // namespace

/**
 * @details Returns the primitive (shortest-path) rings, by index, given a
 *  neighbour list (also by index) and the largest ring size sought.
 *
 *  The construction follows Yuan and Cormack (Comput. Mater. Sci. 24, 343,
 *  2002): every member of a primitive ring sits at a graph distance from any
 *  other member equal to their separation around the ring, so each arc of the
 *  ring is a shortest path from any member taken as a source. A ring of even
 *  size @f$2L@f$ is therefore a source vertex plus two vertex-disjoint
 *  shortest paths to a vertex at level @f$L@f$; a ring of odd size
 *  @f$2L+1@f$ closes instead over an edge joining two vertices both at level
 *  @f$L@f$. Candidates are built only from shortest paths, and each ring is
 *  enumerated exactly once by restricting the search to its lowest-indexed
 *  member, so the output carries no duplicates by construction.
 *
 *  Candidates are then kept or dropped by the same criterion Franzblau's
 *  filter applied: no pair of members may be joined through the graph by
 *  fewer edges than around the ring. On a cycle of length at most 7 that
 *  question is a walk of one or two edges, read from the adjacency rows.
 *  Longer cycles use one bounded neighbourhood per vertex. Pairs that sit
 *  on a shortest path out of the source, and on an even cycle the antipode,
 *  are not tested again.
 *
 *  Each ring starts at its lowest-indexed member and runs toward the lesser
 *  of that member's two ring neighbours, and the rings are sorted
 *  lexicographically: the output depends on the graph alone, not on the
 *  order of the neighbour rows.
 *
 *  primitive::countAllRingsFromIndex and primitive::removeNonSPrings remain
 *  available and produce the same ring set by the generate-then-filter route.
 * @param[in] nList Row-ordered neighbour list by index (and NOT the atom ID)
 * @param[in] maxDepth The largest ring size to search for. Rings larger than
 *  maxDepth are not generated.
 * @return A vector of vectors of the rings; each ring contains the atom
 *  indices of the ring members.
 */
static std::vector<std::vector<int>>
searchRings(const std::vector<std::vector<int>> &nList, int maxDepth,
            const std::vector<char> &sources, bool legacy);

std::vector<std::vector<int>>
primitive::ringNetwork(const std::vector<std::vector<int>> &nList, int maxDepth) {
  return searchRings(nList, maxDepth, std::vector<char>(), false);
}

std::vector<std::vector<int>>
primitive::ringNetwork(const std::vector<int> &offsets, const std::vector<int> &cols,
                       int maxDepth) {
  const int n = offsets.empty() ? 0 : static_cast<int>(offsets.size()) - 1;
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(std::max(n, 0)));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
  for (int i = 0; i < n; ++i) {
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    auto &row = rows[static_cast<std::size_t>(i)];
    row.resize(static_cast<std::size_t>(1 + b - a));
    row[0] = i;
    for (int k = a; k < b; ++k) {
      row[static_cast<std::size_t>(1 + k - a)] = cols[static_cast<std::size_t>(k)];
    }
  }
  return searchRings(rows, maxDepth, std::vector<char>(), false);
}

std::vector<std::vector<int>>
primitive::ringNetworkBalls(const std::vector<std::vector<int>> &nList, int maxDepth) {
  return searchRings(nList, maxDepth, std::vector<char>(), true);
}

/**
 * @details The rings of the whole search whose lowest-indexed member is
 *  flagged in @a sources, in the same order; an empty mask flags every
 *  vertex. Balls are filled only within maxDepth/2 hops of a flagged
 *  vertex, where the members of its rings lie.
 */
static std::vector<std::vector<int>>
searchRings(const std::vector<std::vector<int>> &nList, int maxDepth,
            const std::vector<char> &sources, bool legacy) {
  std::vector<std::vector<int>> rings;
  if (maxDepth < 3 || nList.empty()) {
    return rings;
  }
  const int nVertices = static_cast<int>(nList.size());
  // The fork, the per-source pools, and pausing the team are slower than
  // the search itself on a small graph. Larger graphs keep the team.
  const bool serial = nVertices < 2048;
  if (!serial) {
    seams::adoptThreads();
  }

  // Adjacency by index; the first element of each nList row is the vertex
  // itself, so it is skipped. A repeated neighbour is one edge. A self
  // entry and an index outside the cloud are not edges.
  const auto adjStart = std::chrono::steady_clock::now();
  std::vector<std::vector<int>> uniq(static_cast<std::size_t>(nVertices));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (!serial && !omp_in_parallel())
#endif
  for (int i = 0; i < nVertices; i++) {
    if (nList[static_cast<std::size_t>(i)].size() <= 1) {
      continue;
    }
    auto &row = uniq[static_cast<std::size_t>(i)];
    row.reserve(nList[static_cast<std::size_t>(i)].size() - 1);
    for (std::size_t k = 1; k < nList[static_cast<std::size_t>(i)].size(); k++) {
      const int j = nList[static_cast<std::size_t>(i)][k];
      if (j < 0 || j >= nVertices || j == i) {
        continue;
      }
      if (std::find(row.begin(), row.end(), j) == row.end()) {
        row.push_back(j);
      }
    }
  }
  PackedRows adjacency;
  adjacency.first.assign(nVertices + 1, 0);
  for (int i = 0; i < nVertices; i++) {
    adjacency.first[i + 1] =
        adjacency.first[i] +
        static_cast<int>(uniq[static_cast<std::size_t>(i)].size());
  }
  adjacency.cols.resize(adjacency.first[nVertices]);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (!serial && !omp_in_parallel())
#endif
  for (int i = 0; i < nVertices; i++) {
    const auto &row = uniq[static_cast<std::size_t>(i)];
    if (!row.empty()) {
      std::copy(row.begin(), row.end(),
                adjacency.cols.begin() +
                    static_cast<std::ptrdiff_t>(adjacency.first[i]));
    }
  }

  const int maxLvl = maxDepth / 2;
  const int radius = std::max(maxLvl - 1, 1);
  // Ring distance on a cycle of length <= maxDepth is at most maxLvl.
  // A strict shortcut is then a walk of length <= maxLvl - 1. When that
  // bound is 2 (maxDepth <= 7), adjacency answers every query the
  // primitivity test asks, and the ball table is dead work.
  const bool hopOnly = !legacy && maxLvl <= 3;
  // Rings of size 8, 9 and 10 need distances up to 3 or 4. The occupancy
  // bitset of that ball replaces the sorted pair table. Larger rings keep
  // the table.
  const bool bitBall = !legacy && maxLvl > 3 && maxLvl <= 5;
  std::vector<int> from;
  for (int v = 0; v < nVertices; v++) {
    if (sources.empty() || (v < static_cast<int>(sources.size()) && sources[v])) {
      from.push_back(v);
    }
  }
  BoundedBalls balls;
  BitBalls bits;
  gRingSplit.adj_ms = millisSince(adjStart);
  const auto ballStart = std::chrono::steady_clock::now();
  if (hopOnly) {
    gUsedBalls = 0;
    gRingSplit.ball_count_ms = 0;
    gRingSplit.ball_fill_ms = 0;
    gRingSplit.balls_ms = 0;
  } else if (bitBall) {
    gUsedBalls = 2;
    std::vector<int> reach = from;
    if (static_cast<int>(from.size()) != nVertices) {
      std::vector<int> lvl(nVertices, -1);
      for (const int v : from) {
        lvl[v] = 0;
      }
      for (std::size_t k = 0; k < reach.size(); k++) {
        const int u = reach[k];
        if (lvl[u] == maxLvl) {
          continue;
        }
        for (const int w : adjacency[u]) {
          if (w >= 0 && w < nVertices && lvl[w] == -1) {
            lvl[w] = lvl[u] + 1;
            reach.push_back(w);
          }
        }
      }
    }
    bits.fill(adjacency, radius, reach);
    gRingSplit.ball_count_ms = 0;
    gRingSplit.ball_fill_ms = millisSince(ballStart);
    gRingSplit.balls_ms = gRingSplit.ball_fill_ms;
  } else {
    gUsedBalls = 1;
    if (static_cast<int>(from.size()) == nVertices) {
      balls.fillList(adjacency, radius, from);
    } else {
      std::vector<int> lvl(nVertices, -1);
      std::vector<int> reached = from;
      for (const int v : from) {
        lvl[v] = 0;
      }
      for (std::size_t k = 0; k < reached.size(); k++) {
        const int u = reached[k];
        if (lvl[u] == maxLvl) {
          continue;
        }
        for (const int w : adjacency[u]) {
          if (w >= 0 && w < nVertices && lvl[w] == -1) {
            lvl[w] = lvl[u] + 1;
            reached.push_back(w);
          }
        }
      }
      balls.fillList(adjacency, radius, reached);
    }
    gRingSplit.balls_ms = millisSince(ballStart);
  }
  DistOracle oracle;
  oracle.adj = &adjacency;
  oracle.balls = (hopOnly || bitBall) ? nullptr : &balls;
  oracle.bits = bitBall ? &bits : nullptr;
  oracle.hop2 = hopOnly || bitBall;
  const int nFrom = static_cast<int>(from.size());

#ifdef SEAMS_HAS_OPENMP
  // One thread writes straight into the output vector. The parallel
  // collector gives every source its own ring list, and on one thread
  // that list of lists costs more than the search. The directing step
  // clears levels only inside the ball this source touched, so the
  // search stays linear in the number of vertices.
  const bool parallelSources =
      !serial && !omp_in_parallel() && omp_get_max_threads() > 1;
  if (parallelSources) {
    // Sources are independent given the shared read-only balls; per-source
    // collection keeps the output identical to the serial ascending-source
    // order
    // One cache line per source. The ring vectors of adjacent sources
    // used to share a line, and the threads then waited on each other's
    // size updates. A static share of the sources is even on a crystal
    // and keeps each thread on its own lines.
    const auto enumStart = std::chrono::steady_clock::now();
    std::vector<RingPool> perSource(static_cast<std::size_t>(nFrom));
#pragma omp parallel
    {
      RingScratch scratch;
#pragma omp for schedule(guided)
      for (int k = 0; k < nFrom; k++) {
        enumerateFromSource(adjacency, oracle, from[k], maxDepth, maxLvl,
                            scratch, perSource[static_cast<std::size_t>(k)],
                            nullptr, legacy);
      }
    }
    gRingSplit.enumerate_ms = millisSince(enumStart);
    const auto mergeStart = std::chrono::steady_clock::now();
    // Every ring, source list and ball is a heap block of its own; moving
    // and freeing them on one thread is a sizeable share of a threaded search
    std::vector<std::size_t> offset(static_cast<std::size_t>(nFrom) + 1, 0);
    for (int k = 0; k < nFrom; k++) {
      offset[static_cast<std::size_t>(k + 1)] =
          offset[static_cast<std::size_t>(k)] +
          perSource[static_cast<std::size_t>(k)].size();
    }
    rings.resize(offset[static_cast<std::size_t>(nFrom)]);
#pragma omp parallel
    {
#pragma omp for schedule(static) nowait
      for (int k = 0; k < nFrom; k++) {
        auto &pool = perSource[static_cast<std::size_t>(k)];
        std::move(pool.begin(), pool.end(),
                  rings.begin() +
                      static_cast<std::ptrdiff_t>(offset[static_cast<std::size_t>(k)]));
        pool.release();
      }
    }
    balls.clear();
    std::vector<std::vector<int>>().swap(uniq);
    std::vector<int>().swap(adjacency.cols);
    std::vector<std::size_t>().swap(adjacency.first);
    std::vector<RingPool>().swap(perSource);
    gRingSplit.merge_ms = millisSince(mergeStart);
    return rings;
  }
#endif

  const auto enumStart = std::chrono::steady_clock::now();
  RingScratch scratch;
  for (const int src : from) {
    enumerateFromSource(adjacency, oracle, src, maxDepth, maxLvl, scratch,
                        rings, nullptr, legacy);
  }
  gRingSplit.enumerate_ms = millisSince(enumStart);
  const auto mergeStart = std::chrono::steady_clock::now();
  balls.clear();
  std::vector<std::vector<int>>().swap(uniq);
  std::vector<int>().swap(adjacency.cols);
  std::vector<std::size_t>().swap(adjacency.first);
  gRingSplit.merge_ms = millisSince(mergeStart);

  return rings;
}

std::vector<std::vector<int>>
primitive::ringNetwork(const std::vector<std::vector<int>> &nList, int maxDepth,
                       const std::vector<char> &sources) {
  return searchRings(nList, maxDepth, sources, false);
}

primitive::RingSearchSplit primitive::lastRingSearchSplit() {
  return {gRingSplit.balls_ms, gRingSplit.enumerate_ms, gRingSplit.merge_ms,
          gRingSplit.adj_ms, gRingSplit.ball_count_ms, gRingSplit.ball_fill_ms};
}

int primitive::lastRingSearchUsedBalls() { return gUsedBalls; }

/**
 *  @details Get all possible rings (only atom indices, not IDs). The input
 *   neighbour list is in terms of indices. All possible rings (including
 * non-SP) rings are generated using a recursive backtracking algorithm. This
 * function calls the following functions internally:
 *   - primitive::populateGraphFromIndices (for initializing the Graph object)
 *   - primitive::findRings (for getting all rings, using backtracking)
 *   - primitive::restoreEdgesFromIndices (restores back-up of the edges since
 * some may have been removed)
 *  @param[in] neighHbondList Row-ordered neighbour list by atom index (not ID).
 *  @param[in] maxDepth The maximum depth upto which rings will be searched.
 *   This means that rings larger than maxDepth will not be generated.
 *  @return The Graph object for the current frame.
 */
primitive::Graph
primitive::countAllRingsFromIndex(const std::vector<std::vector<int>> &neighHbondList,
                                  int maxDepth) {
  //
  primitive::Graph fullGraph; // Graph object
  std::vector<int>
      visited; // Contains the indices (NOT IDs) visited for a particular node
  int depth;   // Current depth

  // ------------------------------
  // Init
  // Initialize the graph object with all the information from the neighbour
  // list (of indices only)
  fullGraph = primitive::populateGraphFromIndices(neighHbondList);
  // ------------------------------
  // Loop through every vertex
  for (int iatom = 0; iatom < neighHbondList.size(); iatom++) {
    visited.clear();
    depth = 0;
    // Call the function for getting rings
    primitive::findRings(fullGraph, iatom, visited, maxDepth, depth);
  } // loop through every vertex
  // ------------------------------
  // Restore back-up of the edges (some may have been removed)
  primitive::restoreEdgesFromIndices(fullGraph, neighHbondList);
  // ------------------------------

  return fullGraph;
}

/**
 * @details All possible rings are searched for in this function, which
 * recursively calls itself. The rings are 'grown' from the root node (which is
 * the first vertex) using the backtracking algorithm. When it is first called
 * (before the root node has been assigned), root is a dummy value (which is
 * equal to -1, a value that can never be legitimate).
 *  @param[in, out] fullGraph Graph object containing the vertices (and the
 * neighbour lists). Vertices may be 'removed' from the Graph.
 *  @param[in] v The current vertex being visited or checked. It is added to the
 * list of all vertices visited.
 *  @param[in] visited A vector containing a list of the vertices visited for
 * book-keeping. If the current visited vector fulfills the condition for being
 * a ring, it is added to the rings vector of vector in the Graph.
 *  @param[in] maxDepth The maximum depth upto which rings will be searched.
 * This means that rings larger than maxDepth will not be generated.
 *  @param[in] depth The current depth. When this function is called for the
 * first time from primitive::countAllRingsFromIndex, the depth is initialized
 * to 0. When the depth is greater than or equal to maxDepth, the function
 * exits.
 *  @param[in] root The first vertex, from which the current visited vector
 * (candidate ring) is being grown. This is initialized to a dummy value of -1,
 * when it is called from primitive::countAllRingsFromIndex.
 */
int primitive::findRings(Graph &fullGraph, int v, std::vector<int> &visited,
                         int maxDepth, int depth, int root) {
  //
  int nnumNeighbours; // Number of nearest neighbours for iNode
  int n;              // Node of a neighbour (index)
  // -------------
  // For the first call:
  if (root == -1) {
    root = v;                          // Init the root as the current node
    fullGraph.pts[v].inGraph = false; // Remove the root from the graph
  }                                    // first call
  //
  // Checks
  if (depth >= maxDepth) {
    return 1; // false?
  }           // searched until the maximum length specified
  //
  // Add the current node to the visited vector
  visited.push_back(v);
  // -------------
  // Depth-first search
  // 1. Pick a root node (above)
  // 2. Search all the neighbouring nodes
  // 3. Start a recursion call, from the second nearest neighbours
  // 3(i) If the neighbour equals the root (selected in 1), ring has been found!
  // 3(ii) Or else search for all the unvisited neighbours
  // 4. Remove the root node from the graph
  // 5. Update the vector of vectors containing all rings
  // -------------
  // Start the search
  depth += 1; // Go to the next layer
  nnumNeighbours = fullGraph.pts[v].neighListIndex.size();
  // Loop through the neighbours of iNode
  for (int j = 0; j < nnumNeighbours; j++) {
    n = fullGraph.pts[v].neighListIndex[j]; // Neighbour index
    // Has a ring been found?!
    if (depth > 2 && n == root) {
      // Add the visited vector to the rings vector of vector
      fullGraph.rings.push_back(visited);
    } // A ring has been found!
    // Otherwise search all the neighbours which have not been searched already
    else if (fullGraph.pts[n].inGraph) {
      fullGraph.pts[n].inGraph = false; // Set to false now
      // Recursive call
      primitive::findRings(fullGraph, n, visited, maxDepth, depth, root);
      fullGraph.pts[n].inGraph = true; // Put n back in the graph
    } // Search the other unsearched neighbours
  }   // end of search through all neighbours
  // -------------
  // When the depth is 2, remove just the root from the neighbours of iNode
  if (depth == 2) {
    // Search for root in the neighbour list of v
    for (int j = 0; j < fullGraph.pts[v].neighListIndex.size(); j++) {
      if (root == fullGraph.pts[v].neighListIndex[j]) {
        fullGraph.pts[v].neighListIndex.erase(
            fullGraph.pts[v].neighListIndex.begin() + j);
      } // end of erase
    }   // end of search
  }     // remove root not edges

  //
  visited.pop_back();
  //
  return 0;
}

/**
 * @details Fills a Graph object with information from the PointCloud and the
 *  neighbour list. The indices in the neighbour list in the Vertex object are
 *  NOT the atom IDs (they are the atom indices according to the input
 *  PointCloud). The input neighbour list is by atom ID.
 * @param[in] yCloud The input PointCloud.
 * @param[in] neighHbondList The row-ordered neighbour list, containing atom
 *  IDs, and not the atom indices.
 * @return The Graph object for the current frame.
 */
primitive::Graph primitive::populateGraphFromNListID(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &neighHbondList) {
  //
  primitive::Graph fullGraph; // Contains all the information of the pointCloud
  primitive::Vertex iVertex;  // The vertex corresponding to a particular point
  int nnumNeighbours;         // Number of nearest neighbours for iatom
  std::vector<int> iNeigh;    // Neighbours of the current iatom
  int jatomID;                // Atom ID of the nearest neighbour
  int jatomIndex;             // Atom index of the nearest neighbour
  // ------------------------------
  // Loop through every point in yCloud
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // -----
    // Update the neighbours of iatom (only the indices!)
    iNeigh.clear(); // init
    nnumNeighbours = neighHbondList[iatom].size() -
                     1; // The first element is atomID of iatom
    for (int j = 1; j <= nnumNeighbours; j++) {
      jatomID = neighHbondList[iatom][j]; // Atom ID
      // Get the atom index for the vector nearest neighbour list
      auto it = yCloud.idIndexMap.find(jatomID);
      if (it != yCloud.idIndexMap.end()) {
        jatomIndex = it->second;
      } // found jatomIndex
      else {
        // A missing ID used to return the vertices built so far. Callers
        // then treated that prefix as the whole graph.
        std::cerr << "Skipping neighbour ID " << jatomID
                  << " which is not in the frame.\n";
        continue;
      } // jatomIndex not found
      // Update iNeigh
      iNeigh.push_back(jatomIndex);
    } // end of loop through nearest neighbours
    // -----
    // Update iVertex
    iVertex.atomIndex = iatom;
    iVertex.neighListIndex = iNeigh;
    // Add to the Graph object
    fullGraph.pts.push_back(iVertex);
  } // end of loop through every iatom
  // ------------------------------

  return fullGraph;
}

/**
 * @details Fills a Graph object with information from the PointCloud and the
 *  neighbour list. The indices in the neighbour list in the Vertex object are
 *  NOT the atom IDs (they are the atom indices according to the input
 *  PointCloud). The input neighbour list is by index NOT atom IDs. Otherwise,
 *  this function does the same thing as primitive::populateGraphFromNListID.
 * The only difference is that this function takes the neighbour list BY INDEX.
 * @param[in] nList The row-ordered neighbour list, containing atom
 *  indices (according to the input PointCloud).
 * @return The Graph object for the current frame.
 */
primitive::Graph
primitive::populateGraphFromIndices(const std::vector<std::vector<int>> &nList) {
  //
  primitive::Graph fullGraph; // Contains all the information of the pointCloud
  primitive::Vertex iVertex;  // The vertex corresponding to a particular point
  int nnumNeighbours;         // Number of nearest neighbours for iatom
  int iatom, jatom;           // Atom index being saved
  // ------------------------------
  // Loop through every point in nList
  for (int i = 0; i < nList.size(); i++) {
    iatom = nList[i][0]; // Atom index of i
    // neighListIndex is simply the i^th row of nList
    //
    // Update iVertex
    iVertex.atomIndex = iatom;
    iVertex.neighListIndex = nList[i];
    // Add to the Graph object
    fullGraph.pts.push_back(iVertex);
  } // end of loop through iatom

  return fullGraph;
}

/**
 * @details Re-fills the neighbour lists of a Graph object from a row-ordered
 *  neighbour list (which is BY INDEX not IDs). Some vertices may have been
 *  removed while rings were generated using the backtracking algorithm
 *  (primitive::findRings). Also, the indices in the neighbour list in the
 *  Vertex object are not the atom IDs.
 * @param[in] fullGraph The Graph object for the current frame. The neighbour
 *  lists of component Vertex objects may have been depleted.
 * @param[in] nList The row-ordered neighbour list, containing atom
 *  indices (according to the input PointCloud).
 * @return The Graph object for the current frame.
 */
void
primitive::restoreEdgesFromIndices(Graph &fullGraph,
                                   const std::vector<std::vector<int>> &nList) {
  //
  // ------------------------------
  // Loop through every point in nList
  for (int i = 0; i < nList.size(); i++) {
    // neighListIndex is simply the i^th row of nList
    // Update the neighListIndex list in the graph object
    fullGraph.pts[i].neighListIndex = nList[i];
  } // end of loop through iatom
}

/**
 * @details Removes non-SP rings (according to the Franzblau shortest path
 *  criterion) from the rings vector of vectors member n the Graph object. The
 *  rings vector of vectors contains indices and NOT atom IDs. This function
 *  calls primitive::shortestPath internally to calculate the shortest path.
 * @param[in] fullGraph The Graph object for the current frame. This also
 *  contains the rings vector of vectors, which has all possible rings (possibly
 *  inclding non-SP rings).
 * @return The Graph object for the current frame.
 */
void primitive::removeNonSPrings(primitive::Graph &fullGraph) {
  //
  const int nVertices = fullGraph.pts.size(); // Number of vertices in the graph
  const int nRings = fullGraph.rings.size();  // Number of rings
  std::vector<bool> ringsToRemove; // Vector containing the logical values for
                                   // removal of the current ring index
  std::vector<std::vector<int>>
      emptyTempRings; // Empty vector of vectors to swap
  std::vector<std::vector<int>>
      primitiveRings; // Vector of vectors of rings after removing non SP rings
  // -------------------
  // Make sure all the vertices are in the graph before removing non-SP rings
  for (int iVer = 0; iVer < nVertices; iVer++) {
    fullGraph.pts[iVer].inGraph = true;
  } // end of loop through every vertex
  // -------------------
  // The criterion compares, for each pair of non-adjacent ring members, the
  // separation around the ring against the separation through the graph. Only
  // the latter needs searching, it is bounded by half the longest ring, and it
  // does not depend on which ring the pair was drawn from.
  //
  // Collect every such question first and group them by the vertex they start
  // from. One breadth-first sweep per distinct starting vertex then answers all
  // of its questions at once, in place of one depth-first search per pair per
  // ring. The sweep reaches only the handful of vertices within the bound, so
  // its scratch is cleared by walking what it touched rather than the whole
  // graph.
  ringsToRemove.assign(nRings, false);

  struct Query {
    int target;   // The other ring member
    int ringHops; // Separation around the ring, in edges
    int ring;     // Which ring asked
  };
  std::vector<std::vector<Query>> queries(nVertices);

  int maxHops = 0;
  for (int iRing = 0; iRing < nRings; iRing++) {
    const std::vector<int> &currentRing = fullGraph.rings[iRing];
    const int ringSize = currentRing.size();
    maxHops = std::max(maxHops, ringSize / 2);
    for (int jVer = 0; jVer < ringSize; jVer++) {
      // connect with all other, skip j-j (distance=0) and j-(j+1) (nearest
      // neighbors)
      for (int kVer = jVer + 2; kVer < ringSize; kVer++) {
        const int currentV = currentRing[jVer];
        const int currentN = currentRing[kVer];
        if (currentV < 0 || currentV >= nVertices || currentN < 0 ||
            currentN >= nVertices) {
          continue;
        }
        const int d_jk = std::abs(jVer - kVer);
        queries[currentV].push_back(
            {currentN, std::min(d_jk, std::abs(d_jk - ringSize)), iRing});
      } // end of loop through k^th vertex
    }   // end of loop through j^th vertex
  }     // end of loop through rings

  // Answer each vertex's questions with one bounded sweep
  std::vector<int> distance(nVertices, -1);
  std::vector<int> touched;
  std::vector<int> frontier;
  std::vector<int> next;

  for (int v = 0; v < nVertices; v++) {
    if (queries[v].empty()) {
      continue;
    }

    touched.clear();
    frontier.assign(1, v);
    distance[v] = 0;
    touched.push_back(v);
    for (int depth = 1; depth <= maxHops && !frontier.empty(); depth++) {
      next.clear();
      for (const int u : frontier) {
        for (const int w : fullGraph.pts[u].neighListIndex) {
          if (w >= 0 && w < nVertices && distance[w] == -1) {
            distance[w] = depth;
            touched.push_back(w);
            next.push_back(w);
          }
        }
      }
      frontier.swap(next);
    }

    for (const Query &q : queries[v]) {
      // Separation through the graph, in edges; -1 when beyond the bound
      const int graphHops = distance[q.target];
      if (graphHops >= 0 && graphHops < q.ringHops) {
        ringsToRemove[q.ring] = true;
      } // Decide whether to keep or remove the ring
    }

    // Clear only what this sweep reached
    for (const int w : touched) {
      distance[w] = -1;
    }
  } // end of loop through vertices
  // -------------------
  // Remove all the rings whose indices are given in the ringsToRemove vector
  for (int i = 0; i < ringsToRemove.size(); i++) {
    if (!ringsToRemove[i]) {
      primitiveRings.push_back(fullGraph.rings[i]);
    } // updates new copy
  }   // end of loop through ringsToRemove
  // -------------------
  // Update the graph rings with the primitiveRings
  emptyTempRings.swap(fullGraph.rings);
  fullGraph.rings = std::move(primitiveRings);
  // -------------------
}

/**
 * @details Calculates the shortest path for a particular ring. This function
 * uses recursion.
 * @param[in] fullGraph The Graph object for the current frame.
 * @param[in] v The current vertex being checked.
 * @param[in] goal The first element of the candidate ring being checked (the
 *  root node).
 * @param[in] path The path or length of the visited points (This basically
 *  contains all the indices in the visited vector, excluding the current
 * vertex).
 * @param[in] visited This vector contains the indices of the vertices visited
 *  or checked (for book-keeping).
 * @param[in] maxDepth The maximum depth or maximum length of the rings.
 * @param[in] depth The current depth. When this function is called from
 *  primitive::removeNonSPrings, the depth is initialized as 0.
 * @return The Graph object for the current frame.
 */
int primitive::shortestPath(Graph &fullGraph, int v, int goal,
                            std::vector<int> &path, std::vector<int> &visited,
                            int maxDepth, int depth) {
  int len_path = 0;   // Length of the path
  int nnumNeighbours; // Number of neighbours for a particular v
  int n;              // Index of the nearest neighbour of v
  // Start the search for the shortest path
  if (depth < maxDepth) {
    depth += 1; // One layer below
    visited.push_back(
        v); // Add the current vertex to the path (visited points)
    //
    if (v == goal) {
      len_path = path.size(); // Path of the length of visited points
      // If the current path is shorter OR this is the first path found
      if (depth < len_path || len_path == 0) {
        path = visited;
        maxDepth = depth;
      } // Current path is the shortest
    }   // Goal reached
    // Recursive calls to function
    else {
      nnumNeighbours = fullGraph.pts[v].neighListIndex.size();
      // Search all the neighbours
      for (int j = 0; j < nnumNeighbours; j++) {
        n = fullGraph.pts[v].neighListIndex[j]; // Index of nearest neighbour
        // If n has not already been searched:
        if (fullGraph.pts[n].inGraph == true) {
          fullGraph.pts[n].inGraph = false; // Set to false
          primitive::shortestPath(fullGraph, n, goal, path, visited, maxDepth,
                                  depth);   // Recursive call
          fullGraph.pts[n].inGraph = true; // Put back in the graph
        } // If n is in the graph, call recursively
      }   // End of loop over all neighbours
    }     // Goal not reached
    //
    // Pop the vector
    visited.pop_back();
  } // for depth less than maxDepth
  return 0;
}

/**
 * @details Function for clearing Graph if it is already
 *  filled. This should be called before every frame is read in.
 * @param[out] currentGraph The cleared Graph
 */
primitive::Graph primitive::clearGraph(Graph &currentGraph) {
  //
  std::vector<primitive::Vertex> tempPts;
  std::vector<std::vector<int>> tempRings;
  tempPts.swap(currentGraph.pts);
  tempRings.swap(currentGraph.rings);
  return currentGraph;
}

/**
 * @details Implementation state of primitive::RingUpdater.
 *
 *  The locality bound that makes the update exact. Let C hold the vertices
 *  whose adjacency row differs between the old graph G and the new graph
 *  G', the endpoints of every added or removed edge. A cycle of at most
 *  maxDepth members is a primitive ring of one graph and not of the other
 *  only if it holds a changed edge, whose endpoints are in C, or if a
 *  changed edge lies on a shortcut: a path of l < d edges between two
 *  members d apart around the ring. There d <= maxLvl = maxDepth/2, and one
 *  endpoint of the changed edge lies within (l-1)/2 <= reach = (maxLvl-2)/2
 *  edges of a member along the path. Every ring that differs therefore
 *  meets T, the vertices within reach of C in the union of G and G'; reach
 *  is 0 up to maxDepth 7, where T is C. The update drops the stored rings
 *  that meet T and adds the rings of G' that meet T, each enumerated once,
 *  from its lowest-indexed member in T. Rings are held in the canonical
 *  form and sorted order of ringNetwork, so the result equals a full
 *  recomputation element by element. The balls of a vertex change only when
 *  an edge within ballRadius = maxLvl-1 of it changes, so they are refreshed
 *  around C on the union graph.
 */
struct primitive::RingUpdater::Impl {
  int maxDepth = 0;
  int maxLvl = 0;
  int ballRadius = 1;
  int reach = 0;
  // Two packed, sorted adjacency buffers. `which` is the previous frame.
  // One allocation holds every row, and frame 2 onward reuses it.
  PackedRows graph[2];
  int which = 0;
  bool primed = false;
  std::vector<std::vector<int>> flattened;
  bool ballsStale = true;
  BoundedBalls balls;
  std::vector<RingPool> bySource;
  std::vector<RingScratch> scratches;
  int recomputed = 0;
  int ballsRefreshed = 0;
  int changedBonds = 0;
  double changedBondFraction = 0.0;
  bool bondsPending = false;
  int pendingOld = 0;
  bool slotSorted[2] = {false, false};
  bool hot = false;
  bool hopOnly = false;
  // The default does not use a fraction. After T is built, a full pass
  // runs only when T is the whole vertex set.
  double fullPassFraction = 1.0;
  bool useVertexFraction = false;
  // A negative value leaves the bond-count test off.
  double bondCrossover = -1.0;

  void writeGraph(const std::vector<std::vector<int>> &nList, int slot,
                  bool sortRows) {
    auto &dst = graph[slot];
    const int nVertices = static_cast<int>(nList.size());
    dst.first.resize(static_cast<std::size_t>(nVertices) + 1);
    dst.first[0] = 0;
    for (int v = 0; v < nVertices; v++) {
      const auto n =
          nList[static_cast<std::size_t>(v)].empty()
              ? 0
              : nList[static_cast<std::size_t>(v)].size() - 1;
      dst.first[static_cast<std::size_t>(v) + 1] =
          dst.first[static_cast<std::size_t>(v)] + n;
    }
    dst.cols.resize(dst.first[static_cast<std::size_t>(nVertices)]);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nVertices >= 1024 && !omp_in_parallel())
#endif
    for (int v = 0; v < nVertices; v++) {
      const auto n = dst.first[static_cast<std::size_t>(v) + 1] -
                     dst.first[static_cast<std::size_t>(v)];
      if (n == 0) {
        continue;
      }
      int *row = dst.cols.data() +
                 static_cast<std::ptrdiff_t>(dst.first[static_cast<std::size_t>(v)]);
      std::copy(nList[static_cast<std::size_t>(v)].begin() + 1,
                nList[static_cast<std::size_t>(v)].end(), row);
      if (sortRows && n > 1) {
        std::sort(row, row + static_cast<std::ptrdiff_t>(n));
      }
    }
    slotSorted[slot] = sortRows;
  }

  void sortSlot(int slot) {
    if (slotSorted[slot]) {
      return;
    }
    const int nVertices = static_cast<int>(graph[slot].size());
    auto &g = graph[slot];
    for (int v = 0; v < nVertices; v++) {
      auto row = g[static_cast<std::size_t>(v)];
      if (row.end() - row.begin() > 1) {
        std::sort(const_cast<int *>(row.begin()), const_cast<int *>(row.end()));
      }
    }
    slotSorted[slot] = true;
  }

  void resolveBonds() {
    if (!bondsPending) {
      return;
    }
    bondsPending = false;
    const int nVertices = static_cast<int>(graph[which].size());
    sortSlot(which);
    if (pendingOld == which ||
        graph[pendingOld].size() != graph[which].size()) {
      changedBonds = 0;
      changedBondFraction = 0.0;
      return;
    }
    sortSlot(pendingOld);
    int directed = 0;
    for (int v = 0; v < nVertices; v++) {
      const auto a = graph[pendingOld][static_cast<std::size_t>(v)];
      const auto b = graph[which][static_cast<std::size_t>(v)];
      const int *ab = a.begin();
      const int *ae = a.end();
      const int *bb = b.begin();
      const int *be = b.end();
      if (std::equal(ab, ae, bb, be)) {
        continue;
      }
      while (ab != ae || bb != be) {
        if (bb == be || (ab != ae && *ab < *bb)) {
          ++directed;
          ++ab;
        } else if (ab == ae || *bb < *ab) {
          ++directed;
          ++bb;
        } else {
          ++ab;
          ++bb;
        }
      }
    }
    changedBonds = directed / 2;
    changedBondFraction =
        nVertices > 0 ? static_cast<double>(changedBonds) / nVertices : 0.0;
  }

  // Whole-frame enumeration into persistent per-source pools, then a copy
  // into the flattened buffer that keeps each ring's allocation.
  template <typename Graph> void fastEnumerate(const Graph &adj) {
    const int nVertices = static_cast<int>(adj.size());
    if (hopOnly) {
      ballsRefreshed = 0;
      ballsStale = false;
    } else {
      balls.refillAll(adj, ballRadius);
      ballsStale = false;
      ballsRefreshed = nVertices;
    }
    DistOracle oracle;
    oracle.adj = &adj;
    oracle.balls = hopOnly ? nullptr : &balls;
    oracle.hop2 = hopOnly;
    if (static_cast<int>(bySource.size()) != nVertices) {
      bySource.assign(static_cast<std::size_t>(nVertices), {});
    }
    int nThreads = 1;
#ifdef SEAMS_HAS_OPENMP
    if (!omp_in_parallel()) {
      nThreads = omp_get_max_threads();
    }
#endif
    if (static_cast<int>(scratches.size()) < nThreads) {
      scratches.resize(static_cast<std::size_t>(nThreads));
    }
#ifdef SEAMS_HAS_OPENMP
    if (nThreads > 1) {
#pragma omp parallel
      {
        RingScratch &scratch = scratches[static_cast<std::size_t>(omp_get_thread_num())];
        // Guided chunks: a nucleus and a partial shell put much heavier
        // sources next to empty ones, and a static share leaves the last
        // busy chunk with nothing left to steal.
#pragma omp for schedule(guided)
        for (int src = 0; src < nVertices; src++) {
          bySource[static_cast<std::size_t>(src)].reset();
          enumerateFromSource(adj, oracle, src, maxDepth, maxLvl, scratch,
                              bySource[static_cast<std::size_t>(src)]);
        }
      }
    } else
#endif
    {
      if (scratches.empty()) {
        scratches.emplace_back();
      }
      for (int src = 0; src < nVertices; src++) {
        bySource[static_cast<std::size_t>(src)].reset();
        enumerateFromSource(adj, oracle, src, maxDepth, maxLvl, scratches[0],
                            bySource[static_cast<std::size_t>(src)]);
      }
    }

    std::vector<std::size_t> offset(static_cast<std::size_t>(nVertices) + 1, 0);
    for (int src = 0; src < nVertices; src++) {
      offset[static_cast<std::size_t>(src) + 1] =
          offset[static_cast<std::size_t>(src)] +
          bySource[static_cast<std::size_t>(src)].n;
    }
    const std::size_t total = offset[static_cast<std::size_t>(nVertices)];
    if (flattened.size() < total) {
      flattened.resize(total);
    }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (!omp_in_parallel())
#endif
    for (int src = 0; src < nVertices; src++) {
      RingPool &pool = bySource[static_cast<std::size_t>(src)];
      std::size_t at = offset[static_cast<std::size_t>(src)];
      for (std::size_t k = 0; k < pool.n; k++) {
        // Exchange buffers with the previous frame. The ints are not
        // copied, and the pool keeps a buffer for the next frame.
        flattened[at++].swap(pool.storage[k]);
      }
    }
    if (flattened.size() != total) {
      flattened.resize(total);
    }
    recomputed = nVertices;
  }

  template <typename Graph>
  std::vector<std::vector<int>>
  ringsMeeting(const Graph &adj, const std::vector<int> &T,
               const std::vector<char> &inT) {
    const int n = static_cast<int>(T.size());
    std::vector<std::vector<std::vector<int>>> fresh(T.size());
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (!omp_in_parallel())
#endif
    {
      RingScratch scratch;
      DistOracle oracle;
      oracle.adj = &adj;
      oracle.balls = hopOnly ? nullptr : &balls;
      oracle.hop2 = hopOnly;
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(guided)
#endif
      for (int k = 0; k < n; k++) {
        enumerateFromSource(adj, oracle, T[static_cast<std::size_t>(k)],
                            maxDepth, maxLvl, scratch,
                            fresh[static_cast<std::size_t>(k)], &inT);
      }
    }
    std::size_t total = 0;
    for (const auto &f : fresh) {
      total += f.size();
    }
    std::vector<std::vector<int>> out;
    out.reserve(total);
    for (auto &f : fresh) {
      for (auto &r : f) {
        out.push_back(std::move(r));
      }
    }
    std::sort(out.begin(), out.end());
    return out;
  }

  void splice(const std::vector<char> &inT,
              std::vector<std::vector<int>> added) {
    std::vector<std::vector<int>> merged;
    merged.reserve(flattened.size() + added.size());
    auto a = added.begin();
    for (auto &r : flattened) {
      bool meets = false;
      for (const int v : r) {
        if (inT[static_cast<std::size_t>(v)]) {
          meets = true;
          break;
        }
      }
      if (meets) {
        continue;
      }
      while (a != added.end() && *a < r) {
        merged.push_back(std::move(*a++));
      }
      merged.push_back(std::move(r));
    }
    while (a != added.end()) {
      merged.push_back(std::move(*a++));
    }
    flattened.swap(merged);
  }
};

primitive::RingUpdater::RingUpdater(int maxDepth) : impl_(new Impl) {
  impl_->maxDepth = maxDepth;
  impl_->maxLvl = maxDepth / 2;
  impl_->ballRadius = std::max(impl_->maxLvl - 1, 1);
  impl_->reach = std::max((impl_->maxLvl - 2) / 2, 0);
  impl_->hopOnly = impl_->maxLvl <= 3;
}

primitive::RingUpdater::~RingUpdater() = default;
primitive::RingUpdater::RingUpdater(RingUpdater &&) noexcept = default;
primitive::RingUpdater &
primitive::RingUpdater::operator=(RingUpdater &&) noexcept = default;

int primitive::RingUpdater::lastRecomputedSources() const {
  return impl_->recomputed;
}

int primitive::RingUpdater::lastBallsRefreshed() const {
  return impl_->ballsRefreshed;
}

int primitive::RingUpdater::lastChangedBonds() const {
  impl_->resolveBonds();
  return impl_->changedBonds;
}

double primitive::RingUpdater::lastChangedBondFraction() const {
  impl_->resolveBonds();
  return impl_->changedBondFraction;
}

void primitive::RingUpdater::setFullPassFraction(double fraction) {
  impl_->fullPassFraction = fraction;
  impl_->useVertexFraction = true;
}

void primitive::RingUpdater::setBondCrossover(double changedBondsPerMolecule) {
  impl_->bondCrossover = changedBondsPerMolecule;
}

const std::vector<std::vector<int>> &
primitive::RingUpdater::update(const std::vector<std::vector<int>> &nList) {
  Impl &st = *impl_;
  const int nVertices = static_cast<int>(nList.size());
  st.resolveBonds();
  if (nVertices == 0 || st.maxDepth < 3) {
    st.hot = false;
    st.flattened.clear();
    st.recomputed = 0;
    st.ballsRefreshed = 0;
    st.changedBonds = 0;
    st.changedBondFraction = 0.0;
    st.primed = false;
    return st.flattened;
  }

  const bool sameSize =
      st.primed && st.graph[st.which].size() == static_cast<std::size_t>(nVertices);
  const int slot = sameSize ? 1 - st.which : 0;
  st.hot = false;
  st.writeGraph(nList, slot, true);
  if (!sameSize) {
    // A missing frame or a different molecule count leaves no row to diff.
    // Comparing past the stored prefix walks off the packed table.
    st.hot = false;
    st.which = slot;
    st.primed = true;
    st.changedBonds = 0;
    st.changedBondFraction = 0.0;
    st.fastEnumerate(st.graph[st.which]);
    return st.flattened;
  }

  std::vector<int> changed;
  int directed = 0;
  const auto &prev = st.graph[st.which];
  const auto &next = st.graph[slot];
  auto scanRow = [&](int v, int &localDirected, std::vector<int> &local) {
    const auto a = prev[static_cast<std::size_t>(v)];
    const auto b = next[static_cast<std::size_t>(v)];
    const int *ab = a.begin();
    const int *ae = a.end();
    const int *bb = b.begin();
    const int *be = b.end();
    if (std::equal(ab, ae, bb, be)) {
      return;
    }
    local.push_back(v);
    while (ab != ae || bb != be) {
      if (bb == be || (ab != ae && *ab < *bb)) {
        ++localDirected;
        ++ab;
      } else if (ab == ae || *bb < *ab) {
        ++localDirected;
        ++bb;
      } else {
        ++ab;
        ++bb;
      }
    }
  };
  int nThreads = 1;
#ifdef SEAMS_HAS_OPENMP
  if (!omp_in_parallel()) {
    nThreads = omp_get_max_threads();
  }
  if (nVertices >= 1024 && nThreads > 1) {
    struct alignas(64) Acc {
      std::vector<int> vs;
      int directed = 0;
    };
    std::vector<Acc> acc(static_cast<std::size_t>(nThreads));
#pragma omp parallel
    {
      Acc &me = acc[static_cast<std::size_t>(omp_get_thread_num())];
#pragma omp for schedule(static)
      for (int v = 0; v < nVertices; v++) {
        scanRow(v, me.directed, me.vs);
      }
    }
    for (const Acc &me : acc) {
      directed += me.directed;
      changed.insert(changed.end(), me.vs.begin(), me.vs.end());
    }
  } else
#endif
  {
    for (int v = 0; v < nVertices; v++) {
      scanRow(v, directed, changed);
    }
  }
  // Each undirected bond is an endpoint in two rows.
  st.changedBonds = directed / 2;
  st.changedBondFraction =
      nVertices > 0 ? static_cast<double>(st.changedBonds) / nVertices : 0.0;

  if (changed.empty()) {
    st.recomputed = 0;
    st.ballsRefreshed = 0;
    return st.flattened;
  }

  const auto atomLimit =
      static_cast<std::size_t>(st.fullPassFraction * nVertices);
  const bool overBonds =
      st.bondCrossover >= 0.0 && st.changedBondFraction > st.bondCrossover;
  const bool overAtoms =
      st.useVertexFraction && st.fullPassFraction < 1.0 &&
      (st.fullPassFraction <= 0.0 || changed.size() > atomLimit);
  if (overAtoms || overBonds) {
    st.which = slot;
    st.fastEnumerate(st.graph[st.which]);
    return st.flattened;
  }

  // T: within reach of a changed vertex in the union of both frames' graphs.
  std::vector<char> inT(static_cast<std::size_t>(nVertices), 0);
  std::vector<int> T = changed;
  for (const int v : changed) {
    inT[static_cast<std::size_t>(v)] = 1;
  }
  std::vector<int> frontier = changed;
  std::vector<int> wave;
  for (int depth = 1; depth <= st.reach && !frontier.empty(); depth++) {
    wave.clear();
    for (const int u : frontier) {
      for (const int w : prev[static_cast<std::size_t>(u)]) {
        if (w >= 0 && w < nVertices && !inT[static_cast<std::size_t>(w)]) {
          inT[static_cast<std::size_t>(w)] = 1;
          T.push_back(w);
          wave.push_back(w);
        }
      }
      for (const int w : next[static_cast<std::size_t>(u)]) {
        if (w >= 0 && w < nVertices && !inT[static_cast<std::size_t>(w)]) {
          inT[static_cast<std::size_t>(w)] = 1;
          T.push_back(w);
          wave.push_back(w);
        }
      }
    }
    frontier.swap(wave);
  }

  // Each source in T is enumerated once. A full pass enumerates n sources.
  // Those counts meet when T is the whole vertex set, and only then.
  if (static_cast<int>(T.size()) == nVertices) {
    st.which = slot;
    st.fastEnumerate(st.graph[st.which]);
    return st.flattened;
  }

  if (st.hopOnly) {
    st.ballsRefreshed = 0;
  } else if (st.ballsStale) {
    st.balls.refillAll(next, st.ballRadius);
    st.ballsStale = false;
    st.ballsRefreshed = nVertices;
  } else {
    st.ballsRefreshed =
        st.balls.refreshNear(prev, next, changed, st.ballRadius);
  }
  st.splice(inT, st.ringsMeeting(next, T, inT));
  st.which = slot;
  st.recomputed = static_cast<int>(T.size());
  return st.flattened;
}
