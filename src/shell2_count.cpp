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
#include <atomic>
#include <bit>
#include <cstdint>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

// Seven-ring counts for a shell two words wide. The walk lives in this file
// so the three-word shell, which does not call it, keeps the code it was
// tuned with.

namespace {

struct RowView {
  const std::size_t *first = nullptr;
  const int *cols = nullptr;
  int n = 0;

  struct Row {
    const int *b;
    const int *e;
    const int *begin() const { return b; }
    const int *end() const { return e; }
  };

  Row operator[](std::size_t u) const {
    return {cols + first[u], cols + first[u + 1]};
  }
  std::size_t size() const { return static_cast<std::size_t>(n); }
};

// Hop-2 answers every shortcut a ring of size at most 7 asks. This file is
// only called for that bound.
struct DistOracle {
  const RowView *adj = nullptr;
  bool hop2 = false;

  bool shorter(int u, int v, int ringDist) const {
    if (ringDist <= 1) {
      return false;
    }
    if (hop2 && ringDist <= 3) {
      const auto row = (*adj)[static_cast<std::size_t>(u)];
      for (const int w : row) {
        if (w == v) {
          return true;
        }
      }
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
    return false;
  }
};

bool ringIsPrimitive(const DistOracle &oracle, const std::vector<int> &ring, int skipA,
                     int skipB) {
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

inline bool pushUniqueMember(std::vector<int> &ring, int v) {
  for (const int x : ring) {
    if (x == v) {
      return false;
    }
  }
  ring.push_back(v);
  return true;
}


inline void tallyRing(const int *ring, int k, std::vector<int> *hist, int *perAtom,
                      int perStride) {
  if (perAtom != nullptr) {
    for (int i = 0; i < k; i++) {
      const int a = ring[i];
      if (a >= 0) {
        perAtom[static_cast<std::size_t>(a) * static_cast<std::size_t>(perStride) +
                static_cast<std::size_t>(k)] += 1;
      }
    }
  }
  if (hist != nullptr && k >= 0 && k < static_cast<int>(hist->size())) {
    (*hist)[static_cast<std::size_t>(k)]++;
  }
}

struct alignas(64) RingScratch {
  std::vector<int> lvl;
  std::vector<int> touched, frontier, next;
  std::vector<int> pathFirst, pathCount, pathBuf;
  std::vector<int> ring;
};

struct ChordBits {
  int n = 0;
  int stride = 0;
  std::vector<std::uint64_t> bits;
  std::vector<std::uint64_t> near;
};

std::atomic<int> gChordEpoch{0};

struct NoPool {
  std::size_t n = 0;
  std::size_t size() const { return n; }
  void push_back(const std::vector<int> &) { ++n; }
};

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


template <typename Graph, bool CheckCap>
bool shortestPathsFrom(const Graph &adjacency, int src, RingScratch &scr,
                       int cap) {
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
    if constexpr (CheckCap) {
      for (const int w : adjacency[v]) {
        if (w < 0 || w >= nVertices || scr.lvl[w] != lv - 1) {
          continue;
        }
        count += scr.pathCount[w];
      }
      if (count > cap) {
        return false;
      }
    }
    const int plen = lv - 1;
    for (const int w : adjacency[v]) {
      if (w < 0 || w >= nVertices || scr.lvl[w] != lv - 1) {
        continue;
      }
      if constexpr (!CheckCap) {
        count += scr.pathCount[w];
      }
      for (int k = 0; k < scr.pathCount[w]; k++) {
        const int from = scr.pathFirst[w] + k * plen;
        const int dest = static_cast<int>(scr.pathBuf.size());
        scr.pathBuf.resize(scr.pathBuf.size() + static_cast<std::size_t>(lv));
        if (plen > 0) {
          std::copy_n(scr.pathBuf.begin() + from, plen,
                      scr.pathBuf.begin() + dest);
        }
        scr.pathBuf[static_cast<std::size_t>(dest + plen)] = v;
      }
    }
    scr.pathFirst[v] = first;
    scr.pathCount[v] = count;
  }
  return true;
}


template <typename Graph>
void fillChordBits(const Graph &adjacency, ChordBits &out) {
  const int nVertices = static_cast<int>(adjacency.size());
  const int stride = (nVertices + 63) >> 6;
  out.n = nVertices;
  out.stride = stride;
  out.bits.assign(static_cast<std::size_t>(nVertices) * static_cast<std::size_t>(stride), 0);
  for (int u = 0; u < nVertices; u++) {
    std::uint64_t *row =
        out.bits.data() + static_cast<std::size_t>(u) * static_cast<std::size_t>(stride);
    for (const int v : adjacency[static_cast<std::size_t>(u)]) {
      if (v >= 0 && v < nVertices) {
        row[static_cast<unsigned>(v) >> 6] |= 1ull << (static_cast<unsigned>(v) & 63u);
      }
    }
  }
  out.near.assign(out.bits.size(), 0);
  for (int u = 0; u < nVertices; u++) {
    std::uint64_t *dst =
        out.near.data() + static_cast<std::size_t>(u) * static_cast<std::size_t>(stride);
    const std::uint64_t *self =
        out.bits.data() + static_cast<std::size_t>(u) * static_cast<std::size_t>(stride);
    for (int w = 0; w < stride; w++) {
      dst[w] = self[w];
    }
    for (const int v : adjacency[static_cast<std::size_t>(u)]) {
      if (v < 0 || v >= nVertices) {
        continue;
      }
      const std::uint64_t *row =
          out.bits.data() + static_cast<std::size_t>(v) * static_cast<std::size_t>(stride);
      for (int w = 0; w < stride; w++) {
        dst[w] |= row[w];
      }
    }
  }
}


// A four-ring census. The adjacent pairs among one vertex's arms are a
// popcount of the shell rows, so the pair walk is not repeated.
int countFourShell(const std::vector<int> &L2, const std::uint64_t *shellArms,
                   const std::uint64_t *shellAdj, int lwords) {
  int c4 = 0;
  for (const int b : L2) {
    const std::uint64_t *arms =
        shellArms + static_cast<std::size_t>(b) * static_cast<std::size_t>(lwords);
    int n = 0;
    for (int w = 0; w < lwords; w++) {
      n += std::popcount(arms[w]);
    }
    if (n < 2) {
      continue;
    }
    int cuts = 0;
    for (int w = 0; w < lwords; w++) {
      std::uint64_t word = arms[w];
      const int base = w << 6;
      while (word != 0) {
        const int ra = base + std::countr_zero(word);
        word &= word - 1;
        const std::uint64_t *ba =
            shellAdj + static_cast<std::size_t>(ra) * static_cast<std::size_t>(lwords);
        int hit = 0;
        for (int k = 0; k < lwords; k++) {
          hit += std::popcount(ba[k] & arms[k]);
        }
        cuts += hit;
      }
    }
    c4 += n * (n - 1) / 2 - cuts / 2;
  }
  return c4;
}

template <typename Graph, typename Out, bool PerAtom, int ShellWords = 0>
__attribute__((hot))
void enumerateChordlessPlain(const Graph &adjacency, const DistOracle &oracle,
                        int src, int maxDepth, int minSize, bool legacy,
                        const std::vector<int> &lvl, Out &out,
                        std::vector<int> *hist, int *perAtom, int perStride) {
  const int nVertices = static_cast<int>(adjacency.size());
  if (maxDepth < 3 || nVertices <= 0) {
    return;
  }
  thread_local ChordBits store;
  thread_local int storeEpoch = 0;
  const int epoch = gChordEpoch.load(std::memory_order_relaxed);
  if (epoch == 0 || storeEpoch != epoch || store.n != nVertices) {
    fillChordBits(adjacency, store);
    storeEpoch = epoch;
  }
  const int stride = store.stride;
  const std::uint64_t *bitData = store.bits.data();
  const std::uint64_t *nearData = store.near.data();
  auto adjacent = [&](int u, int v) {
    const std::uint64_t *row =
        bitData + static_cast<std::size_t>(u) * static_cast<std::size_t>(stride);
    return ((row[static_cast<unsigned>(v) >> 6] >> (static_cast<unsigned>(v) & 63u)) & 1ull) != 0ull;
  };
  // Distance at most 2: adjacent, or a shared neighbour. One row per vertex,
  // so each later test is a bit read. The row includes the vertex itself
  // when it has a neighbour; no query asks within2 of a vertex with itself.
  auto within2 = [&](int u, int v) {
    const std::uint64_t *row =
        nearData + static_cast<std::size_t>(u) * static_cast<std::size_t>(stride);
    return ((row[static_cast<unsigned>(v) >> 6] >> (static_cast<unsigned>(v) & 63u)) &
            1ull) != 0ull;
  };
  // The hop-2 oracle answers a shortcut of length 1 or 2. Other oracles
  // (and the legacy pair set) go through the shared predicate.
  const bool fast = !legacy && oracle.hop2;
  auto record = [&](const int *r, int k) {
    if constexpr (PerAtom) {
      tallyRing(r, k, hist, perAtom, perStride);
      return;
    }
    if (hist != nullptr) {
      if (k < static_cast<int>(hist->size())) {
        (*hist)[static_cast<std::size_t>(k)]++;
      }
      return;
    }
    out.push_back(std::vector<int>(r, r + k));
  };
  auto keep = [&](const int *r, int k, int skipB) {
    if (k < minSize || k > maxDepth) {
      return;
    }
    if (!fast) {
      const std::vector<int> ring(r, r + k);
      if (ringIsPrimitive(oracle, ring, legacy ? -1 : 0, legacy ? -1 : skipB)) {
        record(r, k);
      }
      return;
    }
    for (int i = 0; i < k; i++) {
      if (i == 0 || i == skipB) {
        continue;
      }
      for (int j = i + 2; j < k; j++) {
        if (j == skipB) {
          continue;
        }
        const int sep = j - i;
        const int ringDist = std::min(sep, k - sep);
        if (ringDist <= 1) {
          continue;
        }
        const bool cut = ringDist <= 2 ? adjacent(r[i], r[j]) : within2(r[i], r[j]);
        if (cut) {
          return;
        }
      }
    }
    record(r, k);
  };
  std::vector<int> L1, L2, L3;
  L1.reserve(64);
  L2.reserve(64);
  L3.reserve(64);
  for (int v = 0; v < nVertices; v++) {
    const int d = lvl[static_cast<std::size_t>(v)];
    if (d == 1) {
      L1.push_back(v);
    } else if (d == 2) {
      L2.push_back(v);
    } else if (d == 3 && maxDepth >= 6) {
      L3.push_back(v);
    }
  }
  int ring[7];
  ring[0] = src;
  int c5 = 0;
  if (maxDepth >= 3) {
    for (std::size_t i = 0; i < L1.size(); i++) {
      const int a = L1[i];
      for (std::size_t j = i + 1; j < L1.size(); j++) {
        const int b = L1[j];
        if (!adjacent(a, b)) {
          continue;
        }
        ring[1] = a;
        ring[2] = b;
        if (fast) {
          if (minSize <= 3) {
            record(ring, 3);
          }
        } else {
          keep(ring, 3, -1);
        }
      }
    }
  }
  if (maxDepth >= 4 && !L2.empty()) {
    // The shell rows answer both censuses. count4 stays on only when those
    // rows will be built, so skipping the pair walk cannot drop a four-ring.
    const bool countShell =
        !PerAtom && hist != nullptr && fast && minSize <= 5 && !L1.empty();
    const bool count4 = countShell && minSize <= 4;
    const bool count5 = countShell;
    std::vector<std::vector<int>> arms2(static_cast<std::size_t>(nVertices));
    for (const int b : L2) {
      auto &arms = arms2[static_cast<std::size_t>(b)];
      for (const int a : L1) {
        if (adjacent(a, b)) {
          arms.push_back(a);
        }
      }
      if (count4) {
        continue;
      }
      for (std::size_t i = 0; i < arms.size(); i++) {
        for (std::size_t j = i + 1; j < arms.size(); j++) {
          const int a = arms[i];
          const int c = arms[j];
          if (fast) {
            if (adjacent(a, c) || minSize > 4) {
              continue;
            }
          }
          if (a > c) {
            ring[1] = c;
            ring[2] = b;
            ring[3] = a;
          } else {
            ring[1] = a;
            ring[2] = b;
            ring[3] = c;
          }
          if (fast) {
            record(ring, 4);
          } else {
            keep(ring, 4, 2);
          }
        }
      }
    }
    std::vector<int> rank1;
    std::vector<std::uint64_t> shellAdj;
    std::vector<std::uint64_t> shellArms;
    std::vector<std::uint64_t> rightBits;
    std::vector<std::uint64_t> leftBits;
    int lwords = 0;
    if (countShell) {
      // Neighbour bits of the first shell only. A full-frame row is several
      // times wider, and every arm sits in this shell.
      const int nL1 = static_cast<int>(L1.size());
      lwords = (nL1 + 63) >> 6;
      rank1.assign(static_cast<std::size_t>(nVertices), -1);
      for (int i = 0; i < nL1; i++) {
        rank1[static_cast<std::size_t>(L1[static_cast<std::size_t>(i)])] = i;
      }
      shellAdj.assign(static_cast<std::size_t>(nL1) * static_cast<std::size_t>(lwords), 0);
      shellArms.assign(static_cast<std::size_t>(nVertices) * static_cast<std::size_t>(lwords), 0);
      rightBits.assign(static_cast<std::size_t>(lwords), 0);
      leftBits.assign(static_cast<std::size_t>(lwords), 0);
      for (int i = 0; i < nL1; i++) {
        const int a = L1[static_cast<std::size_t>(i)];
        std::uint64_t *row =
            shellAdj.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(lwords);
        for (const int b : adjacency[static_cast<std::size_t>(a)]) {
          if (b < 0 || b >= nVertices) {
            continue;
          }
          const int rb = rank1[static_cast<std::size_t>(b)];
          if (rb >= 0) {
            row[static_cast<unsigned>(rb) >> 6] |= 1ull << (static_cast<unsigned>(rb) & 63u);
          }
        }
      }
      for (const int b : L2) {
        std::uint64_t *row =
            shellArms.data() + static_cast<std::size_t>(b) * static_cast<std::size_t>(lwords);
        for (const int a : arms2[static_cast<std::size_t>(b)]) {
          const int ra = rank1[static_cast<std::size_t>(a)];
          if (ra >= 0) {
            row[static_cast<unsigned>(ra) >> 6] |= 1ull << (static_cast<unsigned>(ra) & 63u);
          }
        }
      }
    }
    if (count4) {
      // A four-ring is two arms that are not adjacent. The arm bit row
      // already lists those arms.
      const int c4 = countFourShell(L2, shellArms.data(), shellAdj.data(), lwords);
      if (c4 != 0 && 4 < static_cast<int>(hist->size())) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp atomic
#endif
        (*hist)[4] += c4;
      }
    }
    if (maxDepth >= 5 && count5) {
      // |{(a,c)}| with the three shortcut tests, without listing rings.
      // a is an arm of p that is not an arm of q, so the bit row already
      // names those arms. Walking the bits skips the arms both vertices share.
      for (std::size_t pi = 0; pi < L2.size(); pi++) {
        const int p = L2[pi];
        if (arms2[static_cast<std::size_t>(p)].empty()) {
          continue;
        }
        const std::uint64_t *pArms =
            shellArms.data() + static_cast<std::size_t>(p) * static_cast<std::size_t>(lwords);
        for (std::size_t qi = pi + 1; qi < L2.size(); qi++) {
          const int q = L2[qi];
          if (!adjacent(p, q)) {
            continue;
          }
          const std::uint64_t *qArms =
              shellArms.data() + static_cast<std::size_t>(q) * static_cast<std::size_t>(lwords);
          int pc = 0;
          int lc = 0;
          for (int w = 0; w < lwords; w++) {
            rightBits[static_cast<std::size_t>(w)] = qArms[w] & ~pArms[w];
            leftBits[static_cast<std::size_t>(w)] = pArms[w] & ~qArms[w];
            pc += std::popcount(rightBits[static_cast<std::size_t>(w)]);
            lc += std::popcount(leftBits[static_cast<std::size_t>(w)]);
          }
          if (pc == 0 || lc == 0) {
            continue;
          }
          // Each exclusive arm of one side is tested against the other side.
          // The two sides are the same count, so the shorter one is enough.
          const std::uint64_t *walk = lc <= pc ? leftBits.data() : rightBits.data();
          const std::uint64_t *other = lc <= pc ? rightBits.data() : leftBits.data();
          const int otherN = lc <= pc ? pc : lc;
          // ShellWords is chosen once per frame from the widest shell. A
          // narrower frame never instantiates this arm, so its loop is the
          // general one.
          if constexpr (ShellWords == 3) {
            if (lwords == 3) {
              const std::uint64_t o0 = other[0];
              const std::uint64_t o1 = other[1];
              const std::uint64_t o2 = other[2];
              for (int w = 0; w < 3; w++) {
                std::uint64_t word = walk[w];
                const int base = w << 6;
                while (word != 0) {
                  const int ra = base + std::countr_zero(word);
                  word &= word - 1;
                  const std::uint64_t *ba =
                      shellAdj.data() + static_cast<std::size_t>(ra) * 3;
                  const int cut = std::popcount(ba[0] & o0) + std::popcount(ba[1] & o1) +
                                  std::popcount(ba[2] & o2);
                  c5 += otherN - cut;
                }
              }
            } else {
              for (int w = 0; w < lwords; w++) {
                std::uint64_t word = walk[w];
                const int base = w << 6;
                while (word != 0) {
                  const int ra = base + std::countr_zero(word);
                  word &= word - 1;
                  const std::uint64_t *ba = shellAdj.data() +
                                            static_cast<std::size_t>(ra) *
                                                static_cast<std::size_t>(lwords);
                  int cut = 0;
                  for (int k = 0; k < lwords; k++) {
                    cut += std::popcount(ba[k] & other[k]);
                  }
                  c5 += otherN - cut;
                }
              }
            }
          } else {
            for (int w = 0; w < lwords; w++) {
              std::uint64_t word = walk[w];
              const int base = w << 6;
              while (word != 0) {
                const int ra = base + std::countr_zero(word);
                word &= word - 1;
                const std::uint64_t *ba =
                    shellAdj.data() + static_cast<std::size_t>(ra) * static_cast<std::size_t>(lwords);
                int cut = 0;
                for (int k = 0; k < lwords; k++) {
                  cut += std::popcount(ba[k] & other[k]);
                }
                c5 += otherN - cut;
              }
            }
          }
        }
      }
    } else if (maxDepth >= 5) {
      for (std::size_t pi = 0; pi < L2.size(); pi++) {
        const int p = L2[pi];
        const auto &ap = arms2[static_cast<std::size_t>(p)];
        if (ap.empty()) {
          continue;
        }
        for (std::size_t qi = pi + 1; qi < L2.size(); qi++) {
          const int q = L2[qi];
          if (!adjacent(p, q)) {
            continue;
          }
          const auto &aq = arms2[static_cast<std::size_t>(q)];
          for (const int a : ap) {
            for (const int c : aq) {
              if (a == c) {
                continue;
              }
              if (fast && (minSize > 5 || adjacent(a, c) || adjacent(a, q) ||
                           adjacent(p, c))) {
                continue;
              }
              if (a < c) {
                ring[1] = a;
                ring[2] = p;
                ring[3] = q;
                ring[4] = c;
              } else {
                ring[1] = c;
                ring[2] = q;
                ring[3] = p;
                ring[4] = a;
              }
              if (fast) {
                record(ring, 5);
              } else {
                keep(ring, 5, -1);
              }
            }
          }
        }
      }
    }
    if (hist != nullptr && c5 != 0 && 5 < static_cast<int>(hist->size())) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp atomic
#endif
      (*hist)[5] += c5;
    }
  }
  if (maxDepth >= 6 && !L3.empty()) {
    std::vector<std::vector<std::pair<int, int>>> arms3(static_cast<std::size_t>(nVertices));
    for (const int p : L3) {
      auto &arms = arms3[static_cast<std::size_t>(p)];
      for (const int b : L2) {
        if (!adjacent(b, p)) {
          continue;
        }
        for (const int a : L1) {
          if (adjacent(a, b)) {
            arms.emplace_back(a, b);
          }
        }
      }
    }
    if (!PerAtom && hist != nullptr && fast) {
      // A six-ring fails for every arm pair that leaves the same second-shell
      // vertex, or whose two second-shell vertices are adjacent. Group the
      // arms by that vertex and drop those pairs before walking them.
      thread_local std::vector<std::vector<int>> byB;
      thread_local std::vector<int> usedB;
      if (static_cast<int>(byB.size()) < nVertices) {
        byB.resize(static_cast<std::size_t>(nVertices));
      }
      int c6 = 0;
      for (const int p : L3) {
        const auto &ap = arms3[static_cast<std::size_t>(p)];
        usedB.clear();
        for (const auto &arm : ap) {
          auto &bucket = byB[static_cast<std::size_t>(arm.second)];
          if (bucket.empty()) {
            usedB.push_back(arm.second);
          }
          bucket.push_back(arm.first);
        }
        for (std::size_t i = 0; i < usedB.size(); i++) {
          const int b = usedB[i];
          const auto &as = byB[static_cast<std::size_t>(b)];
          for (std::size_t j = i + 1; j < usedB.size(); j++) {
            const int b2 = usedB[j];
            if (adjacent(b, b2)) {
              continue;
            }
            const auto &as2 = byB[static_cast<std::size_t>(b2)];
            for (const int a : as) {
              if (within2(a, b2)) {
                continue;
              }
              for (const int a2 : as2) {
                if (a == a2 || adjacent(a, a2) || within2(b, a2)) {
                  continue;
                }
                c6++;
              }
            }
          }
        }
        for (const int b : usedB) {
          byB[static_cast<std::size_t>(b)].clear();
        }
      }
      if (c6 != 0 && 6 < static_cast<int>(hist->size())) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp atomic
#endif
        (*hist)[6] += c6;
      }
    } else {
      for (const int p : L3) {
        const auto &ap = arms3[static_cast<std::size_t>(p)];
        for (std::size_t i = 0; i < ap.size(); i++) {
          for (std::size_t j = i + 1; j < ap.size(); j++) {
            int a = ap[i].first;
            int b = ap[i].second;
            int a2 = ap[j].first;
            int b2 = ap[j].second;
            if (a == a2 || b == b2) {
              continue;
            }
            if (a > a2) {
              std::swap(a, a2);
              std::swap(b, b2);
            }
            ring[1] = a;
            ring[2] = b;
            ring[3] = p;
            ring[4] = b2;
            ring[5] = a2;
            keep(ring, 6, 3);
          }
        }
      }
    }
    if (maxDepth >= 7) {
      // The arms of one level-3 vertex are stored grouped by their near
      // vertex, and the left filter keeps that order. A run shares every
      // test that does not look at the far vertex. The three-word shell
      // keeps the per-arm walk; that instantiation is the dense frame.
      std::vector<std::pair<int, int>> lefts;
      std::vector<std::pair<int, int>> rights;
      std::vector<int> dSlot(static_cast<std::size_t>(nVertices), -1);
      std::vector<int> touchedD;
      std::vector<std::vector<int>> byD;
      int c7 = 0;
      for (std::size_t pi = 0; pi < L3.size(); pi++) {
        const int p = L3[pi];
        const auto &ap = arms3[static_cast<std::size_t>(p)];
        if (ap.empty()) {
          continue;
        }
        for (std::size_t qi = pi + 1; qi < L3.size(); qi++) {
          const int q = L3[qi];
          if (!adjacent(p, q)) {
            continue;
          }
          const auto &aq = arms3[static_cast<std::size_t>(q)];
          // Tests that do not depend on the other arm. Applying them first
          // drops every pair they would reject, and the rings that remain
          // are emitted in the same order.
          if (fast && minSize > 7) {
            continue;
          }
          lefts.clear();
          rights.clear();
          for (const auto &arm : ap) {
            if (fast && (adjacent(arm.first, p) || adjacent(arm.second, q) ||
                         within2(arm.first, q))) {
              continue;
            }
            lefts.push_back(arm);
          }
          for (const auto &arm : aq) {
            if (fast && (adjacent(q, arm.first) || within2(p, arm.first))) {
              continue;
            }
            rights.push_back(arm);
          }
          if (!PerAtom && hist != nullptr && fast) {
            // Arms that share the far vertex fail or pass together on every
            // test that does not look at the near vertex. The hits are a
            // plain total: one add per source, not one atomic per ring.
            if constexpr (ShellWords != 3) {
              // Test each neighbour of the far vertex in place. Staging the
              // ones that clear the near vertex into a side vector is slower
              // than repeating that test.
              thread_local std::vector<int> slotOf;
              thread_local std::vector<int> usedD;
              thread_local std::vector<std::vector<int>> bucket;
              if (static_cast<int>(slotOf.size()) != nVertices) {
                slotOf.assign(static_cast<std::size_t>(nVertices), -1);
              }
              usedD.clear();
              for (const auto &right : rights) {
                const int d = right.second;
                int slot = slotOf[static_cast<std::size_t>(d)];
                if (slot < 0) {
                  slot = static_cast<int>(usedD.size());
                  slotOf[static_cast<std::size_t>(d)] = slot;
                  if (slot == static_cast<int>(bucket.size())) {
                    bucket.emplace_back();
                  }
                  bucket[static_cast<std::size_t>(slot)].clear();
                  usedD.push_back(d);
                }
                bucket[static_cast<std::size_t>(slot)].push_back(right.first);
              }
              const bool record7 = 7 < static_cast<int>(hist->size());
              for (std::size_t i0 = 0; i0 < lefts.size();) {
                const int b = lefts[i0].second;
                std::size_t i1 = i0 + 1;
                while (i1 < lefts.size() && lefts[i1].second == b) {
                  ++i1;
                }
                for (int slot = 0; slot < static_cast<int>(usedD.size()); slot++) {
                  const int d = usedD[static_cast<std::size_t>(slot)];
                  if (d == b || adjacent(p, d) || within2(b, d)) {
                    continue;
                  }
                  const auto &cs = bucket[static_cast<std::size_t>(slot)];
                  for (std::size_t k = i0; k < i1; k++) {
                    const int a = lefts[k].first;
                    if (within2(a, d)) {
                      continue;
                    }
                    for (const int c : cs) {
                      if (c == a || adjacent(a, c) || within2(b, c)) {
                        continue;
                      }
                      if (record7) {
                        c7++;
                      }
                    }
                  }
                }
                i0 = i1;
              }
              for (const int d : usedD) {
                slotOf[static_cast<std::size_t>(d)] = -1;
              }
            } else {
              byD.clear();
              touchedD.clear();
              for (const auto &right : rights) {
                const int d = right.second;
                int slot = dSlot[static_cast<std::size_t>(d)];
                if (slot < 0) {
                  slot = static_cast<int>(byD.size());
                  dSlot[static_cast<std::size_t>(d)] = slot;
                  byD.emplace_back();
                  touchedD.push_back(d);
                }
                byD[static_cast<std::size_t>(slot)].push_back(right.first);
              }
              for (const auto &left : lefts) {
                const int a = left.first;
                const int b = left.second;
                for (int slot = 0; slot < static_cast<int>(byD.size()); slot++) {
                  const int d = touchedD[static_cast<std::size_t>(slot)];
                  if (d == b || adjacent(p, d) || within2(a, d) || within2(b, d)) {
                    continue;
                  }
                  for (const int c : byD[static_cast<std::size_t>(slot)]) {
                    if (c == a || adjacent(a, c) || within2(b, c)) {
                      continue;
                    }
                    if (7 < static_cast<int>(hist->size())) {
                      c7++;
                    }
                  }
                }
              }
              for (const int d : touchedD) {
                dSlot[static_cast<std::size_t>(d)] = -1;
              }
            }
          } else {
            for (const auto &left : lefts) {
              for (const auto &right : rights) {
                const int a = left.first;
                const int b = left.second;
                const int c = right.first;
                const int d = right.second;
                if (a == c || b == d) {
                  continue;
                }
                if (fast && (adjacent(a, c) || adjacent(p, d) || within2(a, d) ||
                             within2(b, d) || within2(b, c))) {
                  continue;
                }
                if (a < c) {
                  ring[1] = a;
                  ring[2] = b;
                  ring[3] = p;
                  ring[4] = q;
                  ring[5] = d;
                  ring[6] = c;
                } else {
                  ring[1] = c;
                  ring[2] = d;
                  ring[3] = q;
                  ring[4] = p;
                  ring[5] = b;
                  ring[6] = a;
                }
                if (fast) {
                  record(ring, 7);
                } else {
                  keep(ring, 7, -1);
                }
              }
            }
          }
        }
      }
      if (hist != nullptr && c7 != 0 && 7 < static_cast<int>(hist->size())) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp atomic
#endif
        (*hist)[7] += c7;
      }
    }
  }
}

void countOne(const RowView &graph, const DistOracle &oracle, RingScratch &scr, int src,
              int nVertices, int maxDepth, int maxLvl, int minSize, int cap,
              std::vector<int> *hist) {
  if (static_cast<int>(scr.lvl.size()) != nVertices) {
    scr.lvl.assign(static_cast<std::size_t>(nVertices), -1);
    scr.touched.clear();
  }
  levelsFrom(graph, src, maxLvl, 0, scr.lvl, scr.touched, scr.frontier, scr.next);
  for (const int v : scr.touched) {
    if (v < src) {
      scr.lvl[v] = -1;
    }
  }
  NoPool out;
  const bool stored = cap > 0 ? shortestPathsFrom<RowView, true>(graph, src, scr, cap)
                              : shortestPathsFrom<RowView, false>(graph, src, scr, 0);
  if (!stored) {
    enumerateChordlessPlain<RowView, NoPool, false, 0>(
        graph, oracle, src, maxDepth, minSize, false, scr.lvl, out, hist, nullptr, 0);
    return;
  }
  auto keep = [&](std::vector<int> &r) {
    if (hist != nullptr && static_cast<int>(r.size()) < static_cast<int>(hist->size())) {
      (*hist)[r.size()]++;
    }
  };
  std::vector<int> &ring = scr.ring;
  for (const int p : scr.touched) {
    const int lp = scr.lvl[p];
    if (lp < 1 || lp > maxLvl) {
      continue;
    }
    const int *pathsP = scr.pathBuf.data() + scr.pathFirst[p];
    const int nP = scr.pathCount[p];
    if (2 * lp >= minSize && 2 * lp <= maxDepth) {
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
          if (ok && ringIsPrimitive(oracle, ring, 0, lp)) {
            keep(ring);
          }
        }
      }
    }
    if (2 * lp + 1 >= minSize && 2 * lp + 1 <= maxDepth) {
      for (const int q : graph[static_cast<std::size_t>(p)]) {
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
            if (ok && ringIsPrimitive(oracle, ring, 0, -1)) {
              keep(ring);
            }
          }
        }
      }
    }
  }
}

} // namespace

void countSourceShell2(const std::size_t *rowFirst, const int *cols, int nVertices, int src,
                       int maxDepth, int maxLvl, int minSize, int cap, std::vector<int> *hist,
                       int epoch) {
  RowView graph{rowFirst, cols, nVertices};
  DistOracle oracle;
  oracle.adj = &graph;
  oracle.hop2 = true;
  gChordEpoch.store(epoch, std::memory_order_relaxed);
  thread_local RingScratch scr;
  countOne(graph, oracle, scr, src, nVertices, maxDepth, maxLvl, minSize, cap, hist);
}
