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

#include <generic.hpp>
#include <par_team.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <span>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

const std::vector<int> kEmptyRings;

/**
 * @details Visiting order for the greedy cage assembly, derived from each
 *  ring's sorted atom content rather than its position in the input. The
 *  assembly claims rings as it accepts cages, so the accepted set depends on
 *  the order rings are tested; keying that order to the rings themselves
 *  makes the result independent of how the ring network was enumerated.
 *  Each atom enters the key as its rank when @a atomRank is given, and as
 *  its index otherwise.
 */
// Lexicographically least reading of a ring under rotation and reversal.
// A least reading starts on a minimum element, so only those starts are
// built. Rings in the cage search are at most a handful of atoms, so the
// key lives on the stack.
struct CycleKey {
  int n = 0;
  std::array<int, 16> v{};

  bool operator<(const CycleKey &other) const {
    if (n != other.n) {
      return n < other.n;
    }
    for (int i = 0; i < n; i++) {
      if (v[static_cast<size_t>(i)] != other.v[static_cast<size_t>(i)]) {
        return v[static_cast<size_t>(i)] < other.v[static_cast<size_t>(i)];
      }
    }
    return false;
  }
};

CycleKey bestCycle(std::span<const int> ringAtoms) {
  CycleKey best;
  const int k = static_cast<int>(ringAtoms.size());
  best.n = k;
  if (k <= 0 || k > static_cast<int>(best.v.size())) {
    return best;
  }
  const int *atoms = ringAtoms.data();
  int minv = atoms[0];
  for (int i = 1; i < k; ++i) {
    if (atoms[i] < minv) {
      minv = atoms[i];
    }
  }
  CycleKey candidate;
  candidate.n = k;
  bool any = false;
  auto take = [&](int start, int step) {
    int idx = start;
    for (int m = 0; m < k; ++m) {
      candidate.v[static_cast<std::size_t>(m)] = atoms[idx];
      idx += step;
      if (idx >= k) {
        idx = 0;
      } else if (idx < 0) {
        idx = k - 1;
      }
    }
    if (!any || candidate < best) {
      best = candidate;
      any = true;
    }
  };
  for (int start = 0; start < k; ++start) {
    if (atoms[start] != minv) {
      continue;
    }
    take(start, 1);
    take(start, -1);
  }
  return best;
}

// Distinct values in xs[0, n). Callers only need the count sort+unique
// would return. A stamp over the value span is linear in n; a wide or
// negative span falls back to a sort of a copy.
int countDistinct(const int *xs, int n) {
  if (n <= 1) {
    return n > 0 ? n : 0;
  }
  int lo = xs[0];
  int hi = xs[0];
  for (int i = 1; i < n; ++i) {
    if (xs[i] < lo) {
      lo = xs[i];
    } else if (xs[i] > hi) {
      hi = xs[i];
    }
  }
  const int span = hi - lo + 1;
  if (lo < 0 || span > 8192) {
    std::vector<int> copy(xs, xs + n);
    std::sort(copy.begin(), copy.end());
    return static_cast<int>(std::unique(copy.begin(), copy.end()) -
                            copy.begin());
  }
  struct Stamp {
    std::vector<int> seen;
    int tick = 0;
  };
  thread_local Stamp stamp;
  if (static_cast<int>(stamp.seen.size()) < span) {
    stamp.seen.assign(static_cast<std::size_t>(span), 0);
  }
  if (++stamp.tick <= 0) {
    std::fill(stamp.seen.begin(), stamp.seen.end(), 0);
    stamp.tick = 1;
  }
  const int mark = stamp.tick;
  int uniq = 0;
  for (int i = 0; i < n; ++i) {
    int &slot = stamp.seen[static_cast<std::size_t>(xs[i] - lo)];
    if (slot != mark) {
      slot = mark;
      ++uniq;
    }
  }
  return uniq;
}

struct AtomKey {
  int n = 0;
  std::array<int, 16> id{};

  bool operator<(const AtomKey &other) const {
    if (n != other.n) {
      return n < other.n;
    }
    for (int i = 0; i < n; i++) {
      if (id[static_cast<size_t>(i)] != other.id[static_cast<size_t>(i)]) {
        return id[static_cast<size_t>(i)] < other.id[static_cast<size_t>(i)];
      }
    }
    return false;
  }
};

std::vector<int> canonicalRingOrder(const std::vector<std::vector<int>> &rings,
                                    const std::vector<int> &atomRank) {
  const int ranked = static_cast<int>(atomRank.size());
  const int n = static_cast<int>(rings.size());
  // Stack keys: one heap vector per ring was the serial cost of this order.
  std::vector<AtomKey> keys(static_cast<size_t>(n));
  std::vector<CycleKey> cycles(static_cast<size_t>(n));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n; i++) {
    const auto &ring = rings[static_cast<size_t>(i)];
    const int len = static_cast<int>(std::min(ring.size(), std::size_t{16}));
    std::array<int, 16> mapped{};
    for (int a = 0; a < len; a++) {
      int atom = ring[static_cast<size_t>(a)];
      if (atom >= 0 && atom < ranked) {
        atom = atomRank[static_cast<size_t>(atom)];
      }
      mapped[static_cast<size_t>(a)] = atom;
    }
    cycles[static_cast<size_t>(i)] =
        bestCycle(std::span<const int>(mapped.data(), static_cast<size_t>(len)));
    AtomKey key;
    key.n = len;
    std::copy_n(mapped.begin(), len, key.id.begin());
    std::sort(key.id.begin(), key.id.begin() + len);
    keys[static_cast<size_t>(i)] = key;
  }
  std::vector<int> order(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    order[static_cast<size_t>(i)] = i;
  }
  // Two rings over the same atoms differ in their cycle. A stable sort
  // keeps the earlier index when the keys and the cycles agree.
  seams::parallelStableSort(order, [&](int a, int b) {
    if (keys[static_cast<size_t>(a)] < keys[static_cast<size_t>(b)]) {
      return true;
    }
    if (keys[static_cast<size_t>(b)] < keys[static_cast<size_t>(a)]) {
      return false;
    }
    if (cycles[static_cast<size_t>(a)] < cycles[static_cast<size_t>(b)]) {
      return true;
    }
    if (cycles[static_cast<size_t>(b)] < cycles[static_cast<size_t>(a)]) {
      return false;
    }
    return a < b;
  });
  return order;
}

static bool basalFromMatrix(const bool bond[6][6]) {
  static constexpr int kPerm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                                      {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  for (int q = 0; q < 2; q++) {
    bool quiet = true;
    for (int t = 0; t < 3 && quiet; t++) {
      const int j = q + 1 + 2 * t;
      for (int i = 0; i < 6; i++) {
        if (bond[i][j % 6]) {
          quiet = false;
          break;
        }
      }
    }
    if (!quiet) {
      continue;
    }
    for (int p = 0; p < 2; p++) {
      for (const auto &s : kPerm) {
        if (bond[p][q + 2 * s[0]] && bond[p + 2][q + 2 * s[1]] &&
            bond[p + 4][q + 2 * s[2]]) {
          return true;
        }
      }
    }
  }
  return false;
}

template <class Bonded>
bool basalConditionsOf(Bonded &&bonded, const std::vector<int> &basal1,
                       const std::vector<int> &basal2) {
  if (basal1.size() != 6 || basal2.size() != 6) {
    return false;
  }
  bool bond[6][6] = {};
  for (int i = 0; i < 6; i++) {
    const int left = basal1[static_cast<size_t>(i)];
    for (int j = 0; j < 6; j++) {
      bond[i][j] = bonded(left, basal2[static_cast<size_t>(j)]);
    }
  }
  return basalFromMatrix(bond);
}

int sharedHex(const int *x, const int *y) {
  int n = 0;
  for (int i = 0; i < 6; i++) {
    const int v = x[i];
    n += (v == y[0]) + (v == y[1]) + (v == y[2]) + (v == y[3]) + (v == y[4]) +
         (v == y[5]);
  }
  return n;
}

int sharedAtoms(const std::vector<int> &a, const std::vector<int> &b) {
  if (a.size() == 6 && b.size() == 6) {
    return sharedHex(a.data(), b.data());
  }
  int n = 0;
  for (const int x : a) {
    n += static_cast<int>(std::find(b.begin(), b.end(), x) != b.end());
  }
  return n;
}

bool commonAtom(const std::vector<int> &a, const std::vector<int> &b,
                const std::vector<int> &c) {
  if (a.size() == 6 && b.size() == 6 && c.size() == 6) {
    const int *x = a.data();
    const int *y = b.data();
    const int *z = c.data();
    for (int i = 0; i < 6; i++) {
      const int v = x[i];
      const bool inB = v == y[0] || v == y[1] || v == y[2] || v == y[3] ||
                       v == y[4] || v == y[5];
      if (!inB) {
        continue;
      }
      if (v == z[0] || v == z[1] || v == z[2] || v == z[3] || v == z[4] ||
          v == z[5]) {
        return true;
      }
    }
    return false;
  }
  for (const int x : a) {
    if (std::find(b.begin(), b.end(), x) != b.end() &&
        std::find(c.begin(), c.end(), x) != c.end()) {
      return true;
    }
  }
  return false;
}

struct RingRow {
  const int *p = nullptr;
  int n = 0;
  [[nodiscard]] int size() const { return n; }
  [[nodiscard]] bool empty() const { return n == 0; }
  [[nodiscard]] int operator[](int i) const { return p[i]; }
  [[nodiscard]] int operator[](std::size_t i) const { return p[static_cast<int>(i)]; }
  [[nodiscard]] const int *begin() const { return p; }
  [[nodiscard]] const int *end() const { return p + n; }
};

RingRow atomRow(const ring::RingSearchIndex &index, int atom) {
  if (atom < 0 || index.ringOff.empty() ||
      static_cast<std::size_t>(atom) + 1 >= index.ringOff.size()) {
    return {};
  }
  const int begin = index.ringOff[static_cast<std::size_t>(atom)];
  const int end = index.ringOff[static_cast<std::size_t>(atom) + 1];
  return {index.ringCol.data() + begin, end - begin};
}

// Six triplet lists per ring stay on the stack. Ice puts one ring in each.
struct IntBuf {
  int stack[16]{};
  int nStack = 0;
  std::vector<int> extra;
  void push_back(int x) {
    if (nStack < 16) {
      stack[nStack++] = x;
    } else {
      extra.push_back(x);
    }
  }
  void clear() {
    nStack = 0;
    extra.clear();
  }
  [[nodiscard]] int size() const {
    return nStack + static_cast<int>(extra.size());
  }
  [[nodiscard]] bool empty() const { return size() == 0; }
  [[nodiscard]] int operator[](int i) const {
    return i < nStack ? stack[i] : extra[static_cast<size_t>(i - nStack)];
  }
  [[nodiscard]] int operator[](std::size_t i) const {
    return operator[](static_cast<int>(i));
  }
};

// Every ring that contains all three atoms, except skipA and skipB.
template <class Out>
void ringsThrough(const ring::RingSearchIndex &index, int a, int b, int c,
                  int skipA, int skipB, Out &out) {
  const auto &A = atomRow(index, a);
  const auto &B = atomRow(index, b);
  const auto &C = atomRow(index, c);
  size_t i = 0;
  size_t j = 0;
  size_t k = 0;
  while (i < A.size() && j < B.size() && k < C.size()) {
    const int x = A[i];
    const int y = B[j];
    const int z = C[k];
    if (x == y && y == z) {
      if (x != skipA && x != skipB) {
        out.push_back(x);
      }
      ++i;
      ++j;
      ++k;
      continue;
    }
    const int lo = std::min(x, std::min(y, z));
    if (x == lo) {
      ++i;
    }
    if (y == lo) {
      ++j;
    }
    if (z == lo) {
      ++k;
    }
  }
}

// A partner ring on an adjacent triplet already shares the atoms of that
// triplet. share>=3 is then a count of the remaining atoms, so the pairs
// are the edges of a sparse graph rather than the full product of the
// three triplet sets. Chains and the assignment chosen from them are the
// same ones the triple loop would accept.
struct DdcChain {
  std::array<int, 3> r;
  unsigned middles;
};

struct PartnerScratch {
  std::vector<std::vector<int>> bucket;
  std::vector<int> used;
  std::vector<int> stamp;
  std::vector<int> hits;
  std::vector<int> seen;
  std::vector<int> tripSeen;
  int tick = 0;
  int tripTick = 1;
  int markTickB = 0;
  int markTickC = 0;

  void prepare(int nAtoms, int nRight) {
    if (static_cast<int>(bucket.size()) < nAtoms) {
      bucket.resize(static_cast<size_t>(nAtoms));
    }
    if (static_cast<int>(tripSeen.size()) < nAtoms) {
      tripSeen.assign(static_cast<size_t>(nAtoms), 0);
    }
    if (static_cast<int>(stamp.size()) < nRight) {
      stamp.resize(static_cast<size_t>(nRight), 0);
      hits.resize(static_cast<size_t>(nRight), 0);
    }
  }

  void clearBuckets() {
    for (const int atom : used) {
      bucket[static_cast<size_t>(atom)].clear();
    }
    used.clear();
  }
};

void partnersOf(const std::vector<std::vector<int>> &rings,
                const IntBuf &left, const IntBuf &right,
                const int *tripRight, int nAtoms, PartnerScratch &scratch,
                std::vector<std::vector<int>> &adj) {
  const int nRight = static_cast<int>(right.size());
  adj.assign(static_cast<std::size_t>(left.size()), {});
  scratch.prepare(nAtoms, nRight);
  scratch.clearBuckets();
  if (++scratch.tripTick <= 0) {
    std::fill(scratch.tripSeen.begin(), scratch.tripSeen.end(), 0);
    scratch.tripTick = 1;
  }
  const int tripMark = scratch.tripTick;
  for (int t = 0; t < 3; ++t) {
    const int atom = tripRight[t];
    if (atom >= 0 && atom < nAtoms) {
      scratch.tripSeen[static_cast<size_t>(atom)] = tripMark;
    }
  }
  for (int j = 0; j < nRight; ++j) {
    const int ring = right[static_cast<size_t>(j)];
    for (const int atom : rings[static_cast<size_t>(ring)]) {
      if (atom < 0 || atom >= nAtoms ||
          scratch.tripSeen[static_cast<size_t>(atom)] == tripMark) {
        continue;
      }
      auto &bucket = scratch.bucket[static_cast<size_t>(atom)];
      if (bucket.empty()) {
        scratch.used.push_back(atom);
      }
      bucket.push_back(j);
    }
  }
  for (int i = 0; i < static_cast<int>(left.size()); ++i) {
    const int ra = left[static_cast<size_t>(i)];
    const auto &atoms = rings[static_cast<size_t>(ra)];
    int always = 0;
    for (const int atom : atoms) {
      if (atom >= 0 && atom < nAtoms &&
          scratch.tripSeen[static_cast<size_t>(atom)] == tripMark) {
        ++always;
      }
    }
    if (always >= 3) {
      auto &out = adj[static_cast<size_t>(i)];
      out.reserve(static_cast<size_t>(nRight));
      for (int j = 0; j < nRight; ++j) {
        if (right[static_cast<size_t>(j)] != ra) {
          out.push_back(j);
        }
      }
      continue;
    }
    if (++scratch.tick <= 0) {
      std::fill(scratch.stamp.begin(), scratch.stamp.end(), 0);
      scratch.tick = 1;
    }
    scratch.seen.clear();
    for (const int atom : atoms) {
      if (atom < 0 || atom >= nAtoms ||
          scratch.tripSeen[static_cast<size_t>(atom)] == tripMark) {
        continue;
      }
      for (const int j : scratch.bucket[static_cast<size_t>(atom)]) {
        if (scratch.stamp[static_cast<size_t>(j)] != scratch.tick) {
          scratch.stamp[static_cast<size_t>(j)] = scratch.tick;
          scratch.hits[static_cast<size_t>(j)] = 0;
          scratch.seen.push_back(j);
        }
        ++scratch.hits[static_cast<size_t>(j)];
      }
    }
    auto &out = adj[static_cast<size_t>(i)];
    for (const int j : scratch.seen) {
      if (always + scratch.hits[static_cast<size_t>(j)] >= 3 &&
          right[static_cast<size_t>(j)] != ra) {
        out.push_back(j);
      }
    }
  }
}

unsigned ddcPartnerMask(unsigned middles) {
  unsigned partners = 0U;
  for (unsigned t = 0; t < 3; ++t) {
    if ((middles >> t) & 1U) {
      partners |= (1U << t) | (1U << ((t + 2U) % 3U));
    }
  }
  return partners;
}

bool ddcDistinct(const DdcChain &e, const DdcChain &o) {
  for (const int x : e.r) {
    for (const int y : o.r) {
      if (x == y) {
        return false;
      }
    }
  }
  return true;
}

// Counts of rings inside one middle-mask of the opposite parity. A stamp
// replaces the hash map: only the rings touched by this mask are written.
struct DdcMaskCount {
  int tick = 1;
  int total = 0;
  std::vector<int> touch;
  std::vector<int> cnt;
  std::vector<std::vector<int>> at;
  std::vector<int> used;
  std::vector<int> chainSeen;
  int chainTick = 1;

  void prepare(int nRings) {
    if (static_cast<int>(touch.size()) < nRings) {
      touch.assign(static_cast<size_t>(nRings), 0);
      cnt.assign(static_cast<size_t>(nRings), 0);
      at.resize(static_cast<size_t>(nRings));
    }
  }

  void reset() {
    for (const int ring : used) {
      at[static_cast<size_t>(ring)].clear();
    }
    used.clear();
    if (++tick <= 0) {
      std::fill(touch.begin(), touch.end(), 0);
      tick = 1;
    }
    total = 0;
  }

  void add(int ring) {
    if (ring < 0 || ring >= static_cast<int>(touch.size())) {
      return;
    }
    if (touch[static_cast<size_t>(ring)] != tick) {
      touch[static_cast<size_t>(ring)] = tick;
      cnt[static_cast<size_t>(ring)] = 0;
    }
    ++cnt[static_cast<size_t>(ring)];
  }

  [[nodiscard]] int get(int ring) const {
    if (ring < 0 || ring >= static_cast<int>(touch.size()) ||
        touch[static_cast<size_t>(ring)] != tick) {
      return 0;
    }
    return cnt[static_cast<size_t>(ring)];
  }

  void lists(const std::vector<int> &ids, const std::vector<DdcChain> &chains) {
    for (const int id : ids) {
      for (const int ring : chains[static_cast<size_t>(id)].r) {
        if (ring < 0 || ring >= static_cast<int>(at.size())) {
          continue;
        }
        auto &row = at[static_cast<size_t>(ring)];
        if (row.empty()) {
          used.push_back(ring);
        }
        row.push_back(id);
      }
    }
  }

  void prepareChains(int nChains) {
    if (static_cast<int>(chainSeen.size()) < nChains) {
      chainSeen.assign(static_cast<size_t>(nChains), 0);
      chainTick = 1;
    }
  }

  // True when some chain of this mask is disjoint from {a, b, c}.
  // Chain ids are indexes into the opposite chain vector.
  [[nodiscard]] bool hasDisjoint(int a, int b, int c) {
    if (++chainTick <= 0) {
      std::fill(chainSeen.begin(), chainSeen.end(), 0);
      chainTick = 1;
    }
    int cover = 0;
    const int ringsHere[3] = {a, b, c};
    for (const int ring : ringsHere) {
      if (ring < 0 || ring >= static_cast<int>(at.size())) {
        continue;
      }
      for (const int id : at[static_cast<size_t>(ring)]) {
        if (id < 0 || id >= static_cast<int>(chainSeen.size())) {
          continue;
        }
        if (chainSeen[static_cast<size_t>(id)] == chainTick) {
          continue;
        }
        chainSeen[static_cast<size_t>(id)] = chainTick;
        ++cover;
      }
    }
    return total > cover;
  }
};

void markDdcLive(const std::vector<DdcChain> &self,
                 const std::vector<DdcChain> &other, int nRings,
                 std::vector<char> &live, DdcMaskCount &count) {
  count.prepare(nRings);
  count.prepareChains(static_cast<int>(other.size()));
  live.assign(self.size(), 0);
  std::array<std::vector<int>, 8> bucket;
  for (int i = 0; i < static_cast<int>(other.size()); ++i) {
    bucket[other[static_cast<size_t>(i)].middles & 7U].push_back(i);
  }
  std::vector<int> slow;
  for (unsigned mask = 1; mask <= 7; ++mask) {
    const std::vector<int> &ids = bucket[mask];
    if (ids.empty()) {
      continue;
    }
    count.reset();
    count.total = static_cast<int>(ids.size());
    for (const int id : ids) {
      for (const int ring : other[static_cast<size_t>(id)].r) {
        count.add(ring);
      }
    }
    slow.clear();
    bool needLists = false;
    for (int i = 0; i < static_cast<int>(self.size()); ++i) {
      if (live[static_cast<size_t>(i)] != 0) {
        continue;
      }
      const unsigned partners =
          ddcPartnerMask(self[static_cast<size_t>(i)].middles);
      if ((partners & mask) == 0U) {
        continue;
      }
      const auto &chain = self[static_cast<size_t>(i)];
      const int coverBound = count.get(chain.r[0]) + count.get(chain.r[1]) +
                             count.get(chain.r[2]);
      if (count.total > coverBound) {
        live[static_cast<size_t>(i)] = 1;
      } else {
        slow.push_back(i);
        needLists = true;
      }
    }
    if (!needLists) {
      continue;
    }
    count.lists(ids, other);
    for (const int i : slow) {
      if (live[static_cast<size_t>(i)] != 0) {
        continue;
      }
      const auto &chain = self[static_cast<size_t>(i)];
      if (count.hasDisjoint(chain.r[0], chain.r[1], chain.r[2])) {
        live[static_cast<size_t>(i)] = 1;
      }
    }
  }
}

void ddcChains(const std::vector<std::vector<int>> &rings,
               const std::vector<int> &equatorial,
               const std::array<IntBuf, 6> &T, int parity, int nAtoms,
               std::vector<DdcChain> &out) {
  thread_local PartnerScratch scratch;
  thread_local std::vector<std::vector<int>> ab;
  thread_local std::vector<std::vector<int>> ac;
  thread_local std::vector<std::vector<int>> bc;
  thread_local std::vector<std::vector<int>> bcFromC;
  thread_local std::vector<int> markB;
  thread_local std::vector<int> markC;
  out.clear();
  const IntBuf &A = T[static_cast<size_t>(parity)];
  const IntBuf &B = T[static_cast<size_t>(parity + 2)];
  const IntBuf &C = T[static_cast<size_t>(parity + 4)];
  int tripB[3];
  int tripC[3];
  for (int t = 0; t < 3; ++t) {
    tripB[t] = equatorial[static_cast<size_t>((parity + 2 + t) % 6)];
    tripC[t] = equatorial[static_cast<size_t>((parity + 4 + t) % 6)];
  }
  partnersOf(rings, A, B, tripB, nAtoms, scratch, ab);
  partnersOf(rings, A, C, tripC, nAtoms, scratch, ac);
  partnersOf(rings, B, C, tripC, nAtoms, scratch, bc);
  for (auto &row : bc) {
    std::sort(row.begin(), row.end());
  }
  if (bcFromC.size() < C.size()) {
    bcFromC.resize(static_cast<size_t>(C.size()));
  }
  for (int ci = 0; ci < C.size(); ++ci) {
    bcFromC[static_cast<size_t>(ci)].clear();
  }
  for (int j = 0; j < static_cast<int>(B.size()); ++j) {
    for (const int kc : bc[static_cast<size_t>(j)]) {
      bcFromC[static_cast<size_t>(kc)].push_back(j);
    }
  }
  if (markB.size() < static_cast<size_t>(B.size())) {
    markB.resize(static_cast<size_t>(B.size()));
  }
  if (markC.size() < static_cast<size_t>(C.size())) {
    markC.resize(static_cast<size_t>(C.size()));
  }
  int tickB = scratch.markTickB;
  int tickC = scratch.markTickC;
  auto fresh = [](int &tick, std::vector<int> &mark) {
    if (++tick <= 0) {
      std::fill(mark.begin(), mark.end(), 0);
      tick = 1;
    }
    return tick;
  };
  for (int i = 0; i < static_cast<int>(A.size()); ++i) {
    const int onB = fresh(tickB, markB);
    const int onC = fresh(tickC, markC);
    for (const int jb : ab[static_cast<size_t>(i)]) {
      markB[static_cast<size_t>(jb)] = onB;
    }
    for (const int kc : ac[static_cast<size_t>(i)]) {
      markC[static_cast<size_t>(kc)] = onC;
    }
    const int ra = A[static_cast<size_t>(i)];
    const auto &atomsA = rings[static_cast<size_t>(ra)];
    for (const int jb : ab[static_cast<size_t>(i)]) {
      const int rb = B[static_cast<size_t>(jb)];
      if (rb == ra) {
        continue;
      }
      const auto &atomsB = rings[static_cast<size_t>(rb)];
      for (const int kc : ac[static_cast<size_t>(i)]) {
        const int rc = C[static_cast<size_t>(kc)];
        if (rc == ra || rc == rb) {
          continue;
        }
        const bool bcOK =
            std::binary_search(bc[static_cast<size_t>(jb)].begin(),
                               bc[static_cast<size_t>(jb)].end(), kc);
        const unsigned middles = 1U | (bcOK ? 2U : 0U) | (bcOK ? 4U : 0U);
        if (!commonAtom(atomsA, atomsB, rings[static_cast<size_t>(rc)])) {
          continue;
        }
        out.push_back(DdcChain{{ra, rb, rc}, middles});
      }
      for (const int kc : bc[static_cast<size_t>(jb)]) {
        if (markC[static_cast<size_t>(kc)] == onC) {
          continue;
        }
        const int rc = C[static_cast<size_t>(kc)];
        if (rc == ra || rc == rb) {
          continue;
        }
        if (!commonAtom(atomsA, atomsB, rings[static_cast<size_t>(rc)])) {
          continue;
        }
        out.push_back(DdcChain{{ra, rb, rc}, 2U});
      }
    }
    for (const int kc : ac[static_cast<size_t>(i)]) {
      const int rc = C[static_cast<size_t>(kc)];
      if (rc == ra) {
        continue;
      }
      const auto &atomsC = rings[static_cast<size_t>(rc)];
      for (const int jb : bcFromC[static_cast<size_t>(kc)]) {
        if (markB[static_cast<size_t>(jb)] == onB) {
          continue;
        }
        const int rb = B[static_cast<size_t>(jb)];
        if (rb == ra || rb == rc) {
          continue;
        }
        if (!commonAtom(atomsA, rings[static_cast<size_t>(rb)], atomsC)) {
          continue;
        }
        out.push_back(DdcChain{{ra, rb, rc}, 4U});
      }
    }
  }
  scratch.markTickB = tickB;
  scratch.markTickC = tickC;
}

// Seven hexagons of a DDC span 14 molecules, and each peripheral shares
// exactly three atoms with the equator. The stamp is the same test as the
// sort in the triple loop.
bool ddcFourteen(const std::vector<std::vector<int>> &rings,
                 const std::vector<int> &equatorial, const DdcChain &e,
                 const DdcChain &o, int nAtoms) {
  thread_local std::vector<int> stamp;
  thread_local int tick = 1;
  if (static_cast<int>(stamp.size()) < nAtoms) {
    stamp.assign(static_cast<size_t>(nAtoms), 0);
  }
  auto fresh = [&]() {
    if (++tick <= 0) {
      std::fill(stamp.begin(), stamp.end(), 0);
      tick = 1;
    }
    return tick;
  };
  const int eqMark = fresh();
  int eqN = 0;
  for (const int atom : equatorial) {
    if (atom < 0 || atom >= nAtoms || stamp[static_cast<size_t>(atom)] == eqMark) {
      return false;
    }
    stamp[static_cast<size_t>(atom)] = eqMark;
    ++eqN;
  }
  if (eqN != 6) {
    return false;
  }
  const int ids[6] = {e.r[0], e.r[1], e.r[2], o.r[0], o.r[1], o.r[2]};
  for (const int ri : ids) {
    const auto &cycle = rings[static_cast<size_t>(ri)];
    if (cycle.size() != 6) {
      return false;
    }
    int share = 0;
    for (const int atom : cycle) {
      if (atom < 0 || atom >= nAtoms) {
        return false;
      }
      share += stamp[static_cast<size_t>(atom)] == eqMark ? 1 : 0;
    }
    if (share != 3) {
      return false;
    }
  }
  const int allMark = fresh();
  int unique = 0;
  auto take = [&](const std::vector<int> &cycle) {
    for (const int atom : cycle) {
      if (stamp[static_cast<size_t>(atom)] != allMark) {
        stamp[static_cast<size_t>(atom)] = allMark;
        ++unique;
      }
    }
  };
  take(equatorial);
  for (const int ri : ids) {
    take(rings[static_cast<size_t>(ri)]);
  }
  return unique == 14;
}

// One chain plus the equator is a subset of the seven rings, so more than
// 14 molecules cannot be completed to a DDC.
bool ddcChainRoom(const std::vector<std::vector<int>> &rings,
                  const std::vector<int> &equatorial, const DdcChain &chain,
                  int nAtoms) {
  thread_local std::vector<int> stamp;
  thread_local int tick = 1;
  if (static_cast<int>(stamp.size()) < nAtoms) {
    stamp.assign(static_cast<size_t>(nAtoms), 0);
  }
  if (++tick <= 0) {
    std::fill(stamp.begin(), stamp.end(), 0);
    tick = 1;
  }
  const int mark = tick;
  int unique = 0;
  auto take = [&](const std::vector<int> &cycle) -> bool {
    for (const int atom : cycle) {
      if (atom < 0 || atom >= nAtoms) {
        return false;
      }
      if (stamp[static_cast<size_t>(atom)] != mark) {
        stamp[static_cast<size_t>(atom)] = mark;
        ++unique;
        if (unique > 14) {
          return false;
        }
      }
    }
    return true;
  };
  if (!take(equatorial)) {
    return false;
  }
  for (const int ri : chain.r) {
    if (!take(rings[static_cast<size_t>(ri)])) {
      return false;
    }
  }
  return true;
}

void ddcKeepRoom(const std::vector<std::vector<int>> &rings,
                 const std::vector<int> &equatorial, int nAtoms,
                 std::vector<DdcChain> &chains) {
  size_t w = 0;
  for (size_t i = 0; i < chains.size(); ++i) {
    if (ddcChainRoom(rings, equatorial, chains[i], nAtoms)) {
      if (w != i) {
        chains[w] = chains[i];
      }
      ++w;
    }
  }
  chains.resize(w);
}

[[gnu::noinline, gnu::cold]] bool equatorialLarge(const std::vector<std::vector<int>> &rings,
                     const std::vector<int> &equatorial,
                     const std::array<IntBuf, 6> &T, int nAtoms,
                     std::vector<int> *peripheralUnion,
                     std::vector<int> *chosen, const std::vector<int> *rank,
                     ring::EquatorChoice *raw) {
  // A peripheral that meets the equator in a fourth atom is not a DDC face.
  // Dropping it here leaves the pairs the triple loop would keep.
  thread_local std::vector<int> eqStamp;
  thread_local int eqTick = 1;
  if (static_cast<int>(eqStamp.size()) < nAtoms) {
    eqStamp.assign(static_cast<size_t>(nAtoms), 0);
  }
  if (++eqTick <= 0) {
    std::fill(eqStamp.begin(), eqStamp.end(), 0);
    eqTick = 1;
  }
  const int eqMark = eqTick;
  for (const int atom : equatorial) {
    if (atom >= 0 && atom < nAtoms) {
      eqStamp[static_cast<size_t>(atom)] = eqMark;
    }
  }
  thread_local std::array<IntBuf, 6> kept;
  for (int k = 0; k < 6; ++k) {
    kept[static_cast<size_t>(k)].clear();
    const IntBuf &src = T[static_cast<size_t>(k)];
    for (int i = 0; i < src.size(); ++i) {
      const int ring = src[static_cast<size_t>(i)];
      const auto &cycle = rings[static_cast<size_t>(ring)];
      if (cycle.size() != 6) {
        continue;
      }
      int share = 0;
      for (const int atom : cycle) {
        if (atom >= 0 && atom < nAtoms &&
            eqStamp[static_cast<size_t>(atom)] == eqMark) {
          ++share;
        }
      }
      if (share == 3) {
        kept[static_cast<size_t>(k)].push_back(ring);
      }
    }
    if (kept[static_cast<size_t>(k)].empty()) {
      return false;
    }
  }
  std::vector<DdcChain> even;
  std::vector<DdcChain> odd;
  ddcChains(rings, equatorial, kept, 0, nAtoms, even);
  if (even.empty()) {
    return false;
  }
  ddcChains(rings, equatorial, kept, 1, nAtoms, odd);
  if (odd.empty()) {
    return false;
  }
  ddcKeepRoom(rings, equatorial, nAtoms, even);
  ddcKeepRoom(rings, equatorial, nAtoms, odd);
  if (even.empty() || odd.empty()) {
    return false;
  }
  thread_local DdcMaskCount maskCount;
  std::vector<char> evenLive;
  std::vector<char> oddLive;
  const int nRings = static_cast<int>(rings.size());
  markDdcLive(even, odd, nRings, evenLive, maskCount);
  bool any = false;
  for (const char flag : evenLive) {
    any = any || flag != 0;
  }
  if (!any) {
    return false;
  }
  markDdcLive(odd, even, nRings, oddLive, maskCount);
  any = false;
  for (const char flag : oddLive) {
    any = any || flag != 0;
  }
  if (!any) {
    return false;
  }
  const bool wantBest =
      rank != nullptr && (chosen != nullptr || raw != nullptr);
  auto pairIsKept = [&](const DdcChain &e, const DdcChain &o) {
    return (ddcPartnerMask(e.middles) & o.middles) != 0U && ddcDistinct(e, o) &&
           ddcFourteen(rings, equatorial, e, o, nAtoms);
  };
  if (!wantBest && peripheralUnion == nullptr) {
    for (size_t i = 0; i < even.size(); ++i) {
      if (evenLive[static_cast<size_t>(i)] == 0) {
        continue;
      }
      for (size_t j = 0; j < odd.size(); ++j) {
        if (oddLive[static_cast<size_t>(j)] == 0) {
          continue;
        }
        if (pairIsKept(even[i], odd[j])) {
          return true;
        }
      }
    }
    return false;
  }
  if (peripheralUnion != nullptr) {
    for (size_t i = 0; i < even.size(); ++i) {
      if (evenLive[static_cast<size_t>(i)] == 0) {
        continue;
      }
      for (size_t j = 0; j < odd.size(); ++j) {
        if (oddLive[static_cast<size_t>(j)] == 0) {
          continue;
        }
        if (!pairIsKept(even[i], odd[j])) {
          continue;
        }
        peripheralUnion->insert(peripheralUnion->end(), even[i].r.begin(),
                                even[i].r.end());
        peripheralUnion->insert(peripheralUnion->end(), odd[j].r.begin(),
                                odd[j].r.end());
      }
    }
    std::sort(peripheralUnion->begin(), peripheralUnion->end());
    peripheralUnion->erase(
        std::unique(peripheralUnion->begin(), peripheralUnion->end()),
        peripheralUnion->end());
    if (peripheralUnion->empty() && !wantBest) {
      return false;
    }
  }
  if (wantBest) {
    int bestRank = std::numeric_limits<int>::max();
    int pivot = -1;
    auto consider = [&](int ring) {
      const int rk = (*rank)[static_cast<size_t>(ring)];
      if (rk < bestRank) {
        bestRank = rk;
        pivot = ring;
      }
    };
    for (size_t i = 0; i < even.size(); ++i) {
      if (evenLive[static_cast<size_t>(i)] == 0) {
        continue;
      }
      for (const int ring : even[static_cast<size_t>(i)].r) {
        consider(ring);
      }
    }
    for (size_t i = 0; i < odd.size(); ++i) {
      if (oddLive[static_cast<size_t>(i)] == 0) {
        continue;
      }
      for (const int ring : odd[static_cast<size_t>(i)].r) {
        consider(ring);
      }
    }
    std::vector<int> bestKey;
    std::vector<int> bestRings;
    auto tryPair = [&](const DdcChain &e, const DdcChain &o) {
      if ((ddcPartnerMask(e.middles) & o.middles) == 0U || !ddcDistinct(e, o)) {
        return;
      }
      int key[6];
      key[0] = (*rank)[static_cast<size_t>(e.r[0])];
      key[1] = (*rank)[static_cast<size_t>(e.r[1])];
      key[2] = (*rank)[static_cast<size_t>(e.r[2])];
      key[3] = (*rank)[static_cast<size_t>(o.r[0])];
      key[4] = (*rank)[static_cast<size_t>(o.r[1])];
      key[5] = (*rank)[static_cast<size_t>(o.r[2])];
      int first = key[0];
      int second = key[1];
      if (second < first) {
        std::swap(first, second);
      }
      for (int t = 2; t < 6; ++t) {
        const int v = key[t];
        if (v < first) {
          second = first;
          first = v;
        } else if (v < second) {
          second = v;
        }
      }
      if (!bestKey.empty() && second > bestKey[1]) {
        return;
      }
      for (int i = 1; i < 6; ++i) {
        const int v = key[i];
        int j = i;
        while (j > 0 && key[j - 1] > v) {
          key[j] = key[j - 1];
          --j;
        }
        key[j] = v;
      }
      if (!bestKey.empty()) {
        for (int i = 0; i < 6; ++i) {
          if (key[i] != bestKey[static_cast<size_t>(i)]) {
            if (key[i] > bestKey[static_cast<size_t>(i)]) {
              return;
            }
            break;
          }
          if (i == 5) {
            return;
          }
        }
      }
      if (!ddcFourteen(rings, equatorial, e, o, nAtoms)) {
        return;
      }
      bestKey.assign(key, key + 6);
      bestRings.clear();
      bestRings.insert(bestRings.end(), e.r.begin(), e.r.end());
      bestRings.insert(bestRings.end(), o.r.begin(), o.r.end());
    };
    auto threeMin = [&](const DdcChain &chain) {
      const int a = (*rank)[static_cast<size_t>(chain.r[0])];
      const int b = (*rank)[static_cast<size_t>(chain.r[1])];
      const int c = (*rank)[static_cast<size_t>(chain.r[2])];
      return std::min(a, std::min(b, c));
    };
    auto companionOf = [&](const DdcChain &chain) {
      int m = std::numeric_limits<int>::max();
      for (const int ring : chain.r) {
        if (ring == pivot) {
          continue;
        }
        m = std::min(m, (*rank)[static_cast<size_t>(ring)]);
      }
      return m;
    };
    std::array<std::vector<int>, 8> evenBucket;
    std::array<std::vector<int>, 8> oddBucket;
    std::vector<int> evenHit;
    std::vector<int> oddHit;
    evenHit.reserve(even.size() / 8 + 1);
    oddHit.reserve(odd.size() / 8 + 1);
    for (int i = 0; i < static_cast<int>(even.size()); ++i) {
      if (evenLive[static_cast<size_t>(i)] == 0) {
        continue;
      }
      evenBucket[even[static_cast<size_t>(i)].middles & 7U].push_back(i);
      for (const int ring : even[static_cast<size_t>(i)].r) {
        if (ring == pivot) {
          evenHit.push_back(i);
          break;
        }
      }
    }
    for (int i = 0; i < static_cast<int>(odd.size()); ++i) {
      if (oddLive[static_cast<size_t>(i)] == 0) {
        continue;
      }
      oddBucket[odd[static_cast<size_t>(i)].middles & 7U].push_back(i);
      for (const int ring : odd[static_cast<size_t>(i)].r) {
        if (ring == pivot) {
          oddHit.push_back(i);
          break;
        }
      }
    }
    std::vector<int> evenMin(even.size());
    std::vector<int> oddMin(odd.size());
    std::vector<int> evenGate(even.size());
    std::vector<int> oddGate(odd.size());
    for (unsigned mask = 1; mask <= 7; ++mask) {
      for (const int i : evenBucket[mask]) {
        evenMin[static_cast<size_t>(i)] =
            threeMin(even[static_cast<size_t>(i)]);
      }
      for (const int i : oddBucket[mask]) {
        oddMin[static_cast<size_t>(i)] = threeMin(odd[static_cast<size_t>(i)]);
      }
      std::sort(evenBucket[mask].begin(), evenBucket[mask].end(),
                [&](int a, int b) {
                  return evenMin[static_cast<size_t>(a)] <
                         evenMin[static_cast<size_t>(b)];
                });
      std::sort(oddBucket[mask].begin(), oddBucket[mask].end(),
                [&](int a, int b) {
                  return oddMin[static_cast<size_t>(a)] <
                         oddMin[static_cast<size_t>(b)];
                });
    }
    for (const int i : evenHit) {
      evenGate[static_cast<size_t>(i)] =
          companionOf(even[static_cast<size_t>(i)]);
    }
    for (const int i : oddHit) {
      oddGate[static_cast<size_t>(i)] = companionOf(odd[static_cast<size_t>(i)]);
    }
    std::sort(evenHit.begin(), evenHit.end(), [&](int a, int b) {
      return evenGate[static_cast<size_t>(a)] < evenGate[static_cast<size_t>(b)];
    });
    std::sort(oddHit.begin(), oddHit.end(), [&](int a, int b) {
      return oddGate[static_cast<size_t>(a)] < oddGate[static_cast<size_t>(b)];
    });
    for (const int i : evenHit) {
      const int gate = evenGate[static_cast<size_t>(i)];
      const unsigned partners =
          ddcPartnerMask(even[static_cast<size_t>(i)].middles);
      for (unsigned mask = 1; mask <= 7; ++mask) {
        if ((partners & mask) == 0U) {
          continue;
        }
        for (const int j : oddBucket[mask]) {
          if (!bestKey.empty() && gate > bestKey[1] &&
              oddMin[static_cast<size_t>(j)] > bestKey[1]) {
            break;
          }
          tryPair(even[static_cast<size_t>(i)], odd[static_cast<size_t>(j)]);
        }
      }
    }
    for (const int j : oddHit) {
      const int gate = oddGate[static_cast<size_t>(j)];
      const unsigned oddMiddles = odd[static_cast<size_t>(j)].middles;
      for (unsigned mask = 1; mask <= 7; ++mask) {
        // tryPair tests partner(even) & odd. The mask that reaches an even
        // chain is the even middles, so the odd chain has to meet
        // partner(that mask), not partner(odd) & mask.
        if ((ddcPartnerMask(mask) & oddMiddles) == 0U) {
          continue;
        }
        for (const int i : evenBucket[mask]) {
          if (!bestKey.empty() && gate > bestKey[1] &&
              evenMin[static_cast<size_t>(i)] > bestKey[1]) {
            break;
          }
          tryPair(even[static_cast<size_t>(i)], odd[static_cast<size_t>(j)]);
        }
      }
    }
    if (bestRings.empty()) {
      return false;
    }
    int ids[6];
    for (int t = 0; t < 6; ++t) {
      ids[t] = bestRings[static_cast<size_t>(t)];
    }
    std::sort(ids, ids + 6, [rank](int a, int b) {
      return (*rank)[static_cast<size_t>(a)] < (*rank)[static_cast<size_t>(b)];
    });
    if (chosen != nullptr) {
      chosen->assign(ids, ids + 6);
    }
    if (raw != nullptr) {
      for (int t = 0; t < 6; ++t) {
        raw->ring[t] = ids[t];
      }
      raw->n = 6;
    }
  }
  return true;
}

} // namespace

// -----------------------------------------------------------------------------------------------------
// BULK RING SEARCH ONLY
// -----------------------------------------------------------------------------------------------------

/**
 * @details Function that loops through the primitive rings (which is a vector
 * of vectors) of all sizes, upto maxDepth (the largest ring size). 
 * - ring::clearRingList (Clears the vector of vectors for rings of a single
 * type, to prevent excessive memory being blocked).
 * - ring::getSingleRingSize (Fill a vector of vectors for rings of a particular
 * ring size).
 * - ring::findPrisms (Now that rings of a particular size have been obtained,
 * prism blocks are found and saved).
 * - topoparam::normHeightPercent (Gets the height% for the prism blocks).
 * - ring::assignPrismType (Assigns a prism type to each atom type).
 * - sout::writePrismNum (Write out the prism information for the current
 * frame).
 * - sout::writeLAMMPSdataAllPrisms (Writes out a LAMMPS data file for the
 * current frame, which can be visualized in OVITO).
 * @param[in] path The string to the output directory, in which files will be
 *  written out.
 * @param[in] rings Row-ordered vector of vectors for rings of a single type.
 * @param[in] nList Row-ordered neighbour list by index.
 * @param[in] yCloud The input PointCloud.
 * @param[in] maxDepth The maximum possible size of the primitive rings.
 * @param[in] firstFrame The first frame to be analyzed
 */
int ring::bulkPolygonRingAnalysis(
    std::string path, const std::vector<std::vector<int>> &rings,
    const std::vector<std::vector<int>> &nList,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, int maxDepth,
    int firstFrame) {
  //
  std::vector<std::vector<int>>
      ringsOneType;           // Vector of vectors of rings of a single size
  int nRings;                 // Number of rings of the current type
  std::vector<int> nRingList; // Vector of the values of the number of rings
                              // for a particular frame
  std::vector<int>
      atomTypes; // contains int values for each ring type considered
  // -------------------------------------------------------------------------------
  // Init
  nRingList.resize(
      static_cast<std::size_t>(maxDepth - 2)); // Has a value for every value of ringSize from 3, upto maxDepth
  // The atomTypes vector is the same size as the pointCloud atoms
  atomTypes.resize(static_cast<std::size_t>(yCloud.nop), 1); // The dummy or unclassified value is 1
  // -------------------------------------------------------------------------------
  // Run this loop for rings of sizes upto maxDepth
  // The smallest possible ring is of size 3
  for (int ringSize = 3; ringSize <= maxDepth; ringSize++) {
    // Clear ringsOneType
    ring::clearRingList(ringsOneType);
    // Get rings of the current ring size
    ringsOneType = ring::getSingleRingSize(rings, ringSize);
    //
    // Continue if there are zero rings of ringSize
    if (ringsOneType.empty()) {
      nRingList[static_cast<std::size_t>(ringSize - 3)] = 0; // Update the number of prisms
      continue;
    } // skip if there are no rings
    //
    // -------------
    // Number of rings with n nodes
    nRings = static_cast<int>(ringsOneType.size());
    // -------------
    // Now that you have rings of a certain size:
    nRingList[static_cast<std::size_t>(ringSize - 3)] = nRings; // Update the number of n-membered rings
    // -------------
  } // end of loop through every possible ringSize

  // Get the atom types for all the ring types
  ring::assignPolygonType(rings, atomTypes, nRingList);

  // Write out the ring information
  sout::writeRingNumBulk(path, yCloud.currentFrame, nRingList, maxDepth, firstFrame);
  // Write out the lammps data file for the particular frame
  sout::writeLAMMPSdataAllRings(yCloud, nList, atomTypes, maxDepth, path, false);

  return 0;
}

// -----------------------------------------------------------------------------------------------------
// DDC / HC ALGORITHMS
// -----------------------------------------------------------------------------------------------------

/**
 * @details Finds out if rings constitute double-diamond cages or hexagonal
 * cages. Requires a neighbour list (by index) and a vector of vectors of
 * primitive rings which should also be by index. This is registered as a Lua
 * function and is accessible to the user. Internally, this function calls the
 * following functions:
 * - ring::getSingleRingSize (Saves rings of a single ring size into a new
 * vector of vectors, which is subsequently used for finding DDCs, HCs etc).
 * - ring::findDDC (Finds the DDCs).
 * - ring::findHC (Finds the HCs).
 * - ring::findMixedRings (Finds the mixed rings, which are shared by DDCs and
 * HCs both).
 * - ring::getStrucNumbers (Gets the number of structures (DDCs, HCs, mixed
 * rings, basal rings, prismatic rings, to be used for write-outs).
 * - sout::writeTopoBulkData (Writes out the numbers and data obtained for the
 * current frame).
 * - ring::getAtomTypesTopoBulk (Gets the atom type for every atom, to be used
 * for printing out the ice types found).
 * - sout::writeLAMMPSdataTopoBulk (Writes out the atoms, with the classified
 * types, into a LAMMPS data file, which can be visualized in OVITO).
 *  @param[in] path The file path of the output directory to which output files
 * will be written.
 *  @param[in] rings Vector of vectors containing the primitive rings. This
 * contains rings of all sizes.
 *  @param[in] nList Row-ordered neighbour list, by index.
 *  @param[in] yCloud The input PointCloud, with respect to which the indices in
 * the rings and nList vector of vectors have been saved.
 *  @param[in] firstFrame First frame to be analyzed
 *  @param[in] onlyTetrahedral Flag for only finding DDCs and HCs (true) or also
 * finding PNCs (false)
 */
int ring::topoBulkAnalysis(
    std::string path, const std::vector<std::vector<int>> &rings,
    const std::vector<std::vector<int>> &nList,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, int firstFrame,
    bool onlyTetrahedral) {
  //
  // Ring IDs of each type will be saved in these vectors
  std::vector<int> listDDC; // Vector for ring indices of DDC
  std::vector<int> listHC;  // Vector for ring indices of HC
  std::vector<int>
      listMixed; // Vector for ring indices of rings that are both DDC and HC
  std::vector<ring::strucType>
      ringType; // This vector will have a value for each ring inside
                // ringList, of type enum strucType in gen.hpp
  // Make a list of all the DDCs and HCs
  std::vector<cage::Cage> cageList;
  std::vector<std::vector<int>>
      ringsOneType;    // Vector of vectors of rings of a single size
  int initRingSize;    // Todo or not: calculate the PNCs or not
  int maxRingSize = 6; // DDCs and HCs are for 6-membered rings
  std::vector<cage::iceType>
      atomTypes; // This vector will have a value for every atom
  // Number of types
  int numHC, numDDC, mixedRings, prismaticRings, basalRings;

  // TODO: handle shapeMatching
  bool doShapeMatching = true;
  // Qualifier for the RMSD per atom:
  std::vector<double> rmsdPerAtom;
  //

  // test
  // std::cout << "rings"
  //           << "\n";
  // for (int i = 0; i < rings.size(); i++) {
  //   for (int j = 0; j < rings[i].size(); j++) {
  //     std::cout << rings[i][j] << " ";
  //   }
  //   std::cout << "\n";
  // }
  // std::cout << "\n";

  if (onlyTetrahedral) {
    initRingSize = 6;
  } else {
    initRingSize = 5;
  }

  // Init the atom type vector
  atomTypes.resize(static_cast<std::size_t>(yCloud.nop)); // Has a value for each atom
  // Init the rmsd per atom (not used yet)
  rmsdPerAtom.resize(static_cast<std::size_t>(yCloud.nop)); // Has a value for each atom

  // ----------------------------------------------
  // Init
  // Get rings of size 5 or 6.
  for (int ringSize = initRingSize; ringSize <= maxRingSize; ringSize++) {
    // Clear ringsOneType
    ring::clearRingList(ringsOneType);
    // Get rings of the current ring size
    ringsOneType = ring::getSingleRingSize(rings, ringSize);
    // Skip for zero rings
    if (ringsOneType.empty()) {
      continue;
    } // skip for no rings of ringSize
    //
    // Init the ringType vector
    ringType.resize(
        ringsOneType.size()); // Has a value for each ring. init to zero.
    // ----------------------------------------------
    if (ringSize == 6) {
      // Get the cages

      // Find HC rings, saving the ring IDs (starting from 0) to listHC
      const auto ringIndex = ring::buildRingSearchIndex(ringsOneType, yCloud);
      listHC = ring::findHC(ringsOneType, ringType, nList, cageList, ringIndex);

      // Find DDC rings, saving the IDs to listDDC
      listDDC =
          ring::findDDC(ringsOneType, ringType, listHC, cageList, ringIndex);

      // Find rings which are both DDCs and HCs (mixed)
      // A dummy value of -10 in the listDDC and listHC vectors for mixed rings
      listMixed =
          ring::findMixedRings(ringsOneType, ringType, listDDC, listHC);

      // Get the number of structures (DDCs, HCs, mixed rings, basal rings,
      // prismatic rings)
      ring::getStrucNumbers(ringType, cageList, numHC, numDDC, mixedRings,
                            prismaticRings, basalRings);

      // Write out to a file
      sout::writeTopoBulkData(path, yCloud.currentFrame, numHC, numDDC,
                              mixedRings, basalRings, prismaticRings,
                              firstFrame);
      // Gets the atom type for every atom, to be used for printing out the ice
      // types found
      ring::getAtomTypesTopoBulk(ringsOneType, ringType, atomTypes);
    }
    // Finding prismatic blocks
    else {
      // Get the prism block classifications
      prism3::findBulkPrisms(ringsOneType, ringType, nList, yCloud,
                             rmsdPerAtom);
      // Gets the atom type for every atom, to be used for printing out the ice
      // types found
      ring::getAtomTypesTopoBulk(ringsOneType, ringType, atomTypes);
    }
  }

  // Print out the lammps data file with the bonds
  sout::writeLAMMPSdataTopoBulk(yCloud, nList, atomTypes, path);
  // To output the bonds between dummy atoms, uncomment the following line
  // sout::writeLAMMPSdataTopoBulk(yCloud, nList, atomTypes, path, true);

  return 0;
}

/**
 * @details Builds the inverted atom-to-rings index. Rows stay in ascending
 *  ring order, which is the order a linear scan over the ring network would
 *  have visited them in, so callers that substitute a lookup for the scan see
 *  the same sequence.
 * @param[in] rings Vector of vectors containing the rings, by atom index.
 * @param[in] numAtoms Number of atoms in the point cloud.
 * @return The index.
 */
ring::RingSearchIndex
ring::buildRingSearchIndex(const std::vector<std::vector<int>> &rings,
                           int numAtoms) {
  ring::RingSearchIndex index;
  if (numAtoms <= 0) {
    return index;
  }
  index.ringsContainingAtom.resize(static_cast<std::size_t>(numAtoms));
  const int nRings = static_cast<int>(rings.size());
  int nThreads = 1;
#ifdef SEAMS_HAS_OPENMP
  if (!omp_in_parallel()) {
    nThreads = std::max(1, omp_get_max_threads());
  }
#endif
  if (nThreads <= 1) {
    std::vector<int> count(static_cast<std::size_t>(numAtoms), 0);
    for (int iring = 0; iring < nRings; iring++) {
      for (const int atom : rings[static_cast<std::size_t>(iring)]) {
        if (atom >= 0 && atom < numAtoms) {
          ++count[static_cast<std::size_t>(atom)];
        }
      }
    }
    index.ringOff.assign(static_cast<std::size_t>(numAtoms) + 1, 0);
    for (int atom = 0; atom < numAtoms; atom++) {
      index.ringOff[static_cast<std::size_t>(atom) + 1] =
          index.ringOff[static_cast<std::size_t>(atom)] + count[static_cast<std::size_t>(atom)];
    }
    index.ringCol.assign(static_cast<std::size_t>(index.ringOff[static_cast<std::size_t>(numAtoms)]),
                         0);
    std::vector<int> cursor(index.ringOff.begin(), index.ringOff.end() - 1);
    for (int iring = 0; iring < nRings; iring++) {
      for (const int atom : rings[static_cast<std::size_t>(iring)]) {
        if (atom >= 0 && atom < numAtoms) {
          index.ringCol[static_cast<std::size_t>(cursor[static_cast<std::size_t>(atom)]++)] =
              iring;
        }
      }
    }
    for (int atom = 0; atom < numAtoms; atom++) {
      const int begin = index.ringOff[static_cast<std::size_t>(atom)];
      const int end = index.ringOff[static_cast<std::size_t>(atom) + 1];
      auto &row = index.ringsContainingAtom[static_cast<std::size_t>(atom)];
      row.assign(index.ringCol.begin() + begin, index.ringCol.begin() + end);
    }
    return index;
  }

  // Each thread counts, then writes, its own ring block. Concatenating those
  // blocks keeps every atom's row in ascending ring order.
  std::vector<std::vector<int>> counts(static_cast<std::size_t>(nThreads));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel
#endif
  {
    const int tid = omp_get_thread_num();
    counts[static_cast<std::size_t>(tid)].assign(static_cast<std::size_t>(numAtoms), 0);
    auto &count = counts[static_cast<std::size_t>(tid)];
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int iring = 0; iring < nRings; iring++) {
      for (const int atom : rings[static_cast<std::size_t>(iring)]) {
        if (atom >= 0 && atom < numAtoms) {
          count[static_cast<std::size_t>(atom)]++;
        }
      }
    }
  }
  std::vector<int> base(static_cast<std::size_t>(numAtoms), 0);
  int total = 0;
  for (int atom = 0; atom < numAtoms; atom++) {
    base[static_cast<std::size_t>(atom)] = total;
    for (int t = 0; t < nThreads; t++) {
      total += counts[static_cast<std::size_t>(t)][static_cast<std::size_t>(atom)];
    }
  }
  std::vector<int> flat(static_cast<std::size_t>(total));
  std::vector<std::vector<int>> head(static_cast<std::size_t>(nThreads));
  for (int t = 0; t < nThreads; t++) {
    head[static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(numAtoms));
  }
  for (int atom = 0; atom < numAtoms; atom++) {
    int cursor = base[static_cast<std::size_t>(atom)];
    for (int t = 0; t < nThreads; t++) {
      head[static_cast<std::size_t>(t)][static_cast<std::size_t>(atom)] = cursor;
      cursor += counts[static_cast<std::size_t>(t)][static_cast<std::size_t>(atom)];
    }
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel
#endif
  {
    const int tid = omp_get_thread_num();
    auto &cursor = head[static_cast<std::size_t>(tid)];
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int iring = 0; iring < nRings; iring++) {
      for (const int atom : rings[static_cast<std::size_t>(iring)]) {
        if (atom >= 0 && atom < numAtoms) {
          flat[static_cast<std::size_t>(cursor[static_cast<std::size_t>(atom)]++)] = iring;
        }
      }
    }
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int atom = 0; atom < numAtoms; atom++) {
    int n = 0;
    for (int t = 0; t < nThreads; t++) {
      n += counts[static_cast<std::size_t>(t)][static_cast<std::size_t>(atom)];
    }
    auto &row = index.ringsContainingAtom[static_cast<std::size_t>(atom)];
    const int begin = base[static_cast<std::size_t>(atom)];
    row.assign(flat.begin() + begin, flat.begin() + begin + n);
  }
  index.ringOff.assign(static_cast<std::size_t>(numAtoms) + 1, 0);
  for (int atom = 0; atom < numAtoms; atom++) {
    index.ringOff[static_cast<std::size_t>(atom)] = base[static_cast<std::size_t>(atom)];
  }
  index.ringOff[static_cast<std::size_t>(numAtoms)] = total;
  index.ringCol = std::move(flat);
  return index;
}

ring::RingSearchIndex ring::buildRingSearchIndex(
    const std::vector<std::vector<int>> &rings,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  auto index =
      ring::buildRingSearchIndex(rings, static_cast<int>(yCloud.pts.size()));
  // The claim key is the ring's atoms and its unwrapped shift, not the
  // absolute position. positionRank stays available for callers that want
  // the coordinate order; it is not the cage order.
  const auto box = gen::makeFracBox(yCloud);
  index.ringShift.assign(rings.size(), {0, 0, 0});
  if (box.ok) {
    const int nAtoms = static_cast<int>(yCloud.pts.size());
    const int nRings = static_cast<int>(rings.size());
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int r = 0; r < nRings; r++) {
      const auto &ring = rings[static_cast<std::size_t>(r)];
      std::array<int, 3> wind{0, 0, 0};
      bool inFrame = true;
      for (size_t i = 0; i < ring.size() && inFrame; i++) {
        const int a = ring[i];
        const int b = ring[(i + 1) % ring.size()];
        if (a < 0 || b < 0 || a >= nAtoms || b >= nAtoms) {
          inFrame = false;
          break;
        }
        const auto &pa = yCloud.pts[static_cast<size_t>(a)];
        const auto &pb = yCloud.pts[static_cast<size_t>(b)];
        const auto mic = gen::fracDelta(box, pb.x, pb.y, pb.z, pa.x, pa.y, pa.z);
        const double raw[3] = {pb.x - pa.x, pb.y - pa.y, pb.z - pa.z};
        const double len[3] = {box.lx, box.ly, box.lz};
        for (int k = 0; k < 3; k++) {
          if (!(len[k] > 0.0)) {
            continue;
          }
          wind[static_cast<size_t>(k)] +=
              static_cast<int>(std::lround((raw[k] - mic[static_cast<size_t>(k)]) / len[k]));
        }
      }
      if (inFrame) {
        index.ringShift[static_cast<std::size_t>(r)] = wind;
      }
    }
  }
  return index;
}

/**
 * @details Throws std::invalid_argument when two atoms share a wrapped
 *  position: no position order separates them.
 */
std::vector<int> ring::positionRank(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  const int n = static_cast<int>(yCloud.pts.size());
  std::vector<std::array<double, 3>> wrapped(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    wrapped[static_cast<std::size_t>(i)] = {p.x, p.y, p.z};
    for (size_t d = 0; d < 3; d++) {
      if (d < yCloud.box.size() && yCloud.box[d] > 0) {
        const double lo = d < yCloud.boxLow.size() ? yCloud.boxLow[d] : 0.0;
        wrapped[static_cast<std::size_t>(i)][d] -=
            yCloud.box[d] * std::floor((wrapped[static_cast<std::size_t>(i)][d] - lo) / yCloud.box[d]);
      }
    }
  }
  std::vector<int> order(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    order[static_cast<std::size_t>(i)] = i;
  }
  std::sort(order.begin(), order.end(), [&wrapped](int a, int b) {
    return wrapped[static_cast<std::size_t>(a)] != wrapped[static_cast<std::size_t>(b)] ? wrapped[static_cast<std::size_t>(a)] < wrapped[static_cast<std::size_t>(b)] : a < b;
  });
  for (int r = 1; r < n; r++) {
    if (wrapped[static_cast<std::size_t>(order[static_cast<std::size_t>(r)])] ==
        wrapped[static_cast<std::size_t>(order[static_cast<std::size_t>(r - 1)])]) {
      throw std::invalid_argument("atoms " + std::to_string(order[static_cast<std::size_t>(r - 1)]) +
                                  " and " + std::to_string(order[static_cast<std::size_t>(r)]) +
                                  " share a wrapped position");
    }
  }
  std::vector<int> rank(static_cast<std::size_t>(n));
  for (int r = 0; r < n; r++) {
    rank[static_cast<std::size_t>(order[static_cast<std::size_t>(r)])] = r;
  }
  return rank;
}

namespace {

// Lattice image of the minimum-image step from a to b. The same wrap as
// gen::fracDelta: orthorhombic axes use floor(d/L + 1/2), a tilt box wraps
// the fractional coordinates with round. Summed around a cycle the images
// are the winding. A contractible ring sums to zero.
std::array<int, 3> bondImage(const gen::FracBox &box, double ax, double ay,
                             double az, double bx, double by, double bz) {
  if (!box.triclinic) {
    auto shift = [](double d, double len) {
      return static_cast<int>(std::floor(d / len + 0.5));
    };
    return {shift(bx - ax, box.lx), shift(by - ay, box.ly),
            shift(bz - az, box.lz)};
  }
  const double szA = (az - box.oz) / box.lz;
  const double syA = (ay - box.oy - box.yz * szA) / box.ly;
  const double sxA = (ax - box.ox - box.xy * syA - box.xz * szA) / box.lx;
  const double szB = (bz - box.oz) / box.lz;
  const double syB = (by - box.oy - box.yz * szB) / box.ly;
  const double sxB = (bx - box.ox - box.xy * syB - box.xz * szB) / box.lx;
  return {static_cast<int>(std::round(sxB - sxA)),
          static_cast<int>(std::round(syB - syA)),
          static_cast<int>(std::round(szB - szA))};
}

bool ringContracts(const gen::FracBox &box,
                   const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                   const std::vector<int> &ring, int nAtoms) {
  const std::size_t n = ring.size();
  std::array<int, 3> sum = {0, 0, 0};
  for (std::size_t i = 0; i < n; i++) {
    const int a = ring[i];
    const int b = ring[(i + 1) % n];
    if (a < 0 || b < 0 || a >= nAtoms || b >= nAtoms) {
      return false;
    }
    const auto &pa = yCloud.pts[static_cast<std::size_t>(a)];
    const auto &pb = yCloud.pts[static_cast<std::size_t>(b)];
    const auto image = bondImage(box, pa.x, pa.y, pa.z, pb.x, pb.y, pb.z);
    sum[0] += image[0];
    sum[1] += image[1];
    sum[2] += image[2];
  }
  return sum[0] == 0 && sum[1] == 0 && sum[2] == 0;
}

} // namespace

/**
 * @details Each bond contributes the integer lattice image of its minimum-image
 *  step. Those images sum to zero exactly when the ring does not wind. A
 *  floating sum against half the cell confuses a short lattice vector with
 *  noise once the cell is thin; the integers do not.
 */
std::vector<std::vector<int>> ring::contractibleRings(
    const std::vector<std::vector<int>> &rings,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  const auto box = gen::makeFracBox(yCloud);
  if (!box.ok) {
    return rings;
  }
  const int nAtoms = static_cast<int>(yCloud.pts.size());
  const int nR = static_cast<int>(rings.size());
  std::vector<char> keep(rings.size(), 0);
  // Each ring reads only its own atoms. The kept list stays in input order.
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
  const bool scale = nR >= 1024 && omp_get_max_threads() > 1 && !omp_in_parallel();
#else
  const bool scale = false;
#endif
  if (scale) {
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
#pragma omp parallel for schedule(static)
    for (int i = 0; i < nR; ++i) {
      keep[static_cast<std::size_t>(i)] =
          ringContracts(box, yCloud, rings[static_cast<std::size_t>(i)], nAtoms) ? 1 : 0;
    }
#endif
  } else {
    for (int i = 0; i < nR; ++i) {
      keep[static_cast<std::size_t>(i)] =
          ringContracts(box, yCloud, rings[static_cast<std::size_t>(i)], nAtoms) ? 1 : 0;
    }
  }
  std::vector<std::vector<int>> kept;
  kept.reserve(rings.size());
  for (int i = 0; i < nR; ++i) {
    if (keep[static_cast<std::size_t>(i)] != 0) {
      kept.push_back(rings[static_cast<std::size_t>(i)]);
    }
  }
  return kept;
}

std::vector<std::vector<int>> ring::cageRings(
    const std::vector<std::vector<int>> &rings,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double cutoff) {
  const auto box = gen::makeFracBox(yCloud);
  // Smith: below half the narrowest face separation a pair has one image
  // inside the cutoff. At or beyond that length the index graph has collapsed
  // two images of one neighbour into one edge, and a cage built from it is
  // not a cage of the structure.
  if (box.ok && !(cutoff + 1e-12 < box.halfMin)) {
    throw std::invalid_argument(
        "cage search refused: the cell is thinner than twice the neighbour "
        "cutoff, so one pair can meet through more than one periodic image");
  }
  return ring::contractibleRings(rings, yCloud);
}

/**
 * @details Convenience overload that builds the index itself. Prefer the
 *  overload taking a prebuilt index when several searches run over one ring
 *  network.
 */
std::vector<int> ring::findDDC(const std::vector<std::vector<int>> &rings,
                               std::vector<ring::strucType> &ringType,
                               const std::vector<int> &listHC,
                               std::vector<cage::Cage> &cageList) {
  int maxAtom = 0;
  for (const auto &r : rings) {
    for (const int atom : r) {
      maxAtom = std::max(maxAtom, atom);
    }
  }
  const auto index = ring::buildRingSearchIndex(rings, maxAtom + 1);
  return ring::findDDC(rings, ringType, listHC, cageList, index);
}

/**
 * @details Determines which hexagonal rings are DDC rings. This function
 * returns a vector which contains the ring IDs of all the rings which are DDC
 * rings. The ring IDs correspond to the index of the rings inside the vector of
 * vector rings, starting from 0. Candidate equatorial rings are visited in the
 * order of canonicalRingOrder (sorted atom indices, then the unwrapped shift,
 * then the cyclic reading). HC rings are never equatorial. Every other ring
 * that passes ring::equatorialConditions is a candidate: every atom of the
 * ring lies in at least three other rings, every triplet of the ring lies in
 * at least one other ring, and six distinct peripheral rings satisfy the
 * alternating-triplet conditions. Of the passing choices of peripherals, the
 * one whose ring keys are least is taken. Two candidates conflict when each
 * equator's ring lies in the other's peripheral set, so at most one of them
 * can be claimed. The reported cages are a maximum non-conflicting set. The
 * equatorial ring and its peripherals are then set to the DDC type inside
 * ringType (bothBasal or bothPrismatic for a peripheral that is already an
 * HC ring), and the cage is appended to cageList.
 *
 * @param[in] rings Vector of vectors containing 6-membered primitive rings
 *  (wherein each ring contains the atom indices of the particles which
 * constitute the ring).
 * @param[in] ringType Vector containing a ring::strucType (structure type) for
 *  each ring.
 * @param[in] listHC Ring IDs of the HC rings, which are never equatorial.
 * @param[in] cageList Vector in which every cage is saved.
 * @param[in] index Atom-to-rings index of rings. ringShift, when set, is the
 *  unwrapped winding of each ring and is part of the ring key.
 * @return A vector of all the ring indices which constitute DDCs.
 */
namespace {
struct DdcSplitTimes {
  double order_ms = 0;
  double equator_ms = 0;
  double pack_ms = 0;
  double conf_ms = 0;
  double split_ms = 0;
  double clique_ms = 0;
  double match_ms = 0;
  double search_ms = 0;
  double write_ms = 0;
};
DdcSplitTimes gDdcSplit;
// Cycle keys built by the last findDDC order. Ideal ice distinguishes rings
// by their atoms, so this stays zero when the comparator can stop early.
std::atomic<std::uint64_t> gDdcCycleEvals{0};
} // namespace

namespace {
std::mutex gEquatorOverflowMu;
std::unordered_map<const ring::EquatorChoice *, std::vector<int>> gEquatorOverflow;
} // namespace

int ring::equatorChoiceOverflow(const ring::EquatorChoice *self, int at) {
  std::lock_guard<std::mutex> lock(gEquatorOverflowMu);
  const auto it = gEquatorOverflow.find(self);
  if (it == gEquatorOverflow.end() || at < 0 ||
      at >= static_cast<int>(it->second.size())) {
    return 0;
  }
  return it->second[static_cast<std::size_t>(at)];
}

void ring::equatorChoiceOverflowPut(ring::EquatorChoice *self, int at, int value) {
  std::lock_guard<std::mutex> lock(gEquatorOverflowMu);
  auto &past = gEquatorOverflow[self];
  if (at >= static_cast<int>(past.size())) {
    past.resize(static_cast<std::size_t>(at) + 1);
  }
  past[static_cast<std::size_t>(at)] = value;
}

void ring::releaseEquatorOverflow() {
  std::lock_guard<std::mutex> lock(gEquatorOverflowMu);
  gEquatorOverflow.clear();
}

extern "C" __attribute__((visibility("default"), used)) std::uint64_t seams_ddc_cycle_evals() {
  return gDdcCycleEvals.load(std::memory_order_relaxed);
}

// Cubic mW has parity product 1. Frame 148 has parity product at most 16.
// The triple loop costs the product, so it stays in use through the square
// of that frame product. The guided team starts at the same count: a shorter
// ring list is no larger than one such loop.
constexpr int kSmallParityProduct = 16;
constexpr int kTripleLoopProduct = kSmallParityProduct * kSmallParityProduct;

std::vector<int> canonicalRingOrderByAtoms(
    const std::vector<std::vector<int>> &rings, const ring::RingSearchIndex &index) {
  const int n = static_cast<int>(rings.size());
  // The rotation key is needed only when two rings share atom ids and a
  // shift. Ideal ice never ties, so that key stays off this record.
  struct Key {
    std::array<int, 6> id{};
    int len = 0;
    std::array<int, 3> shift{};
  };
  std::vector<Key> keys(static_cast<size_t>(n));
  std::vector<CycleKey> cycles;
  std::vector<char> cycleReady;
  std::mutex cycleMu;
  const int nShift = static_cast<int>(index.ringShift.size());
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
  for (int i = 0; i < n; i++) {
    const auto &ring = rings[static_cast<size_t>(i)];
    Key &k = keys[static_cast<size_t>(i)];
    k.len = static_cast<int>(std::min<size_t>(ring.size(), 6));
    for (int j = 0; j < k.len; j++) {
      k.id[static_cast<size_t>(j)] = ring[static_cast<size_t>(j)];
    }
    std::sort(k.id.begin(), k.id.begin() + k.len);
    if (i < nShift) {
      k.shift = index.ringShift[static_cast<size_t>(i)];
    }
  }
  int maxAtom = 0;
  for (int i = 0; i < n; i++) {
    if (keys[static_cast<size_t>(i)].len > 0) {
      maxAtom = std::max(maxAtom, keys[static_cast<size_t>(i)].id[0]);
    }
  }
  // Histogram of the first atom, then an exclusive scan gives each bucket start.
  std::vector<int> bucket(static_cast<size_t>(maxAtom + 3), 0);
  for (int i = 0; i < n; i++) {
    const int b = keys[static_cast<size_t>(i)].len > 0 ? keys[static_cast<size_t>(i)].id[0] + 1 : 0;
    bucket[static_cast<size_t>(b)]++;
  }
  seams::exclusiveScan(bucket);
  std::vector<int> order(static_cast<size_t>(n));
  std::vector<int> place = bucket;
  for (int i = 0; i < n; i++) {
    const int b = keys[static_cast<size_t>(i)].len > 0 ? keys[static_cast<size_t>(i)].id[0] + 1 : 0;
    order[static_cast<size_t>(place[static_cast<size_t>(b)]++)] = i;
  }
  auto before = [&](int a, int b) {
    const Key &ka = keys[static_cast<size_t>(a)];
    const Key &kb = keys[static_cast<size_t>(b)];
    const int m = ka.len < kb.len ? ka.len : kb.len;
    for (int i = 0; i < m; i++) {
      if (ka.id[static_cast<size_t>(i)] != kb.id[static_cast<size_t>(i)]) {
        return ka.id[static_cast<size_t>(i)] < kb.id[static_cast<size_t>(i)];
      }
    }
    if (ka.len != kb.len) {
      return ka.len < kb.len;
    }
    if (ka.shift != kb.shift) {
      return ka.shift < kb.shift;
    }
    // The atom ids already separate almost every pair. The rotation key is
    // built only for a tie, which is the comparison the eager pass used to
    // skip.
    auto ensureCycle = [&](int id) {
      std::lock_guard<std::mutex> lock(cycleMu);
      if (cycles.empty()) {
        cycles.resize(static_cast<std::size_t>(n));
        cycleReady.assign(static_cast<std::size_t>(n), 0);
      }
      if (cycleReady[static_cast<std::size_t>(id)]) {
        return;
      }
      const auto &ring = rings[static_cast<std::size_t>(id)];
      cycles[static_cast<std::size_t>(id)] =
          bestCycle(std::span<const int>(ring.data(), ring.size()));
      cycleReady[static_cast<std::size_t>(id)] = 1;
      gDdcCycleEvals.fetch_add(1, std::memory_order_relaxed);
    };
    ensureCycle(a);
    ensureCycle(b);
    const CycleKey &ca = cycles[static_cast<std::size_t>(a)];
    const CycleKey &cb = cycles[static_cast<std::size_t>(b)];
    if (ca < cb) {
      return true;
    }
    if (cb < ca) {
      return false;
    }
    return a < b;
  };
  const int nBuckets = maxAtom + 2;
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
  for (int b = 0; b < nBuckets; b++) {
    const int lo = bucket[static_cast<size_t>(b)];
    const int hi = b + 1 < static_cast<int>(bucket.size()) ? bucket[static_cast<size_t>(b + 1)]
                                                           : n;
    if (hi - lo > 1) {
      std::sort(order.begin() + lo, order.begin() + hi, before);
    }
  }
  return order;
}

std::vector<int> ring::findDDC(const std::vector<std::vector<int>> &rings,
                               std::vector<ring::strucType> &ringType,
                               const std::vector<int> &listHC,
                               std::vector<cage::Cage> &cageList,
                               const ring::RingSearchIndex &index,
                               ring::PackingCertificate *certificate) {
  std::vector<int> listDDC;
  const int totalRingNum = static_cast<int>(rings.size());
  std::vector<int> peripheralRings;
  std::vector<int> DDCRings; // Equatorial ring first, then the peripherals
  // true if the ring is prismatic, basal or peripheral
  std::vector<char> notEquatorial(static_cast<size_t>(totalRingNum), 0);
  for (const int currentRingIndex : listHC) {
    notEquatorial[static_cast<std::size_t>(currentRingIndex)] = true;
  }

  // The peripheral tie-break and the emitted order use the ring key
  // (sorted atom ids, unwrapped shift). Which equators are kept is a
  // maximum non-conflicting set, so the count does not follow a claim order.
  gDdcCycleEvals.store(0, std::memory_order_relaxed);
  const auto td0 = std::chrono::steady_clock::now();
  const std::vector<int> visitOrder = canonicalRingOrderByAtoms(rings, index);
  const auto td1 = std::chrono::steady_clock::now();
  std::vector<int> rank(static_cast<size_t>(totalRingNum));
  for (int r = 0; r < totalRingNum; r++) {
    rank[static_cast<size_t>(visitOrder[static_cast<size_t>(r)])] = r;
  }
  std::vector<ring::EquatorChoice> chosenPeriph(static_cast<size_t>(totalRingNum));
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel for schedule(static, 64) if (totalRingNum >= 256 && omp_get_max_threads() > 1)
#endif
  for (int iring = 0; iring < totalRingNum; iring++) {
    ring::EquatorChoice &local = chosenPeriph[static_cast<size_t>(iring)];
    if (notEquatorial[static_cast<size_t>(iring)]) {
      continue;
    }
    if (!ring::equatorialConditions(rings, iring, index, nullptr, nullptr, &rank, &local)) {
      continue;
    }
  }
  const auto td2 = std::chrono::steady_clock::now();
  const std::vector<int> equators =
      ring::maximumEquators(chosenPeriph, visitOrder, certificate);
  const auto td3 = std::chrono::steady_clock::now();
  const ring::PackParts parts = ring::lastPackParts();
  gDdcSplit.order_ms = std::chrono::duration<double, std::milli>(td1 - td0).count();
  gDdcSplit.equator_ms = std::chrono::duration<double, std::milli>(td2 - td1).count();
  gDdcSplit.pack_ms = std::chrono::duration<double, std::milli>(td3 - td2).count();
  gDdcSplit.conf_ms = parts.conf_ms;
  gDdcSplit.split_ms = parts.split_ms;
  gDdcSplit.clique_ms = parts.clique_ms;
  gDdcSplit.match_ms = parts.match_ms;
  gDdcSplit.search_ms = parts.search_ms;
  if (std::getenv("DSEAMS_CAGE_TIME") != nullptr) {
    std::fprintf(stderr, "ddc order %.3f eq %.3f pack %.3f kept %zu\n", gDdcSplit.order_ms,
                 gDdcSplit.equator_ms, gDdcSplit.pack_ms, equators.size());
  }
  for (const int iring : equators) {
    const ring::EquatorChoice &choice = chosenPeriph[static_cast<size_t>(iring)];
    peripheralRings.resize(static_cast<std::size_t>(choice.n));
    for (int k = 0; k < choice.n; k++) {
      peripheralRings[static_cast<std::size_t>(k)] = choice.idAt(k);
    }
    if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::unclassified) {
      ringType[static_cast<std::size_t>(iring)] = ring::strucType::DDC;
      listDDC.push_back(iring);
    }
    for (const int jring : peripheralRings) {
      if (ringType[static_cast<std::size_t>(jring)] == ring::strucType::unclassified) {
        ringType[static_cast<std::size_t>(jring)] = ring::strucType::DDC;
        listDDC.push_back(jring);
      } else if (ringType[static_cast<std::size_t>(jring)] == ring::strucType::HCbasal) {
        ringType[static_cast<std::size_t>(jring)] = ring::strucType::bothBasal;
        listDDC.push_back(jring);
      } else if (ringType[static_cast<std::size_t>(jring)] == ring::strucType::HCprismatic) {
        ringType[static_cast<std::size_t>(jring)] = ring::strucType::bothPrismatic;
        listDDC.push_back(jring);
      }
    }
    DDCRings.clear();
    DDCRings.push_back(iring);
    DDCRings.insert(DDCRings.end(), peripheralRings.begin(),
                    peripheralRings.end());
    cageList.push_back({cage::cageType::DoubleDiaC, DDCRings});
  }
  ring::releaseEquatorOverflow();
  std::vector<ring::EquatorChoice>().swap(chosenPeriph);
  gDdcSplit.write_ms =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - td3).count();

  return listDDC;
}

ring::DdcSplit ring::lastDdcSplit() {
  return {gDdcSplit.order_ms,     gDdcSplit.equator_ms, gDdcSplit.pack_ms,
          gDdcSplit.conf_ms,      gDdcSplit.split_ms,   gDdcSplit.clique_ms,
          gDdcSplit.match_ms,     gDdcSplit.search_ms,  gDdcSplit.write_ms};
}

/**
 * @details The triplet T_k of ring iring is (m_k, m_{k+1}, m_{k+2}). A
 *  rotation or reversal of the ring maps the triplets onto each other and
 *  the even ones onto the even or the odd ones, so the conditions are
 *  tested for every starting atom and direction at once: a chain of three
 *  alternate rings records which of its rings can sit in the middle (sharing
 *  three atoms with both others), and an even and an odd chain combine when
 *  two of their middles sit on adjacent triplets. The choice of P_k is
 *  existential as well, with the six rings distinct.
 */
bool ring::equatorialConditions(const std::vector<std::vector<int>> &rings,
                                int iring, const ring::RingSearchIndex &index,
                                std::vector<int> *peripheralUnion,
                                std::vector<int> *chosen,
                                const std::vector<int> *rank,
                                ring::EquatorChoice *raw) {
  if (peripheralUnion != nullptr) {
    peripheralUnion->clear();
  }
  if (chosen != nullptr) {
    chosen->clear();
  }
  if (raw != nullptr) {
    raw->n = 0;
  }
  const std::vector<int> &R = rings[static_cast<std::size_t>(iring)];
  if (R.size() != 6) {
    return false;
  }
  for (const int atom : R) {
    const auto row = atomRow(index, atom);
    int others = 0;
    for (const int r : row) {
      others += static_cast<int>(r != iring);
    }
    if (others < 3) {
      return false;
    }
  }
  std::array<IntBuf, 6> T;
  for (int k = 0; k < 6; k++) {
    ringsThrough(index, R[static_cast<std::size_t>(k)], R[static_cast<std::size_t>((k + 1) % 6)], R[static_cast<std::size_t>((k + 2) % 6)], iring, -1, T[static_cast<std::size_t>(k)]);
    if (T[static_cast<std::size_t>(k)].empty()) {
      return false;
    }
  }
  // One ring on every triplet is the ice network. That assignment is the
  // only candidate, so the chain search is a single 14-molecule check.
  bool oneEach = true;
  for (int k = 0; k < 6; k++) {
    if (T[static_cast<std::size_t>(k)].size() != 1) {
      oneEach = false;
      break;
    }
  }
  if (oneEach) {
    int r[6];
    bool distinct = true;
    for (int k = 0; k < 6; k++) {
      r[k] = T[static_cast<std::size_t>(k)][0];
      for (int j = 0; j < k; j++) {
        distinct = distinct && r[j] != r[k];
      }
    }
    auto share6 = [](const std::vector<int> &a, const std::vector<int> &b) {
      const int *x = a.data();
      const int *y = b.data();
      int n = 0;
      for (int i = 0; i < 6; i++) {
        const int v = x[i];
        n += (v == y[0]) + (v == y[1]) + (v == y[2]) + (v == y[3]) + (v == y[4]) + (v == y[5]);
      }
      return n;
    };
    auto chainBits = [&](int a, int b, int c) -> unsigned {
      const auto &A = rings[static_cast<size_t>(a)];
      const auto &B = rings[static_cast<size_t>(b)];
      const auto &C = rings[static_cast<size_t>(c)];
      if (!commonAtom(A, B, C)) {
        return 0U;
      }
      const bool ab = share6(A, B) >= 3;
      const bool bc = share6(B, C) >= 3;
      const bool ca = share6(C, A) >= 3;
      return (ca && ab ? 1U : 0U) | (ab && bc ? 2U : 0U) | (bc && ca ? 4U : 0U);
    };
    const unsigned evenBits = chainBits(r[0], r[2], r[4]);
    const unsigned oddBits = chainBits(r[1], r[3], r[5]);
    unsigned partners = 0U;
    for (unsigned t = 0; t < 3; t++) {
      if ((evenBits >> t) & 1U) {
        partners |= (1U << t) | (1U << ((t + 2) % 3));
      }
    }
    if (!(distinct && evenBits != 0U && (partners & oddBits) != 0U)) {
      return false;
    }
    bool hex = true;
    for (int k = 0; k < 6; k++) {
      hex = hex && rings[static_cast<size_t>(r[k])].size() == 6;
    }
    if (!hex) {
      return false;
    }
    int atoms[42];
    int nAtomsHere = 0;
    auto take = [&](const std::vector<int> &cycle) {
      for (const int a : cycle) {
        atoms[nAtomsHere++] = a;
      }
    };
    take(R);
    for (int k = 0; k < 6; k++) {
      take(rings[static_cast<size_t>(r[k])]);
    }
    // Only the number of distinct molecules matters. A stamp over the frame
    // counts them without sorting the forty-two ids. A wide or negative id
    // falls back to the same count sort+unique would return.
    const int nAtoms = static_cast<int>(index.ringsContainingAtom.size());
    int nUnique = 0;
    bool counted = nAtoms > 0;
    if (counted) {
      thread_local std::vector<int> stampStore;
      thread_local int era = 0;
      std::vector<int> *stamp = &stampStore;
      asm volatile("" : "+r"(stamp));
      if (static_cast<int>(stamp->size()) < nAtoms) {
        stamp->assign(static_cast<std::size_t>(nAtoms), 0);
      }
      if (++era == 0) {
        std::fill(stamp->begin(), stamp->end(), 0);
        era = 1;
      }
      int *seen = stamp->data();
      for (int i = 0; i < nAtomsHere; i++) {
        const int atom = atoms[i];
        if (atom < 0 || atom >= nAtoms) {
          counted = false;
          break;
        }
        if (seen[atom] != era) {
          seen[atom] = era;
          nUnique++;
        }
      }
    }
    if (!counted) {
      nUnique = countDistinct(atoms, nAtomsHere);
    }
    bool folded = nUnique != 14;
    for (int k = 0; k < 6 && !folded; k++) {
      folded = share6(rings[static_cast<size_t>(r[k])], R) != 3;
    }
    if (folded) {
      return false;
    }
    if (peripheralUnion != nullptr) {
      peripheralUnion->assign(r, r + 6);
      std::sort(peripheralUnion->begin(), peripheralUnion->end());
    }
    if (rank != nullptr && (chosen != nullptr || raw != nullptr)) {
      int ord[6] = {r[0], r[1], r[2], r[3], r[4], r[5]};
      std::sort(ord, ord + 6, [rank](int a, int b) {
        return (*rank)[static_cast<size_t>(a)] < (*rank)[static_cast<size_t>(b)];
      });
      if (chosen != nullptr) {
        chosen->assign(ord, ord + 6);
      }
      if (raw != nullptr) {
        for (int k = 0; k < 6; k++) {
          raw->ring[k] = ord[k];
        }
        raw->n = 6;
      }
    }
    return true;
  }
  const std::int64_t prodEven = static_cast<std::int64_t>(T[0].size()) *
                                static_cast<std::int64_t>(T[2].size()) *
                                static_cast<std::int64_t>(T[4].size());
  const std::int64_t prodOdd = static_cast<std::int64_t>(T[1].size()) *
                               static_cast<std::int64_t>(T[3].size()) *
                               static_cast<std::int64_t>(T[5].size());
  // Products at or below kTripleLoopProduct stay on the triple loop. Above
  // that the share graph enumerates the same chains, and a pair is kept only
  // when the seven rings span 14 molecules.
  if (prodEven > kTripleLoopProduct || prodOdd > kTripleLoopProduct) {
    return equatorialLarge(rings, R, T,
                           static_cast<int>(index.ringsContainingAtom.size()),
                           peripheralUnion, chosen, rank, raw);
  }
  struct Chain {
    std::array<int, 3> r;
    unsigned middles; // bit t: the ring on triplet parity + 2t is a middle
  };
  auto chainsOf = [&](int parity) {
    std::vector<Chain> out;
    const IntBuf &Ta = T[static_cast<size_t>(parity)];
    const IntBuf &Tb = T[static_cast<size_t>(parity + 2)];
    const IntBuf &Tc = T[static_cast<size_t>(parity + 4)];
    const int na = Ta.size();
    const int nb = Tb.size();
    const int nc = Tc.size();
    // An empty leg has no chain. The share of a pair does not depend on the
    // third ring, and the product of the three lengths is what this loop
    // used to pay for three share tests. The pair tables are at most the
    // product, because this path is only taken when that product is small.
    if (na == 0 || nb == 0 || nc == 0) {
      return out;
    }
    std::vector<char> shareAB(static_cast<std::size_t>(na) * static_cast<std::size_t>(nb));
    std::vector<char> shareBC(static_cast<std::size_t>(nb) * static_cast<std::size_t>(nc));
    std::vector<char> shareCA(static_cast<std::size_t>(nc) * static_cast<std::size_t>(na));
    for (int ia = 0; ia < na; ia++) {
      const int a = Ta[ia];
      const auto &A = rings[static_cast<std::size_t>(a)];
      for (int ib = 0; ib < nb; ib++) {
        const int b = Tb[ib];
        shareAB[static_cast<std::size_t>(ia) * static_cast<std::size_t>(nb) +
                static_cast<std::size_t>(ib)] =
            static_cast<char>(b != a && sharedAtoms(A, rings[static_cast<std::size_t>(b)]) >= 3);
      }
      for (int ic = 0; ic < nc; ic++) {
        const int c = Tc[ic];
        shareCA[static_cast<std::size_t>(ic) * static_cast<std::size_t>(na) +
                static_cast<std::size_t>(ia)] =
            static_cast<char>(c != a && sharedAtoms(rings[static_cast<std::size_t>(c)], A) >= 3);
      }
    }
    for (int ib = 0; ib < nb; ib++) {
      const int b = Tb[ib];
      const auto &B = rings[static_cast<std::size_t>(b)];
      for (int ic = 0; ic < nc; ic++) {
        const int c = Tc[ic];
        shareBC[static_cast<std::size_t>(ib) * static_cast<std::size_t>(nc) +
                static_cast<std::size_t>(ic)] =
            static_cast<char>(c != b && sharedAtoms(B, rings[static_cast<std::size_t>(c)]) >= 3);
      }
    }
    for (int ia = 0; ia < na; ia++) {
      const int a = Ta[ia];
      for (int ib = 0; ib < nb; ib++) {
        const int b = Tb[ib];
        if (b == a) {
          continue;
        }
        const bool ab =
            shareAB[static_cast<std::size_t>(ia) * static_cast<std::size_t>(nb) +
                    static_cast<std::size_t>(ib)] != 0;
        for (int ic = 0; ic < nc; ic++) {
          const int c = Tc[ic];
          if (c == a || c == b) {
            continue;
          }
          const bool bc =
              shareBC[static_cast<std::size_t>(ib) * static_cast<std::size_t>(nc) +
                      static_cast<std::size_t>(ic)] != 0;
          const bool ca =
              shareCA[static_cast<std::size_t>(ic) * static_cast<std::size_t>(na) +
                      static_cast<std::size_t>(ia)] != 0;
          const unsigned middles = (ca && ab ? 1U : 0U) |
                                   (ab && bc ? 2U : 0U) | (bc && ca ? 4U : 0U);
          if (middles == 0U) {
            continue;
          }
          if (!commonAtom(rings[static_cast<std::size_t>(a)], rings[static_cast<std::size_t>(b)], rings[static_cast<std::size_t>(c)])) {
            continue;
          }
          out.push_back({{a, b, c}, middles});
        }
      }
    }
    return out;
  };
  const std::vector<Chain> even = chainsOf(0);
  if (even.empty()) {
    return false;
  }
  const std::vector<Chain> odd = chainsOf(1);
  if (odd.empty()) {
    return false;
  }
  bool pass = false;
  std::vector<int> best;
  std::vector<int> candidate(6);
  int bestIds[6] = {};
  bool haveBest = false;
  for (const Chain &e : even) {
    // Even middle 2t sits next to the odd middles 2t+1 and 2t-1
    unsigned partners = 0U;
    for (unsigned t = 0; t < 3; t++) {
      if ((e.middles >> t) & 1U) {
        partners |= (1U << t) | (1U << ((t + 2) % 3));
      }
    }
    for (const Chain &o : odd) {
      if ((partners & o.middles) == 0U) {
        continue;
      }
      bool distinct = true;
      for (const int x : e.r) {
        for (const int y : o.r) {
          distinct = distinct && x != y;
        }
      }
      if (!distinct) {
        continue;
      }
      // Seven hexagons of a DDC span 14 molecules: the equatorial ring,
      // three first neighbours and an apex on each staggered triplet.
      // A peripheral that shares a fourth equatorial atom is not that cage.
      int atoms[42];
      int nAtomsHere = 0;
      bool hexagons = true;
      auto take = [&](const std::vector<int> &cycle) {
        if (cycle.size() != 6 || nAtomsHere + 6 > 42) {
          hexagons = false;
          return;
        }
        for (const int a : cycle) {
          atoms[nAtomsHere++] = a;
        }
      };
      take(R);
      for (const int ri : e.r) {
        take(rings[static_cast<std::size_t>(ri)]);
      }
      for (const int ri : o.r) {
        take(rings[static_cast<std::size_t>(ri)]);
      }
      if (!hexagons) {
        continue;
      }
      const int nUnique = countDistinct(atoms, nAtomsHere);
      bool folded = false;
      for (const int ri : e.r) {
        folded = folded || sharedAtoms(rings[static_cast<std::size_t>(ri)], R) != 3;
      }
      for (const int ri : o.r) {
        folded = folded || sharedAtoms(rings[static_cast<std::size_t>(ri)], R) != 3;
      }
      if (folded || nUnique != 14) {
        continue;
      }
      pass = true;
      if (peripheralUnion != nullptr) {
        peripheralUnion->insert(peripheralUnion->end(), e.r.begin(), e.r.end());
        peripheralUnion->insert(peripheralUnion->end(), o.r.begin(), o.r.end());
      }
      if (rank != nullptr && (chosen != nullptr || raw != nullptr)) {
        for (int t = 0; t < 3; t++) {
          candidate[static_cast<std::size_t>(t)] = (*rank)[static_cast<size_t>(e.r[static_cast<size_t>(t)])];
          candidate[static_cast<std::size_t>(3 + t)] = (*rank)[static_cast<size_t>(o.r[static_cast<size_t>(t)])];
        }
        std::sort(candidate.begin(), candidate.end());
        if (!haveBest || candidate < best) {
          best = candidate;
          haveBest = true;
          for (int t = 0; t < 3; t++) {
            bestIds[t] = e.r[static_cast<size_t>(t)];
            bestIds[3 + t] = o.r[static_cast<size_t>(t)];
          }
          if (chosen != nullptr) {
            chosen->assign(bestIds, bestIds + 6);
          }
        }
      } else if (peripheralUnion == nullptr && raw == nullptr) {
        return true;
      }
    }
  }
  if (!pass) {
    return false;
  }
  if (peripheralUnion != nullptr) {
    std::sort(peripheralUnion->begin(), peripheralUnion->end());
    peripheralUnion->erase(
        std::unique(peripheralUnion->begin(), peripheralUnion->end()),
        peripheralUnion->end());
  }
  if (haveBest && rank != nullptr) {
    std::sort(bestIds, bestIds + 6, [rank](int a, int b) {
      return (*rank)[static_cast<size_t>(a)] < (*rank)[static_cast<size_t>(b)];
    });
    if (chosen != nullptr) {
      chosen->assign(bestIds, bestIds + 6);
    }
    if (raw != nullptr) {
      for (int t = 0; t < 6; t++) {
        raw->ring[t] = bestIds[t];
      }
      raw->n = 6;
    }
  }
  return true;
}

/**
 * @details For a given ring, which is being tested as the equatorial ring,
 * this function tests if each element of the ring \f$ (m_k) \f$ is present
 * in at least three other rings or not. Returns false if this is not true.
 * The ring indices of rings that have the common element are saved inside
 * the periperal ring vector (peripheralRings) as potential peripheral rings,
 * which is then passed as an input to the subsequent functions testing
 * conditions.
 * @param[in] rings Vector of vectors containing the 6-membered primitive
 *  rings.
 * @param[in] peripheralRings Vector containing the indices of rings which are
 *  potential peripheral rings.
 * @param[in] iring Index of the ring being tested as equatorial.
 * @return A bool; true if ring being tested as equatorial, iring, satisfies
 *  the above condition for being an equatorial ring, and false otherwise.
 */
bool ring::conditionOneDDC(const std::vector<std::vector<int>> &rings,
                           std::vector<int> &peripheralRings, int iring) {
  int maxAtom = 0;
  for (const auto &r : rings) {
    for (const int atom : r) {
      maxAtom = std::max(maxAtom, atom);
    }
  }
  const auto index = ring::buildRingSearchIndex(rings, maxAtom + 1);
  return ring::conditionOneDDC(rings, peripheralRings, iring, index);
}

/**
 * @details As the overload above, but answering "which rings contain this
 *  atom" from a prebuilt index rather than by scanning the whole ring network
 *  once per element. The scan made the search quadratic in the ring count for
 *  no gain: the answer is fixed for a given ring network.
 * @param[in] rings Vector of vectors containing the 6-membered primitive
 *  rings.
 * @param[in] peripheralRings Vector containing the indices of rings which are
 *  potential peripheral rings.
 * @param[in] iring Index of the ring being tested as equatorial.
 * @param[in] index Inverted atom-to-rings index over @a rings.
 * @return A bool; true if iring satisfies the condition for being an
 *  equatorial ring, and false otherwise.
 */
bool ring::conditionOneDDC(const std::vector<std::vector<int>> &rings,
                           std::vector<int> &peripheralRings, int iring,
                           const ring::RingSearchIndex &index) {
  // Loop through each element of iring for finding matches
  for (int m = 0; m < 6; m++) {
    const int atom = rings[static_cast<std::size_t>(iring)][static_cast<std::size_t>(m)]; // Atom index to be matched with
    int noOfCommonRings = 0;          // init to zero.

    if (atom < 0 ||
        static_cast<size_t>(atom) >= index.ringsContainingAtom.size()) {
      return false;
    }

    // Every ring holding this atom, in ascending ring order
    for (const int jring : index.ringsContainingAtom[static_cast<std::size_t>(atom)]) {
      if (jring == iring) {
        continue;
      } // Skip for iring
      noOfCommonRings++;
      peripheralRings.push_back(jring);
    } // end of loop through all rings holding atom

    // If less than 3 rings have been found for each element, then this is
    // not an equatorial ring
    if (noOfCommonRings < 3) {
      return false;
    } // end of check for common ring number per element in iring
  }   // end of loop through each element of iring

  // iring is an equatorial ring. The duplicate ring IDs inside
  // peripheralRings should be removed
  std::vector<int>::iterator ip; // Iterator to find the last element upto
                                 // which unique elements are present
  // Duplicate IDs must be removed
  int numElements =
      static_cast<int>(peripheralRings.size()); // number of elements in peripheralRings
  sort(peripheralRings.begin(), peripheralRings.end());
  ip = std::unique(peripheralRings.begin(),
                   peripheralRings.begin() + numElements);
  // Resize peripheral rings to remove undefined terms
  peripheralRings.resize(static_cast<std::size_t>(std::distance(peripheralRings.begin(), ip)));

  return true;
}

/**
 * @details Determines which hexagonal rings are HCs. This function
 * returns a vector which contains the ring IDs of all the rings which are HC
 * rings. The ring IDs correspond to the index of the rings inside the vector of
 * vector rings, starting from 0. HC rings can be found using a three-step
 * procedure, in which first two basal rings are found. Prismatic rings are
 * simply rings which share every face made by upper and lower triplets of the
 * basal rings The neighbour list is also required as an input, which is a
 * vector of vectors, containing atom IDs. The first element of the neighbour
 * list is the atomID of the atom for which the other elements are nearest
 * neighbours. The function calls the following functions internally:
 * - ring::hasCommonElements (Step one: Checks to see if basal1 and basal2, two
 * candidate basal rings, have common elements or not. If they don't, then they
 * cannot be basal rings).
 * - ring::basalConditions (Step two and three: One of the elements of basal2
 * must be the nearest neighbour of either the first (index0; \f$ l_1 \f$) or
 * second (index1; \f$ l_2 \f$) element of basal1. If \f$ m_k \f$ is the nearest
 * neighbour of \f$ l_1 \f$, \f$ m_{k+2} \f$ and \f$ m_{k+4} \f$ must be
 * neighbours of \f$ l_3 \f$ and \f$ l_5 \f$(\f$ l_5 \f$ or \f$ l_3 \f$).
 * Similarly modified for \f$ l_2 \f$).
 * @param[in] rings Vector of vectors containing 6-membered primitive rings
 *  (wherein each ring contains the atom indices of the particles which
 * constitute the ring).
 * @param[in] ringType Vector containing a ring::strucType (structure type) for
 *  each ring.
 * @param[in] nList Row-ordered vector of vectors of the neighbour list (by
 *  index).
 * @param[in] cageList Vector in which every cage is saved.
 * @return A vector of all the ring indices which constitute HCs.
 */
std::vector<int> ring::findHC(const std::vector<std::vector<int>> &rings,
                              std::vector<ring::strucType> &ringType,
                              const std::vector<std::vector<int>> &nList,
                              std::vector<cage::Cage> &cageList) {
  const auto index =
      ring::buildRingSearchIndex(rings, static_cast<int>(nList.size()));
  return ring::findHC(rings, ringType, nList, cageList, index);
}

/**
 * @details As the overload above, but drawing candidate partner rings from a
 *  prebuilt index instead of testing every pair.
 *
 *  ring::basalConditions can only succeed when the second ring contains an
 *  atom that neighbours the first or second atom of the first ring, and a pair
 *  it rejects leaves no state behind. Only rings holding one of those
 *  neighbours are therefore worth testing, which is a small set rather than
 *  the whole network. Candidates are visited in the canonical order of
 *  @a index, so the rings are claimed in the same sequence as an exhaustive
 *  search in that order.
 * @param[in] rings Vector of vectors containing the 6-membered primitive
 *  rings.
 * @param[in,out] ringType Vector describing the type of each ring.
 * @param[in] nList Row-ordered neighbour list, by index.
 * @param[in,out] cageList Vector of cages.
 * @param[in] index Inverted atom-to-rings index over @a rings.
 * @return Vector of ring indices which are HC rings.
 */
std::vector<int> ring::findHC(const std::vector<std::vector<int>> &rings,
                              std::vector<ring::strucType> &ringType,
                              const std::vector<std::vector<int>> &nList,
                              std::vector<cage::Cage> &cageList,
                              const ring::RingSearchIndex &index) {
  std::vector<int> listHC;
  int totalRingNum = static_cast<int>(rings.size());
  std::vector<int>
      HCRings; // Indices of rings which constitute a single HC, with the basal
               // rings first, followed by prismatic rings
  std::vector<int> prismaticRings; // Ring indices of prismatic rings
  std::vector<bool>
      isPrismatic; // Flag for checking if the ring is prismatic (true) or not
                   // (false), since the basal rings are checked
  isPrismatic.resize(static_cast<size_t>(totalRingNum)); // Initialized to false

  // Rings are claimed as prismatic as cages are accepted, so visit them in
  // the canonical order and test each unordered pair once, in both
  // directions. The pair test does not depend on earlier claims; the claim
  // does. Several threads can therefore test pairs, and one thread then
  // accepts them in visit order, which is the serial result at every
  // thread count.
  const std::vector<int> visitOrder = canonicalRingOrder(rings, index.atomRank);
  std::vector<int> rank(static_cast<size_t>(totalRingNum));
  for (int r = 0; r < totalRingNum; r++) {
    rank[static_cast<size_t>(visitOrder[static_cast<size_t>(r)])] = r;
  }
  ring::NeighBond bonds;
  bonds.build(nList);

  // basalConditions(iring, jring) and the swapped call both need jring to
  // hold a neighbour of each atom of some alternate triplet of iring. Rings
  // that fail that test cannot be the second basal ring, so they are not
  // gathered. The surviving rings are still tested with basalConditions and
  // claimed in rank order.
  const int nAtoms = static_cast<int>(nList.size());
  struct GatherScratch {
    std::vector<int> onRing;
    std::vector<int> arm1;
    std::vector<int> arm2;
    std::vector<int> seen;
    int tickOn = 0;
    int tick1 = 0;
    int tick2 = 0;
    int tickSeen = 0;
    void prepare(int atoms, int nRings) {
      if (static_cast<int>(onRing.size()) < atoms) {
        onRing.assign(static_cast<size_t>(atoms), 0);
        arm1.assign(static_cast<size_t>(atoms), 0);
        arm2.assign(static_cast<size_t>(atoms), 0);
      }
      if (static_cast<int>(seen.size()) < nRings) {
        seen.assign(static_cast<size_t>(nRings), 0);
      }
    }
    int next(int &tick, std::vector<int> &stamp) {
      if (++tick <= 0) {
        std::fill(stamp.begin(), stamp.end(), 0);
        tick = 1;
      }
      return tick;
    }
  };

  auto gatherCandidates = [&](int iring, std::vector<int> &candidates,
                              GatherScratch &scratch) {
    candidates.clear();
    const std::vector<int> &basal1 = rings[static_cast<size_t>(iring)];
    if (basal1.size() != 6) {
      return;
    }
    scratch.prepare(nAtoms, totalRingNum);
    const int on = scratch.next(scratch.tickOn, scratch.onRing);
    for (const int atom : basal1) {
      if (atom >= 0 && atom < nAtoms) {
        scratch.onRing[static_cast<size_t>(atom)] = on;
      }
    }
    auto neighbourSpan = [&](int atom) {
      int span = 0;
      if (atom < 0 || atom >= nAtoms) {
        return span;
      }
      const auto &row = nList[static_cast<size_t>(atom)];
      for (size_t n = 1; n < row.size(); n++) {
        const int nb = row[n];
        if (nb < 0 || nb >= nAtoms ||
            scratch.onRing[static_cast<size_t>(nb)] == on) {
          continue;
        }
        span += atomRow(index, nb).size();
      }
      return span;
    };
    auto stampArm = [&](int atom, std::vector<int> &stamp, int tick) {
      if (atom < 0 || atom >= nAtoms) {
        return;
      }
      const auto &row = nList[static_cast<size_t>(atom)];
      for (size_t n = 1; n < row.size(); n++) {
        const int nb = row[n];
        if (nb >= 0 && nb < nAtoms && scratch.onRing[static_cast<size_t>(nb)] != on) {
          stamp[static_cast<size_t>(nb)] = tick;
        }
      }
    };
    for (int parity = 0; parity < 2; parity++) {
      int spans[3];
      int arm = 0;
      int best = 0;
      for (int t = 0; t < 3; t++) {
        spans[t] = neighbourSpan(basal1[static_cast<size_t>((parity + 2 * t) % 6)]);
        if (t == 0 || spans[t] < best) {
          best = spans[t];
          arm = t;
        }
      }
      const int tick1 = scratch.next(scratch.tick1, scratch.arm1);
      const int tick2 = scratch.next(scratch.tick2, scratch.arm2);
      stampArm(basal1[static_cast<size_t>((parity + 2 * ((arm + 1) % 3)) % 6)],
               scratch.arm1, tick1);
      stampArm(basal1[static_cast<size_t>((parity + 2 * ((arm + 2) % 3)) % 6)],
               scratch.arm2, tick2);
      const int seenTick = scratch.next(scratch.tickSeen, scratch.seen);
      const int pivot =
          basal1[static_cast<size_t>((parity + 2 * arm) % 6)];
      if (pivot < 0 || pivot >= nAtoms) {
        continue;
      }
      const auto &row = nList[static_cast<size_t>(pivot)];
      for (size_t n = 1; n < row.size(); n++) {
        const int nb = row[n];
        if (nb < 0 || nb >= nAtoms ||
            scratch.onRing[static_cast<size_t>(nb)] == on) {
          continue;
        }
        const auto ringsOf = atomRow(index, nb);
        for (int t = 0; t < ringsOf.size(); ++t) {
          const int candidate = ringsOf[t];
          if (scratch.seen[static_cast<size_t>(candidate)] == seenTick) {
            continue;
          }
          scratch.seen[static_cast<size_t>(candidate)] = seenTick;
          if (rank[static_cast<size_t>(candidate)] <=
              rank[static_cast<size_t>(iring)]) {
            continue;
          }
          bool share = false;
          bool hit1 = false;
          bool hit2 = false;
          for (const int atom : rings[static_cast<size_t>(candidate)]) {
            if (atom < 0 || atom >= nAtoms) {
              continue;
            }
            if (scratch.onRing[static_cast<size_t>(atom)] == on) {
              share = true;
            }
            if (scratch.arm1[static_cast<size_t>(atom)] == tick1) {
              hit1 = true;
            }
            if (scratch.arm2[static_cast<size_t>(atom)] == tick2) {
              hit2 = true;
            }
          }
          if (!share && hit1 && hit2) {
            candidates.push_back(candidate);
          }
        }
      }
    }
    std::sort(candidates.begin(), candidates.end(),
              [&rank](int a, int b) {
                return rank[static_cast<size_t>(a)] < rank[static_cast<size_t>(b)];
              });
    candidates.erase(std::unique(candidates.begin(), candidates.end()),
                     candidates.end());
  };

  auto claimPair = [&](int iring, int jring, bool cond2, bool reverse) {
    if (ringType[static_cast<size_t>(iring)] == ring::strucType::unclassified) {
      ringType[static_cast<size_t>(iring)] = ring::strucType::HCbasal;
      listHC.push_back(iring);
    } else if (ringType[static_cast<size_t>(iring)] == ring::strucType::DDC) {
      ringType[static_cast<size_t>(iring)] = ring::strucType::bothBasal;
      listHC.push_back(iring);
    }
    if (ringType[static_cast<size_t>(jring)] == ring::strucType::unclassified) {
      ringType[static_cast<size_t>(jring)] = ring::strucType::HCbasal;
      listHC.push_back(jring);
    } else if (ringType[static_cast<size_t>(jring)] == ring::strucType::DDC) {
      ringType[static_cast<size_t>(jring)] = ring::strucType::bothBasal;
      listHC.push_back(jring);
    }
    prismaticRings.clear();
    if (cond2) {
      ring::findPrismatic(rings, listHC, ringType, iring, jring, prismaticRings,
                          index);
    }
    if (reverse) {
      ring::findPrismatic(rings, listHC, ringType, jring, iring, prismaticRings,
                          index);
    }
    std::sort(prismaticRings.begin(), prismaticRings.end(),
              [&rank](int a, int b) {
                return rank[static_cast<size_t>(a)] < rank[static_cast<size_t>(b)];
              });
    prismaticRings.erase(
        std::unique(prismaticRings.begin(), prismaticRings.end()),
        prismaticRings.end());
    for (const int kring : prismaticRings) {
      if (ringType[static_cast<size_t>(kring)] == ring::strucType::unclassified) {
        ringType[static_cast<size_t>(kring)] = ring::strucType::HCprismatic;
        listHC.push_back(kring);
      } else if (ringType[static_cast<size_t>(kring)] == ring::strucType::DDC) {
        ringType[static_cast<size_t>(kring)] = ring::strucType::bothPrismatic;
        listHC.push_back(kring);
      }
      isPrismatic[static_cast<size_t>(kring)] = true;
    }
    HCRings.clear();
    HCRings.push_back(iring);
    HCRings.push_back(jring);
    HCRings.insert(HCRings.end(), prismaticRings.begin(), prismaticRings.end());
    cageList.push_back({cage::cageType::HexC, HCRings});
  };

  struct Passing {
    int jring;
    bool cond2;
    bool reverse;
  };
  const int nThreads =
#if defined(SEAMS_HAS_OPENMP)
      omp_in_parallel() ? 1 : omp_get_max_threads();
#else
      1;
#endif
  const bool parallel = nThreads > 1;
  if (parallel) {
    // Hits stay on the thread that found them. A vector per ring shares
    // cache lines between threads and is built before any search starts.
    struct Hit {
      int iring;
      int jring;
      bool cond2;
      bool reverse;
    };
    std::vector<std::vector<Hit>> perThread(static_cast<std::size_t>(nThreads));
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel
#endif
    {
#ifdef SEAMS_HAS_OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      auto &local = perThread[static_cast<std::size_t>(tid)];
      std::vector<int> candidates;
      GatherScratch scratch;
#if defined(SEAMS_HAS_OPENMP)
#pragma omp for schedule(guided)
#endif
      for (int t = 0; t < totalRingNum; t++) {
        const int iring = visitOrder[static_cast<size_t>(t)];
        gatherCandidates(iring, candidates, scratch);
        for (const int jring : candidates) {
          const std::vector<int> &basal1 = rings[static_cast<size_t>(iring)];
          const std::vector<int> &basal2 = rings[static_cast<size_t>(jring)];
          if (ring::hasCommonElements(basal1, basal2)) {
            continue;
          }
          const bool cond2 = ring::basalConditions(bonds, basal1, basal2);
          const bool reverse = ring::basalConditions(bonds, basal2, basal1);
          if (cond2 || reverse) {
            local.push_back({iring, jring, cond2, reverse});
          }
        }
      }
    }
    std::vector<Hit> hits;
    for (auto &local : perThread) {
      hits.insert(hits.end(), local.begin(), local.end());
    }
    std::sort(hits.begin(), hits.end(), [&](const Hit &a, const Hit &b) {
      const int ra = rank[static_cast<size_t>(a.iring)];
      const int rb = rank[static_cast<size_t>(b.iring)];
      if (ra != rb) {
        return ra < rb;
      }
      return rank[static_cast<size_t>(a.jring)] < rank[static_cast<size_t>(b.jring)];
    });
    for (const Hit &hit : hits) {
      if (isPrismatic[static_cast<size_t>(hit.iring)] ||
          isPrismatic[static_cast<size_t>(hit.jring)]) {
        continue;
      }
      claimPair(hit.iring, hit.jring, hit.cond2, hit.reverse);
    }
  } else {
    std::vector<int> candidates;
    GatherScratch scratch;
    for (const int iring : visitOrder) {
      if (isPrismatic[static_cast<size_t>(iring)]) {
        continue;
      }
      gatherCandidates(iring, candidates, scratch);
      const std::vector<int> &basal1 = rings[static_cast<size_t>(iring)];
      for (const int jring : candidates) {
        if (isPrismatic[static_cast<size_t>(jring)]) {
          continue;
        }
        const std::vector<int> &basal2 = rings[static_cast<size_t>(jring)];
        if (ring::hasCommonElements(basal1, basal2)) {
          continue;
        }
        const bool cond2 = ring::basalConditions(bonds, basal1, basal2);
        const bool reverse = ring::basalConditions(bonds, basal2, basal1);
        if (!cond2 && !reverse) {
          continue;
        }
        claimPair(iring, jring, cond2, reverse);
      }
    }
  }

  sort(listHC.begin(), listHC.end());
  auto ip = std::unique(listHC.begin(), listHC.end());
  // Resize peripheral rings to remove undefined terms
  listHC.resize(static_cast<std::size_t>(std::distance(listHC.begin(), ip)));

  return listHC;
}

/**
 * @details Check to see if two basal rings are HCs or not, using the neighbour
 * list information. The neighbour list nList is a vector of vectors, containing
 * atom indices (not atom IDs!). The first element of each subvector in nList is
 * the atom index of the particle for which the other elements are the nearest
 * neighbours.
 * Internally, the following functions are called:
 * - ring::basalNeighbours (CONDITION1: @f$ m_{k+2} @f$ and @f$ m_{k+4} @f$ must
 * be bonded to @f$ l_3 @f$ and @f$ l_5 @f$ (if @f$ l_1 @f$ is a neighbour) or
 * @f$ m_{k+2} @f$ and @f$ m_{k+4} @f$ must be bonded to @f$ l_4 @f$ and @f$ l_6
 * @f$ (if @f$ l_2 @f$ is a neighbour)).
 * - ring::notNeighboursOfRing (CONDITION2: @f$ m_{k+1} @f$, @f$ m_{k+3} @f$ and
 * @f$ m_{k+5}@f$ must NOT be bonded to any element in basal1).
 * @param[in] nList Vector of vectors of the neighbour list (by index).
 * @param[in] basal1 Vector containing the first candidate basal ring.
 * @param[in] basal2 Vector containing the second candidate basal ring.
 * @return A bool; true if the basal rings being tested fulfill this condition
 *  for being the basal rings of an HC.
 */
void ring::NeighBond::build(const std::vector<std::vector<int>> &nList) {
  n_ = static_cast<int>(nList.size());
  off_.assign(static_cast<size_t>(n_) + 1, 0);
  cols_.clear();
  std::vector<int> row;
  for (int i = 0; i < n_; i++) {
    row.clear();
    if (static_cast<size_t>(i) < nList.size()) {
      const auto &nbrs = nList[static_cast<size_t>(i)];
      row.reserve(nbrs.size());
      for (size_t k = 1; k < nbrs.size(); k++) {
        const int j = nbrs[k];
        if (j >= 0 && j < n_) {
          row.push_back(j);
        }
      }
    }
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
    off_[static_cast<size_t>(i)] = static_cast<int>(cols_.size());
    cols_.insert(cols_.end(), row.begin(), row.end());
  }
  off_[static_cast<size_t>(n_)] = static_cast<int>(cols_.size());
}

bool ring::basalConditions(const std::vector<std::vector<int>> &nList,
                           const std::vector<int> &basal1,
                           const std::vector<int> &basal2) {
  // An empty neighbour row has no bonds, which fails the one-to-one test
  return basalConditionsOf(
      [&](int left, int right) {
        if (left < 0 || left >= static_cast<int>(nList.size()) ||
            nList[static_cast<size_t>(left)].empty()) {
          return false;
        }
        const auto &row = nList[static_cast<size_t>(left)];
        return std::find(row.begin() + 1, row.end(), right) != row.end();
      },
      basal1, basal2);
}

bool ring::basalConditions(const ring::NeighBond &bonds,
                           const std::vector<int> &basal1,
                           const std::vector<int> &basal2) {
  if (basal1.size() != 6 || basal2.size() != 6) {
    return false;
  }
  // One sequential pass of each neighbour row against the six atoms of
  // basal2. Thirty-six binary searches asked the same question and jumped
  // through the row each time.
  int atom[6];
  int pos[6];
  for (int j = 0; j < 6; j++) {
    atom[j] = basal2[static_cast<std::size_t>(j)];
    pos[j] = j;
  }
  for (int a = 1; a < 6; a++) {
    const int key = atom[a];
    const int at = pos[a];
    int b = a;
    while (b > 0 && atom[b - 1] > key) {
      atom[b] = atom[b - 1];
      pos[b] = pos[b - 1];
      b--;
    }
    atom[b] = key;
    pos[b] = at;
  }
  bool bond[6][6] = {};
  for (int i = 0; i < 6; i++) {
    const int *p = nullptr;
    const int *end = nullptr;
    bonds.row(basal1[static_cast<std::size_t>(i)], p, end);
    if (p == nullptr) {
      continue;
    }
    int t = 0;
    while (p != end && t < 6) {
      if (*p == atom[t]) {
        bond[i][pos[t]] = true;
        ++t;
        if (t >= 6 || atom[t] != *p) {
          ++p;
        }
      } else if (*p < atom[t]) {
        ++p;
      } else {
        ++t;
      }
    }
  }
  return basalFromMatrix(bond);
}

/**
 * @details Tests whether the last two elements of a triplet are neighbours of
 *two atom indices which have been passed in as inputs.
 * @param[in] nList Vector of vectors of the neighbour list (by index).
 * @param[in] triplet Vector containing the current triplet being tested.
 * @param[in] atomOne Index of the first atom.
 * @param[in] atomTwo Index of the second atom.
 * @return A bool; true if the condition is met and false otherwise.
 */
bool ring::basalNeighbours(const std::vector<std::vector<int>> &nList,
                           std::vector<int> &triplet, int atomOne,
                           int atomTwo) {
  // Search for needles in a haystack :)
  int needle1 = triplet[1];
  int needle2 = triplet[2];
  bool neighbourFound = false;
  bool neighOne = false,
       neighTwo = false; // Neighbour for atomOne, or neighbour for atomTwo
  // ----------------------------
  // For first element needle1, which must belong to either atomOne's or
  // atomTwo's neighbour list Search atomOne's neighbours
  auto it =
      std::find(nList[static_cast<std::size_t>(atomOne)].begin() + 1, nList[static_cast<std::size_t>(atomOne)].end(), needle1);
  if (it != nList[static_cast<std::size_t>(atomOne)].end()) {
    neighbourFound = true;
    neighOne = true;
  } // atomOne's neighbour
  // If it is not atomOne's neighbour, it might be atomTwo's neighbour
  if (!neighOne) {
    it = std::find(nList[static_cast<std::size_t>(atomTwo)].begin() + 1, nList[static_cast<std::size_t>(atomTwo)].end(), needle1);
    if (it != nList[static_cast<std::size_t>(atomTwo)].end()) {
      neighbourFound = true;
      neighTwo = true;
    } // end of check to see if neighbour was found
  }   // End of check to see if needle1 is atomTwo's neighbour
  // ----------------------------
  // If needle1 is not a neighbour of atomOne or atomTwo, return false
  if (!neighbourFound) {
    return false;
  }

  // Check to see if needle2 is a neighbour of either atomOne or atomTwo
  // ===============
  // if atomOne was a neighbour of needle1, needle2 must be a neighbour of
  // atomTwo
  if (neighOne) {
    it = std::find(nList[static_cast<std::size_t>(atomTwo)].begin() + 1, nList[static_cast<std::size_t>(atomTwo)].end(), needle2);
    // It is a neighbour
    if (it != nList[static_cast<std::size_t>(atomTwo)].end()) {
      return true;
    }
    // It is not a neighbour
    else {
      return false;
    }
  } // End of check for neighbour of atomTwo
  // ===============
  // if atomTwo was a neighbour of needle1, needle2 must be a neighbour of
  // atomOne
  else {
    it = std::find(nList[static_cast<std::size_t>(atomOne)].begin() + 1, nList[static_cast<std::size_t>(atomOne)].end(), needle2);
    // It is a neighbour
    if (it != nList[static_cast<std::size_t>(atomOne)].end()) {
      return true;
    }
    // It is not a neighbour
    else {
      return false;
    }
  }
  // ===============
}

/**
 * @details Checks to make sure that the elements of the triplet are NOT
 * neighbours of any elements inside a vector (ring) passed in (false)
 * If any of them are neighbours, this function returns false.
 * @param[in] nList Vector of vectors of the neighbour list (by index).
 * @param[in] triplet Vector containing the current triplet being tested.
 * @param[in] ring Ring passed in.
 * @return A bool; true if the condition is met and false otherwise.
 */
bool ring::notNeighboursOfRing(const std::vector<std::vector<int>> &nList,
                               std::vector<int> &triplet,
                               const std::vector<int> &ring) {
  int iatom; // AtomID of the atom to be searched for inside the neighbour
             // lists
  int jatom; // AtomID of in whose neighbour list iatom will be searched for
  std::vector<int>::const_iterator it;

  for (int i = 0; i < triplet.size(); i++) {
    iatom = triplet[static_cast<std::size_t>(i)]; // AtomID to be searched for
    // iatom must be searched for in the neighbour lists of all elements
    // of the ring vector
    for (int j = 0; j < ring.size(); j++) {
      jatom = ring[static_cast<std::size_t>(j)];
      // ------------------
      // Search for iatom in the neighbour list of jatom
      it = std::find(nList[static_cast<std::size_t>(jatom)].begin() + 1, nList[static_cast<std::size_t>(jatom)].end(), iatom);
      // It is a neighbour!
      if (it != nList[static_cast<std::size_t>(jatom)].end()) {
        return false;
      }
      // ------------------
    } // end of loop through every element of ring
  }   // end of loop through all triplet elements

  return true;
}

/**
 * @details Finding prismatic rings when passed the information in the
 *ringType input vector.
 * @param[in] rings The 6-membered primitive rings.
 * @param[in] listHC Vector containing the ring indices of rings which are part
 *  of HCs.
 * @param[in] ringType Contains a structure type for each ring.
 * @param[in] iring Index of the \f$ i^{th} \f$ ring.
 * @param[in] jring Index of the \f$ j^{th} \f$ ring.
 * @param[in] prismaticRings Vector containing the indices of rings which are
 *  prismatic.
 */
int ring::findPrismatic(const std::vector<std::vector<int>> &rings,
                        std::vector<int> &listHC,
                        std::vector<ring::strucType> &ringType, int iring,
                        int jring, std::vector<int> &prismaticRings) {
  int maxAtom = 0;
  for (const auto &r : rings) {
    for (const int atom : r) {
      maxAtom = std::max(maxAtom, atom);
    }
  }
  const auto index = ring::buildRingSearchIndex(rings, maxAtom + 1);
  return ring::findPrismatic(rings, listHC, ringType, iring, jring,
                             prismaticRings, index);
}

int ring::findPrismatic(const std::vector<std::vector<int>> &rings,
                        std::vector<int> &listHC,
                        std::vector<ring::strucType> &ringType, int iring,
                        int jring, std::vector<int> &prismaticRings,
                        const ring::RingSearchIndex &index) {
  int iIndex;                     // Used for making rings to be searched
  int ringSize = static_cast<int>(rings[0].size()); // This is 6 for hexagons
  std::vector<int> iTriplet;      // triplet formed from iring
  std::vector<int> jTriplet;      // triplet formed from jring
  std::vector<int> common;        // Common elements
  std::vector<int> candidates;

  // Make all possible triplets out of iring
  for (int i = 0; i < ringSize; i++) {
    // init
    iTriplet.clear();
    // ------
    // Get a triplet
    for (int m = 0; m < 3; m++) {
      iIndex = i + m;
      if (iIndex >= ringSize) {
        iIndex -= ringSize;
      }
      iTriplet.push_back(rings[static_cast<std::size_t>(iring)][static_cast<std::size_t>(iIndex)]);
    } // end of getting a triplet from iring

    candidates.clear();
    ringsThrough(index, iTriplet[0], iTriplet[1], iTriplet[2], iring, jring,
                 candidates);
    for (const int kring : candidates) {

      // -----------------
      // If kring does have the triplet, check to see if at least three other
      // elements of kring are shared by jring
      jTriplet.clear(); // init
      // Make jTriplet
      for (int j = 0; j < ringSize; j++) {
        int katom = rings[static_cast<std::size_t>(kring)][static_cast<std::size_t>(j)];
        auto it = std::find(iTriplet.begin(), iTriplet.end(), katom);

        // If not found, add it to jTriplet
        if (it == iTriplet.end()) {
          jTriplet.push_back(katom);
        } // update jTriplet
      }   // end of making jTriplet out of kring
      // -----------------

      // Now search for jTriplet inside jring
      common = ring::findsCommonElements(jTriplet, rings[static_cast<std::size_t>(jring)]);

      // Update the prismatic rings
      if (common.size() == 3) {
        //
        listHC.push_back(kring);         // Update listHC vector
        prismaticRings.push_back(kring); // Update prismatic rings
        // Update the type inside ringType
        // If the ring is already a DDC ring, it is a mixed ring
        if (ringType[static_cast<std::size_t>(kring)] == ring::strucType::DDC) {
          ringType[static_cast<std::size_t>(kring)] = ring::strucType::bothPrismatic;
        }
        // If it is unclassified, it is just a prismatic ring
        if (ringType[static_cast<std::size_t>(kring)] == ring::strucType::unclassified) {
          ringType[static_cast<std::size_t>(kring)] = ring::strucType::HCprismatic;
        } // end ring update
      }   // add kring to the list of prismatic rings
    }     // end of searching through rings for kring
    // -------------------------------------------
  } // end of getting pairs of triplets

  return 0;
}

/**
 * @details Determines which hexagonal rings are both DDCs and HCs. This
 * function returns a vector which contains the ring IDs of all the rings which
 * are both. The ring IDs correspond to the index of the rings inside the vector
 * of vector rings, starting from 0. Rings which are both are of enum type
 * bothBasal OR bothPrismatic. Reassign rings which are mixed in listDDC and
 * listHC as the dummy value -10.
 * @param[in] rings The 6-membered primitive rings.
 * @param[in] ringType Structure type for every ring.
 * @param[in] listDDC Contains a vector of ring indices of DDCs.
 * @param[in] listHC Contains a vector of ring indices of HCs.
 * @return A vector containing the indices of rings which are mixed.
 */
std::vector<int> ring::findMixedRings(const std::vector<std::vector<int>> &rings,
                                      std::vector<strucType> &ringType,
                                      std::vector<int> &listDDC,
                                      std::vector<int> &listHC) {
  std::vector<int> listMixed;
  int dummyValue = -10;

  // Loop through all rings in the ringType and
  // adds the ring Indices of all rings which are both DDCs and HCs
  for (int iring = 0; iring < ringType.size(); iring++) {
    // If iring is of mixed type, add it to the listMixed vector
    if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::bothBasal ||
        ringType[static_cast<std::size_t>(iring)] == ring::strucType::bothPrismatic) {
      listMixed.push_back(iring);

      //-----------------
      // Search for iring in listDDC
      std::sort(listDDC.begin(), listDDC.end());
      auto iter = std::find(listDDC.begin(), listDDC.end(), iring);
      if (iter != listDDC.end()) {
        *iter = dummyValue; // Assign dummy value to the mixed ring
      }                     // found in listDDC
      //-----------------
      //-----------------
      // Search for iring in listHC
      std::sort(listHC.begin(), listHC.end());
      auto itr = std::find(listHC.begin(), listHC.end(), iring);
      if (itr != listHC.end()) {
        *itr = dummyValue; // Assign dummy value to the mixed ring
      }                    // found in listHC
      //-----------------

    } // end of check for type
  }   // end of loop through all every ring

  return listMixed;
}

/**
 * @details Assigns a type (cage::iceType) to each atom, according to the
 * previous classification of every ring in ringType. Each subring or vector
 * inside the vector of vector rings, is by index, meaning that the atoms are
 * saved by their indices starting from 0 in the PointCloud.
 * @param[in] rings Vector of vectors of 6-membered rings.
 * @param[in] ringType Vector containing the structural type for each ring.
 * @param[in] atomTypes Structural type for each atom.
 */
int ring::getAtomTypesTopoBulk(const std::vector<std::vector<int>> &rings,
                               const std::vector<ring::strucType> &ringType,
                               std::vector<cage::iceType> &atomTypes) {
  cage::iceType iRingType;        // Type of the current ring
  int iatom;                      // Current ring
  int ringSize = static_cast<int>(rings[0].size()); // Size of the ring

  // Loop through every ring in ringType
  for (int iring = 0; iring < ringType.size(); iring++) {
    //
    // Skip if the ring is unclassified
    if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::unclassified) {
      continue;
    } // skip for unclassified rings

    // ------------
    // Get the current ring type
    // DDC
    if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::DDC) {
      iRingType = cage::iceType::ddc;
    } // DDC atoms
    //
    // HC
    else if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::HCbasal ||
             ringType[static_cast<std::size_t>(iring)] == ring::strucType::HCprismatic) {
      iRingType = cage::iceType::hc;
    } // HC atoms
    //
    // Mixed
    else if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::bothBasal ||
             ringType[static_cast<std::size_t>(iring)] == ring::strucType::bothPrismatic) {
      iRingType = cage::iceType::mixed;
    } // HC atoms
    // Prism
    else if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::Prism ||
             ringType[static_cast<std::size_t>(iring)] == ring::strucType::deformedPrism ||
             ringType[static_cast<std::size_t>(iring)] == ring::strucType::mixedPrismRing) {
      iRingType = cage::iceType::pnc; // 5 membered pnc
    }                        // prism
    // Should never go here
    else {
      continue;
    } //
    // ------------
    // Otherwise, loop through every inside the ring and assign atomTypes the
    // iRingType
    for (int i = 0; i < ringSize; i++) {
      iatom = rings[static_cast<std::size_t>(iring)][static_cast<std::size_t>(i)]; // Atom index in ring
      if (atomTypes[static_cast<std::size_t>(iatom)] == cage::iceType::mixed ||
          atomTypes[static_cast<std::size_t>(iatom)] == cage::iceType::mixed2) {
        continue;
      } // Don't reassign
      // For atoms shared by PNCs and DDCs/HCs
      if (ringSize == 6) {
        if (atomTypes[static_cast<std::size_t>(iatom)] == cage::iceType::pnc) {
          atomTypes[static_cast<std::size_t>(iatom)] = cage::iceType::mixed2;
        } else {
          atomTypes[static_cast<std::size_t>(iatom)] = iRingType;
        }
      } else {
        atomTypes[static_cast<std::size_t>(iatom)] = iRingType;
      }
    } // end of loop thorugh the current ring
  }   // end of loop through every ring

  return 0;
}

/**
 * @details Determines the number of HCs, DDCs from the cageList vector,
 * containing a list of cages.
 * The number of mixed rings, prismatic rings and basal rings are obtained
 * from the ringType vector.
 * @param[in] ringType Vector containing the structural type for each ring.
 * @param[in] cageList Contains all the cages (DDCs and HCs).
 * @param[in] numHC The number of HCs.
 * @param[in] numDDC The number of DDCs.
 * @param[in] mixedRings The number of mixed rings.
 * @param[in] prismaticRings The number of prismatic rings (of HCs/mixed
 *  rings).
 * @param[in] basalRings TThe number of basal rings (of HCs/mixed
 *  rings).
 */
int ring::getStrucNumbers(const std::vector<ring::strucType> &ringType,
                          const std::vector<cage::Cage> &cageList, int &numHC,
                          int &numDDC, int &mixedRings, int &prismaticRings,
                          int &basalRings) {
  // Determines the number of HCs, DDCs, Mixed rings, prismatic and basal rings
  // Init
  numHC = 0;
  numDDC = 0;
  mixedRings = 0;
  prismaticRings = 0;
  basalRings = 0;
  // ------------------------------------
  // GETTING THE CAGES (DDCs and HCs)
  // Loop through cages
  for (int icage = 0; icage < cageList.size(); icage++) {
    //
    // HC
    if (cageList[static_cast<std::size_t>(icage)].type == cage::cageType::HexC) {
      numHC += 1;
    } // end of updating HC number
    //
    // DDC
    if (cageList[static_cast<std::size_t>(icage)].type == cage::cageType::DoubleDiaC) {
      numDDC += 1;
    } // end of updating DDC number
  }   // end of loop through cages
  // ------------------------------------
  // GETTING THE RINGSS (Mixed, Prismatic and Basal rings)
  // Loop through the rings
  for (int iring = 0; iring < ringType.size(); iring++) {
    // Mixed
    if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::bothBasal ||
        ringType[static_cast<std::size_t>(iring)] == ring::strucType::bothPrismatic) {
      mixedRings += 1;
      // Also update basal rings
      if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::bothBasal) {
        basalRings += 1;
      } // mixed basal rings
      // Also update prismatic rings
      if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::bothPrismatic) {
        prismaticRings += 1;
      } // mixed prismatic rings
    }   // end of updating mixed
    //
    // HCs
    if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::HCprismatic) {
      prismaticRings += 1;
    } // HC prismatic
    // basal HCs
    if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::HCbasal) {
      basalRings += 1;
    } // HC basal
  }   // end of loop through every ring
  // ------------------------------------

  return 0;
} // end of function

/**
 * @details Determines which rings are n-sided prisms. This function
 * returns a vector which contains the ring indices of all the rings which are
 * prisms. The ring indices correspond to the index of the rings inside the
 * vector of vector rings, starting from 0. Prism rings can be found using a
 * three-step procedure, in which first two basal rings are found. Prismatic
 * rings are simply rings which share every face made by upper and lower
 * triplets of the basal rings The neighbour list is also required as an input,
 * which is a vector of vectors, containing atom IDs. The first element of the
 * neighbour list is the atom index of
 * the atom for which the other elements are nearest neighbours.\
 * @param[in] rings The input vector of vectors containing the primitive rings
 *  of a single ring size (number of nodes).
 * @param[in] ringType A vector containing a ring::strucType value (a
 *  classification type) for each ring.
 * @param[in] nPrisms The number of prism blocks identified for the number of
 *  nodes.
 * @param[in] nList The row-ordered neighbour list (by atom index).
 * @param[in] yCloud The input PointCloud.
 * @return A vector containing the ring indices of all the rings which have
 *  been classified as prisms. The indices are with respect to the input rings
 *  vector of vectors.
 */
int prism3::findBulkPrisms(
    const std::vector<std::vector<int>> &rings, std::vector<ring::strucType> &ringType,
    const std::vector<std::vector<int>> &nList,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    std::vector<double> &rmsdPerAtom, double heightCutoff) {
  int totalRingNum = static_cast<int>(rings.size()); // Total number of rings
  std::vector<int> basal1;         // First basal ring
  std::vector<int> basal2;         // Second basal ring
  bool cond1, cond2; // Conditions for rings to be basal (true) or not (false)
  int ringSize = static_cast<int>(rings[0].size()); // Number of nodes in each ring
  // Matrix for the reference ring for a given ringSize.
  Eigen::MatrixXd refPointSet(ringSize, 3);

  int axialDim = 2; // Default=z
  refPointSet = pntToPnt::getPointSetRefRing(ringSize, axialDim);
  //

  // Two loops through all the rings are required to find pairs of basal rings
  for (int iring = 0; iring < totalRingNum - 1; iring++) {
    cond1 = false;
    cond2 = false;
    basal1 = rings[static_cast<std::size_t>(iring)]; // Assign iring to basal1
    // Loop through the other rings to get a pair
    for (int jring = iring + 1; jring < totalRingNum; jring++) {
      basal2 = rings[static_cast<std::size_t>(jring)]; // Assign jring to basal2
      // ------------
      // Put extra check for axial basal rings if shapeMatching is being done
      // ------------
      // Step one: Check to see if basal1 and basal2 have common
      // elements or not. If they don't, then they cannot be basal rings
      cond1 = ring::hasCommonElements(basal1, basal2);
      if (cond1) {
        continue;
      }

      // ------------
      bool smallDist =
          prism3::basalRingsSeparation(yCloud, basal1, basal2, heightCutoff);
      if (!smallDist) {
        continue;
      } // the basal rings are too far apart

      // Otherwise
      // Do shape matching here
      bool isPrism = match::matchUntetheredPrism(yCloud, nList, refPointSet,
                                                 basal1, basal2, rmsdPerAtom);

      // Success! The rings are basal rings of a prism!
      if (isPrism) {
        //
        // Update iring
        if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::unclassified) {
          ringType[static_cast<std::size_t>(iring)] = ring::strucType::Prism;
        }
        // Update jring
        if (ringType[static_cast<std::size_t>(jring)] == ring::strucType::unclassified) {
          ringType[static_cast<std::size_t>(jring)] = ring::strucType::Prism;
        }
      } // end of reduced criteria
      // Strict criteria
      else {
        cond2 = prism3::basalPrismConditions(nList, basal1, basal2);
        // If the condition is false then the strict criterion has not been met
        if (!cond2) {
          continue;
        }
        // Update iring
        if (ringType[static_cast<std::size_t>(iring)] == ring::strucType::unclassified) {
          ringType[static_cast<std::size_t>(iring)] = ring::strucType::Prism;
        }
        // Update jring
        if (ringType[static_cast<std::size_t>(jring)] == ring::strucType::unclassified) {
          ringType[static_cast<std::size_t>(jring)] = ring::strucType::Prism;
        }
        //
        // Shape-matching to get the RMSD (if shape-matching is desired)

        // bool isKnownPrism = match::matchPrism(
        //     yCloud, nList, refPointSet, basal1, basal2, rmsdPerAtom, true);

        // -----------
      } // end of strict criteria

    } // end of loop through rest of the rings to get the second basal ring
  }   // end of loop through all rings for first basal ring

  return 0;
}

/**
 * @details A function that checks to see if two basal rings are basal rings of
 * a prism block or not, using the neighbour list information. The neighbour
 * list nList is a row-ordered vector of vectors, containing atom indices (not
 * atom IDs!). The first element of each subvector in nList is the atom index of
 * the particle for which the other elements are the nearest neighbours.
 * @param[in] nList Row-ordered neighbour list by atom index.
 * @param[in] basal1 The vector for one of the basal rings.
 * @param[in] basal2 The vector for the other basal ring.
 * @return A value that is true if the basal rings constitute a prism block,
 *  and false if they do not make up a prism block.
 */
bool prism3::basalPrismConditions(const std::vector<std::vector<int>> &nList,
                                  std::vector<int> &basal1,
                                  std::vector<int> &basal2) {
  if (basal1.empty() || basal2.empty() || basal1.size() != basal2.size()) {
    return false;
  }
  int l1 = basal1[0]; // first element of basal1 ring
  int ringSize =
      static_cast<int>(basal1.size()); // Size of the ring; each ring contains n elements
  int m_k;              // Atom ID of element in basal2
  bool l1_neighbour = false; // m_k is a neighbour of l1(true) or not (false)

  // isNeighbour is initialized to false for all basal2 elements; indication if
  // basal2 elements are neighbours of basal1
  std::vector<bool> isNeighbour(static_cast<std::size_t>(ringSize), false);
  int kIndex = 0; // m_k index
  int lAtomID; // atomID of the current element of basal1
  int kAtomID; // atomID of the current element of basal2

  // ---------------------------------------------
  // COMPARISON OF basal2 ELEMENTS WITH l1
  for (int k = 0; k < ringSize; k++) {
    l1_neighbour = false;
    m_k = basal2[static_cast<std::size_t>(k)];
    // =================================
    // Checking to seee if m_k is be a neighbour of l1
    // Find m_k inside l1 neighbour list
    auto it = std::find(nList[static_cast<std::size_t>(l1)].begin() + 1, nList[static_cast<std::size_t>(l1)].end(), m_k);

    // If the element has been found, for l1
    if (it != nList[static_cast<std::size_t>(l1)].end()) {
      l1_neighbour = true;
      kIndex = k;
      break;
    }
  } // l1 is a neighbour of m_k

  // If there is no nearest neighbour, then the two rings are not part of the
  // prism
  if (!l1_neighbour) {
    return false;
  }

  // ---------------------------------------------
  // NEIGHBOURS of basal1 in basal2
  isNeighbour[static_cast<std::size_t>(kIndex)] = true;

  // All elements of basal1 must be neighbours of basal2
  for (int i = 1; i < ringSize; i++) {
    lAtomID = basal1[static_cast<std::size_t>(i)]; // element of basal1 ring
    for (int k = 0; k < ringSize; k++) {
      // Skip if already a neighbour
      if (isNeighbour[static_cast<std::size_t>(k)]) {
        continue;
      }
      // Get the comparison basal2 element
      kAtomID = basal2[static_cast<std::size_t>(k)]; // element of basal2 ring;

      // Checking to see if kAtomID is a neighbour of lAtomID
      // Find kAtomID inside lAtomID neighbour list
      auto it1 =
          std::find(nList[static_cast<std::size_t>(lAtomID)].begin() + 1, nList[static_cast<std::size_t>(lAtomID)].end(), kAtomID);

      // If the element has been found, for l1
      if (it1 != nList[static_cast<std::size_t>(lAtomID)].end()) {
        isNeighbour[static_cast<std::size_t>(k)] = true;
      }
    } // Loop through basal2
  }   // end of check for neighbours of basal1

  // ---------------------------------------------

  // They should all be neighbours
  for (int k = 0; k < ringSize; k++) {
    // Check to see if any element is false
    if (!isNeighbour[static_cast<std::size_t>(k)]) {
      return false;
    }
  }

  // Everything works out!
  return true;
}

/**
 * @details Relaxed criteria for deformed
 * prism blocks: at least one bond
 * should exist between the basal
 * rings.
 */
bool prism3::relaxedPrismConditions(const std::vector<std::vector<int>> &nList,
                                    std::vector<int> &basal1,
                                    std::vector<int> &basal2) {
  int ringSize =
      static_cast<int>(basal1.size());     // Size of the ring; each ring contains n elements
  int m_k;                  // Atom ID of element in basal2
  bool isNeighbour = false; // This is true if there is at least one bond
                            // between the basal rings
  int l_k;                  // Atom ID of element in basal1

  // ---------------------------------------------
  // COMPARISON OF basal2 ELEMENTS (m_k) WITH basal1 ELEMENTS (l_k)
  // Loop through all the elements of basal1
  for (int l = 0; l < ringSize; l++) {
    l_k = basal1[static_cast<std::size_t>(l)];
    // Search for the nearest neighbour of l_k in basal2
    // Loop through basal2 elements
    for (int m = 0; m < ringSize; m++) {
      m_k = basal2[static_cast<std::size_t>(m)];
      // Find m_k inside l_k neighbour list
      auto it = std::find(nList[static_cast<std::size_t>(l_k)].begin() + 1, nList[static_cast<std::size_t>(l_k)].end(), m_k);

      // If the element has been found, for l1
      if (it != nList[static_cast<std::size_t>(l_k)].end()) {
        isNeighbour = true;
        break;
      } // found element
    }   // end of loop through all the elements of basal2

    // If a neighbour has been found then
    if (isNeighbour) {
      return true;
    }
  } // end of loop through all the elements of basal1

  // If a neighbour has not been found, return false
  return false;
}

//! Check to see that candidate basal prisms are not really far from each other
//! Return true if the basal rings are within the heightCutoff
bool prism3::basalRingsSeparation(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &basal1, const std::vector<int> &basal2, double heightCutoff) {
  //
  int ringSize = static_cast<int>(basal1.size());
  int l_k, m_k; // Atom indices
  double infHugeNumber = 100000;
  double leastDist = infHugeNumber;
  int index = -1; // starting index
  // For the first element of basal1:

  l_k = basal1[0]; // This is the atom particle C++ index

  // Search for the nearest neighbour of l_k in basal2
  // Loop through basal2 elements
  for (int m = 0; m < ringSize; m++) {
    m_k = basal2[static_cast<std::size_t>(m)]; // Atom index to find in the neighbour list of iatom

    // Calculate the distance
    double dist = gen::periodicDist(yCloud, l_k, m_k);

    // Update the least distance
    if (leastDist > dist) {
      leastDist = dist; // This is the new least distance
      index = m;
    } // end of update of the least distance

  } // found element

  // If the element has been found, for l1
  if (leastDist < heightCutoff) {
    return true;
  } // end of check
  else {
    return false;
  }
}
