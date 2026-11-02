#include <topo_bulk.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

struct Row {
  int v[8]{};
  int n = 0;
  bool push(int x) {
    if (n >= 8) {
      return false;
    }
    v[n++] = x;
    return true;
  }
  void sortUnique() {
    std::sort(v, v + n);
    n = static_cast<int>(std::unique(v, v + n) - v);
  }
  [[nodiscard]] bool has(int x) const {
    for (int i = 0; i < n; i++) {
      if (v[i] == x) {
        return true;
      }
    }
    return false;
  }
};

struct Bits {
  static constexpr int kWords = 4;
  static constexpr int kCap = 256;
  std::uint64_t w[kWords]{};

  void set(int i) { w[i >> 6] |= 1ULL << (i & 63); }
  void reset(int i) { w[i >> 6] &= ~(1ULL << (i & 63)); }
  [[nodiscard]] bool test(int i) const { return (w[i >> 6] >> (i & 63)) & 1ULL; }
  [[nodiscard]] int count() const {
    int n = 0;
    for (int i = 0; i < kWords; i++) {
      n += __builtin_popcountll(w[i]);
    }
    return n;
  }
  [[nodiscard]] bool empty() const {
    for (int i = 0; i < kWords; i++) {
      if (w[i]) {
        return false;
      }
    }
    return true;
  }
  Bits operator&(const Bits &o) const {
    Bits r;
    for (int i = 0; i < kWords; i++) {
      r.w[i] = w[i] & o.w[i];
    }
    return r;
  }
  Bits operator-(const Bits &o) const {
    Bits r;
    for (int i = 0; i < kWords; i++) {
      r.w[i] = w[i] & ~o.w[i];
    }
    return r;
  }
  [[nodiscard]] int first() const {
    for (int i = 0; i < kWords; i++) {
      if (w[i]) {
        return i * 64 + __builtin_ctzll(w[i]);
      }
    }
    return -1;
  }
};

struct MisOut {
  std::vector<int> ids;
  int cover = 0;
  bool solved = false;
};

// Exact maximum independent set on one component. Local indices follow the
// canonical rank order of mem. A clique partition is an upper bound; the
// search stops when that bound cannot beat the incumbent.
MisOut exactMis(const std::vector<Row> &conf, const std::vector<int> &mem) {
  MisOut out;
  const int n = static_cast<int>(mem.size());
  if (n == 0) {
    out.solved = true;
    return out;
  }
  if (n > Bits::kCap) {
    return out;
  }
  std::vector<int> at(conf.size(), -1);
  for (int i = 0; i < n; i++) {
    at[static_cast<size_t>(mem[static_cast<size_t>(i)])] = i;
  }
  std::vector<Bits> nbr(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    nbr[static_cast<size_t>(i)].set(i);
    const int g = mem[static_cast<size_t>(i)];
    const Row &row = conf[static_cast<size_t>(g)];
    for (int k = 0; k < row.n; k++) {
      const int j = at[static_cast<size_t>(row.v[k])];
      if (j >= 0) {
        nbr[static_cast<size_t>(i)].set(j);
      }
    }
  }
  std::vector<int> cliqueOf(static_cast<size_t>(n), -1);
  std::vector<char> used(static_cast<size_t>(n), 0);
  std::vector<int> grow(static_cast<size_t>(n));
  int nCliques = 0;
  std::vector<int> seq(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    seq[static_cast<size_t>(i)] = i;
  }
  std::sort(seq.begin(), seq.end(), [&](int a, int b) {
    const int da = nbr[static_cast<size_t>(a)].count();
    const int db = nbr[static_cast<size_t>(b)].count();
    if (da != db) {
      return da > db;
    }
    return a < b;
  });
  for (const int seed : seq) {
    if (used[static_cast<size_t>(seed)]) {
      continue;
    }
    int nMem = 0;
    grow[static_cast<size_t>(nMem++)] = seed;
    used[static_cast<size_t>(seed)] = 1;
    for (int u = 0; u < n; u++) {
      if (used[static_cast<size_t>(u)]) {
        continue;
      }
      bool ok = true;
      for (int t = 0; t < nMem; t++) {
        if (!nbr[static_cast<size_t>(grow[static_cast<size_t>(t)])].test(u)) {
          ok = false;
          break;
        }
      }
      if (ok) {
        grow[static_cast<size_t>(nMem++)] = u;
        used[static_cast<size_t>(u)] = 1;
      }
    }
    for (int t = 0; t < nMem; t++) {
      cliqueOf[static_cast<size_t>(grow[static_cast<size_t>(t)])] = nCliques;
    }
    nCliques++;
  }
  out.cover = nCliques;
  std::vector<int> stamp(static_cast<size_t>(std::max(nCliques, 1)), 0);
  int era = 0;
  auto bound = [&](Bits avail) {
    if (++era == 0) {
      std::fill(stamp.begin(), stamp.end(), 0);
      era = 1;
    }
    int b = 0;
    for (int v = 0; v < n; v++) {
      if (!avail.test(v)) {
        continue;
      }
      const int c = cliqueOf[static_cast<size_t>(v)];
      if (stamp[static_cast<size_t>(c)] != era) {
        stamp[static_cast<size_t>(c)] = era;
        b++;
      }
    }
    return b;
  };
  Bits all;
  for (int i = 0; i < n; i++) {
    all.set(i);
  }
  std::vector<int> path;
  std::vector<int> bestPath;
  int best = 0;
  {
    Bits left = all;
    while (!left.empty()) {
      const int v = left.first();
      path.push_back(v);
      left = left - nbr[static_cast<size_t>(v)];
    }
    best = static_cast<int>(path.size());
    bestPath = path;
    path.clear();
  }
  const auto rec = [&](auto &&self, Bits avail, int taken) -> void {
    const int mark = static_cast<int>(path.size());
    for (;;) {
      int leaf = -1;
      int leafN = -1;
      int zero = -1;
      for (int v = 0; v < n; v++) {
        if (!avail.test(v)) {
          continue;
        }
        const int deg = (nbr[static_cast<size_t>(v)] & avail).count() - 1;
        if (deg <= 0) {
          zero = v;
          break;
        }
        if (deg == 1 && leaf < 0) {
          Bits rest = nbr[static_cast<size_t>(v)] & avail;
          rest.reset(v);
          leaf = v;
          leafN = rest.first();
        }
      }
      if (zero >= 0) {
        avail.reset(zero);
        path.push_back(zero);
        taken++;
        continue;
      }
      if (leaf >= 0) {
        avail.reset(leaf);
        avail.reset(leafN);
        path.push_back(leaf);
        taken++;
        avail = avail - nbr[static_cast<size_t>(leafN)];
        continue;
      }
      break;
    }
    if (taken + bound(avail) <= best) {
      path.resize(static_cast<size_t>(mark));
      return;
    }
    if (avail.empty()) {
      if (taken > best) {
        best = taken;
        bestPath = path;
      }
      path.resize(static_cast<size_t>(mark));
      return;
    }
    int branch = -1;
    int bestDeg = -1;
    for (int v = 0; v < n; v++) {
      if (!avail.test(v)) {
        continue;
      }
      const int deg = (nbr[static_cast<size_t>(v)] & avail).count() - 1;
      if (deg > bestDeg) {
        bestDeg = deg;
        branch = v;
      }
    }
    path.push_back(branch);
    self(self, avail - nbr[static_cast<size_t>(branch)], taken + 1);
    path.pop_back();
    avail.reset(branch);
    self(self, avail, taken);
    path.resize(static_cast<size_t>(mark));
  };
  if (out.cover > best) {
    rec(rec, all, 0);
  }
  out.solved = true;
  out.ids.reserve(static_cast<size_t>(best));
  for (const int i : bestPath) {
    out.ids.push_back(mem[static_cast<size_t>(i)]);
  }
  return out;
}

struct PKey {
  std::uint64_t k;
  int ring;
  int slot;
};

std::uint64_t packRanks(const int *id, int n, const std::vector<int> &rankOf) {
  if (n <= 0 || n > 4) {
    return ~0ULL;
  }
  int mem[4];
  for (int i = 0; i < n; i++) {
    mem[i] = rankOf[static_cast<size_t>(id[i])];
    if (mem[i] < 0 || mem[i] >= 65536) {
      return ~0ULL;
    }
  }
  for (int i = 1; i < n; i++) {
    const int x = mem[i];
    int j = i;
    while (j > 0 && mem[j - 1] > x) {
      mem[j] = mem[j - 1];
      j--;
    }
    mem[j] = x;
  }
  std::uint64_t k = 0;
  for (int i = 0; i < n; i++) {
    k |= static_cast<std::uint64_t>(mem[i]) << (16 * i);
  }
  return k;
}

struct AdjE {
  int right;
  int ring;
};

// Every other edge of an Euler tour of a 2k-regular bipartite graph is a
// k-factor. One more pass on the resulting 2-factor is a perfect matching.
bool eulerPerfect(int nL, int nR, const std::vector<int> &off, const std::vector<AdjE> &adj,
                  std::vector<int> &matchL, std::vector<int> &matchR, std::vector<int> &matchRing) {
  const int m = static_cast<int>(adj.size());
  if (nL == 0 || nL != nR || m == 0) {
    return false;
  }
  const int d0 = off[1] - off[0];
  if (d0 <= 0 || (d0 % 2) != 0) {
    return false;
  }
  for (int u = 0; u < nL; u++) {
    if (off[static_cast<size_t>(u + 1)] - off[static_cast<size_t>(u)] != d0) {
      return false;
    }
  }
  // One scratch per thread: the fuzz harness calls findDDC concurrently.
  thread_local std::vector<int> rCount, rOff, rEid, rPut, owner, lIt, rIt, tour, cuts;
  thread_local std::vector<int> succ, succ2, succE, succE2;
  thread_local std::vector<char> used, vis;
  auto fit = [](std::vector<int> &v, int n) {
    if (static_cast<int>(v.size()) < n) {
      v.resize(static_cast<size_t>(n));
    }
  };
  fit(rCount, nR);
  std::fill(rCount.begin(), rCount.begin() + nR, 0);
  for (const AdjE &e : adj) {
    rCount[static_cast<size_t>(e.right)]++;
  }
  for (int r = 0; r < nR; r++) {
    if (rCount[static_cast<size_t>(r)] != d0) {
      return false;
    }
  }
  fit(rOff, nR + 1);
  rOff[0] = 0;
  for (int r = 0; r < nR; r++) {
    rOff[static_cast<size_t>(r + 1)] = rOff[static_cast<size_t>(r)] + rCount[static_cast<size_t>(r)];
  }
  fit(rEid, m);
  fit(rPut, nR + 1);
  std::copy(rOff.begin(), rOff.begin() + nR + 1, rPut.begin());
  fit(owner, m);
  for (int u = 0; u < nL; u++) {
    for (int e = off[static_cast<size_t>(u)]; e < off[static_cast<size_t>(u + 1)]; e++) {
      owner[static_cast<size_t>(e)] = u;
      rEid[static_cast<size_t>(rPut[static_cast<size_t>(adj[static_cast<size_t>(e)].right)]++)] = e;
    }
  }
  fit(lIt, nL);
  if (static_cast<int>(used.size()) < m) {
    used.resize(static_cast<size_t>(m));
  }
  std::fill(used.begin(), used.begin() + m, 0);
  for (int u = 0; u < nL; u++) {
    lIt[static_cast<size_t>(u)] = off[static_cast<size_t>(u)];
  }
  fit(rIt, nR + 1);
  std::copy(rOff.begin(), rOff.begin() + nR + 1, rIt.begin());
  tour.clear();
  if (tour.capacity() < static_cast<size_t>(m)) {
    tour.reserve(static_cast<size_t>(m));
  }
  cuts.clear();
  cuts.push_back(0);
  struct Fr {
    int u;
    int inEdge;
  };
  std::vector<Fr> st;
  const int nv = nL + nR;
  auto advance = [&](int u, int &eid, int &nxt) -> bool {
    if (u < nL) {
      int &p = lIt[static_cast<size_t>(u)];
      const int end = off[static_cast<size_t>(u + 1)];
      while (p < end && used[static_cast<size_t>(p)]) {
        p++;
      }
      if (p == end) {
        return false;
      }
      eid = p;
      nxt = nL + adj[static_cast<size_t>(p)].right;
      used[static_cast<size_t>(p)] = 1;
      p++;
      return true;
    }
    const int r = u - nL;
    int &p = rIt[static_cast<size_t>(r)];
    const int end = rOff[static_cast<size_t>(r + 1)];
    while (p < end && used[static_cast<size_t>(rEid[static_cast<size_t>(p)])]) {
      p++;
    }
    if (p == end) {
      return false;
    }
    eid = rEid[static_cast<size_t>(p)];
    nxt = owner[static_cast<size_t>(eid)];
    used[static_cast<size_t>(eid)] = 1;
    p++;
    return true;
  };
  for (int start = 0; start < nv; start++) {
    const bool pending = start < nL
                             ? lIt[static_cast<size_t>(start)] < off[static_cast<size_t>(start + 1)]
                             : rIt[static_cast<size_t>(start - nL)] < rOff[static_cast<size_t>(start - nL + 1)];
    if (!pending) {
      continue;
    }
    st.clear();
    st.push_back({start, -1});
    while (!st.empty()) {
      const int u = st.back().u;
      int eid = 0;
      int nxt = 0;
      if (!advance(u, eid, nxt)) {
        if (st.back().inEdge >= 0) {
          tour.push_back(st.back().inEdge);
        }
        st.pop_back();
      } else {
        st.push_back({nxt, eid});
      }
    }
    cuts.push_back(static_cast<int>(tour.size()));
  }
  if (static_cast<int>(tour.size()) != m) {
    return false;
  }
  // Alternate edges of each circuit: a 2-factor when every degree is 4.
  fit(succ, nv);
  fit(succ2, nv);
  fit(succE, nv);
  fit(succE2, nv);
  std::fill(succ.begin(), succ.begin() + nv, -1);
  std::fill(succ2.begin(), succ2.begin() + nv, -1);
  std::fill(succE.begin(), succE.begin() + nv, -1);
  std::fill(succE2.begin(), succE2.begin() + nv, -1);
  auto link = [&](int u, int w, int id) {
    if (succ[static_cast<size_t>(u)] < 0) {
      succ[static_cast<size_t>(u)] = w;
      succE[static_cast<size_t>(u)] = id;
    } else {
      succ2[static_cast<size_t>(u)] = w;
      succE2[static_cast<size_t>(u)] = id;
    }
  };
  for (size_t c = 0; c + 1 < cuts.size(); c++) {
    for (int i = cuts[c]; i < cuts[c + 1]; i += 2) {
      const int id = tour[static_cast<size_t>(i)];
      const int L = owner[static_cast<size_t>(id)];
      const int R = nL + adj[static_cast<size_t>(id)].right;
      if (succ2[static_cast<size_t>(L)] >= 0 || succ2[static_cast<size_t>(R)] >= 0) {
        return false;
      }
      link(L, R, id);
      link(R, L, id);
    }
  }
  if (static_cast<int>(vis.size()) < nv) {
    vis.resize(static_cast<size_t>(nv));
  }
  std::fill(vis.begin(), vis.begin() + nv, 0);
  int got = 0;
  for (int s = 0; s < nv; s++) {
    if (vis[static_cast<size_t>(s)] || succ[static_cast<size_t>(s)] < 0 ||
        succ2[static_cast<size_t>(s)] < 0) {
      continue;
    }
    int u = s;
    int prev = -1;
    int step = 0;
    while (step < nv + 2) {
      if (vis[static_cast<size_t>(u)] && u == s && step > 0) {
        break;
      }
      vis[static_cast<size_t>(u)] = 1;
      const bool first = succ[static_cast<size_t>(u)] != prev;
      const int w = first ? succ[static_cast<size_t>(u)] : succ2[static_cast<size_t>(u)];
      const int id = first ? succE[static_cast<size_t>(u)] : succE2[static_cast<size_t>(u)];
      if (w < 0 || id < 0) {
        return false;
      }
      if ((step % 2) == 0) {
        const int L = u < nL ? u : w;
        const int R = u < nL ? w - nL : u - nL;
        if (L < 0 || L >= nL || R < 0 || R >= nR) {
          return false;
        }
        if (matchL[static_cast<size_t>(L)] < 0 && matchR[static_cast<size_t>(R)] < 0) {
          matchL[static_cast<size_t>(L)] = R;
          matchR[static_cast<size_t>(R)] = L;
          matchRing[static_cast<size_t>(L)] = adj[static_cast<size_t>(id)].ring;
          got++;
        }
      }
      prev = u;
      u = w;
      step++;
    }
  }
  return got == nL;
}

// Line graph of a bipartite graph on this component. Neighbourhoods split
// into the two endpoint cliques; a maximum set of rings is a maximum matching.
bool linePack(const std::vector<Row> &conf, const std::vector<int> &mem,
              const std::vector<int> &rankOf, std::vector<int> &out, int &matched,
              int &cover) {
  const int nAll = static_cast<int>(conf.size());
  const int m = static_cast<int>(mem.size());
  if (m == 0) {
    matched = 0;
    cover = 0;
    return true;
  }
  const bool timing = std::getenv("DSEAMS_CAGE_TIME") != nullptr;
  auto clockNow = std::chrono::steady_clock::now();
  auto mark = [&](const char *name) {
    if (!timing) {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    std::fprintf(stderr, "    %s %.3f\n", name,
                 std::chrono::duration<double, std::milli>(now - clockNow).count());
    clockNow = now;
  };
  auto linked = [&](int x, int y) { return conf[static_cast<size_t>(x)].has(y); };
  std::vector<PKey> keys;
  bool packed = true;
  bool all6 = true;
  for (const int v : mem) {
    if (conf[static_cast<size_t>(v)].n != 6) {
      all6 = false;
      break;
    }
  }
  if (all6) {
    keys.resize(static_cast<size_t>(m) * 2);
    std::vector<char> bad(static_cast<size_t>(m), 0);
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel for schedule(static) if (m >= 256)
#endif
    for (int i = 0; i < m; i++) {
      const int v = mem[static_cast<size_t>(i)];
      const Row &nb = conf[static_cast<size_t>(v)];
      auto near = [&](int x, int y) {
        const Row &row = conf[static_cast<size_t>(x)];
        for (int t = 0; t < row.n; t++) {
          if (row.v[t] == y) {
            return true;
          }
        }
        return false;
      };
      int A[3];
      int B[3];
      int nA = 1;
      int nB = 0;
      A[0] = nb.v[0];
      bool ok = true;
      for (int t = 1; t < 6; t++) {
        if (near(nb.v[0], nb.v[t])) {
          if (nA >= 3) {
            ok = false;
            break;
          }
          A[nA++] = nb.v[t];
        } else if (nB >= 3) {
          ok = false;
          break;
        } else {
          B[nB++] = nb.v[t];
        }
      }
      if (!ok || nA != 3 || nB != 3 || !near(A[1], A[2]) || !near(B[0], B[1]) ||
          !near(B[0], B[2]) || !near(B[1], B[2])) {
        bad[static_cast<size_t>(i)] = 1;
        continue;
      }
      for (int a = 0; a < 3 && ok; a++) {
        for (int b = 0; b < 3; b++) {
          if (near(A[a], B[b])) {
            ok = false;
            break;
          }
        }
      }
      if (!ok) {
        bad[static_cast<size_t>(i)] = 1;
        continue;
      }
      const int left[4] = {v, A[0], A[1], A[2]};
      const int right[4] = {v, B[0], B[1], B[2]};
      std::uint64_t k0 = packRanks(left, 4, rankOf);
      std::uint64_t k1 = packRanks(right, 4, rankOf);
      if (k0 == ~0ULL || k1 == ~0ULL || k0 == k1) {
        bad[static_cast<size_t>(i)] = 1;
        continue;
      }
      // Smaller key first, so clique ids follow ranks rather than storage order.
      if (k0 > k1) {
        const std::uint64_t tmp = k0;
        k0 = k1;
        k1 = tmp;
      }
      keys[static_cast<size_t>(2 * i)] = {k0, v, 0};
      keys[static_cast<size_t>(2 * i + 1)] = {k1, v, 1};
    }
    for (int i = 0; i < m; i++) {
      packed = packed && bad[static_cast<size_t>(i)] == 0;
    }
  } else {
  keys.reserve(static_cast<size_t>(m) * 2);
  for (const int v : mem) {
    const Row &nb = conf[static_cast<size_t>(v)];
    if (nb.n == 0 || nb.n > 6) {
      packed = false;
      break;
    }
    int A[6];
    int B[6];
    int nA = 0;
    int nB = 0;
    if (nb.n == 6) {
      nA = 1;
      A[0] = nb.v[0];
      bool ok = true;
      for (int i = 1; i < 6; i++) {
        if (linked(nb.v[0], nb.v[i])) {
          if (nA >= 3) {
            ok = false;
            break;
          }
          A[nA++] = nb.v[i];
        } else if (nB >= 3) {
          ok = false;
          break;
        } else {
          B[nB++] = nb.v[i];
        }
      }
      if (!ok || nA != 3 || nB != 3 || !linked(A[1], A[2]) || !linked(B[0], B[1]) ||
          !linked(B[0], B[2]) || !linked(B[1], B[2])) {
        packed = false;
        break;
      }
      for (int i = 0; i < 3 && packed; i++) {
        for (int j = 0; j < 3; j++) {
          if (linked(A[i], B[j])) {
            packed = false;
            break;
          }
        }
      }
      if (!packed) {
        break;
      }
    } else {
      // Neighbourhood splits into two cliques, found from the closed
      // neighbourhood of the first neighbour when the two sides are disjoint.
      std::uint32_t bit[8] = {};
      for (int a = 0; a < nb.n; a++) {
        for (int b = a + 1; b < nb.n; b++) {
          if (linked(nb.v[a], nb.v[b])) {
            bit[a] |= 1u << b;
            bit[b] |= 1u << a;
          }
        }
      }
      const unsigned all = (1u << nb.n) - 1u;
      const unsigned mask = bit[0] | 1u;
      const unsigned other = all ^ mask;
      auto closed = [&](unsigned side) {
        if (side == 0) {
          return false;
        }
        for (int a = 0; a < nb.n; a++) {
          if (((side >> a) & 1u) == 0) {
            continue;
          }
          if (bit[a] != (side ^ (1u << a))) {
            return false;
          }
        }
        return true;
      };
      auto clique = [&](unsigned side) {
        for (int a = 0; a < nb.n; a++) {
          if (((side >> a) & 1u) == 0) {
            continue;
          }
          for (int b = a + 1; b < nb.n; b++) {
            if (((side >> b) & 1u) && (bit[a] & (1u << b)) == 0) {
              return false;
            }
          }
        }
        return true;
      };
      unsigned sideA = 0;
      unsigned sideB = 0;
      if (other != 0 && closed(mask) && closed(other)) {
        sideA = mask;
        sideB = other;
      } else if (clique(all)) {
        sideA = all;
        sideB = 0;
      } else if (nb.n <= 10) {
        bool found = false;
        for (unsigned s = 1; s < all; s++) {
          const unsigned o = all ^ s;
          if (o < s) {
            continue;
          }
          if (clique(s) && clique(o)) {
            sideA = s;
            sideB = o;
            found = true;
            break;
          }
        }
        if (!found) {
          packed = false;
          break;
        }
      } else {
        packed = false;
        break;
      }
      nA = 0;
      nB = 0;
      for (int a = 0; a < nb.n; a++) {
        if ((sideA >> a) & 1u) {
          A[nA++] = nb.v[a];
        } else if ((sideB >> a) & 1u) {
          B[nB++] = nb.v[a];
        }
      }
    }
    int left[8];
    int right[8];
    int nL = 0;
    int nR = 0;
    left[nL++] = v;
    for (int i = 0; i < nA; i++) {
      if (nL >= 4) {
        packed = false;
        break;
      }
      left[nL++] = A[i];
    }
    right[nR++] = v;
    for (int i = 0; i < nB; i++) {
      if (nR >= 4) {
        packed = false;
        break;
      }
      right[nR++] = B[i];
    }
    if (!packed) {
      break;
    }
    std::uint64_t k0 = packRanks(left, nL, rankOf);
    std::uint64_t k1 = packRanks(right, nR, rankOf);
    if (k0 == ~0ULL || k1 == ~0ULL || k0 == k1) {
      packed = false;
      break;
    }
    if (k0 > k1) {
      const std::uint64_t tmp = k0;
      k0 = k1;
      k1 = tmp;
    }
    keys.push_back({k0, v, 0});
    keys.push_back({k1, v, 1});
  }
  }
  if (!packed || keys.size() != static_cast<size_t>(m) * 2) {
    return false;
  }
  mark("split");
  // Clique ids follow the first time each key appears in canonical ring order,
  // so the numbering does not depend on a sort of the keys.
  std::vector<std::array<int, 2>> ends(static_cast<size_t>(nAll), std::array<int, 2>{-1, -1});
  const size_t nTab = 1u << 14;
  struct Ent {
    std::uint64_t k;
    int id;
  };
  std::vector<Ent> tab(nTab, Ent{0, -1});
  std::vector<int> cliqueCount;
  int nC = 0;
  for (const PKey &key : keys) {
    size_t h = static_cast<size_t>(key.k * 0x9E3779B185EBCA87ULL) & (nTab - 1);
    while (tab[h].id >= 0 && tab[h].k != key.k) {
      h = (h + 1) & (nTab - 1);
    }
    if (tab[h].id < 0) {
      tab[h] = Ent{key.k, nC++};
      cliqueCount.push_back(0);
    }
    cliqueCount[static_cast<size_t>(tab[h].id)]++;
    ends[static_cast<size_t>(key.ring)][static_cast<size_t>(key.slot)] = tab[h].id;
  }
  bool groupsOk = all6;
  if (all6) {
    for (const int c : cliqueCount) {
      if (c != 4) {
        groupsOk = false;
        break;
      }
    }
  } else {
    groupsOk = false;
  }
  mark("radix");
  std::vector<int> cOff(static_cast<size_t>(nC + 1), 0);
  std::vector<int> gOff(static_cast<size_t>(nC + 1), 0);
  for (const int v : mem) {
    const int a = ends[static_cast<size_t>(v)][0];
    const int b = ends[static_cast<size_t>(v)][1];
    if (a < 0 || b < 0 || a == b) {
      return false;
    }
    cOff[static_cast<size_t>(a + 1)]++;
    cOff[static_cast<size_t>(b + 1)]++;
    gOff[static_cast<size_t>(a + 1)]++;
    gOff[static_cast<size_t>(b + 1)]++;
  }
  for (int c = 0; c < nC; c++) {
    cOff[static_cast<size_t>(c + 1)] += cOff[static_cast<size_t>(c)];
    gOff[static_cast<size_t>(c + 1)] += gOff[static_cast<size_t>(c)];
  }
  std::vector<int> cMem(static_cast<size_t>(cOff[static_cast<size_t>(nC)]));
  std::vector<int> gTo(static_cast<size_t>(gOff[static_cast<size_t>(nC)]));
  std::vector<int> cPut(cOff.begin(), cOff.end() - 1);
  std::vector<int> gPut(gOff.begin(), gOff.end() - 1);
  for (const int v : mem) {
    const int a = ends[static_cast<size_t>(v)][0];
    const int b = ends[static_cast<size_t>(v)][1];
    cMem[static_cast<size_t>(cPut[static_cast<size_t>(a)]++)] = v;
    cMem[static_cast<size_t>(cPut[static_cast<size_t>(b)]++)] = v;
    gTo[static_cast<size_t>(gPut[static_cast<size_t>(a)]++)] = b;
    gTo[static_cast<size_t>(gPut[static_cast<size_t>(b)]++)] = a;
  }
  mark("csr");
  if (!groupsOk) {
  std::vector<int> stampV(static_cast<size_t>(nAll), 0);
  int era = 1;
  for (const int v : mem) {
    const int a = ends[static_cast<size_t>(v)][0];
    const int b = ends[static_cast<size_t>(v)][1];
    if (++era == 0) {
      std::fill(stampV.begin(), stampV.end(), 0);
      era = 1;
    }
    int ne = 0;
    for (const int c : {a, b}) {
      for (int p = cOff[static_cast<size_t>(c)]; p < cOff[static_cast<size_t>(c + 1)]; p++) {
        const int u = cMem[static_cast<size_t>(p)];
        if (u != v && stampV[static_cast<size_t>(u)] != era) {
          stampV[static_cast<size_t>(u)] = era;
          ne++;
        }
      }
    }
    const Row &row = conf[static_cast<size_t>(v)];
    if (ne != row.n) {
      return false;
    }
    for (int i = 0; i < row.n; i++) {
      if (stampV[static_cast<size_t>(row.v[i])] != era) {
        return false;
      }
    }
  }
  }
  mark("equal");
  std::vector<int> color(static_cast<size_t>(nC), -1);
  std::vector<int> q(static_cast<size_t>(nC));
  for (int s = 0; s < nC; s++) {
    if (color[static_cast<size_t>(s)] >= 0) {
      continue;
    }
    color[static_cast<size_t>(s)] = 0;
    int qh = 0;
    int qt = 0;
    q[static_cast<size_t>(qt++)] = s;
    while (qh < qt) {
      const int u = q[static_cast<size_t>(qh++)];
      for (int e = gOff[static_cast<size_t>(u)]; e < gOff[static_cast<size_t>(u + 1)]; e++) {
        const int w = gTo[static_cast<size_t>(e)];
        if (color[static_cast<size_t>(w)] < 0) {
          color[static_cast<size_t>(w)] = 1 - color[static_cast<size_t>(u)];
          q[static_cast<size_t>(qt++)] = w;
        } else if (color[static_cast<size_t>(w)] == color[static_cast<size_t>(u)]) {
          return false;
        }
      }
    }
  }
  std::vector<int> leftOf(static_cast<size_t>(nC), -1);
  std::vector<int> rightOf(static_cast<size_t>(nC), -1);
  int nL = 0;
  int nR = 0;
  for (int c = 0; c < nC; c++) {
    if (color[static_cast<size_t>(c)] == 0) {
      leftOf[static_cast<size_t>(c)] = nL++;
    } else {
      rightOf[static_cast<size_t>(c)] = nR++;
    }
  }
  std::vector<int> off(static_cast<size_t>(nL + 1), 0);
  std::vector<AdjE> adj(static_cast<size_t>(m));
  for (const int v : mem) {
    const int a = ends[static_cast<size_t>(v)][0];
    const int b = ends[static_cast<size_t>(v)][1];
    int L = leftOf[static_cast<size_t>(a)];
    int R = rightOf[static_cast<size_t>(b)];
    if (L < 0 || R < 0) {
      L = leftOf[static_cast<size_t>(b)];
      R = rightOf[static_cast<size_t>(a)];
    }
    if (L < 0 || R < 0) {
      return false;
    }
    off[static_cast<size_t>(L + 1)]++;
    (void)R;
  }
  // Second pass writes the adjacency in ring-rank order inside each left vertex.
  for (int u = 0; u < nL; u++) {
    off[static_cast<size_t>(u + 1)] += off[static_cast<size_t>(u)];
  }
  std::vector<int> put = off;
  for (const int v : mem) {
    const int a = ends[static_cast<size_t>(v)][0];
    const int b = ends[static_cast<size_t>(v)][1];
    int L = leftOf[static_cast<size_t>(a)];
    int R = rightOf[static_cast<size_t>(b)];
    if (L < 0 || R < 0) {
      L = leftOf[static_cast<size_t>(b)];
      R = rightOf[static_cast<size_t>(a)];
    }
    adj[static_cast<size_t>(put[static_cast<size_t>(L)]++)] = {R, v};
  }
  for (int u = 0; u < nL; u++) {
    std::sort(adj.begin() + off[static_cast<size_t>(u)],
              adj.begin() + off[static_cast<size_t>(u + 1)],
              [&](const AdjE &x, const AdjE &y) {
                return rankOf[static_cast<size_t>(x.ring)] < rankOf[static_cast<size_t>(y.ring)];
              });
  }
  std::vector<int> matchL(static_cast<size_t>(nL), -1);
  std::vector<int> matchR(static_cast<size_t>(nR), -1);
  std::vector<int> matchRing(static_cast<size_t>(nL), -1);
  const bool usedEuler = eulerPerfect(nL, nR, off, adj, matchL, matchR, matchRing);
  mark("euler");
  if (!usedEuler) {
  std::vector<int> byRank = mem;
  std::sort(byRank.begin(), byRank.end(), [&](int a, int b) {
    return rankOf[static_cast<size_t>(a)] < rankOf[static_cast<size_t>(b)];
  });
  mark("color");
  int pre = 0;
  for (const int v : byRank) {
    const int a = ends[static_cast<size_t>(v)][0];
    const int b = ends[static_cast<size_t>(v)][1];
    int L = leftOf[static_cast<size_t>(a)];
    int R = rightOf[static_cast<size_t>(b)];
    if (L < 0 || R < 0) {
      L = leftOf[static_cast<size_t>(b)];
      R = rightOf[static_cast<size_t>(a)];
    }
    if (matchL[static_cast<size_t>(L)] < 0 && matchR[static_cast<size_t>(R)] < 0) {
      matchL[static_cast<size_t>(L)] = R;
      matchR[static_cast<size_t>(R)] = L;
      matchRing[static_cast<size_t>(L)] = v;
      pre++;
    }
  }
  if (timing) {
    std::fprintf(stderr, "    greedy %d / %d\n", pre, nL);
  }
  mark("greedy");
  std::vector<int> seen(static_cast<size_t>(nR), 0);
  int stamp = 0;
  struct Step {
    int u;
    int e;
    int viaR;
    int viaRing;
  };
  std::vector<Step> st;
  st.reserve(static_cast<size_t>(nL));
  for (int start = 0; start < nL; start++) {
    if (matchL[static_cast<size_t>(start)] >= 0) {
      continue;
    }
    if (++stamp == 0) {
      std::fill(seen.begin(), seen.end(), 0);
      stamp = 1;
    }
    st.clear();
    st.push_back({start, off[static_cast<size_t>(start)], -1, -1});
    bool found = false;
    while (!st.empty()) {
      Step &fr = st.back();
      if (fr.e >= off[static_cast<size_t>(fr.u + 1)]) {
        st.pop_back();
        continue;
      }
      const AdjE edge = adj[static_cast<size_t>(fr.e++)];
      const int r = edge.right;
      if (seen[static_cast<size_t>(r)] == stamp) {
        continue;
      }
      seen[static_cast<size_t>(r)] = stamp;
      const int w = matchR[static_cast<size_t>(r)];
      if (w < 0) {
        int curR = r;
        int curRing = edge.ring;
        for (int k = static_cast<int>(st.size()) - 1; k >= 0; --k) {
          const int left = st[static_cast<size_t>(k)].u;
          const int oldR = matchL[static_cast<size_t>(left)];
          matchL[static_cast<size_t>(left)] = curR;
          matchR[static_cast<size_t>(curR)] = left;
          matchRing[static_cast<size_t>(left)] = curRing;
          if (oldR < 0) {
            break;
          }
          curR = st[static_cast<size_t>(k)].viaR;
          curRing = st[static_cast<size_t>(k)].viaRing;
        }
        found = true;
        break;
      }
      st.push_back({w, off[static_cast<size_t>(w)], r, edge.ring});
    }
    (void)found;
  }
  }
  mark("kuhn");
  // König cover from the free left vertices. An augmenting path would mean
  // the search above stopped short. A perfect matching of the left shore is
  // already a cover of that shore.
  int nFree = 0;
  for (int u = 0; u < nL; u++) {
    nFree += matchL[static_cast<size_t>(u)] < 0 ? 1 : 0;
  }
  std::vector<char> zL(static_cast<size_t>(nL), 0);
  std::vector<char> zR(static_cast<size_t>(nR), 0);
  int qh = 0;
  int qt = 0;
  bool augmenting = false;
  if (nFree != 0) {
  for (int u = 0; u < nL; u++) {
    if (matchL[static_cast<size_t>(u)] < 0) {
      zL[static_cast<size_t>(u)] = 1;
      q[static_cast<size_t>(qt++)] = u;
    }
  }
  // q was sized nC, which is nL+nR. Free-left BFS stays inside that.
  while (qh < qt) {
    const int u = q[static_cast<size_t>(qh++)];
    for (int e = off[static_cast<size_t>(u)]; e < off[static_cast<size_t>(u + 1)]; e++) {
      const int r = adj[static_cast<size_t>(e)].right;
      if (r == matchL[static_cast<size_t>(u)] || zR[static_cast<size_t>(r)]) {
        continue;
      }
      zR[static_cast<size_t>(r)] = 1;
      const int w = matchR[static_cast<size_t>(r)];
      if (w < 0) {
        augmenting = true;
      } else if (!zL[static_cast<size_t>(w)]) {
        zL[static_cast<size_t>(w)] = 1;
        q[static_cast<size_t>(qt++)] = w;
      }
    }
  }
  }
  if (augmenting || qt > nC) {
    return false;
  }
  int nMatch = 0;
  int nCover = 0;
  for (int u = 0; u < nL; u++) {
    nMatch += matchL[static_cast<size_t>(u)] >= 0 ? 1 : 0;
    nCover += zL[static_cast<size_t>(u)] ? 0 : 1;
  }
  for (int r = 0; r < nR; r++) {
    nCover += zR[static_cast<size_t>(r)] ? 1 : 0;
  }
  if (nMatch != nCover) {
    return false;
  }
  out.clear();
  out.reserve(static_cast<size_t>(nMatch));
  for (int u = 0; u < nL; u++) {
    if (matchRing[static_cast<size_t>(u)] >= 0) {
      out.push_back(matchRing[static_cast<size_t>(u)]);
    }
  }
  matched = nMatch;
  cover = nCover;
  return true;
}

} // namespace

int bruteMis(const std::vector<Row> &conf, const std::vector<int> &mem) {
  const int n = static_cast<int>(mem.size());
  int best = 0;
  const int limit = 1 << n;
  for (int mask = 0; mask < limit; mask++) {
    bool ok = true;
    int cnt = 0;
    for (int i = 0; i < n && ok; i++) {
      if (((mask >> i) & 1) == 0) {
        continue;
      }
      cnt++;
      const Row &row = conf[static_cast<size_t>(mem[static_cast<size_t>(i)])];
      for (int k = 0; k < row.n; k++) {
        for (int j = 0; j < n; j++) {
          if (mem[static_cast<size_t>(j)] == row.v[k] && ((mask >> j) & 1)) {
            ok = false;
          }
        }
      }
    }
    if (ok) {
      best = std::max(best, cnt);
    }
  }
  return best;
}

void misSelfTest() {
  auto check = [](const std::vector<std::pair<int, int>> &edges, int n, int expect) {
    std::vector<Row> conf(static_cast<size_t>(n));
    for (const auto &e : edges) {
      conf[static_cast<size_t>(e.first)].push(e.second);
      conf[static_cast<size_t>(e.second)].push(e.first);
    }
    std::vector<int> mem(static_cast<size_t>(n));
    for (int i = 0; i < n; i++) {
      mem[static_cast<size_t>(i)] = i;
    }
    const MisOut got = exactMis(conf, mem);
    const int brute = bruteMis(conf, mem);
    if (!got.solved || static_cast<int>(got.ids.size()) != brute || brute != expect ||
        got.cover < brute) {
      std::fprintf(stderr, "mis-selftest FAIL n %d got %zu brute %d expect %d cover %d\n", n,
                   got.ids.size(), brute, expect, got.cover);
      std::exit(1);
    }
  };
  check({}, 1, 1);
  check({{0, 1}}, 2, 1);
  check({{0, 1}, {1, 2}, {2, 0}}, 3, 1);
  check({{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, 5, 2);
  check({{0, 1}, {2, 3}, {4, 5}}, 6, 3);
  std::fprintf(stderr, "mis-selftest ok\n");
}

std::vector<int> ring::maximumEquators(const std::vector<ring::EquatorChoice> &periph,
                                       const std::vector<int> &order,
                                       ring::PackingCertificate *certificate) {
  static int selfTested = 0;
  if (selfTested == 0 && std::getenv("DSEAMS_MIS_SELFTEST") != nullptr) {
    selfTested = 1;
    misSelfTest();
  }
  const auto t0 = std::chrono::steady_clock::now();
  const int n = static_cast<int>(periph.size());
  std::vector<Row> conf(static_cast<size_t>(n));
  int nEdges = 0;
  int nCand = 0;
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel for schedule(static) reduction(+ : nEdges, nCand) if (n >= 256)
#endif
  for (int i = 0; i < n; i++) {
    if (periph[static_cast<size_t>(i)].n == 0) {
      continue;
    }
    nCand++;
    const ring::EquatorChoice &row = periph[static_cast<size_t>(i)];
    for (int k = 0; k < row.n; k++) {
      const int f = row.ring[k];
      if (f < 0 || f >= n || f == i || periph[static_cast<size_t>(f)].n == 0) {
        continue;
      }
      bool back = false;
      const ring::EquatorChoice &other = periph[static_cast<size_t>(f)];
      for (int t = 0; t < other.n; t++) {
        if (other.ring[t] == i) {
          back = true;
          break;
        }
      }
      if (!back) {
        continue;
      }
      conf[static_cast<size_t>(i)].push(f);
    }
    conf[static_cast<size_t>(i)].sortUnique();
    nEdges += conf[static_cast<size_t>(i)].n;
  }
  nEdges /= 2;
  std::vector<int> rankOf(static_cast<size_t>(n));
  for (int r = 0; r < static_cast<int>(order.size()) && r < n; r++) {
    rankOf[static_cast<size_t>(order[static_cast<size_t>(r)])] = r;
  }
  const auto tConf = std::chrono::steady_clock::now();
  const auto tComp = tConf;

  std::vector<int> seen(static_cast<size_t>(n), 0);
  std::vector<int> take;
  int nComp = 0;
  int nIso = 0;
  int matchSum = 0;
  int coverSum = 0;
  bool anySearch = false;
  bool allSolved = true;
  struct Comp {
    int minRank;
    std::vector<int> mem;
  };
  std::vector<Comp> comps;
  std::vector<int> stack;
  for (int s = 0; s < n; s++) {
    if (periph[static_cast<size_t>(s)].n == 0 || seen[static_cast<size_t>(s)]) {
      continue;
    }
    stack.clear();
    stack.push_back(s);
    seen[static_cast<size_t>(s)] = 1;
    Comp comp;
    comp.minRank = rankOf[static_cast<size_t>(s)];
    while (!stack.empty()) {
      const int u = stack.back();
      stack.pop_back();
      comp.mem.push_back(u);
      comp.minRank = std::min(comp.minRank, rankOf[static_cast<size_t>(u)]);
      const Row &row = conf[static_cast<size_t>(u)];
      for (int k = 0; k < row.n; k++) {
        const int w = row.v[k];
        if (!seen[static_cast<size_t>(w)]) {
          seen[static_cast<size_t>(w)] = 1;
          stack.push_back(w);
        }
      }
    }
    nComp++;
    if (comp.mem.size() == 1) {
      take.push_back(comp.mem[0]);
      nIso++;
      matchSum++;
      coverSum++;
      continue;
    }
    comps.push_back(std::move(comp));
  }
  std::sort(comps.begin(), comps.end(), [](const Comp &a, const Comp &b) {
    if (a.minRank != b.minRank) {
      return a.minRank < b.minRank;
    }
    return a.mem[0] < b.mem[0];
  });
  for (Comp &comp : comps) {
    std::sort(comp.mem.begin(), comp.mem.end(), [&](int a, int b) {
      return rankOf[static_cast<size_t>(a)] < rankOf[static_cast<size_t>(b)];
    });
    std::vector<int> picked;
    int partMatch = 0;
    int partCover = 0;
    if (linePack(conf, comp.mem, rankOf, picked, partMatch, partCover)) {
      take.insert(take.end(), picked.begin(), picked.end());
      matchSum += partMatch;
      coverSum += partCover;
      continue;
    }
    anySearch = true;
    const MisOut solved = exactMis(conf, comp.mem);
    if (!solved.solved) {
      allSolved = false;
      continue;
    }
    take.insert(take.end(), solved.ids.begin(), solved.ids.end());
    matchSum += static_cast<int>(solved.ids.size());
    coverSum += solved.cover;
  }
  std::sort(take.begin(), take.end(), [&](int a, int b) {
    return rankOf[static_cast<size_t>(a)] < rankOf[static_cast<size_t>(b)];
  });
  if (certificate != nullptr) {
    const int cages = static_cast<int>(take.size());
    certificate->cages = cages;
    certificate->matching = matchSum;
    certificate->cover = coverSum;
    certificate->upper = allSolved ? cages : coverSum;
    certificate->isolates = nIso;
    certificate->vertices = nCand;
    certificate->edges = nEdges;
    certificate->components = nComp;
    certificate->method = anySearch ? 2 : 1;
    certificate->tight = allSolved && coverSum == cages;
    certificate->milliseconds =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  }
  if (std::getenv("DSEAMS_CAGE_TIME")) {
    const auto now = std::chrono::steady_clock::now();
    std::fprintf(stderr, "  pack conf %.3f solve %.3f kept %zu comps %d\n",
                 std::chrono::duration<double, std::milli>(tConf - t0).count(),
                 std::chrono::duration<double, std::milli>(now - tComp).count(), take.size(),
                 nComp);
  }
  return take;
}
