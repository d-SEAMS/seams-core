#include <par_team.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <numeric>
#include <unordered_map>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

std::uint64_t gDdcMatchSteps = 0;
std::uint64_t gDdcRowCopies = 0;
std::uint64_t gDdcKeyDigits = 0;
std::uint64_t gDdcBucketSteps = 0;
thread_local std::vector<std::uint64_t> *gPackedLive = nullptr;

// One conflict row per candidate. Eight neighbours sit in the row, which is
// the whole row on ice (degree six) and matches the dense layout of that
// search. A wider neighbourhood spills into a side table so the inline row
// stays a plain block of ints.
struct SpillBlock {
  int *p = nullptr;
  int cap = 0;
};
std::mutex gRowSpillMu;
std::unordered_map<const struct Row *, SpillBlock> gRowSpill;

struct Row {
  static constexpr int kInline = 8;
  int head[kInline]{};
  int n = 0;

  [[nodiscard]] int operator[](int i) const {
    if (i < kInline) {
      return head[i];
    }
    std::lock_guard<std::mutex> lock(gRowSpillMu);
    return gRowSpill.find(this)->second.p[i - kInline];
  }

  bool push(int x) {
    if (n < kInline) {
      head[n++] = x;
      return true;
    }
    std::lock_guard<std::mutex> lock(gRowSpillMu);
    SpillBlock &sp = gRowSpill[this];
    const int at = n - kInline;
    if (at == sp.cap) {
      const int cap = sp.cap == 0 ? 4 : sp.cap * 2;
      int *grown = new int[static_cast<std::size_t>(cap)];
      if (sp.cap > 0) {
        std::copy(sp.p, sp.p + sp.cap, grown);
      }
      delete[] sp.p;
      sp.p = grown;
      sp.cap = cap;
    }
    sp.p[at] = x;
    ++n;
    return true;
  }

  void sortUnique() {
    if (n <= 1) {
      return;
    }
    if (n <= kInline) {
      for (int i = 1; i < n; i++) {
        const int value = head[i];
        int j = i;
        while (j > 0 && head[j - 1] > value) {
          head[j] = head[j - 1];
          --j;
        }
        head[j] = value;
      }
      int kept = 1;
      for (int i = 1; i < n; i++) {
        if (head[i] != head[kept - 1]) {
          head[kept++] = head[i];
        }
      }
      n = kept;
      return;
    }
    gDdcRowCopies += static_cast<std::uint64_t>(n);
    int stack[64];
    int *buf = stack;
    std::vector<int> heap;
    if (n > 64) {
      heap.resize(static_cast<std::size_t>(n));
      buf = heap.data();
    }
    const int filled = n;
    for (int i = 0; i < filled; i++) {
      buf[i] = (*this)[i];
    }
    std::sort(buf, buf + filled);
    const int kept = static_cast<int>(std::unique(buf, buf + filled) - buf);
    {
      std::lock_guard<std::mutex> lock(gRowSpillMu);
      const auto it = gRowSpill.find(this);
      if (it != gRowSpill.end()) {
        delete[] it->second.p;
        gRowSpill.erase(it);
      }
    }
    n = 0;
    for (int i = 0; i < kept; i++) {
      push(buf[i]);
    }
  }

  [[nodiscard]] bool has(int x) const {
    for (int i = 0; i < n && i < kInline; i++) {
      if (head[i] == x) {
        return true;
      }
    }
    for (int i = kInline; i < n; i++) {
      if ((*this)[i] == x) {
        return true;
      }
    }
    return false;
  }
};

void releaseRowSpills() {
  std::lock_guard<std::mutex> lock(gRowSpillMu);
  for (auto &kv : gRowSpill) {
    delete[] kv.second.p;
  }
  gRowSpill.clear();
}

struct RowSpillGuard {
  ~RowSpillGuard() { releaseRowSpills(); }
};

// One bit per vertex of the component under search. The word count is
// the component order, so a larger crystal is not discarded.
struct Bits {
  std::vector<std::uint64_t> w;

  Bits() = default;
  explicit Bits(int nbits) {
    if (nbits > 0) {
      w.assign(static_cast<std::size_t>((nbits + 63) >> 6), 0);
    }
  }

  void set(int i) { w[static_cast<std::size_t>(i >> 6)] |= 1ULL << (i & 63); }
  void reset(int i) { w[static_cast<std::size_t>(i >> 6)] &= ~(1ULL << (i & 63)); }
  [[nodiscard]] bool test(int i) const {
    return (w[static_cast<std::size_t>(i >> 6)] >> (i & 63)) & 1ULL;
  }
  [[nodiscard]] int count() const {
    int n = 0;
    for (const std::uint64_t word : w) {
      n += __builtin_popcountll(word);
    }
    return n;
  }
  [[nodiscard]] bool empty() const {
    for (const std::uint64_t word : w) {
      if (word) {
        return false;
      }
    }
    return true;
  }
  Bits operator&(const Bits &o) const {
    Bits r;
    r.w.resize(w.size());
    for (std::size_t i = 0; i < w.size(); i++) {
      r.w[i] = w[i] & o.w[i];
    }
    return r;
  }
  Bits operator-(const Bits &o) const {
    Bits r;
    r.w.resize(w.size());
    for (std::size_t i = 0; i < w.size(); i++) {
      r.w[i] = w[i] & ~o.w[i];
    }
    return r;
  }
  [[nodiscard]] int first() const {
    for (std::size_t i = 0; i < w.size(); i++) {
      if (w[i]) {
        return static_cast<int>(i) * 64 + __builtin_ctzll(w[i]);
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
  std::vector<int> at(conf.size(), -1);
  for (int i = 0; i < n; i++) {
    at[static_cast<size_t>(mem[static_cast<size_t>(i)])] = i;
  }
  std::vector<Bits> nbr(static_cast<size_t>(n), Bits(n));
  for (int i = 0; i < n; i++) {
    nbr[static_cast<size_t>(i)].set(i);
    const int g = mem[static_cast<size_t>(i)];
    const Row &row = conf[static_cast<size_t>(g)];
    for (int k = 0; k < row.n; k++) {
      const int j = at[static_cast<size_t>(row[k])];
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
  Bits all(n);
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

// Sorted canonical ranks of one clique. The integers are the ranks
// themselves, so a frame with more than 2^16 rings still names its cliques.
struct RankKey {
  int id[4]{};
  int n = 0;
  bool ok = false;

  [[nodiscard]] bool operator==(const RankKey &o) const {
    if (!ok || !o.ok || n != o.n) {
      return false;
    }
    for (int i = 0; i < n; i++) {
      if (id[i] != o.id[i]) {
        return false;
      }
    }
    return true;
  }
  [[nodiscard]] bool operator<(const RankKey &o) const {
    const int m = n < o.n ? n : o.n;
    for (int i = 0; i < m; i++) {
      if (id[i] != o.id[i]) {
        return id[i] < o.id[i];
      }
    }
    return n < o.n;
  }
};

struct PKey {
  RankKey k;
  int ring;
  int slot;
};

RankKey packRanks(const int *id, int n, const std::vector<int> &rankOf) {
  RankKey key;
  if (n <= 0 || n > 4) {
    return key;
  }
  int mem[4];
  for (int i = 0; i < n; i++) {
    mem[i] = rankOf[static_cast<size_t>(id[i])];
    if (mem[i] < 0) {
      return key;
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
  key.n = n;
  for (int i = 0; i < n; i++) {
    key.id[i] = mem[i];
  }
  key.ok = true;
  return key;
}

struct AdjE {
  int right;
  int ring;
};

struct HalfEdge {
  int left;
  int right;
  int ring;
};

// Pair the edges at each vertex, then walk each alternating cycle once.
// Each split halves a positive even degree. The last split is a perfect
// matching. The representative of a cycle is its lowest arc, so the kept
// edges are the ones the pointer-jumping split kept, without the log passes.
struct EulerScratch {
  std::vector<HalfEdge> edges;
  std::vector<HalfEdge> next;
  std::vector<int> leftCount, rightCount, leftOff, rightOff;
  std::vector<int> leftAt, rightAt, leftPut, rightPut;
  std::vector<int> partnerL, partnerR, succ;
  std::vector<int> keptAt, keptRight, rep, arcRank, cycle;
  std::vector<char> seen, keep;
  std::vector<int> gotL, gotR, gotRing;

  template <class T> static void fit(std::vector<T> &v, int n) {
    if (static_cast<int>(v.size()) < n) {
      v.resize(static_cast<std::size_t>(n));
    }
  }
};

bool parallelEulerMatch(int nL, int nR, const std::vector<int> &off, const std::vector<AdjE> &adj,
                        std::vector<int> &matchL, std::vector<int> &matchR,
                        std::vector<int> &matchRing) {
  // Pin the scratch. Reloading the thread_local address on each field was
  // a million instructions at a few thousand molecules.
  thread_local EulerScratch storage;
  EulerScratch *sc = &storage;
  asm volatile("" : "+r"(sc));
  const int m0 = static_cast<int>(adj.size());
  if (nL <= 0 || nL != nR || m0 <= 0) {
    return false;
  }
  EulerScratch::fit(sc->edges, m0);
  for (int u = 0; u < nL; u++) {
    for (int e = off[static_cast<std::size_t>(u)]; e < off[static_cast<std::size_t>(u + 1)]; e++) {
      sc->edges[static_cast<std::size_t>(e)] = {u, adj[static_cast<std::size_t>(e)].right,
                                               adj[static_cast<std::size_t>(e)].ring};
    }
  }
  int m = m0;
  int degree = off[1] - off[0];
  while (degree > 1) {
    if (degree <= 0 || (degree % 2) != 0) {
      return false;
    }
    EulerScratch::fit(sc->leftCount, nL);
    EulerScratch::fit(sc->rightCount, nR);
    std::fill(sc->leftCount.begin(), sc->leftCount.begin() + nL, 0);
    std::fill(sc->rightCount.begin(), sc->rightCount.begin() + nR, 0);
    for (int e = 0; e < m; e++) {
      const HalfEdge &edge = sc->edges[static_cast<std::size_t>(e)];
      if (edge.left < 0 || edge.left >= nL || edge.right < 0 || edge.right >= nR) {
        return false;
      }
      sc->leftCount[static_cast<std::size_t>(edge.left)]++;
      sc->rightCount[static_cast<std::size_t>(edge.right)]++;
    }
    for (int u = 0; u < nL; u++) {
      if (sc->leftCount[static_cast<std::size_t>(u)] != degree) {
        return false;
      }
    }
    for (int r = 0; r < nR; r++) {
      if (sc->rightCount[static_cast<std::size_t>(r)] != degree) {
        return false;
      }
    }
    EulerScratch::fit(sc->leftOff, nL + 1);
    EulerScratch::fit(sc->rightOff, nR + 1);
    sc->leftOff[0] = 0;
    sc->rightOff[0] = 0;
    for (int u = 0; u < nL; u++) {
      sc->leftOff[static_cast<std::size_t>(u + 1)] =
          sc->leftOff[static_cast<std::size_t>(u)] + sc->leftCount[static_cast<std::size_t>(u)];
    }
    for (int r = 0; r < nR; r++) {
      sc->rightOff[static_cast<std::size_t>(r + 1)] =
          sc->rightOff[static_cast<std::size_t>(r)] + sc->rightCount[static_cast<std::size_t>(r)];
    }
    EulerScratch::fit(sc->leftAt, m);
    EulerScratch::fit(sc->rightAt, m);
    EulerScratch::fit(sc->leftPut, nL + 1);
    EulerScratch::fit(sc->rightPut, nR + 1);
    std::copy(sc->leftOff.begin(), sc->leftOff.begin() + nL + 1, sc->leftPut.begin());
    std::copy(sc->rightOff.begin(), sc->rightOff.begin() + nR + 1, sc->rightPut.begin());
    for (int e = 0; e < m; e++) {
      sc->leftAt[static_cast<std::size_t>(
          sc->leftPut[static_cast<std::size_t>(sc->edges[static_cast<std::size_t>(e)].left)]++)] = e;
      sc->rightAt[static_cast<std::size_t>(
          sc->rightPut[static_cast<std::size_t>(sc->edges[static_cast<std::size_t>(e)].right)]++)] = e;
    }
    auto byRing = [&](int a, int b) {
      const int ra = sc->edges[static_cast<std::size_t>(a)].ring;
      const int rb = sc->edges[static_cast<std::size_t>(b)].ring;
      if (ra != rb) {
        return ra < rb;
      }
      return a < b;
    };
    // Degree is a handful (four on ice). Insertion sort of that run is the
    // same order as std::sort and does not pay for an introsort frame.
    auto sortRun = [&](std::vector<int> &at, int lo, int hi) {
      if (hi - lo > 16) {
        std::sort(at.begin() + lo, at.begin() + hi, byRing);
        return;
      }
      for (int i = lo + 1; i < hi; i++) {
        const int v = at[static_cast<std::size_t>(i)];
        int j = i;
        while (j > lo && byRing(v, at[static_cast<std::size_t>(j - 1)])) {
          at[static_cast<std::size_t>(j)] = at[static_cast<std::size_t>(j - 1)];
          --j;
        }
        at[static_cast<std::size_t>(j)] = v;
      }
    };
    for (int u = 0; u < nL; u++) {
      sortRun(sc->leftAt, sc->leftOff[static_cast<std::size_t>(u)],
              sc->leftOff[static_cast<std::size_t>(u + 1)]);
    }
    for (int r = 0; r < nR; r++) {
      sortRun(sc->rightAt, sc->rightOff[static_cast<std::size_t>(r)],
              sc->rightOff[static_cast<std::size_t>(r + 1)]);
    }
    EulerScratch::fit(sc->partnerL, m);
    EulerScratch::fit(sc->partnerR, m);
    auto pairShore = [&](int nShore, const std::vector<int> &shoreOff, const std::vector<int> &shoreAt,
                         std::vector<int> &partner) {
      for (int u = 0; u < nShore; u++) {
        const int lo = shoreOff[static_cast<std::size_t>(u)];
        const int hi = shoreOff[static_cast<std::size_t>(u + 1)];
        if (((hi - lo) % 2) != 0) {
          return false;
        }
        for (int p = lo; p < hi; p += 2) {
          const int a = shoreAt[static_cast<std::size_t>(p)];
          const int b = shoreAt[static_cast<std::size_t>(p + 1)];
          partner[static_cast<std::size_t>(a)] = b;
          partner[static_cast<std::size_t>(b)] = a;
        }
      }
      return true;
    };
    if (!pairShore(nL, sc->leftOff, sc->leftAt, sc->partnerL) ||
        !pairShore(nR, sc->rightOff, sc->rightAt, sc->partnerR)) {
      return false;
    }
    const int narcs = m * 2;
    EulerScratch::fit(sc->succ, narcs);
    for (int e = 0; e < m; e++) {
      const int goRight = sc->partnerR[static_cast<std::size_t>(e)];
      const int goLeft = sc->partnerL[static_cast<std::size_t>(e)];
      if (goRight < 0 || goLeft < 0) {
        return false;
      }
      sc->succ[static_cast<std::size_t>(2 * e)] = 2 * goRight + 1;
      sc->succ[static_cast<std::size_t>(2 * e + 1)] = 2 * goLeft;
    }
    EulerScratch::fit(sc->seen, narcs);
    std::fill(sc->seen.begin(), sc->seen.begin() + narcs, 0);
    EulerScratch::fit(sc->rep, narcs);
    EulerScratch::fit(sc->arcRank, narcs);
    sc->cycle.clear();
    if (sc->cycle.capacity() < static_cast<std::size_t>(narcs)) {
      sc->cycle.reserve(static_cast<std::size_t>(narcs));
    }
    for (int start = 0; start < narcs; start++) {
      if (sc->seen[static_cast<std::size_t>(start)]) {
        continue;
      }
      sc->cycle.clear();
      int a = start;
      int guard = 0;
      do {
        if (a < 0 || a >= narcs || sc->seen[static_cast<std::size_t>(a)]) {
          return false;
        }
        sc->seen[static_cast<std::size_t>(a)] = 1;
        sc->cycle.push_back(a);
        a = sc->succ[static_cast<std::size_t>(a)];
        gDdcMatchSteps++;
        if (++guard > narcs) {
          return false;
        }
      } while (a != start);
      int minArc = sc->cycle[0];
      int minAt = 0;
      for (int i = 1; i < static_cast<int>(sc->cycle.size()); i++) {
        if (sc->cycle[static_cast<std::size_t>(i)] < minArc) {
          minArc = sc->cycle[static_cast<std::size_t>(i)];
          minAt = i;
        }
      }
      const int len = static_cast<int>(sc->cycle.size());
      for (int i = 0; i < len; i++) {
        const int arc = sc->cycle[static_cast<std::size_t>((minAt + i) % len)];
        sc->rep[static_cast<std::size_t>(arc)] = minArc;
        sc->arcRank[static_cast<std::size_t>(arc)] = i;
      }
    }
    EulerScratch::fit(sc->keep, m);
    std::fill(sc->keep.begin(), sc->keep.begin() + m, 0);
    EulerScratch::fit(sc->keptAt, nL);
    EulerScratch::fit(sc->keptRight, nR);
    std::fill(sc->keptAt.begin(), sc->keptAt.begin() + nL, 0);
    std::fill(sc->keptRight.begin(), sc->keptRight.begin() + nR, 0);
    for (int e = 0; e < m; e++) {
      const int a0 = 2 * e;
      const int a1 = 2 * e + 1;
      const int chosen = sc->rep[static_cast<std::size_t>(a0)] <= sc->rep[static_cast<std::size_t>(a1)] ? a0 : a1;
      if ((sc->arcRank[static_cast<std::size_t>(chosen)] % 2) == 0) {
        sc->keep[static_cast<std::size_t>(e)] = 1;
        sc->keptAt[static_cast<std::size_t>(sc->edges[static_cast<std::size_t>(e)].left)]++;
        sc->keptRight[static_cast<std::size_t>(sc->edges[static_cast<std::size_t>(e)].right)]++;
      }
    }
    const int half = degree / 2;
    for (int u = 0; u < nL; u++) {
      if (sc->keptAt[static_cast<std::size_t>(u)] != half) {
        return false;
      }
    }
    for (int r = 0; r < nR; r++) {
      if (sc->keptRight[static_cast<std::size_t>(r)] != half) {
        return false;
      }
    }
    sc->next.clear();
    if (sc->next.capacity() < static_cast<std::size_t>(m / 2)) {
      sc->next.reserve(static_cast<std::size_t>(m / 2));
    }
    for (int e = 0; e < m; e++) {
      if (sc->keep[static_cast<std::size_t>(e)]) {
        sc->next.push_back(sc->edges[static_cast<std::size_t>(e)]);
      }
    }
    sc->edges.swap(sc->next);
    m = static_cast<int>(sc->edges.size());
    degree = half;
  }
  if (m != nL) {
    return false;
  }
  EulerScratch::fit(sc->gotL, nL);
  EulerScratch::fit(sc->gotR, nR);
  EulerScratch::fit(sc->gotRing, nL);
  std::fill(sc->gotL.begin(), sc->gotL.begin() + nL, -1);
  std::fill(sc->gotR.begin(), sc->gotR.begin() + nR, -1);
  for (int e = 0; e < m; e++) {
    const HalfEdge &edge = sc->edges[static_cast<std::size_t>(e)];
    if (sc->gotL[static_cast<std::size_t>(edge.left)] >= 0 ||
        sc->gotR[static_cast<std::size_t>(edge.right)] >= 0) {
      return false;
    }
    sc->gotL[static_cast<std::size_t>(edge.left)] = edge.right;
    sc->gotR[static_cast<std::size_t>(edge.right)] = edge.left;
    sc->gotRing[static_cast<std::size_t>(edge.left)] = edge.ring;
  }
  matchL.assign(sc->gotL.begin(), sc->gotL.begin() + nL);
  matchR.assign(sc->gotR.begin(), sc->gotR.begin() + nR);
  matchRing.assign(sc->gotRing.begin(), sc->gotRing.begin() + nL);
  return true;
}

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
  if (parallelEulerMatch(nL, nR, off, adj, matchL, matchR, matchRing)) {
    return true;
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

ring::PackParts gPackParts;

void resetPackParts() { gPackParts = {}; }


// Same order as sorting the rank keys. Digits are the key ints, least
// significant first, so equal keys stay in input order and the run ids match
// a comparison sort. Ranks above 16 bits fall back to that sort.
// Ids are the sorted order of the packed rank keys. Equal keys share one
// id. A hash collects the distinct words, then those few words are numbered
// by sorting them, which is the same numbering as a full radix of every copy.
bool assignPackedIds(const std::vector<PKey> &keys, std::vector<std::array<int, 2>> &ends,
                     std::vector<int> &cliqueCount) {
  const int nKeys = static_cast<int>(keys.size());
  if (nKeys == 0) {
    cliqueCount.clear();
    return true;
  }
  thread_local std::vector<std::uint64_t> wordStore;
  std::vector<std::uint64_t> *words = &wordStore;
  asm volatile("" : "+r"(words));
  words->resize(static_cast<std::size_t>(nKeys));
  std::uint64_t *word = words->data();
  for (int i = 0; i < nKeys; i++) {
    const RankKey &key = keys[static_cast<std::size_t>(i)].k;
    if (!key.ok || key.n != 4) {
      return false;
    }
    std::uint64_t packed = 0;
    for (int d = 0; d < 4; d++) {
      if (key.id[d] < 0 || key.id[d] > 65535) {
        return false;
      }
      packed |= static_cast<std::uint64_t>(key.id[3 - d]) << (16 * d);
    }
    word[i] = packed;
  }
  int nTab = 1;
  while (nTab < nKeys) {
    nTab <<= 1;
  }
  // A generation stamp marks live slots, so a repeat search does not clear
  // the whole table. The ids are still the sorted order of the distinct words.
  struct Slot {
    std::uint64_t k = 0;
    int id = -1;
    std::uint32_t gen = 0;
  };
  thread_local std::vector<Slot> tabStore;
  thread_local std::uint32_t tabGen = 0;
  thread_local std::vector<std::uint64_t> uniqStore;
  thread_local std::vector<int> tmpStore;
  thread_local std::vector<int> orderStore;
  thread_local std::vector<int> rankStore;
  std::vector<Slot> *tabP = &tabStore;
  asm volatile("" : "+r"(tabP));
  if (static_cast<int>(tabP->size()) < nTab) {
    tabP->assign(static_cast<std::size_t>(nTab), Slot{});
  }
  if (++tabGen == 0) {
    for (Slot &slot : *tabP) {
      slot.gen = 0;
    }
    tabGen = 1;
  }
  const std::uint32_t gen = tabGen;
  Slot *tab = tabP->data();
  const std::size_t mask = static_cast<std::size_t>(nTab - 1);
  uniqStore.clear();
  if (uniqStore.capacity() < static_cast<std::size_t>(nKeys / 4 + 1)) {
    uniqStore.reserve(static_cast<std::size_t>(nKeys / 4 + 1));
  }
  tmpStore.resize(static_cast<std::size_t>(nKeys));
  for (int i = 0; i < nKeys; i++) {
    const std::uint64_t w = word[i];
    std::size_t h = static_cast<std::size_t>(w * 0x9E3779B185EBCA87ULL) & mask;
    int probes = 0;
    while (tab[h].gen == gen && tab[h].k != w) {
      h = (h + 1) & mask;
      if (++probes > nTab) {
        return false;
      }
    }
    if (tab[h].gen != gen) {
      tab[h].k = w;
      tab[h].id = static_cast<int>(uniqStore.size());
      tab[h].gen = gen;
      uniqStore.push_back(w);
    }
    tmpStore[static_cast<std::size_t>(i)] = tab[h].id;
  }
  const int nU = static_cast<int>(uniqStore.size());
  orderStore.resize(static_cast<std::size_t>(nU));
  std::iota(orderStore.begin(), orderStore.end(), 0);
  const std::uint64_t *uniq = uniqStore.data();
  std::sort(orderStore.begin(), orderStore.end(),
            [&](int a, int b) { return uniq[a] < uniq[b]; });
  rankStore.resize(static_cast<std::size_t>(nU));
  for (int i = 0; i < nU; i++) {
    rankStore[static_cast<std::size_t>(orderStore[static_cast<std::size_t>(i)])] = i;
  }
  cliqueCount.assign(static_cast<std::size_t>(std::max(nU, 1)), 0);
  for (int i = 0; i < nKeys; i++) {
    const int id = rankStore[static_cast<std::size_t>(tmpStore[static_cast<std::size_t>(i)])];
    const PKey &key = keys[static_cast<std::size_t>(i)];
    cliqueCount[static_cast<std::size_t>(id)]++;
    ends[static_cast<std::size_t>(key.ring)][static_cast<std::size_t>(key.slot)] = id;
  }
  gDdcKeyDigits += static_cast<std::uint64_t>(nKeys) * 4;
  return true;
}

void orderCliqueKeys(std::vector<int> &ord, const std::vector<PKey> &keys) {
  const int n = static_cast<int>(ord.size());
  auto fallback = [&]() {
    seams::parallelStableSort(ord, [&](int a, int b) {
      const PKey &ka = keys[static_cast<std::size_t>(a)];
      const PKey &kb = keys[static_cast<std::size_t>(b)];
      if (ka.k < kb.k) {
        return true;
      }
      if (kb.k < ka.k) {
        return false;
      }
      if (ka.ring != kb.ring) {
        return ka.ring < kb.ring;
      }
      return ka.slot < kb.slot;
    });
  };
  if (n < 2) {
    return;
  }
  int maxId = 0;
  for (const PKey &key : keys) {
    if (!key.k.ok || key.k.n < 0 || key.k.n > 4) {
      fallback();
      return;
    }
    for (int d = 0; d < key.k.n; d++) {
      if (key.k.id[d] < 0) {
        fallback();
        return;
      }
      maxId = std::max(maxId, key.k.id[d]);
    }
  }
  if (maxId > 65535) {
    fallback();
    return;
  }
  // Four 16-bit digits, id[3] in the low bits so a radix matches RankKey
  // order when every key has length 4. One scratch, reused.
  thread_local std::vector<std::uint64_t> packedStore;
  thread_local std::vector<int> cntStore;
  thread_local std::vector<int> bufStore;
  std::vector<std::uint64_t> *packedP = &packedStore;
  std::vector<int> *cntP = &cntStore;
  std::vector<int> *bufP = &bufStore;
  asm volatile("" : "+r"(packedP), "+r"(cntP), "+r"(bufP));
  packedP->resize(static_cast<std::size_t>(n));
  bool packedOk = true;
  for (int i = 0; i < n; i++) {
    const RankKey &key = keys[static_cast<std::size_t>(ord[static_cast<std::size_t>(i)])].k;
    if (key.n != 4) {
      packedOk = false;
      break;
    }
    std::uint64_t word = 0;
    for (int d = 0; d < 4; d++) {
      word |= static_cast<std::uint64_t>(key.id[3 - d]) << (16 * d);
    }
    (*packedP)[static_cast<std::size_t>(i)] = word;
  }
  if (!packedOk) {
    const int buckets = std::max(maxId, 4) + 1;
    cntP->assign(static_cast<std::size_t>(buckets), 0);
    bufP->resize(static_cast<std::size_t>(n));
    auto pass = [&](auto digit, int nBuckets) {
      std::fill(cntP->begin(), cntP->begin() + nBuckets, 0);
      for (int i = 0; i < n; i++) {
        (*cntP)[static_cast<std::size_t>(digit(ord[static_cast<std::size_t>(i)]))]++;
      }
      int sum = 0;
      for (int b = 0; b < nBuckets; b++) {
        const int take = (*cntP)[static_cast<std::size_t>(b)];
        (*cntP)[static_cast<std::size_t>(b)] = sum;
        sum += take;
      }
      for (int i = 0; i < n; i++) {
        const int id = ord[static_cast<std::size_t>(i)];
        (*bufP)[static_cast<std::size_t>((*cntP)[static_cast<std::size_t>(digit(id))]++)] = id;
      }
      ord.swap(*bufP);
    };
    pass([&](int i) { return keys[static_cast<std::size_t>(i)].k.n; }, 5);
    for (int d = 3; d >= 0; --d) {
      pass(
          [&](int i) {
            const RankKey &k = keys[static_cast<std::size_t>(i)].k;
            return d < k.n ? k.id[d] : 0;
          },
          maxId + 1);
    }
    gDdcKeyDigits += static_cast<std::uint64_t>(n) * 5;
    return;
  }
  const int buckets = maxId + 1;
  cntP->assign(static_cast<std::size_t>(buckets), 0);
  bufP->resize(static_cast<std::size_t>(n));
  std::uint64_t *word = packedP->data();
  int *cnt = cntP->data();
  int *buf = bufP->data();
  int *src = ord.data();
  int *dst = buf;
  for (int shift = 0; shift < 64; shift += 16) {
    std::fill(cnt, cnt + buckets, 0);
    for (int i = 0; i < n; i++) {
      cnt[(word[src[i]] >> shift) & 65535U]++;
    }
    int sum = 0;
    for (int b = 0; b < buckets; b++) {
      const int take = cnt[b];
      cnt[b] = sum;
      sum += take;
    }
    for (int i = 0; i < n; i++) {
      const int id = src[i];
      dst[cnt[(word[id] >> shift) & 65535U]++] = id;
    }
    std::swap(src, dst);
  }
  if (src != ord.data()) {
    std::copy(src, src + n, ord.begin());
  }
  gDdcKeyDigits += static_cast<std::uint64_t>(n) * 4;
  gDdcBucketSteps += static_cast<std::uint64_t>(buckets) * 4U;
  gPackedLive = packedP;
}

// Line graph of a bipartite graph on this component. Neighbourhoods split
// into the two endpoint cliques; a maximum set of rings is a maximum matching.
bool linePack(const std::vector<Row> &conf, const std::vector<int> &mem,
              const std::vector<int> &rankOf, std::vector<int> &out, int &matched,
              int &cover) {
  const auto splitStart = std::chrono::steady_clock::now();
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
#pragma omp parallel for schedule(static) if (!omp_in_parallel())
#endif
    for (int i = 0; i < m; i++) {
      const int v = mem[static_cast<size_t>(i)];
      const Row &nb = conf[static_cast<size_t>(v)];
      // Degree is six here, so the neighbour test reads the inline block.
      auto near = [&](int x, int y) {
        const int *head = conf[static_cast<size_t>(x)].head;
        return head[0] == y || head[1] == y || head[2] == y || head[3] == y || head[4] == y ||
               head[5] == y;
      };
      int A[3];
      int B[3];
      int nA = 1;
      int nB = 0;
      A[0] = nb.head[0];
      bool ok = true;
      for (int t = 1; t < 6; t++) {
        if (near(nb.head[0], nb.head[t])) {
          if (nA >= 3) {
            ok = false;
            break;
          }
          A[nA++] = nb.head[t];
        } else if (nB >= 3) {
          ok = false;
          break;
        } else {
          B[nB++] = nb.head[t];
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
      RankKey k0 = packRanks(left, 4, rankOf);
      RankKey k1 = packRanks(right, 4, rankOf);
      if (!k0.ok || !k1.ok || k0 == k1) {
        bad[static_cast<size_t>(i)] = 1;
        continue;
      }
      // Smaller key first, so clique ids follow ranks rather than storage order.
      if (k1 < k0) {
        const RankKey tmp = k0;
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
      A[0] = nb[0];
      bool ok = true;
      for (int i = 1; i < 6; i++) {
        if (linked(nb[0], nb[i])) {
          if (nA >= 3) {
            ok = false;
            break;
          }
          A[nA++] = nb[i];
        } else if (nB >= 3) {
          ok = false;
          break;
        } else {
          B[nB++] = nb[i];
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
          if (linked(nb[a], nb[b])) {
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
          A[nA++] = nb[a];
        } else if ((sideB >> a) & 1u) {
          B[nB++] = nb[a];
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
    RankKey k0 = packRanks(left, nL, rankOf);
    RankKey k1 = packRanks(right, nR, rankOf);
    if (!k0.ok || !k1.ok || k0 == k1) {
      packed = false;
      break;
    }
    if (k1 < k0) {
      const RankKey tmp = k0;
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
  gPackParts.split_ms += std::chrono::duration<double, std::milli>(
                             std::chrono::steady_clock::now() - splitStart)
                             .count();
  const auto cliqueStart = std::chrono::steady_clock::now();
  // Clique ids are the run index of the rank keys in sorted order. A prefix
  // sum of the run starts assigns those ids, and equal keys keep one id.
  std::vector<std::array<int, 2>> ends(static_cast<size_t>(nAll), std::array<int, 2>{-1, -1});
  const int nKeys = static_cast<int>(keys.size());
  std::vector<int> ord;
  std::vector<int> cliqueCount;
  if (!assignPackedIds(keys, ends, cliqueCount)) {
  ord.resize(static_cast<std::size_t>(nKeys));
  std::iota(ord.begin(), ord.end(), 0);
  gPackedLive = nullptr;
  orderCliqueKeys(ord, keys);
  if (gPackedLive != nullptr && static_cast<int>(gPackedLive->size()) == nKeys) {
    const std::uint64_t *word = gPackedLive->data();
    int runs = 0;
    std::uint64_t prev = 0;
    for (int i = 0; i < nKeys; i++) {
      const std::uint64_t w = word[ord[static_cast<std::size_t>(i)]];
      if (i == 0 || w != prev) {
        ++runs;
        prev = w;
      }
    }
    cliqueCount.assign(static_cast<std::size_t>(std::max(runs, 1)), 0);
    int id = -1;
    prev = 0;
    for (int i = 0; i < nKeys; i++) {
      const std::uint64_t w = word[ord[static_cast<std::size_t>(i)]];
      if (i == 0 || w != prev) {
        ++id;
        prev = w;
      }
      const PKey &key = keys[static_cast<std::size_t>(ord[static_cast<std::size_t>(i)])];
      cliqueCount[static_cast<std::size_t>(id)]++;
      ends[static_cast<std::size_t>(key.ring)][static_cast<std::size_t>(key.slot)] = id;
    }
  } else {
  std::vector<int> fresh(static_cast<std::size_t>(std::max(nKeys, 1)), 0);
  for (int i = 1; i < nKeys; i++) {
    const RankKey &prev = keys[static_cast<std::size_t>(ord[static_cast<std::size_t>(i - 1)])].k;
    const RankKey &cur = keys[static_cast<std::size_t>(ord[static_cast<std::size_t>(i)])].k;
    fresh[static_cast<std::size_t>(i)] = cur == prev ? 0 : 1;
  }
  std::vector<int> runId = fresh;
  seams::exclusiveScan(runId);
  for (int i = 0; i < nKeys; i++) {
    runId[static_cast<std::size_t>(i)] += fresh[static_cast<std::size_t>(i)];
  }
  const int nC = nKeys == 0 ? 0 : runId[static_cast<std::size_t>(nKeys - 1)] + 1;
  cliqueCount.assign(static_cast<std::size_t>(std::max(nC, 1)), 0);
  for (int i = 0; i < nKeys; i++) {
    const PKey &key = keys[static_cast<std::size_t>(ord[static_cast<std::size_t>(i)])];
    const int id = runId[static_cast<std::size_t>(i)];
    cliqueCount[static_cast<std::size_t>(id)]++;
    ends[static_cast<std::size_t>(key.ring)][static_cast<std::size_t>(key.slot)] = id;
  }
  }
  }
  gPackParts.clique_ms += std::chrono::duration<double, std::milli>(
                              std::chrono::steady_clock::now() - cliqueStart)
                              .count();
  struct MatchClock {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    ~MatchClock() {
      gPackParts.match_ms += std::chrono::duration<double, std::milli>(
                                 std::chrono::steady_clock::now() - t0)
                                 .count();
    }
  } matchClock;
  struct FreeKeys {
    std::vector<PKey> &keys;
    std::vector<std::array<int, 2>> &ends;
    std::vector<int> &ord;
    std::vector<int> &cliqueCount;
    ~FreeKeys() {
      std::vector<PKey>().swap(keys);
      std::vector<std::array<int, 2>>().swap(ends);
      std::vector<int>().swap(ord);
      std::vector<int>().swap(cliqueCount);
    }
  } freeKeys{keys, ends, ord, cliqueCount};
  const int nC = static_cast<int>(cliqueCount.size());
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
      if (stampV[static_cast<size_t>(row[i])] != era) {
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
  // Each free left vertex proposes its lowest-rank free neighbour. Each right
  // vertex accepts the lowest-rank proposal. The rounds are a maximal
  // matching, and the augmenting search below grows it to a maximum one.
  bool grew = true;
  int guard = 0;
  while (grew && guard++ <= nL) {
    grew = false;
    std::vector<int> offerR(static_cast<std::size_t>(nL), -1);
    std::vector<int> offerRing(static_cast<std::size_t>(nL), -1);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (!omp_in_parallel())
#endif
    for (int u = 0; u < nL; u++) {
      if (matchL[static_cast<std::size_t>(u)] >= 0) {
        continue;
      }
      int bestR = -1;
      int bestRing = -1;
      int bestRank = 0;
      for (int e = off[static_cast<std::size_t>(u)]; e < off[static_cast<std::size_t>(u + 1)]; e++) {
        const int r = adj[static_cast<std::size_t>(e)].right;
        const int ring = adj[static_cast<std::size_t>(e)].ring;
        if (matchR[static_cast<std::size_t>(r)] >= 0) {
          continue;
        }
        const int rank = rankOf[static_cast<std::size_t>(ring)];
        if (bestR < 0 || rank < bestRank || (rank == bestRank && ring < bestRing)) {
          bestR = r;
          bestRing = ring;
          bestRank = rank;
        }
      }
      offerR[static_cast<std::size_t>(u)] = bestR;
      offerRing[static_cast<std::size_t>(u)] = bestRing;
    }
    std::vector<int> bestL(static_cast<std::size_t>(nR), -1);
    std::vector<int> bestRank(static_cast<std::size_t>(nR), 0);
    std::vector<int> bestRing(static_cast<std::size_t>(nR), -1);
    for (int u = 0; u < nL; u++) {
      const int r = offerR[static_cast<std::size_t>(u)];
      if (r < 0 || matchR[static_cast<std::size_t>(r)] >= 0 ||
          matchL[static_cast<std::size_t>(u)] >= 0) {
        continue;
      }
      const int ring = offerRing[static_cast<std::size_t>(u)];
      const int rank = rankOf[static_cast<std::size_t>(ring)];
      const int held = bestL[static_cast<std::size_t>(r)];
      if (held < 0 || rank < bestRank[static_cast<std::size_t>(r)] ||
          (rank == bestRank[static_cast<std::size_t>(r)] && u < held)) {
        bestL[static_cast<std::size_t>(r)] = u;
        bestRank[static_cast<std::size_t>(r)] = rank;
        bestRing[static_cast<std::size_t>(r)] = ring;
      }
    }
    for (int r = 0; r < nR; r++) {
      const int u = bestL[static_cast<std::size_t>(r)];
      if (u < 0) {
        continue;
      }
      matchL[static_cast<std::size_t>(u)] = r;
      matchR[static_cast<std::size_t>(r)] = u;
      matchRing[static_cast<std::size_t>(u)] = bestRing[static_cast<std::size_t>(r)];
      grew = true;
    }
  }
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
          if (mem[static_cast<size_t>(j)] == row[k] && ((mask >> j) & 1)) {
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
  RowSpillGuard spillGuard;
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
  // A component past the old fixed word count must still be solved.
  {
    constexpr int wide = 300;
    std::vector<Row> conf(static_cast<size_t>(wide));
    std::vector<int> mem(static_cast<size_t>(wide));
    for (int i = 0; i < wide; i++) {
      mem[static_cast<size_t>(i)] = i;
      if (i + 1 < wide) {
        conf[static_cast<size_t>(i)].push(i + 1);
        conf[static_cast<size_t>(i + 1)].push(i);
      }
    }
    const MisOut path = exactMis(conf, mem);
    if (!path.solved || static_cast<int>(path.ids.size()) != wide / 2) {
      std::fprintf(stderr, "mis-selftest FAIL path %d solved %d got %zu\n", wide,
                   path.solved ? 1 : 0, path.ids.size());
      std::exit(1);
    }
  }
  std::fprintf(stderr, "mis-selftest ok\n");
}

std::vector<int> ring::maximumEquators(const std::vector<ring::EquatorChoice> &periph,
                                       const std::vector<int> &order,
                                       ring::PackingCertificate *certificate) {
  RowSpillGuard spillGuard;
  static int selfTested = 0;
  if (selfTested == 0 && std::getenv("DSEAMS_MIS_SELFTEST") != nullptr) {
    selfTested = 1;
    misSelfTest();
  }
  const auto t0 = std::chrono::steady_clock::now();
  resetPackParts();
  gDdcMatchSteps = 0;
  gDdcRowCopies = 0;
  gDdcKeyDigits = 0;
  gDdcBucketSteps = 0;
  const int n = static_cast<int>(periph.size());
  std::vector<Row> conf(static_cast<size_t>(n));
  int nEdges = 0;
  int nCand = 0;
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel for schedule(static) reduction(+ : nEdges, nCand) if (!omp_in_parallel())
#endif
  for (int i = 0; i < n; i++) {
    if (periph[static_cast<size_t>(i)].n == 0) {
      continue;
    }
    nCand++;
    const ring::EquatorChoice &row = periph[static_cast<size_t>(i)];
    for (int k = 0; k < row.n; k++) {
      const int f = row.idAt(k);
      if (f < 0 || f >= n || f == i || periph[static_cast<size_t>(f)].n == 0) {
        continue;
      }
      bool back = false;
      const ring::EquatorChoice &other = periph[static_cast<size_t>(f)];
      for (int t = 0; t < other.n; t++) {
        if (other.idAt(t) == i) {
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
  int nDropped = 0;
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
        const int w = row[k];
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
  gPackParts.conf_ms =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  for (Comp &comp : comps) {
    const auto orderStart = std::chrono::steady_clock::now();
    std::sort(comp.mem.begin(), comp.mem.end(), [&](int a, int b) {
      return rankOf[static_cast<size_t>(a)] < rankOf[static_cast<size_t>(b)];
    });
    gPackParts.conf_ms += std::chrono::duration<double, std::milli>(
                              std::chrono::steady_clock::now() - orderStart)
                              .count();
    std::vector<int> picked;
    int partMatch = 0;
    int partCover = 0;
    if (linePack(conf, comp.mem, rankOf, picked, partMatch, partCover)) {
      const auto keepStart = std::chrono::steady_clock::now();
      take.insert(take.end(), picked.begin(), picked.end());
      gPackParts.conf_ms += std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - keepStart)
                                .count();
      matchSum += partMatch;
      coverSum += partCover;
      continue;
    }
    anySearch = true;
    const auto searchStart = std::chrono::steady_clock::now();
    const MisOut solved = exactMis(conf, comp.mem);
    gPackParts.search_ms += std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - searchStart)
                                .count();
    if (!solved.solved) {
      allSolved = false;
      nDropped++;
      continue;
    }
    take.insert(take.end(), solved.ids.begin(), solved.ids.end());
    matchSum += static_cast<int>(solved.ids.size());
    coverSum += solved.cover;
  }
  const auto tailStart = std::chrono::steady_clock::now();
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
    certificate->dropped = nDropped;
    certificate->method = anySearch ? 2 : 1;
    certificate->tight = allSolved && coverSum == cages;
    certificate->milliseconds =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  }
  gPackParts.conf_ms +=
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tailStart)
          .count();
  // The named blocks are the clique split, the id assignment, the matching
  // and the search. The conflict block is the rest of this call: the row
  // build, the component walk and the bookkeeping between those blocks.
  std::vector<Row>().swap(conf);
  std::vector<Comp>().swap(comps);
  std::vector<int>().swap(seen);
  std::vector<int>().swap(rankOf);
  const double packWall =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  const double named = gPackParts.split_ms + gPackParts.clique_ms + gPackParts.match_ms +
                       gPackParts.search_ms;
  if (packWall > named) {
    gPackParts.conf_ms = packWall - named;
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

ring::PackParts ring::lastPackParts() { return gPackParts; }

extern "C" __attribute__((visibility("default"), used)) std::uint64_t seams_ddc_match_steps() {
  return gDdcMatchSteps;
}

extern "C" __attribute__((visibility("default"), used)) std::uint64_t seams_ddc_row_copies() {
  return gDdcRowCopies;
}

extern "C" __attribute__((visibility("default"), used)) std::uint64_t seams_ddc_key_digits() {
  return gDdcKeyDigits;
}

extern "C" __attribute__((visibility("default"), used)) std::uint64_t seams_ddc_bucket_steps() {
  return gDdcBucketSteps;
}

