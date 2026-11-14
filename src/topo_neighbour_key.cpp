//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------
// Keep this translation unit out of the cold partition. A cold clone
// lands before the ring search and moves it.
#pragma GCC optimize("no-reorder-blocks-and-partition")

#include <topo_fingerprint.hpp>

#include <algorithm>
#include <cstring>
#include <limits>

namespace topo {

struct LocalGraph {
  std::vector<int> atoms;
  std::vector<std::vector<int>> adjacency;
  std::vector<int> localColors;
  std::vector<std::pair<int, int>> edges;
};

LocalKey localKeyFromGraph(LocalGraph &g, int hops, std::uint64_t *refinement);
void fillLocalGraph(LocalGraph &g, const Rows &rows, int atom, int hops,
                    const std::vector<int> &colors);

namespace {

// One flat image of a small rooted ball. Ice is seventeen vertices; a
// vector per vertex costs more than the walk, and the same image repeats.
constexpr int kFlatVerts = 64;
constexpr int kFlatNbrs = kFlatVerts * 16;

struct NbrCsr {
  std::vector<int> off;
  std::vector<int> col;
};

struct Image {
  int n = 0;
  int nNbr = 0;
  int nCol = 0;
  int off[kFlatVerts + 1]{};
  int nbr[kFlatNbrs]{};
  int colors[kFlatVerts]{};

  bool same(const Image &o) const {
    return n == o.n && nNbr == o.nNbr && nCol == o.nCol &&
           std::memcmp(off, o.off, static_cast<std::size_t>(n + 1) * sizeof(int)) == 0 &&
           std::memcmp(nbr, o.nbr, static_cast<std::size_t>(nNbr) * sizeof(int)) == 0 &&
           (nCol == 0 ||
            std::memcmp(colors, o.colors, static_cast<std::size_t>(nCol) * sizeof(int)) == 0);
  }
};

void buildNbrCsr(const Rows &rows, NbrCsr &csr) {
  const int n = static_cast<int>(rows.size());
  csr.off.resize(static_cast<std::size_t>(n) + 1);
  csr.off[0] = 0;
  int count = 0;
  for (int i = 0; i < n; i++) {
    const auto &row = rows[static_cast<std::size_t>(i)];
    count += row.empty() ? 0 : static_cast<int>(row.size()) - 1;
    csr.off[static_cast<std::size_t>(i) + 1] = count;
  }
  csr.col.resize(static_cast<std::size_t>(count));
  for (int i = 0; i < n; i++) {
    const auto &row = rows[static_cast<std::size_t>(i)];
    int at = csr.off[static_cast<std::size_t>(i)];
    for (std::size_t m = 1; m < row.size(); m++) {
      csr.col[static_cast<std::size_t>(at++)] = row[m];
    }
  }
}

struct Seen {
  std::vector<int> seen;
  std::vector<int> local;
  std::vector<int> depth;
  int epoch = 0;

  void prepare(int nRows) {
    if (static_cast<int>(seen.size()) < nRows) {
      seen.assign(static_cast<std::size_t>(nRows), 0);
      local.assign(static_cast<std::size_t>(nRows), 0);
      depth.assign(static_cast<std::size_t>(nRows), 0);
      epoch = 0;
    }
    if (epoch == std::numeric_limits<int>::max()) {
      std::fill(seen.begin(), seen.end(), 0);
      epoch = 0;
    }
    ++epoch;
  }
};

bool hopImage(const NbrCsr &csr, int atom, int hops, int n, Seen &vis, int *atoms, int &count) {
  count = 0;
  if (atom < 0 || atom >= n) {
    return true;
  }
  vis.prepare(n);
  vis.depth[static_cast<std::size_t>(atom)] = 0;
  vis.seen[static_cast<std::size_t>(atom)] = vis.epoch;
  atoms[count++] = atom;
  int head = 0;
  while (head < count) {
    const int v = atoms[head++];
    const int d = vis.depth[static_cast<std::size_t>(v)];
    if (d == hops) {
      continue;
    }
    const int begin = csr.off[static_cast<std::size_t>(v)];
    const int end = csr.off[static_cast<std::size_t>(v) + 1];
    for (int m = begin; m < end; m++) {
      const int u = csr.col[static_cast<std::size_t>(m)];
      if (u < 0 || u >= n || vis.seen[static_cast<std::size_t>(u)] == vis.epoch) {
        continue;
      }
      if (count == kFlatVerts) {
        return false;
      }
      vis.depth[static_cast<std::size_t>(u)] = d + 1;
      vis.seen[static_cast<std::size_t>(u)] = vis.epoch;
      atoms[count++] = u;
    }
  }
  return true;
}

bool fillImage(Image &img, const int *atoms, int count, const NbrCsr &csr, int nRows,
               const std::vector<int> &colors, Seen &vis) {
  img.n = count;
  img.nNbr = 0;
  img.nCol = 0;
  img.off[0] = 0;
  vis.prepare(nRows);
  for (int i = 0; i < count; i++) {
    const int a = atoms[i];
    vis.seen[static_cast<std::size_t>(a)] = vis.epoch;
    vis.local[static_cast<std::size_t>(a)] = i;
  }
  for (int i = 0; i < count; i++) {
    const int v = atoms[i];
    const int begin = csr.off[static_cast<std::size_t>(v)];
    const int end = csr.off[static_cast<std::size_t>(v) + 1];
    int at = img.nNbr;
    for (int m = begin; m < end; m++) {
      const int u = csr.col[static_cast<std::size_t>(m)];
      if (u < 0 || u >= nRows || vis.seen[static_cast<std::size_t>(u)] != vis.epoch) {
        continue;
      }
      const int j = vis.local[static_cast<std::size_t>(u)];
      if (j == i) {
        continue;
      }
      if (at == kFlatNbrs) {
        return false;
      }
      img.nbr[at++] = j;
    }
    for (int a = img.nNbr + 1; a < at; a++) {
      const int val = img.nbr[a];
      int b = a;
      while (b > img.nNbr && img.nbr[b - 1] > val) {
        img.nbr[b] = img.nbr[b - 1];
        --b;
      }
      img.nbr[b] = val;
    }
    img.nNbr = at;
    img.off[i + 1] = at;
  }
  if (colors.size() == static_cast<std::size_t>(nRows) && nRows > 0) {
    img.nCol = count;
    for (int i = 0; i < count; i++) {
      img.colors[i] = colors[static_cast<std::size_t>(atoms[i])];
    }
  }
  return true;
}

void materialize(LocalGraph &g, const int *atoms, const Image &img) {
  g.atoms.assign(atoms, atoms + img.n);
  g.edges.clear();
  g.adjacency.resize(static_cast<std::size_t>(img.n));
  for (int i = 0; i < img.n; i++) {
    const int begin = img.off[i];
    const int end = img.off[i + 1];
    g.adjacency[static_cast<std::size_t>(i)].assign(img.nbr + begin, img.nbr + end);
  }
  if (img.nCol == 0) {
    g.localColors.clear();
  } else {
    g.localColors.assign(img.colors, img.colors + img.nCol);
  }
}

} // namespace

LocalKey flatNeighbourKey(const Rows &rows, int atom, int hops, const std::vector<int> &colors,
                         std::uint64_t *refined) {
  // One thread. The packed rows are rebuilt when the frame changes, which
  // for a full fingerprint is the first atom. A team never calls this.
  const int n = static_cast<int>(rows.size());
  thread_local const void *cached = nullptr;
  thread_local int cachedN = -1;
  thread_local NbrCsr csr;
  if (atom == 0 || cached != rows.data() || cachedN != n) {
    buildNbrCsr(rows, csr);
    cached = rows.data();
    cachedN = n;
  }
  thread_local Seen vis;
  thread_local Image scratch;
  thread_local Image memo;
  thread_local LocalGraph g;
  thread_local LocalKey memoKey;
  thread_local std::uint64_t memoRefined = 0;
  thread_local int memoHops = -1;
  thread_local bool memoReady = false;
  int atomsBuf[kFlatVerts];
  int count = 0;
  std::uint64_t got = 0;
  const bool fitted = hopImage(csr, atom, hops, n, vis, atomsBuf, count);
  if (fitted && count > 0 && fillImage(scratch, atomsBuf, count, csr, n, colors, vis)) {
    if (memoReady && memoHops == hops && memo.same(scratch)) {
      if (refined != nullptr) {
        *refined = memoRefined;
      }
      return memoKey;
    }
    materialize(g, atomsBuf, scratch);
    LocalKey lk = localKeyFromGraph(g, hops, &got);
    memo = scratch;
    memoKey = lk;
    memoRefined = got;
    memoHops = hops;
    memoReady = true;
    if (refined != nullptr) {
      *refined = got;
    }
    return lk;
  }
  fillLocalGraph(g, rows, atom, hops, colors);
  return localKeyFromGraph(g, hops, refined);
}

} // namespace topo
