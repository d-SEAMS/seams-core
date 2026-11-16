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

#include <chordless.hpp>

#include <algorithm>
#include <cstdint>
#include <stdexcept>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

// A cycle of length L through its lowest vertex s splits at the vertex c
// opposite s into two induced paths from s: a short half of p = L / 2
// steps and a long half of q = L - p steps. The cycle is chordless exactly
// when both halves are induced, their interiors are disjoint, and no edge
// joins the two interiors. Walking every path of L - 1 steps from s costs
// about d^(L-1) visits at degree d; the halves cost about d^q, and the
// pairing touches only halves that end on the same vertex.

namespace {

struct Scratch {
  std::vector<unsigned char> near; // see walk()
  std::vector<int> head;           // endpoint -> first long half, or -1
  std::vector<std::uint32_t> stamp;
  std::uint32_t token = 0;
  std::vector<int> touched;
  std::vector<int> shortHalves; // p vertices each, s implied; unused when p == q
  std::vector<int> longHalves;  // q vertices each, s implied
  std::vector<int> link;        // next long half with the same endpoint
  std::vector<int> local;       // this source's cycles
  int path[chordless::kMaxSize];
};

struct Search {
  const int *off;
  const int *cols;
  int size;
  int p;
  int q;
};

// near[w] counts the path vertices before the last one that w touches.
// A candidate with a nonzero count would close a chord, so each step reads
// one byte per candidate instead of scanning the candidate's row. A vertex
// already on the path touches its predecessor and so has a nonzero count,
// except s, which v > s excludes.
void walk(const Search &g, Scratch &scr, int s, int k) {
  if (k >= 2) {
    const int prev = scr.path[k - 2];
    for (int b = g.off[prev]; b < g.off[prev + 1]; ++b) {
      ++scr.near[static_cast<std::size_t>(g.cols[b])];
    }
  }
  const int last = scr.path[k - 1];
  for (int a = g.off[last]; a < g.off[last + 1]; ++a) {
    const int v = g.cols[a];
    if (v <= s || scr.near[static_cast<std::size_t>(v)] != 0) {
      continue;
    }
    scr.path[k] = v;
    if (k == g.p && g.p != g.q) {
      scr.shortHalves.insert(scr.shortHalves.end(), scr.path + 1, scr.path + 1 + g.p);
    }
    if (k == g.q) {
      scr.longHalves.insert(scr.longHalves.end(), scr.path + 1, scr.path + 1 + g.q);
      continue;
    }
    walk(g, scr, s, k + 1);
  }
  if (k >= 2) {
    const int prev = scr.path[k - 2];
    for (int b = g.off[prev]; b < g.off[prev + 1]; ++b) {
      --scr.near[static_cast<std::size_t>(g.cols[b])];
    }
  }
}

// Appends the cycles of source s to scr.local, size vertices each.
void fromSource(const Search &g, Scratch &scr, int s) {
  scr.local.clear();
  if (g.size == 3) {
    for (int a = g.off[s]; a < g.off[s + 1]; ++a) {
      const int u = g.cols[a];
      if (u <= s) {
        continue;
      }
      for (int b = g.off[u]; b < g.off[u + 1]; ++b) {
        const int w = g.cols[b];
        if (w <= u) {
          continue;
        }
        for (int c = g.off[s]; c < g.off[s + 1]; ++c) {
          if (g.cols[c] == w) {
            scr.local.insert(scr.local.end(), {s, u, w});
            break;
          }
        }
      }
    }
    return;
  }
  scr.shortHalves.clear();
  scr.longHalves.clear();
  scr.path[0] = s;
  walk(g, scr, s, 1);

  const int p = g.p;
  const int q = g.q;
  const int nLong = static_cast<int>(scr.longHalves.size()) / q;
  scr.link.resize(static_cast<std::size_t>(nLong));
  for (int j = 0; j < nLong; ++j) {
    const int c = scr.longHalves[static_cast<std::size_t>(j * q + q - 1)];
    int &first = scr.head[static_cast<std::size_t>(c)];
    if (first < 0) {
      scr.touched.push_back(c);
    }
    scr.link[static_cast<std::size_t>(j)] = first;
    first = j;
  }
  const std::vector<int> &shortHalves = p == q ? scr.longHalves : scr.shortHalves;
  const int nShort = static_cast<int>(shortHalves.size()) / p;
  for (int i = 0; i < nShort; ++i) {
    const int *hp = shortHalves.data() + static_cast<std::size_t>(i) * p;
    const int c = hp[p - 1];
    for (int j = scr.head[static_cast<std::size_t>(c)]; j >= 0;
         j = scr.link[static_cast<std::size_t>(j)]) {
      const int *hq = scr.longHalves.data() + static_cast<std::size_t>(j) * q;
      if (hp[0] >= hq[0]) {
        continue;
      }
      if (++scr.token == 0) {
        std::fill(scr.stamp.begin(), scr.stamp.end(), 0u);
        scr.token = 1;
      }
      const std::uint32_t tok = scr.token;
      for (int k = 0; k + 1 < q; ++k) {
        scr.stamp[static_cast<std::size_t>(hq[k])] = tok;
      }
      bool ok = true;
      for (int k = 0; k + 1 < p && ok; ++k) {
        const int u = hp[k];
        if (scr.stamp[static_cast<std::size_t>(u)] == tok) {
          ok = false;
          break;
        }
        for (int b = g.off[u]; b < g.off[u + 1]; ++b) {
          if (scr.stamp[static_cast<std::size_t>(g.cols[b])] == tok) {
            ok = false;
            break;
          }
        }
      }
      if (!ok) {
        continue;
      }
      scr.local.push_back(s);
      scr.local.insert(scr.local.end(), hp, hp + p);
      for (int k = q - 2; k >= 0; --k) {
        scr.local.push_back(hq[k]);
      }
    }
  }
  for (int c : scr.touched) {
    scr.head[static_cast<std::size_t>(c)] = -1;
  }
  scr.touched.clear();
}

void sortLocal(std::vector<int> &local, int size, std::vector<int> &tmp) {
  const std::size_t count = local.size() / static_cast<std::size_t>(size);
  if (count < 2) {
    return;
  }
  std::vector<int> order(count);
  for (std::size_t i = 0; i < count; ++i) {
    order[i] = static_cast<int>(i);
  }
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    return std::lexicographical_compare(
        local.begin() + static_cast<std::ptrdiff_t>(a) * size,
        local.begin() + static_cast<std::ptrdiff_t>(a + 1) * size,
        local.begin() + static_cast<std::ptrdiff_t>(b) * size,
        local.begin() + static_cast<std::ptrdiff_t>(b + 1) * size);
  });
  tmp.resize(local.size());
  for (std::size_t i = 0; i < count; ++i) {
    std::copy_n(local.begin() + static_cast<std::ptrdiff_t>(order[i]) * size, size,
                tmp.begin() + static_cast<std::ptrdiff_t>(i) * size);
  }
  local.swap(tmp);
}

// Sources per thread below which a team costs more than it saves.
constexpr int kSourcesPerThread = 512;

} // namespace

std::size_t chordless::cycles(const std::vector<int> &offsets, const std::vector<int> &cols,
                              int size, std::vector<int> *rings, int *per_atom) {
  if (size < 3 || size > kMaxSize) {
    throw std::invalid_argument("chordless cycle size must be from 3 to 12");
  }
  const int n = offsets.empty() ? 0 : static_cast<int>(offsets.size()) - 1;
  if (per_atom != nullptr) {
    std::fill(per_atom, per_atom + n, 0);
  }
  if (rings != nullptr) {
    rings->clear();
  }
  if (n <= 0) {
    return 0;
  }
  const Search g{offsets.data(), cols.data(), size, size / 2, size - size / 2};
  const bool keep = rings != nullptr;
  std::vector<int> perSource(keep ? static_cast<std::size_t>(n) : 0, 0);
  std::size_t total = 0;
  int team = 1;
#ifdef SEAMS_HAS_OPENMP
  team = std::max(1, std::min(omp_get_max_threads(), n / kSourcesPerThread));
#endif
  std::vector<std::vector<int>> buffers(static_cast<std::size_t>(team));
  std::vector<std::vector<int>> blockSource(static_cast<std::size_t>(team));

#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel num_threads(team) reduction(+ : total)
#endif
  {
    int tid = 0;
#ifdef SEAMS_HAS_OPENMP
    tid = omp_get_thread_num();
#endif
    Scratch scr;
    scr.near.assign(static_cast<std::size_t>(n), 0);
    scr.head.assign(static_cast<std::size_t>(n), -1);
    scr.stamp.assign(static_cast<std::size_t>(n), 0u);
    std::vector<int> tmp;
    auto &buf = buffers[static_cast<std::size_t>(tid)];
    auto &mine = blockSource[static_cast<std::size_t>(tid)];
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
    for (int s = 0; s < n; ++s) {
      fromSource(g, scr, s);
      const std::size_t found = scr.local.size() / static_cast<std::size_t>(size);
      if (found == 0) {
        continue;
      }
      total += found;
      if (per_atom != nullptr) {
        for (int v : scr.local) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp atomic
#endif
          per_atom[v] += 1;
        }
      }
      if (keep) {
        sortLocal(scr.local, size, tmp);
        perSource[static_cast<std::size_t>(s)] = static_cast<int>(found);
        mine.push_back(s);
        buf.insert(buf.end(), scr.local.begin(), scr.local.end());
      }
    }
  }

  if (keep) {
    std::vector<std::size_t> start(static_cast<std::size_t>(n) + 1, 0);
    for (int s = 0; s < n; ++s) {
      start[static_cast<std::size_t>(s) + 1] =
          start[static_cast<std::size_t>(s)] +
          static_cast<std::size_t>(perSource[static_cast<std::size_t>(s)]);
    }
    rings->resize(total * static_cast<std::size_t>(size));
    for (int t = 0; t < team; ++t) {
      const auto &buf = buffers[static_cast<std::size_t>(t)];
      std::size_t at = 0;
      for (int s : blockSource[static_cast<std::size_t>(t)]) {
        const std::size_t len = static_cast<std::size_t>(perSource[static_cast<std::size_t>(s)]) *
                                static_cast<std::size_t>(size);
        std::copy_n(buf.begin() + static_cast<std::ptrdiff_t>(at), len,
                    rings->begin() +
                        static_cast<std::ptrdiff_t>(start[static_cast<std::size_t>(s)] *
                                                    static_cast<std::size_t>(size)));
        at += len;
      }
    }
  }
  return total;
}
