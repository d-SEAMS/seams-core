//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <cage_canon.hpp>

#include <algorithm>
#include <cstdint>
#include <map>
#include <sstream>
#include <vector>

#ifdef SEAMS_HAS_NAUTY
#define MAXN 64
#if __has_include(<nauty/nauty.h>)
#include <nauty/nauty.h>
#else
#include <nauty.h>
#endif
#endif

namespace {

void collectAtoms(const std::vector<std::vector<int>> &rings,
                  std::vector<int> &atoms, std::map<int, int> &toLocal) {
  atoms.clear();
  toLocal.clear();
  for (const auto &ring : rings) {
    for (int a : ring) {
      if (toLocal.find(a) == toLocal.end()) {
        toLocal[a] = static_cast<int>(atoms.size());
        atoms.push_back(a);
      }
    }
  }
}

#ifdef SEAMS_HAS_NAUTY
std::string densenautyCertificate(int n, const std::vector<std::pair<int, int>> &edges,
                                  int root = -1) {
  if (n <= 0 || n > MAXN) {
    return {};
  }
  DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;
  const int m = SETWORDSNEEDED(n);
  options.getcanon = TRUE;
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
  std::vector<graph> g(static_cast<size_t>(n * m), 0);
  std::vector<graph> cg(static_cast<size_t>(n * m), 0);
  std::vector<int> lab(static_cast<size_t>(n), 0);
  std::vector<int> ptn(static_cast<size_t>(n), 0);
  std::vector<int> orbits(static_cast<size_t>(n), 0);
  if (root >= 0 && root < n) {
    // two colour cells: {root} and the rest; ptn marks the end of a cell
    // with 0 and the inside of a cell with 1
    options.defaultptn = FALSE;
    int pos = 0;
    lab[static_cast<size_t>(pos)] = root;
    ptn[static_cast<size_t>(pos)] = 0;
    ++pos;
    for (int v = 0; v < n; v++) {
      if (v == root) {
        continue;
      }
      lab[static_cast<size_t>(pos)] = v;
      ptn[static_cast<size_t>(pos)] = 1;
      ++pos;
    }
    ptn[static_cast<size_t>(n - 1)] = 0;
  }
  EMPTYGRAPH(g.data(), m, n);
  for (const auto &e : edges) {
    if (e.first == e.second) {
      continue;
    }
    ADDONEEDGE(g.data(), e.first, e.second, m);
  }
  densenauty(g.data(), lab.data(), ptn.data(), orbits.data(), &options, &stats,
             m, n, cg.data());
  std::ostringstream out;
  out << std::hex;
  for (int i = 0; i < n * m; i++) {
    out << cg[static_cast<size_t>(i)];
    if (i + 1 < n * m) {
      out << '-';
    }
  }
  return out.str();
}
#endif

std::vector<std::pair<int, int>>
ringEdges(const std::vector<std::vector<int>> &rings,
          const std::map<int, int> &toLocal) {
  std::vector<std::pair<int, int>> edges;
  for (const auto &ring : rings) {
    if (ring.size() < 2) {
      continue;
    }
    for (size_t i = 0; i < ring.size(); i++) {
      const int a = toLocal.at(ring[i]);
      const int b = toLocal.at(ring[(i + 1) % ring.size()]);
      if (a == b) {
        continue;
      }
      edges.emplace_back(std::min(a, b), std::max(a, b));
    }
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  return edges;
}

#ifdef SEAMS_HAS_NAUTY
std::string hexPrismCertificate() {
  std::vector<std::pair<int, int>> edges;
  for (int i = 0; i < 6; i++) {
    edges.emplace_back(i, (i + 1) % 6);
    edges.emplace_back(6 + i, 6 + (i + 1) % 6);
    edges.emplace_back(i, 6 + i);
  }
  return densenautyCertificate(12, edges);
}
#endif

} // namespace

namespace cage {

bool nautyAvailable() {
#ifdef SEAMS_HAS_NAUTY
  return true;
#else
  return false;
#endif
}

std::string canonicalCertificate(const std::vector<std::vector<int>> &rings) {
#ifndef SEAMS_HAS_NAUTY
  (void)rings;
  return {};
#else
  std::vector<int> atoms;
  std::map<int, int> toLocal;
  collectAtoms(rings, atoms, toLocal);
  if (atoms.empty()) {
    return {};
  }
  return densenautyCertificate(static_cast<int>(atoms.size()),
                               ringEdges(rings, toLocal));
#endif
}

bool isHexagonalPrism(const std::vector<std::vector<int>> &rings) {
#ifndef SEAMS_HAS_NAUTY
  (void)rings;
  return false;
#else
  const std::string got = canonicalCertificate(rings);
  if (got.empty()) {
    return false;
  }
  return got == hexPrismCertificate();
#endif
}

bool sameCertificate(const std::vector<std::vector<int>> &a,
                     const std::vector<std::vector<int>> &b) {
  const std::string ca = canonicalCertificate(a);
  const std::string cb = canonicalCertificate(b);
  return !ca.empty() && ca == cb;
}

std::string canonicalCertificateRooted(int n,
                                       const std::vector<std::pair<int, int>> &edges,
                                       int root) {
#ifndef SEAMS_HAS_NAUTY
  (void)n;
  (void)edges;
  (void)root;
  return {};
#else
  std::vector<std::pair<int, int>> local;
  local.reserve(edges.size());
  for (const auto &e : edges) {
    if (e.first < 0 || e.second < 0 || e.first >= n || e.second >= n || e.first == e.second) {
      continue;
    }
    local.emplace_back(std::min(e.first, e.second), std::max(e.first, e.second));
  }
  std::sort(local.begin(), local.end());
  local.erase(std::unique(local.begin(), local.end()), local.end());
  return densenautyCertificate(n, local, root);
#endif
}

} // namespace cage
