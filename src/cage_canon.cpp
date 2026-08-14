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
#include <nauty.h>
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
std::string densenautyCertificate(int n, const std::vector<std::pair<int, int>> &edges) {
  if (n <= 0 || n > MAXN) {
    return {};
  }
  DYNALLSTAT(graph, g, g_sz);
  DYNALLSTAT(graph, cg, cg_sz);
  DYNALLSTAT(int, lab, lab_sz);
  DYNALLSTAT(int, ptn, ptn_sz);
  DYNALLSTAT(int, orbits, orbits_sz);
  DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;
  const int m = SETWORDSNEEDED(n);
  options.getcanon = TRUE;
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
  DYNALLOC2(graph, g, g_sz, n, m, "densenauty");
  DYNALLOC2(graph, cg, cg_sz, n, m, "densenauty");
  DYNALLOC1(int, lab, lab_sz, n, "densenauty");
  DYNALLOC1(int, ptn, ptn_sz, n, "densenauty");
  DYNALLOC1(int, orbits, orbits_sz, n, "densenauty");
  EMPTYGRAPH(g, m, n);
  for (const auto &e : edges) {
    if (e.first == e.second) {
      continue;
    }
    ADDONEEDGE(g, e.first, e.second, m);
  }
  densenauty(g, lab, ptn, orbits, &options, &stats, m, n, cg);
  std::ostringstream out;
  out << std::hex;
  const int words = m * n;
  for (int i = 0; i < words; i++) {
    out << cg[i];
    if (i + 1 < words) {
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

} // namespace cage
