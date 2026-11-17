#ifndef SEAMS_EDGE_SET_H_
#define SEAMS_EDGE_SET_H_

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace seams {

// Undirected edges, self entries dropped, each pair stored once as (lo, hi)
// and sorted. The hash is FNV-1a over that list. Reuse is allowed only when
// the hashes match and the lists then compare equal. A hash match alone is
// not enough.
struct EdgeSet {
  std::uint64_t hash = 1469598103934665603ull;
  std::vector<std::pair<int, int>> edges;
};

namespace {

void sortEdgePairs(std::vector<std::pair<int, int>> &edges) {
#if !defined(SEAMS_HAS_OPENMP)
  std::sort(edges.begin(), edges.end());
#else
  const std::size_t n = edges.size();
  const int nT = omp_get_max_threads();
  if (n < 8192 || nT < 2 || omp_in_parallel()) {
    std::sort(edges.begin(), edges.end());
    return;
  }
  const std::size_t chunk =
      (n + static_cast<std::size_t>(nT) - 1) / static_cast<std::size_t>(nT);
#pragma omp parallel for schedule(static)
  for (int t = 0; t < nT; ++t) {
    const std::size_t begin = static_cast<std::size_t>(t) * chunk;
    if (begin >= n) {
      continue;
    }
    const std::size_t end = std::min(n, begin + chunk);
    std::sort(edges.begin() + static_cast<std::ptrdiff_t>(begin),
              edges.begin() + static_cast<std::ptrdiff_t>(end));
  }
  for (std::size_t width = chunk; width < n; width *= 2) {
    const std::size_t span = width * 2;
    const int blocks =
        static_cast<int>((n + span - 1) / span);
#pragma omp parallel for schedule(static)
    for (int b = 0; b < blocks; ++b) {
      const std::size_t begin = static_cast<std::size_t>(b) * span;
      const std::size_t mid = std::min(n, begin + width);
      const std::size_t end = std::min(n, begin + span);
      if (begin < mid && mid < end) {
        std::inplace_merge(edges.begin() + static_cast<std::ptrdiff_t>(begin),
                           edges.begin() + static_cast<std::ptrdiff_t>(mid),
                           edges.begin() + static_cast<std::ptrdiff_t>(end));
      }
    }
  }
#endif
}

void hashEdges(EdgeSet &out) {
  std::uint64_t h = 1469598103934665603ull;
  for (const auto &e : out.edges) {
    h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(e.first));
    h *= 1099511628211ull;
    h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(e.second));
    h *= 1099511628211ull;
  }
  out.hash = h;
}

} // namespace

inline EdgeSet edgeSetOf(const std::vector<std::vector<int>> &rows) {
  EdgeSet out;
  const int nRows = static_cast<int>(rows.size());
  bool parallel = false;
#ifdef SEAMS_HAS_OPENMP
  parallel = nRows >= 4096 && omp_get_max_threads() > 1 && !omp_in_parallel();
#endif
  if (parallel) {
#ifdef SEAMS_HAS_OPENMP
    const int team = omp_get_max_threads();
    std::vector<std::vector<std::pair<int, int>>> parts(
        static_cast<std::size_t>(team));
#pragma omp parallel
    {
      std::vector<std::pair<int, int>> local;
#pragma omp for schedule(static)
      for (int i = 0; i < nRows; ++i) {
        const auto &row = rows[static_cast<std::size_t>(i)];
        for (std::size_t k = 1; k < row.size(); ++k) {
          const int j = row[k];
          if (j < 0 || j == i) {
            continue;
          }
          local.emplace_back(std::min(i, j), std::max(i, j));
        }
      }
      parts[static_cast<std::size_t>(omp_get_thread_num())] = std::move(local);
    }
    std::size_t total = 0;
    for (const auto &part : parts) {
      total += part.size();
    }
    out.edges.reserve(total);
    for (auto &part : parts) {
      out.edges.insert(out.edges.end(), part.begin(), part.end());
    }
#endif
  } else {
    for (int i = 0; i < nRows; ++i) {
      const auto &row = rows[static_cast<std::size_t>(i)];
      for (std::size_t k = 1; k < row.size(); ++k) {
        const int j = row[k];
        if (j < 0 || j == i) {
          continue;
        }
        out.edges.emplace_back(std::min(i, j), std::max(i, j));
      }
    }
  }
  sortEdgePairs(out.edges);
  out.edges.erase(std::unique(out.edges.begin(), out.edges.end()), out.edges.end());
  hashEdges(out);
  return out;
}

inline bool provenEqual(const EdgeSet &a, const EdgeSet &b) {
  if (a.hash != b.hash) {
    return false;
  }
  return a.edges == b.edges;
}

} // namespace seams

#endif
