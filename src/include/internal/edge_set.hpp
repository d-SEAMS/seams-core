#ifndef SEAMS_EDGE_SET_H_
#define SEAMS_EDGE_SET_H_

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

namespace seams {

// Undirected edges, self entries dropped, each pair stored once as (lo, hi)
// and sorted. The hash is FNV-1a over that list. Reuse is allowed only when
// the hashes match and the lists then compare equal. A hash match alone is
// not enough.
struct EdgeSet {
  std::uint64_t hash = 1469598103934665603ull;
  std::vector<std::pair<int, int>> edges;
};

inline EdgeSet edgeSetOf(const std::vector<std::vector<int>> &rows) {
  EdgeSet out;
  for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
    const auto &row = rows[static_cast<std::size_t>(i)];
    for (std::size_t k = 1; k < row.size(); ++k) {
      const int j = row[k];
      if (j < 0 || j == i) {
        continue;
      }
      out.edges.emplace_back(std::min(i, j), std::max(i, j));
    }
  }
  std::sort(out.edges.begin(), out.edges.end());
  out.edges.erase(std::unique(out.edges.begin(), out.edges.end()),
                  out.edges.end());
  std::uint64_t h = 1469598103934665603ull;
  for (const auto &e : out.edges) {
    h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(e.first));
    h *= 1099511628211ull;
    h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(e.second));
    h *= 1099511628211ull;
  }
  out.hash = h;
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
