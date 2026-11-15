//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------
#include <topo_fingerprint.hpp>

#include <cage_canon.hpp>
#include <franzblau.hpp>

namespace primitive {
std::vector<int> ringCensusCounts(const std::vector<std::vector<int>> &nList, int maxDepth);
}

#include <algorithm>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

constexpr std::uint64_t kOffset = 1469598103934665603ULL;
constexpr std::uint64_t kPrime = 1099511628211ULL;

std::uint64_t mix(std::uint64_t h, std::uint64_t v) {
  // FNV-1a over the eight bytes of v, then a final avalanche
  for (int b = 0; b < 8; b++) {
    h ^= (v >> (8 * b)) & 0xffULL;
    h *= kPrime;
  }
  h ^= h >> 29;
  h *= 0xbf58476d1ce4e5b9ULL;
  h ^= h >> 32;
  return h;
}

void fillRingCensus(std::vector<int> &census,
                    const std::vector<std::vector<int>> &rows, int maxRingSize) {
  census.assign(static_cast<std::size_t>(maxRingSize) + 1, 0);
  if (rows.empty() || maxRingSize < 3) {
    return;
  }
  // The frame key reads the census, not the ring list. Count while
  // searching so a sparse frame does not store every ring either.
  census = primitive::ringCensusCounts(rows, maxRingSize);
  if (static_cast<int>(census.size()) < maxRingSize + 1) {
    census.resize(static_cast<std::size_t>(maxRingSize) + 1, 0);
  }
}

std::uint64_t hashString(const std::string &s) {
  std::uint64_t h = kOffset;
  for (unsigned char c : s) {
    h ^= c;
    h *= kPrime;
  }
  return h;
}

} // namespace

namespace topo {

std::string hex(std::uint64_t value) {
  static constexpr char kDigits[] = "0123456789abcdef";
  std::string out(16, '0');
  for (int i = 15; i >= 0; --i) {
    out[static_cast<std::size_t>(i)] = kDigits[value & 0xfu];
    value >>= 4;
  }
  return out;
}

void hopNeighbourhoodInto(const Rows &rows, int atom, int hops, std::vector<int> &order) {
  order.clear();
  if (atom < 0 || atom >= static_cast<int>(rows.size())) {
    return;
  }
  // One stamp for the walk. The visit order is the neighbour-row order,
  // the same order a map of first sightings produces.
  thread_local std::vector<int> seen;
  thread_local std::vector<int> depth;
  thread_local int epoch = 0;
  const int n = static_cast<int>(rows.size());
  if (static_cast<int>(seen.size()) < n) {
    seen.assign(static_cast<std::size_t>(n), 0);
    depth.assign(static_cast<std::size_t>(n), 0);
    epoch = 0;
  }
  if (epoch == std::numeric_limits<int>::max()) {
    std::fill(seen.begin(), seen.end(), 0);
    epoch = 0;
  }
  ++epoch;
  depth[static_cast<std::size_t>(atom)] = 0;
  seen[static_cast<std::size_t>(atom)] = epoch;
  order.push_back(atom);
  std::size_t head = 0;
  while (head < order.size()) {
    const int v = order[head++];
    const int d = depth[static_cast<std::size_t>(v)];
    if (d == hops) {
      continue;
    }
    const auto &row = rows[static_cast<std::size_t>(v)];
    for (std::size_t m = 1; m < row.size(); m++) {
      const int u = row[m];
      if (u < 0 || u >= n || seen[static_cast<std::size_t>(u)] == epoch) {
        continue;
      }
      depth[static_cast<std::size_t>(u)] = d + 1;
      seen[static_cast<std::size_t>(u)] = epoch;
      order.push_back(u);
    }
  }
}

std::vector<int> hopNeighbourhood(const Rows &rows, int atom, int hops) {
  std::vector<int> order;
  hopNeighbourhoodInto(rows, atom, hops, order);
  return order;
}

std::uint64_t wlHash(const std::vector<std::vector<int>> &adjacency, int root,
                     int rounds, const std::vector<int> &colors) {
  const int n = static_cast<int>(adjacency.size());
  const bool colored = static_cast<int>(colors.size()) == n;
  // The same four buffers serve every neighbourhood. Each call overwrites
  // the n colors it reads, so a later graph does not see a previous one.
  thread_local std::vector<std::uint64_t> color;
  thread_local std::vector<std::uint64_t> next;
  thread_local std::vector<std::uint64_t> bag;
  thread_local std::vector<std::uint64_t> all;
  color.resize(static_cast<std::size_t>(n));
  next.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    const std::uint64_t deg = static_cast<std::uint64_t>(adjacency[static_cast<std::size_t>(i)].size());
    std::uint64_t h = mix(kOffset, deg + (i == root ? 0x9e3779b97f4a7c15ULL : 0ULL));
    if (colored) {
      h = mix(h, static_cast<std::uint64_t>(static_cast<std::int64_t>(colors[static_cast<std::size_t>(i)])) + 1ULL);
    }
    color[static_cast<std::size_t>(i)] = h;
  }
  for (int r = 0; r < rounds; r++) {
    for (int i = 0; i < n; i++) {
      bag.clear();
      for (int j : adjacency[static_cast<std::size_t>(i)]) {
        bag.push_back(color[static_cast<std::size_t>(j)]);
      }
      // Neighbourhoods on a sparse frame are a handful of colors. An
      // insertion sort of that bag is the refinement; a long bag still
      // uses the general sort.
      if (bag.size() <= 16) {
        for (std::size_t a = 1; a < bag.size(); a++) {
          const std::uint64_t v = bag[a];
          std::size_t b = a;
          while (b > 0 && bag[b - 1] > v) {
            bag[b] = bag[b - 1];
            --b;
          }
          bag[b] = v;
        }
      } else {
        std::sort(bag.begin(), bag.end());
      }
      std::uint64_t h = mix(kOffset, color[static_cast<std::size_t>(i)]);
      for (std::uint64_t c : bag) {
        h = mix(h, c);
      }
      next[static_cast<std::size_t>(i)] = h;
    }
    color.swap(next);
  }
  all.assign(color.begin(), color.begin() + n);
  std::sort(all.begin(), all.end());
  std::uint64_t h = mix(kOffset, static_cast<std::uint64_t>(n));
  if (root >= 0 && root < n) {
    h = mix(h, color[static_cast<std::size_t>(root)]);
  }
  for (std::uint64_t c : all) {
    h = mix(h, c);
  }
  return h;
}

// Same refinement as wlHash. A wide ball repeats neighbour colors, so one
// sort names them and each vertex counts copies instead of sorting its list.
std::uint64_t wlHashWide(const std::vector<std::vector<int>> &adjacency, int root,
                         int rounds, const std::vector<int> &colors) {
  const int n = static_cast<int>(adjacency.size());
  const bool colored = static_cast<int>(colors.size()) == n;
  thread_local std::vector<std::uint64_t> color;
  thread_local std::vector<std::uint64_t> next;
  thread_local std::vector<std::uint64_t> all;
  thread_local std::vector<std::uint64_t> palette;
  thread_local std::vector<int> colorId;
  thread_local std::vector<int> counts;
  thread_local std::vector<std::pair<std::uint64_t, int>> order;
  color.resize(static_cast<std::size_t>(n));
  next.resize(static_cast<std::size_t>(n));
  bool uniform = n >= 2 && root >= 0 && root < n;
  const std::size_t deg0 = adjacency.empty() ? 0 : adjacency[0].size();
  const int paint0 = colored ? colors[0] : 0;
  for (int i = 0; i < n; i++) {
    const std::size_t deg = adjacency[static_cast<std::size_t>(i)].size();
    if (deg != deg0) {
      uniform = false;
    }
    if (colored && colors[static_cast<std::size_t>(i)] != paint0) {
      uniform = false;
    }
    std::uint64_t h = mix(kOffset, static_cast<std::uint64_t>(deg) +
                                       (i == root ? 0x9e3779b97f4a7c15ULL : 0ULL));
    if (colored) {
      h = mix(h, static_cast<std::uint64_t>(
                     static_cast<std::int64_t>(colors[static_cast<std::size_t>(i)])) +
                     1ULL);
    }
    color[static_cast<std::size_t>(i)] = h;
  }
  for (int r = 0; r < rounds; r++) {
    // Every non-root starts on one color, so the first round is the root's
    // neighbour list: one hash for the root, one for its neighbours, one for
    // the rest. Later rounds still walk the edges.
    if (r == 0 && uniform) {
      const std::uint64_t common = color[static_cast<std::size_t>(root == 0 ? 1 : 0)];
      const std::uint64_t rooted = color[static_cast<std::size_t>(root)];
      thread_local std::vector<int> hits;
      hits.assign(static_cast<std::size_t>(n), 0);
      bool simple = true;
      for (const int j : adjacency[static_cast<std::size_t>(root)]) {
        if (j < 0 || j >= n || ++hits[static_cast<std::size_t>(j)] > 1) {
          simple = false;
          break;
        }
      }
      if (simple) {
        const int deg = static_cast<int>(deg0);
        auto sprinkle = [](std::uint64_t h, std::uint64_t col, int copies) {
          for (int t = 0; t < copies; t++) {
            h = mix(h, col);
          }
          return h;
        };
        const std::uint64_t hFar = sprinkle(mix(kOffset, common), common, deg);
        std::uint64_t hTouch = mix(kOffset, common);
        if (rooted < common) {
          hTouch = sprinkle(mix(hTouch, rooted), common, deg - 1);
        } else if (common < rooted) {
          hTouch = sprinkle(hTouch, common, deg - 1);
          hTouch = mix(hTouch, rooted);
        } else {
          hTouch = hFar;
        }
        const std::uint64_t hRoot = sprinkle(mix(kOffset, rooted), common, deg);
        for (int i = 0; i < n; i++) {
          if (i == root) {
            next[static_cast<std::size_t>(i)] = hRoot;
          } else if (hits[static_cast<std::size_t>(i)] != 0) {
            next[static_cast<std::size_t>(i)] = hTouch;
          } else {
            next[static_cast<std::size_t>(i)] = hFar;
          }
        }
        color.swap(next);
        continue;
      }
    }
    order.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
      order[static_cast<std::size_t>(i)] = {color[static_cast<std::size_t>(i)], i};
    }
    std::sort(order.begin(), order.end());
    palette.clear();
    colorId.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
      if (palette.empty() || palette.back() != order[static_cast<std::size_t>(i)].first) {
        palette.push_back(order[static_cast<std::size_t>(i)].first);
      }
      colorId[static_cast<std::size_t>(order[static_cast<std::size_t>(i)].second)] =
          static_cast<int>(palette.size()) - 1;
    }
    const int nColors = static_cast<int>(palette.size());
    if (nColors > 64) {
      thread_local std::vector<std::uint64_t> bag;
      for (int i = 0; i < n; i++) {
        bag.clear();
        for (const int j : adjacency[static_cast<std::size_t>(i)]) {
          bag.push_back(color[static_cast<std::size_t>(j)]);
        }
        std::sort(bag.begin(), bag.end());
        std::uint64_t h = mix(kOffset, color[static_cast<std::size_t>(i)]);
        for (const std::uint64_t c : bag) {
          h = mix(h, c);
        }
        next[static_cast<std::size_t>(i)] = h;
      }
    } else if (n >= 400) {
      // The mix of a vertex is fixed by its own color and the neighbour
      // counts. A dense crystal repeats that signature, so the long mix
      // runs once per signature. A smaller ball does not repay the grouping.
      thread_local std::vector<int> grid;
      thread_local std::vector<int> ord;
      const int cols = std::max(nColors, 1);
      grid.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(cols), 0);
      if (nColors > 0) {
        for (int i = 0; i < n; i++) {
          int *row = grid.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(cols);
          for (const int j : adjacency[static_cast<std::size_t>(i)]) {
            row[colorId[static_cast<std::size_t>(j)]]++;
          }
        }
      }
      ord.resize(static_cast<std::size_t>(n));
      std::iota(ord.begin(), ord.end(), 0);
      std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        const int ca = colorId[static_cast<std::size_t>(a)];
        const int cb = colorId[static_cast<std::size_t>(b)];
        if (ca != cb) {
          return ca < cb;
        }
        const int *ra = grid.data() + static_cast<std::size_t>(a) * static_cast<std::size_t>(cols);
        const int *rb = grid.data() + static_cast<std::size_t>(b) * static_cast<std::size_t>(cols);
        return std::lexicographical_compare(ra, ra + cols, rb, rb + cols);
      });
      auto same = [&](int a, int b) {
        if (colorId[static_cast<std::size_t>(a)] != colorId[static_cast<std::size_t>(b)]) {
          return false;
        }
        const int *ra = grid.data() + static_cast<std::size_t>(a) * static_cast<std::size_t>(cols);
        const int *rb = grid.data() + static_cast<std::size_t>(b) * static_cast<std::size_t>(cols);
        return std::equal(ra, ra + cols, rb);
      };
      for (int k = 0; k < n; k++) {
        const int i = ord[static_cast<std::size_t>(k)];
        if (k > 0 && same(i, ord[static_cast<std::size_t>(k - 1)])) {
          next[static_cast<std::size_t>(i)] =
              next[static_cast<std::size_t>(ord[static_cast<std::size_t>(k - 1)])];
          continue;
        }
        std::uint64_t h = mix(kOffset, color[static_cast<std::size_t>(i)]);
        if (nColors > 0) {
          const int *row = grid.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(cols);
          for (int c = 0; c < nColors; c++) {
            const std::uint64_t col = palette[static_cast<std::size_t>(c)];
            const int copies = row[c];
            for (int t = 0; t < copies; t++) {
              h = mix(h, col);
            }
          }
        }
        next[static_cast<std::size_t>(i)] = h;
      }
    } else {
      counts.assign(static_cast<std::size_t>(std::max(nColors, 0)), 0);
      for (int i = 0; i < n; i++) {
        std::fill(counts.begin(), counts.end(), 0);
        for (const int j : adjacency[static_cast<std::size_t>(i)]) {
          counts[static_cast<std::size_t>(colorId[static_cast<std::size_t>(j)])]++;
        }
        std::uint64_t h = mix(kOffset, color[static_cast<std::size_t>(i)]);
        for (int c = 0; c < nColors; c++) {
          const std::uint64_t col = palette[static_cast<std::size_t>(c)];
          const int copies = counts[static_cast<std::size_t>(c)];
          for (int t = 0; t < copies; t++) {
            h = mix(h, col);
          }
        }
        next[static_cast<std::size_t>(i)] = h;
      }
    }
    color.swap(next);
  }
  all.assign(color.begin(), color.begin() + n);
  std::sort(all.begin(), all.end());
  std::uint64_t h = mix(kOffset, static_cast<std::uint64_t>(n));
  if (root >= 0 && root < n) {
    h = mix(h, color[static_cast<std::size_t>(root)]);
  }
  for (const std::uint64_t c : all) {
    h = mix(h, c);
  }
  return h;
}

// The ball under 400 vertices. A wide palette is counted, not sorted:
// sorting each neighbour list costs more than the count. The hash matches
// the bag sort.
std::uint64_t wlHashCount(const std::vector<std::vector<int>> &adjacency, int root,
                         int rounds, const std::vector<int> &colors) {
  const int n = static_cast<int>(adjacency.size());
  const bool colored = static_cast<int>(colors.size()) == n;
  thread_local std::vector<std::uint64_t> color;
  thread_local std::vector<std::uint64_t> next;
  thread_local std::vector<std::uint64_t> all;
  thread_local std::vector<std::uint64_t> palette;
  thread_local std::vector<int> colorId;
  thread_local std::vector<int> counts;
  thread_local std::vector<std::pair<std::uint64_t, int>> order;
  color.resize(static_cast<std::size_t>(n));
  next.resize(static_cast<std::size_t>(n));
  bool uniform = n >= 2 && root >= 0 && root < n;
  const std::size_t deg0 = adjacency.empty() ? 0 : adjacency[0].size();
  const int paint0 = colored ? colors[0] : 0;
  for (int i = 0; i < n; i++) {
    const std::size_t deg = adjacency[static_cast<std::size_t>(i)].size();
    if (deg != deg0) {
      uniform = false;
    }
    if (colored && colors[static_cast<std::size_t>(i)] != paint0) {
      uniform = false;
    }
    std::uint64_t h = mix(kOffset, static_cast<std::uint64_t>(deg) +
                                       (i == root ? 0x9e3779b97f4a7c15ULL : 0ULL));
    if (colored) {
      h = mix(h, static_cast<std::uint64_t>(
                     static_cast<std::int64_t>(colors[static_cast<std::size_t>(i)])) +
                     1ULL);
    }
    color[static_cast<std::size_t>(i)] = h;
  }
  for (int r = 0; r < rounds; r++) {
    // Every non-root starts on one color, so the first round is the root's
    // neighbour list: one hash for the root, one for its neighbours, one for
    // the rest. Later rounds still walk the edges.
    if (r == 0 && uniform) {
      const std::uint64_t common = color[static_cast<std::size_t>(root == 0 ? 1 : 0)];
      const std::uint64_t rooted = color[static_cast<std::size_t>(root)];
      thread_local std::vector<int> hits;
      hits.assign(static_cast<std::size_t>(n), 0);
      bool simple = true;
      for (const int j : adjacency[static_cast<std::size_t>(root)]) {
        if (j < 0 || j >= n || ++hits[static_cast<std::size_t>(j)] > 1) {
          simple = false;
          break;
        }
      }
      if (simple) {
        const int deg = static_cast<int>(deg0);
        auto sprinkle = [](std::uint64_t h, std::uint64_t col, int copies) {
          for (int t = 0; t < copies; t++) {
            h = mix(h, col);
          }
          return h;
        };
        const std::uint64_t hFar = sprinkle(mix(kOffset, common), common, deg);
        std::uint64_t hTouch = mix(kOffset, common);
        if (rooted < common) {
          hTouch = sprinkle(mix(hTouch, rooted), common, deg - 1);
        } else if (common < rooted) {
          hTouch = sprinkle(hTouch, common, deg - 1);
          hTouch = mix(hTouch, rooted);
        } else {
          hTouch = hFar;
        }
        const std::uint64_t hRoot = sprinkle(mix(kOffset, rooted), common, deg);
        for (int i = 0; i < n; i++) {
          if (i == root) {
            next[static_cast<std::size_t>(i)] = hRoot;
          } else if (hits[static_cast<std::size_t>(i)] != 0) {
            next[static_cast<std::size_t>(i)] = hTouch;
          } else {
            next[static_cast<std::size_t>(i)] = hFar;
          }
        }
        color.swap(next);
        continue;
      }
    }
    order.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
      order[static_cast<std::size_t>(i)] = {color[static_cast<std::size_t>(i)], i};
    }
    std::sort(order.begin(), order.end());
    palette.clear();
    colorId.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
      if (palette.empty() || palette.back() != order[static_cast<std::size_t>(i)].first) {
        palette.push_back(order[static_cast<std::size_t>(i)].first);
      }
      colorId[static_cast<std::size_t>(order[static_cast<std::size_t>(i)].second)] =
          static_cast<int>(palette.size()) - 1;
    }
    const int nColors = static_cast<int>(palette.size());
    if (nColors > 64 && n >= 400) {
      thread_local std::vector<std::uint64_t> bag;
      for (int i = 0; i < n; i++) {
        bag.clear();
        for (const int j : adjacency[static_cast<std::size_t>(i)]) {
          bag.push_back(color[static_cast<std::size_t>(j)]);
        }
        std::sort(bag.begin(), bag.end());
        std::uint64_t h = mix(kOffset, color[static_cast<std::size_t>(i)]);
        for (const std::uint64_t c : bag) {
          h = mix(h, c);
        }
        next[static_cast<std::size_t>(i)] = h;
      }
    } else if (n >= 400) {
      // The mix of a vertex is fixed by its own color and the neighbour
      // counts. A dense crystal repeats that signature, so the long mix
      // runs once per signature. A smaller ball does not repay the grouping.
      thread_local std::vector<int> grid;
      thread_local std::vector<int> ord;
      const int cols = std::max(nColors, 1);
      grid.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(cols), 0);
      if (nColors > 0) {
        for (int i = 0; i < n; i++) {
          int *row = grid.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(cols);
          for (const int j : adjacency[static_cast<std::size_t>(i)]) {
            row[colorId[static_cast<std::size_t>(j)]]++;
          }
        }
      }
      ord.resize(static_cast<std::size_t>(n));
      std::iota(ord.begin(), ord.end(), 0);
      std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        const int ca = colorId[static_cast<std::size_t>(a)];
        const int cb = colorId[static_cast<std::size_t>(b)];
        if (ca != cb) {
          return ca < cb;
        }
        const int *ra = grid.data() + static_cast<std::size_t>(a) * static_cast<std::size_t>(cols);
        const int *rb = grid.data() + static_cast<std::size_t>(b) * static_cast<std::size_t>(cols);
        return std::lexicographical_compare(ra, ra + cols, rb, rb + cols);
      });
      auto same = [&](int a, int b) {
        if (colorId[static_cast<std::size_t>(a)] != colorId[static_cast<std::size_t>(b)]) {
          return false;
        }
        const int *ra = grid.data() + static_cast<std::size_t>(a) * static_cast<std::size_t>(cols);
        const int *rb = grid.data() + static_cast<std::size_t>(b) * static_cast<std::size_t>(cols);
        return std::equal(ra, ra + cols, rb);
      };
      for (int k = 0; k < n; k++) {
        const int i = ord[static_cast<std::size_t>(k)];
        if (k > 0 && same(i, ord[static_cast<std::size_t>(k - 1)])) {
          next[static_cast<std::size_t>(i)] =
              next[static_cast<std::size_t>(ord[static_cast<std::size_t>(k - 1)])];
          continue;
        }
        std::uint64_t h = mix(kOffset, color[static_cast<std::size_t>(i)]);
        if (nColors > 0) {
          const int *row = grid.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(cols);
          for (int c = 0; c < nColors; c++) {
            const std::uint64_t col = palette[static_cast<std::size_t>(c)];
            const int copies = row[c];
            for (int t = 0; t < copies; t++) {
              h = mix(h, col);
            }
          }
        }
        next[static_cast<std::size_t>(i)] = h;
      }
    } else {
      counts.assign(static_cast<std::size_t>(std::max(nColors, 0)), 0);
      for (int i = 0; i < n; i++) {
        std::fill(counts.begin(), counts.end(), 0);
        for (const int j : adjacency[static_cast<std::size_t>(i)]) {
          counts[static_cast<std::size_t>(colorId[static_cast<std::size_t>(j)])]++;
        }
        std::uint64_t h = mix(kOffset, color[static_cast<std::size_t>(i)]);
        for (int c = 0; c < nColors; c++) {
          const std::uint64_t col = palette[static_cast<std::size_t>(c)];
          const int copies = counts[static_cast<std::size_t>(c)];
          for (int t = 0; t < copies; t++) {
            h = mix(h, col);
          }
        }
        next[static_cast<std::size_t>(i)] = h;
      }
    }
    color.swap(next);
  }
  all.assign(color.begin(), color.begin() + n);
  std::sort(all.begin(), all.end());
  std::uint64_t h = mix(kOffset, static_cast<std::uint64_t>(n));
  if (root >= 0 && root < n) {
    h = mix(h, color[static_cast<std::size_t>(root)]);
  }
  for (const std::uint64_t c : all) {
    h = mix(h, c);
  }
  return h;
}

struct LocalGraph {
  std::vector<int> atoms;
  std::vector<std::vector<int>> adjacency;
  std::vector<int> localColors;
  std::vector<std::pair<int, int>> edges;
};

void fillLocalGraph(LocalGraph &g, const Rows &rows, int atom, int hops,
                    const std::vector<int> &colors) {
  g.edges.clear();
  g.localColors.clear();
  hopNeighbourhoodInto(rows, atom, hops, g.atoms);
  if (g.atoms.empty()) {
    g.adjacency.clear();
    return;
  }
  // Local index of each atom in this ball. The stamp is the call, so a
  // neighbour outside the ball is simply unmarked.
  thread_local std::vector<int> seen;
  thread_local std::vector<int> local;
  thread_local int epoch = 0;
  const int nRows = static_cast<int>(rows.size());
  if (static_cast<int>(seen.size()) < nRows) {
    seen.assign(static_cast<std::size_t>(nRows), 0);
    local.assign(static_cast<std::size_t>(nRows), 0);
    epoch = 0;
  }
  if (epoch == std::numeric_limits<int>::max()) {
    std::fill(seen.begin(), seen.end(), 0);
    epoch = 0;
  }
  ++epoch;
  for (std::size_t i = 0; i < g.atoms.size(); i++) {
    const int a = g.atoms[i];
    seen[static_cast<std::size_t>(a)] = epoch;
    local[static_cast<std::size_t>(a)] = static_cast<int>(i);
  }
  const int n = static_cast<int>(g.atoms.size());
  g.adjacency.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    const auto &row = rows[static_cast<std::size_t>(g.atoms[static_cast<std::size_t>(i)])];
    auto &adj = g.adjacency[static_cast<std::size_t>(i)];
    adj.clear();
    for (std::size_t m = 1; m < row.size(); m++) {
      const int u = row[m];
      if (u < 0 || u >= nRows || seen[static_cast<std::size_t>(u)] != epoch) {
        continue;
      }
      const int j = local[static_cast<std::size_t>(u)];
      if (j == i) {
        continue;
      }
      adj.push_back(j);
    }
    // A repeated crystal ball matches only when the rows are ordered. A long
    // row in a disordered ball does not repeat, and sorting it costs more
    // than the match it would find.
    if (adj.size() <= 16) {
      std::sort(adj.begin(), adj.end());
    }
  }
  const bool colored = colors.size() == rows.size() && !rows.empty();
  if (colored) {
    g.localColors.reserve(g.atoms.size());
    for (int a : g.atoms) {
      g.localColors.push_back(colors[static_cast<std::size_t>(a)]);
    }
  }
}

LocalKey localKeyFromGraph(LocalGraph &g, int hops, std::uint64_t *refinement) {
  LocalKey out;
  if (g.atoms.empty()) {
    out.method = "wl";
    out.key = hex(0);
    if (refinement != nullptr) {
      *refinement = 0;
    }
    return out;
  }
  // A crystal repeats one rooted ball. The refinement depends only on that
  // adjacency and the vertex colors, so a repeat is the key already hashed.
  struct Memo {
    int hops = -1;
    std::vector<std::vector<int>> adjacency;
    std::vector<int> colors;
    LocalKey key;
    std::uint64_t refined = 0;
    bool ready = false;
  };
  thread_local Memo memo;
  if (memo.ready && memo.hops == hops && memo.colors == g.localColors &&
      memo.adjacency == g.adjacency) {
    if (refinement != nullptr) {
      *refinement = memo.refined;
    }
    return memo.key;
  }
  if (g.edges.empty()) {
    for (int i = 0; i < static_cast<int>(g.adjacency.size()); i++) {
      for (const int j : g.adjacency[static_cast<std::size_t>(i)]) {
        if (i < j) {
          g.edges.emplace_back(i, j);
        }
      }
    }
  }
  std::sort(g.edges.begin(), g.edges.end());
  g.edges.erase(std::unique(g.edges.begin(), g.edges.end()), g.edges.end());
  out.vertices = static_cast<int>(g.atoms.size());
  out.edges = static_cast<int>(g.edges.size());
  // The frame key is this refinement even when the local key is a
  // certificate, so it is hashed once and handed back to the caller.
  // The wide refinement counts repeated colors instead of sorting
  // every neighbour list. It is the same hash; a dense ball has a
  // long list at every vertex.
  // Under 400 vertices a wide palette is counted. A larger ball keeps
  // the bag sort used for a dense shared frame.
  const int ballN = static_cast<int>(g.adjacency.size());
  const std::uint64_t refined =
      ballN < 400 ? wlHashCount(g.adjacency, 0, hops + 2, g.localColors)
                  : wlHashWide(g.adjacency, 0, hops + 2, g.localColors);
  if (refinement != nullptr) {
    *refinement = refined;
  }
  const bool colored = !g.localColors.empty();
  const int n = static_cast<int>(g.atoms.size());
  const std::string cert =
      colored ? cage::canonicalCertificateColored(n, g.edges, g.localColors, 0)
               : cage::canonicalCertificateRooted(n, g.edges, 0);
  if (!cert.empty()) {
    out.method = "nauty";
    std::string tagged = cert;
    if (colored) {
      std::vector<int> sorted(g.localColors);
      std::sort(sorted.begin(), sorted.end());
      tagged += "|";
      for (std::size_t i = 0; i < sorted.size(); i++) {
        if (i == 0 || sorted[i] != sorted[i - 1]) {
          tagged += std::to_string(sorted[i]) + ":";
        }
        tagged += ".";
      }
    }
    out.key = hex(hashString(tagged));
  } else {
    out.method = "wl";
    out.key = hex(refined);
  }
  memo.ready = true;
  memo.hops = hops;
  memo.adjacency = g.adjacency;
  memo.colors = g.localColors;
  memo.key = out;
  memo.refined = refined;
  return out;
}

LocalKey flatNeighbourKey(const Rows &rows, int atom, int hops, const std::vector<int> &colors,
                         std::uint64_t *refined);

LocalKey localKeyRefined(const Rows &rows, int atom, int hops,
                        const std::vector<int> &colors, std::uint64_t *refined) {
  // A large crystal on one thread repeats one small ball. The flat image
  // for that ball lives in its own unit so this search, and the team that
  // shares it, stay the code they were tuned as. A team does not build the
  // packed rows: that walk would finish before the team starts.
#ifdef SEAMS_HAS_OPENMP
  const bool flat = static_cast<int>(rows.size()) > 4096 && omp_get_max_threads() <= 1 &&
                    !omp_in_parallel();
#else
  const bool flat = static_cast<int>(rows.size()) > 4096;
#endif
  if (flat) {
    return flatNeighbourKey(rows, atom, hops, colors, refined);
  }
  thread_local LocalGraph g;
  fillLocalGraph(g, rows, atom, hops, colors);
  return localKeyFromGraph(g, hops, refined);
}

LocalKey localKey(const Rows &rows, int atom, int hops, const std::vector<int> &colors) {
  return localKeyRefined(rows, atom, hops, colors, nullptr);
}

bool ballCoversAll(const Rows &rows, int hops) {
  const int n = static_cast<int>(rows.size());
  if (n <= 1 || hops < 1) {
    return n <= 1;
  }
  std::vector<int> seen(static_cast<std::size_t>(n), 0);
  std::vector<int> depth(static_cast<std::size_t>(n), 0);
  std::vector<int> order;
  order.reserve(static_cast<std::size_t>(n));
  for (int atom = 0; atom < n; atom++) {
    std::fill(seen.begin(), seen.end(), 0);
    order.clear();
    seen[static_cast<std::size_t>(atom)] = 1;
    depth[static_cast<std::size_t>(atom)] = 0;
    order.push_back(atom);
    for (std::size_t head = 0; head < order.size(); head++) {
      const int v = order[head];
      if (depth[static_cast<std::size_t>(v)] == hops) {
        continue;
      }
      const auto &row = rows[static_cast<std::size_t>(v)];
      for (std::size_t m = 1; m < row.size(); m++) {
        const int u = row[m];
        if (u < 0 || u >= n || seen[static_cast<std::size_t>(u)] != 0) {
          continue;
        }
        depth[static_cast<std::size_t>(u)] = depth[static_cast<std::size_t>(v)] + 1;
        seen[static_cast<std::size_t>(u)] = 1;
        order.push_back(u);
      }
    }
    if (static_cast<int>(order.size()) != n) {
      return false;
    }
  }
  return true;
}

std::vector<std::vector<int>> adjacencyInBall(const Rows &rows) {
  const int n = static_cast<int>(rows.size());
  std::vector<std::vector<int>> adj(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    const auto &row = rows[static_cast<std::size_t>(i)];
    auto &out = adj[static_cast<std::size_t>(i)];
    out.reserve(row.empty() ? 0 : row.size() - 1);
    for (std::size_t m = 1; m < row.size(); m++) {
      const int u = row[m];
      if (u < 0 || u >= n || u == i) {
        continue;
      }
      out.push_back(u);
    }
  }
  return adj;
}

FrameFingerprint fingerprintChoose(const Rows &rows, int hops, int maxRingSize,
                                   const std::vector<int> &colors, int which) {
  FrameFingerprint out;
  out.hops = hops;
  const int n = static_cast<int>(rows.size());
  const bool colored = colors.size() == rows.size() && !rows.empty();
  out.colored = colored;
  out.atomKeys.reserve(static_cast<std::size_t>(n));
  // The frame key is built from the refinement hash in every build so that
  // a nauty host and a plain host name the same frame the same way; the
  // per-atom keys carry the exact certificate when it exists.
  std::vector<std::uint64_t> wl(static_cast<std::size_t>(n));
  // A hop ball that is the whole frame is the same graph from every atom.
  // The refinement only moves the distinguished root. Ice is far larger
  // than this and keeps the per-atom ball.
  // which: 0 auto, 1 shared ball, 2 one ball per atom.
#ifndef SEAMS_HAS_NAUTY
  int maxDeg = 0;
  if (which != 2 && n > 1 && n <= 4096) {
    for (const auto &row : rows) {
      maxDeg = std::max(maxDeg, static_cast<int>(row.size()) - 1);
    }
  }
  // Only a very dense frame can have every hop ball equal the whole graph.
  const bool shared =
      which == 1 || (which != 2 && maxDeg > 100 && ballCoversAll(rows, hops));
#else
  const bool shared = which == 1;
#endif
  if (shared) {
    const auto adj = adjacencyInBall(rows);
    const std::vector<int> none;
    const std::vector<int> &paint = colored ? colors : none;
    out.atomKeys.resize(static_cast<std::size_t>(n));
    out.method = "wl";
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 128 && !omp_in_parallel())
#endif
    for (int i = 0; i < n; i++) {
      const std::uint64_t refined = wlHashWide(adj, i, hops + 2, paint);
      wl[static_cast<std::size_t>(i)] = refined;
      out.atomKeys[static_cast<std::size_t>(i)] = hex(refined);
    }
    for (const auto &key : out.atomKeys) {
      out.classes[key] += 1;
    }
  } else {
    out.atomKeys.resize(static_cast<std::size_t>(n));
    // Each atom's ball is independent. The class histogram is filled
    // afterwards so the counts stay in atom order, and the recorded
    // method is the last atom's, as a serial walk leaves it.
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256 && !omp_in_parallel() && \
                                              omp_get_max_threads() > 1)
#endif
    for (int i = 0; i < n; i++) {
      std::uint64_t refined = 0;
      LocalKey lk = localKeyRefined(rows, i, hops, colors, &refined);
      if (i + 1 == n) {
        out.method = lk.method;
      }
      out.atomKeys[static_cast<std::size_t>(i)] = std::move(lk.key);
      wl[static_cast<std::size_t>(i)] = refined;
    }
    for (const auto &key : out.atomKeys) {
      out.classes[key] += 1;
    }
  }
  out.wlAtom = wl;
  std::sort(wl.begin(), wl.end());
  std::uint64_t h = mix(kOffset, static_cast<std::uint64_t>(n));
  for (std::uint64_t v : wl) {
    h = mix(h, v);
  }
  fillRingCensus(out.ringCensus, rows, maxRingSize);
  for (int s = 0; s <= maxRingSize; s++) {
    h = mix(h, static_cast<std::uint64_t>(out.ringCensus[static_cast<std::size_t>(s)]));
  }
  out.key = hex(h);
  return out;
}

FrameFingerprint fingerprint(const Rows &rows, int hops, int maxRingSize,
                             const std::vector<int> &colors) {
  return fingerprintChoose(rows, hops, maxRingSize, colors, 0);
}

FrameFingerprint fingerprintPaths(const Rows &rows, int hops, int maxRingSize,
                                  const std::vector<int> &colors, int which) {
  return fingerprintChoose(rows, hops, maxRingSize, colors, which);
}

namespace {

void rebuildFrameKey(FrameFingerprint &fp, int maxRingSize) {
  std::vector<std::uint64_t> wl = fp.wlAtom;
  std::sort(wl.begin(), wl.end());
  std::uint64_t h = mix(kOffset, static_cast<std::uint64_t>(fp.atomKeys.size()));
  for (std::uint64_t v : wl) {
    h = mix(h, v);
  }
  for (int s = 0; s <= maxRingSize && s < static_cast<int>(fp.ringCensus.size());
       s++) {
    h = mix(h, static_cast<std::uint64_t>(fp.ringCensus[static_cast<std::size_t>(s)]));
  }
  fp.key = hex(h);
}

} // namespace

FrameFingerprint incrementalFingerprint(const FrameFingerprint &prev,
                                        const Rows &rows,
                                        const std::vector<int> &dirtyAtoms,
                                        int hops, int maxRingSize,
                                        const std::vector<int> &colors) {
  const int n = static_cast<int>(rows.size());
  if (n == 0 || static_cast<int>(prev.atomKeys.size()) != n ||
      prev.hops != hops || prev.wlAtom.size() != rows.size()) {
    return fingerprint(rows, hops, maxRingSize, colors);
  }
  std::unordered_set<int> dirty;
  for (int a : dirtyAtoms) {
    if (a < 0 || a >= n) {
      continue;
    }
    for (int b : hopNeighbourhood(rows, a, hops)) {
      dirty.insert(b);
    }
  }
  FrameFingerprint out = prev;
  out.hops = hops;
  out.colored = colors.size() == rows.size() && !rows.empty();
  out.classes.clear();
  for (int i : dirty) {
    std::uint64_t refined = 0;
    LocalKey lk = localKeyRefined(rows, i, hops, colors, &refined);
    out.method = lk.method;
    out.atomKeys[static_cast<std::size_t>(i)] = std::move(lk.key);
    out.wlAtom[static_cast<std::size_t>(i)] = refined;
  }
  for (const auto &k : out.atomKeys) {
    out.classes[k] += 1;
  }
  fillRingCensus(out.ringCensus, rows, maxRingSize);
  rebuildFrameKey(out, maxRingSize);
  return out;
}

void addToLibrary(KeyLibrary &lib, const FrameFingerprint &fp, const std::string &label) {
  if (lib.labelOf.empty()) {
    lib.method = fp.method;
    lib.hops = fp.hops;
    lib.colored = fp.colored;
  } else if (lib.method != fp.method || lib.hops != fp.hops || lib.colored != fp.colored) {
    throw std::invalid_argument("key library mixes methods, hop counts or colorings");
  }
  for (const auto &kv : fp.classes) {
    auto it = lib.labelOf.find(kv.first);
    if (it == lib.labelOf.end()) {
      lib.labelOf.emplace(kv.first, label);
      continue;
    }
    // a key seen under several labels names all of them, sorted and
    // joined by '|', so the census says which references share an
    // environment
    std::vector<std::string> parts;
    std::string cur;
    for (char ch : it->second) {
      if (ch == '|') {
        parts.push_back(cur);
        cur.clear();
      } else {
        cur.push_back(ch);
      }
    }
    parts.push_back(cur);
    if (std::find(parts.begin(), parts.end(), label) != parts.end()) {
      continue;
    }
    parts.push_back(label);
    std::sort(parts.begin(), parts.end());
    std::string joined;
    for (std::size_t i = 0; i < parts.size(); i++) {
      joined += (i ? "|" : "") + parts[i];
    }
    it->second = joined;
  }
}

std::string writeLibrary(const KeyLibrary &lib) {
  std::ostringstream out;
  out << "# method " << lib.method << " hops " << lib.hops << " colors "
      << (lib.colored ? 1 : 0) << "\n";
  for (const auto &kv : lib.labelOf) {
    out << kv.first << " " << kv.second << "\n";
  }
  return out.str();
}

KeyLibrary readLibrary(const std::string &text) {
  KeyLibrary lib;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    if (line[0] == '#') {
      std::istringstream head(line.substr(1));
      std::string word;
      while (head >> word) {
        if (word == "method") {
          head >> lib.method;
        } else if (word == "hops") {
          head >> lib.hops;
        } else if (word == "colors" || word == "colours") {
          // `colours` is the 2.13 spelling. The library has no
          // deprecation-warning channel, so the alias is silent.
          int flag = 0;
          head >> flag;
          lib.colored = flag != 0;
        }
      }
      continue;
    }
    std::istringstream row(line);
    std::string key;
    std::string label;
    if (row >> key >> label) {
      lib.labelOf[key] = label;
    }
  }
  return lib;
}

LibraryMatch matchLibrary(const FrameFingerprint &fp, const KeyLibrary &lib) {
  if (fp.method != lib.method || fp.hops != lib.hops || fp.colored != lib.colored) {
    throw std::invalid_argument(
        "fingerprint and library differ in method, hops or coloring");
  }
  LibraryMatch out;
  out.labels.reserve(fp.atomKeys.size());
  out.depth.reserve(fp.atomKeys.size());
  for (const auto &key : fp.atomKeys) {
    const auto it = lib.labelOf.find(key);
    const std::string label = it == lib.labelOf.end() ? std::string() : it->second;
    out.counts[label] += 1;
    if (!label.empty()) {
      ++out.matched;
    }
    out.depth.push_back(label.empty() ? 0 : lib.hops);
    out.labels.push_back(label);
  }
  return out;
}

LibraryMatch matchLibraries(const Rows &rows, const std::vector<KeyLibrary> &libs,
                            int maxRingSize, const std::vector<int> &colors) {
  if (libs.empty()) {
    throw std::invalid_argument("no key library to match against");
  }
  // deepest first, so the most specific name wins
  std::vector<const KeyLibrary *> order;
  for (const auto &lib : libs) {
    order.push_back(&lib);
  }
  std::sort(order.begin(), order.end(),
            [](const KeyLibrary *a, const KeyLibrary *b) { return a->hops > b->hops; });
  for (std::size_t i = 1; i < order.size(); ++i) {
    if (order[i - 1]->hops == order[i]->hops) {
      throw std::invalid_argument("two key libraries share a hop count");
    }
  }
  LibraryMatch out;
  out.labels.assign(rows.size(), std::string());
  out.depth.assign(rows.size(), 0);
  for (const KeyLibrary *lib : order) {
    if (lib->colored != (colors.size() == rows.size() && !rows.empty())) {
      throw std::invalid_argument(
          "fingerprint and library differ in method or coloring");
    }
    for (std::size_t a = 0; a < rows.size(); ++a) {
      if (!out.labels[a].empty()) {
        continue;
      }
      const LocalKey lk = localKey(rows, static_cast<int>(a), lib->hops, colors);
      if (lk.method != lib->method) {
        throw std::invalid_argument(
            "fingerprint and library differ in method or coloring");
      }
      const auto it = lib->labelOf.find(lk.key);
      if (it != lib->labelOf.end()) {
        out.labels[a] = it->second;
        out.depth[a] = lib->hops;
      }
    }
  }
  for (const auto &label : out.labels) {
    out.counts[label] += 1;
    if (!label.empty()) {
      ++out.matched;
    }
  }
  return out;
}

} // namespace topo
