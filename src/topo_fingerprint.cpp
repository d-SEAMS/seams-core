//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------
#include <topo_fingerprint.hpp>

#include <cage_canon.hpp>
#include <franzblau.hpp>

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <unordered_map>

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
  std::ostringstream out;
  out << std::hex << std::setw(16) << std::setfill('0') << value;
  return out.str();
}

std::vector<int> hopNeighbourhood(const Rows &rows, int atom, int hops) {
  std::vector<int> order;
  if (atom < 0 || atom >= static_cast<int>(rows.size())) {
    return order;
  }
  std::unordered_map<int, int> depth;
  depth[atom] = 0;
  order.push_back(atom);
  std::size_t head = 0;
  while (head < order.size()) {
    const int v = order[head++];
    const int d = depth[v];
    if (d == hops) {
      continue;
    }
    const auto &row = rows[static_cast<std::size_t>(v)];
    for (std::size_t m = 1; m < row.size(); m++) {
      const int u = row[m];
      if (u < 0 || u >= static_cast<int>(rows.size()) || depth.count(u)) {
        continue;
      }
      depth[u] = d + 1;
      order.push_back(u);
    }
  }
  return order;
}

std::uint64_t wlHash(const std::vector<std::vector<int>> &adjacency, int root,
                     int rounds, const std::vector<int> &colours) {
  const int n = static_cast<int>(adjacency.size());
  const bool coloured = static_cast<int>(colours.size()) == n;
  std::vector<std::uint64_t> colour(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    const std::uint64_t deg = static_cast<std::uint64_t>(adjacency[static_cast<std::size_t>(i)].size());
    std::uint64_t h = mix(kOffset, deg + (i == root ? 0x9e3779b97f4a7c15ULL : 0ULL));
    if (coloured) {
      h = mix(h, static_cast<std::uint64_t>(static_cast<std::int64_t>(colours[static_cast<std::size_t>(i)])) + 1ULL);
    }
    colour[static_cast<std::size_t>(i)] = h;
  }
  std::vector<std::uint64_t> next(static_cast<std::size_t>(n));
  std::vector<std::uint64_t> bag;
  for (int r = 0; r < rounds; r++) {
    for (int i = 0; i < n; i++) {
      bag.clear();
      for (int j : adjacency[static_cast<std::size_t>(i)]) {
        bag.push_back(colour[static_cast<std::size_t>(j)]);
      }
      std::sort(bag.begin(), bag.end());
      std::uint64_t h = mix(kOffset, colour[static_cast<std::size_t>(i)]);
      for (std::uint64_t c : bag) {
        h = mix(h, c);
      }
      next[static_cast<std::size_t>(i)] = h;
    }
    colour.swap(next);
  }
  std::vector<std::uint64_t> all(colour.begin(), colour.end());
  std::sort(all.begin(), all.end());
  std::uint64_t h = mix(kOffset, static_cast<std::uint64_t>(n));
  if (root >= 0 && root < n) {
    h = mix(h, colour[static_cast<std::size_t>(root)]);
  }
  for (std::uint64_t c : all) {
    h = mix(h, c);
  }
  return h;
}

LocalKey localKey(const Rows &rows, int atom, int hops, const std::vector<int> &colours) {
  LocalKey out;
  const bool coloured = colours.size() == rows.size() && !rows.empty();
  const std::vector<int> atoms = hopNeighbourhood(rows, atom, hops);
  if (atoms.empty()) {
    out.method = "wl";
    out.key = hex(0);
    return out;
  }
  std::unordered_map<int, int> local;
  for (std::size_t i = 0; i < atoms.size(); i++) {
    local[atoms[i]] = static_cast<int>(i);
  }
  const int n = static_cast<int>(atoms.size());
  std::vector<std::vector<int>> adjacency(static_cast<std::size_t>(n));
  std::vector<std::pair<int, int>> edges;
  for (int i = 0; i < n; i++) {
    const auto &row = rows[static_cast<std::size_t>(atoms[static_cast<std::size_t>(i)])];
    for (std::size_t m = 1; m < row.size(); m++) {
      const auto it = local.find(row[m]);
      if (it == local.end() || it->second == i) {
        continue;
      }
      adjacency[static_cast<std::size_t>(i)].push_back(it->second);
      if (i < it->second) {
        edges.emplace_back(i, it->second);
      }
    }
    std::sort(adjacency[static_cast<std::size_t>(i)].begin(),
              adjacency[static_cast<std::size_t>(i)].end());
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  out.vertices = n;
  out.edges = static_cast<int>(edges.size());
  std::vector<int> localColours;
  if (coloured) {
    localColours.reserve(atoms.size());
    for (int a : atoms) {
      localColours.push_back(colours[static_cast<std::size_t>(a)]);
    }
  }
  const std::string cert = coloured
                               ? cage::canonicalCertificateColoured(n, edges, localColours, 0)
                               : cage::canonicalCertificateRooted(n, edges, 0);
  if (!cert.empty()) {
    out.method = "nauty";
    // the certificate is canonical for the coloured graph given the colour
    // cells; the cells themselves (colour value and size, sorted) complete
    // the class and do not depend on the input numbering
    std::string tagged = cert;
    if (coloured) {
      std::vector<int> sorted(localColours);
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
    return out;
  }
  out.method = "wl";
  out.key = hex(wlHash(adjacency, 0, hops + 2, localColours));
  return out;
}

FrameFingerprint fingerprint(const Rows &rows, int hops, int maxRingSize,
                             const std::vector<int> &colours) {
  FrameFingerprint out;
  out.hops = hops;
  const int n = static_cast<int>(rows.size());
  const bool coloured = colours.size() == rows.size() && !rows.empty();
  out.atomKeys.reserve(static_cast<std::size_t>(n));
  // The frame key is built from the refinement hash in every build so that
  // a nauty host and a plain host name the same frame the same way; the
  // per-atom keys carry the exact certificate when it exists.
  std::vector<std::uint64_t> wl(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    LocalKey lk = localKey(rows, i, hops, colours);
    out.method = lk.method;
    out.classes[lk.key] += 1;
    out.atomKeys.push_back(std::move(lk.key));
    const std::vector<int> atoms = hopNeighbourhood(rows, i, hops);
    std::unordered_map<int, int> local;
    for (std::size_t k = 0; k < atoms.size(); k++) {
      local[atoms[k]] = static_cast<int>(k);
    }
    std::vector<std::vector<int>> adjacency(atoms.size());
    for (std::size_t k = 0; k < atoms.size(); k++) {
      const auto &row = rows[static_cast<std::size_t>(atoms[k])];
      for (std::size_t m = 1; m < row.size(); m++) {
        const auto it = local.find(row[m]);
        if (it != local.end() && it->second != static_cast<int>(k)) {
          adjacency[k].push_back(it->second);
        }
      }
      std::sort(adjacency[k].begin(), adjacency[k].end());
    }
    std::vector<int> localColours;
    if (coloured) {
      for (int a : atoms) {
        localColours.push_back(colours[static_cast<std::size_t>(a)]);
      }
    }
    wl[static_cast<std::size_t>(i)] =
        atoms.empty() ? 0 : wlHash(adjacency, 0, hops + 2, localColours);
  }
  std::sort(wl.begin(), wl.end());
  std::uint64_t h = mix(kOffset, static_cast<std::uint64_t>(n));
  for (std::uint64_t v : wl) {
    h = mix(h, v);
  }
  out.ringCensus.assign(static_cast<std::size_t>(maxRingSize) + 1, 0);
  if (n > 0 && maxRingSize >= 3) {
    for (const auto &ring : primitive::ringNetwork(rows, maxRingSize)) {
      if (ring.size() <= static_cast<std::size_t>(maxRingSize)) {
        out.ringCensus[ring.size()] += 1;
      }
    }
  }
  for (int s = 0; s <= maxRingSize; s++) {
    h = mix(h, static_cast<std::uint64_t>(out.ringCensus[static_cast<std::size_t>(s)]));
  }
  out.key = hex(h);
  return out;
}

} // namespace topo
