//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <cage_canon.hpp>
#include <cage_enum.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

std::string trimmed(std::string_view in) {
  size_t a = 0;
  size_t b = in.size();
  while (a < b && std::isspace(static_cast<unsigned char>(in[a]))) {
    ++a;
  }
  while (b > a && std::isspace(static_cast<unsigned char>(in[b - 1]))) {
    --b;
  }
  return std::string(in.substr(a, b - a));
}

std::string lowerCopy(std::string s) {
  for (char &c : s) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return s;
}

const std::unordered_map<std::string, std::string> &namedTable() {
  static const std::unordered_map<std::string, std::string> table = {
      {"sodalite", "4:6,6:8"},
      {"alpha", "4:6,6:8,8:6"},
      {"512", "5:12"},
      {"51262", "5:12,6:2"},
      {"51264", "5:12,6:4"},
      {"51268", "5:12,6:8"},
      {"sh", "4:3,5:6,6:3"},
      {"sH", "4:3,5:6,6:3"},
      {"hc", "4:6,6:2"},
      {"ddc", "6:7"},
  };
  return table;
}

using Edge = std::uint64_t;

Edge packEdge(int a, int b) {
  const auto lo = static_cast<std::uint32_t>(std::min(a, b));
  const auto hi = static_cast<std::uint32_t>(std::max(a, b));
  return (static_cast<std::uint64_t>(lo) << 32) | hi;
}

std::vector<Edge> edgesOfRing(const std::vector<int> &ring) {
  std::vector<Edge> edges;
  const size_t n = ring.size();
  if (n < 2) {
    return edges;
  }
  edges.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    const int a = ring[i];
    const int b = ring[(i + 1) % n];
    if (a != b) {
      edges.push_back(packEdge(a, b));
    }
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  return edges;
}

std::vector<int> uniqueVertices(const std::vector<std::vector<int>> &rings,
                                const std::vector<int> &faces) {
  std::set<int> verts;
  for (const int fi : faces) {
    if (fi < 0 || fi >= static_cast<int>(rings.size())) {
      continue;
    }
    for (const int a : rings[static_cast<size_t>(fi)]) {
      verts.insert(a);
    }
  }
  return {verts.begin(), verts.end()};
}

bool countsEqual(const std::map<int, int> &have, const cage::Signature &sig) {
  if (have.size() != sig.counts.size()) {
    return false;
  }
  return have == sig.counts;
}

bool wouldOverflowSize(const std::map<int, int> &have, const cage::Signature &sig,
                       int size) {
  const auto need = sig.counts.find(size);
  if (need == sig.counts.end()) {
    return true;
  }
  const auto got = have.find(size);
  const int n = got == have.end() ? 0 : got->second;
  return n + 1 > need->second;
}

struct Search {
  const std::vector<std::vector<int>> &rings;
  const cage::Signature &sig;
  std::vector<int> cand;
  std::vector<std::vector<Edge>> edges;
  std::unordered_map<Edge, std::vector<int>> ringsOf;
  std::vector<char> used;
  std::vector<int> chosen;
  std::map<int, int> have;
  std::unordered_map<Edge, int> edgeUse;
  std::set<std::vector<int>> seenVerts;
  std::vector<cage::FoundCage> found;
  int seed = 0;
  bool allowIncomplete = false;
  int minFaces = 1;

  explicit Search(const std::vector<std::vector<int>> &allRings,
                  const cage::Signature &signature)
      : rings(allRings), sig(signature) {
    const int nRings = static_cast<int>(rings.size());
    used.assign(static_cast<size_t>(nRings), 0);
    edges.resize(static_cast<size_t>(nRings));
    for (int i = 0; i < nRings; ++i) {
      const int sz = static_cast<int>(rings[static_cast<size_t>(i)].size());
      if (!sig.containsSize(sz)) {
        continue;
      }
      cand.push_back(i);
      edges[static_cast<size_t>(i)] = edgesOfRing(rings[static_cast<size_t>(i)]);
      for (const Edge e : edges[static_cast<size_t>(i)]) {
        ringsOf[e].push_back(i);
      }
    }
  }

  bool edgeWouldOverflow(int ring) const {
    for (const Edge e : edges[static_cast<size_t>(ring)]) {
      const auto it = edgeUse.find(e);
      const int n = it == edgeUse.end() ? 0 : it->second;
      if (n >= 2) {
        return true;
      }
    }
    return false;
  }

  void addRing(int ring) {
    used[static_cast<size_t>(ring)] = 1;
    chosen.push_back(ring);
    const int sz = static_cast<int>(rings[static_cast<size_t>(ring)].size());
    have[sz] += 1;
    for (const Edge e : edges[static_cast<size_t>(ring)]) {
      edgeUse[e] += 1;
    }
  }

  void undoRing() {
    const int ring = chosen.back();
    chosen.pop_back();
    used[static_cast<size_t>(ring)] = 0;
    const int sz = static_cast<int>(rings[static_cast<size_t>(ring)].size());
    auto hit = have.find(sz);
    if (hit != have.end()) {
      hit->second -= 1;
      if (hit->second <= 0) {
        have.erase(hit);
      }
    }
    for (const Edge e : edges[static_cast<size_t>(ring)]) {
      auto eit = edgeUse.find(e);
      if (eit == edgeUse.end()) {
        continue;
      }
      eit->second -= 1;
      if (eit->second <= 0) {
        edgeUse.erase(eit);
      }
    }
  }

  bool canAdd(int ring) const {
    if (ring < seed || used[static_cast<size_t>(ring)]) {
      return false;
    }
    const int sz = static_cast<int>(rings[static_cast<size_t>(ring)].size());
    if (wouldOverflowSize(have, sig, sz)) {
      return false;
    }
    return !edgeWouldOverflow(ring);
  }

  // Returns 0 on a dead branch, 1 when a forced ring was added, 2 when
  // every unsaturated edge still has a choice.
  int addForced() {
    for (const auto &kv : edgeUse) {
      if (kv.second != 1) {
        continue;
      }
      const auto rit = ringsOf.find(kv.first);
      if (rit == ringsOf.end()) {
        if (allowIncomplete) {
          continue;
        }
        return 0;
      }
      int only = -1;
      int nOpt = 0;
      for (const int r : rit->second) {
        if (!canAdd(r)) {
          continue;
        }
        ++nOpt;
        only = r;
        if (nOpt > 1) {
          break;
        }
      }
      if (nOpt == 0) {
        if (allowIncomplete) {
          continue;
        }
        return 0;
      }
      if (nOpt == 1) {
        addRing(only);
        return 1;
      }
    }
    return 2;
  }

  bool propagate() {
    for (;;) {
      const int st = addForced();
      if (st == 0) {
        return false;
      }
      if (st == 2) {
        return true;
      }
    }
  }

  bool allEdgesPaired() const {
    if (edgeUse.empty()) {
      return false;
    }
    for (const auto &kv : edgeUse) {
      if (kv.second != 2) {
        return false;
      }
    }
    return true;
  }

  void accept(bool closed) {
    std::vector<int> faces = chosen;
    std::sort(faces.begin(), faces.end());
    std::vector<int> verts = uniqueVertices(rings, faces);
    if (!seenVerts.insert(verts).second) {
      return;
    }
    int dangling = 0;
    for (const auto &kv : edgeUse) {
      if (kv.second == 1) {
        ++dangling;
      }
    }
    if (!closed && dangling == 0) {
      return;
    }
    std::vector<std::vector<int>> faceRings;
    faceRings.reserve(faces.size());
    for (const int fi : faces) {
      faceRings.push_back(rings[static_cast<size_t>(fi)]);
    }
    cage::FoundCage cage;
    cage.signature = sig;
    cage.faces = std::move(faces);
    cage.vertices = std::move(verts);
    cage.certificate = cage::canonicalCertificate(faceRings);
    cage.closed = closed;
    cage.danglingEdges = dangling;
    found.push_back(std::move(cage));
  }

  std::pair<Edge, std::vector<int>> branchEdge() const {
    Edge best = 0;
    std::vector<int> opts;
    bool haveBest = false;
    for (const auto &kv : edgeUse) {
      if (kv.second != 1) {
        continue;
      }
      const auto rit = ringsOf.find(kv.first);
      if (rit == ringsOf.end()) {
        continue;
      }
      std::vector<int> cur;
      for (const int r : rit->second) {
        if (canAdd(r)) {
          cur.push_back(r);
        }
      }
      if (!haveBest || cur.size() < opts.size()) {
        best = kv.first;
        opts = std::move(cur);
        haveBest = true;
        if (opts.size() <= 1) {
          break;
        }
      }
    }
    return {best, opts};
  }

  void restore(size_t mark) {
    while (chosen.size() > mark) {
      undoRing();
    }
  }

  void maybeAcceptIncomplete() {
    if (allowIncomplete && static_cast<int>(chosen.size()) >= minFaces &&
        !allEdgesPaired()) {
      accept(false);
    }
  }

  void search() {
    const size_t mark = chosen.size();
    if (!propagate()) {
      maybeAcceptIncomplete();
      restore(mark);
      return;
    }
    if (countsEqual(have, sig)) {
      if (allEdgesPaired()) {
        accept(true);
      }
      restore(mark);
      return;
    }
    const auto br = branchEdge();
    if (br.second.empty()) {
      maybeAcceptIncomplete();
      restore(mark);
      return;
    }
    for (const int r : br.second) {
      if (!canAdd(r)) {
        continue;
      }
      addRing(r);
      search();
      undoRing();
    }
    restore(mark);
  }

  void run() {
    for (const int r : cand) {
      seed = r;
      if (!allowIncomplete) {
        bool seedOk = true;
        for (const Edge e : edges[static_cast<size_t>(r)]) {
          const auto it = ringsOf.find(e);
          if (it == ringsOf.end() || it->second.size() < 2) {
            seedOk = false;
            break;
          }
        }
        if (!seedOk) {
          continue;
        }
      }
      addRing(r);
      search();
      undoRing();
    }
  }
};

} // namespace

namespace cage {

Signature Signature::parse(std::string_view spec) {
  const std::string raw = trimmed(spec);
  if (raw.empty()) {
    throw std::invalid_argument("empty cage signature");
  }
  std::string body = raw;
  Signature out;
  const std::string key = lowerCopy(raw);
  if (key == "hc") {
    out.kind = Signature::Kind::HexC;
    body = "4:6,6:2";
  } else if (key == "ddc") {
    out.kind = Signature::Kind::DoubleDiaC;
    body = "6:7";
  } else {
    const auto named = namedTable().find(key);
    if (named != namedTable().end()) {
      body = named->second;
    }
  }
  std::string token;
  std::istringstream in(body);
  while (std::getline(in, token, ',')) {
    const std::string part = trimmed(token);
    if (part.empty()) {
      continue;
    }
    const auto colon = part.find(':');
    if (colon == std::string::npos || colon == 0 || colon + 1 == part.size()) {
      throw std::invalid_argument("bad cage signature '" + raw + "'");
    }
    int size = 0;
    int count = 0;
    try {
      size = std::stoi(part.substr(0, colon));
      count = std::stoi(part.substr(colon + 1));
    } catch (const std::exception &) {
      throw std::invalid_argument("bad cage signature '" + raw + "'");
    }
    if (size <= 0 || count <= 0) {
      throw std::invalid_argument("cage signature sizes and counts must be positive");
    }
    out.counts[size] += count;
  }
  if (out.counts.empty()) {
    throw std::invalid_argument("empty cage signature");
  }
  return out;
}

std::string Signature::str() const {
  std::ostringstream out;
  bool first = true;
  for (const auto &kv : counts) {
    if (!first) {
      out << ',';
    }
    first = false;
    out << kv.first << ':' << kv.second;
  }
  return out.str();
}

int Signature::faceCount() const {
  int n = 0;
  for (const auto &kv : counts) {
    n += kv.second;
  }
  return n;
}

int Signature::maxRingSize() const {
  if (counts.empty()) {
    return 0;
  }
  return counts.rbegin()->first;
}

bool Signature::containsSize(int size) const {
  return counts.find(size) != counts.end();
}

bool isClosedPolyhedron(const std::vector<std::vector<int>> &rings,
                        const std::vector<int> &faces) {
  if (faces.empty()) {
    return false;
  }
  std::unordered_map<Edge, int> use;
  for (const int fi : faces) {
    if (fi < 0 || fi >= static_cast<int>(rings.size())) {
      return false;
    }
    const auto edges = edgesOfRing(rings[static_cast<size_t>(fi)]);
    if (edges.empty()) {
      return false;
    }
    for (const Edge e : edges) {
      use[e] += 1;
    }
  }
  for (const auto &kv : use) {
    if (kv.second != 2) {
      return false;
    }
  }
  return true;
}

std::vector<FoundCage> findBySignature(const std::vector<std::vector<int>> &rings,
                                       const Signature &signature) {
  if (signature.counts.empty() || signature.faceCount() <= 0) {
    return {};
  }
  Search search(rings, signature);
  search.run();
  return search.found;
}

std::vector<FoundCage>
findIncompleteBySignature(const std::vector<std::vector<int>> &rings,
                          const Signature &signature, int minFaces) {
  if (signature.counts.empty() || signature.faceCount() <= 0) {
    return {};
  }
  const int floor = minFaces > 0 ? minFaces
                                 : std::max(1, signature.faceCount() / 2);
  Search search(rings, signature);
  search.allowIncomplete = true;
  search.minFaces = floor;
  search.run();
  const auto closed = findBySignature(rings, signature);
  std::set<int> closedFaces;
  for (const auto &cl : closed) {
    closedFaces.insert(cl.faces.begin(), cl.faces.end());
  }
  std::vector<FoundCage> out;
  out.reserve(search.found.size());
  for (auto &c : search.found) {
    if (c.closed) {
      continue;
    }
    bool allInClosed = !closedFaces.empty();
    if (allInClosed) {
      for (const int f : c.faces) {
        if (closedFaces.find(f) == closedFaces.end()) {
          allInClosed = false;
          break;
        }
      }
    }
    if (!allInClosed) {
      out.push_back(std::move(c));
    }
  }
  return out;
}

namespace {

std::vector<FoundCage>
tumByType(const std::vector<std::vector<int>> &rings,
          const std::vector<std::vector<int>> &nList, cageType want,
          const Signature &signature) {
  std::vector<std::vector<int>> six;
  std::vector<int> origin;
  six.reserve(rings.size());
  origin.reserve(rings.size());
  for (int i = 0; i < static_cast<int>(rings.size()); ++i) {
    if (rings[static_cast<size_t>(i)].size() == 6) {
      six.push_back(rings[static_cast<size_t>(i)]);
      origin.push_back(i);
    }
  }
  std::vector<ring::strucType> ringType(six.size());
  std::vector<Cage> cageList;
  const auto listHC = ring::findHC(six, ringType, nList, cageList);
  if (want == cageType::DoubleDiaC) {
    ring::findDDC(six, ringType, listHC, cageList);
  }
  std::vector<FoundCage> out;
  for (const auto &c : cageList) {
    if (c.type != want) {
      continue;
    }
    FoundCage found;
    found.signature = signature;
    found.faces.reserve(c.rings.size());
    for (const int si : c.rings) {
      if (si >= 0 && si < static_cast<int>(origin.size())) {
        found.faces.push_back(origin[static_cast<size_t>(si)]);
      }
    }
    found.vertices = uniqueVertices(rings, found.faces);
    std::vector<std::vector<int>> faceRings;
    faceRings.reserve(found.faces.size());
    for (const int fi : found.faces) {
      faceRings.push_back(rings[static_cast<size_t>(fi)]);
    }
    found.certificate = canonicalCertificate(faceRings);
    out.push_back(std::move(found));
  }
  return out;
}

} // namespace

std::vector<FoundCage>
findBySignature(const std::vector<std::vector<int>> &rings,
                const std::vector<std::vector<int>> &nList,
                const Signature &signature) {
  if (signature.kind == Signature::Kind::HexC) {
    return tumByType(rings, nList, cageType::HexC, signature);
  }
  if (signature.kind == Signature::Kind::DoubleDiaC) {
    return tumByType(rings, nList, cageType::DoubleDiaC, signature);
  }
  return findBySignature(rings, signature);
}

} // namespace cage
