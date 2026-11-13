//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <cage_canon.hpp>
#include <cage_enum.hpp>
#include <ira_sofi.hpp>
#include <number_parse.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

// `lo` and `hi` are inclusive and use the coordinate units of `xyz`.
// The reference oxygen-oxygen contacts sit near 2.74–2.81 when those
// units are angstroms (mW 2.757, sI 2.771–2.805, sII 2.737–2.771).
bool waterFace(const std::vector<int> &ring,
               const std::vector<std::array<double, 3>> &xyz,
               const std::vector<double> *box, double lo, double hi) {
  const int n = static_cast<int>(ring.size());
  double len[3] = {0, 0, 0};
  const bool ortho = box != nullptr && box->size() >= 3;
  if (ortho) {
    len[0] = (*box)[0];
    len[1] = (*box)[1];
    len[2] = (*box)[2];
  }
  for (int i = 0; i < n; ++i) {
    const int ia = ring[static_cast<size_t>(i)];
    const int ib = ring[static_cast<size_t>((i + 1) % n)];
    if (ia < 0 || ib < 0 || static_cast<size_t>(ia) >= xyz.size() ||
        static_cast<size_t>(ib) >= xyz.size()) {
      return false;
    }
    double d[3];
    for (int k = 0; k < 3; ++k) {
      d[k] = xyz[static_cast<size_t>(ia)][static_cast<size_t>(k)] -
             xyz[static_cast<size_t>(ib)][static_cast<size_t>(k)];
      if (ortho && len[k] > 0) {
        d[k] -= len[k] * std::round(d[k] / len[k]);
      }
    }
    const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (r < lo || r > hi) {
      return false;
    }
  }
  return true;
}

struct FaceSetHash {
  size_t operator()(const std::vector<int> &v) const {
    size_t h = v.size();
    for (const int x : v) {
      h ^= static_cast<size_t>(x) + 0x9e3779b9u + (h << 6) + (h >> 2);
    }
    return h;
  }
};

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

bool facesTileSphere(const std::vector<std::vector<int>> &rings,
                     const std::vector<int> &faces) {
  std::unordered_map<int, std::vector<std::pair<int, int>>> wedges;
  std::set<Edge> edgeSet;
  for (const int fi : faces) {
    const auto &ring = rings[static_cast<size_t>(fi)];
    const size_t n = ring.size();
    for (size_t i = 0; i < n; ++i) {
      const int prev = ring[(i + n - 1) % n];
      const int next = ring[(i + 1) % n];
      wedges[ring[i]].emplace_back(prev, next);
      edgeSet.insert(packEdge(ring[i], next));
    }
  }
  const long v = static_cast<long>(wedges.size());
  const long e = static_cast<long>(edgeSet.size());
  const long f = static_cast<long>(faces.size());
  if (v - e + f != 2) {
    return false;
  }
  // With paired edges every neighbour of a vertex closes exactly two of its
  // wedges, so the wedges form cycles; one cycle must hold them all.
  for (const auto &kv : wedges) {
    const auto &w = kv.second;
    if (w.size() < 3) {
      return false;
    }
    std::vector<char> seen(w.size(), 0);
    size_t at = 0;
    int tail = w[0].second;
    size_t walked = 1;
    seen[0] = 1;
    for (;;) {
      size_t nextWedge = w.size();
      for (size_t k = 0; k < w.size(); ++k) {
        if (k != at && !seen[k] && (w[k].first == tail || w[k].second == tail)) {
          nextWedge = k;
          break;
        }
      }
      if (nextWedge == w.size()) {
        break;
      }
      seen[nextWedge] = 1;
      tail = w[nextWedge].first == tail ? w[nextWedge].second : w[nextWedge].first;
      at = nextWedge;
      ++walked;
    }
    if (walked != w.size()) {
      return false;
    }
  }
  return true;
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

  const std::vector<char> *admit = nullptr;

  explicit Search(const std::vector<std::vector<int>> &allRings,
                  const cage::Signature &signature,
                  const std::vector<char> *admitMask = nullptr)
      : rings(allRings), sig(signature), admit(admitMask) {
    const int nRings = static_cast<int>(rings.size());
    used.assign(static_cast<size_t>(nRings), 0);
    edges.resize(static_cast<size_t>(nRings));
    for (int i = 0; i < nRings; ++i) {
      if (admit != nullptr &&
          (static_cast<size_t>(i) >= admit->size() ||
           (*admit)[static_cast<size_t>(i)] == 0)) {
        continue;
      }
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
      if (allEdgesPaired() && facesTileSphere(rings, chosen)) {
        accept(true);
      } else {
        maybeAcceptIncomplete();
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

// Closed cages only, with Search's forced moves and branching on dense
// edge ids. Two cuts follow from facesTileSphere and drop no cage: a ring
// that shares two consecutive edges with a chosen ring leaves their common
// vertex in two faces, and a ring through a vertex whose chosen faces
// already close around it gives that vertex a second cycle. Up to 32
// signatures share one search: a face set grows while its census fits in
// one of them, and a closed set is a cage of the signature it equals.
// A ring is a simple cycle when its vertices are unique. Sorting a handful
// of ids is the whole test; a node-based set made the cage index allocate
// once per face.
bool distinctVertices(const std::vector<int> &ring) {
  const int sz = static_cast<int>(ring.size());
  if (sz <= 1) {
    return true;
  }
  int stack[32];
  int *tmp = stack;
  std::vector<int> heap;
  if (sz > 32) {
    heap = ring;
    tmp = heap.data();
  } else {
    std::copy(ring.begin(), ring.end(), tmp);
  }
  std::sort(tmp, tmp + sz);
  for (int i = 1; i < sz; ++i) {
    if (tmp[i] == tmp[i - 1]) {
      return false;
    }
  }
  return true;
}

// One allocation for every row. A vector per ring was most of the cage
// index, and the search only walks the rows.
struct IdRows {
  std::vector<int> off;
  std::vector<int> col;
  struct View {
    const int *b = nullptr;
    const int *e = nullptr;
    const int *begin() const { return b; }
    const int *end() const { return e; }
    std::size_t size() const { return static_cast<std::size_t>(e - b); }
    bool empty() const { return b == e; }
  };
  View operator[](std::size_t i) const {
    return {col.data() + off[i], col.data() + off[i + 1]};
  }
  void reset(int nRows) {
    off.assign(static_cast<std::size_t>(nRows) + 1, 0);
    col.clear();
  }
};

struct ClosedSearch {
  // Edge and ring tables are read by every worker. Sharing them is what
  // lets a team search seeds without copying the index onto each thread.
  struct Spine {
    std::vector<std::vector<int>> need;
    std::vector<int> facesOf;
    std::vector<int> maxVerts;
    std::vector<int> maxEdges;
    std::vector<int> cand;
    std::vector<int> ringSize;
    IdRows ringEdges;
    IdRows ringWedges;
    IdRows ringVerts;
    IdRows ringsOfEdge;
    std::vector<std::array<int, 2>> edgeEnds;
  };
  std::shared_ptr<Spine> spine;
  const std::vector<std::vector<int>> &rings;
  const std::vector<cage::Signature> &sigs;
  std::vector<std::vector<int>> &need;
  std::vector<int> &facesOf;
  std::vector<int> &maxVerts;
  std::vector<int> &maxEdges;
  std::vector<int> &cand;
  std::vector<int> &ringSize;
  IdRows &ringEdges;
  IdRows &ringWedges;
  IdRows &ringVerts;
  IdRows &ringsOfEdge;
  std::vector<std::array<int, 2>> &edgeEnds;
  int maxFaces = 0;
  std::vector<int> have;
  // room[s]: the signatures with space for one more ring of size s.
  // fits[d]: the signatures whose census holds the first d chosen rings.
  std::vector<unsigned> room;
  std::vector<unsigned> fits;
  std::vector<unsigned char> edgeUse;
  std::vector<unsigned char> wedgeUse;
  std::vector<int> half;
  std::vector<int> halfPos;
  std::vector<int> openAt;
  std::vector<int> facesAt;
  std::vector<char> used;
  std::vector<int> chosen;
  std::vector<std::vector<int>> optsAt;
  std::vector<std::vector<cage::FoundCage>> found;
  std::vector<int> branchOverride;
  // Face sets already expanded from this seed. The edges and vertices of a
  // partial cage are a function of the set, so a second route to the same
  // faces has the same continuations.
  std::unordered_set<std::vector<int>, FaceSetHash> seenFaces;
  bool trivalent = true;
  int seed = 0;
  int nVerts = 0;
  int nEdges = 0;

  ClosedSearch(const std::vector<std::vector<int>> &allRings,
               const std::vector<cage::Signature> &signatures,
               const std::vector<char> *admit = nullptr,
               const std::vector<std::array<double, 3>> *xyz = nullptr,
               const std::vector<double> *box = nullptr,
               const cage::EdgeWindow &window = {})
      : spine(std::make_shared<Spine>()), rings(allRings), sigs(signatures),
        need(spine->need), facesOf(spine->facesOf), maxVerts(spine->maxVerts),
        maxEdges(spine->maxEdges), cand(spine->cand), ringSize(spine->ringSize),
        ringEdges(spine->ringEdges), ringWedges(spine->ringWedges),
        ringVerts(spine->ringVerts), ringsOfEdge(spine->ringsOfEdge),
        edgeEnds(spine->edgeEnds) {
    const size_t nSig = sigs.size();
    int maxSize = 0;
    for (const auto &s : sigs) {
      maxSize = std::max(maxSize, s.maxRingSize());
    }
    need.assign(nSig, std::vector<int>(static_cast<size_t>(maxSize) + 1, 0));
    facesOf.assign(nSig, 0);
    maxVerts.assign(nSig, 0);
    maxEdges.assign(nSig, 0);
    std::vector<char> inSome(static_cast<size_t>(maxSize) + 1, 0);
    for (size_t k = 0; k < nSig; ++k) {
      int sides = 0;
      for (const auto &kv : sigs[k].counts) {
        need[k][static_cast<size_t>(kv.first)] = kv.second;
        inSome[static_cast<size_t>(kv.first)] = 1;
        sides += kv.first * kv.second;
      }
      facesOf[k] = sigs[k].faceCount();
      maxFaces = std::max(maxFaces, facesOf[k]);
      // A sphere with every edge in two faces has E = sides / 2 and
      // V = 2 - F + E.
      maxEdges[k] = sides / 2;
      maxVerts[k] = 2 - facesOf[k] + maxEdges[k];
      if (2 * maxEdges[k] != 3 * maxVerts[k]) {
        trivalent = false;
      }
    }
    have.assign(static_cast<size_t>(maxSize) + 1, 0);
    room.assign(have.size(), 0);
    for (int s = 0; s <= maxSize; ++s) {
      updateRoom(s);
    }
    fits.assign(static_cast<size_t>(maxFaces) + 1, 0);
    fits[0] = nSig >= 32 ? ~0u : (1u << nSig) - 1u;
    const int nRings = static_cast<int>(rings.size());
    used.assign(static_cast<size_t>(nRings), 0);
    ringSize.assign(static_cast<size_t>(nRings), 0);
    std::vector<char> keep(static_cast<size_t>(nRings), 0);
    int indexTeam = 1;
#ifdef SEAMS_HAS_OPENMP
    if (!omp_in_parallel()) {
      indexTeam = std::max(1, omp_get_max_threads());
    }
#endif
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel for schedule(static) if (nRings >= 128 && indexTeam > 1)
#endif
    for (int i = 0; i < nRings; ++i) {
      const auto &ring = rings[static_cast<size_t>(i)];
      const int sz = static_cast<int>(ring.size());
      if (admit != nullptr &&
          (static_cast<size_t>(i) >= admit->size() ||
           (*admit)[static_cast<size_t>(i)] == 0)) {
        continue;
      }
      if (sz < 3 || sz > maxSize || inSome[static_cast<size_t>(sz)] == 0) {
        continue;
      }
      if (window.enabled && xyz != nullptr &&
          !waterFace(ring, *xyz, box, window.lo, window.hi)) {
        continue;
      }
      if (!distinctVertices(ring)) {
        continue;
      }
      keep[static_cast<size_t>(i)] = 1;
    }
    int slots = 0;
    int maxV = -1;
    cand.reserve(static_cast<size_t>(nRings));
    for (int i = 0; i < nRings; ++i) {
      if (keep[static_cast<size_t>(i)] == 0) {
        continue;
      }
      cand.push_back(i);
      const int sz = static_cast<int>(rings[static_cast<size_t>(i)].size());
      ringSize[static_cast<size_t>(i)] = sz;
      slots += sz;
      for (const int v : rings[static_cast<size_t>(i)]) {
        maxV = std::max(maxV, v);
      }
    }
    auto lay = [&](IdRows &rows) {
      rows.off.assign(static_cast<size_t>(nRings) + 1, 0);
      int cursor = 0;
      for (int i = 0; i < nRings; ++i) {
        rows.off[static_cast<size_t>(i)] = cursor;
        if (keep[static_cast<size_t>(i)] != 0) {
          cursor += ringSize[static_cast<size_t>(i)];
        }
      }
      rows.off[static_cast<size_t>(nRings)] = cursor;
      rows.col.assign(static_cast<size_t>(cursor), 0);
    };
    lay(ringEdges);
    lay(ringWedges);
    lay(ringVerts);

    // Open-addressed ids in first-seen order. Node-based maps allocated once
    // per edge and dominated the index on a few thousand faces.
    struct FlatId {
      std::vector<std::uint64_t> keys;
      std::vector<int> ids;
      std::size_t mask = 0;
      int n = 0;
      void reset(int expect) {
        std::size_t cap = 16;
        const std::size_t need = static_cast<std::size_t>(std::max(expect, 1)) * 2;
        while (cap < need) {
          cap *= 2;
        }
        keys.assign(cap, ~std::uint64_t{0});
        ids.assign(cap, -1);
        mask = cap - 1;
        n = 0;
      }
      int intern(std::uint64_t key) {
        std::size_t h = static_cast<std::size_t>(key * 11400714819323198485ull) & mask;
        while (keys[h] != ~std::uint64_t{0} && keys[h] != key) {
          h = (h + 1) & mask;
        }
        if (keys[h] == ~std::uint64_t{0}) {
          keys[h] = key;
          const int id = n++;
          ids[h] = id;
          return id;
        }
        return ids[h];
      }
      int find(std::uint64_t key) const {
        std::size_t h = static_cast<std::size_t>(key * 11400714819323198485ull) & mask;
        while (keys[h] != ~std::uint64_t{0} && keys[h] != key) {
          h = (h + 1) & mask;
        }
        return keys[h] == key ? ids[h] : -1;
      }
    };
    FlatId edgeId;
    FlatId wedgeId;
    edgeId.reset(std::max(slots, 1));
    wedgeId.reset(std::max(slots, 1));
    const bool serialIndex =
        indexTeam <= 1 || static_cast<int>(cand.size()) < 256 || maxV < 0;
    std::vector<std::array<int, 2>> edgeAtoms;
    std::vector<int> vertexOf;
    int nV = 0;
    if (!serialIndex) {
      // Each thread records the edges and vertices first seen in its
      // slice of the ring list. Merging those slices in ring order
      // assigns the same ids as one thread walking the rings.
      struct FreshEdge {
        std::uint64_t key;
        int a;
        int b;
      };
      struct Bucket {
        std::vector<FreshEdge> edges;
        std::vector<int> verts;
        std::vector<std::uint64_t> wedges;
      };
      const int nCand = static_cast<int>(cand.size());
      std::vector<Bucket> buckets(static_cast<std::size_t>(indexTeam));
      std::vector<std::vector<int>> edgeMap(static_cast<std::size_t>(indexTeam));
      std::vector<std::vector<int>> wedgeMap(static_cast<std::size_t>(indexTeam));
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel num_threads(indexTeam)
#endif
      {
#if defined(SEAMS_HAS_OPENMP)
        const int tid = omp_get_thread_num();
        const int team = omp_get_num_threads();
#else
        const int tid = 0;
        const int team = 1;
#endif
        const int c0 = static_cast<int>((static_cast<long long>(nCand) * tid) / team);
        const int c1 = static_cast<int>((static_cast<long long>(nCand) * (tid + 1)) / team);
        FlatId local;
        local.reset(std::max(c1 - c0, 1) * 8);
        std::vector<int> stamp(static_cast<std::size_t>(maxV) + 1, 0);
        auto &bucket = buckets[static_cast<std::size_t>(tid)];
        bucket.edges.reserve(static_cast<std::size_t>(std::max(c1 - c0, 0)));
        bucket.verts.reserve(static_cast<std::size_t>(std::max(c1 - c0, 0)));
        for (int ci = c0; ci < c1; ++ci) {
          const int i = cand[static_cast<std::size_t>(ci)];
          const auto &ring = rings[static_cast<std::size_t>(i)];
          const int sz = ringSize[static_cast<std::size_t>(i)];
          int *eCol = ringEdges.col.data() + ringEdges.off[static_cast<std::size_t>(i)];
          for (int k = 0; k < sz; ++k) {
            const int a = ring[static_cast<std::size_t>(k)];
            const int b = ring[static_cast<std::size_t>((k + 1) % sz)];
            const std::uint64_t key = packEdge(a, b);
            const int before = local.n;
            const int id = local.intern(key);
            eCol[k] = id;
            if (id == before) {
              bucket.edges.push_back({key, std::min(a, b), std::max(a, b)});
            }
            if (a >= 0 && a <= maxV && stamp[static_cast<std::size_t>(a)] == 0) {
              stamp[static_cast<std::size_t>(a)] = 1;
              bucket.verts.push_back(a);
            }
          }
        }
      }
      vertexOf.assign(static_cast<std::size_t>(maxV) + 1, -1);
      for (int t = 0; t < indexTeam; ++t) {
        const auto &bucket = buckets[static_cast<std::size_t>(t)];
        auto &map = edgeMap[static_cast<std::size_t>(t)];
        map.resize(bucket.edges.size());
        for (std::size_t localId = 0; localId < bucket.edges.size(); ++localId) {
          const auto &edge = bucket.edges[localId];
          const int before = edgeId.n;
          const int id = edgeId.intern(edge.key);
          map[localId] = id;
          if (id == before) {
            edgeAtoms.push_back({edge.a, edge.b});
          }
        }
        for (const int atom : bucket.verts) {
          if (atom >= 0 && atom <= maxV &&
              vertexOf[static_cast<std::size_t>(atom)] < 0) {
            vertexOf[static_cast<std::size_t>(atom)] = nV++;
          }
        }
      }
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel for schedule(static) num_threads(indexTeam)
#endif
      for (int ci = 0; ci < nCand; ++ci) {
        const int i = cand[static_cast<std::size_t>(ci)];
        const int sz = ringSize[static_cast<std::size_t>(i)];
        const auto &ring = rings[static_cast<std::size_t>(i)];
        const int tid = static_cast<int>((static_cast<long long>(ci) * indexTeam) / nCand);
        const auto &map = edgeMap[static_cast<std::size_t>(tid)];
        int *eCol = ringEdges.col.data() + ringEdges.off[static_cast<std::size_t>(i)];
        int *vCol = ringVerts.col.data() + ringVerts.off[static_cast<std::size_t>(i)];
        for (int k = 0; k < sz; ++k) {
          eCol[k] = map[static_cast<std::size_t>(eCol[k])];
          const int a = ring[static_cast<std::size_t>(k)];
          vCol[k] = (a >= 0 && a <= maxV) ? vertexOf[static_cast<std::size_t>(a)] : -1;
        }
      }
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel num_threads(indexTeam)
#endif
      {
#if defined(SEAMS_HAS_OPENMP)
        const int tid = omp_get_thread_num();
        const int team = omp_get_num_threads();
#else
        const int tid = 0;
        const int team = 1;
#endif
        const int c0 = static_cast<int>((static_cast<long long>(nCand) * tid) / team);
        const int c1 = static_cast<int>((static_cast<long long>(nCand) * (tid + 1)) / team);
        FlatId local;
        local.reset(std::max(c1 - c0, 1) * 8);
        auto &bucket = buckets[static_cast<std::size_t>(tid)];
        bucket.wedges.clear();
        bucket.wedges.reserve(static_cast<std::size_t>(std::max(c1 - c0, 0)) * 6);
        for (int ci = c0; ci < c1; ++ci) {
          const int i = cand[static_cast<std::size_t>(ci)];
          const int sz = ringSize[static_cast<std::size_t>(i)];
          const int *eCol =
              ringEdges.col.data() + ringEdges.off[static_cast<std::size_t>(i)];
          int *wCol = ringWedges.col.data() + ringWedges.off[static_cast<std::size_t>(i)];
          for (int k = 0; k < sz; ++k) {
            const auto e1 = static_cast<std::uint64_t>(eCol[k]);
            const auto e2 = static_cast<std::uint64_t>(eCol[(k + 1) % sz]);
            const std::uint64_t key = (std::min(e1, e2) << 32) | std::max(e1, e2);
            const int before = local.n;
            const int id = local.intern(key);
            wCol[k] = id;
            if (id == before) {
              bucket.wedges.push_back(key);
            }
          }
        }
      }
      for (int t = 0; t < indexTeam; ++t) {
        const auto &bucket = buckets[static_cast<std::size_t>(t)];
        auto &map = wedgeMap[static_cast<std::size_t>(t)];
        map.resize(bucket.wedges.size());
        for (std::size_t localId = 0; localId < bucket.wedges.size(); ++localId) {
          map[localId] = wedgeId.intern(bucket.wedges[localId]);
        }
      }
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel for schedule(static) num_threads(indexTeam)
#endif
      for (int ci = 0; ci < nCand; ++ci) {
        const int i = cand[static_cast<std::size_t>(ci)];
        const int sz = ringSize[static_cast<std::size_t>(i)];
        int *wCol = ringWedges.col.data() + ringWedges.off[static_cast<std::size_t>(i)];
        const int tid = static_cast<int>((static_cast<long long>(ci) * indexTeam) / nCand);
        const auto &map = wedgeMap[static_cast<std::size_t>(tid)];
        for (int k = 0; k < sz; ++k) {
          wCol[k] = map[static_cast<std::size_t>(wCol[k])];
        }
      }
    } else {
    if (maxV >= 0) {
      vertexOf.assign(static_cast<size_t>(maxV) + 1, -1);
    }
    for (const int i : cand) {
      const auto &ring = rings[static_cast<size_t>(i)];
      const int sz = ringSize[static_cast<size_t>(i)];
      int *eCol = ringEdges.col.data() + ringEdges.off[static_cast<size_t>(i)];
      int *vCol = ringVerts.col.data() + ringVerts.off[static_cast<size_t>(i)];
      for (int k = 0; k < sz; ++k) {
        const int a = ring[static_cast<size_t>(k)];
        const int b = ring[static_cast<size_t>((k + 1) % sz)];
        const int before = edgeId.n;
        const int id = edgeId.intern(packEdge(a, b));
        if (id == before) {
          edgeAtoms.push_back({std::min(a, b), std::max(a, b)});
        }
        eCol[k] = id;
        int vid = -1;
        if (a >= 0 && a < static_cast<int>(vertexOf.size())) {
          vid = vertexOf[static_cast<size_t>(a)];
          if (vid < 0) {
            vid = nV++;
            vertexOf[static_cast<size_t>(a)] = vid;
          }
        }
        vCol[k] = vid;
      }
      int *wCol = ringWedges.col.data() + ringWedges.off[static_cast<size_t>(i)];
      for (int k = 0; k < sz; ++k) {
        const auto e1 = static_cast<std::uint64_t>(eCol[k]);
        const auto e2 = static_cast<std::uint64_t>(eCol[(k + 1) % sz]);
        const std::uint64_t key = (std::min(e1, e2) << 32) | std::max(e1, e2);
        wCol[k] = wedgeId.intern(key);
      }
    }
    }
    edgeEnds.resize(edgeAtoms.size());
    for (std::size_t e = 0; e < edgeAtoms.size(); ++e) {
      const int a = edgeAtoms[e][0];
      const int b = edgeAtoms[e][1];
      const int va = (a >= 0 && a < static_cast<int>(vertexOf.size()))
                         ? vertexOf[static_cast<size_t>(a)]
                         : -1;
      const int vb = (b >= 0 && b < static_cast<int>(vertexOf.size()))
                         ? vertexOf[static_cast<size_t>(b)]
                         : -1;
      edgeEnds[e] = {va, vb};
    }
    ringsOfEdge.off.assign(edgeEnds.size() + 1, 0);
    const int nEdge = static_cast<int>(edgeEnds.size());
    const bool parallelRows = !serialIndex && nEdge > 0 && indexTeam > 1;
    if (!parallelRows) {
    for (const int i : cand) {
      const auto row = ringEdges[static_cast<size_t>(i)];
      for (const int e : row) {
        ringsOfEdge.off[static_cast<size_t>(e) + 1]++;
      }
    }
    for (std::size_t e = 0; e < edgeEnds.size(); ++e) {
      ringsOfEdge.off[e + 1] += ringsOfEdge.off[e];
    }
    ringsOfEdge.col.assign(ringsOfEdge.off.back(), 0);
    std::vector<int> at(ringsOfEdge.off.begin(), ringsOfEdge.off.end() - 1);
    for (const int i : cand) {
      for (const int e : ringEdges[static_cast<size_t>(i)]) {
        ringsOfEdge.col[static_cast<size_t>(at[static_cast<size_t>(e)]++)] = i;
      }
    }
    } else {
      const int nCandRows = static_cast<int>(cand.size());
      std::vector<std::vector<int>> part(
          static_cast<std::size_t>(indexTeam),
          std::vector<int>(static_cast<std::size_t>(nEdge), 0));
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel num_threads(indexTeam)
#endif
      {
#if defined(SEAMS_HAS_OPENMP)
        const int tid = omp_get_thread_num();
        const int team = omp_get_num_threads();
#else
        const int tid = 0;
        const int team = 1;
#endif
        const int c0 = static_cast<int>((static_cast<long long>(nCandRows) * tid) / team);
        const int c1 =
            static_cast<int>((static_cast<long long>(nCandRows) * (tid + 1)) / team);
        auto &mine = part[static_cast<std::size_t>(tid)];
        for (int ci = c0; ci < c1; ++ci) {
          const int i = cand[static_cast<std::size_t>(ci)];
          for (const int e : ringEdges[static_cast<std::size_t>(i)]) {
            if (e >= 0 && e < nEdge) {
              ++mine[static_cast<std::size_t>(e)];
            }
          }
        }
      }
      for (int e = 0; e < nEdge; ++e) {
        int sum = 0;
        for (int t = 0; t < indexTeam; ++t) {
          sum += part[static_cast<std::size_t>(t)][static_cast<std::size_t>(e)];
        }
        ringsOfEdge.off[static_cast<std::size_t>(e) + 1] =
            ringsOfEdge.off[static_cast<std::size_t>(e)] + sum;
      }
      ringsOfEdge.col.assign(ringsOfEdge.off.back(), 0);
      std::vector<std::vector<int>> start(
          static_cast<std::size_t>(indexTeam),
          std::vector<int>(static_cast<std::size_t>(nEdge), 0));
      for (int e = 0; e < nEdge; ++e) {
        int cursor = ringsOfEdge.off[static_cast<std::size_t>(e)];
        for (int t = 0; t < indexTeam; ++t) {
          start[static_cast<std::size_t>(t)][static_cast<std::size_t>(e)] = cursor;
          cursor += part[static_cast<std::size_t>(t)][static_cast<std::size_t>(e)];
        }
      }
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel num_threads(indexTeam)
#endif
      {
#if defined(SEAMS_HAS_OPENMP)
        const int tid = omp_get_thread_num();
        const int team = omp_get_num_threads();
#else
        const int tid = 0;
        const int team = 1;
#endif
        const int c0 = static_cast<int>((static_cast<long long>(nCandRows) * tid) / team);
        const int c1 =
            static_cast<int>((static_cast<long long>(nCandRows) * (tid + 1)) / team);
        auto &cursor = start[static_cast<std::size_t>(tid)];
        for (int ci = c0; ci < c1; ++ci) {
          const int i = cand[static_cast<std::size_t>(ci)];
          for (const int e : ringEdges[static_cast<std::size_t>(i)]) {
            ringsOfEdge.col[static_cast<std::size_t>(cursor[static_cast<std::size_t>(e)]++)] = i;
          }
        }
      }
    }
    edgeUse.assign(edgeEnds.size(), 0);
    halfPos.assign(edgeEnds.size(), -1);
    wedgeUse.assign(static_cast<size_t>(wedgeId.n), 0);
    openAt.assign(static_cast<size_t>(nV), 0);
    facesAt.assign(static_cast<size_t>(nV), 0);
    optsAt.resize(static_cast<size_t>(maxFaces) + 1);
    found.resize(nSig);
  }

  // Workers share the spine and keep a private stack. The stack is empty
  // here, so the copy is the occupancy vectors, not the ring index.
  ClosedSearch(const ClosedSearch &o)
      : spine(o.spine), rings(o.rings), sigs(o.sigs), need(spine->need),
        facesOf(spine->facesOf), maxVerts(spine->maxVerts),
        maxEdges(spine->maxEdges), cand(spine->cand), ringSize(spine->ringSize),
        ringEdges(spine->ringEdges), ringWedges(spine->ringWedges),
        ringVerts(spine->ringVerts), ringsOfEdge(spine->ringsOfEdge),
        edgeEnds(spine->edgeEnds), maxFaces(o.maxFaces), have(o.have),
        room(o.room), fits(o.fits), edgeUse(o.edgeUse), wedgeUse(o.wedgeUse),
        half(o.half), halfPos(o.halfPos), openAt(o.openAt), facesAt(o.facesAt),
        used(o.used), chosen(o.chosen), optsAt(o.optsAt), found(o.found),
        branchOverride(o.branchOverride), seenFaces(o.seenFaces),
        trivalent(o.trivalent), seed(o.seed), nVerts(o.nVerts),
        nEdges(o.nEdges) {}

  void updateRoom(int s) {
    unsigned m = 0;
    for (size_t k = 0; k < sigs.size(); ++k) {
      if (have[static_cast<size_t>(s)] < need[k][static_cast<size_t>(s)]) {
        m |= 1u << k;
      }
    }
    room[static_cast<size_t>(s)] = m;
  }

  void halfAdd(int e) {
    halfPos[static_cast<size_t>(e)] = static_cast<int>(half.size());
    half.push_back(e);
  }

  void halfDrop(int e) {
    const int at = halfPos[static_cast<size_t>(e)];
    const int last = half.back();
    half[static_cast<size_t>(at)] = last;
    halfPos[static_cast<size_t>(last)] = at;
    half.pop_back();
    halfPos[static_cast<size_t>(e)] = -1;
  }

  bool canAdd(int ring) const {
    if (ring < seed || used[static_cast<size_t>(ring)]) {
      return false;
    }
    const int sz = ringSize[static_cast<size_t>(ring)];
    const unsigned m = fits[chosen.size()] & room[static_cast<size_t>(sz)];
    if (m == 0) {
      return false;
    }
    for (const int e : ringEdges[static_cast<size_t>(ring)]) {
      if (edgeUse[static_cast<size_t>(e)] >= 2) {
        return false;
      }
    }
    for (const int w : ringWedges[static_cast<size_t>(ring)]) {
      if (wedgeUse[static_cast<size_t>(w)] != 0) {
        return false;
      }
    }
    int newVerts = 0;
    for (const int v : ringVerts[static_cast<size_t>(ring)]) {
      if (facesAt[static_cast<size_t>(v)] == 0) {
        ++newVerts;
      } else if (openAt[static_cast<size_t>(v)] == 0) {
        return false;
      }
    }
    int newEdges = 0;
    for (const int e : ringEdges[static_cast<size_t>(ring)]) {
      newEdges += edgeUse[static_cast<size_t>(e)] == 0 ? 1 : 0;
    }
    for (unsigned b = m; b != 0; b &= b - 1) {
      const auto k = static_cast<size_t>(std::countr_zero(b));
      if (nVerts + newVerts <= maxVerts[k] && nEdges + newEdges <= maxEdges[k]) {
        return true;
      }
    }
    return false;
  }

  void addRing(int ring) {
    const int sz = ringSize[static_cast<size_t>(ring)];
    const unsigned m = fits[chosen.size()] & room[static_cast<size_t>(sz)];
    used[static_cast<size_t>(ring)] = 1;
    chosen.push_back(ring);
    fits[chosen.size()] = m;
    have[static_cast<size_t>(sz)] += 1;
    updateRoom(sz);
    for (const int e : ringEdges[static_cast<size_t>(ring)]) {
      const auto &ends = edgeEnds[static_cast<size_t>(e)];
      auto &u = edgeUse[static_cast<size_t>(e)];
      if (u == 0) {
        ++nEdges;
        halfAdd(e);
        ++openAt[static_cast<size_t>(ends[0])];
        ++openAt[static_cast<size_t>(ends[1])];
      } else {
        halfDrop(e);
        --openAt[static_cast<size_t>(ends[0])];
        --openAt[static_cast<size_t>(ends[1])];
      }
      ++u;
    }
    for (const int w : ringWedges[static_cast<size_t>(ring)]) {
      ++wedgeUse[static_cast<size_t>(w)];
    }
    for (const int v : ringVerts[static_cast<size_t>(ring)]) {
      nVerts += facesAt[static_cast<size_t>(v)] == 0 ? 1 : 0;
      ++facesAt[static_cast<size_t>(v)];
    }
  }

  void undoRing() {
    const int ring = chosen.back();
    const int sz = ringSize[static_cast<size_t>(ring)];
    chosen.pop_back();
    used[static_cast<size_t>(ring)] = 0;
    have[static_cast<size_t>(sz)] -= 1;
    updateRoom(sz);
    for (const int e : ringEdges[static_cast<size_t>(ring)]) {
      const auto &ends = edgeEnds[static_cast<size_t>(e)];
      auto &u = edgeUse[static_cast<size_t>(e)];
      --u;
      if (u == 0) {
        --nEdges;
        halfDrop(e);
        --openAt[static_cast<size_t>(ends[0])];
        --openAt[static_cast<size_t>(ends[1])];
      } else {
        halfAdd(e);
        ++openAt[static_cast<size_t>(ends[0])];
        ++openAt[static_cast<size_t>(ends[1])];
      }
    }
    for (const int w : ringWedges[static_cast<size_t>(ring)]) {
      --wedgeUse[static_cast<size_t>(w)];
    }
    for (const int v : ringVerts[static_cast<size_t>(ring)]) {
      --facesAt[static_cast<size_t>(v)];
      nVerts -= facesAt[static_cast<size_t>(v)] == 0 ? 1 : 0;
    }
  }

  void restore(size_t mark) {
    while (chosen.size() > mark) {
      undoRing();
    }
  }

  // These three clathrate signatures are trivalent, and so are the hydrate
  // cages the tests ask for when every requested signature has 2E = 3V.
  // A vertex then never takes a fourth face, three faces leave it closed,
  // and two faces leave exactly one open wedge.
  [[nodiscard]] bool shortVertex(int ring) const {
    if (!trivalent) {
      return false;
    }
    for (const int v : ringVerts[static_cast<size_t>(ring)]) {
      const int faces = facesAt[static_cast<size_t>(v)];
      const int open = openAt[static_cast<size_t>(v)];
      if (faces > 3) {
        return true;
      }
      if (faces == 3 && open != 0) {
        return true;
      }
      if (faces >= 2 && open > 2) {
        return true;
      }
      if (faces == 2 && open != 2) {
        return true;
      }
      if (faces > 0 && faces < 3 && open == 0) {
        return true;
      }
    }
    return false;
  }

  // Each open edge still needs one face-edge. Faces not yet placed have only
  // so many edges between them.
  [[nodiscard]] bool boundaryCloseable() const {
    const int facesNow = static_cast<int>(chosen.size());
    if (facesNow > maxFaces) {
      return false;
    }
    const int openEdges = static_cast<int>(half.size());
    if (openEdges == 0) {
      return true;
    }
    const unsigned mask = fits[static_cast<size_t>(facesNow)];
    for (size_t k = 0; k < sigs.size(); ++k) {
      if (((mask >> k) & 1u) == 0u) {
        continue;
      }
      int cap = 0;
      for (size_t s = 0; s < need[k].size(); ++s) {
        const int left = need[k][s] - have[static_cast<size_t>(s)];
        if (left > 0) {
          cap += left * static_cast<int>(s);
        }
      }
      if (openEdges <= cap) {
        return true;
      }
    }
    return false;
  }

  // A vertex with two faces in a trivalent cage has two open edges, and the
  // third face is a ring through both. That ring is tried before a free
  // choice on a single edge.
  bool forceWedge(int &edgeOut) {
    branchOverride.clear();
    if (!trivalent) {
      return false;
    }
    for (const int e : half) {
      const auto &ends = edgeEnds[static_cast<size_t>(e)];
      for (int t = 0; t < 2; ++t) {
        const int v = ends[static_cast<size_t>(t)];
        if (facesAt[static_cast<size_t>(v)] != 2 ||
            openAt[static_cast<size_t>(v)] != 2) {
          continue;
        }
        int e1 = -1;
        int e2 = -1;
        for (const int h : half) {
          const auto &he = edgeEnds[static_cast<size_t>(h)];
          if (he[0] != v && he[1] != v) {
            continue;
          }
          if (e1 < 0) {
            e1 = h;
          } else {
            e2 = h;
          }
        }
        if (e1 < 0 || e2 < 0) {
          return false;
        }
        const auto &a = ringsOfEdge[static_cast<size_t>(e1)];
        const auto &b = ringsOfEdge[static_cast<size_t>(e2)];
        auto ia = std::lower_bound(a.begin(), a.end(), seed);
        auto ib = std::lower_bound(b.begin(), b.end(), seed);
        while (ia != a.end() && ib != b.end()) {
          if (*ia == *ib) {
            if (canAdd(*ia)) {
              branchOverride.push_back(*ia);
            }
            ++ia;
            ++ib;
          } else if (*ia < *ib) {
            ++ia;
          } else {
            ++ib;
          }
        }
        edgeOut = e1;
        return true;
      }
    }
    return false;
  }

  int narrow() {
    for (;;) {
      int wedgeEdge = -1;
      if (forceWedge(wedgeEdge)) {
        if (branchOverride.empty()) {
          return -2;
        }
        if (branchOverride.size() == 1) {
          const int only = branchOverride[0];
          branchOverride.clear();
          addRing(only);
          if (shortVertex(only) || !boundaryCloseable()) {
            return -2;
          }
          continue;
        }
        return wedgeEdge;
      }
      int bestEdge = -1;
      size_t best = 0;
      bool forced = false;
      for (const int e : half) {
        size_t n = 0;
        int only = -1;
        const auto &row = ringsOfEdge[static_cast<size_t>(e)];
        for (auto it = std::lower_bound(row.begin(), row.end(), seed);
             it != row.end(); ++it) {
          if (!canAdd(*it)) {
            continue;
          }
          only = *it;
          if (++n > 1 && bestEdge >= 0 && n >= best) {
            break;
          }
        }
        if (n == 0) {
          return -2;
        }
        if (n == 1) {
          addRing(only);
          if (shortVertex(only) || !boundaryCloseable()) {
            return -2;
          }
          forced = true;
          break;
        }
        if (bestEdge < 0 || n < best) {
          best = n;
          bestEdge = e;
        }
      }
      if (!forced) {
        return bestEdge;
      }
    }
  }

  void accept() {
    const unsigned m = fits[chosen.size()];
    for (size_t k = 0; k < sigs.size(); ++k) {
      if (((m >> k) & 1u) == 0 ||
          static_cast<int>(chosen.size()) != facesOf[k]) {
        continue;
      }
      if (!facesTileSphere(rings, chosen)) {
        return;
      }
      cage::FoundCage cage;
      cage.signature = sigs[k];
      cage.faces = chosen;
      std::sort(cage.faces.begin(), cage.faces.end());
      cage.vertices = uniqueVertices(rings, cage.faces);
      found[k].push_back(std::move(cage));
      return;
    }
  }

  void search() {
    const size_t mark = chosen.size();
    const int edge = narrow();
    if (edge == -2) {
      restore(mark);
      return;
    }
    if (half.empty()) {
      accept();
      restore(mark);
      return;
    }
    if (edge < 0 || !boundaryCloseable()) {
      restore(mark);
      return;
    }
    std::vector<int> faceKey = chosen;
    std::sort(faceKey.begin(), faceKey.end());
    if (!seenFaces.insert(std::move(faceKey)).second) {
      restore(mark);
      return;
    }
    auto &opts = optsAt[chosen.size()];
    opts.clear();
    if (!branchOverride.empty()) {
      opts.swap(branchOverride);
    } else {
      const auto &row = ringsOfEdge[static_cast<size_t>(edge)];
      for (auto it = std::lower_bound(row.begin(), row.end(), seed);
           it != row.end(); ++it) {
        if (canAdd(*it)) {
          opts.push_back(*it);
        }
      }
    }
    for (size_t k = 0; k < opts.size(); ++k) {
      const int r = opts[k];
      if (!canAdd(r)) {
        continue;
      }
      addRing(r);
      if (shortVertex(r) || !boundaryCloseable()) {
        undoRing();
        continue;
      }
      search();
      undoRing();
    }
    restore(mark);
  }

  // A signature whose ring-size census exceeds the rings on hand cannot
  // match. Ice, which has no pentagons, then skips the clathrate search.
  bool feasible() const {
    if (sigs.empty()) {
      return false;
    }
    std::vector<int> got(room.size(), 0);
    for (const int r : cand) {
      const int s = ringSize[static_cast<size_t>(r)];
      if (s >= 0 && s < static_cast<int>(got.size())) {
        ++got[static_cast<size_t>(s)];
      }
    }
    for (size_t k = 0; k < sigs.size(); ++k) {
      bool ok = true;
      for (size_t s = 0; s < got.size(); ++s) {
        if (need[k][s] > got[s]) {
          ok = false;
          break;
        }
      }
      if (ok) {
        return true;
      }
    }
    return false;
  }

  void searchOne(int r) {
    seed = r;
    seenFaces.clear();
    bool seedOk = true;
    for (const int e : ringEdges[static_cast<size_t>(r)]) {
      if (ringsOfEdge[static_cast<size_t>(e)].size() < 2) {
        seedOk = false;
        break;
      }
    }
    if (!seedOk) {
      return;
    }
    addRing(r);
    search();
    undoRing();
  }

  void searchSeeds(const std::vector<int> &seeds) {
    for (const int r : seeds) {
      searchOne(r);
    }
  }

  void publish() {
    for (auto &list : found) {
      std::sort(list.begin(), list.end(),
                [](const cage::FoundCage &a, const cage::FoundCage &b) {
                  return a.vertices != b.vertices ? a.vertices < b.vertices
                                                  : a.faces < b.faces;
                });
      list.erase(std::unique(list.begin(), list.end(),
                             [](const cage::FoundCage &a,
                                const cage::FoundCage &b) {
                               return a.vertices == b.vertices;
                             }),
                 list.end());
      std::sort(list.begin(), list.end(),
                [](const cage::FoundCage &a, const cage::FoundCage &b) {
                  return a.faces < b.faces;
                });
      std::vector<std::vector<int>> faceRings;
      for (auto &c : list) {
        faceRings.clear();
        for (const int fi : c.faces) {
          faceRings.push_back(rings[static_cast<size_t>(fi)]);
        }
        c.certificate = cage::canonicalCertificate(faceRings);
      }
    }
  }

  // Every face set is reached once, from its lowest ring. Face sets that
  // share their vertices keep the one with the lowest faces, so the result
  // does not depend on the search order or on the other signatures.
  // Threads search disjoint seeds and the lists are merged before that
  // sort, so the cages do not depend on the thread count.
  void run() {
    if (!feasible()) {
      publish();
      return;
    }
    int nThreads = 1;
#if defined(SEAMS_HAS_OPENMP)
    if (!omp_in_parallel()) {
      nThreads = omp_get_max_threads();
    }
#endif
    if (nThreads <= 1 || cand.size() < 96) {
      searchSeeds(cand);
    } else {
      std::vector<std::vector<std::vector<cage::FoundCage>>> partial(
          static_cast<size_t>(nThreads));
#if defined(SEAMS_HAS_OPENMP)
#pragma omp parallel
#endif
      {
#if defined(SEAMS_HAS_OPENMP)
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        ClosedSearch local(*this);
        for (auto &list : local.found) {
          list.clear();
        }
        // Low ring indices own more completions, so a static stride leaves
        // the first threads holding the long seeds. Dynamic chunks hand
        // those seeds out as threads finish.
        const int nCand = static_cast<int>(cand.size());
#if defined(SEAMS_HAS_OPENMP)
#pragma omp for schedule(dynamic, 32)
#endif
        for (int i = 0; i < nCand; ++i) {
          local.searchOne(cand[static_cast<size_t>(i)]);
        }
        if (tid < nThreads) {
          partial[static_cast<size_t>(tid)] = std::move(local.found);
        }
      }
      for (int t = 0; t < nThreads; t++) {
        if (partial[static_cast<size_t>(t)].size() != found.size()) {
          continue;
        }
        for (size_t k = 0; k < found.size(); k++) {
          auto &src = partial[static_cast<size_t>(t)][k];
          found[k].insert(found[k].end(), src.begin(), src.end());
        }
      }
    }
    publish();
  }
};

} // namespace

namespace cage {

namespace {

bool parseBound(std::string_view text, double &out) {
  while (!text.empty() &&
         std::isspace(static_cast<unsigned char>(text.front()))) {
    text.remove_prefix(1);
  }
  while (!text.empty() &&
         std::isspace(static_cast<unsigned char>(text.back()))) {
    text.remove_suffix(1);
  }
  if (text.empty()) {
    return false;
  }
  double value = 0.0;
  if (!seams::parseExactDouble(text, value) || !std::isfinite(value)) {
    return false;
  }
  out = value;
  return true;
}

} // namespace

EdgeWindow parseEdgeWindow(std::string_view spec) {
  const std::string key = lowerCopy(trimmed(spec));
  if (key.empty()) {
    return {};
  }
  if (key == "none" || key == "off" || key == "false") {
    EdgeWindow window;
    window.enabled = false;
    return window;
  }
  const auto sep = key.find_first_of(",:/");
  if (sep == std::string::npos || sep == 0 || sep + 1 >= key.size() ||
      key.find_first_of(",:/", sep + 1) != std::string::npos) {
    throw std::invalid_argument(
        "bad cage edge window (want LO,HI or none)");
  }
  EdgeWindow window;
  if (!parseBound(key.substr(0, sep), window.lo) ||
      !parseBound(key.substr(sep + 1), window.hi) || !(window.lo >= 0.0) ||
      !(window.hi > window.lo)) {
    throw std::invalid_argument(
        "bad cage edge window (want 0 <= lo < hi, or none)");
  }
  window.enabled = true;
  return window;
}

std::string edgeWindowSpec(const EdgeWindow &window) {
  if (!window.enabled) {
    return "none";
  }
  std::ostringstream os;
  os << window.lo << "," << window.hi;
  return os.str();
}

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

bool speciesCycleMatches(const std::vector<int> &cycle,
                         const std::vector<int> &pattern) {
  const size_t n = cycle.size();
  if (n == 0 || n != pattern.size()) {
    return false;
  }
  for (size_t shift = 0; shift < n; ++shift) {
    bool forward = true;
    bool reverse = true;
    for (size_t i = 0; i < n; ++i) {
      if (cycle[(shift + i) % n] != pattern[i]) {
        forward = false;
      }
      if (cycle[(shift + n - i) % n] != pattern[i]) {
        reverse = false;
      }
    }
    if (forward || reverse) {
      return true;
    }
  }
  return false;
}

namespace {

std::vector<int> speciesCycle(const std::vector<int> &ring,
                              const std::vector<int> &species) {
  std::vector<int> cycle;
  cycle.reserve(ring.size());
  for (const int atom : ring) {
    if (atom < 0 || static_cast<size_t>(atom) >= species.size()) {
      return {};
    }
    cycle.push_back(species[static_cast<size_t>(atom)]);
  }
  return cycle;
}

std::vector<char> admitBySpecies(const std::vector<std::vector<int>> &rings,
                                 const std::vector<int> &species,
                                 const std::vector<std::vector<int>> &patterns) {
  std::vector<char> admit(rings.size(), 0);
  for (size_t i = 0; i < rings.size(); ++i) {
    const auto cycle = speciesCycle(rings[i], species);
    if (cycle.empty()) {
      continue;
    }
    for (const auto &pattern : patterns) {
      if (speciesCycleMatches(cycle, pattern)) {
        admit[i] = 1;
        break;
      }
    }
  }
  return admit;
}

} // namespace

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

bool isSpherePolyhedron(const std::vector<std::vector<int>> &rings,
                        const std::vector<int> &faces) {
  return isClosedPolyhedron(rings, faces) && facesTileSphere(rings, faces);
}

std::vector<FoundCage> findBySignature(const std::vector<std::vector<int>> &rings,
                                       const Signature &signature) {
  if (signature.counts.empty() || signature.faceCount() <= 0) {
    return {};
  }
  const std::vector<Signature> one = {signature};
  ClosedSearch search(rings, one);
  search.run();
  return std::move(search.found[0]);
}

std::vector<FoundCage>
findBySignature(const std::vector<std::vector<int>> &rings,
                const Signature &signature, const std::vector<int> &species,
                const std::vector<std::vector<int>> &patterns) {
  if (patterns.empty()) {
    return findBySignature(rings, signature);
  }
  if (signature.counts.empty() || signature.faceCount() <= 0) {
    return {};
  }
  const auto admit = admitBySpecies(rings, species, patterns);
  const std::vector<Signature> one = {signature};
  ClosedSearch search(rings, one, &admit);
  search.run();
  return std::move(search.found[0]);
}

namespace {

Eigen::MatrixXd asMatrix(const std::vector<std::array<double, 3>> &xyz) {
  Eigen::MatrixXd m(static_cast<Eigen::Index>(xyz.size()), 3);
  for (size_t i = 0; i < xyz.size(); ++i) {
    m(static_cast<Eigen::Index>(i), 0) = xyz[i][0];
    m(static_cast<Eigen::Index>(i), 1) = xyz[i][1];
    m(static_cast<Eigen::Index>(i), 2) = xyz[i][2];
  }
  return m;
}

} // namespace

std::vector<std::array<double, 3>>
coordsOfVertices(const std::vector<std::array<double, 3>> &all,
                 const std::vector<int> &vertices) {
  std::vector<std::array<double, 3>> out;
  out.reserve(vertices.size());
  for (const int v : vertices) {
    if (v < 0 || static_cast<size_t>(v) >= all.size()) {
      continue;
    }
    out.push_back(all[static_cast<size_t>(v)]);
  }
  return out;
}

CageShape shapeOfVertices(const std::vector<std::array<double, 3>> &cageXyz) {
  CageShape shape;
  shape.nVertices = static_cast<int>(cageXyz.size());
  if (cageXyz.empty() || !ira::available()) {
    return shape;
  }
  ira::PointGroup group;
  shape.status = ira::pointGroup(asMatrix(cageXyz), group);
  if (shape.status == 0) {
    shape.pointGroup = group.symbol;
  }
  return shape;
}

std::vector<FormerRow>
formerRows(const std::vector<std::vector<int>> &nList,
           const std::vector<std::vector<int>> &rings,
           const std::vector<int> &species, int formerSpecies) {
  const int n = static_cast<int>(nList.size());
  std::vector<std::map<int, int>> through(static_cast<size_t>(n));
  for (const auto &ring : rings) {
    const int size = static_cast<int>(ring.size());
    for (const int atom : ring) {
      if (atom < 0 || atom >= n) {
        continue;
      }
      through[static_cast<size_t>(atom)][size] += 1;
    }
  }
  std::vector<FormerRow> out;
  for (int i = 0; i < n; ++i) {
    const int sp = static_cast<size_t>(i) < species.size()
                       ? species[static_cast<size_t>(i)]
                       : 0;
    if (formerSpecies >= 0 && sp != formerSpecies) {
      continue;
    }
    FormerRow row;
    row.index = i;
    row.species = sp;
    for (const int nb : nList[static_cast<size_t>(i)]) {
      if (nb == i) {
        continue;
      }
      row.coord += 1;
      if (nb < 0 || static_cast<size_t>(nb) >= species.size()) {
        continue;
      }
      if (species[static_cast<size_t>(nb)] == sp) {
        row.homopolar += 1;
      }
    }
    row.rings = through[static_cast<size_t>(i)];
    out.push_back(row);
  }
  return out;
}

bool sameNetwork(const std::vector<FormerRow> &early,
                 const std::vector<FormerRow> &late) {
  if (early.size() != late.size()) {
    return false;
  }
  for (size_t i = 0; i < early.size(); ++i) {
    const auto &a = early[i];
    const auto &b = late[i];
    if (a.index != b.index || a.species != b.species || a.coord != b.coord ||
        a.homopolar != b.homopolar || a.rings != b.rings) {
      return false;
    }
  }
  return true;
}

CageShape overlayVertices(const std::vector<std::array<double, 3>> &ref,
                          const std::vector<std::array<double, 3>> &cageXyz) {
  CageShape shape;
  shape.nVertices = static_cast<int>(cageXyz.size());
  if (cageXyz.empty() || ref.size() != cageXyz.size() || !ira::available()) {
    return shape;
  }
  ira::Match hit;
  shape.status = ira::match(asMatrix(ref), asMatrix(cageXyz), hit);
  if (shape.status == 0) {
    shape.rmsd = hit.rmsd;
  }
  return shape;
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
          const Signature &signature, const std::vector<int> &atomRank) {
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
  auto index = ring::buildRingSearchIndex(six, static_cast<int>(nList.size()));
  index.atomRank = atomRank;
  const auto listHC = ring::findHC(six, ringType, nList, cageList, index);
  if (want == cageType::DoubleDiaC) {
    ring::findDDC(six, ringType, listHC, cageList, index);
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
  return findBySignature(rings, nList, signature, {});
}

std::vector<FoundCage>
findBySignature(const std::vector<std::vector<int>> &rings,
                const std::vector<std::vector<int>> &nList,
                const Signature &signature, const std::vector<int> &atomRank) {
  if (signature.kind == Signature::Kind::HexC) {
    return tumByType(rings, nList, cageType::HexC, signature, atomRank);
  }
  if (signature.kind == Signature::Kind::DoubleDiaC) {
    return tumByType(rings, nList, cageType::DoubleDiaC, signature, atomRank);
  }
  return findBySignature(rings, signature);
}

std::vector<std::vector<FoundCage>>
findBySignatures(const std::vector<std::vector<int>> &rings,
                 const std::vector<std::vector<int>> &nList,
                 const std::vector<Signature> &signatures,
                 const std::vector<int> &atomRank,
                 const std::vector<std::array<double, 3>> *xyz,
                 const std::vector<double> *box, const EdgeWindow &window) {
  std::vector<std::vector<FoundCage>> out(signatures.size());
  std::vector<Signature> census;
  std::vector<std::vector<size_t>> askedAt;
  for (size_t i = 0; i < signatures.size(); ++i) {
    const auto &s = signatures[i];
    if (s.kind != Signature::Kind::Census || s.counts.empty() ||
        s.faceCount() <= 0) {
      out[i] = findBySignature(rings, nList, s, atomRank);
      continue;
    }
    const auto at = std::find(census.begin(), census.end(), s);
    if (at == census.end()) {
      census.push_back(s);
      askedAt.push_back({i});
    } else {
      askedAt[static_cast<size_t>(at - census.begin())].push_back(i);
    }
  }
  constexpr size_t kPerSearch = 32;
  for (size_t from = 0; from < census.size(); from += kPerSearch) {
    const size_t to = std::min(census.size(), from + kPerSearch);
    const std::vector<Signature> chunk(census.begin() + static_cast<std::ptrdiff_t>(from),
                                       census.begin() + static_cast<std::ptrdiff_t>(to));
    int maxSize = 0;
    for (const auto &s : chunk) {
      maxSize = std::max(maxSize, s.maxRingSize());
    }
    std::vector<int> got(static_cast<size_t>(maxSize) + 1, 0);
    for (const auto &ring : rings) {
      const int sz = static_cast<int>(ring.size());
      if (sz >= 0 && sz <= maxSize) {
        ++got[static_cast<size_t>(sz)];
      }
    }
    bool any = false;
    for (const auto &s : chunk) {
      bool ok = true;
      for (const auto &kv : s.counts) {
        if (kv.first > maxSize ||
            got[static_cast<size_t>(kv.first)] < kv.second) {
          ok = false;
          break;
        }
      }
      any = any || ok;
    }
    if (!any) {
      continue;
    }
    ClosedSearch search(rings, chunk, nullptr, xyz, box, window);
    search.run();
    for (size_t j = 0; j < chunk.size(); ++j) {
      for (const size_t i : askedAt[from + j]) {
        out[i] = search.found[j];
      }
    }
  }
  return out;
}

} // namespace cage
