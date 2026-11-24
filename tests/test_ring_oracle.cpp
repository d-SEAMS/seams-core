// Copyright (c) d-SEAMS developers
// The primitive ring search against a brute-force enumerator.
//
// tests/ring_oracle.hpp lists every simple cycle and keeps the ones with no
// shortcut. Each case hands one graph to every entry that returns rings:
// ringNetwork over its row builders, the CSR overload, the ball table, the
// source mask, the size floor, the size and per-atom counts, the
// generate-then-filter route, RingUpdater, the chordless cycles, and the C
// API on coordinates. The ring lists must be equal element by element.

#include "ice_ih.hpp"
#include "ring_oracle.hpp"
#include "seams.h"

#include <catch2/catch_test_macros.hpp>

#include <chordless.hpp>
#include <franzblau.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include "hydrate_sites.inc"
#include "ring_oracle_cells.inc"

namespace {

using oracle::Adjacency;
using oracle::Rings;

// Rows as the library takes them: the vertex first, then its partners in
// any order. Repeats and self entries are allowed in a row and are not
// edges.
std::vector<std::vector<int>> rowsOf(const Adjacency &adj, std::mt19937 *rng = nullptr) {
  std::vector<std::vector<int>> rows(adj.size());
  for (std::size_t i = 0; i < adj.size(); ++i) {
    auto &row = rows[i];
    row.push_back(static_cast<int>(i));
    row.insert(row.end(), adj[i].begin(), adj[i].end());
    if (rng != nullptr) {
      std::uniform_int_distribution<int> coin(0, 3);
      if (!adj[i].empty() && coin(*rng) == 0) {
        row.push_back(adj[i][static_cast<std::size_t>(coin(*rng)) % adj[i].size()]);
      }
      if (coin(*rng) == 0) {
        row.push_back(static_cast<int>(i));
      }
      std::shuffle(row.begin() + 1, row.end(), *rng);
    }
  }
  return rows;
}

void csrOf(const Adjacency &adj, std::vector<int> &offsets, std::vector<int> &cols) {
  offsets.assign(1, 0);
  cols.clear();
  for (const auto &row : adj) {
    cols.insert(cols.end(), row.rbegin(), row.rend());
    offsets.push_back(static_cast<int>(cols.size()));
  }
}

std::string show(const Rings &rings) {
  std::string s;
  for (const auto &r : rings) {
    s += "(";
    for (std::size_t k = 0; k < r.size(); ++k) {
      s += (k ? " " : "") + std::to_string(r[k]);
    }
    s += ") ";
  }
  return s;
}

std::string showGraph(const Adjacency &adj) {
  std::string s;
  for (std::size_t i = 0; i < adj.size(); ++i) {
    for (int j : adj[i]) {
      if (static_cast<int>(i) < j) {
        s += std::to_string(i) + "-" + std::to_string(j) + " ";
      }
    }
  }
  return s;
}

#define REQUIRE_RINGS(got, want, adj, depth, what)                                         \
  do {                                                                                     \
    if ((got) != (want)) {                                                                 \
      INFO(what << " depth " << (depth) << "\ngraph " << showGraph(adj) << "\nwant "        \
                << show(want) << "\ngot  " << show(got));                                  \
      REQUIRE((got) == (want));                                                            \
    }                                                                                      \
  } while (0)

std::vector<int> histogramOf(const Rings &rings, int len) {
  std::vector<int> h(static_cast<std::size_t>(len), 0);
  for (const auto &r : rings) {
    if (static_cast<int>(r.size()) < len) {
      ++h[r.size()];
    }
  }
  return h;
}

// Every in-process entry on one graph against the oracle.
void checkAllEntries(const Adjacency &adj, int depth, std::mt19937 &rng) {
  const Rings want = oracle::primitiveRings(adj, depth);
  const auto rows = rowsOf(adj, &rng);
  const int n = static_cast<int>(adj.size());

  REQUIRE_RINGS(primitive::ringNetwork(rows, depth), want, adj, depth, "ringNetwork");
  REQUIRE_RINGS(primitive::ringNetwork(rows, depth, primitive::RingRows::Vectors), want, adj,
                depth, "RingRows::Vectors");
  REQUIRE_RINGS(primitive::ringNetwork(rows, depth, primitive::RingRows::Packed), want, adj,
                depth, "RingRows::Packed");
  REQUIRE_RINGS(primitive::ringNetworkBalls(rows, depth), want, adj, depth, "ringNetworkBalls");

  std::vector<int> offsets;
  std::vector<int> cols;
  csrOf(adj, offsets, cols);
  if (n > 0) {
    REQUIRE_RINGS(primitive::ringNetwork(offsets, cols, depth), want, adj, depth,
                  "CSR ringNetwork");
    std::vector<int> counts(16, -1);
    primitive::ringSizeCounts(offsets, cols, depth, counts.data(), 16);
    REQUIRE(counts == histogramOf(want, 16));
  }

  std::uniform_int_distribution<int> pick(3, std::max(depth, 3));
  const int floor = pick(rng);
  Rings tall;
  for (const auto &r : want) {
    if (static_cast<int>(r.size()) >= floor) {
      tall.push_back(r);
    }
  }
  REQUIRE_RINGS(primitive::ringNetworkMinSize(rows, depth, floor), tall, adj, depth,
                "ringNetworkMinSize " + std::to_string(floor));

  std::vector<char> mask(static_cast<std::size_t>(n), 0);
  std::bernoulli_distribution half(0.5);
  for (auto &m : mask) {
    m = half(rng) ? 1 : 0;
  }
  Rings masked;
  for (const auto &r : want) {
    if (mask[static_cast<std::size_t>(r.front())] != 0) {
      masked.push_back(r);
    }
  }
  if (n > 0) {
    REQUIRE_RINGS(primitive::ringNetwork(rows, depth, mask), masked, adj, depth,
                  "source mask");
    REQUIRE_RINGS(primitive::ringNetwork(offsets, cols, depth, mask), masked, adj, depth,
                  "CSR source mask");
  }

  const auto perAtom = primitive::atomRingSizeCounts(rows, depth);
  const int stride = depth + 1;
  std::vector<int> wantPerAtom(static_cast<std::size_t>(n) * static_cast<std::size_t>(stride),
                               0);
  for (const auto &r : want) {
    for (int v : r) {
      ++wantPerAtom[static_cast<std::size_t>(v) * static_cast<std::size_t>(stride) + r.size()];
    }
  }
  REQUIRE(perAtom == wantPerAtom);

  if (depth <= 8) {
    auto graph = primitive::countAllRingsFromIndex(rows, depth);
    primitive::removeNonSPrings(graph);
    REQUIRE_RINGS(oracle::canonicalSorted(graph.rings), want, adj, depth,
                  "countAllRingsFromIndex + removeNonSPrings");
  }

  for (int size = 3; size <= std::min(depth, chordless::kMaxSize); ++size) {
    if (n == 0) {
      break;
    }
    std::vector<int> flat;
    const std::size_t count = chordless::cycles(offsets, cols, size, &flat, nullptr);
    Rings got(count);
    for (std::size_t c = 0; c < count; ++c) {
      got[c].assign(flat.begin() + static_cast<std::ptrdiff_t>(c * size),
                    flat.begin() + static_cast<std::ptrdiff_t>((c + 1) * size));
    }
    REQUIRE_RINGS(got, oracle::chordlessCycles(adj, size), adj, size, "chordless");
  }
}

Adjacency randomGraph(std::mt19937 &rng, int n, double meanDegree) {
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
  if (n > 1) {
    std::bernoulli_distribution edge(std::min(1.0, meanDegree / (n - 1)));
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (edge(rng)) {
          rows[static_cast<std::size_t>(i)].push_back(j);
        }
      }
    }
  }
  return oracle::simpleGraph(rows, n);
}

Adjacency cycleGraph(int k) {
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(k));
  for (int i = 0; i < k; ++i) {
    rows[static_cast<std::size_t>(i)].push_back((i + 1) % k);
  }
  return oracle::simpleGraph(rows, k);
}

struct Frame {
  std::vector<double> xyz;
  std::array<double, 9> box9{};
};

Frame replicate(const double (*xyz)[3], int n, const double *box9, int na, int nb, int nc) {
  Frame f;
  const auto h = oracle::lattice(box9);
  for (int a = 0; a < na; ++a) {
    for (int b = 0; b < nb; ++b) {
      for (int c = 0; c < nc; ++c) {
        for (int i = 0; i < n; ++i) {
          for (int k = 0; k < 3; ++k) {
            f.xyz.push_back(xyz[i][k] + a * h[0][k] + b * h[1][k] + c * h[2][k]);
          }
        }
      }
    }
  }
  f.box9 = {0, 0, 0, box9[3] * na, box9[4] * nb, box9[5] * nc,
            box9[6] * nb, box9[7] * nc, box9[8] * nc};
  return f;
}

Frame cubic(const double (*xyz)[3], int n, double side) {
  const double box9[9] = {0, 0, 0, side, side, side, 0, 0, 0};
  return replicate(xyz, n, box9, 1, 1, 1);
}

Frame iceIh(int nx, int ny, int nz) {
  const IceIhCell cell = iceIhOrtho(nx, ny, nz);
  Frame f;
  f.xyz = cell.xyz;
  std::copy(std::begin(cell.box9), std::end(cell.box9), f.box9.begin());
  return f;
}

Frame iceIc(int side) {
  const double a = 4.0 * 2.75 / std::sqrt(3.0);
  const double basis[8][3] = {{0.0, 0.0, 0.0},    {0.5, 0.5, 0.0},    {0.5, 0.0, 0.5},
                              {0.0, 0.5, 0.5},    {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  double xyz[8][3];
  for (int b = 0; b < 8; ++b) {
    for (int k = 0; k < 3; ++k) {
      xyz[b][k] = basis[b][k] * a;
    }
  }
  const double box9[9] = {0, 0, 0, a, a, a, 0, 0, 0};
  return replicate(xyz, 8, box9, side, side, side);
}

// A periodic column of stacked k-gons along z: each vertex bonds to its two
// ring partners and to the vertices above and below it. The prism search
// takes the longest box edge as the tube axis, so z is the longest.
Frame prismColumn(int k, int layers, double bond) {
  Frame f;
  const double pi = std::acos(-1.0);
  const double radius = bond / (2.0 * std::sin(pi / k));
  const double side = std::min(16.0, 0.9 * layers * bond);
  for (int l = 0; l < layers; ++l) {
    for (int i = 0; i < k; ++i) {
      f.xyz.push_back(side / 2 + radius * std::cos(2 * pi * i / k));
      f.xyz.push_back(side / 2 + radius * std::sin(2 * pi * i / k));
      f.xyz.push_back((l + 0.5) * bond);
    }
  }
  f.box9 = {0, 0, 0, side, side, layers * bond, 0, 0, 0};
  return f;
}

// The regular dodecahedron with edge `bond`, alone in a large box: one
// 5^12 cage of twelve pentagons.
Frame dodecahedron(double bond) {
  const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
  const double scale = bond / (2.0 / phi);
  std::vector<std::array<double, 3>> v;
  for (int sx : {-1, 1}) {
    for (int sy : {-1, 1}) {
      for (int sz : {-1, 1}) {
        v.push_back({1.0 * sx, 1.0 * sy, 1.0 * sz});
      }
      v.push_back({0.0, sx / phi, sy * phi});
      v.push_back({sx / phi, sy * phi, 0.0});
      v.push_back({sx * phi, 0.0, sy / phi});
    }
  }
  Frame f;
  const double side = 30.0;
  for (const auto &p : v) {
    for (int k = 0; k < 3; ++k) {
      f.xyz.push_back(side / 2 + scale * p[static_cast<std::size_t>(k)]);
    }
  }
  f.box9 = {0, 0, 0, side, side, side, 0, 0, 0};
  return f;
}

seams_system *systemOf(const Frame &f) {
  seams_system *sys = seams_system_create();
  REQUIRE(sys != nullptr);
  REQUIRE(seams_system_set_atoms(sys, f.xyz.data(), f.xyz.size() / 3, f.box9.data(), 9,
                                 nullptr, nullptr) == SEAMS_OK);
  return sys;
}

Rings ringsOf(const seams_rings *rings) {
  const int n = seams_rings_count(rings);
  Rings out(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    std::size_t len = 0;
    REQUIRE(seams_rings_vertices(rings, i, nullptr, 0, &len) == SEAMS_OK);
    out[static_cast<std::size_t>(i)].resize(len);
    REQUIRE(seams_rings_vertices(rings, i, out[static_cast<std::size_t>(i)].data(), len,
                                 &len) == SEAMS_OK);
  }
  return out;
}

Adjacency neighboursOf(seams_system *sys, double cutoff) {
  seams_neighbours *nb = nullptr;
  REQUIRE(seams_neighbours_build(sys, cutoff, 0, SEAMS_GRAPH_CUTOFF, &nb) == SEAMS_OK);
  std::size_t nOff = 0;
  std::size_t nnz = 0;
  REQUIRE(seams_neighbours_copy(nb, nullptr, 0, &nOff, nullptr, 0, &nnz) == SEAMS_OK);
  std::vector<int> off(nOff);
  std::vector<int> cols(nnz);
  REQUIRE(seams_neighbours_copy(nb, off.data(), nOff, &nOff, cols.data(), nnz, &nnz) ==
          SEAMS_OK);
  std::vector<std::vector<int>> rows(nOff > 0 ? nOff - 1 : 0);
  for (std::size_t i = 0; i + 1 < nOff; ++i) {
    rows[i].assign(cols.begin() + off[i], cols.begin() + off[i + 1]);
  }
  seams_neighbours_destroy(nb);
  return oracle::simpleGraph(rows, static_cast<int>(rows.size()));
}

// seams_rings_find and seams_rings_chordless on a frame against the oracle
// on the cutoff graph of the same coordinates.
void checkFrame(const Frame &f, double cutoff, int depth, const char *what) {
  double gap = 0;
  const Adjacency adj = oracle::cutoffGraph(f.xyz, f.box9.data(), cutoff, &gap);
  INFO(what << " cutoff " << cutoff << " depth " << depth);
  REQUIRE(gap > 1e-3);
  seams_system *sys = systemOf(f);
  REQUIRE(neighboursOf(sys, cutoff) == adj);
  seams_rings *rings = nullptr;
  REQUIRE(seams_rings_find(sys, cutoff, 0, SEAMS_GRAPH_CUTOFF, depth, &rings) == SEAMS_OK);
  const Rings want = oracle::primitiveRings(adj, depth);
  REQUIRE_RINGS(ringsOf(rings), want, adj, depth, std::string("seams_rings_find ") + what);
  std::vector<int> hist(16, -1);
  REQUIRE(seams_rings_histogram(rings, hist.data(), 16) == SEAMS_OK);
  REQUIRE(hist == histogramOf(want, 16));
  seams_rings_destroy(rings);

  seams_neighbours *nb = nullptr;
  REQUIRE(seams_neighbours_build(sys, cutoff, 0, SEAMS_GRAPH_CUTOFF, &nb) == SEAMS_OK);
  for (int size = 3; size <= depth; ++size) {
    const Rings cycles = oracle::chordlessCycles(adj, size);
    std::size_t count = 0;
    REQUIRE(seams_rings_chordless(nb, size, nullptr, 0, &count, nullptr, 0) == SEAMS_OK);
    REQUIRE(count == cycles.size());
    std::vector<int> flat(count * static_cast<std::size_t>(size));
    if (count > 0) {
      REQUIRE(seams_rings_chordless(nb, size, flat.data(), flat.size(), &count, nullptr, 0) ==
              SEAMS_OK);
    }
    Rings got(count);
    for (std::size_t c = 0; c < count; ++c) {
      got[c].assign(flat.begin() + static_cast<std::ptrdiff_t>(c * size),
                    flat.begin() + static_cast<std::ptrdiff_t>((c + 1) * size));
    }
    REQUIRE_RINGS(got, cycles, adj, size, std::string("seams_rings_chordless ") + what);
  }
  seams_neighbours_destroy(nb);
  seams_system_destroy(sys);
}

} // namespace

TEST_CASE("primitive rings of random graphs match the brute-force oracle", "[ring_oracle]") {
  for (int threads : {1, 4}) {
    REQUIRE(seams_set_num_threads(threads) == SEAMS_OK);
    std::mt19937 rng(20261121u + static_cast<unsigned>(threads));
    for (int trial = 0; trial < 1500; ++trial) {
      std::uniform_int_distribution<int> depthOf(3, 12);
      const int depth = depthOf(rng);
      const int nMax = depth <= 8 ? 40 : 26;
      std::uniform_int_distribution<int> nOf(0, nMax);
      std::uniform_real_distribution<double> degreeOf(0.5, depth <= 8 ? 4.5 : 3.4);
      const Adjacency adj = randomGraph(rng, nOf(rng), degreeOf(rng));
      INFO("threads " << threads << " trial " << trial);
      checkAllEntries(adj, depth, rng);
    }
  }
  REQUIRE(seams_set_num_threads(0) == SEAMS_OK);
}

TEST_CASE("primitive rings at the edges of the input", "[ring_oracle]") {
  std::mt19937 rng(7u);
  SECTION("no vertices, no edges, and lone atoms") {
    for (int depth : {3, 6, 12}) {
      checkAllEntries(Adjacency{}, depth, rng);
      checkAllEntries(Adjacency(1), depth, rng);
      checkAllEntries(Adjacency(9), depth, rng);
    }
  }
  SECTION("a k-cycle needs a depth of at least k") {
    for (int k = 3; k <= 13; ++k) {
      const Adjacency ring = cycleGraph(k);
      for (int depth = 3; depth <= 12; ++depth) {
        INFO("cycle " << k);
        const Rings want = oracle::primitiveRings(ring, depth);
        REQUIRE(want.size() == (depth >= k ? 1u : 0u));
        checkAllEntries(ring, depth, rng);
      }
    }
  }
  SECTION("disconnected components with isolated atoms between them") {
    std::vector<std::vector<int>> rows(20);
    rows[1] = {2, 4};
    rows[2] = {4};
    rows[6] = {7};
    rows[7] = {8};
    rows[8] = {9};
    rows[9] = {6};
    rows[12] = {13, 15, 16, 17};
    rows[13] = {14};
    rows[14] = {15};
    rows[16] = {17};
    rows[18] = {19};
    const Adjacency adj = oracle::simpleGraph(rows, 20);
    for (int depth = 3; depth <= 12; ++depth) {
      checkAllEntries(adj, depth, rng);
    }
  }
  SECTION("hub vertices past the 48-entry and 256-entry row limits") {
    for (int hubDegree : {49, 120, 300}) {
      std::mt19937 local(static_cast<unsigned>(hubDegree));
      Adjacency base = randomGraph(local, 400, 2.2);
      std::vector<std::vector<int>> rows(base.begin(), base.end());
      rows.emplace_back();
      std::vector<int> order(400);
      for (int i = 0; i < 400; ++i) {
        order[static_cast<std::size_t>(i)] = i;
      }
      std::shuffle(order.begin(), order.end(), local);
      for (int k = 0; k < hubDegree; ++k) {
        rows[400].push_back(order[static_cast<std::size_t>(k)]);
      }
      const Adjacency adj = oracle::simpleGraph(rows, 401);
      for (int depth : {5, 6, 8}) {
        INFO("hub degree " << hubDegree);
        checkAllEntries(adj, depth, rng);
      }
    }
  }
  SECTION("a sparse random graph of more than 4096 vertices") {
    std::mt19937 local(4097u);
    const Adjacency adj = randomGraph(local, 5000, 3.0);
    for (int depth : {7, 8, 10}) {
      checkAllEntries(adj, depth, rng);
    }
  }
  SECTION("theta graphs, where two rings tie with a third") {
    for (int a = 1; a <= 6; ++a) {
      for (int b = a; b <= 6; ++b) {
        for (int c = b; c <= 6; ++c) {
          // Three paths of a, b and c edges between vertices 0 and 1.
          std::vector<std::vector<int>> rows(2);
          for (int len : {a, b, c}) {
            int prev = 0;
            for (int s = 1; s < len; ++s) {
              rows.emplace_back();
              const int v = static_cast<int>(rows.size()) - 1;
              rows[static_cast<std::size_t>(prev)].push_back(v);
              prev = v;
            }
            rows[static_cast<std::size_t>(prev)].push_back(1);
          }
          const Adjacency adj = oracle::simpleGraph(rows, static_cast<int>(rows.size()));
          for (int depth : {6, 9, 12}) {
            INFO("theta " << a << " " << b << " " << c);
            checkAllEntries(adj, depth, rng);
          }
        }
      }
    }
  }
}

TEST_CASE("RingUpdater follows the oracle across changing frames", "[ring_oracle]") {
  std::mt19937 rng(99u);
  for (int run = 0; run < 12; ++run) {
    std::uniform_int_distribution<int> depthOf(4, 10);
    const int depth = depthOf(rng);
    const int n = 36;
    std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
    Adjacency adj = randomGraph(rng, n, 3.0);
    primitive::RingUpdater updater(depth);
    for (int frame = 0; frame < 8; ++frame) {
      INFO("run " << run << " frame " << frame << " depth " << depth);
      const Rings want = oracle::primitiveRings(adj, depth);
      REQUIRE_RINGS(updater.update(rowsOf(adj)), want, adj, depth, "RingUpdater");
      // Flip a few random pairs.
      std::uniform_int_distribution<int> vertex(0, n - 1);
      std::vector<std::vector<int>> next(adj.begin(), adj.end());
      for (int flip = 0; flip < 3; ++flip) {
        const int i = vertex(rng);
        const int j = vertex(rng);
        if (i == j) {
          continue;
        }
        auto &ri = next[static_cast<std::size_t>(i)];
        auto &rj = next[static_cast<std::size_t>(j)];
        const auto it = std::find(ri.begin(), ri.end(), j);
        if (it != ri.end()) {
          ri.erase(it);
          rj.erase(std::find(rj.begin(), rj.end(), i));
        } else {
          ri.push_back(j);
        }
      }
      adj = oracle::simpleGraph(next, n);
    }
  }
}

TEST_CASE("RingUpdater follows the oracle on changing ice frames", "[ring_oracle]") {
  std::mt19937 rng(1234u);
  for (const Frame &f : {iceIh(4, 2, 2), iceIh(9, 8, 8)}) {
    const Adjacency lattice = oracle::cutoffGraph(f.xyz, f.box9.data(), 3.2);
    const int n = static_cast<int>(lattice.size());
    for (int depth : {6, 8}) {
      primitive::RingUpdater updater(depth);
      std::vector<std::vector<int>> rows(lattice.begin(), lattice.end());
      for (int frame = 0; frame < 5; ++frame) {
        const Adjacency adj = oracle::simpleGraph(rows, n);
        INFO("molecules " << n << " depth " << depth << " frame " << frame);
        REQUIRE_RINGS(updater.update(rowsOf(adj)), oracle::primitiveRings(adj, depth), adj,
                      depth, "RingUpdater on ice");
        // Break a few bonds and add a few second-shell ones.
        std::uniform_int_distribution<int> vertex(0, n - 1);
        rows.assign(adj.begin(), adj.end());
        for (int change = 0; change < 6; ++change) {
          const int i = vertex(rng);
          auto &ri = rows[static_cast<std::size_t>(i)];
          if (change % 2 == 0 && !ri.empty()) {
            const int j = ri.front();
            ri.erase(ri.begin());
            auto &rj = rows[static_cast<std::size_t>(j)];
            rj.erase(std::find(rj.begin(), rj.end(), i));
          } else if (!ri.empty()) {
            const auto &mid = adj[static_cast<std::size_t>(ri.back())];
            const int j = mid[static_cast<std::size_t>(vertex(rng)) % mid.size()];
            if (j != i) {
              ri.push_back(j);
            }
          }
        }
      }
    }
  }
}

TEST_CASE("ring entries on ice lattices match the oracle", "[ring_oracle]") {
  std::mt19937 rng(3u);
  SECTION("ice Ih and Ic graphs at every depth") {
    for (const Frame &f : {iceIh(2, 2, 2), iceIc(2)}) {
      const Adjacency adj = oracle::cutoffGraph(f.xyz, f.box9.data(), 3.2);
      for (int depth = 3; depth <= 12; ++depth) {
        checkAllEntries(adj, depth, rng);
      }
    }
  }
  SECTION("ice II and ice III graphs at every depth") {
    const Frame two = replicate(kIceIIXyz, kIceIIN, kIceIIBox9, 2, 2, 2);
    const Frame three = replicate(kIceIIIXyz, kIceIIIN, kIceIIIBox9, 2, 2, 2);
    for (const Frame *f : {&two, &three}) {
      const Adjacency adj = oracle::cutoffGraph(f->xyz, f->box9.data(), 3.2);
      for (const auto &row : adj) {
        REQUIRE(row.size() == 4);
      }
      for (int depth = 3; depth <= 10; ++depth) {
        checkAllEntries(adj, depth, rng);
      }
    }
  }
  SECTION("more than 4096 vertices") {
    const Frame f = iceIh(9, 8, 8);
    const Adjacency adj = oracle::cutoffGraph(f.xyz, f.box9.data(), 3.2);
    for (int depth : {6, 7, 8}) {
      checkAllEntries(adj, depth, rng);
    }
  }
}

namespace {

std::vector<int> symmetryCensus(const Frame &f, double cutoff, int depth) {
  seams_system *sys = systemOf(f);
  seams_rings *rings = nullptr;
  REQUIRE(seams_rings_find(sys, cutoff, 0, SEAMS_GRAPH_CUTOFF, depth, &rings) == SEAMS_OK);
  const int n = seams_rings_count(rings);
  std::vector<int> sym(static_cast<std::size_t>(n));
  REQUIRE(seams_rings_labels(rings, nullptr, sym.data()) == SEAMS_OK);
  std::vector<int> census(12, 0);
  for (int code : sym) {
    REQUIRE(code >= 0);
    REQUIRE(code < 12);
    ++census[static_cast<std::size_t>(code)];
  }
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  return census;
}

// Hexagonal only, double-diamond only, neither, both.
std::array<int, 4> cageCounts(const Frame &f, int graph) {
  seams_system *sys = systemOf(f);
  seams_cage_set *cages = nullptr;
  REQUIRE(seams_cages_find(sys, 3.2, 0, graph, 0, 0, &cages) == SEAMS_OK);
  std::array<int, 4> counts{};
  REQUIRE(seams_cages_counts(cages, &counts[0], &counts[1], &counts[2], &counts[3]) ==
          SEAMS_OK);
  seams_cages_destroy(cages);
  seams_system_destroy(sys);
  return counts;
}

} // namespace

TEST_CASE("every six-ring of ideal ice I is a cage face", "[ring_oracle][c_api]") {
  // Ice Ih: each molecule sits in one puckered basal layer, whose hexagons
  // number half its molecules; the other 3N/2 of the 2N six-rings are
  // prismatic. Ice Ic: every six-ring is a double-diamond face.
  for (int threads : {1, 4}) {
    REQUIRE(seams_set_num_threads(threads) == SEAMS_OK);
    for (const auto dims : {std::array<int, 3>{4, 2, 2}, std::array<int, 3>{4, 4, 4}}) {
      const Frame ih = iceIh(dims[0], dims[1], dims[2]);
      const int n = static_cast<int>(ih.xyz.size() / 3);
      INFO("ice Ih " << n << " threads " << threads);
      const auto census = symmetryCensus(ih, 3.2, 6);
      CHECK(census[SEAMS_RING_HC_BASAL] == n / 2);
      CHECK(census[SEAMS_RING_HC_PRISMATIC] == 3 * n / 2);
      CHECK(census[SEAMS_RING_HEXAGON] == 0);
      CHECK(cageCounts(ih, SEAMS_GRAPH_CUTOFF) == std::array<int, 4>{n, 0, 0, 0});
      CHECK(cageCounts(ih, SEAMS_GRAPH_MUTUAL) == std::array<int, 4>{n, 0, 0, 0});
    }
    for (int side : {2, 3}) {
      const Frame ic = iceIc(side);
      const int n = static_cast<int>(ic.xyz.size() / 3);
      INFO("ice Ic " << n << " threads " << threads);
      const auto census = symmetryCensus(ic, 3.2, 6);
      CHECK(census[SEAMS_RING_DDC] == 2 * n);
      CHECK(census[SEAMS_RING_HEXAGON] == 0);
      CHECK(cageCounts(ic, SEAMS_GRAPH_CUTOFF) == std::array<int, 4>{0, n, 0, 0});
      CHECK(cageCounts(ic, SEAMS_GRAPH_MUTUAL) == std::array<int, 4>{0, n, 0, 0});
    }
  }
  REQUIRE(seams_set_num_threads(0) == SEAMS_OK);
}

TEST_CASE("prism and clathrate cages on hand-built frames", "[ring_oracle][c_api]") {
  SECTION("a periodic column of m stacked k-gons holds m perfect k-prisms") {
    // Seven layers or more, so the axial chain closes no ring of six or fewer.
    for (int k = 4; k <= 6; ++k) {
      for (int layers : {7, 8, 10}) {
        INFO(k << "-gon column of " << layers << " layers");
        const Frame f = prismColumn(k, layers, 2.76);
        seams_system *sys = systemOf(f);
        int perfect[8] = {};
        int imperfect[8] = {};
        REQUIRE(seams_prisms(sys, 3.2, 0, perfect, imperfect, 8) == SEAMS_OK);
        for (int s = 3; s < 8; ++s) {
          INFO("size " << s);
          CHECK(perfect[s] == (s == k ? layers : 0));
          CHECK(imperfect[s] == 0);
        }
        const auto census = symmetryCensus(f, 3.2, 8);
        CHECK(census[SEAMS_RING_PRISM] == layers);
        seams_system_destroy(sys);
      }
    }
  }
  SECTION("closed cages of the sI and sII cells and of one dodecahedron") {
    struct Case {
      Frame f;
      const char *spec;
      std::size_t cages;
      std::size_t atoms;
    };
    // Four cages meet at every vertex. In sI the six 6c waters sit on no
    // dodecahedron; in sII the eight 8a waters sit on dodecahedra only.
    const Frame sI = cubic(kSIxyz, kSI, kSIcell);
    const Frame sII = cubic(kSIIxyz, kSII, kSIIcell);
    for (const Case &c : {Case{sI, "512", 2, 40}, Case{sI, "51262", 6, 46},
                          Case{sII, "512", 16, 136}, Case{sII, "51264", 8, 128},
                          Case{dodecahedron(2.76), "512", 1, 20}}) {
      INFO(c.spec << " on " << c.f.xyz.size() / 3 << " sites");
      seams_system *sys = systemOf(c.f);
      std::size_t cages = 0;
      std::size_t atoms = 0;
      REQUIRE(seams_signature_cages(sys, 3.2, 0, SEAMS_GRAPH_CUTOFF, c.spec, &cages, &atoms) ==
              SEAMS_OK);
      CHECK(cages == c.cages);
      CHECK(atoms == c.atoms);
      seams_system_destroy(sys);
    }
  }
}

TEST_CASE("seams_rings_find matches the oracle on periodic cells", "[ring_oracle][c_api]") {
  for (int threads : {1, 4}) {
    REQUIRE(seams_set_num_threads(threads) == SEAMS_OK);
    for (int depth : {6, 8, 10, 12}) {
      checkFrame(iceIh(2, 2, 2), 3.2, depth, "ice Ih 2x2x2");
      checkFrame(iceIc(2), 3.2, depth, "ice Ic 2x2x2");
      checkFrame(replicate(kIceIIXyz, kIceIIN, kIceIIBox9, 2, 2, 2), 3.2, depth,
                 "ice II 2x2x2 triclinic");
      checkFrame(replicate(kIceIIIXyz, kIceIIIN, kIceIIIBox9, 2, 2, 2), 3.2, depth,
                 "ice III 2x2x2");
      checkFrame(cubic(kSIxyz, kSI, kSIcell), 3.2, depth, "sI");
    }
    // One cell across: rings close through the periodic boundary, and an
    // ice III or ice II cell is narrower than twice the cutoff.
    for (int depth : {6, 8, 10}) {
      checkFrame(iceIh(1, 1, 1), 3.2, depth, "ice Ih one cell");
      checkFrame(iceIc(1), 3.2, depth, "ice Ic one cell");
      checkFrame(replicate(kIceIIIXyz, kIceIIIN, kIceIIIBox9, 1, 1, 1), 3.2, depth,
                 "ice III one cell");
      checkFrame(replicate(kIceIIXyz, kIceIIN, kIceIIBox9, 1, 1, 1), 3.2, depth,
                 "ice II one cell");
    }
    checkFrame(prismColumn(4, 3, 2.76), 3.2, 8, "square column of three layers");
    checkFrame(prismColumn(5, 9, 2.76), 3.2, 10, "pentagonal column");
    checkFrame(prismColumn(6, 9, 2.76), 3.2, 10, "hexagonal column");
    checkFrame(dodecahedron(2.76), 3.2, 8, "dodecahedron");
  }
  REQUIRE(seams_set_num_threads(0) == SEAMS_OK);
}

TEST_CASE("seams_rings_find keeps frame indices under a species filter and on the mutual graph",
          "[ring_oracle][c_api]") {
  // Type-2 decoys at every bond midpoint would bond to both ends if the
  // filter let them in. On ideal ice I the four nearest neighbours are the
  // only partners inside the cutoff, so the mutual graph is the cutoff graph.
  for (const auto &[f, name] : {std::pair<Frame, const char *>{iceIh(2, 2, 2), "ice Ih"},
                                std::pair<Frame, const char *>{iceIc(2), "ice Ic"},
                                std::pair<Frame, const char *>{
                                    replicate(kIceIIXyz, kIceIIN, kIceIIBox9, 2, 2, 2),
                                    "ice II"}}) {
    INFO(name);
    const Adjacency adj = oracle::cutoffGraph(f.xyz, f.box9.data(), 3.2);
    const int n = static_cast<int>(adj.size());
    const auto h = oracle::lattice(f.box9.data());
    Frame mixed;
    mixed.box9 = f.box9;
    std::vector<int> types;
    std::vector<int> oxygenAt;
    for (int i = 0; i < n; ++i) {
      oxygenAt.push_back(static_cast<int>(types.size()));
      types.push_back(1);
      mixed.xyz.insert(mixed.xyz.end(), f.xyz.begin() + 3 * i, f.xyz.begin() + 3 * i + 3);
      for (int j : adj[static_cast<std::size_t>(i)]) {
        if (j < i) {
          continue;
        }
        // The image of j nearest to i.
        double best = 1e300;
        std::array<double, 3> mid{};
        for (int a = -1; a <= 1; ++a) {
          for (int b = -1; b <= 1; ++b) {
            for (int c = -1; c <= 1; ++c) {
              std::array<double, 3> d{};
              double d2 = 0;
              for (int k = 0; k < 3; ++k) {
                d[static_cast<std::size_t>(k)] =
                    f.xyz[static_cast<std::size_t>(3 * j + k)] -
                    f.xyz[static_cast<std::size_t>(3 * i + k)] + a * h[0][k] + b * h[1][k] +
                    c * h[2][k];
                d2 += d[static_cast<std::size_t>(k)] * d[static_cast<std::size_t>(k)];
              }
              if (d2 < best) {
                best = d2;
                for (int k = 0; k < 3; ++k) {
                  mid[static_cast<std::size_t>(k)] =
                      f.xyz[static_cast<std::size_t>(3 * i + k)] +
                      0.5 * d[static_cast<std::size_t>(k)];
                }
              }
            }
          }
        }
        types.push_back(2);
        mixed.xyz.insert(mixed.xyz.end(), mid.begin(), mid.end());
      }
    }
    Rings want = oracle::primitiveRings(adj, 8);
    for (auto &r : want) {
      for (int &v : r) {
        v = oxygenAt[static_cast<std::size_t>(v)];
      }
    }
    want = oracle::canonicalSorted(want);
    seams_system *sys = seams_system_create();
    REQUIRE(seams_system_set_atoms(sys, mixed.xyz.data(), mixed.xyz.size() / 3,
                                   mixed.box9.data(), 9, types.data(), nullptr) == SEAMS_OK);
    seams_rings *rings = nullptr;
    REQUIRE(seams_rings_find(sys, 3.2, 1, SEAMS_GRAPH_CUTOFF, 8, &rings) == SEAMS_OK);
    REQUIRE_RINGS(ringsOf(rings), want, adj, 8, "species filter");
    seams_rings_destroy(rings);
    seams_system_destroy(sys);

    if (std::string(name) != "ice II") {
      seams_system *pure = systemOf(f);
      REQUIRE(seams_rings_find(pure, 3.2, 0, SEAMS_GRAPH_MUTUAL, 8, &rings) == SEAMS_OK);
      REQUIRE_RINGS(ringsOf(rings), oracle::primitiveRings(adj, 8), adj, 8, "mutual graph");
      seams_rings_destroy(rings);
      seams_system_destroy(pure);
    }
  }
}
