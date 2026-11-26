// Copyright (c) d-SEAMS developers
//
// Adversarial bands for the Euclidean minimum image, scored against the
// binary128 reference in oracle.hpp.
//
// usage: mic_oracle [compact-cells] [seed] [-v]
//
// meson test draws 800 cells on each compact band and scores every
// displacement shape on that cell. Needle and cliff bands draw 8 cells,
// except the long-edge needle, which draws 48 cells and one displacement
// each: that is the sample that still reaches a descent landing outside
// the Voronoi cell and the basis swap in the coset search. The basis is
// built once per cell. A coverage build compiles the minimum-image
// headers at -O1 without inlining and with atomic arc counters, and that
// first call is the runtime. Later shapes on the same cell reuse the
// cached basis. A fresh cell per pair repeats edges this draw already
// took. `mic_oracle 100000` raises the compact-band cell count.
// needle-1e10-selling-cap is not part of that draw. It is one fixed cell
// with the origin at zero: lx 1e10, ly 0.94595037089926481,
// lz 0.41367195128549233, xy 0.22894346172539071,
// xz 0.18843458615172756, yz 0.27052973178254824. Selling reduction of
// that cell keeps replacing an acute pair with a shorter one, so the
// obtuse cleanup in snapToLattice runs out its iteration cap. The other
// bands do not reach that exit. The band then scores four fixed
// displacements, the first of them (1.3, 0.2, -0.1), through the same
// binary128 closest and the same scoreOne limit as every other band.
// A gate rebuilds the binary at other flags from the build's own
// commands: the mic_oracle entry of build/compile_commands.json with, for
// example, -O3 -march=native appended, linked as `ninja -t commands
// tests/mic_oracle` links it. The minimum-image code is inline in the
// headers, so it is compiled at those flags.
//
// Each pair scores gen::mic::euclideanIfShorter, gen::refineEuclidean and
// gen::exactDelta. The unit is 2^-52 times the pair's evaluation scale: the
// largest component of the displacement, its fractional wrap, the exact
// image, the translation between them, and the cell rows times the
// fractional coordinates of the displacement and of its wrap
// (translationScale and fractionalScale in oracle.hpp). The fractional wrap
// that every function can start from adds up numbers that large, so
// binary64 places no image closer. A result is a
// violation when its length differs from the exact minimum by more than
// kUnits units, it is not a lattice image of the displacement to kUnits
// units (offLattice), or, unless another image is within kUnits units of
// the minimum length, a component is more than kUnits units from the exact
// image. -v prints the first violations of each band and function in full.
//
// Throws are counted. A throw on an in-range band is a violation unless the
// wrapped length is at least 2^50 times the shortest lattice vector. On the
// cliff band (lx 2e16, 1e17 and 1e160 with ly 1, lz 0.4 and tilts 0.3)
// euclideanIfShorter and refineEuclidean must throw mic::RangeError for
// every pair, and exactDelta may only answer inside the Smith ball. The
// off-lattice band passes euclideanIfShorter a finite candidate off the
// lattice, which must never come back. With minimage linked,
// refineEuclidean starts from minimage's image through
// mic::euclideanFromImage. The minimage band covers needle cells where
// minimage 0.1.4 returns images off the lattice, and the tilt-1e-14 band
// cells where its image drops the tilt and is off the lattice by it.
#include <generic.hpp>
#include <neighbours.hpp>

#include "oracle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <functional>
#include <string>
#include <vector>

namespace {

using micoracle::Cell;
using micoracle::Q;

struct Rng {
  std::uint64_t s;
  std::uint64_t next() {
    std::uint64_t z = (s += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
  }
  double uniform() { return (double)(next() >> 11) * 0x1p-53; }
  double in(double lo, double hi) { return lo + (hi - lo) * uniform(); }
  int pick(int lo, int hi) { return lo + (int)(next() % (std::uint64_t)(hi - lo + 1)); }
  double sign() { return (next() & 1) ? 1.0 : -1.0; }
};

Cell legal(Rng &g, double scale) {
  Cell c{};
  c.lx = scale * g.in(1, 10);
  c.ly = scale * g.in(1, 10);
  c.lz = scale * g.in(1, 10);
  c.xy = g.in(-0.5, 0.5) * c.lx;
  c.xz = g.in(-0.5, 0.5) * c.lx;
  c.yz = g.in(-0.5, 0.5) * c.ly;
  return c;
}

Cell needle(double lx, Rng &g, bool jitter) {
  const double j = jitter ? 1.0 : 0.0;
  return {lx, 1.0 + j * g.in(-0.1, 0.1), 0.4 + j * g.in(-0.05, 0.05),
          0.3 + j * g.in(-0.05, 0.05), 0.3 + j * g.in(-0.05, 0.05),
          0.3 + j * g.in(-0.05, 0.05)};
}

struct Band {
  const char *name;
  std::function<Cell(Rng &)> cell;
  enum { Mixed, Window } sampler;
  enum { InRange, OffLattice, Cliff } kind;
};

std::vector<Band> bands() {
  std::vector<Band> out;
  out.push_back({"random-legal", [](Rng &g) { return legal(g, 1.0); }, Band::Mixed,
                 Band::InRange});
  out.push_back({"tilt-corners",
                 [](Rng &g) {
                   Cell c = legal(g, 1.0);
                   c.xy = g.sign() * 0.5 * c.lx;
                   c.xz = g.sign() * 0.5 * c.lx;
                   c.yz = g.sign() * 0.5 * c.ly;
                   return c;
                 },
                 Band::Mixed, Band::InRange});
  out.push_back({"aspect-1e3",
                 [](Rng &g) {
                   Cell c = legal(g, 1.0);
                   const int axis = g.pick(0, 2);
                   (axis == 0 ? c.lx : axis == 1 ? c.ly : c.lz) *= 1e3;
                   c.xy = g.in(-0.5, 0.5) * c.lx;
                   c.xz = g.in(-0.5, 0.5) * c.lx;
                   c.yz = g.in(-0.5, 0.5) * c.ly;
                   return c;
                 },
                 Band::Mixed, Band::InRange});
  out.push_back({"hex-prism",
                 [](Rng &g) {
                   const double a = g.in(5, 40);
                   Cell c{a, a * std::sqrt(3.0) / 2.0, a * g.in(0.5, 5), a / 2.0,
                          g.in(-0.01, 0.01) * a, g.in(-0.01, 0.01) * a};
                   return c;
                 },
                 Band::Mixed, Band::InRange});
  out.push_back({"cof-skew",
                 [](Rng &g) {
                   Cell c{g.in(10, 40), g.in(10, 40), g.in(1, 4), 0, 0, 0};
                   c.xy = g.in(-0.5, 0.5) * c.lx;
                   c.xz = g.in(-0.5, 0.5) * c.lx;
                   c.yz = g.in(-0.5, 0.5) * c.ly;
                   return c;
                 },
                 Band::Mixed, Band::InRange});
  out.push_back({"needle-1e8..1e15",
                 [](Rng &g) {
                   Cell c = needle(std::pow(10.0, g.pick(8, 15)), g, true);
                   c.xy = g.in(-0.5, 0.5) * c.ly;
                   c.xz = g.in(-0.5, 0.5) * c.ly;
                   return c;
                 },
                 Band::Mixed, Band::InRange});
  out.push_back({"thin-ly",
                 [](Rng &g) {
                   Cell c = legal(g, 1.0);
                   c.ly = 1e-3 * c.lx;
                   c.yz = g.in(-0.5, 0.5) * c.ly;
                   return c;
                 },
                 Band::Mixed, Band::InRange});
  out.push_back({"shear-x100",
                 [](Rng &g) {
                   Cell c = legal(g, 1.0);
                   c.xy = g.in(-50, 50) * c.lx;
                   c.xz = g.in(-50, 50) * c.lx;
                   c.yz = g.in(-50, 50) * c.ly;
                   return c;
                 },
                 Band::Mixed, Band::InRange});
  out.push_back({"scale-1e-10", [](Rng &g) { return legal(g, 1e-10); }, Band::Mixed,
                 Band::InRange});
  out.push_back({"scale-1e+10", [](Rng &g) { return legal(g, 1e10); }, Band::Mixed,
                 Band::InRange});
  out.push_back({"tilt-1e-14",
                 [](Rng &g) {
                   Cell c = legal(g, 1.0);
                   c.xy = g.sign() * 1e-14 * c.lx;
                   c.xz = g.sign() * 1e-14 * c.lx;
                   c.yz = g.sign() * 1e-14 * c.ly;
                   return c;
                 },
                 Band::Mixed, Band::InRange});
  out.push_back({"near-orthogonal",
                 [](Rng &g) {
                   Cell c = legal(g, 1.0);
                   c.xy = g.in(-1e-3, 1e-3) * c.lx;
                   c.xz = g.in(-1e-3, 1e-3) * c.lx;
                   c.yz = g.in(-1e-3, 1e-3) * c.ly;
                   return c;
                 },
                 Band::Mixed, Band::InRange});
  out.push_back({"minimage-needle",
                 [](Rng &g) {
                   const double lx[3] = {1e6, 1e12, 1e15};
                   return needle(lx[g.pick(0, 2)], g, false);
                 },
                 Band::Window, Band::InRange});
  out.push_back({"off-lattice-candidate", [](Rng &g) { return legal(g, 1.0); }, Band::Mixed,
                 Band::OffLattice});
  out.push_back({"cliff-2e16-1e17-1e160",
                 [](Rng &g) {
                   const double lx[3] = {2e16, 1e17, 1e160};
                   return needle(lx[g.pick(0, 2)], g, false);
                 },
                 Band::Window, Band::Cliff});
  // Appended so the earlier bands keep their seeds. One cell, not a draw:
  // these literals are the cell whose obtuse cleanup runs out the iteration
  // cap in snapToLattice.
  out.push_back({"needle-1e10-selling-cap",
                 [](Rng &) {
                   return Cell{1e10,
                               0.94595037089926481,
                               0.41367195128549233,
                               0.22894346172539071,
                               0.18843458615172756,
                               0.27052973178254824};
                 },
                 Band::Mixed, Band::InRange});
  return out;
}

std::array<double, 3> cellTimes(const Cell &c, double s0, double s1, double s2) {
  return {c.lx * s0 + c.xy * s1 + c.xz * s2, c.ly * s1 + c.yz * s2, c.lz * s2};
}

// Displacements for needle-1e10-selling-cap. The cell is what exhausts the
// obtuse cleanup; these four stay near the origin and are scored like any
// other pair.
const double kSellingCapRaw[][3] = {
    {1.3, 0.2, -0.1},
    {0.7, -0.35, 0.2},
    {-1.1, 0.15, -0.45},
    {2.4, 0.55, 0.3},
};
constexpr int kSellingCapDisp = 4;

// Four displacement shapes in turn: uniform in a few cells, at the wrap's
// cell faces, a short image under a long lattice vector, and half a lattice
// vector, where two images tie.
std::array<double, 3> mixed(const Cell &c, double shortest, Rng &g, int which) {
  switch (which % 4) {
  case 0:
    return cellTimes(c, g.in(-0.5, 0.5) + g.pick(-2, 2), g.in(-0.5, 0.5) + g.pick(-2, 2),
                     g.in(-0.5, 0.5) + g.pick(-2, 2));
  case 1: {
    double s[3];
    for (double &v : s) {
      v = (g.next() & 1) ? g.sign() * 0.5 + g.in(-1e-12, 1e-12) : g.in(-0.5, 0.5);
    }
    return cellTimes(c, s[0], s[1], s[2]);
  }
  case 2: {
    auto v = cellTimes(c, g.pick(-3, 3), g.pick(-3, 3), g.pick(-3, 3));
    for (double &x : v) {
      x += shortest * g.in(-0.5, 0.5);
    }
    return v;
  }
  default: {
    auto v = cellTimes(c, 0.5 * g.pick(-3, 3), 0.5 * g.pick(-3, 3), 0.5 * g.pick(-3, 3));
    for (double &x : v) {
      x += shortest * g.in(-1e-9, 1e-9);
    }
    return v;
  }
  }
}

// The difference of two points spread over b and c and a strip 3 long in x.
std::array<double, 3> window(const Cell &c, Rng &g) {
  const auto p = cellTimes(c, 0, g.uniform(), g.uniform());
  const auto q = cellTimes(c, 0, g.uniform(), g.uniform());
  return {(p[0] + g.in(0, 3)) - (q[0] + g.in(0, 3)), p[1] - q[1], p[2] - q[2]};
}

// Each function rounds a handful of times at the evaluation scale, the
// wrap's products and sums and then a few lattice steps, each time by at
// most half a unit.
constexpr double kUnits = 8;

bool verbose = false;

struct Score {
  long pairs = 0, viol = 0, offLattice = 0, throws = 0, spurious = 0, missedThrow = 0, other = 0;
  double worstLen = 0, worstOff = 0, worstVec = 0;
  std::string first;
};

void note(Score &s, const std::string &what, const std::string &detail) {
  s.viol++;
  if (s.first.empty()) {
    s.first = what;
  }
  if (verbose && s.viol <= 3) {
    std::printf("  %s\n%s", what.c_str(), detail.c_str());
  }
}

std::string triple(const char *name, double x, double y, double z) {
  char buf[128];
  std::snprintf(buf, sizeof buf, "    %-5s %.17g %.17g %.17g\n", name, x, y, z);
  return buf;
}

using Fn = std::function<std::array<double, 3>(const gen::FracBox &, const std::array<double, 3> &,
                                               const std::array<double, 3> &)>;

struct Pair {
  Cell cell;
  gen::FracBox box;
  std::array<double, 3> raw;
  std::array<double, 3> wrap;
  micoracle::Image exact;
  double shortest;
  double unit;
  bool tie;
};

void scoreOne(Score &s, const Band &band, const Pair &p, const Fn &fn, bool mustThrow,
              const std::array<double, 3> *candidate) {
  s.pairs++;
  const std::array<double, 3> &given = candidate != nullptr ? *candidate : p.wrap;
  const double wlen = std::sqrt(p.wrap[0] * p.wrap[0] + p.wrap[1] * p.wrap[1] + p.wrap[2] * p.wrap[2]);
  std::array<double, 3> d{};
  char where[160];
  std::snprintf(where, sizeof where, "%s lx %.17g ly %.17g raw (%.17g, %.17g, %.17g)", band.name,
                p.cell.lx, p.cell.ly, p.raw[0], p.raw[1], p.raw[2]);
  char cellText[200];
  std::snprintf(cellText, sizeof cellText, "    cell  %.17g %.17g %.17g %.17g %.17g %.17g\n",
                p.cell.lx, p.cell.ly, p.cell.lz, p.cell.xy, p.cell.xz, p.cell.yz);
  std::string detail = std::string(cellText) + triple("raw", p.raw[0], p.raw[1], p.raw[2]) +
                       triple("given", given[0], given[1], given[2]) +
                       triple("exact", (double)p.exact.d.x, (double)p.exact.d.y,
                              (double)p.exact.d.z);
  try {
    d = fn(p.box, p.raw, given);
  } catch (const gen::mic::RangeError &) {
    s.throws++;
    if (!mustThrow && band.kind != Band::Cliff && wlen < 0x1p50 * p.shortest) {
      s.spurious++;
      note(s, std::string("spurious RangeError: ") + where, detail);
    }
    return;
  } catch (const std::exception &ex) {
    s.other++;
    note(s, std::string("unexpected exception '") + ex.what() + "': " + where, detail);
    return;
  }
  detail += triple("got", d[0], d[1], d[2]);
  if (mustThrow) {
    s.missedThrow++;
    note(s, std::string("returned instead of RangeError: ") + where, detail);
    return;
  }
  const double raw[3] = {p.raw[0], p.raw[1], p.raw[2]};
  const double got[3] = {d[0], d[1], d[2]};
  const double off = micoracle::offLattice(p.cell, raw, got);
  s.worstOff = std::max(s.worstOff, off);
  if (!(off <= kUnits)) {
    s.offLattice++;
    note(s, std::string("not a lattice image: ") + where, detail);
    return;
  }
  const Q got2 = (Q)d[0] * d[0] + (Q)d[1] * d[1] + (Q)d[2] * d[2];
  const double glen = (double)sqrtq(got2);
  const double elen = (double)sqrtq(p.exact.n2);
  const double eu = (glen - elen) / p.unit;
  s.worstLen = std::max(s.worstLen, std::fabs(eu));
  if (!(std::fabs(eu) <= kUnits)) {
    note(s, std::string("length off the exact minimum: ") + where, detail);
    return;
  }
  if (!p.tie) {
    const double ex[3] = {(double)p.exact.d.x, (double)p.exact.d.y, (double)p.exact.d.z};
    double vec = 0;
    for (int k = 0; k < 3; k++) {
      vec = std::max(vec, std::fabs(d[k] - ex[k]) / p.unit);
    }
    s.worstVec = std::max(s.worstVec, vec);
    if (!(vec <= kUnits)) {
      note(s, std::string("image vector off the exact image: ") + where, detail);
    }
  }
}

int selftest() {
  int bad = 0;
  Rng g{12345};
  for (int t = 0; t < 2000; t++) {
    const Cell c = legal(g, 1.0);
    const auto r = micoracle::reduce(c);
    const auto raw = cellTimes(c, g.in(-3, 3), g.in(-3, 3), g.in(-3, 3));
    const micoracle::V w{(Q)raw[0], (Q)raw[1], (Q)raw[2]};
    const auto img = micoracle::closest(c, r, w);
    Q brute = (Q)INFINITY;
    for (int a = -8; a <= 8; a++) {
      for (int b = -8; b <= 8; b++) {
        for (int k = -8; k <= 8; k++) {
          const micoracle::I n[3] = {a, b, k};
          const auto d = w - micoracle::combo(c, n);
          brute = std::min(brute, micoracle::dot(d, d));
        }
      }
    }
    if (!img.ok || img.n2 != brute) {
      bad++;
    }
    const double rawd[3] = {raw[0], raw[1], raw[2]};
    const double imgd[3] = {(double)img.d.x, (double)img.d.y, (double)img.d.z};
    if (!(micoracle::offLattice(c, rawd, imgd) <= 2)) {
      bad++;
    }
    const double moved[3] = {imgd[0] + 1e-6 * c.lz, imgd[1], imgd[2]};
    if (!(micoracle::offLattice(c, rawd, moved) > 1e3)) {
      bad++;
    }
  }
  for (int t = 0; t < 500; t++) {
    Cell c = legal(g, 1.0);
    c.xy = c.xz = c.yz = 0;
    const auto r = micoracle::reduce(c);
    const double raw[3] = {g.in(-20, 20), g.in(-20, 20), g.in(-20, 20)};
    const auto img = micoracle::closest(c, r, micoracle::fromDouble(raw));
    const double len[3] = {c.lx, c.ly, c.lz};
    Q want = 0;
    for (int k = 0; k < 3; k++) {
      const Q v = (Q)raw[k] - (Q)len[k] * roundq((Q)raw[k] / (Q)len[k]);
      want += v * v;
    }
    if (!img.ok || fabsq(img.n2 - want) > (Q)1e-28 * want) {
      bad++;
    }
  }
  std::printf("selftest: %d failures\n", bad);
  return bad;
}

}  // namespace

int main(int argc, char **argv) {
  // Compact bands saturate their reduction and search edges by a few
  // hundred cells once every displacement shape is scored on the cell.
  // Needle and cliff cells are the slow basis builds, and eight of each
  // already cover those bands. argv replaces the compact-band count.
  constexpr long kCompactCells = 800;
  constexpr long kHeavyCells = 8;
  constexpr long kNeedleCells = 48;
  constexpr int kCompactDisp = 8;
  constexpr int kHeavyDisp = 4;
  const long perBand = argc > 1 ? std::strtol(argv[1], nullptr, 10) : kCompactCells;
  const std::uint64_t seed = argc > 2 ? std::strtoull(argv[2], nullptr, 10) : 20261125;
  verbose = argc > 3 && std::string(argv[3]) == "-v";
  if (perBand <= 0 || argc > 4 || (argc > 3 && !verbose)) {
    std::fprintf(stderr, "usage: mic_oracle [compact-cells] [seed] [-v]\n");
    return 2;
  }
#ifdef SEAMS_HAS_MINIMAGE
  std::printf("minimage linked; refineEuclidean starts from its image\n");
#else
  std::printf("minimage not linked\n");
#endif
  int failed = selftest() != 0;

  const char *names[3] = {"euclideanIfShorter", "refineEuclidean", "exactDelta"};
  const Fn fns[3] = {
      [](const gen::FracBox &b, const std::array<double, 3> &raw, const std::array<double, 3> &w) {
        return gen::mic::euclideanIfShorter(b.lx, b.ly, b.lz, b.xy, b.xz, b.yz, raw, w);
      },
      [](const gen::FracBox &b, const std::array<double, 3> &raw, const std::array<double, 3> &w) {
        return gen::refineEuclidean(b.lx, b.ly, b.lz, b.xy, b.xz, b.yz, raw, w);
      },
      [](const gen::FracBox &b, const std::array<double, 3> &raw, const std::array<double, 3> &) {
        return gen::exactDelta(b, raw[0], raw[1], raw[2], 0.0, 0.0, 0.0);
      },
  };
  long total = 0, violations = 0, oracleFail = 0;
  std::printf("units of 2^-52 times the evaluation scale; a violation is above %g\n", kUnits);
  std::printf("%-24s %-19s %8s %6s %6s %7s %8s %9s %9s %9s %9s\n", "band", "function", "pairs",
              "viol", "offlat", "throws", "spurious", "worstLen", "worstOff", "worstVec",
              "missThrow");
  const auto all = bands();
  for (std::size_t bi = 0; bi < all.size(); bi++) {
    const Band &band = all[bi];
    Rng cells{seed + 0x1000 * bi};
    Rng dispRng{seed + 0x8000 * bi + 17};
    Score score[3];
    const bool needle = std::string(band.name) == "needle-1e8..1e15";
    const bool sellingCap = std::string(band.name) == "needle-1e10-selling-cap";
    const bool heavy = !sellingCap && (band.sampler == Band::Window ||
                                       std::string(band.name).find("needle") != std::string::npos);
    // The long-edge needle shares one generator between the cell and its
    // single displacement, matching the draw that leaves the Voronoi cell.
    Rng &disp = needle ? cells : dispRng;
    const long nCells = sellingCap ? 1
                                   : needle ? std::min(perBand, kNeedleCells)
                                            : heavy ? std::min(perBand, kHeavyCells) : perBand;
    const int nDisp = sellingCap ? kSellingCapDisp : needle ? 1 : heavy ? kHeavyDisp : kCompactDisp;
    for (long ci = 0; ci < nCells; ci++) {
      Pair held{};
      const Cell want = band.cell(cells);
      molSys::PointCloud<molSys::Point<double>, double> cloud;
      if (!nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, 0.0, 0.0, 0.0, want.lx, want.ly,
                                         want.lz, want.xy, want.xz, want.yz)) {
        oracleFail++;
        continue;
      }
      held.box = gen::makeFracBox(cloud);
      held.cell = {held.box.lx, held.box.ly, held.box.lz, held.box.xy, held.box.xz, held.box.yz};
      const auto basis = micoracle::reduce(held.cell);
      if (!held.box.ok || !basis.ok) {
        oracleFail++;
        continue;
      }
      held.shortest = (double)sqrtq(basis.gn[0]);
      for (int shape = 0; shape < nDisp; shape++) {
        Pair p = held;
        // The long-edge needle draws one displacement per cell and cycles
        // the four shapes with the cell index, the same sequence as a
        // pair-per-cell run. The selling-cap cell uses its fixed list.
        const int which = needle ? static_cast<int>(ci) : shape;
        if (sellingCap) {
          p.raw = {kSellingCapRaw[shape][0], kSellingCapRaw[shape][1], kSellingCapRaw[shape][2]};
        } else {
          p.raw = band.sampler == Band::Window ? window(p.cell, disp)
                                              : mixed(p.cell, p.shortest, disp, which);
        }
        p.wrap = gen::fracDelta(p.box, p.raw[0], p.raw[1], p.raw[2], 0.0, 0.0, 0.0);
        p.exact = micoracle::closest(p.cell, basis, micoracle::fromDouble(p.raw.data()));
        if (!p.exact.ok) {
          oracleFail++;
          continue;
        }
        {
          const double ex[3] = {(double)p.exact.d.x, (double)p.exact.d.y, (double)p.exact.d.z};
          const auto t = micoracle::translationScale(p.cell, p.exact.n);
          const auto fr = micoracle::fractionalScale(p.cell, p.raw.data());
          const auto fw = micoracle::fractionalScale(p.cell, p.wrap.data());
          double scale = 0;
          for (int k = 0; k < 3; k++) {
            scale = std::max({scale, std::fabs(p.raw[k]), std::fabs(p.wrap[k]), std::fabs(ex[k]),
                              t[k], fr[k], fw[k]});
          }
          p.unit = scale > 0 ? scale * 0x1p-52 : 0x1p-1074;
        }
        // Another image whose length is within kUnits of the minimum is a tie:
        // binary64 cannot tell the two apart, so only the length is scored.
        {
          const Q elen = sqrtq(p.exact.n2);
          const Q tol = (Q)kUnits * (Q)p.unit;
          p.tie = micoracle::nearTie(basis, micoracle::fromDouble(p.raw.data()), p.exact,
                                     2 * elen * tol + tol * tol);
        }
        std::array<double, 3> moved = p.wrap;
        const std::array<double, 3> *candidate = nullptr;
        if (band.kind == Band::OffLattice) {
          for (double &x : moved) {
            x += p.shortest * disp.in(-1e-3, 1e-3);
          }
          candidate = &moved;
        }
        for (int f = 0; f < 3; f++) {
          const bool smith = gen::smithInside(p.box, p.wrap[0] * p.wrap[0] + p.wrap[1] * p.wrap[1] +
                                                         p.wrap[2] * p.wrap[2]);
          const bool mustThrow = band.kind == Band::Cliff && !(f == 2 && smith);
          if (band.kind == Band::OffLattice && f != 0) {
            continue;
          }
          scoreOne(score[f], band, p, fns[f], mustThrow, candidate);
        }
      }
    }
    for (int f = 0; f < 3; f++) {
      const Score &s = score[f];
      if (s.pairs == 0) {
        continue;
      }
      std::printf("%-24s %-19s %8ld %6ld %6ld %7ld %8ld %9.3g %9.3g %9.3g %9ld\n", band.name,
                  names[f], s.pairs, s.viol, s.offLattice, s.throws, s.spurious, s.worstLen,
                  s.worstOff, s.worstVec, s.missedThrow);
      if (!s.first.empty()) {
        std::printf("  first: %s\n", s.first.c_str());
      }
      total += s.pairs;
      violations += s.viol;
    }
  }
  std::printf("total %ld scored, %ld violations, %ld pairs the oracle could not set up\n", total,
              violations, oracleFail);
  failed |= violations != 0 || oracleFail != 0;
  std::printf(failed ? "FAIL\n" : "ok\n");
  return failed ? 1 : 0;
}
