// Copyright (c) d-SEAMS developers
// Smooth collective variables: gradients, continuity, and the 0/1 limit.
//
// Every gradient is checked against a central difference of the value with
// step h = 1e-5 A. The truncation error is h^2 |f'''| / 6. The third
// derivative of the switch used here (r0 = 0.3 A, n = 6, m = 12, taper over
// 0.5 A) stays below 1.2e3 A^-3, which bounds the error of ds/dr by 2e-8;
// a variable adds such terms over the few dozen bonds an atom reaches
// through its rings, each scaled by weights at most one. The rounding error
// is about 1e-16 |f| / h, below 1e-8 for |f| < 1e3. A bound of
// 1e-6 (1 + |g|) is above both, and an analytic gradient that is wrong in
// one term is off by the size of that term, 1e-2 or more here.

#include "ice_ih.hpp"
#include "seams.h"

#include <catch2/catch_test_macros.hpp>

#include <franzblau.hpp>
#include <smooth_cv.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace {

constexpr double kStep = 1e-5;

struct Frame {
  std::vector<double> xyz;
  double box9[9]{};
  std::size_t n = 0;
};

double lengthOf(const Frame &f, int d) { return f.box9[3 + d] - f.box9[d]; }

void wrap(Frame &f) {
  for (std::size_t i = 0; i < f.n; ++i) {
    for (int d = 0; d < 3; ++d) {
      double &x = f.xyz[3 * i + static_cast<std::size_t>(d)];
      const double l = lengthOf(f, d);
      x -= std::floor((x - f.box9[d]) / l) * l;
    }
  }
}

// Ideal ice Ih with every coordinate moved by up to amp, so that bond
// lengths spread over the soft part of the switch.
Frame jitteredIce(int nx, int ny, int nz, double amp, unsigned seed) {
  const IceIhCell cell = iceIhOrtho(nx, ny, nz);
  Frame f;
  f.xyz = cell.xyz;
  std::copy(cell.box9, cell.box9 + 9, f.box9);
  f.n = static_cast<std::size_t>(cell.n);
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(-amp, amp);
  for (double &x : f.xyz) {
    x += u(rng);
  }
  wrap(f);
  return f;
}

// The ice atoms within radius of the cell centre, in a cubic box of side
// box, centred on `at` so that the cluster straddles a periodic boundary.
Frame iceCluster(double radius, double box, std::array<double, 3> at, double amp,
                 unsigned seed) {
  const Frame ice = jitteredIce(3, 3, 3, amp, seed);
  std::array<double, 3> mid{};
  for (int d = 0; d < 3; ++d) {
    mid[static_cast<std::size_t>(d)] = 0.5 * lengthOf(ice, d);
  }
  Frame f;
  for (std::size_t i = 0; i < ice.n; ++i) {
    double r2 = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
      const double x = ice.xyz[3 * i + d] - mid[d];
      r2 += x * x;
    }
    if (r2 < radius * radius) {
      for (std::size_t d = 0; d < 3; ++d) {
        f.xyz.push_back(ice.xyz[3 * i + d] - mid[d] + at[d]);
      }
      ++f.n;
    }
  }
  f.box9[3] = f.box9[4] = f.box9[5] = box;
  wrap(f);
  return f;
}

struct Api {
  seams_system *sys = nullptr;
  seams_cv *cv = nullptr;
  Frame frame;

  Api(const Frame &f, double r0, double d0, int n, int m, double dOn, double dMax,
      int maxRing)
      : frame(f) {
    sys = seams_system_create();
    REQUIRE(sys != nullptr);
    REQUIRE(seams_system_set_atoms(sys, frame.xyz.data(), frame.n, frame.box9, 9, nullptr,
                                   nullptr) == SEAMS_OK);
    REQUIRE(seams_cv_create(r0, d0, n, m, dOn, dMax, maxRing, &cv) == SEAMS_OK);
  }
  ~Api() {
    seams_cv_destroy(cv);
    seams_system_destroy(sys);
  }
  Api(const Api &) = delete;
  Api &operator=(const Api &) = delete;

  void rebuild(double skin) {
    REQUIRE(seams_cv_rebuild(cv, sys, 0, nullptr, 0, skin) == SEAMS_OK);
  }
  void place(const std::vector<double> &xyz) {
    REQUIRE(seams_system_set_positions(sys, xyz.data(), frame.n, nullptr, 0) == SEAMS_OK);
  }
  double value(int kind, int arg, std::vector<double> *grad = nullptr) {
    double v = 0.0;
    if (grad != nullptr) {
      grad->assign(3 * frame.n, 0.0);
    }
    REQUIRE(seams_cv_value(cv, sys, kind, arg, &v, grad ? grad->data() : nullptr,
                           grad ? grad->size() : 0) == SEAMS_OK);
    return v;
  }
  std::vector<double> states(std::vector<double> *grad = nullptr) {
    std::size_t count = 0;
    REQUIRE(seams_cv_states(cv, sys, nullptr, 0, &count, nullptr, 0) == SEAMS_OK);
    std::vector<double> chi(count, 0.0);
    if (grad != nullptr) {
      grad->assign(count * 3 * frame.n, 0.0);
    }
    REQUIRE(seams_cv_states(cv, sys, chi.data(), chi.size(), &count,
                            grad ? grad->data() : nullptr, grad ? grad->size() : 0) ==
            SEAMS_OK);
    return chi;
  }
};

// Central differences of f over every coordinate, with the frame restored
// after each pair of evaluations.
template <class F> std::vector<double> centralDifference(Api &api, F &&f) {
  std::vector<double> xyz = api.frame.xyz;
  std::vector<double> out(xyz.size(), 0.0);
  for (std::size_t c = 0; c < xyz.size(); ++c) {
    const double x = xyz[c];
    xyz[c] = x + kStep;
    api.place(xyz);
    const double up = f();
    xyz[c] = x - kStep;
    api.place(xyz);
    const double down = f();
    xyz[c] = x;
    out[c] = (up - down) / (2.0 * kStep);
  }
  api.place(xyz);
  return out;
}

void requireClose(const std::vector<double> &analytic, const std::vector<double> &numeric,
                  const std::string &what) {
  REQUIRE(analytic.size() == numeric.size());
  double worst = 0.0;
  double norm = 0.0;
  for (std::size_t c = 0; c < analytic.size(); ++c) {
    const double err = std::abs(analytic[c] - numeric[c]) / (1.0 + std::abs(analytic[c]));
    worst = std::max(worst, err);
    norm = std::max(norm, std::abs(analytic[c]));
  }
  INFO(what << ": worst scaled error " << worst << ", largest entry " << norm);
  REQUIRE(norm > 1e-3);
  REQUIRE(worst < 1e-6);
}

void checkGradient(Api &api, int kind, int arg, const std::string &what) {
  std::vector<double> grad;
  api.value(kind, arg, &grad);
  const auto fd = centralDifference(api, [&] { return api.value(kind, arg); });
  requireClose(grad, fd, what);
}

// Bonds whose raw coordinate difference spans more than half the cell, so
// they are measured through a periodic image.
int bondsAcrossBoundary(const Frame &f, double cutoff) {
  int count = 0;
  for (std::size_t i = 0; i < f.n; ++i) {
    for (std::size_t j = i + 1; j < f.n; ++j) {
      double r2 = 0.0;
      bool wrapped = false;
      for (int d = 0; d < 3; ++d) {
        double x = f.xyz[3 * j + static_cast<std::size_t>(d)] -
                   f.xyz[3 * i + static_cast<std::size_t>(d)];
        const double l = lengthOf(f, d);
        if (std::abs(x) > 0.5 * l) {
          wrapped = true;
          x -= std::round(x / l) * l;
        }
        r2 += x * x;
      }
      count += static_cast<int>(wrapped && r2 < cutoff * cutoff);
    }
  }
  return count;
}

// Switch shared by the soft cases: half weight near 2.85 A, taper from 3.0
// to 3.5 A, so the jittered ice bonds sit on the slope.
constexpr double kR0 = 0.3;
constexpr double kD0 = 2.6;
constexpr int kN = 6;
constexpr int kM = 12;
constexpr double kOn = 3.0;
constexpr double kMax = 3.5;

} // namespace

TEST_CASE("rational switch matches its derivative and tapers to zero", "[smooth_cv]") {
  double s = 0.0;
  double ds = 0.0;
  // x = 1 is the removable singularity of (1 - x^n) / (1 - x^m).
  REQUIRE(seams_switch_rational(kD0 + kR0, kR0, kD0, kN, kM, kOn, kMax, &s, &ds) == SEAMS_OK);
  REQUIRE(std::abs(s - static_cast<double>(kN) / kM) < 1e-14);
  REQUIRE(seams_switch_rational(1.0, kR0, kD0, kN, kM, kOn, kMax, &s, &ds) == SEAMS_OK);
  REQUIRE(s == 1.0);
  REQUIRE(ds == 0.0);
  REQUIRE(seams_switch_rational(kMax, kR0, kD0, kN, kM, kOn, kMax, &s, &ds) == SEAMS_OK);
  REQUIRE(s == 0.0);
  REQUIRE(ds == 0.0);
  for (double r = kD0 + 0.01; r < kMax; r += 0.0173) {
    double up = 0.0;
    double down = 0.0;
    double unused = 0.0;
    REQUIRE(seams_switch_rational(r, kR0, kD0, kN, kM, kOn, kMax, &s, &ds) == SEAMS_OK);
    REQUIRE(seams_switch_rational(r + kStep, kR0, kD0, kN, kM, kOn, kMax, &up, &unused) ==
            SEAMS_OK);
    REQUIRE(seams_switch_rational(r - kStep, kR0, kD0, kN, kM, kOn, kMax, &down, &unused) ==
            SEAMS_OK);
    INFO("r = " << r);
    REQUIRE(std::abs((up - down) / (2.0 * kStep) - ds) < 1e-7);
    REQUIRE(s > 0.0);
    REQUIRE(s <= 1.0);
  }
  // Near d_max the taper is (1 - u)^3 (1 + 3u + 6u^2) <= 10 (1 - u)^3 and
  // its slope 30 u^2 (1 - u)^2 / width, and the rational factor is at most
  // one, so s and ds/dr vanish as the cube and the square of the distance.
  const double width = kMax - kOn;
  for (double delta : {1e-2, 1e-3, 1e-4}) {
    REQUIRE(seams_switch_rational(kMax - delta, kR0, kD0, kN, kM, kOn, kMax, &s, &ds) ==
            SEAMS_OK);
    const double u = delta / width;
    INFO("delta = " << delta);
    REQUIRE(s <= 10.0 * u * u * u);
    REQUIRE(std::abs(ds) <= 31.0 * u * u / width);
  }
  // m = 0 selects 2n, as in PLUMED.
  double s0 = 0.0;
  REQUIRE(seams_switch_rational(2.9, kR0, kD0, kN, 0, kOn, kMax, &s0, nullptr) == SEAMS_OK);
  REQUIRE(seams_switch_rational(2.9, kR0, kD0, kN, 2 * kN, kOn, kMax, &s, nullptr) == SEAMS_OK);
  REQUIRE(s0 == s);
  REQUIRE(seams_switch_rational(2.9, 0.0, kD0, kN, kM, kOn, kMax, &s, nullptr) ==
          SEAMS_ERR_RANGE);
  REQUIRE(seams_switch_rational(2.9, kR0, kD0, kN, kN, kOn, kMax, &s, nullptr) ==
          SEAMS_ERR_RANGE);
  REQUIRE(seams_switch_rational(2.9, kR0, kD0, kN, kM, kMax, kOn, &s, nullptr) ==
          SEAMS_ERR_RANGE);
  REQUIRE(seams_switch_rational(2.9, kR0, kOn + 0.1, kN, kM, kOn, kMax, &s, nullptr) ==
          SEAMS_ERR_RANGE);
  REQUIRE(seams_switch_rational(-1.0, kR0, kD0, kN, kM, kOn, kMax, &s, nullptr) ==
          SEAMS_ERR_RANGE);
}

TEST_CASE("weights of one reproduce the primitive rings and the HC basal test",
          "[smooth_cv]") {
  seams_system *sys = seams_system_create();
  REQUIRE(sys != nullptr);
  REQUIRE(seams_system_read(sys, "traj/ice1h.lammpstrj", 1, -1) == SEAMS_OK);
  seams_neighbours *nb = nullptr;
  REQUIRE(seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) == SEAMS_OK);
  std::size_t nOff = 0;
  std::size_t nnz = 0;
  REQUIRE(seams_neighbours_copy(nb, nullptr, 0, &nOff, nullptr, 0, &nnz) == SEAMS_OK);
  std::vector<int> off(nOff);
  std::vector<int> cols(nnz);
  REQUIRE(seams_neighbours_copy(nb, off.data(), nOff, &nOff, cols.data(), nnz, &nnz) ==
          SEAMS_OK);
  seams_neighbours_destroy(nb);
  const int n = static_cast<int>(nOff) - 1;
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    rows[static_cast<std::size_t>(i)].push_back(i);
    for (int k = off[static_cast<std::size_t>(i)]; k < off[static_cast<std::size_t>(i) + 1];
         ++k) {
      rows[static_cast<std::size_t>(i)].push_back(cols[static_cast<std::size_t>(k)]);
    }
  }

  const smooth::Topology topo = smooth::buildTopology(rows, 6);
  const std::vector<double> ones(topo.edges.size(), 1.0);
  const std::vector<double> w = smooth::ringWeights(topo, ones);
  std::set<std::vector<int>> soft;
  for (int r = 0; r < topo.ringCount(); ++r) {
    const double wr = w[static_cast<std::size_t>(r)];
    REQUIRE((wr == 0.0 || wr == 1.0));
    if (wr == 1.0) {
      std::vector<int> ring(topo.ringAtoms.begin() + topo.ringOff[static_cast<std::size_t>(r)],
                            topo.ringAtoms.begin() +
                                topo.ringOff[static_cast<std::size_t>(r) + 1]);
      std::sort(ring.begin(), ring.end());
      soft.insert(ring);
    }
  }

  seams_rings *rings = nullptr;
  REQUIRE(seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK);
  std::set<std::vector<int>> hard;
  std::vector<std::vector<int>> hexes;
  std::vector<int> atoms(16);
  for (int k = 0; k < seams_rings_count(rings); ++k) {
    std::size_t len = 0;
    REQUIRE(seams_rings_vertices(rings, k, atoms.data(), atoms.size(), &len) == SEAMS_OK);
    std::vector<int> ring(atoms.begin(), atoms.begin() + static_cast<std::ptrdiff_t>(len));
    if (len == 6) {
      hexes.push_back(ring);
    }
    std::sort(ring.begin(), ring.end());
    hard.insert(ring);
  }
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  REQUIRE(hard.size() == 256);
  REQUIRE(soft == hard);

  // The oracle runs ring::basalConditions both ways on every pair of
  // primitive hexagons; the soft count must equal it at weights 0 and 1.
  int oracle = 0;
  for (std::size_t a = 0; a < hexes.size(); ++a) {
    for (std::size_t b = a + 1; b < hexes.size(); ++b) {
      oracle += static_cast<int>(ring::basalConditions(rows, hexes[a], hexes[b]) ||
                                 ring::basalConditions(rows, hexes[b], hexes[a]));
    }
  }
  REQUIRE(oracle > 0);
  const double pairs = smooth::hexBasalPairs(topo, ones, w, nullptr, nullptr);
  REQUIRE(pairs == static_cast<double>(oracle));
}

TEST_CASE("ring and hexagon-pair gradients match central differences", "[smooth_cv]") {
  Frame ice = jitteredIce(3, 2, 2, 0.15, 7);
  // Shift the crystal so the cell faces cut through rings.
  for (std::size_t i = 0; i < ice.n; ++i) {
    ice.xyz[3 * i] += 1.3;
    ice.xyz[3 * i + 1] += 0.7;
  }
  wrap(ice);
  REQUIRE(bondsAcrossBoundary(ice, kMax) > 10);
  Api api(ice, kR0, kD0, kN, kM, kOn, kMax, 6);
  api.rebuild(0.4);
  REQUIRE(api.value(SEAMS_CV_RING_COUNT, 6) > 1.0);
  checkGradient(api, SEAMS_CV_RING_COUNT, 0, "ring count, every size");
  checkGradient(api, SEAMS_CV_RING_COUNT, 6, "ring count, six");
  checkGradient(api, SEAMS_CV_RING_ATOMS, 6, "ring atoms, six");
  checkGradient(api, SEAMS_CV_HEX_BASAL_PAIRS, 0, "hexagon basal pairs");
}

TEST_CASE("gyration gradients match central differences across the boundary", "[smooth_cv]") {
  const Frame cluster = iceCluster(6.5, 30.0, {0.5, 29.6, 15.0}, 0.15, 11);
  REQUIRE(cluster.n > 30);
  REQUIRE(bondsAcrossBoundary(cluster, kMax) > 5);
  Api api(cluster, kR0, kD0, kN, kM, kOn, kMax, 6);
  api.rebuild(0.4);
  const double rg2 = api.value(SEAMS_CV_GYRATION_RG2, -1);
  // The unwrapped cluster has the radius of gyration of a ball of radius
  // about 6.5 A, not one stretched over the 30 A cell.
  REQUIRE(rg2 > 9.0);
  REQUIRE(rg2 < 30.0);
  REQUIRE(api.value(SEAMS_CV_GYRATION_VOLUME, 6) > 0.0);
  checkGradient(api, SEAMS_CV_GYRATION_RG2, -1, "rg2, uniform");
  checkGradient(api, SEAMS_CV_GYRATION_RG2, 6, "rg2, six-ring weights");
  checkGradient(api, SEAMS_CV_GYRATION_VOLUME, -1, "volume, uniform");
  checkGradient(api, SEAMS_CV_GYRATION_VOLUME, 0, "volume, ring weights");
}

TEST_CASE("translation through the boundary leaves values and gradients alone",
          "[smooth_cv]") {
  const Frame cluster = iceCluster(6.5, 30.0, {15.0, 15.0, 15.0}, 0.15, 11);
  Frame moved = cluster;
  for (std::size_t i = 0; i < moved.n; ++i) {
    moved.xyz[3 * i] += 14.6;
    moved.xyz[3 * i + 2] -= 14.8;
  }
  wrap(moved);
  REQUIRE(bondsAcrossBoundary(moved, kMax) > bondsAcrossBoundary(cluster, kMax));
  Api a(cluster, kR0, kD0, kN, kM, kOn, kMax, 6);
  Api b(moved, kR0, kD0, kN, kM, kOn, kMax, 6);
  a.rebuild(0.4);
  b.rebuild(0.4);
  const std::array<std::pair<int, int>, 6> cvs{{{SEAMS_CV_RING_COUNT, 6},
                                                {SEAMS_CV_RING_ATOMS, 0},
                                                {SEAMS_CV_HEX_BASAL_PAIRS, 0},
                                                {SEAMS_CV_GYRATION_RG2, -1},
                                                {SEAMS_CV_GYRATION_VOLUME, 6},
                                                {SEAMS_CV_GYRATION_RG2, 0}}};
  for (const auto &[kind, arg] : cvs) {
    std::vector<double> ga;
    std::vector<double> gb;
    const double va = a.value(kind, arg, &ga);
    const double vb = b.value(kind, arg, &gb);
    INFO("kind " << kind << " arg " << arg);
    REQUIRE(std::abs(va - vb) <= 1e-9 * (1.0 + std::abs(va)));
    for (std::size_t c = 0; c < ga.size(); ++c) {
      REQUIRE(std::abs(ga[c] - gb[c]) <= 1e-9 * (1.0 + std::abs(ga[c])));
    }
  }
}

TEST_CASE("pairs beyond d_max in the list change nothing", "[smooth_cv]") {
  // Skin 0 lists pairs to d_max; skin 0.7 lists the near edge of the next
  // shell as well, which adds candidate rings, shortcut paths and hexagon
  // pairs that hold a zero bond weight. The polynomial is the same, so the
  // values and gradients agree to rounding.
  const Frame ice = jitteredIce(3, 2, 2, 0.15, 5);
  Api narrow(ice, kR0, kD0, kN, kM, kOn, kMax, 7);
  Api wide(ice, kR0, kD0, kN, kM, kOn, kMax, 7);
  narrow.rebuild(0.0);
  wide.rebuild(0.7);
  std::size_t edgesN = 0;
  std::size_t ringsN = 0;
  std::size_t edgesW = 0;
  std::size_t ringsW = 0;
  REQUIRE(seams_cv_topology(narrow.cv, &edgesN, &ringsN, nullptr) == SEAMS_OK);
  REQUIRE(seams_cv_topology(wide.cv, &edgesW, &ringsW, nullptr) == SEAMS_OK);
  REQUIRE(edgesW > edgesN);
  REQUIRE(ringsW > ringsN);
  for (const int kind : {SEAMS_CV_RING_COUNT, SEAMS_CV_RING_ATOMS}) {
    for (const int arg : {0, 3, 4, 5, 6, 7}) {
      std::vector<double> gn;
      std::vector<double> gw;
      const double vn = narrow.value(kind, arg, &gn);
      const double vw = wide.value(kind, arg, &gw);
      INFO("kind " << kind << " arg " << arg);
      REQUIRE(std::abs(vn - vw) <= 1e-12 * (1.0 + std::abs(vn)));
      for (std::size_t c = 0; c < gn.size(); ++c) {
        REQUIRE(std::abs(gn[c] - gw[c]) <= 1e-12 * (1.0 + std::abs(gn[c])));
      }
    }
  }
  std::vector<double> gn;
  std::vector<double> gw;
  const double vn = narrow.value(SEAMS_CV_HEX_BASAL_PAIRS, 0, &gn);
  const double vw = wide.value(SEAMS_CV_HEX_BASAL_PAIRS, 0, &gw);
  REQUIRE(vn > 0.0);
  REQUIRE(std::abs(vn - vw) <= 1e-12 * (1.0 + vn));
  for (std::size_t c = 0; c < gn.size(); ++c) {
    REQUIRE(std::abs(gn[c] - gw[c]) <= 1e-12 * (1.0 + std::abs(gn[c])));
  }
}

TEST_CASE("a chord crossing d_max turns a square into two triangles smoothly",
          "[smooth_cv]") {
  // A rhombus of side 2.8 A whose short diagonal sweeps from 3.8 A to
  // 2.9 A, through d_max = 3.5 A. Inside d_max the chord is a shortcut of
  // the square and an edge of two triangles. A count that kept only the
  // hard primitive rings would drop the square by its full weight, about
  // 0.7, where the chord enters; the soft count follows one polynomial on
  // both sides and must agree with a topology that holds the chord
  // throughout.
  const double side = 2.8;
  const double boxLen = 20.0;
  auto rhombus = [&](double chord) {
    Frame f;
    const double half = 0.5 * chord;
    const double h = std::sqrt(side * side - half * half);
    // Centred on x = 0, so the square straddles the periodic boundary.
    f.xyz = {-half, 10.0, 10.0, 0.0, 10.0 - h, 10.0, half, 10.0, 10.0, 0.0, 10.0 + h, 10.0};
    f.n = 4;
    f.box9[3] = f.box9[4] = f.box9[5] = boxLen;
    wrap(f);
    return f;
  };
  Api ref(rhombus(2.9), kR0, kD0, kN, kM, kOn, kMax, 6);
  ref.rebuild(0.0);
  std::size_t edges = 0;
  REQUIRE(seams_cv_topology(ref.cv, &edges, nullptr, nullptr) == SEAMS_OK);
  REQUIRE(edges == 5);
  // Every atom stays within 0.25 A of its place at a chord of 3.35 A, so
  // skin 0.5 keeps this topology valid over the sweep, and the long
  // diagonal, at least 4.1 A, stays outside the 4.0 A list.
  Api fixed(rhombus(3.35), kR0, kD0, kN, kM, kOn, kMax, 6);
  fixed.rebuild(0.5);
  REQUIRE(seams_cv_topology(fixed.cv, &edges, nullptr, nullptr) == SEAMS_OK);
  REQUIRE(edges == 5);
  double before = 0.0;
  double slopeBefore = 0.0;
  double far = 0.0;
  bool first = true;
  const double step = 2e-3;
  for (int k = 0; k <= 450; ++k) {
    const double chord = 3.8 - step * k;
    const Frame f = rhombus(chord);
    Api fresh(f, kR0, kD0, kN, kM, kOn, kMax, 6);
    fresh.rebuild(0.0);
    fixed.place(f.xyz);
    std::vector<double> grad;
    const double total = fresh.value(SEAMS_CV_RING_COUNT, 0, &grad);
    const double four = fresh.value(SEAMS_CV_RING_COUNT, 4);
    const double three = fresh.value(SEAMS_CV_RING_COUNT, 3);
    INFO("chord " << chord);
    REQUIRE(std::abs(total - four - three) < 1e-14);
    REQUIRE(std::abs(fixed.value(SEAMS_CV_RING_COUNT, 0) - total) < 1e-14);
    // d total / d chord: atoms 0 and 2 move apart along x and atoms 1 and 3
    // along y, keeping the side length.
    const double h = std::sqrt(side * side - 0.25 * chord * chord);
    const double dh = -0.25 * chord / h;
    const double slope = -0.5 * grad[0] + 0.5 * grad[6] - dh * grad[4] + dh * grad[10];
    if (first) {
      far = four;
      REQUIRE(three == 0.0);
      REQUIRE(far > 0.6);
    } else {
      // With the sides fixed the count is a + b s_c, s_c the chord weight
      // and b = 2 s^2 - s^4 < 1, so it moves by at most max |s'| step and
      // its slope by max |s''| step. Over 2.9 .. 3.8 A the switch has
      // |s'| <= 5 and |s''| <= 36 per A, far below the 0.7 per step a
      // dropped square costs.
      REQUIRE(std::abs(total - before) < 5.0 * step);
      REQUIRE(std::abs(slope - slopeBefore) < 36.0 * step);
    }
    before = total;
    slopeBefore = slope;
    first = false;
    if (k == 450) {
      // At 2.9 A the chord has weight one half: the square keeps half its
      // weight and each triangle has the chord's weight times two sides.
      REQUIRE(std::abs(four - 0.5 * far) < 1e-12);
      REQUIRE(three > 0.5);
    }
  }
}

TEST_CASE("soft states compose and their gradients match central differences",
          "[smooth_cv]") {
  const Frame cluster = iceCluster(6.5, 30.0, {29.8, 15.0, 0.4}, 0.15, 13);
  Api api(cluster, kR0, kD0, kN, kM, kOn, kMax, 6);
  api.rebuild(0.4);
  const double six = api.value(SEAMS_CV_RING_COUNT, 6);
  const double rg2 = api.value(SEAMS_CV_GYRATION_RG2, -1);
  const double pairs = api.value(SEAMS_CV_HEX_BASAL_PAIRS, 0);
  // State 0: many six-rings and a compact cluster. State 1: not state 0.
  // State 2: state 0 and some basal pairs.
  REQUIRE(seams_cv_state_clause(api.cv, 0, SEAMS_CV_RING_COUNT, 6, six - 0.5, 1.0, 1) ==
          SEAMS_OK);
  REQUIRE(seams_cv_state_clause(api.cv, 0, SEAMS_CV_GYRATION_RG2, -1, rg2 + 0.3, 0.5, 0) ==
          SEAMS_OK);
  REQUIRE(seams_cv_state_clause(api.cv, 1, SEAMS_CV_STATE, 0, 0.0, 0.0, 0) == SEAMS_OK);
  REQUIRE(seams_cv_state_clause(api.cv, 2, SEAMS_CV_STATE, 0, 0.0, 0.0, 1) == SEAMS_OK);
  REQUIRE(seams_cv_state_clause(api.cv, 2, SEAMS_CV_HEX_BASAL_PAIRS, 0, pairs - 0.2, 0.4, 1) ==
          SEAMS_OK);
  REQUIRE(seams_cv_state_clause(api.cv, 4, SEAMS_CV_STATE, 0, 0.0, 0.0, 1) ==
          SEAMS_ERR_RANGE);
  REQUIRE(seams_cv_state_clause(api.cv, 3, SEAMS_CV_STATE, 3, 0.0, 0.0, 1) ==
          SEAMS_ERR_RANGE);
  REQUIRE(seams_cv_state_clause(api.cv, 3, SEAMS_CV_RING_COUNT, 6, 1.0, 0.0, 1) ==
          SEAMS_ERR_RANGE);

  std::vector<double> grad;
  const std::vector<double> chi = api.states(&grad);
  REQUIRE(chi.size() == 3);
  const double logistic0 = 1.0 / (1.0 + std::exp(-0.5));
  const double logistic1 = 1.0 - 1.0 / (1.0 + std::exp(-(-0.3 / 0.5)));
  REQUIRE(std::abs(chi[0] - logistic0 * logistic1) < 1e-12);
  REQUIRE(std::abs(chi[0] + chi[1] - 1.0) < 1e-14);
  REQUIRE(std::abs(chi[2] - chi[0] / (1.0 + std::exp(-0.5))) < 1e-12);
  const std::size_t dim = 3 * cluster.n;
  for (std::size_t k = 0; k < chi.size(); ++k) {
    std::vector<double> row(grad.begin() + static_cast<std::ptrdiff_t>(k * dim),
                            grad.begin() + static_cast<std::ptrdiff_t>((k + 1) * dim));
    std::vector<double> single;
    REQUIRE(std::abs(api.value(SEAMS_CV_STATE, static_cast<int>(k), &single) - chi[k]) <
            1e-14);
    for (std::size_t c = 0; c < dim; ++c) {
      REQUIRE(std::abs(single[c] - row[c]) < 1e-14);
    }
    const auto fd = centralDifference(api, [&] { return api.states()[k]; });
    requireClose(row, fd, "state " + std::to_string(k));
  }
}

TEST_CASE("the list goes stale after half the skin and selections are honoured",
          "[smooth_cv]") {
  Frame ice = jitteredIce(3, 2, 2, 0.1, 3);
  Api api(ice, kR0, kD0, kN, kM, kOn, kMax, 6);
  double v = 0.0;
  REQUIRE(seams_cv_value(api.cv, api.sys, SEAMS_CV_RING_COUNT, 6, &v, nullptr, 0) ==
          SEAMS_ERR_STATE);
  api.rebuild(0.4);
  int stale = -1;
  REQUIRE(seams_cv_needs_rebuild(api.cv, api.sys, &stale) == SEAMS_OK);
  REQUIRE(stale == 0);
  std::vector<double> xyz = ice.xyz;
  xyz[3 * 5] += 0.19;
  api.place(xyz);
  REQUIRE(seams_cv_needs_rebuild(api.cv, api.sys, &stale) == SEAMS_OK);
  REQUIRE(stale == 0);
  xyz[3 * 5] += 0.02;
  api.place(xyz);
  REQUIRE(seams_cv_needs_rebuild(api.cv, api.sys, &stale) == SEAMS_OK);
  REQUIRE(stale == 1);
  REQUIRE(seams_cv_value(api.cv, api.sys, SEAMS_CV_RING_COUNT, 6, &v, nullptr, 0) ==
          SEAMS_ERR_STATE);
  api.rebuild(0.4);
  REQUIRE(seams_cv_value(api.cv, api.sys, SEAMS_CV_RING_COUNT, 6, &v, nullptr, 0) == SEAMS_OK);

  // A selection drops the other atoms from the graph and their gradient.
  std::vector<int> half;
  for (int i = 0; i < static_cast<int>(ice.n) / 2; ++i) {
    half.push_back(i);
  }
  REQUIRE(seams_cv_rebuild(api.cv, api.sys, 0, half.data(), half.size(), 0.4) == SEAMS_OK);
  std::vector<double> grad;
  api.value(SEAMS_CV_RING_ATOMS, 0, &grad);
  for (std::size_t i = half.size(); i < ice.n; ++i) {
    for (std::size_t d = 0; d < 3; ++d) {
      REQUIRE(grad[3 * i + d] == 0.0);
    }
  }
  const std::vector<int> twice{1, 2, 1};
  REQUIRE(seams_cv_rebuild(api.cv, api.sys, 0, twice.data(), twice.size(), 0.4) ==
          SEAMS_ERR_ARG);
  const std::vector<int> outside{1, static_cast<int>(ice.n)};
  REQUIRE(seams_cv_rebuild(api.cv, api.sys, 0, outside.data(), outside.size(), 0.4) ==
          SEAMS_ERR_RANGE);
  REQUIRE(seams_cv_rebuild(api.cv, api.sys, 0, nullptr, 0, 5.0) == SEAMS_ERR_BOX);
}
