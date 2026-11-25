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
#include <limits>
#include <random>
#include <set>
#include <string>
#include <utility>
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

// The ice atoms within radius of the cell centre, centred on `at`, not
// wrapped and with no cell.
Frame iceBall(double radius, std::array<double, 3> at, double amp, unsigned seed) {
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
  return f;
}

// iceBall in a cubic box of side box, so that the cluster straddles a
// periodic boundary.
Frame iceCluster(double radius, double box, std::array<double, 3> at, double amp,
                 unsigned seed) {
  Frame f = iceBall(radius, at, amp, seed);
  f.box9[3] = f.box9[4] = f.box9[5] = box;
  wrap(f);
  return f;
}

// The tilted cell of the triclinic BOP test: edges 7 A, every tilt 3.4 A.
// Its faces are 6.14 A apart, so pairs from 3.07 A to 3.5 A have a shorter
// Euclidean image than the wrapped fractional difference.
constexpr double kTilted[9] = {0.0, 0.0, 0.0, 7.0, 7.0, 7.0, 3.4, 3.4, 3.4};

// Wraps the fractional coordinates of every atom of a cell given as
// ox, oy, oz, lx, ly, lz, xy, xz, yz into [0, 1).
void wrapTilted(Frame &f) {
  const double lx = f.box9[3];
  const double ly = f.box9[4];
  const double lz = f.box9[5];
  const double xy = f.box9[6];
  const double xz = f.box9[7];
  const double yz = f.box9[8];
  for (std::size_t i = 0; i < f.n; ++i) {
    double *p = f.xyz.data() + 3 * i;
    double sz = (p[2] - f.box9[2]) / lz;
    double sy = (p[1] - f.box9[1] - yz * sz) / ly;
    double sx = (p[0] - f.box9[0] - xy * sy - xz * sz) / lx;
    sx -= std::floor(sx);
    sy -= std::floor(sy);
    sz -= std::floor(sz);
    p[0] = f.box9[0] + lx * sx + xy * sy + xz * sz;
    p[1] = f.box9[1] + ly * sy + yz * sz;
    p[2] = f.box9[2] + lz * sz;
  }
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

TEST_CASE("the switch slope at d0 is the right derivative", "[smooth_cv]") {
  // Past d0, R = 1 - x^n + O(x^(n+1)) with m = 2n, so the right derivative
  // at d0 is -1 / r0 for n = 1 and 0 for n >= 2, and the forward difference
  // is within h / r0^2 = 1.2e-4 / A of it. At d0 = 0 that is the only side.
  for (const double d0 : {0.0, kD0}) {
    for (const int n : {1, 2, 3}) {
      double s = 0.0;
      double ds = 0.0;
      double up = 0.0;
      double unused = 0.0;
      REQUIRE(seams_switch_rational(d0, kR0, d0, n, 0, kOn, kMax, &s, &ds) == SEAMS_OK);
      REQUIRE(seams_switch_rational(d0 + kStep, kR0, d0, n, 0, kOn, kMax, &up, &unused) ==
              SEAMS_OK);
      const double right = n == 1 ? -1.0 / kR0 : 0.0;
      INFO("d0 = " << d0 << ", n = " << n << ", ds/dr = " << ds);
      REQUIRE(s == 1.0);
      REQUIRE(std::abs(ds - right) <= 1e-12 / kR0);
      REQUIRE(std::abs((up - s) / kStep - right) < 2e-4);
    }
  }
}

TEST_CASE("large switch exponents stay finite past x = 1", "[smooth_cv]") {
  // x^(m-1) overflows for these exponents. With m = 2n, R = 1 / (1 + x^n),
  // within a relative 1e-300 of y^n for y = 1 / x, and ds/dr is
  // -n R / (r - d0) T + R T'.
  double s = 0.0;
  double ds = 0.0;
  const double r = kD0 + 2.0 * kR0;
  REQUIRE(seams_switch_rational(r, kR0, kD0, 1000, 0, kOn, kMax, &s, &ds) == SEAMS_OK);
  const double y = kR0 / (r - kD0);
  const double u = (r - kOn) / (kMax - kOn);
  const double taper = 1.0 - u * u * u * (10.0 - 15.0 * u + 6.0 * u * u);
  const double dTaper = -30.0 * u * u * (1.0 - u) * (1.0 - u) / (kMax - kOn);
  const double rational = std::pow(y, 1000);
  INFO("s = " << s << ", ds/dr = " << ds);
  REQUIRE(rational > 1e-305);
  REQUIRE(std::abs(s - rational * taper) < 1e-10 * rational * taper);
  const double slope = -1000.0 * rational / (r - kD0) * taper + rational * dTaper;
  REQUIRE(std::abs(ds - slope) < 1e-10 * std::abs(slope));
  double up = 0.0;
  double down = 0.0;
  REQUIRE(seams_switch_rational(r + kStep, kR0, kD0, 1000, 0, kOn, kMax, &up, nullptr) ==
          SEAMS_OK);
  REQUIRE(seams_switch_rational(r - kStep, kR0, kD0, 1000, 0, kOn, kMax, &down, nullptr) ==
          SEAMS_OK);
  REQUIRE(std::abs((up - down) / (2.0 * kStep) - ds) < 1e-4 * std::abs(ds));
  // x = 40 with n = 200: R is about 1e-320, below the normal doubles.
  REQUIRE(seams_switch_rational(3.0, 0.01, kD0, 200, 400, kOn, kMax, &s, &ds) == SEAMS_OK);
  REQUIRE(std::isfinite(s));
  REQUIRE(std::isfinite(ds));
  REQUIRE(s >= 0.0);
  REQUIRE(s < 1e-300);
}

TEST_CASE("every finite switch argument gives a finite result or SEAMS_ERR_RANGE",
          "[smooth_cv]") {
  constexpr double kTiny = 1e-310;
  constexpr int kIntMax = std::numeric_limits<int>::max();
  double s = -1.0;
  double ds = -1.0;
  // A slope of order 1 / r0 or 1 / (d_max - d_on) beyond the largest double
  // is an error and leaves both outputs alone; s alone is still available.
  REQUIRE(seams_switch_rational(0.5 * kTiny, kTiny, 0.0, kN, kM, kOn, kMax, &s, &ds) ==
          SEAMS_ERR_RANGE);
  REQUIRE(s == -1.0);
  REQUIRE(ds == -1.0);
  REQUIRE(seams_switch_rational(0.5 * kTiny, kTiny, 0.0, kN, kM, kOn, kMax, &s, nullptr) ==
          SEAMS_OK);
  REQUIRE(s > 0.9);
  REQUIRE(s <= 1.0);
  s = -1.0;
  REQUIRE(seams_switch_rational(0.5 * kTiny, kR0, 0.0, kN, kM, 0.0, kTiny, &s, &ds) ==
          SEAMS_ERR_RANGE);
  REQUIRE(s == -1.0);
  REQUIRE(seams_switch_rational(0.5 * kTiny, kR0, 0.0, kN, kM, 0.0, kTiny, &s, nullptr) ==
          SEAMS_OK);
  REQUIRE(std::isfinite(s));
  // m = 0 doubles n, which must fit in an int. Exponents above 4096 are
  // rejected even when 2n still fits.
  REQUIRE(seams_switch_rational(1.0, kR0, kD0, kIntMax, 0, kOn, kMax, &s, nullptr) ==
          SEAMS_ERR_RANGE);
  REQUIRE(seams_switch_rational(1.0, kR0, kD0, kIntMax / 2 + 1, 0, kOn, kMax, &s, nullptr) ==
          SEAMS_ERR_RANGE);
  s = -1.0;
  REQUIRE(seams_switch_rational(1.0, kR0, kD0, kIntMax / 2, 0, kOn, kMax, &s, nullptr) ==
          SEAMS_ERR_RANGE);
  REQUIRE(s == -1.0);
  REQUIRE(std::string(seams_last_error_message()).find("4096") != std::string::npos);

  const double huge = std::numeric_limits<double>::max();
  const std::array<std::pair<int, int>, 6> exponents{
      {{1, 0}, {2, 0}, {6, 12}, {50, 0}, {1000, 0}, {999, 1000}}};
  std::size_t checked = 0;
  std::size_t refused = 0;
  for (const double r0 : {kTiny, 1e-200, 1e-6, kR0, 1e6, 1e300, huge}) {
    for (const double d0 : {0.0, kD0}) {
      const std::array<std::array<double, 2>, 3> tapers{
          {{d0, d0 + 1e-12}, {kOn, kMax}, {1e300, huge}}};
      for (const auto &[n, m] : exponents) {
        for (const auto &[on, max] : tapers) {
          for (const double r : {0.0, d0, d0 + r0, d0 + 0.5 * r0, 0.5 * (on + max), 2.99, 3.2,
                                 std::nextafter(max, 0.0), 1e299, 1e301}) {
            if (!std::isfinite(r)) {
              continue;
            }
            s = -1.0;
            ds = -1.0;
            const seams_status status = seams_switch_rational(r, r0, d0, n, m, on, max, &s, &ds);
            INFO("r = " << r << ", r0 = " << r0 << ", d0 = " << d0 << ", n = " << n
                        << ", m = " << m << ", d_on = " << on << ", d_max = " << max);
            REQUIRE((status == SEAMS_OK || status == SEAMS_ERR_RANGE));
            if (status == SEAMS_OK) {
              REQUIRE(std::isfinite(s));
              REQUIRE(std::isfinite(ds));
            } else {
              REQUIRE(s == -1.0);
              REQUIRE(ds == -1.0);
              ++refused;
            }
            REQUIRE(seams_switch_rational(r, r0, d0, n, m, on, max, &s, nullptr) == SEAMS_OK);
            REQUIRE(std::isfinite(s));
            REQUIRE(s >= 0.0);
            REQUIRE(s <= 1.0 + 1e-15);
            ++checked;
          }
        }
      }
    }
  }
  INFO(checked << " arguments, " << refused << " slopes out of range");
  REQUIRE(checked > 2000);
  REQUIRE(refused > 0);
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

TEST_CASE("the gyration volume has a kink where the cluster becomes flat", "[smooth_cv]") {
  // A planar hexagon has det G = 0, so V = 0 with a zero gradient. Lifting
  // one atom by h leaves the in-plane block of G alone and makes det G
  // proportional to h^2, so V = c |h|, with slopes of opposite sign on the
  // two sides of the plane.
  const double pi = std::acos(-1.0);
  Frame hex;
  hex.n = 6;
  for (int k = 0; k < 6; ++k) {
    const double t = k * pi / 3.0;
    hex.xyz.insert(hex.xyz.end(),
                   {10.0 + 2.75 * std::cos(t), 10.0 + 2.75 * std::sin(t), 10.0});
  }
  hex.box9[3] = hex.box9[4] = hex.box9[5] = 20.0;
  Api api(hex, kR0, kD0, kN, kM, kOn, kMax, 6);
  api.rebuild(0.3);
  std::vector<double> grad;
  REQUIRE(api.value(SEAMS_CV_GYRATION_VOLUME, -1, &grad) == 0.0);
  REQUIRE(std::all_of(grad.begin(), grad.end(), [](double g) { return g == 0.0; }));
  REQUIRE(api.value(SEAMS_CV_GYRATION_RG2, -1) > 0.0);
  double slope = 0.0;
  for (const double h : {1e-2, 1e-3, -1e-3, -1e-2}) {
    std::vector<double> xyz = hex.xyz;
    xyz[2] += h;
    api.place(xyz);
    const double v = api.value(SEAMS_CV_GYRATION_VOLUME, -1, &grad);
    if (slope == 0.0) {
      slope = v / std::abs(h);
    }
    INFO("h = " << h << ", V = " << v << ", dV/dz0 = " << grad[2]);
    REQUIRE(slope > 1.0);
    REQUIRE(std::abs(v / std::abs(h) - slope) < 1e-9 * slope);
    REQUIRE(std::abs(grad[2] - std::copysign(slope, h)) < 1e-6 * slope);
  }
}

TEST_CASE("ring-weighted gyration drops to 0 when the last ring breaks", "[smooth_cv]") {
  // The ring-weighted values are ratios over the total participation W, so
  // an equilateral triangle keeps Rg^2 = side^2 / 3 while W fades and has
  // Rg^2 = 0 once its bonds pass d_max. With r0 = 1 A the ring weight at a
  // side of 3.49 A is still about 1e-13. Unit weights stay continuous.
  const double pi = std::acos(-1.0);
  auto triangle = [pi](double side) {
    std::vector<double> xyz;
    const double radius = side / std::sqrt(3.0);
    for (int k = 0; k < 3; ++k) {
      const double t = 2.0 * k * pi / 3.0;
      xyz.insert(xyz.end(), {10.0 + radius * std::cos(t), 10.0 + radius * std::sin(t), 10.0});
    }
    return xyz;
  };
  Frame f;
  f.n = 3;
  f.xyz = triangle(3.4);
  f.box9[3] = f.box9[4] = f.box9[5] = 20.0;
  Api api(f, 1.0, kD0, kN, kM, kOn, kMax, 3);
  api.rebuild(0.3);
  for (const double side : {3.4, 3.49, 3.51}) {
    api.place(triangle(side));
    std::vector<double> grad;
    const double ring = api.value(SEAMS_CV_RING_COUNT, 3);
    const double rg2 = api.value(SEAMS_CV_GYRATION_RG2, 3, &grad);
    const double plain = api.value(SEAMS_CV_GYRATION_RG2, -1);
    INFO("side " << side << ": ring weight " << ring << ", ring Rg2 " << rg2);
    REQUIRE(std::abs(plain - side * side / 3.0) < 1e-12);
    if (side < kMax) {
      REQUIRE(ring > 0.0);
      REQUIRE(std::abs(rg2 - side * side / 3.0) < 1e-9);
    } else {
      REQUIRE(ring == 0.0);
      REQUIRE(rg2 == 0.0);
      REQUIRE(std::all_of(grad.begin(), grad.end(), [](double g) { return g == 0.0; }));
    }
  }
}

TEST_CASE("triclinic gyration places an unbonded part by its Euclidean minimum image",
          "[smooth_cv][triclinic]") {
  // The pair of the triclinic BOP test is 3.330796 A apart through its
  // nearest image, while the wrapped fractional difference is 5.399463 A.
  // No bond joins them, so the second atom is placed from the first.
  Frame f;
  f.xyz = {0.0, 0.0, 0.0, 4.79, 2.26, 1.05};
  f.n = 2;
  std::copy(kTilted, kTilted + 9, f.box9);
  Api api(f, 0.2, 2.5, 6, 0, 2.7, 2.9, 6);
  api.rebuild(0.16);
  std::size_t edges = 1;
  REQUIRE(seams_cv_topology(api.cv, &edges, nullptr, nullptr) == SEAMS_OK);
  REQUIRE(edges == 0);
  const std::array<double, 3> d{4.79 - 7.0, 2.26, 1.05};
  const double d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  REQUIRE(std::abs(std::sqrt(d2) - 3.330796) < 1e-6);
  std::vector<double> grad;
  const double rg2 = api.value(SEAMS_CV_GYRATION_RG2, -1, &grad);
  INFO("rg2 " << rg2);
  REQUIRE(std::abs(rg2 - 2.773550) < 1e-12);
  // Two atoms of weight one: Rg^2 = |d|^2 / 4 and its gradient is -/+ d / 2.
  for (std::size_t k = 0; k < 3; ++k) {
    REQUIRE(std::abs(grad[k] + 0.5 * d[k]) < 1e-12);
    REQUIRE(std::abs(grad[3 + k] - 0.5 * d[k]) < 1e-12);
  }
  checkGradient(api, SEAMS_CV_GYRATION_RG2, -1, "rg2, unbonded pair in the tilted cell");
}

TEST_CASE("triclinic gyration unwraps a bond stretched past the Smith ball",
          "[smooth_cv][triclinic]") {
  // A pair 3.00 A apart along the normal of the narrowest faces is listed
  // with d_max + skin = 3.06 A. Moving each atom 0.07 A, inside half the
  // skin, stretches it to 3.14 A, past half the face separation (3.07 A),
  // where the wrapped fractional difference is 4.50 A.
  // The first row of the inverse cell matrix, b x c up to scale.
  std::array<double, 3> normal{kTilted[4] * kTilted[5], -kTilted[6] * kTilted[5],
                               kTilted[6] * kTilted[8] - kTilted[4] * kTilted[7]};
  const double norm = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                                normal[2] * normal[2]);
  for (double &x : normal) {
    x /= norm;
  }
  const std::array<double, 3> mid{6.9, 5.2, 3.5};
  auto pair = [&](double dist) {
    Frame f;
    f.n = 2;
    std::copy(kTilted, kTilted + 9, f.box9);
    for (const double sign : {-0.5, 0.5}) {
      for (std::size_t k = 0; k < 3; ++k) {
        f.xyz.push_back(mid[k] + sign * dist * normal[k]);
      }
    }
    wrapTilted(f);
    return f;
  };
  Api api(pair(3.0), 0.2, 2.5, 6, 0, 2.7, 2.9, 6);
  api.rebuild(0.16);
  std::size_t edges = 0;
  REQUIRE(seams_cv_topology(api.cv, &edges, nullptr, nullptr) == SEAMS_OK);
  REQUIRE(edges == 1);
  REQUIRE(std::abs(api.value(SEAMS_CV_GYRATION_RG2, -1) - 2.25) < 1e-12);
  const Frame stretched = pair(3.14);
  api.place(stretched.xyz);
  int stale = 1;
  REQUIRE(seams_cv_needs_rebuild(api.cv, api.sys, &stale) == SEAMS_OK);
  REQUIRE(stale == 0);
  api.frame = stretched;
  const double rg2 = api.value(SEAMS_CV_GYRATION_RG2, -1);
  INFO("rg2 " << rg2);
  REQUIRE(std::abs(rg2 - 0.25 * 3.14 * 3.14) < 1e-12);
  checkGradient(api, SEAMS_CV_GYRATION_RG2, -1, "rg2, stretched bond in the tilted cell");
}

TEST_CASE("a cluster in a triclinic cell has the values and gradients of the cubic cell",
          "[smooth_cv][triclinic]") {
  // The same ice cluster, cut through the faces of a cubic and of a tilted
  // 30 A cell. Its bonds and atoms are the same minimum-image geometry in
  // both, so every variable and gradient agrees to rounding, and the tilted
  // gradients match central differences.
  // Centred on the face centre at fractional (0, 1/2, 1/2) of the tilted cell.
  const double cell[9] = {0.0, 0.0, 0.0, 30.0, 30.0, 30.0, 12.0, -9.0, 10.5};
  const std::array<double, 3> face{0.5 * (cell[6] + cell[7]), 0.5 * (cell[4] + cell[8]),
                                   0.5 * cell[5]};
  const Frame ball = iceBall(6.5, face, 0.15, 11);
  const Frame cubic = iceCluster(6.5, 30.0, face, 0.15, 11);
  Frame tilted = ball;
  std::copy(cell, cell + 9, tilted.box9);
  wrapTilted(tilted);
  std::size_t moved = 0;
  for (std::size_t i = 0; i < ball.n; ++i) {
    moved += static_cast<std::size_t>(tilted.xyz[3 * i] != ball.xyz[3 * i] ||
                                      tilted.xyz[3 * i + 1] != ball.xyz[3 * i + 1] ||
                                      tilted.xyz[3 * i + 2] != ball.xyz[3 * i + 2]);
  }
  INFO("atoms wrapped through a tilted face: " << moved << " of " << ball.n);
  REQUIRE(moved > 5);
  REQUIRE(moved + 5 < ball.n);
  Api a(cubic, kR0, kD0, kN, kM, kOn, kMax, 6);
  Api b(tilted, kR0, kD0, kN, kM, kOn, kMax, 6);
  a.rebuild(0.4);
  b.rebuild(0.4);
  const std::array<std::pair<int, int>, 7> cvs{{{SEAMS_CV_RING_COUNT, 6},
                                                {SEAMS_CV_RING_ATOMS, 0},
                                                {SEAMS_CV_HEX_BASAL_PAIRS, 0},
                                                {SEAMS_CV_GYRATION_RG2, -1},
                                                {SEAMS_CV_GYRATION_RG2, 6},
                                                {SEAMS_CV_GYRATION_VOLUME, -1},
                                                {SEAMS_CV_GYRATION_VOLUME, 0}}};
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
  REQUIRE(b.value(SEAMS_CV_RING_COUNT, 6) > 1.0);
  checkGradient(b, SEAMS_CV_RING_COUNT, 6, "ring count, tilted cell");
  checkGradient(b, SEAMS_CV_HEX_BASAL_PAIRS, 0, "hexagon basal pairs, tilted cell");
  checkGradient(b, SEAMS_CV_GYRATION_RG2, 6, "rg2, six-ring weights, tilted cell");
  checkGradient(b, SEAMS_CV_GYRATION_VOLUME, -1, "volume, uniform, tilted cell");
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

TEST_CASE("collective variables are finite or SEAMS_ERR_RANGE for every finite switch",
          "[smooth_cv]") {
  Frame square;
  square.n = 4;
  square.xyz = {5.0, 5.0, 5.0, 7.0, 5.0, 5.0, 7.0, 7.0, 5.0, 5.0, 7.0, 5.0};
  square.box9[3] = square.box9[4] = square.box9[5] = 20.0;
  {
    // Sides of 2 A at x = 200 with n = 200: R = 200^-200 underflows to 0,
    // and x^399 overflows. The square is present with weight 0.
    Api api(square, 0.002, 1.6, 200, 0, 3.0, 3.5, 4);
    api.rebuild(0.1);
    std::vector<double> grad;
    REQUIRE(api.value(SEAMS_CV_RING_COUNT, 4, &grad) == 0.0);
    REQUIRE(std::all_of(grad.begin(), grad.end(), [](double g) { return std::isfinite(g); }));
    REQUIRE(api.value(SEAMS_CV_GYRATION_RG2, 4, &grad) == 0.0);
  }
  // Sides of exactly d0 = 2 A with n = 1 and r0 = 1e-310 A: the slope there
  // is -1 / r0, beyond the largest double, and the diagonals sit at
  // x = 2.8e309. The square has weight one; its gradient is out of range.
  Api api(square, 1e-310, 2.0, 1, 0, 2.5, 3.0, 4);
  api.rebuild(0.1);
  std::size_t edges = 0;
  REQUIRE(seams_cv_topology(api.cv, &edges, nullptr, nullptr) == SEAMS_OK);
  REQUIRE(edges == 6);
  double v = -1.0;
  REQUIRE(seams_cv_value(api.cv, api.sys, SEAMS_CV_RING_COUNT, 4, &v, nullptr, 0) == SEAMS_OK);
  REQUIRE(std::abs(v - 1.0) < 1e-12);
  v = -1.0;
  std::vector<double> grad(12, -1.0);
  REQUIRE(seams_cv_value(api.cv, api.sys, SEAMS_CV_RING_COUNT, 4, &v, grad.data(),
                         grad.size()) == SEAMS_ERR_RANGE);
  REQUIRE(v == -1.0);
  REQUIRE(std::all_of(grad.begin(), grad.end(), [](double g) { return g == -1.0; }));
  REQUIRE(seams_cv_state_clause(api.cv, 0, SEAMS_CV_RING_COUNT, 4, 0.5, 0.1, 1) == SEAMS_OK);
  double chi = -1.0;
  std::size_t count = 0;
  REQUIRE(seams_cv_states(api.cv, api.sys, &chi, 1, &count, nullptr, 0) == SEAMS_OK);
  REQUIRE(std::isfinite(chi));
  REQUIRE(chi > 0.99);
  chi = -1.0;
  REQUIRE(seams_cv_states(api.cv, api.sys, &chi, 1, &count, grad.data(), grad.size()) ==
          SEAMS_ERR_RANGE);
  REQUIRE(chi == -1.0);
  REQUIRE(std::all_of(grad.begin(), grad.end(), [](double g) { return g == -1.0; }));
}

TEST_CASE("switch exponents stop at 4096", "[smooth_cv][limits]") {
  // At the limit, m = 0 selects 2n = 4096 and the value is the closed form
  // y^n / (1 + y^n) for y = 0.9. One past the limit is an error that names
  // 4096 and leaves the outputs untouched.
  const double r = kD0 + kR0 / 0.9;
  REQUIRE(r < kOn);
  const double y = 0.9;
  const double p = std::pow(y, 2048);
  const double expect = p / (1.0 + p);
  double s = 0.0;
  double ds = 0.0;
  REQUIRE(seams_switch_rational(r, kR0, kD0, 2048, 0, kOn, kMax, &s, &ds) == SEAMS_OK);
  INFO("s = " << s << ", expected " << expect << ", ds/dr = " << ds);
  REQUIRE(expect > 1e-300);
  REQUIRE(std::abs(s - expect) <= 1e-8 * expect);
  REQUIRE(std::isfinite(ds));
  REQUIRE(ds < 0.0);
  double up = 0.0;
  double down = 0.0;
  REQUIRE(seams_switch_rational(r + kStep, kR0, kD0, 2048, 0, kOn, kMax, &up, nullptr) ==
          SEAMS_OK);
  REQUIRE(seams_switch_rational(r - kStep, kR0, kD0, 2048, 0, kOn, kMax, &down, nullptr) ==
          SEAMS_OK);
  REQUIRE(std::abs((up - down) / (2.0 * kStep) - ds) <= 1e-2 * std::abs(ds));

  s = 1.25;
  ds = 2.5;
  REQUIRE(seams_switch_rational(r, kR0, kD0, 2049, 0, kOn, kMax, &s, &ds) == SEAMS_ERR_RANGE);
  REQUIRE(s == 1.25);
  REQUIRE(ds == 2.5);
  REQUIRE(std::string(seams_last_error_message()).find("4096") != std::string::npos);
  s = 1.25;
  REQUIRE(seams_switch_rational(r, kR0, kD0, 100, 4097, kOn, kMax, &s, nullptr) ==
          SEAMS_ERR_RANGE);
  REQUIRE(s == 1.25);
  REQUIRE(std::string(seams_last_error_message()).find("4096") != std::string::npos);
}

TEST_CASE("an atom type change invalidates a type-filtered selection", "[smooth_cv][limits]") {
  const double pi = std::acos(-1.0);
  Frame f;
  f.n = 3;
  const double side = 2.0;
  const double radius = side / std::sqrt(3.0);
  for (int k = 0; k < 3; ++k) {
    const double t = 2.0 * k * pi / 3.0;
    f.xyz.insert(f.xyz.end(), {10.0 + radius * std::cos(t), 10.0 + radius * std::sin(t), 10.0});
  }
  f.box9[3] = f.box9[4] = f.box9[5] = 30.0;
  Api api(f, 0.4, 2.5, 6, 12, 3.0, 3.5, 3);
  int types[3] = {1, 1, 1};
  REQUIRE(seams_system_set_atoms(api.sys, f.xyz.data(), f.n, f.box9, 9, types, nullptr) ==
          SEAMS_OK);
  REQUIRE(seams_cv_rebuild(api.cv, api.sys, 1, nullptr, 0, 0.2) == SEAMS_OK);
  REQUIRE(std::abs(api.value(SEAMS_CV_RING_COUNT, 3) - 1.0) < 1e-12);

  types[2] = 2;
  REQUIRE(seams_system_set_atoms(api.sys, f.xyz.data(), f.n, f.box9, 9, types, nullptr) ==
          SEAMS_OK);
  int stale = 0;
  REQUIRE(seams_cv_needs_rebuild(api.cv, api.sys, &stale) == SEAMS_OK);
  REQUIRE(stale == 1);
  REQUIRE(seams_cv_rebuild(api.cv, api.sys, 1, nullptr, 0, 0.2) == SEAMS_OK);
  const double after = api.value(SEAMS_CV_RING_COUNT, 3);

  Api fresh(f, 0.4, 2.5, 6, 12, 3.0, 3.5, 3);
  REQUIRE(seams_system_set_atoms(fresh.sys, f.xyz.data(), f.n, f.box9, 9, types, nullptr) ==
          SEAMS_OK);
  REQUIRE(seams_cv_rebuild(fresh.cv, fresh.sys, 1, nullptr, 0, 0.2) == SEAMS_OK);
  const double other = fresh.value(SEAMS_CV_RING_COUNT, 3);
  REQUIRE(after == other);
  REQUIRE(after == 0.0);
}

TEST_CASE("exactly and nearly planar clusters have zero volume and a matching gradient",
          "[smooth_cv][limits]") {
  // The plane z = 0.3 x + 0.2 y is not a coordinate plane, so the
  // determinant of G is a cancellation. A lift of 1e-8 along the normal
  // is coplanar on the scale of the cluster and must stay at volume 0
  // with a zero gradient. Central differences use a step that remains
  // inside that flat neighbourhood.
  const double pi = std::acos(-1.0);
  const double nx = -0.3;
  const double ny = -0.2;
  const double nz = 1.0;
  const double nn = std::sqrt(nx * nx + ny * ny + nz * nz);
  auto cloud = [&](double lift) {
    Frame f;
    f.n = 6;
    for (int k = 0; k < 6; ++k) {
      const double t = k * pi / 3.0;
      const double x = 10.0 + 2.75 * std::cos(t);
      const double y = 10.0 + 2.75 * std::sin(t);
      f.xyz.insert(f.xyz.end(), {x, y, 0.3 * x + 0.2 * y});
    }
    f.xyz[0] += lift * nx / nn;
    f.xyz[1] += lift * ny / nn;
    f.xyz[2] += lift * nz / nn;
    f.box9[3] = f.box9[4] = f.box9[5] = 40.0;
    return f;
  };
  for (const double lift : {0.0, 1e-8}) {
    const Frame f = cloud(lift);
    Api api(f, kR0, kD0, kN, kM, 6.0, 8.0, 6);
    api.rebuild(0.4);
    std::vector<double> grad;
    const double v = api.value(SEAMS_CV_GYRATION_VOLUME, -1, &grad);
    INFO("lift " << lift << ", V " << v);
    REQUIRE(v == 0.0);
    REQUIRE(std::all_of(grad.begin(), grad.end(), [](double g) { return g == 0.0; }));
    const auto fd = centralDifference(api, [&] { return api.value(SEAMS_CV_GYRATION_VOLUME, -1); });
    REQUIRE(fd.size() == grad.size());
    for (std::size_t c = 0; c < fd.size(); ++c) {
      INFO("component " << c << " fd " << fd[c]);
      REQUIRE(std::abs(fd[c]) < 1e-8);
    }
  }
}

TEST_CASE("a dense neighbour graph stops at the ring-search limit", "[smooth_cv][limits]") {
  // Eight mutually bonded atoms are inside the budget. Every triangle has
  // weight 1 and every longer cycle has a chord, so the ring count is
  // C(8, 3) = 56. Thirty-two mutually bonded atoms with max_ring 7 would
  // visit more than 20000000 ring-search steps; the rebuild names that
  // limit and does not return a partial count.
  auto clique = [](int n, double spacing) {
    Frame f;
    f.n = static_cast<std::size_t>(n);
    const int side = 4;
    for (int i = 0; i < n; ++i) {
      f.xyz.insert(f.xyz.end(), {10.0 + spacing * static_cast<double>(i % side),
                                 10.0 + spacing * static_cast<double>((i / side) % side),
                                 10.0 + spacing * static_cast<double>(i / (side * side))});
    }
    f.box9[3] = f.box9[4] = f.box9[5] = 40.0;
    return f;
  };
  {
    const Frame small = clique(8, 0.4);
    Api api(small, 0.5, 3.0, 6, 12, 4.0, 5.0, 5);
    api.rebuild(0.2);
    REQUIRE(std::abs(api.value(SEAMS_CV_RING_COUNT, 3) - 56.0) < 1e-9);
    REQUIRE(api.value(SEAMS_CV_RING_COUNT, 4) == 0.0);
    REQUIRE(api.value(SEAMS_CV_RING_COUNT, 5) == 0.0);
    REQUIRE(std::abs(api.value(SEAMS_CV_RING_COUNT, 0) - 56.0) < 1e-9);
  }
  const Frame dense = clique(32, 0.45);
  Api api(dense, 0.5, 4.0, 6, 12, 5.0, 6.0, 7);
  double value = 1.0;
  const seams_status rc =
      seams_cv_rebuild(api.cv, api.sys, 0, nullptr, 0, 0.2);
  INFO(seams_last_error_message());
  REQUIRE(rc == SEAMS_ERR_RANGE);
  REQUIRE(value == 1.0);
  const std::string msg = seams_last_error_message();
  REQUIRE(msg.find("20000000") != std::string::npos);
  REQUIRE(msg.find("steps") != std::string::npos);
  double kept = -1.0;
  REQUIRE(seams_cv_value(api.cv, api.sys, SEAMS_CV_RING_COUNT, 0, &kept, nullptr, 0) ==
          SEAMS_ERR_STATE);
}
