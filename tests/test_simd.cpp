/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** SIMD verification tests: compare BatchPeriodicDistSq against a scalar
** reference for various input sizes, including edge cases for the SIMD
** remainder loop.
*/

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <simd_distance.hpp>

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

namespace {

/// Scalar reference: squared periodic distance for one atom pair.
double scalarPeriodicDistSq(double dx, double dy, double dz,
                            double bx, double by, double bz) {
  double ddx = std::fabs(dx);
  double ddy = std::fabs(dy);
  double ddz = std::fabs(dz);
  ddx -= bx * std::round(ddx / bx);
  ddy -= by * std::round(ddy / by);
  ddz -= bz * std::round(ddz / bz);
  return ddx * ddx + ddy * ddy + ddz * ddz;
}

/// Fill a vector with uniformly distributed random doubles in [lo, hi].
void fillRandom(std::vector<double> &v, std::mt19937_64 &rng,
                double lo, double hi) {
  std::uniform_real_distribution<double> dist(lo, hi);
  for (auto &x : v) x = dist(rng);
}

/// Compare SIMD output against scalar reference for a given size.
void checkSize(size_t n) {
  std::mt19937_64 rng(42 + n);  // deterministic per size
  double bx = 25.0, by = 30.0, bz = 35.0;

  std::vector<double> dx(n), dy(n), dz(n);
  fillRandom(dx, rng, -50.0, 50.0);
  fillRandom(dy, rng, -50.0, 50.0);
  fillRandom(dz, rng, -50.0, 50.0);

  std::vector<double> simd_out(n, -1.0);
  seams::BatchPeriodicDistSq(dx.data(), dy.data(), dz.data(),
                             bx, by, bz, simd_out.data(), n);

  for (size_t i = 0; i < n; i++) {
    double ref = scalarPeriodicDistSq(dx[i], dy[i], dz[i], bx, by, bz);
    REQUIRE_THAT(simd_out[i],
                 Catch::Matchers::WithinRel(ref, 1e-12) ||
                 Catch::Matchers::WithinAbs(ref, 1e-15));
  }
}

}  // namespace

TEST_CASE("SIMD matches scalar for size 1 (pure remainder)", "[simd]") {
  checkSize(1);
}

TEST_CASE("SIMD matches scalar for size 4", "[simd]") {
  checkSize(4);
}

TEST_CASE("SIMD matches scalar for size 7 (odd, partial vector)", "[simd]") {
  checkSize(7);
}

TEST_CASE("SIMD matches scalar for size 8", "[simd]") {
  checkSize(8);
}

TEST_CASE("SIMD matches scalar for size 15", "[simd]") {
  checkSize(15);
}

TEST_CASE("SIMD matches scalar for size 16", "[simd]") {
  checkSize(16);
}

TEST_CASE("SIMD matches scalar for size 100", "[simd]") {
  checkSize(100);
}

TEST_CASE("SIMD matches scalar for size 1000", "[simd]") {
  checkSize(1000);
}

TEST_CASE("SIMD handles zero-length input", "[simd]") {
  double dummy;
  // Should not crash or write anything
  seams::BatchPeriodicDistSq(nullptr, nullptr, nullptr,
                             10.0, 10.0, 10.0, &dummy, 0);
}

TEST_CASE("SIMD handles zero displacements", "[simd]") {
  constexpr size_t n = 8;
  std::vector<double> dx(n, 0.0), dy(n, 0.0), dz(n, 0.0);
  std::vector<double> out(n, -1.0);
  seams::BatchPeriodicDistSq(dx.data(), dy.data(), dz.data(),
                             10.0, 10.0, 10.0, out.data(), n);
  for (size_t i = 0; i < n; i++) {
    REQUIRE_THAT(out[i], Catch::Matchers::WithinAbs(0.0, 1e-15));
  }
}

TEST_CASE("SIMD handles exact half-box displacements", "[simd]") {
  // At exactly half the box, the wrapped distance should be ~0 or bx/2
  // depending on rounding.  Verify consistency with scalar.
  constexpr size_t n = 3;
  double bx = 10.0, by = 10.0, bz = 10.0;
  std::vector<double> dx = {5.0, -5.0, 2.5};
  std::vector<double> dy = {5.0,  0.0, 5.0};
  std::vector<double> dz = {0.0,  5.0, 5.0};
  std::vector<double> out(n);

  seams::BatchPeriodicDistSq(dx.data(), dy.data(), dz.data(),
                             bx, by, bz, out.data(), n);
  for (size_t i = 0; i < n; i++) {
    double ref = scalarPeriodicDistSq(dx[i], dy[i], dz[i], bx, by, bz);
    REQUIRE_THAT(out[i], Catch::Matchers::WithinAbs(ref, 1e-15));
  }
}

TEST_CASE("ortho SIMD wrap on dump bound spans misses the a-image", "[simd]") {
  // Sheared pair (0.2, 0.1, 1) and (9.7, 0.1, 1) on dump box
  // {15, 8.66, 10, xy=5}. H-image r^2 is 0.25. Independent wrap
  // on the bound span 15 yields 5.5^2 = 30.25.
  constexpr size_t n = 1;
  const double dx[n] = {0.2 - 9.7};
  const double dy[n] = {0.1 - 0.1};
  const double dz[n] = {1.0 - 1.0};
  double out[n] = {-1.0};
  seams::BatchPeriodicDistSq(dx, dy, dz, 15.0, 8.660254037844386, 10.0, out, n);
  REQUIRE_THAT(out[0], Catch::Matchers::WithinAbs(30.25, 1e-12));
  REQUIRE(std::abs(out[0] - 0.25) > 1.0);
}

TEST_CASE("SIMD with negative displacements matches scalar", "[simd]") {
  constexpr size_t n = 16;
  std::mt19937_64 rng(999);
  double bx = 20.0, by = 20.0, bz = 20.0;

  std::vector<double> dx(n), dy(n), dz(n);
  // All negative
  std::uniform_real_distribution<double> dist(-40.0, -0.1);
  for (size_t i = 0; i < n; i++) {
    dx[i] = dist(rng);
    dy[i] = dist(rng);
    dz[i] = dist(rng);
  }

  std::vector<double> out(n);
  seams::BatchPeriodicDistSq(dx.data(), dy.data(), dz.data(),
                             bx, by, bz, out.data(), n);
  for (size_t i = 0; i < n; i++) {
    double ref = scalarPeriodicDistSq(dx[i], dy[i], dz[i], bx, by, bz);
    REQUIRE_THAT(out[i],
                 Catch::Matchers::WithinRel(ref, 1e-12) ||
                 Catch::Matchers::WithinAbs(ref, 1e-15));
  }
}
