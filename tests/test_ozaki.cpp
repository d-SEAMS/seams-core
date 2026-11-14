#include <ozaki.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <mpfr.h>

#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <vector>

namespace {

double mpfrDot(const double *x, const double *y, int n) {
  mpfr_t acc, a, b, p;
  mpfr_init2(acc, 4096);
  mpfr_init2(a, 4096);
  mpfr_init2(b, 4096);
  mpfr_init2(p, 4096);
  mpfr_set_zero(acc, 1);
  for (int i = 0; i < n; i++) {
    mpfr_set_d(a, x[i], MPFR_RNDN);
    mpfr_set_d(b, y[i], MPFR_RNDN);
    mpfr_mul(p, a, b, MPFR_RNDN);
    mpfr_add(acc, acc, p, MPFR_RNDN);
  }
  const double out = mpfr_get_d(acc, MPFR_RNDN);
  mpfr_clears(acc, a, b, p, static_cast<mpfr_ptr>(nullptr));
  return out;
}

double mpfrAbsErr(const double *x, const double *y, int n, double value) {
  mpfr_t acc, a, b, p, got, err;
  mpfr_init2(acc, 4096);
  mpfr_init2(a, 4096);
  mpfr_init2(b, 4096);
  mpfr_init2(p, 4096);
  mpfr_init2(got, 4096);
  mpfr_init2(err, 4096);
  mpfr_set_zero(acc, 1);
  for (int i = 0; i < n; i++) {
    mpfr_set_d(a, x[i], MPFR_RNDN);
    mpfr_set_d(b, y[i], MPFR_RNDN);
    mpfr_mul(p, a, b, MPFR_RNDN);
    mpfr_add(acc, acc, p, MPFR_RNDN);
  }
  mpfr_set_d(got, value, MPFR_RNDN);
  mpfr_sub(err, got, acc, MPFR_RNDN);
  mpfr_abs(err, err, MPFR_RNDN);
  const double out = mpfr_get_d(err, MPFR_RNDU);
  mpfr_clears(acc, a, b, p, got, err, static_cast<mpfr_ptr>(nullptr));
  return out;
}

void fill(std::vector<double> &x, unsigned seed, double scale) {
  std::srand(seed);
  for (double &v : x) {
    v = scale * ((std::rand() % 20001) / 10000.0 - 1.0);
  }
}

} // namespace

TEST_CASE("TwoProd matches the exact product at 4096 bits", "[ozaki]") {
  const double samples[][2] = {{1.25, -0.5}, {3.141592653589793, 2.718281828459045},
                               {1e-8, 1e-6}, {-1024.5, 7.0}};
  for (const auto &s : samples) {
    double hi = 0.0;
    double lo = 0.0;
    seams::product::twoProd(s[0], s[1], hi, lo);
    mpfr_t a, b, p, h, l, sum;
    mpfr_init2(a, 256);
    mpfr_init2(b, 256);
    mpfr_init2(p, 256);
    mpfr_init2(h, 256);
    mpfr_init2(l, 256);
    mpfr_init2(sum, 256);
    mpfr_set_d(a, s[0], MPFR_RNDN);
    mpfr_set_d(b, s[1], MPFR_RNDN);
    mpfr_mul(p, a, b, MPFR_RNDN);
    mpfr_set_d(h, hi, MPFR_RNDN);
    mpfr_set_d(l, lo, MPFR_RNDN);
    mpfr_add(sum, h, l, MPFR_RNDN);
    REQUIRE(mpfr_equal_p(sum, p));
    mpfr_clears(a, b, p, h, l, sum, static_cast<mpfr_ptr>(nullptr));
  }
}

TEST_CASE("certified dots stay inside their MPFR bounds", "[ozaki]") {
  std::vector<double> x(14);
  std::vector<double> y(14);
  fill(x, 11, 1.0);
  fill(y, 29, 1.0);
  for (auto mode : {seams::product::Mode::Dot2, seams::product::Mode::OzakiI,
                    seams::product::Mode::OzakiII}) {
    const auto d = seams::product::realDot(x.data(), y.data(), 14, mode);
    REQUIRE(d.certified);
    const double err = mpfrAbsErr(x.data(), y.data(), 14, d.value);
    REQUIRE(err <= d.bound);
  }
  std::vector<double> wide(8);
  std::vector<double> wideY(8);
  fill(wide, 3, 1.0);
  fill(wideY, 4, 1e-6);
  const auto wideDot =
      seams::product::realDot(wide.data(), wideY.data(), 8, seams::product::Mode::OzakiI);
  REQUIRE(wideDot.certified);
  REQUIRE(mpfrAbsErr(wide.data(), wideY.data(), 8, wideDot.value) <= wideDot.bound);
}

TEST_CASE("pairwise sum does not depend on a left fold and is repeatable", "[ozaki]") {
  const double x[] = {1.0, 1.0, 1e20, -1e20};
  const double left = seams::product::leftSum(x, 4);
  const double pair = seams::product::pairwiseSum(x, 4);
  REQUIRE(pair == seams::product::pairwiseSum(x, 4));
  REQUIRE(pair != left);
  mpfr_t acc;
  mpfr_init2(acc, 256);
  mpfr_set_zero(acc, 1);
  for (double v : x) {
    mpfr_add_d(acc, acc, v, MPFR_RNDN);
  }
  const double exact = mpfr_get_d(acc, MPFR_RNDN);
  mpfr_clear(acc);
  REQUIRE(pair == Catch::Approx(exact));
  REQUIRE(std::fabs(pair - exact) <= seams::product::pairwiseBound(4, 2e20L + 2.0L));
}

TEST_CASE("INT8 tiles match the scalar product", "[ozaki]") {
  constexpr int m = 32;
  constexpr int n = 32;
  constexpr int k = 64;
  std::vector<signed char> A(m * k);
  std::vector<signed char> B(k * n);
  std::srand(7);
  for (auto &v : A) {
    v = static_cast<signed char>((std::rand() % 255) - 128);
  }
  for (auto &v : B) {
    v = static_cast<signed char>((std::rand() % 255) - 128);
  }
  std::vector<int> scalar(m * n, 0);
  std::vector<int> vnni(m * n, 0);
  std::vector<int> amx(m * n, 0);
  REQUIRE(seams::product::int8Gemm(m, n, k, A.data(), k, B.data(), n, scalar.data(), n,
                                   seams::product::Int8Backend::Scalar));
  if (seams::product::vnniAvailable()) {
    REQUIRE(seams::product::int8Gemm(m, n, k, A.data(), k, B.data(), n, vnni.data(), n,
                                     seams::product::Int8Backend::Vnni));
    REQUIRE(vnni == scalar);
  }
  if (seams::product::amxAvailable()) {
    REQUIRE(seams::product::int8Gemm(m, n, k, A.data(), k, B.data(), n, amx.data(), n,
                                     seams::product::Int8Backend::Amx));
    REQUIRE(amx == scalar);
  }
}

TEST_CASE("Ozaki-II GEMM matches MPFR inside the certificate", "[ozaki]") {
  constexpr int m = 6;
  constexpr int n = 5;
  constexpr int k = 9;
  std::vector<double> A(m * k);
  std::vector<double> B(k * n);
  fill(A, 19, 0.5);
  fill(B, 23, 1.5);
  std::vector<double> C(m * n);
  std::vector<double> bound(m * n);
  const auto info =
      seams::product::gemm(m, n, k, A.data(), k, B.data(), n, C.data(), n,
                           seams::product::Mode::OzakiII, bound.data(),
                           seams::product::Int8Backend::Scalar);
  REQUIRE(info.ok);
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      std::vector<double> row(k);
      std::vector<double> col(k);
      for (std::size_t t = 0; t < k; t++) {
        row[t] = A[i * k + t];
        col[t] = B[t * n + j];
      }
      const double err = mpfrAbsErr(row.data(), col.data(), k, C[i * n + j]);
      REQUIRE(err <= bound[i * n + j]);
    }
  }
  if (seams::product::amxAvailable()) {
    std::vector<double> C2(m * n);
    std::vector<double> bound2(m * n);
    const auto amx = seams::product::gemm(m, n, k, A.data(), k, B.data(), n, C2.data(), n,
                                          seams::product::Mode::OzakiII, bound2.data(),
                                          seams::product::Int8Backend::Amx);
    REQUIRE(amx.ok);
    for (int i = 0; i < m * n; i++) {
      REQUIRE(std::fabs(C2[static_cast<std::size_t>(i)] - C[static_cast<std::size_t>(i)]) <=
              bound[static_cast<std::size_t>(i)] + bound2[static_cast<std::size_t>(i)]);
    }
  }
  (void)mpfrDot;
}

TEST_CASE("bond correlation certificate holds for an l=3 pair", "[ozaki]") {
  std::vector<std::complex<double>> a(7);
  std::vector<std::complex<double>> b(7);
  std::srand(5);
  for (int i = 0; i < 7; i++) {
    a[static_cast<std::size_t>(i)] = {std::rand() / 1e5, std::rand() / 1e5};
    b[static_cast<std::size_t>(i)] = {std::rand() / 1e5, -std::rand() / 1e5};
  }
  const auto cij =
      seams::product::bondCorrelation(a.data(), b.data(), 7, seams::product::Mode::OzakiI);
  REQUIRE(cij.certified);
  std::vector<double> ar(7), ai(7), br(7), bi(7);
  for (int i = 0; i < 7; i++) {
    ar[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i)].real();
    ai[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i)].imag();
    br[static_cast<std::size_t>(i)] = b[static_cast<std::size_t>(i)].real();
    bi[static_cast<std::size_t>(i)] = b[static_cast<std::size_t>(i)].imag();
  }
  mpfr_t acc, t, u, p;
  mpfr_init2(acc, 4096);
  mpfr_init2(t, 4096);
  mpfr_init2(u, 4096);
  mpfr_init2(p, 4096);
  mpfr_set_zero(acc, 1);
  auto addDot = [&](const std::vector<double> &x, const std::vector<double> &y, int sign) {
    for (int i = 0; i < 7; i++) {
      mpfr_set_d(t, x[static_cast<std::size_t>(i)], MPFR_RNDN);
      mpfr_set_d(u, y[static_cast<std::size_t>(i)], MPFR_RNDN);
      mpfr_mul(p, t, u, MPFR_RNDN);
      if (sign < 0) {
        mpfr_neg(p, p, MPFR_RNDN);
      }
      mpfr_add(acc, acc, p, MPFR_RNDN);
    }
  };
  addDot(ar, br, 1);
  addDot(ai, bi, 1);
  mpfr_t nA, nB;
  mpfr_init2(nA, 4096);
  mpfr_init2(nB, 4096);
  mpfr_set_zero(nA, 1);
  mpfr_set_zero(nB, 1);
  for (int i = 0; i < 7; i++) {
    mpfr_set_d(t, ar[static_cast<std::size_t>(i)], MPFR_RNDN);
    mpfr_sqr(p, t, MPFR_RNDN);
    mpfr_add(nA, nA, p, MPFR_RNDN);
    mpfr_set_d(t, ai[static_cast<std::size_t>(i)], MPFR_RNDN);
    mpfr_sqr(p, t, MPFR_RNDN);
    mpfr_add(nA, nA, p, MPFR_RNDN);
    mpfr_set_d(t, br[static_cast<std::size_t>(i)], MPFR_RNDN);
    mpfr_sqr(p, t, MPFR_RNDN);
    mpfr_add(nB, nB, p, MPFR_RNDN);
    mpfr_set_d(t, bi[static_cast<std::size_t>(i)], MPFR_RNDN);
    mpfr_sqr(p, t, MPFR_RNDN);
    mpfr_add(nB, nB, p, MPFR_RNDN);
  }
  mpfr_mul(p, nA, nB, MPFR_RNDN);
  mpfr_sqrt(p, p, MPFR_RNDN);
  mpfr_div(acc, acc, p, MPFR_RNDN);
  mpfr_t got, err;
  mpfr_init2(got, 4096);
  mpfr_init2(err, 4096);
  mpfr_set_d(got, cij.value, MPFR_RNDN);
  mpfr_sub(err, got, acc, MPFR_RNDN);
  mpfr_abs(err, err, MPFR_RNDN);
  const double e = mpfr_get_d(err, MPFR_RNDU);
  REQUIRE(e <= cij.bound);
  mpfr_clears(acc, t, u, p, nA, nB, got, err, static_cast<mpfr_ptr>(nullptr));
}
