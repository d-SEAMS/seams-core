//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <ylm_cartesian.hpp>

#include <immintrin.h>

#include <algorithm>
#include <cstddef>

namespace seams {
namespace ylm {

namespace {

struct V {
  __m512d v;
  V() = default;
  explicit V(__m512d x) : v(x) {}
  V(double s) : v(_mm512_set1_pd(s)) {}
  friend V operator+(V a, V b) { return V(_mm512_add_pd(a.v, b.v)); }
  friend V operator-(V a, V b) { return V(_mm512_sub_pd(a.v, b.v)); }
  friend V operator*(V a, V b) { return V(_mm512_mul_pd(a.v, b.v)); }
  friend V operator-(V a) { return V(_mm512_sub_pd(_mm512_setzero_pd(), a.v)); }
};

void storeLanes(const V *reV, const V *imV, int nComp, int b, int lanes,
                double *re, double *im) {
  alignas(64) double rt[8];
  alignas(64) double it[8];
  for (int m = 0; m < nComp; ++m) {
    _mm512_store_pd(rt, reV[m].v);
    _mm512_store_pd(it, imV[m].v);
    for (int lane = 0; lane < lanes; ++lane) {
      const std::size_t at =
          (static_cast<std::size_t>(b + lane) * static_cast<std::size_t>(nComp)) +
          static_cast<std::size_t>(m);
      re[at] = rt[lane];
      im[at] = it[lane];
    }
  }
}

} // namespace

void evalSolidLanes(int orderL, __m512d x0, __m512d y0, __m512d z0, __mmask8 live,
                    V *reV, V *imV) {
  const __m512d r2 = _mm512_add_pd(
      _mm512_add_pd(_mm512_mul_pd(x0, x0), _mm512_mul_pd(y0, y0)),
      _mm512_mul_pd(z0, z0));
  const __m512d zero = _mm512_setzero_pd();
  const __m512d one = _mm512_set1_pd(1.0);
  __m512d inv;
  if (orderL == 6) {
    const __m512d r4 = _mm512_mul_pd(r2, r2);
    inv = _mm512_mask_div_pd(zero, live, one, _mm512_mul_pd(r4, r2));
  } else if (orderL == 4) {
    inv = _mm512_mask_div_pd(zero, live, one, _mm512_mul_pd(r2, r2));
  } else {
    inv = _mm512_mask_div_pd(zero, live, one, _mm512_mul_pd(r2, _mm512_sqrt_pd(r2)));
  }
  const double *norm = (orderL == 3) ? N3 : (orderL == 4) ? N4 : N6;
  solid::eval(orderL, V(x0), V(y0), V(z0), V(r2), V(inv), norm, reV, imV);
  const int nComp = 2 * orderL + 1;
  for (int m = 0; m < nComp; ++m) {
    reV[m].v = _mm512_mask_blend_pd(live, zero, reV[m].v);
    imV[m].v = _mm512_mask_blend_pd(live, zero, imV[m].v);
  }
}

double hsum4(__m256d v) {
  __m128d lo = _mm256_castpd256_pd128(v);
  __m128d hi = _mm256_extractf128_pd(v, 1);
  lo = _mm_add_pd(lo, hi);
  lo = _mm_add_sd(lo, _mm_unpackhi_pd(lo, lo));
  return _mm_cvtsd_f64(lo);
}

void ylmAdd8(int orderL, const double *dx, const double *dy, const double *dz,
             double *accA, double *accB) {
  if (orderL != 4 && orderL != 6) {
    return;
  }
  const int nComp = 2 * orderL + 1;
  const __m512d x0 = _mm512_loadu_pd(dx);
  const __m512d y0 = _mm512_loadu_pd(dy);
  const __m512d z0 = _mm512_loadu_pd(dz);
  V reV[13];
  V imV[13];
  evalSolidLanes(orderL, x0, y0, z0, static_cast<__mmask8>(0xff), reV, imV);
  for (int m = 0; m < nComp; ++m) {
    const __m256d reLo = _mm512_castpd512_pd256(reV[m].v);
    const __m256d reHi = _mm512_extractf64x4_pd(reV[m].v, 1);
    const __m256d imLo = _mm512_castpd512_pd256(imV[m].v);
    const __m256d imHi = _mm512_extractf64x4_pd(imV[m].v, 1);
    accA[2 * m] += hsum4(reLo);
    accA[2 * m + 1] += hsum4(imLo);
    accB[2 * m] += hsum4(reHi);
    accB[2 * m + 1] += hsum4(imHi);
  }
}

void ylmBatchAvx512(int orderL, const double *dx, const double *dy,
                    const double *dz, int n, double *re, double *im) {
  const int nComp = 2 * orderL + 1;
  int b = 0;
  for (; b + 8 <= n; b += 8) {
    const __m512d x0 = _mm512_loadu_pd(dx + b);
    const __m512d y0 = _mm512_loadu_pd(dy + b);
    const __m512d z0 = _mm512_loadu_pd(dz + b);
    const __m512d r2 = _mm512_add_pd(
        _mm512_add_pd(_mm512_mul_pd(x0, x0), _mm512_mul_pd(y0, y0)),
        _mm512_mul_pd(z0, z0));
    const __mmask8 live = _mm512_cmp_pd_mask(r2, _mm512_setzero_pd(), _CMP_GT_OQ);
    V reV[13];
    V imV[13];
    evalSolidLanes(orderL, x0, y0, z0, live, reV, imV);
    storeLanes(reV, imV, nComp, b, 8, re, im);
  }
  if (b < n) {
    const int tail = n - b;
    const __mmask8 live = static_cast<__mmask8>((1u << tail) - 1u);
    const __m512d x0 = _mm512_maskz_loadu_pd(live, dx + b);
    const __m512d y0 = _mm512_maskz_loadu_pd(live, dy + b);
    const __m512d z0 = _mm512_maskz_loadu_pd(live, dz + b);
    const __m512d r2 = _mm512_add_pd(
        _mm512_add_pd(_mm512_mul_pd(x0, x0), _mm512_mul_pd(y0, y0)),
        _mm512_mul_pd(z0, z0));
    const __mmask8 ok = _mm512_mask_cmp_pd_mask(live, r2, _mm512_setzero_pd(), _CMP_GT_OQ);
    V reV[13];
    V imV[13];
    evalSolidLanes(orderL, x0, y0, z0, ok, reV, imV);
    storeLanes(reV, imV, nComp, b, tail, re, im);
  }
}

} // namespace ylm
} // namespace seams
