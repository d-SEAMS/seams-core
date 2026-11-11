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
    const __m512d inv = _mm512_mask_div_pd(_mm512_setzero_pd(), live, _mm512_set1_pd(1.0),
                                           _mm512_sqrt_pd(r2));
    V x(_mm512_mul_pd(x0, inv));
    V y(_mm512_mul_pd(y0, inv));
    V z(_mm512_mul_pd(z0, inv));
    V reV[13];
    V imV[13];
    evalUnit(orderL, x, y, z, reV, imV);
    // Dead lanes (zero displacement) stay zero, matching ylmReal.
    for (int m = 0; m < nComp; ++m) {
      reV[m].v = _mm512_mask_blend_pd(live, _mm512_setzero_pd(), reV[m].v);
      imV[m].v = _mm512_mask_blend_pd(live, _mm512_setzero_pd(), imV[m].v);
    }
    storeLanes(reV, imV, nComp, b, 8, re, im);
  }
  for (; b < n; ++b) {
    ylmReal(orderL, dx[b], dy[b], dz[b],
            re + static_cast<std::size_t>(b) * static_cast<std::size_t>(nComp),
            im + static_cast<std::size_t>(b) * static_cast<std::size_t>(nComp));
  }
}

} // namespace ylm
} // namespace seams
