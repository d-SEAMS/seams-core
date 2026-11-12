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

void transpose8(__m512d in[8], __m512d out[8]) {
  // in[component] holds eight atoms. Unpack pairs of components, then
  // gather the 128-bit lane that belongs to one atom.
  __m512d pair[8];
  pair[0] = _mm512_unpacklo_pd(in[0], in[1]);
  pair[1] = _mm512_unpackhi_pd(in[0], in[1]);
  pair[2] = _mm512_unpacklo_pd(in[2], in[3]);
  pair[3] = _mm512_unpackhi_pd(in[2], in[3]);
  pair[4] = _mm512_unpacklo_pd(in[4], in[5]);
  pair[5] = _mm512_unpackhi_pd(in[4], in[5]);
  pair[6] = _mm512_unpacklo_pd(in[6], in[7]);
  pair[7] = _mm512_unpackhi_pd(in[6], in[7]);
  for (int atom = 0; atom < 8; ++atom) {
    const int odd = atom & 1;
    const int lane = atom >> 1;
    // The extract immediate is a compile-time lane index.
    const auto laneOf = [](__m512d v, int which) {
      switch (which) {
      case 0:
        return _mm512_extractf64x2_pd(v, 0);
      case 1:
        return _mm512_extractf64x2_pd(v, 1);
      case 2:
        return _mm512_extractf64x2_pd(v, 2);
      default:
        return _mm512_extractf64x2_pd(v, 3);
      }
    };
    __m512d row = _mm512_setzero_pd();
    row = _mm512_insertf64x2(row, laneOf(pair[odd], lane), 0);
    row = _mm512_insertf64x2(row, laneOf(pair[2 + odd], lane), 1);
    row = _mm512_insertf64x2(row, laneOf(pair[4 + odd], lane), 2);
    row = _mm512_insertf64x2(row, laneOf(pair[6 + odd], lane), 3);
    out[atom] = row;
  }
}

void storePairs(double *dst, __m512d re, __m512d im, int nComp) {
  const __m512d u0 = _mm512_unpacklo_pd(re, im);
  const __m512d u1 = _mm512_unpackhi_pd(re, im);
  const __m512i idx = _mm512_set_epi64(7, 6, 3, 2, 5, 4, 1, 0);
  const __m512d low =
      _mm512_permutexvar_pd(idx, _mm512_shuffle_f64x2(u0, u1, 0x44));
  const __m512d high =
      _mm512_permutexvar_pd(idx, _mm512_shuffle_f64x2(u0, u1, 0xEE));
  const int lowN = std::min(8, nComp * 2);
  if (lowN >= 8) {
    _mm512_storeu_pd(dst, low);
  } else if (lowN > 0) {
    _mm512_mask_storeu_pd(dst, static_cast<__mmask8>((1u << lowN) - 1u), low);
  }
  if (nComp > 4) {
    const int highN = (nComp - 4) * 2;
    if (highN >= 8) {
      _mm512_storeu_pd(dst + 8, high);
    } else if (highN > 0) {
      _mm512_mask_storeu_pd(dst + 8, static_cast<__mmask8>((1u << highN) - 1u),
                            high);
    }
  }
}

bool ylmSumRows4(int orderL, const double *dx, const double *dy,
                 const double *dz, double *rows, int stride) {
  if ((orderL != 4 && orderL != 6) || stride < 2 * (2 * orderL + 1)) {
    return false;
  }
  const int nComp = 2 * orderL + 1;
  V accRe[13];
  V accIm[13];
  const V zero(0.0);
  for (int m = 0; m < 13; ++m) {
    accRe[m] = zero;
    accIm[m] = zero;
  }
  for (int slot = 0; slot < 4; ++slot) {
    const __m512d x0 = _mm512_loadu_pd(dx + slot * 8);
    const __m512d y0 = _mm512_loadu_pd(dy + slot * 8);
    const __m512d z0 = _mm512_loadu_pd(dz + slot * 8);
    const __m512d r2 = _mm512_add_pd(
        _mm512_add_pd(_mm512_mul_pd(x0, x0), _mm512_mul_pd(y0, y0)),
        _mm512_mul_pd(z0, z0));
    const __mmask8 live = _mm512_cmp_pd_mask(r2, _mm512_setzero_pd(), _CMP_GT_OQ);
    if (live != static_cast<__mmask8>(0xff)) {
      return false;
    }
    V reV[13];
    V imV[13];
    evalSolidLanes(orderL, x0, y0, z0, live, reV, imV);
    for (int m = 0; m < nComp; ++m) {
      accRe[m] = accRe[m] + reV[m];
      accIm[m] = accIm[m] + imV[m];
    }
  }
  const __m512d zreg = _mm512_setzero_pd();
  for (int block = 0; block < nComp; block += 8) {
    const int have = std::min(8, nComp - block);
    __m512d reIn[8];
    __m512d imIn[8];
    __m512d reOut[8];
    __m512d imOut[8];
    for (int k = 0; k < 8; ++k) {
      const int m = block + k;
      reIn[k] = (m < nComp) ? accRe[m].v : zreg;
      imIn[k] = (m < nComp) ? accIm[m].v : zreg;
    }
    transpose8(reIn, reOut);
    transpose8(imIn, imOut);
    for (int atom = 0; atom < 8; ++atom) {
      storePairs(rows + static_cast<std::size_t>(atom) * static_cast<std::size_t>(stride) +
                     static_cast<std::size_t>(2 * block),
                 reOut[atom], imOut[atom], have);
    }
  }
  return true;
}

bool ylmLocalQlAvx(int orderL, const double *dx, const double *dy, const double *dz,
                   int nbonds, const int *nUsed, double prefactor, double *ql) {
  if ((orderL != 4 && orderL != 6) || nbonds < 0) {
    return false;
  }
  const int nComp = 2 * orderL + 1;
  V accRe[13];
  V accIm[13];
  const V zero(0.0);
  for (int m = 0; m < 13; ++m) {
    accRe[m] = zero;
    accIm[m] = zero;
  }
  for (int bond = 0; bond < nbonds; ++bond) {
    unsigned bits = 0;
    for (int a = 0; a < 8; ++a) {
      if (bond < nUsed[a]) {
        bits |= 1u << a;
      }
    }
    if (bits == 0) {
      continue;
    }
    const __m512d x0 = _mm512_loadu_pd(dx + bond * 8);
    const __m512d y0 = _mm512_loadu_pd(dy + bond * 8);
    const __m512d z0 = _mm512_loadu_pd(dz + bond * 8);
    V reV[13];
    V imV[13];
    evalSolidLanes(orderL, x0, y0, z0, static_cast<__mmask8>(bits), reV, imV);
    for (int m = 0; m < nComp; ++m) {
      accRe[m] = accRe[m] + reV[m];
      accIm[m] = accIm[m] + imV[m];
    }
  }
  alignas(64) double reT[13][8];
  alignas(64) double imT[13][8];
  for (int m = 0; m < nComp; ++m) {
    _mm512_store_pd(reT[m], accRe[m].v);
    _mm512_store_pd(imT[m], accIm[m].v);
  }
  for (int a = 0; a < 8; ++a) {
    const int used = nUsed[a];
    if (used <= 0) {
      ql[a] = 0.0;
      continue;
    }
    const double inv = 1.0 / static_cast<double>(used);
    double sum = 0.0;
    for (int m = 0; m < nComp; ++m) {
      const double r = reT[m][a] * inv;
      const double i = imT[m][a] * inv;
      sum += r * r + i * i;
    }
    ql[a] = std::sqrt(prefactor * sum);
  }
  return true;
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
