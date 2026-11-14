//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <cpu_features.hpp>
#include <ylm_batch.hpp>

#include <cstddef>

namespace seams {
namespace ylm {

namespace {

void ylmBatchScalar(int orderL, const double *dx, const double *dy,
                    const double *dz, int n, double *re, double *im) {
  const int nComp = 2 * orderL + 1;
  for (int b = 0; b < n; ++b) {
    ylmReal(orderL, dx[b], dy[b], dz[b], re + static_cast<std::size_t>(b) * nComp,
            im + static_cast<std::size_t>(b) * nComp);
  }
}

} // namespace

#if defined(SEAMS_YLM_AVX512)
void ylmBatchAvx512(int orderL, const double *dx, const double *dy,
                    const double *dz, int n, double *re, double *im);
#endif

#if defined(SEAMS_YLM_AVX512)
bool ylmLocalQlAvx(int orderL, const double *dx, const double *dy, const double *dz,
                   int nbonds, const int *nUsed, double prefactor, double *ql);
#endif

bool ylmLocalQl(int orderL, const double *dx, const double *dy, const double *dz,
                int nbonds, const int *nUsed, double prefactor, double *ql) {
  if ((orderL != 4 && orderL != 6) || nbonds < 0 || nUsed == nullptr || ql == nullptr) {
    return false;
  }
#if defined(SEAMS_YLM_AVX512)
  if (seams::cpu::features().avx512dq) {
    return ylmLocalQlAvx(orderL, dx, dy, dz, nbonds, nUsed, prefactor, ql);
  }
#endif
  const int nComp = 2 * orderL + 1;
  for (int a = 0; a < 8; ++a) {
    double accRe[13] = {};
    double accIm[13] = {};
    double re[13], im[13];
    const int used = nUsed[a];
    for (int bond = 0; bond < used && bond < nbonds; ++bond) {
      const std::size_t at = static_cast<std::size_t>(bond) * 8u + static_cast<std::size_t>(a);
      if (!ylmReal(orderL, dx[at], dy[at], dz[at], re, im)) {
        continue;
      }
      for (int m = 0; m < nComp; ++m) {
        accRe[m] += re[m];
        accIm[m] += im[m];
      }
    }
    if (used <= 0) {
      ql[a] = 0.0;
      continue;
    }
    const double inv = 1.0 / static_cast<double>(used);
    double sum = 0.0;
    for (int m = 0; m < nComp; ++m) {
      const double r = accRe[m] * inv;
      const double i = accIm[m] * inv;
      sum += r * r + i * i;
    }
    ql[a] = std::sqrt(prefactor * sum);
  }
  return true;
}

void ylmBatch(int orderL, const double *dx, const double *dy, const double *dz,
              int n, double *re, double *im) {
  if (n <= 0 || (orderL != 3 && orderL != 4 && orderL != 6)) {
    return;
  }
#if defined(SEAMS_YLM_AVX512)
  if (seams::cpu::features().avx512dq) {
    ylmBatchAvx512(orderL, dx, dy, dz, n, re, im);
    return;
  }
#endif
  ylmBatchScalar(orderL, dx, dy, dz, n, re, im);
}

} // namespace ylm
} // namespace seams
