//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

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

void ylmBatch(int orderL, const double *dx, const double *dy, const double *dz,
              int n, double *re, double *im) {
  if (n <= 0 || (orderL != 3 && orderL != 4 && orderL != 6)) {
    return;
  }
#if defined(SEAMS_YLM_AVX512)
  if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq")) {
    ylmBatchAvx512(orderL, dx, dy, dz, n, re, im);
    return;
  }
#endif
  ylmBatchScalar(orderL, dx, dy, dz, n, re, im);
}

} // namespace ylm
} // namespace seams
