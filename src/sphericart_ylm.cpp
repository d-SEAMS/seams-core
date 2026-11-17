//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <sphericart_ylm.hpp>

#ifdef SEAMS_HAS_SPHERICART
#include <sphericart.hpp>
#include <vector>
#endif

namespace seams {
namespace sphericart_ylm {

bool available() {
#ifdef SEAMS_HAS_SPHERICART
  return true;
#else
  return false;
#endif
}

int ylmCartesian(int orderL, const double *xyz, int nVec, double *ylmOut) {
  if (nVec == 0) {
    return 0;
  }
#ifndef SEAMS_HAS_SPHERICART
  (void)orderL;
  (void)xyz;
  (void)ylmOut;
  return 1;
#else
  if ((orderL != 3 && orderL != 4 && orderL != 6 && orderL != 8 &&
       orderL != 12) ||
      nVec < 0) {
    return 1;
  }
  const size_t l = static_cast<size_t>(orderL);
  const size_t rows = static_cast<size_t>(nVec);
  const size_t nComp = 2 * l + 1;
  std::vector<double> cart(xyz, xyz + rows * 3);
  // The calculator owns mutable scratch sized for the calling OpenMP team.
  // Frame workers and calls with different team sizes need independent storage.
  sphericart::SphericalHarmonics<double> calc(l);
  std::vector<double> sph;
  calc.compute(cart, sph);
  const size_t base = l * l;
  for (size_t i = 0; i < rows; i++) {
    const double *row = sph.data() + i * (l + 1) * (l + 1);
    double *out = ylmOut + i * nComp * 2;
    for (size_t k = 0; k < nComp; k++) {
      out[2 * k] = row[base + k];
      out[2 * k + 1] = 0.0;
    }
  }
  return 0;
#endif
}

} // namespace sphericart_ylm
} // namespace seams
