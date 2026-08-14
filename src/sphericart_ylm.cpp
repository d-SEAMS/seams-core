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
  if ((orderL != 3 && orderL != 4 && orderL != 6) || nVec < 0) {
    return 1;
  }
  const int nComp = 2 * orderL + 1;
  std::vector<double> cart(static_cast<size_t>(nVec) * 3);
  for (int i = 0; i < nVec * 3; i++) {
    cart[static_cast<size_t>(i)] = xyz[i];
  }
  static sphericart::SphericalHarmonics<double> calc3(3);
  static sphericart::SphericalHarmonics<double> calc4(4);
  static sphericart::SphericalHarmonics<double> calc6(6);
  sphericart::SphericalHarmonics<double> *calc = &calc6;
  if (orderL == 3) {
    calc = &calc3;
  } else if (orderL == 4) {
    calc = &calc4;
  }
  std::vector<double> sph;
  calc->compute(cart, sph);
  const int base = orderL * orderL;
  for (int i = 0; i < nVec; i++) {
    const double *row = sph.data() + static_cast<size_t>(i) * (orderL + 1) * (orderL + 1);
    double *out = ylmOut + static_cast<size_t>(i) * nComp * 2;
    for (int k = 0; k < nComp; k++) {
      out[2 * k] = row[base + k];
      out[2 * k + 1] = 0.0;
    }
  }
  return 0;
#endif
}

} // namespace sphericart_ylm
} // namespace seams
