#ifndef SEAMS_YLM_CARTESIAN_H_
#define SEAMS_YLM_CARTESIAN_H_

// Closed-form complex Y_lm for l = 3 and l = 6 as real polynomials in the
// unit-vector components. The phase convention matches harmonicPair and
// ylmAllTrig: phi = atan2(dx, dy), so e^{i phi} = (y + i x) / rho and
// sin(theta)^|m| e^{i m phi} = (y + i x)^|m| with no inverse sine or
// cosine. Norms are the Condon-Shortley constants used by normQ3/normQ6.

#include <cmath>

namespace seams {
namespace ylm {

inline constexpr double N3[4] = {
    0x1.7e21efc0768b0p-2, // 0.25 * sqrt(7/pi)
    0x1.4aefbf05f1305p-2, // 0.125 * sqrt(21/pi)
    0x1.05a0d7164496dp+0, // 0.25 * sqrt(105/(2 pi))
    0x1.ab3cb8d675a30p-2  // 0.125 * sqrt(35/pi)
};

inline constexpr double N6[7] = {
    0x1.046123cc5b3b1p-4,  // 0.03125 * sqrt(13/pi)
    0x1.a5dce8e510806p-2,  // 0.0625 * sqrt(273/(2 pi))
    0x1.4d830dc3746c6p-2,  // 0.015625 * sqrt(1365/pi)
    0x1.4d830dc3746c6p-1,  // 0.03125 * sqrt(1365/pi)
    0x1.6d58113f224e2p-2,  // 0.09375 * sqrt(91/(2 pi))
    0x1.ac67615af3632p+0,  // 0.09375 * sqrt(1001/pi)
    0x1.eead9a142c943p-2   // 0.015625 * sqrt(3003/pi)
};

// Associated-Legendre amplitude without the azimuthal factor, in the
// hard-coded table convention (the sin^|m| factor is carried by w^m).
template <class T>
inline T legendrePoly(int orderL, int absM, const T *zpow) {
  if (orderL == 3) {
    switch (absM) {
    case 0:
      return T(5) * zpow[3] - T(3) * zpow[1];
    case 1:
      return T(5) * zpow[2] - T(1);
    case 2:
      return zpow[1];
    default:
      return T(1);
    }
  }
  switch (absM) {
  case 0:
    return T(231) * zpow[6] - T(315) * zpow[4] + T(105) * zpow[2] - T(5);
  case 1:
    return T(33) * zpow[5] - T(30) * zpow[3] + T(5) * zpow[1];
  case 2:
    return T(33) * zpow[4] - T(18) * zpow[2] + T(1);
  case 3:
    return T(11) * zpow[3] - T(3) * zpow[1];
  case 4:
    return T(11) * zpow[2] - T(1);
  case 5:
    return zpow[1];
  default:
    return T(1);
  }
}

// x, y, z are components of a unit vector. re/im hold 2*l+1 values,
// m = -l .. +l. l is 3 or 6.
template <class T>
inline void evalUnit(int orderL, T x, T y, T z, T *re, T *im) {
  T wre[7];
  T wim[7];
  T zpow[7];
  wre[0] = T(1);
  wim[0] = T(0);
  zpow[0] = T(1);
  for (int m = 1; m <= orderL; ++m) {
    const T prevRe = wre[m - 1];
    const T prevIm = wim[m - 1];
    wre[m] = prevRe * y - prevIm * x;
    wim[m] = prevRe * x + prevIm * y;
    zpow[m] = zpow[m - 1] * z;
  }
  const double *norm = (orderL == 3) ? N3 : N6;
  for (int m = 0; m <= orderL; ++m) {
    const T amp = T(norm[m]) * legendrePoly(orderL, m, zpow);
    const int ineg = orderL - m;
    const int ipos = orderL + m;
    re[ineg] = amp * wre[m];
    im[ineg] = -amp * wim[m];
    const T signedAmp = (m & 1) ? -amp : amp;
    re[ipos] = signedAmp * wre[m];
    im[ipos] = signedAmp * wim[m];
  }
}

// Displacement, not necessarily a unit vector. A zero vector writes
// zeros and returns false so the caller can drop the bond.
inline bool ylmReal(int orderL, double dx, double dy, double dz, double *re,
                    double *im) {
  const int nComp = 2 * orderL + 1;
  const double r2 = dx * dx + dy * dy + dz * dz;
  if (!(r2 > 0.0) || (orderL != 3 && orderL != 6)) {
    for (int m = 0; m < nComp; ++m) {
      re[m] = 0.0;
      im[m] = 0.0;
    }
    return false;
  }
  const double inv = 1.0 / std::sqrt(r2);
  evalUnit(orderL, dx * inv, dy * inv, dz * inv, re, im);
  return true;
}

} // namespace ylm
} // namespace seams

#endif
