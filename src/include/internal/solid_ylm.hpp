#ifndef SEAMS_SOLID_YLM_H_
#define SEAMS_SOLID_YLM_H_

// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT

// Solid harmonics for l = 3, 4 and 6. A bond (x, y, z) of squared length
// r2 contributes Y_lm = S_lm(x, y, z) / r^l, and S_lm is a polynomial.
// Even l divides by a power of r2, so the harmonic has no square root.
// l = 3 keeps one factor 1/(r2 * sqrt(r2)).
//
// phi = atan2(x, y), so (y + i x)^m carries sin(theta)^m exp(i m phi).
// The caller passes inv = 1/r^l. A zero bond passes inv = 0 and a zero
// vector; the polynomials then contribute nothing.

namespace seams {
namespace solid {

template <class T>
inline void eval(int orderL, T x, T y, T z, T r2, T inv, const double *norm,
                 T *re, T *im) {
  T q0 = T(0);
  T q1 = T(0);
  T q2 = T(0);
  T q3 = T(0);
  T q4 = T(0);
  T q5 = T(0);
  T q6 = T(0);
  const T z2 = z * z;
  if (orderL == 3) {
    q0 = z * (T(5) * z2 - T(3) * r2);
    q1 = T(5) * z2 - r2;
    q2 = z;
    q3 = T(1);
  } else if (orderL == 4) {
    const T z3 = z2 * z;
    const T z4 = z2 * z2;
    const T r4 = r2 * r2;
    q0 = T(35) * z4 - T(30) * z2 * r2 + T(3) * r4;
    q1 = T(7) * z3 - T(3) * z * r2;
    q2 = T(7) * z2 - r2;
    q3 = z;
    q4 = T(1);
  } else {
    const T z3 = z2 * z;
    const T z4 = z2 * z2;
    const T z5 = z4 * z;
    const T z6 = z4 * z2;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    q0 = T(231) * z6 - T(315) * z4 * r2 + T(105) * z2 * r4 - T(5) * r6;
    q1 = T(33) * z5 - T(30) * z3 * r2 + T(5) * z * r4;
    q2 = T(33) * z4 - T(18) * z2 * r2 + r4;
    q3 = T(11) * z3 - T(3) * z * r2;
    q4 = T(11) * z2 - r2;
    q5 = z;
    q6 = T(1);
  }

  T q[7];
  q[0] = q0;
  q[1] = q1;
  q[2] = q2;
  q[3] = q3;
  if (orderL >= 4) {
    q[4] = q4;
  }
  if (orderL == 6) {
    q[5] = q5;
    q[6] = q6;
  }

  T pr = T(1);
  T pi = T(0);
  const T a0 = T(norm[0]) * q[0] * inv;
  re[orderL] = a0;
  im[orderL] = T(0);
  for (int m = 1; m <= orderL; ++m) {
    const T npr = pr * y - pi * x;
    const T npi = pr * x + pi * y;
    pr = npr;
    pi = npi;
    const T amp = T(norm[m]) * q[m] * inv;
    const int ineg = orderL - m;
    const int ipos = orderL + m;
    re[ineg] = amp * pr;
    im[ineg] = -amp * pi;
    const T sign = (m & 1) ? T(-1) : T(1);
    re[ipos] = sign * amp * pr;
    im[ipos] = sign * amp * pi;
  }
}

} // namespace solid
} // namespace seams

#endif
