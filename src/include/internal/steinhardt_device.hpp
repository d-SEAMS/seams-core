#ifndef SEAMS_STEINHARDT_DEVICE_H_
#define SEAMS_STEINHARDT_DEVICE_H_

// Device-safe Steinhardt pieces: no STL containers, no function-local
// statics, no hash maps. Host OpenMP, MPI ranks, and OpenMP target
// regions all call the same functions. The Y_lm pair matches
// Sph::harmonicPair: Y_{l,m} = (-1)^m Y_{l,-m}^* from one amplitude
// (Steinhardt, Nelson and Ronchetti, Phys. Rev. B 28, 784 (1983);
//  the averaged qlBar is Lechner and Dellago, J. Chem. Phys. 129,
//  114707 (2008)).

#include <cmath>

#ifdef SEAMS_HAS_OFFLOAD
#pragma omp declare target
#endif

namespace seams {
namespace steinhardt {

inline void minImage(double &dx, double &dy, double &dz, double bx, double by,
                     double bz) {
  if (dx < -0.5 * bx) {
    dx += bx;
  }
  if (dx >= 0.5 * bx) {
    dx -= bx;
  }
  if (dy < -0.5 * by) {
    dy += by;
  }
  if (dy >= 0.5 * by) {
    dy -= by;
  }
  if (dz < -0.5 * bz) {
    dz += bz;
  }
  if (dz >= 0.5 * bz) {
    dz -= bz;
  }
}

inline double normLM(int orderL, int absM) {
  constexpr double pi = 3.14159265358979323846;
  if (orderL == 3) {
    switch (absM) {
    case 0:
      return 0.25 * std::sqrt(7.0 / pi);
    case 1:
      return 0.125 * std::sqrt(21.0 / pi);
    case 2:
      return 0.25 * std::sqrt(105.0 / (2.0 * pi));
    case 3:
      return 0.125 * std::sqrt(35.0 / pi);
    default:
      return 0.0;
    }
  }
  if (orderL == 4) {
    switch (absM) {
    case 0:
      return 0.1875 * std::sqrt(1.0 / pi);
    case 1:
      return 0.375 * std::sqrt(5.0 / pi);
    case 2:
      return 0.375 * std::sqrt(5.0 / (2.0 * pi));
    case 3:
      return 0.375 * std::sqrt(35.0 / pi);
    case 4:
      return 0.1875 * std::sqrt(35.0 / (2.0 * pi));
    default:
      return 0.0;
    }
  }
  if (orderL == 6) {
    switch (absM) {
    case 0:
      return 0.03125 * std::sqrt(13.0 / pi);
    case 1:
      return 0.0625 * std::sqrt(273.0 / (2.0 * pi));
    case 2:
      return 0.015625 * std::sqrt(1365.0 / pi);
    case 3:
      return 0.03125 * std::sqrt(1365.0 / pi);
    case 4:
      return 0.09375 * std::sqrt(91.0 / (2.0 * pi));
    case 5:
      return 0.09375 * std::sqrt(1001.0 / pi);
    case 6:
      return 0.015625 * std::sqrt(3003.0 / pi);
    default:
      return 0.0;
    }
  }
  return 0.0;
}

inline double legendreAmp(int orderL, int absM, const double *s,
                          const double *c) {
  if (orderL == 3) {
    switch (absM) {
    case 0:
      return 5.0 * c[3] - 3.0 * c[1];
    case 1:
      return s[1] * (5.0 * c[2] - 1.0);
    case 2:
      return s[2] * c[1];
    case 3:
      return s[3];
    default:
      return 0.0;
    }
  }
  if (orderL == 4) {
    switch (absM) {
    case 0:
      return 35.0 * c[4] - 30.0 * c[2] + 3.0;
    case 1:
      return s[1] * (7.0 * c[3] - 3.0 * c[1]);
    case 2:
      return s[2] * (7.0 * c[2] - 1.0);
    case 3:
      return s[3] * c[1];
    case 4:
      return s[4];
    default:
      return 0.0;
    }
  }
  if (orderL == 6) {
    switch (absM) {
    case 0:
      return 231.0 * c[6] - 315.0 * c[4] + 105.0 * c[2] - 5.0;
    case 1:
      return s[1] * (33.0 * c[5] - 30.0 * c[3] + 5.0 * c[1]);
    case 2:
      return s[2] * (33.0 * c[4] - 18.0 * c[2] + 1.0);
    case 3:
      return s[3] * (11.0 * c[3] - 3.0 * c[1]);
    case 4:
      return s[4] * (11.0 * c[2] - 1.0);
    case 5:
      return s[5] * c[1];
    case 6:
      return s[6];
    default:
      return 0.0;
    }
  }
  return 0.0;
}

// Writes 2*orderL+1 complex Y_lm as interleaved re,im starting at out.
inline void ylmAll(int orderL, double theta, double phi, double *out) {
  const int nComp = 2 * orderL + 1;
  for (int k = 0; k < nComp; k++) {
    out[2 * k] = 0.0;
    out[2 * k + 1] = 0.0;
  }

  double s[7];
  double c[7];
  double pr[7];
  double pi[7];
  const double sinT = std::sin(theta);
  const double cosT = std::cos(theta);
  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);
  s[0] = 1.0;
  c[0] = 1.0;
  pr[0] = 1.0;
  pi[0] = 0.0;
  for (int k = 1; k <= orderL; k++) {
    s[k] = s[k - 1] * sinT;
    c[k] = c[k - 1] * cosT;
    pr[k] = pr[k - 1] * cphi - pi[k - 1] * sphi;
    pi[k] = pr[k - 1] * sphi + pi[k - 1] * cphi;
  }

  for (int absM = 0; absM <= orderL; absM++) {
    const double amp = normLM(orderL, absM) * legendreAmp(orderL, absM, s, c);
    // Y_{l,-m} = amp * conj(e^{i m phi})
    const double nre = amp * pr[absM];
    const double nim = -amp * pi[absM];
    const int ineg = orderL - absM;
    out[2 * ineg] = nre;
    out[2 * ineg + 1] = nim;
    // Y_{l,+m} = (-1)^m amp * e^{i m phi}
    const double sign = (absM % 2 == 0) ? 1.0 : -1.0;
    const int ipos = orderL + absM;
    out[2 * ipos] = sign * amp * pr[absM];
    out[2 * ipos + 1] = sign * amp * pi[absM];
  }
}

// One particle, first pass: average Y_lm over CSR neighbours.
inline void qlmOneAtom(int iatom, int orderL, const double *xyz,
                       const int *offsets, const int *cols, double bx,
                       double by, double bz, double *qlmInterleaved) {
  const int nComp = 2 * orderL + 1;
  const int iOff = 3 * iatom;
  const int row = iatom * nComp;
  for (int m = 0; m < nComp; m++) {
    qlmInterleaved[2 * (row + m)] = 0.0;
    qlmInterleaved[2 * (row + m) + 1] = 0.0;
  }
  const int j0 = offsets[iatom];
  const int j1 = offsets[iatom + 1];
  int nUsed = 0;
  double ylm[26];
  for (int p = j0; p < j1; p++) {
    const int jatom = cols[p];
    double dx = xyz[iOff] - xyz[3 * jatom];
    double dy = xyz[iOff + 1] - xyz[3 * jatom + 1];
    double dz = xyz[iOff + 2] - xyz[3 * jatom + 2];
    minImage(dx, dy, dz, bx, by, bz);
    const double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 == 0.0) {
      continue;
    }
    const double r = std::sqrt(r2);
    const double phi = std::atan2(dx, dy);
    const double theta = std::acos(dz / r);
    ylmAll(orderL, theta, phi, ylm);
    for (int m = 0; m < nComp; m++) {
      qlmInterleaved[2 * (row + m)] += ylm[2 * m];
      qlmInterleaved[2 * (row + m) + 1] += ylm[2 * m + 1];
    }
    nUsed++;
  }
  if (nUsed == 0) {
    return;
  }
  const double inv = 1.0 / static_cast<double>(nUsed);
  for (int m = 0; m < nComp; m++) {
    qlmInterleaved[2 * (row + m)] *= inv;
    qlmInterleaved[2 * (row + m) + 1] *= inv;
  }
}

inline void qlOneAtom(int iatom, int orderL, const double *qlmInterleaved,
                      const int *offsets, const int *cols, double *ql,
                      double *qlBar) {
  constexpr double pi = 3.14159265358979323846;
  const int nComp = 2 * orderL + 1;
  const double prefactor = 4.0 * pi / (2.0 * orderL + 1.0);
  const int row = iatom * nComp;
  double sumLocal = 0.0;
  for (int m = 0; m < nComp; m++) {
    const double re = qlmInterleaved[2 * (row + m)];
    const double im = qlmInterleaved[2 * (row + m) + 1];
    sumLocal += re * re + im * im;
  }
  ql[iatom] = std::sqrt(prefactor * sumLocal);

  const int j0 = offsets[iatom];
  const int j1 = offsets[iatom + 1];
  if (j0 == j1) {
    qlBar[iatom] = ql[iatom];
    return;
  }
  double barRe[13];
  double barIm[13];
  for (int m = 0; m < nComp; m++) {
    barRe[m] = qlmInterleaved[2 * (row + m)];
    barIm[m] = qlmInterleaved[2 * (row + m) + 1];
  }
  int nContrib = 1;
  for (int p = j0; p < j1; p++) {
    const int jatom = cols[p];
    const int jRow = jatom * nComp;
    for (int m = 0; m < nComp; m++) {
      barRe[m] += qlmInterleaved[2 * (jRow + m)];
      barIm[m] += qlmInterleaved[2 * (jRow + m) + 1];
    }
    nContrib++;
  }
  const double inv = 1.0 / static_cast<double>(nContrib);
  double sumBar = 0.0;
  for (int m = 0; m < nComp; m++) {
    const double re = barRe[m] * inv;
    const double im = barIm[m] * inv;
    sumBar += re * re + im * im;
  }
  qlBar[iatom] = std::sqrt(prefactor * sumBar);
}

} // namespace steinhardt
} // namespace seams

#ifdef SEAMS_HAS_OFFLOAD
#pragma omp end declare target
#endif

#endif
