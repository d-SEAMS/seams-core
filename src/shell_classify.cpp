//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Per-atom ice label. The neighbour stencil is the orthorhombic cell list
// used by the rest of d-SEAMS: cells no wider than the 5.6 A shell, and the
// 27 cells around each atom. A skewed or tiny box keeps the finer 2.8 A
// grid. The quaternion is the largest eigenvector of the 4x4 matrix K from
// the hand-rolled Jacobi solver below. Eigen is not included in this file,
// so the translation unit can be built with AVX2.

#include <shell_hot.hpp>

#include <algorithm>
#include <cmath>
#include <new>
#include <cstddef>
#include <memory>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace symm {
namespace {

constexpr double kBond = 3.5;
constexpr double kWide = 4.05;
constexpr double kCentre = 5.6;
constexpr double kCentreHex = 6.2;
constexpr double kGrid = 2.8;
constexpr int kMaxBond = 16;
constexpr int kMaxWide = 4;

struct Bond {
  int j;
  double dx;
  double dy;
  double dz;
};

// Uninitialized storage. new T[] value-initializes and the serial zero of the
// bond table is larger than the parallel fill that overwrites it.
template <class T>
struct RawBuf {
  T *p;
  explicit RawBuf(size_t n) : p(static_cast<T *>(::operator new(n * sizeof(T)))) {}
  ~RawBuf() { ::operator delete(p); }
  RawBuf(const RawBuf &) = delete;
  RawBuf &operator=(const RawBuf &) = delete;
  T *get() const { return p; }
  T &operator[](size_t i) const { return p[i]; }
};
constexpr double kCleanRms = 0.02;
constexpr double kCleanShift = 0.05;

struct Image {
  int cell;
  int sx;
  int sy;
  int sz;
};

struct Packed {
  double x, y, z;
  int id;
};

double dist2(const double a[3], const double b[3]) {
  const double dx = a[0] - b[0];
  const double dy = a[1] - b[1];
  const double dz = a[2] - b[2];
  return dx * dx + dy * dy + dz * dz;
}

void invert3(const double M[3][3], double inv[3][3]) {
  const double a = M[0][0], b = M[0][1], c = M[0][2];
  const double d = M[1][0], e = M[1][1], f = M[1][2];
  const double g = M[2][0], h = M[2][1], i = M[2][2];
  const double det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
  const double s = 1.0 / det;
  inv[0][0] = s * (e * i - f * h);
  inv[0][1] = s * (c * h - b * i);
  inv[0][2] = s * (b * f - c * e);
  inv[1][0] = s * (f * g - d * i);
  inv[1][1] = s * (a * i - c * g);
  inv[1][2] = s * (c * d - a * f);
  inv[2][0] = s * (d * h - e * g);
  inv[2][1] = s * (b * g - a * h);
  inv[2][2] = s * (a * e - b * d);
}

// Largest eigenvector of a 4x4 symmetric matrix. One Jacobi sweep is a
// closed rotation; a handful of sweeps is the quaternion form of Kabsch
// (Coutsias, Seok, Dill 2004; Theobald 2005) without a 3x3 SVD.
void jacobi4(double A[4][4], double Q[4][4], int sweeps) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      Q[i][j] = i == j ? 1.0 : 0.0;
    }
  }
  for (int sweep = 0; sweep < sweeps; sweep++) {
    for (int p = 0; p < 3; p++) {
      for (int q = p + 1; q < 4; q++) {
        const double apq = A[p][q];
        if (std::abs(apq) < 1e-14) {
          continue;
        }
        const double app = A[p][p];
        const double aqq = A[q][q];
        const double tau = (aqq - app) / (2.0 * apq);
        const double t = std::copysign(1.0, tau) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double c = 1.0 / std::sqrt(1.0 + t * t);
        const double s = t * c;
        A[p][p] = app - t * apq;
        A[q][q] = aqq + t * apq;
        A[p][q] = 0.0;
        A[q][p] = 0.0;
        for (int r = 0; r < 4; r++) {
          if (r == p || r == q) {
            continue;
          }
          const double arp = A[r][p];
          const double arq = A[r][q];
          A[r][p] = A[p][r] = c * arp - s * arq;
          A[r][q] = A[q][r] = s * arp + c * arq;
        }
        for (int r = 0; r < 4; r++) {
          const double qp = Q[r][p];
          const double qq = Q[r][q];
          Q[r][p] = c * qp - s * qq;
          Q[r][q] = s * qp + c * qq;
        }
      }
    }
  }
}

void quatToR(const double q[4], double R[3][3]) {
  const double w = q[0];
  const double x = q[1];
  const double y = q[2];
  const double z = q[3];
  R[0][0] = 1.0 - 2.0 * (y * y + z * z);
  R[0][1] = 2.0 * (x * y - z * w);
  R[0][2] = 2.0 * (x * z + y * w);
  R[1][0] = 2.0 * (x * y + z * w);
  R[1][1] = 1.0 - 2.0 * (x * x + z * z);
  R[1][2] = 2.0 * (y * z - x * w);
  R[2][0] = 2.0 * (x * z - y * w);
  R[2][1] = 2.0 * (y * z + x * w);
  R[2][2] = 1.0 - 2.0 * (x * x + y * y);
}

void kabsch17(const Shell17 &fixed, const Shell17 &mobile, int n, double R[3][3], double &rmsd,
              double *scaleFree) {
  double cf[3] = {0, 0, 0};
  double cm[3] = {0, 0, 0};
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < 3; k++) {
      cf[k] += fixed.p[i][k];
      cm[k] += mobile.p[i][k];
    }
  }
  const double inv = 1.0 / static_cast<double>(std::max(n, 1));
  for (int k = 0; k < 3; k++) {
    cf[k] *= inv;
    cm[k] *= inv;
  }
  double H[3][3] = {};
  for (int i = 0; i < n; i++) {
    const double p[3] = {fixed.p[i][0] - cf[0], fixed.p[i][1] - cf[1], fixed.p[i][2] - cf[2]};
    const double q[3] = {mobile.p[i][0] - cm[0], mobile.p[i][1] - cm[1], mobile.p[i][2] - cm[2]};
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        H[a][b] += p[a] * q[b];
      }
    }
  }
  double K[4][4];
  K[0][0] = H[0][0] + H[1][1] + H[2][2];
  K[0][1] = H[1][2] - H[2][1];
  K[0][2] = H[2][0] - H[0][2];
  K[0][3] = H[0][1] - H[1][0];
  K[1][1] = H[0][0] - H[1][1] - H[2][2];
  K[1][2] = H[0][1] + H[1][0];
  K[1][3] = H[2][0] + H[0][2];
  K[2][2] = -H[0][0] + H[1][1] - H[2][2];
  K[2][3] = H[1][2] + H[2][1];
  K[3][3] = -H[0][0] - H[1][1] + H[2][2];
  K[1][0] = K[0][1];
  K[2][0] = K[0][2];
  K[2][1] = K[1][2];
  K[3][0] = K[0][3];
  K[3][1] = K[1][3];
  K[3][2] = K[2][3];
  double Q[4][4];
  jacobi4(K, Q, n <= 4 ? 2 : 4);
  int best = 0;
  for (int i = 1; i < 4; i++) {
    if (K[i][i] > K[best][best]) {
      best = i;
    }
  }
  double q[4] = {Q[0][best], -Q[1][best], -Q[2][best], -Q[3][best]};
  const double qn = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  if (qn > 1e-15) {
    for (int i = 0; i < 4; i++) {
      q[i] /= qn;
    }
  }
  quatToR(q, R);
  double gp = 0.0;
  double gq = 0.0;
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < 3; k++) {
      const double dp = fixed.p[i][k] - cf[k];
      const double dq = mobile.p[i][k] - cm[k];
      gp += dp * dp;
      gq += dq * dq;
    }
  }
  const double lam = K[best][best];
  const double acc = gp + gq - 2.0 * lam;
  rmsd = std::sqrt(std::max(0.0, acc) / static_cast<double>(n));
  if (scaleFree != nullptr) {
    double frac = 1.0;
    if (gp > 1e-18 && gq > 1e-18) {
      frac = std::sqrt(std::max(0.0, 1.0 - (lam * lam) / (gp * gq)));
    }
    *scaleFree = frac;
  }
}

void kabsch4(const double fixed[][3], const double mobile[][3], double R[3][3], double &rmsd) {
  Shell17 a{};
  Shell17 b{};
  for (int i = 0; i < 4; i++) {
    for (int k = 0; k < 3; k++) {
      a.p[i][k] = fixed[i][k];
      b.p[i][k] = mobile[i][k];
    }
  }
  kabsch17(a, b, 4, R, rmsd, nullptr);
}

double tetraChirality(const double v[4][3]) {
  const double a[3] = {v[1][0] - v[0][0], v[1][1] - v[0][1], v[1][2] - v[0][2]};
  const double b[3] = {v[2][0] - v[0][0], v[2][1] - v[0][1], v[2][2] - v[0][2]};
  const double c[3] = {v[3][0] - v[0][0], v[3][1] - v[0][1], v[3][2] - v[0][2]};
  const double cx = b[1] * c[2] - b[2] * c[1];
  const double cy = b[2] * c[0] - b[0] * c[2];
  const double cz = b[0] * c[1] - b[1] * c[0];
  return a[0] * cx + a[1] * cy + a[2] * cz;
}

double bondAngle(const double v[3], const double e1[3], const double e2[3]) {
  return std::atan2(v[0] * e2[0] + v[1] * e2[1] + v[2] * e2[2],
                    v[0] * e1[0] + v[1] * e1[1] + v[2] * e1[2]);
}

void bondFrame(const double axis[3], const double ref[3], double e1[3], double e2[3]) {
  double an = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  double u[3] = {axis[0] / an, axis[1] / an, axis[2] / an};
  const double dot = ref[0] * u[0] + ref[1] * u[1] + ref[2] * u[2];
  e1[0] = ref[0] - dot * u[0];
  e1[1] = ref[1] - dot * u[1];
  e1[2] = ref[2] - dot * u[2];
  const double e1n = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  if (e1n > 1e-8) {
    e1[0] /= e1n;
    e1[1] /= e1n;
    e1[2] /= e1n;
  }
  e2[0] = u[1] * e1[2] - u[2] * e1[1];
  e2[1] = u[2] * e1[0] - u[0] * e1[2];
  e2[2] = u[0] * e1[1] - u[1] * e1[0];
}

// True when at least one bond has a staggered second neighbour (cosine at
// most 0.75). Out of line and cold so the ideal-ice fit keeps its layout.
bool __attribute__((noinline, cold)) hasStaggeredBond(const double envF[4][3], const double envS[12][3],
                                                      int nS) {
  for (int a = 0; a < 4; a++) {
    double u[3] = {envF[a][0], envF[a][1], envF[a][2]};
    const double un = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    if (!(un > 1e-8)) {
      return false;
    }
    u[0] /= un;
    u[1] /= un;
    u[2] /= un;
    double mx = -1.0;
    for (int b = 0; b < 4 && mx < 0.75; b++) {
      if (b == a) {
        continue;
      }
      double ref[3] = {envF[b][0], envF[b][1], envF[b][2]};
      const double rd = ref[0] * u[0] + ref[1] * u[1] + ref[2] * u[2];
      ref[0] -= rd * u[0];
      ref[1] -= rd * u[1];
      ref[2] -= rd * u[2];
      const double rn2 = ref[0] * ref[0] + ref[1] * ref[1] + ref[2] * ref[2];
      if (rn2 < 1e-12) {
        continue;
      }
      const double invR = 1.0 / std::sqrt(rn2);
      for (int e = 0; e < nS; e++) {
        double v[3] = {envS[e][0] - envF[a][0], envS[e][1] - envF[a][1], envS[e][2] - envF[a][2]};
        const double bond = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        if (bond < 2.35 * 2.35 || bond > kBond * kBond) {
          continue;
        }
        const double vd = v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
        v[0] -= vd * u[0];
        v[1] -= vd * u[1];
        v[2] -= vd * u[2];
        const double vn2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        if (vn2 < 1e-12) {
          continue;
        }
        const double c = (ref[0] * v[0] + ref[1] * v[1] + ref[2] * v[2]) * invR / std::sqrt(vn2);
        if (c > mx) {
          mx = c;
        }
      }
    }
    if (mx <= 0.75) {
      return true;
    }
  }
  return false;
}

// Second-shell diamond test on the four neighbours, the same split OVITO's
// identifier and pyscal3 use: the twelve second neighbours of a cubic
// diamond are an FCC shell (every bond staggered) and those of a hexagonal
// diamond are an HCP shell (one bond eclipsed). A clean shell is a single
// radius. Anything wider falls through to the full fit.
int diamondSecondShell(const double first[4][3], const double *sx, const double *sy, const double *sz) {
  double rmin = 1e300;
  double rmax = 0.0;
  for (int i = 0; i < 12; i++) {
    const double r2 = sx[i] * sx[i] + sy[i] * sy[i] + sz[i] * sz[i];
    if (r2 < rmin) {
      rmin = r2;
    }
    if (r2 > rmax) {
      rmax = r2;
    }
  }
  if (!(rmin > 1.0) || rmax > rmin * 1.0404) {
    return -1;
  }
  int nEclipsed = 0;
  int nStaggered = 0;
  for (int a = 0; a < 4; a++) {
    double u[3] = {first[a][0], first[a][1], first[a][2]};
    const double un2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    if (!(un2 > 1e-8)) {
      return -1;
    }
    const double un = std::sqrt(un2);
    u[0] /= un;
    u[1] /= un;
    u[2] /= un;
    double mx = -1.0;
    for (int b = 0; b < 4; b++) {
      if (b == a) {
        continue;
      }
      double ref[3] = {first[b][0], first[b][1], first[b][2]};
      const double rd = ref[0] * u[0] + ref[1] * u[1] + ref[2] * u[2];
      ref[0] -= rd * u[0];
      ref[1] -= rd * u[1];
      ref[2] -= rd * u[2];
      const double rn2 = ref[0] * ref[0] + ref[1] * ref[1] + ref[2] * ref[2];
      if (rn2 < 1e-12) {
        continue;
      }
      const double invR = 1.0 / std::sqrt(rn2);
      for (int e = 0; e < 3; e++) {
        const int s = 3 * a + e;
        double v[3] = {sx[s] - first[a][0], sy[s] - first[a][1], sz[s] - first[a][2]};
        const double vd = v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
        v[0] -= vd * u[0];
        v[1] -= vd * u[1];
        v[2] -= vd * u[2];
        const double vn2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        if (vn2 < 1e-12) {
          continue;
        }
        const double c = (ref[0] * v[0] + ref[1] * v[1] + ref[2] * v[2]) * invR / std::sqrt(vn2);
        if (c > mx) {
          mx = c;
        }
      }
    }
    if (mx >= 0.95) {
      nEclipsed++;
    } else if (mx <= 0.65) {
      nStaggered++;
    } else {
      return -1;
    }
  }
  if (nEclipsed == 0 && nStaggered == 4) {
    return 1;
  }
  if (nEclipsed == 1 && nStaggered == 3) {
    return 0;
  }
  return -1;
}

int shellFit(const double center[3], const double nb[][3], int nNb, const CanonShell &canon, double gate,
             bool recover = false) {
  const int m = std::min(nNb, 16);
  int order[16];
  double dd[16];
  for (int i = 0; i < m; i++) {
    order[i] = i;
    dd[i] = dist2(center, nb[i]);
  }
  for (int a = 0; a < 4 && a < m; a++) {
    int best = a;
    for (int b = a + 1; b < m; b++) {
      if (dd[b] < dd[best]) {
        best = b;
      }
    }
    std::swap(dd[a], dd[best]);
    std::swap(order[a], order[best]);
  }
  int nFirst = 0;
  int first[4];
  for (int a = 0; a < 4 && a < m; a++) {
    if (dd[a] >= 2.15 * 2.15 && dd[a] <= 3.85 * 3.85) {
      first[nFirst++] = order[a];
    }
  }
  if (nFirst < 4) {
    return -1;
  }
  double envF[4][3];
  bool used[16] = {};
  for (int k = 0; k < 4; k++) {
    used[first[k]] = true;
    for (int c = 0; c < 3; c++) {
      envF[k][c] = nb[first[k]][c] - center[c];
    }
  }
  double envS[12][3];
  int nS = 0;
  for (int i = 0; i < m && nS < 12; i++) {
    if (used[i]) {
      continue;
    }
    for (int c = 0; c < 3; c++) {
      envS[nS][c] = nb[i][c] - center[c];
    }
    nS++;
  }
  if (nS < 12) {
    return -1;
  }
  // Ideal ice has four equal bonds and skips this. A tetrahedron whose four
  // bonds are all eclipsed is the liquid pattern and does not enter the
  // permutation search.
  double lenMin = 1e300;
  double lenMax = 0.0;
  for (int k = 0; k < 4; k++) {
    const double len2 = envF[k][0] * envF[k][0] + envF[k][1] * envF[k][1] + envF[k][2] * envF[k][2];
    if (len2 < lenMin) {
      lenMin = len2;
    }
    if (len2 > lenMax) {
      lenMax = len2;
    }
  }
  // 1.02 squared. Ideal ice is under it and does not call the eclipsed scan.
  // An all-eclipsed shell next to a crystal is ice Ih, so that pass does not return here.
  const bool eclipsed = lenMax >= lenMin * 1.0404 && !hasStaggeredBond(envF, envS, nS);
  if (eclipsed && !recover) {
    return -1;
  }
  double best = 1e300;
  int bestK = -1;
  const int kinds[3] = {1, 0, 2};
  double closest = 1e300;
  for (int ko = 0; ko < 3 && !(best <= gate); ko++) {
    if (ko > 0 && closest > 40.0) {
      break;
    }
    const int kind = kinds[ko];
    const double templChi = tetraChirality(canon.first[kind]);
    int pidx[4] = {0, 1, 2, 3};
    int tried = 0;
    do {
      double fixed[4][3];
      double mobile[4][3];
      for (int k = 0; k < 4; k++) {
        for (int c = 0; c < 3; c++) {
          fixed[k][c] = canon.first[kind][k][c];
          mobile[k][c] = envF[pidx[k]][c];
        }
      }
      if (templChi * tetraChirality(mobile) < 0.0) {
        continue;
      }
      tried++;
      if (tried > 12 || (!recover && tried > 3 && closest > 40.0 && closest < 1.0e200)) {
        break;
      }
      double R[3][3];
      double rms = 0.0;
      kabsch4(fixed, mobile, R, rms);
      if (rms > 0.8) {
        continue;
      }
      double gT = 0.0;
      double gE = 0.0;
      for (int k = 0; k < 4; k++) {
        for (int c = 0; c < 3; c++) {
          gT += fixed[k][c] * fixed[k][c];
          gE += mobile[k][c] * mobile[k][c];
        }
      }
      const double shellScale = gT > 1e-18 ? std::sqrt(gE / gT) : 1.0;
      Shell17 fixed17{};
      Shell17 mobile17{};
      for (int k = 0; k < 4; k++) {
        for (int c = 0; c < 3; c++) {
          fixed17.p[k + 1][c] = canon.first[kind][k][c];
          mobile17.p[k + 1][c] = envF[pidx[k]][c];
        }
      }
      bool taken[16] = {};
      bool paired = true;
      int slot = 5;
      double shiftSum = 0.0;
      for (int a = 0; a < 4 && paired; a++) {
        const int envArm = pidx[a];
        int cand[3];
        double candD[3];
        int nC = 0;
        for (int e = 0; e < nS; e++) {
          if (taken[e]) {
            continue;
          }
          const double dx = envS[e][0] - envF[envArm][0];
          const double dy = envS[e][1] - envF[envArm][1];
          const double dz = envS[e][2] - envF[envArm][2];
          const double d2 = dx * dx + dy * dy + dz * dz;
          if (nC < 3) {
            cand[nC] = e;
            candD[nC] = d2;
            nC++;
          } else {
            int worst = 0;
            for (int c = 1; c < 3; c++) {
              if (candD[c] > candD[worst]) {
                worst = c;
              }
            }
            if (d2 < candD[worst]) {
              cand[worst] = e;
              candD[worst] = d2;
            }
          }
        }
        if (nC < 3 || canon.nArm[kind][a] < 3) {
          paired = false;
          break;
        }
        double e1[3];
        double e2[3];
        bondFrame(envF[envArm], envF[pidx[(a + 1) % 4]], e1, e2);
        double ang[3];
        for (int i = 0; i < 3; i++) {
          const double outer[3] = {envS[cand[i]][0] - envF[envArm][0],
                                   envS[cand[i]][1] - envF[envArm][1],
                                   envS[cand[i]][2] - envF[envArm][2]};
          ang[i] = bondAngle(outer, e1, e2);
        }
        for (int i = 0; i < 3; i++) {
          for (int j = i + 1; j < 3; j++) {
            if (ang[j] < ang[i]) {
              std::swap(ang[i], ang[j]);
              std::swap(cand[i], cand[j]);
            }
          }
        }
        int seq[2][3];
        for (int i = 0; i < 3; i++) {
          seq[0][i] = cand[i];
          seq[1][i] = cand[2 - i];
        }
        // The three template arms rotate once. The six pairings only permute them.
        double pred[3][3];
        for (int i = 0; i < 3; i++) {
          const int s = canon.arm[kind][a][i];
          const double q0 = canon.second[kind][s][0];
          const double q1 = canon.second[kind][s][1];
          const double q2 = canon.second[kind][s][2];
          pred[i][0] = shellScale * (R[0][0] * q0 + R[1][0] * q1 + R[2][0] * q2);
          pred[i][1] = shellScale * (R[0][1] * q0 + R[1][1] * q1 + R[2][1] * q2);
          pred[i][2] = shellScale * (R[0][2] * q0 + R[1][2] * q1 + R[2][2] * q2);
        }
        int bestDir = 0;
        int bestShift = 0;
        double bestShiftErr = 1e300;
        for (int dir = 0; dir < 2; dir++) {
          for (int shift = 0; shift < 3; shift++) {
            double err = 0.0;
            for (int i = 0; i < 3; i++) {
              const int e = seq[dir][i];
              const double dx = pred[(i + shift) % 3][0] - envS[e][0];
              const double dy = pred[(i + shift) % 3][1] - envS[e][1];
              const double dz = pred[(i + shift) % 3][2] - envS[e][2];
              err += dx * dx + dy * dy + dz * dz;
            }
            if (err < bestShiftErr) {
              bestShiftErr = err;
              bestShift = shift;
              bestDir = dir;
            }
          }
        }
        for (int i = 0; i < 3; i++) {
          const int s = canon.arm[kind][a][(i + bestShift) % 3];
          const int e = seq[bestDir][i];
          taken[e] = true;
          for (int c = 0; c < 3; c++) {
            fixed17.p[slot][c] = canon.second[kind][s][c];
            mobile17.p[slot][c] = envS[e][c];
          }
          slot++;
        }
        shiftSum += bestShiftErr;
      }
      if (paired && slot == 17 && shiftSum < closest) {
        closest = shiftSum;
      }
      if (!paired || slot != 17 || shiftSum > 20.0) {
        continue;
      }
      // A first shell that already sits on the template does not need the
      // 17-point scale fit. The arm error has named the same template.
      if (rms <= kCleanRms && shiftSum <= kCleanShift) {
        best = 0.0;
        bestK = kind;
        break;
      }
      double R17[3][3];
      double rms17 = 0.0;
      double frac = 1.0;
      kabsch17(fixed17, mobile17, 17, R17, rms17, &frac);
      (void)rms17;
      if (frac < best) {
        best = frac;
        bestK = kind;
      }
      if (best <= gate) {
        break;
      }
    } while (std::next_permutation(pidx, pidx + 4));
  }
  if (bestK < 0 || !(best <= gate)) {
    return -1;
  }
  return bestK;
}

} // namespace

int bondsBeyond(int atom, const Packed *packed, const int *cellStart, const Image *images, int nImage,
               const int *bin, const double lat[3][3], const double *px, const double *py, const double *pz,
               Bond *out);

void recoverHex(int n, const Bond *bonds, const int *nbon, const Packed *packed, const int *cellStart,
                const Image *images, int nImage, const int *bin, const double lat[3][3], const double *px,
                const double *py, const double *pz, int *kind, const CanonShell &canon, double gate);

void classifyPrepared(const double *x, const double *y, const double *z, int n,
                      const double lat[3][3], const CanonShell &canon, double gate, int *kind) {
  if (n <= 0) {
    return;
  }
  double invT[3][3];
  double M[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      M[i][j] = lat[j][i];
    }
  }
  invert3(M, invT);
  auto length = [](const double v[3]) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  };
  const double la = length(lat[0]);
  const double lb = length(lat[1]);
  const double lc = length(lat[2]);
  const double skew = std::abs(lat[1][0]) + std::abs(lat[2][0]) + std::abs(lat[2][1]);
  const bool ortho = skew < 1e-8;
  int nx = 0;
  int ny = 0;
  int nz = 0;
  bool fast = false;
  if (ortho && la >= 3.0 * kBond && lb >= 3.0 * kBond && lc >= 3.0 * kBond) {
    // The 3.5 A bond list. Cells are at most one bond wide, so the 27-cell
    // stencil holds every bonded neighbour. The second shell is read from
    // those neighbours' own rows, not from a 5.6 A search.
    nx = std::max(3, static_cast<int>(std::ceil(la / kBond)));
    ny = std::max(3, static_cast<int>(std::ceil(lb / kBond)));
    nz = std::max(3, static_cast<int>(std::ceil(lc / kBond)));
    fast = true;
  } else {
    nx = std::max(1, static_cast<int>(la / kGrid));
    ny = std::max(1, static_cast<int>(lb / kGrid));
    nz = std::max(1, static_cast<int>(lc / kGrid));
  }
  const int ncell = nx * ny * nz;
  const int reach = fast ? 1 : (skew > 1.0 ? 4 : 2);
  RawBuf<double> px(static_cast<size_t>(n));
  RawBuf<double> py(static_cast<size_t>(n));
  RawBuf<double> pz(static_cast<size_t>(n));
  RawBuf<int> bin(static_cast<size_t>(n));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 4096)
#endif
  for (int i = 0; i < n; i++) {
    const double c0 = x[i];
    const double c1 = y[i];
    const double c2 = z[i];
    double f[3] = {0, 0, 0};
    f[0] = invT[0][0] * c0 + invT[0][1] * c1 + invT[0][2] * c2;
    f[1] = invT[1][0] * c0 + invT[1][1] * c1 + invT[1][2] * c2;
    f[2] = invT[2][0] * c0 + invT[2][1] * c1 + invT[2][2] * c2;
    for (int a = 0; a < 3; a++) {
      f[a] -= std::floor(f[a]);
    }
    double wrapped[3] = {0, 0, 0};
    for (int c = 0; c < 3; c++) {
      for (int a = 0; a < 3; a++) {
        wrapped[c] += f[a] * lat[a][c];
      }
    }
    px[static_cast<size_t>(i)] = wrapped[0];
    py[static_cast<size_t>(i)] = wrapped[1];
    pz[static_cast<size_t>(i)] = wrapped[2];
    int ix = std::min(nx - 1, static_cast<int>(f[0] * nx));
    int iy = std::min(ny - 1, static_cast<int>(f[1] * ny));
    int iz = std::min(nz - 1, static_cast<int>(f[2] * nz));
    if (ix < 0) {
      ix = 0;
    }
    if (iy < 0) {
      iy = 0;
    }
    if (iz < 0) {
      iz = 0;
    }
    const int cell = (iz * ny + iy) * nx + ix;
    bin[static_cast<size_t>(i)] = cell;
  }
  // Atoms packed by cell, which is the layout of the orthorhombic neighbour list.
  std::vector<int> cellStart(static_cast<size_t>(ncell) + 1, 0);
  for (int i = 0; i < n; i++) {
    ++cellStart[static_cast<size_t>(bin[static_cast<size_t>(i)]) + 1];
  }
  for (int c = 0; c < ncell; c++) {
    cellStart[static_cast<size_t>(c) + 1] += cellStart[static_cast<size_t>(c)];
  }
  std::vector<int> fill(cellStart.begin(), cellStart.end() - 1);
  RawBuf<Packed> packed(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    const int slot = fill[static_cast<size_t>(bin[static_cast<size_t>(i)])]++;
    packed[static_cast<size_t>(slot)] = Packed{px[static_cast<size_t>(i)], py[static_cast<size_t>(i)],
                                               pz[static_cast<size_t>(i)], i};
  }
  const int span = 2 * reach + 1;
  const int nImage = span * span * span;
  RawBuf<Image> images(static_cast<size_t>(ncell) * static_cast<size_t>(nImage));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (ncell >= 256)
#endif
  for (int cell = 0; cell < ncell; cell++) {
    const int ix0 = cell % nx;
    const int iy0 = (cell / nx) % ny;
    const int iz0 = cell / (nx * ny);
    int slot = 0;
    for (int dz = -reach; dz <= reach; dz++) {
      int iz = iz0 + dz;
      int sz = 0;
      while (iz < 0) {
        iz += nz;
        sz--;
      }
      while (iz >= nz) {
        iz -= nz;
        sz++;
      }
      for (int dy = -reach; dy <= reach; dy++) {
        int iy = iy0 + dy;
        int sy = 0;
        while (iy < 0) {
          iy += ny;
          sy--;
        }
        while (iy >= ny) {
          iy -= ny;
          sy++;
        }
        for (int dx = -reach; dx <= reach; dx++) {
          int ix = ix0 + dx;
          int sx = 0;
          while (ix < 0) {
            ix += nx;
            sx--;
          }
          while (ix >= nx) {
            ix -= nx;
            sx++;
          }
          images[static_cast<size_t>(cell) * static_cast<size_t>(nImage) + static_cast<size_t>(slot)] =
              Image{(iz * ny + iy) * nx + ix, sx, sy, sz};
          slot++;
        }
      }
    }
  }

  RawBuf<Bond> bonds(static_cast<size_t>(n) * static_cast<size_t>(kMaxBond));
  std::vector<int> nbon(static_cast<size_t>(n), 0);
#ifdef SEAMS_HAS_OPENMP
  const int threads = omp_get_max_threads();
#pragma omp parallel for schedule(static) if (n >= 64 && threads > 1)
#endif
  for (int cell = 0; cell < ncell; cell++) {
    const Image *im =
        &images[static_cast<size_t>(cell) * static_cast<size_t>(nImage)];
    const int here0 = cellStart[static_cast<size_t>(cell)];
    const int here1 = cellStart[static_cast<size_t>(cell) + 1];
    for (int mine = here0; mine < here1; mine++) {
    const int i = packed[static_cast<size_t>(mine)].id;
    const double ci0 = packed[static_cast<size_t>(mine)].x;
    const double ci1 = packed[static_cast<size_t>(mine)].y;
    const double ci2 = packed[static_cast<size_t>(mine)].z;
    Bond row[kMaxBond];
    double rd[kMaxBond];
    int nb = 0;
    for (int t = 0; t < nImage; t++) {
      const int sx = im[t].sx;
      const int sy = im[t].sy;
      const int sz = im[t].sz;
      const int begin = cellStart[static_cast<size_t>(im[t].cell)];
      const int end = cellStart[static_cast<size_t>(im[t].cell) + 1];
      for (int m = begin; m < end; m++) {
        const Packed &nbr = packed[static_cast<size_t>(m)];
        if (nbr.id == i && sx == 0 && sy == 0 && sz == 0) {
          continue;
        }
        const double dx = nbr.x + sx * lat[0][0] + sy * lat[1][0] + sz * lat[2][0] - ci0;
        const double dy = nbr.y + sx * lat[0][1] + sy * lat[1][1] + sz * lat[2][1] - ci1;
        const double dz = nbr.z + sx * lat[0][2] + sy * lat[1][2] + sz * lat[2][2] - ci2;
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 > kBond * kBond || d2 < 2.15 * 2.15) {
          continue;
        }
        if (nb < kMaxBond) {
          row[nb] = Bond{nbr.id, dx, dy, dz};
          rd[nb] = d2;
          nb++;
        } else {
          int worst = 0;
          for (int c = 1; c < kMaxBond; c++) {
            if (rd[c] > rd[worst]) {
              worst = c;
            }
          }
          if (d2 < rd[worst]) {
            row[worst] = Bond{nbr.id, dx, dy, dz};
            rd[worst] = d2;
          }
        }
      }
    }
    for (int a = 1; a < nb; a++) {
      const Bond b = row[a];
      const double d2 = rd[a];
      int at = a;
      while (at > 0 && d2 < rd[at - 1]) {
        row[at] = row[at - 1];
        rd[at] = rd[at - 1];
        at--;
      }
      row[at] = b;
      rd[at] = d2;
    }
    nbon[static_cast<size_t>(i)] = nb;
    Bond *dest = &bonds[static_cast<size_t>(i) * static_cast<size_t>(kMaxBond)];
    for (int a = 0; a < nb; a++) {
      dest[a] = row[a];
    }
    }
  }

#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 64 && threads > 1)
#endif
  for (int i = 0; i < n; i++) {
    const int nb = nbon[static_cast<size_t>(i)];
    const Bond *mine = &bonds[static_cast<size_t>(i) * static_cast<size_t>(kMaxBond)];
    int firstId[4];
    double nbxyz[16][3];
    const double ci0 = px[static_cast<size_t>(i)];
    const double ci1 = py[static_cast<size_t>(i)];
    const double ci2 = pz[static_cast<size_t>(i)];
    if (nb < 4) {
      kind[i] = -1;
      continue;
    }
    double uvec[4][3];
    double acc = 0.0;
    bool tetra = true;
    for (int k = 0; k < 4; k++) {
      const double nrm = std::sqrt(mine[k].dx * mine[k].dx + mine[k].dy * mine[k].dy +
                                   mine[k].dz * mine[k].dz);
      if (!(nrm > 1e-8)) {
        tetra = false;
        break;
      }
      uvec[k][0] = mine[k].dx / nrm;
      uvec[k][1] = mine[k].dy / nrm;
      uvec[k][2] = mine[k].dz / nrm;
    }
    if (tetra) {
      for (int a = 0; a < 4; a++) {
        for (int b = a + 1; b < 4; b++) {
          const double c = uvec[a][0] * uvec[b][0] + uvec[a][1] * uvec[b][1] + uvec[a][2] * uvec[b][2];
          const double d = c + (1.0 / 3.0);
          acc += d * d;
        }
      }
    }
    // q = 1 - (3/8) sum (cos + 1/3)^2. Below 0.88 the four neighbours are not
    // a template tetrahedron, so the second shell is not walked.
    if (!tetra || (1.0 - 0.375 * acc) < 0.88) {
      kind[i] = -1;
      continue;
    }
    // Cheap diamond test, before the full shell fit. The twelve second
    // neighbours are the three outward bonds of each of the four neighbours.
    // A perfect FCC or HCP shell is labelled here and does not enter Kabsch.
    {
      double first[4][3];
      for (int k = 0; k < 4; k++) {
        first[k][0] = mine[k].dx;
        first[k][1] = mine[k].dy;
        first[k][2] = mine[k].dz;
      }
      double fx[12];
      double fy[12];
      double fz[12];
      int nf = 0;
      bool closed = true;
      for (int k = 0; k < 4 && closed; k++) {
        const int f = mine[k].j;
        const int na = nbon[static_cast<size_t>(f)];
        const Bond *arm = &bonds[static_cast<size_t>(f) * static_cast<size_t>(kMaxBond)];
        int taken = 0;
        for (int t = 0; t < na && taken < 3; t++) {
          const double x = mine[k].dx + arm[t].dx;
          const double y = mine[k].dy + arm[t].dy;
          const double z = mine[k].dz + arm[t].dz;
          const double origin = x * x + y * y + z * z;
          if (origin < 0.25) {
            continue;
          }
          bool isFirst = false;
          for (int q = 0; q < 4; q++) {
            const double dx = x - mine[q].dx;
            const double dy = y - mine[q].dy;
            const double dz = z - mine[q].dz;
            if (dx * dx + dy * dy + dz * dz < 0.25) {
              isFirst = true;
            }
          }
          if (isFirst) {
            continue;
          }
          if (nf >= 12) {
            closed = false;
            break;
          }
          fx[nf] = x;
          fy[nf] = y;
          fz[nf] = z;
          nf++;
          taken++;
        }
        if (taken < 3) {
          closed = false;
        }
      }
      if (closed && nf == 12) {
        const int fast = diamondSecondShell(first, fx, fy, fz);
        if (fast >= 0) {
          kind[i] = fast;
          continue;
        }
      }
    }
    for (int k = 0; k < 4; k++) {
      firstId[k] = mine[k].j;
      nbxyz[k][0] = ci0 + mine[k].dx;
      nbxyz[k][1] = ci1 + mine[k].dy;
      nbxyz[k][2] = ci2 + mine[k].dz;
    }
    double sx[48];
    double sy[48];
    double sz[48];
    double sb[48];
    int ns = 0;
    for (int k = 0; k < 4 && nb >= 4; k++) {
      const int f = firstId[k];
      const int na = nbon[static_cast<size_t>(f)];
      const Bond *arm = &bonds[static_cast<size_t>(f) * static_cast<size_t>(kMaxBond)];
      for (int t = 0; t < na; t++) {
        const double x = mine[k].dx + arm[t].dx;
        const double y = mine[k].dy + arm[t].dy;
        const double z = mine[k].dz + arm[t].dz;
        // A periodic image is a different shell site. Reject the centre and
        // the four first-shell positions, not every copy of their atom id.
        const double origin = x * x + y * y + z * z;
        if (origin < 0.25 || origin > kCentre * kCentre) {
          continue;
        }
        bool isFirst = false;
        for (int q = 0; q < 4; q++) {
          const double dx = x - mine[q].dx;
          const double dy = y - mine[q].dy;
          const double dz = z - mine[q].dz;
          if (dx * dx + dy * dy + dz * dz < 0.25) {
            isFirst = true;
          }
        }
        if (isFirst) {
          continue;
        }
        const double bond = arm[t].dx * arm[t].dx + arm[t].dy * arm[t].dy + arm[t].dz * arm[t].dz;
        int found = -1;
        for (int s = 0; s < ns; s++) {
          const double dx = x - sx[s];
          const double dy = y - sy[s];
          const double dz = z - sz[s];
          if (dx * dx + dy * dy + dz * dz < 0.25) {
            found = s;
          }
        }
        if (found >= 0) {
          if (bond < sb[found]) {
            sx[found] = x;
            sy[found] = y;
            sz[found] = z;
            sb[found] = bond;
          }
          continue;
        }
        if (ns < 48) {
          sx[ns] = x;
          sy[ns] = y;
          sz[ns] = z;
          sb[ns] = bond;
          ns++;
        }
      }
    }
    // An incomplete 3.5 A shell is not ice. Widening the bond pulled every
    // liquid atom into a second search.
    if (ns < 12 || nb < 4) {
      kind[i] = -1;
      continue;
    }

    for (int a = 0; a < 12; a++) {
      int best = a;
      for (int b = a + 1; b < ns; b++) {
        if (sb[b] < sb[best]) {
          best = b;
        }
      }
      std::swap(sb[a], sb[best]);
      std::swap(sx[a], sx[best]);
      std::swap(sy[a], sy[best]);
      std::swap(sz[a], sz[best]);
      nbxyz[4 + a][0] = ci0 + sx[a];
      nbxyz[4 + a][1] = ci1 + sy[a];
      nbxyz[4 + a][2] = ci2 + sz[a];
    }
    const double ci[3] = {ci0, ci1, ci2};
    kind[i] = shellFit(ci, nbxyz, 16, canon, gate);
  }
  // One crystal neighbour among the four bonds. Isolated liquid labels go.
  // The test reads the labels from the fit, not the filtered ones.
  std::vector<int> raw(kind, kind + n);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 64 && threads > 1)
#endif
  for (int i = 0; i < n; i++) {
    if (raw[static_cast<size_t>(i)] < 0) {
      continue;
    }
    const int nb = nbon[static_cast<size_t>(i)];
    const Bond *mine = &bonds[static_cast<size_t>(i) * static_cast<size_t>(kMaxBond)];
    int nCrystal = 0;
    const int lim = nb < 4 ? nb : 4;
    for (int k = 0; k < lim; k++) {
      if (raw[static_cast<size_t>(mine[k].j)] >= 0) {
        nCrystal++;
      }
    }
    if (nCrystal < 1) {
      kind[i] = -1;
    }
  }
  recoverHex(n, bonds.get(), nbon.data(), packed.get(), cellStart.data(), images.get(), nImage, bin.get(),
             lat, px.get(), py.get(), pz.get(), kind, canon, gate);
}

int bondsBeyond(int atom, const Packed *packed, const int *cellStart, const Image *images, int nImage,
               const int *bin, const double lat[3][3], const double *px, const double *py, const double *pz,
               Bond *out) {
  Bond row[kMaxWide];
  double rd[kMaxWide];
  int nw = 0;
  const int cell = bin[atom];
  const Image *im = images + static_cast<size_t>(cell) * static_cast<size_t>(nImage);
  const double ci0 = px[atom];
  const double ci1 = py[atom];
  const double ci2 = pz[atom];
  for (int t = 0; t < nImage; t++) {
    const int begin = cellStart[im[t].cell];
    const int end = cellStart[im[t].cell + 1];
    for (int m = begin; m < end; m++) {
      const Packed &nbr = packed[m];
      if (nbr.id == atom && im[t].sx == 0 && im[t].sy == 0 && im[t].sz == 0) {
        continue;
      }
      const double dx = nbr.x + im[t].sx * lat[0][0] + im[t].sy * lat[1][0] + im[t].sz * lat[2][0] - ci0;
      const double dy = nbr.y + im[t].sx * lat[0][1] + im[t].sy * lat[1][1] + im[t].sz * lat[2][1] - ci1;
      const double dz = nbr.z + im[t].sx * lat[0][2] + im[t].sy * lat[1][2] + im[t].sz * lat[2][2] - ci2;
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 <= kBond * kBond || d2 > kWide * kWide) {
        continue;
      }
      if (nw < kMaxWide) {
        row[nw] = Bond{nbr.id, dx, dy, dz};
        rd[nw] = d2;
        nw++;
      } else {
        int worst = 0;
        for (int c = 1; c < kMaxWide; c++) {
          if (rd[c] > rd[worst]) {
            worst = c;
          }
        }
        if (d2 < rd[worst]) {
          row[worst] = Bond{nbr.id, dx, dy, dz};
          rd[worst] = d2;
        }
      }
    }
  }
  for (int a = 1; a < nw; a++) {
    const Bond b = row[a];
    const double d2 = rd[a];
    int at = a;
    while (at > 0 && d2 < rd[at - 1]) {
      row[at] = row[at - 1];
      rd[at] = rd[at - 1];
      at--;
    }
    row[at] = b;
    rd[at] = d2;
  }
  for (int a = 0; a < nw; a++) {
    out[a] = row[a];
  }
  return nw;
}

void recoverHex(int n, const Bond *bonds, const int *nbon, const Packed *packed, const int *cellStart,
                const Image *images, int nImage, const int *bin, const double lat[3][3], const double *px,
                const double *py, const double *pz, int *kind, const CanonShell &canon, double gate) {
  int unlabeled = 0;
  for (int i = 0; i < n; i++) {
    if (kind[i] < 0) {
      unlabeled++;
    }
  }
  if (unlabeled == 0) {
    return;
  }
  std::vector<char> crystal(static_cast<size_t>(n), 0);
  std::vector<char> tried(static_cast<size_t>(n), 0);
  for (int pass = 0; pass < 2; pass++) {
  for (int i = 0; i < n; i++) {
    crystal[static_cast<size_t>(i)] = kind[i] >= 0 ? 1 : 0;
  }
#ifdef SEAMS_HAS_OPENMP
  const int threads = omp_get_max_threads();
#pragma omp parallel for schedule(static) if (n >= 64 && threads > 1)
#endif
  for (int i = 0; i < n; i++) {
    if (crystal[static_cast<size_t>(i)] || tried[static_cast<size_t>(i)]) {
      continue;
    }
    const int nb = nbon[i];
    const Bond *mine = bonds + static_cast<size_t>(i) * kMaxBond;
    Bond mineWbuf[kMaxWide];
    const int nw = bondsBeyond(i, packed, cellStart, images, nImage, bin, lat, px, py, pz, mineWbuf);
    const Bond *mineW = mineWbuf;
    if (nb + nw < 4) {
      continue;
    }
    int nCrystal = 0;
    const int lim = nb < 4 ? nb : 4;
    for (int k = 0; k < lim; k++) {
      if (crystal[static_cast<size_t>(mine[k].j)]) {
        nCrystal++;
      }
    }
    if (nCrystal < 1 && nw > 0 && nb >= 3) {
      for (int k = 0; k < nb && k < 3; k++) {
        if (crystal[static_cast<size_t>(mine[k].j)]) {
          nCrystal++;
        }
      }
    }
    if (nCrystal < 3) {
      continue;
    }
    tried[static_cast<size_t>(i)] = 1;
    Bond first[4];
    int nf = 0;
    for (int k = 0; k < nb && nf < 4; k++) {
      first[nf++] = mine[k];
    }
    if (nf < 4 && nw > 0) {
      const double d4 = mineW[0].dx * mineW[0].dx + mineW[0].dy * mineW[0].dy + mineW[0].dz * mineW[0].dz;
      double rNext = kWide;
      if (nw >= 2) {
        const double d5 = mineW[1].dx * mineW[1].dx + mineW[1].dy * mineW[1].dy + mineW[1].dz * mineW[1].dz;
        rNext = std::sqrt(d5);
      }
      const double r4 = std::sqrt(d4);
      if (rNext - r4 >= 0.12) {
        first[nf++] = mineW[0];
      }
    }
    if (nf < 4) {
      continue;
    }
    const double ci0 = px[i];
    const double ci1 = py[i];
    const double ci2 = pz[i];
    double nbxyz[16][3];
    for (int k = 0; k < 4; k++) {
      nbxyz[k][0] = ci0 + first[k].dx;
      nbxyz[k][1] = ci1 + first[k].dy;
      nbxyz[k][2] = ci2 + first[k].dz;
    }
    double sx[48], sy[48], sz[48], sb[48];
    int ns = 0;
    for (int k = 0; k < 4; k++) {
      const int f = first[k].j;
      const Bond *arm = bonds + static_cast<size_t>(f) * kMaxBond;
      const int na = nbon[f];
      Bond armWbuf[kMaxWide];
      const int naw = bondsBeyond(f, packed, cellStart, images, nImage, bin, lat, px, py, pz, armWbuf);
      const Bond *armW = armWbuf;
      int seen = 0;
      auto push = [&](const Bond &b) {
        const double x = first[k].dx + b.dx;
        const double y = first[k].dy + b.dy;
        const double z = first[k].dz + b.dz;
        const double origin = x * x + y * y + z * z;
        if (origin < 0.25 || origin > kCentreHex * kCentreHex) {
          return;
        }
        for (int q = 0; q < 4; q++) {
          const double dx = x - first[q].dx;
          const double dy = y - first[q].dy;
          const double dz = z - first[q].dz;
          if (dx * dx + dy * dy + dz * dz < 0.25) {
            return;
          }
        }
        const double bond = b.dx * b.dx + b.dy * b.dy + b.dz * b.dz;
        for (int s = 0; s < ns; s++) {
          const double dx = x - sx[s];
          const double dy = y - sy[s];
          const double dz = z - sz[s];
          if (dx * dx + dy * dy + dz * dz < 0.25) {
            if (bond < sb[s]) {
              sx[s] = x;
              sy[s] = y;
              sz[s] = z;
              sb[s] = bond;
            }
            return;
          }
        }
        if (ns < 48) {
          sx[ns] = x;
          sy[ns] = y;
          sz[ns] = z;
          sb[ns] = bond;
          ns++;
        }
      };
      for (int t = 0; t < na && seen < 4; t++) {
        seen++;
        push(arm[t]);
      }
      if (seen < 4 && naw > 0) {
        const double d4 = armW[0].dx * armW[0].dx + armW[0].dy * armW[0].dy + armW[0].dz * armW[0].dz;
        double rNext = kWide;
        if (naw >= 2) {
          const double d5 = armW[1].dx * armW[1].dx + armW[1].dy * armW[1].dy + armW[1].dz * armW[1].dz;
          rNext = std::sqrt(d5);
        }
        if (rNext - std::sqrt(d4) >= 0.12) {
          push(armW[0]);
        }
      }
    }
    if (ns < 12) {
      continue;
    }
    for (int a = 0; a < 12; a++) {
      int best = a;
      for (int b = a + 1; b < ns; b++) {
        if (sb[b] < sb[best]) {
          best = b;
        }
      }
      std::swap(sb[a], sb[best]);
      std::swap(sx[a], sx[best]);
      std::swap(sy[a], sy[best]);
      std::swap(sz[a], sz[best]);
      nbxyz[4 + a][0] = ci0 + sx[a];
      nbxyz[4 + a][1] = ci1 + sy[a];
      nbxyz[4 + a][2] = ci2 + sz[a];
    }
    const double ci[3] = {ci0, ci1, ci2};
    kind[i] = shellFit(ci, nbxyz, 16, canon, gate, true);
  }
  }
}

} // namespace symm
