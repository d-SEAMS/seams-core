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
#include <cstddef>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace symm {
namespace {

constexpr double kShell = 5.6;
constexpr double kGrid = 2.8;
constexpr double kCleanRms = 0.02;
constexpr double kCleanShift = 0.05;

struct Image {
  int cell;
  int sx;
  int sy;
  int sz;
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
void jacobi4(double A[4][4], double Q[4][4]) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      Q[i][j] = i == j ? 1.0 : 0.0;
    }
  }
  for (int sweep = 0; sweep < 4; sweep++) {
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
  jacobi4(K, Q);
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

void frameOf(const double a[3], const double b[3], double e[3][3]) {
  double an = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  e[0][0] = a[0] / an;
  e[0][1] = a[1] / an;
  e[0][2] = a[2] / an;
  const double dot = b[0] * e[0][0] + b[1] * e[0][1] + b[2] * e[0][2];
  e[1][0] = b[0] - dot * e[0][0];
  e[1][1] = b[1] - dot * e[0][1];
  e[1][2] = b[2] - dot * e[0][2];
  double bn = std::sqrt(e[1][0] * e[1][0] + e[1][1] * e[1][1] + e[1][2] * e[1][2]);
  if (bn < 1e-12) {
    bn = 1.0;
  }
  e[1][0] /= bn;
  e[1][1] /= bn;
  e[1][2] /= bn;
  e[2][0] = e[0][1] * e[1][2] - e[0][2] * e[1][1];
  e[2][1] = e[0][2] * e[1][0] - e[0][0] * e[1][2];
  e[2][2] = e[0][0] * e[1][1] - e[0][1] * e[1][0];
}

// Identity pairing of a gapped tetrahedron. The rotation is the bond frame,
// not the 4x4 eigenproblem. Returns -2 when the shell is not that clean, so
// the caller keeps the full fit and the label does not change.
int cleanLabel(const double center[3], const double nb[][3], const CanonShell &canon) {
  double r[5];
  for (int k = 0; k < 5; k++) {
    r[k] = std::sqrt(dist2(center, nb[k]));
  }
  if (!(r[3] - r[0] < 0.2 && r[4] - r[3] > 0.8)) {
    return -2;
  }
  double env[4][3];
  for (int k = 0; k < 4; k++) {
    for (int c = 0; c < 3; c++) {
      env[k][c] = nb[k][c] - center[c];
    }
  }
  double envS[12][3];
  for (int s = 0; s < 12; s++) {
    for (int c = 0; c < 3; c++) {
      envS[s][c] = nb[4 + s][c] - center[c];
    }
  }
  const int kinds[3] = {1, 0, 2};
  for (int ko = 0; ko < 3; ko++) {
    const int kind = kinds[ko];
    if (tetraChirality(canon.first[kind]) * tetraChirality(env) < 0.0) {
      continue;
    }
    double eF[3][3];
    double tF[3][3];
    frameOf(env[0], env[1], eF);
    frameOf(canon.first[kind][0], canon.first[kind][1], tF);
    // R maps the environment frame onto the template frame.
    double R[3][3];
    for (int row = 0; row < 3; row++) {
      for (int col = 0; col < 3; col++) {
        R[row][col] = tF[0][row] * eF[0][col] + tF[1][row] * eF[1][col] + tF[2][row] * eF[2][col];
      }
    }
    double gT = 0.0;
    double gE = 0.0;
    for (int k = 0; k < 4; k++) {
      for (int c = 0; c < 3; c++) {
        gT += canon.first[kind][k][c] * canon.first[kind][k][c];
        gE += env[k][c] * env[k][c];
      }
    }
    const double shellScale = gT > 1e-18 ? std::sqrt(gE / gT) : 1.0;
    bool taken[12] = {};
    double shiftSum = 0.0;
    bool paired = true;
    for (int a = 0; a < 4 && paired; a++) {
      int cand[3];
      double candD[3];
      int nC = 0;
      for (int e = 0; e < 12; e++) {
        if (taken[e]) {
          continue;
        }
        const double dx = envS[e][0] - env[a][0];
        const double dy = envS[e][1] - env[a][1];
        const double dz = envS[e][2] - env[a][2];
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
      if (nC < 3) {
        paired = false;
        break;
      }
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
      double bestShiftErr = 1e300;
      int bestShift = 0;
      int bestDir = 0;
      int seq[2][3];
      // Order candidates by azimuth so the shift matches the eigenproblem path.
      double e1[3];
      double e2[3];
      bondFrame(env[a], env[(a + 1) % 4], e1, e2);
      double ang[3];
      for (int i = 0; i < 3; i++) {
        const double outer[3] = {envS[cand[i]][0] - env[a][0], envS[cand[i]][1] - env[a][1],
                                 envS[cand[i]][2] - env[a][2]};
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
      for (int i = 0; i < 3; i++) {
        seq[0][i] = cand[i];
        seq[1][i] = cand[2 - i];
      }
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
      (void)bestDir;
      (void)bestShift;
      for (int i = 0; i < 3; i++) {
        taken[seq[bestDir][i]] = true;
      }
      shiftSum += bestShiftErr;
    }
    if (paired && shiftSum <= kCleanShift) {
      return kind;
    }
  }
  return -2;
}

int shellFit(const double center[3], const double nb[][3], int nNb, const CanonShell &canon,
             double gate) {
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
      if (tried > 12 || (tried > 3 && closest > 40.0 && closest < 1.0e200)) {
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

void classifyPrepared(const double *xyz, int n, const double lat[3][3], const CanonShell &canon,
                      double gate, int *kind) {
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
  if (ortho && la >= 3.0 * kShell && lb >= 3.0 * kShell && lc >= 3.0 * kShell) {
    // Cutoff-sized cells, the orthorhombic neighbour list. The 27-cell stencil
    // is complete because each cell is at most the shell cutoff wide.
    nx = std::max(3, static_cast<int>(std::ceil(la / kShell)));
    ny = std::max(3, static_cast<int>(std::ceil(lb / kShell)));
    nz = std::max(3, static_cast<int>(std::ceil(lc / kShell)));
    fast = true;
  } else {
    nx = std::max(1, static_cast<int>(la / kGrid));
    ny = std::max(1, static_cast<int>(lb / kGrid));
    nz = std::max(1, static_cast<int>(lc / kGrid));
  }
  const int ncell = nx * ny * nz;
  const int reach = fast ? 1 : (skew > 0.5 * kShell ? 4 : 2);
  std::vector<double> px(static_cast<size_t>(n));
  std::vector<double> py(static_cast<size_t>(n));
  std::vector<double> pz(static_cast<size_t>(n));
  std::vector<int> bin(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    const double c0 = xyz[3 * i];
    const double c1 = xyz[3 * i + 1];
    const double c2 = xyz[3 * i + 2];
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
  struct Packed {
    double x, y, z;
    int id;
  };
  std::vector<int> cellStart(static_cast<size_t>(ncell) + 1, 0);
  for (int i = 0; i < n; i++) {
    ++cellStart[static_cast<size_t>(bin[static_cast<size_t>(i)]) + 1];
  }
  for (int c = 0; c < ncell; c++) {
    cellStart[static_cast<size_t>(c) + 1] += cellStart[static_cast<size_t>(c)];
  }
  std::vector<int> fill(cellStart.begin(), cellStart.end() - 1);
  std::vector<Packed> packed(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    const int slot = fill[static_cast<size_t>(bin[static_cast<size_t>(i)])]++;
    packed[static_cast<size_t>(slot)] = Packed{px[static_cast<size_t>(i)], py[static_cast<size_t>(i)],
                                               pz[static_cast<size_t>(i)], i};
  }
  const int span = 2 * reach + 1;
  const int nImage = span * span * span;
  std::vector<Image> images(static_cast<size_t>(ncell) * static_cast<size_t>(nImage));
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

#ifdef SEAMS_HAS_OPENMP
  const int threads = omp_get_max_threads();
#pragma omp parallel for schedule(static) if (n >= 64 && threads > 1)
#endif
  for (int i = 0; i < n; i++) {
    const double ci[3] = {px[static_cast<size_t>(i)], py[static_cast<size_t>(i)],
                          pz[static_cast<size_t>(i)]};
    double kd[48];
    double kx[48];
    double ky[48];
    double kz[48];
    int nkeep = 0;
    const Image *im =
        &images[static_cast<size_t>(bin[static_cast<size_t>(i)]) * static_cast<size_t>(nImage)];
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
        const double x = nbr.x + sx * lat[0][0] + sy * lat[1][0] + sz * lat[2][0];
        const double y = nbr.y + sx * lat[0][1] + sy * lat[1][1] + sz * lat[2][1];
        const double z = nbr.z + sx * lat[0][2] + sy * lat[1][2] + sz * lat[2][2];
        const double ddx = x - ci[0];
        const double ddy = y - ci[1];
        const double ddz = z - ci[2];
        const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
        if (d2 > kShell * kShell || d2 < 2.15 * 2.15) {
          continue;
        }
        if (nkeep < 48) {
          kd[nkeep] = d2;
          kx[nkeep] = x;
          ky[nkeep] = y;
          kz[nkeep] = z;
          nkeep++;
        }
      }
    }
    for (int a = 1; a < nkeep; a++) {
      const double d2 = kd[a];
      const double x = kx[a];
      const double y = ky[a];
      const double z = kz[a];
      int at = a;
      while (at > 0 && d2 < kd[at - 1]) {
        kd[at] = kd[at - 1];
        kx[at] = kx[at - 1];
        ky[at] = ky[at - 1];
        kz[at] = kz[at - 1];
        at--;
      }
      kd[at] = d2;
      kx[at] = x;
      ky[at] = y;
      kz[at] = z;
    }
    if (nkeep < 4) {
      kind[i] = -1;
      continue;
    }
    const int pool = std::min(nkeep, 8);
    int pick[4] = {0, 1, 2, 3};
    double pickScore = 0.0;
    const bool outlier = std::sqrt(kd[std::min(3, nkeep - 1)]) - std::sqrt(kd[0]) > 0.45;
    if (outlier) {
      pickScore = 1e300;
    }
    for (int a = 0; a < pool && outlier; a++) {
      for (int b = a + 1; b < pool; b++) {
        for (int c = b + 1; c < pool; c++) {
          for (int d = c + 1; d < pool; d++) {
            const int s[4] = {a, b, c, d};
            double rad[4];
            double rsum = 0.0;
            for (int t = 0; t < 4; t++) {
              rad[t] = std::sqrt(kd[s[t]]);
              rsum += rad[t];
            }
            const double rmean = rsum / 4.0;
            if (rmean < 2.3 || rmean > 3.5) {
              continue;
            }
            double rvar = 0.0;
            for (int t = 0; t < 4; t++) {
              const double dr = rad[t] - rmean;
              rvar += dr * dr;
            }
            double edges[6];
            int ne = 0;
            double esum = 0.0;
            for (int t = 0; t < 4; t++) {
              for (int u = t + 1; u < 4; u++) {
                const double dx = kx[s[t]] - kx[s[u]];
                const double dy = ky[s[t]] - ky[s[u]];
                const double dz = kz[s[t]] - kz[s[u]];
                edges[ne] = std::sqrt(dx * dx + dy * dy + dz * dz);
                esum += edges[ne];
                ne++;
              }
            }
            const double emean = esum / 6.0;
            double evar = 0.0;
            for (int t = 0; t < 6; t++) {
              const double de = edges[t] - emean;
              evar += de * de;
            }
            const double score = rvar + 0.25 * evar;
            if (score < pickScore) {
              pickScore = score;
              for (int t = 0; t < 4; t++) {
                pick[t] = s[t];
              }
            }
          }
        }
      }
    }
    if (!(pickScore < 1.0e200)) {
      kind[i] = -1;
      continue;
    }
    double nb[16][3];
    for (int k = 0; k < 4; k++) {
      nb[k][0] = kx[pick[k]];
      nb[k][1] = ky[pick[k]];
      nb[k][2] = kz[pick[k]];
    }
    int bondAt[48];
    double bondTo[48];
    int nBond = 0;
    for (int k = 0; k < nkeep; k++) {
      bool taken = false;
      for (int f = 0; f < 4; f++) {
        if (k == pick[f]) {
          taken = true;
        }
      }
      if (taken) {
        continue;
      }
      double best = 1e300;
      for (int f = 0; f < 4; f++) {
        const double dx = kx[k] - nb[f][0];
        const double dy = ky[k] - nb[f][1];
        const double dz = kz[k] - nb[f][2];
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < best) {
          best = d2;
        }
      }
      if (best < 4.2 * 4.2 && nBond < 48) {
        bondAt[nBond] = k;
        bondTo[nBond] = best;
        nBond++;
      }
    }
    for (int a = 0; a < nBond; a++) {
      for (int b = a + 1; b < nBond; b++) {
        if (bondTo[b] < bondTo[a]) {
          std::swap(bondTo[a], bondTo[b]);
          std::swap(bondAt[a], bondAt[b]);
        }
      }
    }
    int nSecond = 0;
    for (int a = 0; a < nBond && nSecond < 12; a++) {
      const int k = bondAt[a];
      nb[4 + nSecond][0] = kx[k];
      nb[4 + nSecond][1] = ky[k];
      nb[4 + nSecond][2] = kz[k];
      nSecond++;
    }
    if (nSecond < 12) {
      kind[i] = -1;
      continue;
    }
    const int quick = cleanLabel(ci, nb, canon);
    kind[i] = quick == -2 ? shellFit(ci, nb, 16, canon, gate) : quick;
  }
}

} // namespace symm
