//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Fixed-size shells. The ice label is the scale-free residual of the
// seventeen-atom shell after the optimal rotation and scale. The second
// shell is the atoms bonded to the tetrahedron. The rotation is the
// quaternion form of Kabsch used by polyhedral template matching
// (Kabsch 1976; Theobald 2005; Larsen, Schmidt, Schiotz 2016). The
// continuous shape measure is the scaled residual of Alvarez et al. 2005
// after that alignment (Pinsky and Avnir 1998).

#include <local_symmetry.hpp>
#include <shell_hot.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

constexpr int kMaxN = 32;
constexpr double kBondLo = 2.35;
constexpr double kBondHi = 3.25;
constexpr double kShell = 5.6;
constexpr double kGrid = 2.8;
constexpr double kEclipsed = 0.90;

// Rotate mobile (env) onto fixed (template). Proper rotation.
double rmsdOf(const double fixed[][3], const double mobile[][3], int n, const double R[3][3]);

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

bool kabschQuat(const double fixed[][3], const double mobile[][3], int n, double R[3][3],
                double &rmsd, double *scaleFree = nullptr) {
  double cf[3] = {0, 0, 0};
  double cm[3] = {0, 0, 0};
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < 3; k++) {
      cf[k] += fixed[i][k];
      cm[k] += mobile[i][k];
    }
  }
  const double inv = 1.0 / static_cast<double>(std::max(n, 1));
  for (int k = 0; k < 3; k++) {
    cf[k] *= inv;
    cm[k] *= inv;
  }
  double H[3][3] = {};
  for (int i = 0; i < n; i++) {
    const double p[3] = {fixed[i][0] - cf[0], fixed[i][1] - cf[1], fixed[i][2] - cf[2]};
    const double q[3] = {mobile[i][0] - cm[0], mobile[i][1] - cm[1], mobile[i][2] - cm[2]};
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        H[a][b] += p[a] * q[b];
      }
    }
  }
  // K such that its largest eigenvector is the quaternion taking q onto p.
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
  // The eigenvector is the quaternion of the inverse map. Conjugating it
  // takes the mobile point onto the fixed point.
  double q[4] = {Q[0][best], -Q[1][best], -Q[2][best], -Q[3][best]};
  const double qn = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  if (qn > 1e-15) {
    for (int i = 0; i < 4; i++) {
      q[i] /= qn;
    }
  }
  quatToR(q, R);
  // Theobald: λ_max of K is the maximum of sum p·(R q). With both clouds
  // centered, N * RMSD^2 = |p|^2 + |q|^2 - 2 λ_max.
  double gp = 0.0;
  double gq = 0.0;
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < 3; k++) {
      const double dp = fixed[i][k] - cf[k];
      const double dq = mobile[i][k] - cm[k];
      gp += dp * dp;
      gq += dq * dq;
    }
  }
  const double lam = K[best][best];
  const double acc = gp + gq - 2.0 * lam;
  rmsd = std::sqrt(std::max(0.0, acc) / static_cast<double>(n));
  if (scaleFree != nullptr) {
    // Optimal scale s = λ / |q|^2 removes the lattice constant. The
    // unexplained fraction of the template is then 1 - λ^2 / (|p|^2 |q|^2)
    // (Horn 1987; Larsen, Schmidt, Schiotz 2016).
    double frac = 1.0;
    if (gp > 1e-18 && gq > 1e-18) {
      frac = std::sqrt(std::max(0.0, 1.0 - (lam * lam) / (gp * gq)));
    }
    *scaleFree = frac;
  }
  return true;
}

void applyR(const double R[3][3], const double v[3], double out[3]) {
  for (int i = 0; i < 3; i++) {
    out[i] = R[i][0] * v[0] + R[i][1] * v[1] + R[i][2] * v[2];
  }
}

double rmsdOf(const double fixed[][3], const double mobile[][3], int n, const double R[3][3]) {
  double cf[3] = {0, 0, 0};
  double cm[3] = {0, 0, 0};
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < 3; k++) {
      cf[k] += fixed[i][k];
      cm[k] += mobile[i][k];
    }
  }
  const double inv = 1.0 / static_cast<double>(n);
  for (int k = 0; k < 3; k++) {
    cf[k] *= inv;
    cm[k] *= inv;
  }
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    double y[3] = {mobile[i][0] - cm[0], mobile[i][1] - cm[1], mobile[i][2] - cm[2]};
    double ry[3];
    applyR(R, y, ry);
    const double dx = (fixed[i][0] - cf[0]) - ry[0];
    const double dy = (fixed[i][1] - cf[1]) - ry[1];
    const double dz = (fixed[i][2] - cf[2]) - ry[2];
    sum += dx * dx + dy * dy + dz * dz;
  }
  return std::sqrt(sum / static_cast<double>(n));
}

// Proper and improper. R maps mobile onto fixed: fixed ≈ R * (mobile - cm) + cf.
bool kabschEither(const double fixed[][3], const double mobile[][3], int n, double R[3][3],
                  double &rmsd) {
  double bestR[3][3];
  double best = 1e300;
  bool ok = false;
  double trial[kMaxN][3];
  for (int mirror = 0; mirror < 2; mirror++) {
    for (int i = 0; i < n; i++) {
      trial[i][0] = mobile[i][0];
      trial[i][1] = mobile[i][1];
      trial[i][2] = mirror ? -mobile[i][2] : mobile[i][2];
    }
    double Rm[3][3];
    double rms = 0.0;
    if (!kabschQuat(fixed, trial, n, Rm, rms)) {
      continue;
    }
    const double measured = rmsdOf(fixed, trial, n, Rm);
    if (measured < best) {
      best = measured;
      ok = true;
      if (mirror) {
        for (int i = 0; i < 3; i++) {
          bestR[i][0] = Rm[i][0];
          bestR[i][1] = Rm[i][1];
          bestR[i][2] = -Rm[i][2];
        }
      } else {
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            bestR[i][j] = Rm[i][j];
          }
        }
      }
    }
  }
  if (!ok) {
    return false;
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      R[i][j] = bestR[i][j];
    }
  }
  rmsd = best;
  return true;
}

void translationOf(const double fixed[][3], const double mobile[][3], int n, const double R[3][3],
                   double t[3]) {
  double cf[3] = {0, 0, 0};
  double cm[3] = {0, 0, 0};
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < 3; k++) {
      cf[k] += fixed[i][k];
      cm[k] += mobile[i][k];
    }
  }
  const double inv = 1.0 / static_cast<double>(n);
  for (int k = 0; k < 3; k++) {
    cf[k] *= inv;
    cm[k] *= inv;
  }
  double rcm[3];
  applyR(R, cm, rcm);
  for (int k = 0; k < 3; k++) {
    t[k] = cf[k] - rcm[k];
  }
}

double dist2(const double a[3], const double b[3]) {
  const double dx = a[0] - b[0];
  const double dy = a[1] - b[1];
  const double dz = a[2] - b[2];
  return dx * dx + dy * dy + dz * dz;
}

int hubOf(const double p[][3], int n) {
  int hub = 0;
  double best = 1e300;
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      sum += std::sqrt(dist2(p[i], p[j]));
    }
    if (sum < best) {
      best = sum;
      hub = i;
    }
  }
  return hub;
}

void nearest4(const double p[][3], int n, int hub, int out[4]) {
  int idx[kMaxN];
  double d[kMaxN];
  int m = 0;
  for (int j = 0; j < n; j++) {
    if (j == hub) {
      continue;
    }
    idx[m] = j;
    d[m] = dist2(p[hub], p[j]);
    m++;
  }
  for (int a = 0; a < 4 && a < m; a++) {
    int best = a;
    for (int b = a + 1; b < m; b++) {
      if (d[b] < d[best]) {
        best = b;
      }
    }
    std::swap(d[a], d[best]);
    std::swap(idx[a], idx[best]);
    out[a] = idx[a];
  }
}

int eclipsedCountPts(const double p[][3], int n, int hub, const int first[4],
                     double cosine, double *stored) {
  int count = 0;
  for (int a = 0; a < 4; a++) {
    const int ia = first[a];
    double u[3] = {p[ia][0] - p[hub][0], p[ia][1] - p[hub][1], p[ia][2] - p[hub][2]};
    const double un = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    if (!(un > 1e-8)) {
      continue;
    }
    u[0] /= un;
    u[1] /= un;
    u[2] /= un;
    double w[3][3];
    int nw = 0;
    for (int b = 0; b < 4; b++) {
      if (b == a) {
        continue;
      }
      const int ib = first[b];
      double v[3] = {p[ib][0] - p[hub][0], p[ib][1] - p[hub][1], p[ib][2] - p[hub][2]};
      const double dot = v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
      v[0] -= dot * u[0];
      v[1] -= dot * u[1];
      v[2] -= dot * u[2];
      const double vn = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      if (vn < 1e-8) {
        continue;
      }
      w[nw][0] = v[0] / vn;
      w[nw][1] = v[1] / vn;
      w[nw][2] = v[2] / vn;
      nw++;
    }
    double outer[6][3];
    int no = 0;
    for (int j = 0; j < n && no < 6; j++) {
      if (j == hub || j == ia) {
        continue;
      }
      const double bond = dist2(p[ia], p[j]);
      if (bond < kBondLo * kBondLo || bond > kBondHi * kBondHi) {
        continue;
      }
      double v[3] = {p[j][0] - p[ia][0], p[j][1] - p[ia][1], p[j][2] - p[ia][2]};
      const double dot = v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
      v[0] -= dot * u[0];
      v[1] -= dot * u[1];
      v[2] -= dot * u[2];
      const double vn = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      if (vn < 1e-8) {
        continue;
      }
      outer[no][0] = v[0] / vn;
      outer[no][1] = v[1] / vn;
      outer[no][2] = v[2] / vn;
      no++;
    }
    double mx = -1.0;
    for (int i = 0; i < nw; i++) {
      for (int j = 0; j < no; j++) {
        const double c =
            w[i][0] * outer[j][0] + w[i][1] * outer[j][1] + w[i][2] * outer[j][2];
        if (c > mx) {
          mx = c;
        }
      }
    }
    if (stored != nullptr) {
      stored[a] = mx;
    }
    if (mx > cosine) {
      count++;
    }
  }
  return count;
}

symm::Kind iceKindFromEclipsed(int nEcl) {
  if (nEcl == 0) {
    return symm::Kind::iceIc;
  }
  if (nEcl == 1) {
    return symm::Kind::iceIh;
  }
  return symm::Kind::iceIsd;
}

bool loadPoints(const Eigen::MatrixXd &m, double p[][3], int &n) {
  if (m.cols() != 3 || m.rows() < 1 || m.rows() > kMaxN) {
    return false;
  }
  n = static_cast<int>(m.rows());
  for (int i = 0; i < n; i++) {
    p[i][0] = m(i, 0);
    p[i][1] = m(i, 1);
    p[i][2] = m(i, 2);
  }
  return true;
}

void loadEigen(const Eigen::MatrixXd &m, double p[][3]) {
  const int n = static_cast<int>(m.rows());
  for (int i = 0; i < n; i++) {
    p[i][0] = m(i, 0);
    p[i][1] = m(i, 1);
    p[i][2] = m(i, 2);
  }
}

bool matchRing(const double env[][3], int n, const double templ[][3], int perm[kMaxN], double R[3][3],
               double t[3], double &rmsd) {
  double cf[3] = {0, 0, 0};
  for (int i = 0; i < n; i++) {
    cf[0] += env[i][0];
    cf[1] += env[i][1];
    cf[2] += env[i][2];
  }
  for (int k = 0; k < 3; k++) {
    cf[k] /= static_cast<double>(n);
  }
  double normal[3] = {0, 0, 0};
  for (int i = 0; i < n; i++) {
    const int j = (i + 1) % n;
    // Newell on the raw order is not cyclic-invariant. Use the two longest
    // axes of the covariance instead.
    (void)j;
  }
  double cov[3][3] = {};
  for (int i = 0; i < n; i++) {
    const double d[3] = {env[i][0] - cf[0], env[i][1] - cf[1], env[i][2] - cf[2]};
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        cov[a][b] += d[a] * d[b];
      }
    }
  }
  // Smallest eigenvector of cov is the plane normal. Power-iterate on the
  // adjugate (cofactors) which points along the smallest axis for a planar set.
  normal[0] = cov[1][1] * cov[2][2] - cov[1][2] * cov[2][1];
  normal[1] = cov[0][2] * cov[2][1] - cov[0][0] * cov[2][2];
  normal[2] = cov[0][0] * cov[1][1] - cov[0][1] * cov[1][0];
  // That is only one cofactor column. Use the cross of the two largest
  // centered vectors instead, which is stable for a polygon.
  int i0 = 0;
  int i1 = 1;
  double far = -1.0;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      const double d = dist2(env[i], env[j]);
      if (d > far) {
        far = d;
        i0 = i;
        i1 = j;
      }
    }
  }
  double bestArea = -1.0;
  double a[3] = {0, 0, 0};
  double b[3] = {0, 0, 0};
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      const double ai[3] = {env[i][0] - cf[0], env[i][1] - cf[1], env[i][2] - cf[2]};
      const double bj[3] = {env[j][0] - cf[0], env[j][1] - cf[1], env[j][2] - cf[2]};
      const double cx = ai[1] * bj[2] - ai[2] * bj[1];
      const double cy = ai[2] * bj[0] - ai[0] * bj[2];
      const double cz = ai[0] * bj[1] - ai[1] * bj[0];
      const double area = cx * cx + cy * cy + cz * cz;
      if (area > bestArea) {
        bestArea = area;
        a[0] = ai[0];
        a[1] = ai[1];
        a[2] = ai[2];
        b[0] = bj[0];
        b[1] = bj[1];
        b[2] = bj[2];
        normal[0] = cx;
        normal[1] = cy;
        normal[2] = cz;
      }
    }
  }
  (void)a;
  (void)b;
  const double nn = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
  if (!(nn > 1e-12)) {
    return false;
  }
  normal[0] /= nn;
  normal[1] /= nn;
  normal[2] /= nn;
  double axis[3] = {1, 0, 0};
  if (std::abs(normal[0]) > 0.9) {
    axis[0] = 0;
    axis[1] = 1;
  }
  double e1[3] = {axis[1] * normal[2] - axis[2] * normal[1],
                  axis[2] * normal[0] - axis[0] * normal[2],
                  axis[0] * normal[1] - axis[1] * normal[0]};
  const double e1n = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  e1[0] /= e1n;
  e1[1] /= e1n;
  e1[2] /= e1n;
  double e2[3] = {normal[1] * e1[2] - normal[2] * e1[1], normal[2] * e1[0] - normal[0] * e1[2],
                  normal[0] * e1[1] - normal[1] * e1[0]};
  struct Ang {
    double a;
    int i;
  } order[kMaxN];
  for (int i = 0; i < n; i++) {
    const double d[3] = {env[i][0] - cf[0], env[i][1] - cf[1], env[i][2] - cf[2]};
    const double x = d[0] * e1[0] + d[1] * e1[1] + d[2] * e1[2];
    const double y = d[0] * e2[0] + d[1] * e2[1] + d[2] * e2[2];
    order[i].a = std::atan2(y, x);
    order[i].i = i;
  }
  std::sort(order, order + n, [](const Ang &p, const Ang &q) { return p.a < q.a; });
  double best = 1e300;
  int bestPerm[kMaxN];
  double bestR[3][3];
  for (int dir = 0; dir < 2; dir++) {
    for (int shift = 0; shift < n; shift++) {
      int permTry[kMaxN];
      double mob[kMaxN][3];
      for (int step = 0; step < n; step++) {
        const int slot = dir == 0 ? (shift + step) % n : (shift + n - step) % n;
        permTry[step] = order[slot].i;
        mob[step][0] = env[permTry[step]][0];
        mob[step][1] = env[permTry[step]][1];
        mob[step][2] = env[permTry[step]][2];
      }
      double Rtry[3][3];
      double rms = 0.0;
      if (!kabschEither(templ, mob, n, Rtry, rms)) {
        continue;
      }
      const double measured = rmsdOf(templ, mob, n, Rtry);
      if (measured < best) {
        best = measured;
        rmsd = measured;
        for (int step = 0; step < n; step++) {
          bestPerm[step] = permTry[step];
        }
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            bestR[i][j] = Rtry[i][j];
          }
        }
      }
      // A ring has one cyclic class. The first overlay inside the noise is enough.
      if (best < 0.06) {
        dir = 2;
        break;
      }
    }
  }
  if (!(best < 1.0)) {
    return false;
  }
  for (int step = 0; step < n; step++) {
    perm[step] = bestPerm[step];
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      R[i][j] = bestR[i][j];
    }
  }
  double mob[kMaxN][3];
  for (int step = 0; step < n; step++) {
    mob[step][0] = env[perm[step]][0];
    mob[step][1] = env[perm[step]][1];
    mob[step][2] = env[perm[step]][2];
  }
  translationOf(templ, mob, n, R, t);
  rmsd = best;
  return true;
}

bool matchIce(const double env[][3], int n, const double templ[][3], int perm[kMaxN], double R[3][3],
              double t[3], double &rmsd, symm::Kind &kind, bool gateEclipsed = true) {
  const int hub = hubOf(env, n);
  int first[4];
  nearest4(env, n, hub, first);
  const int nEcl = eclipsedCountPts(env, n, hub, first, kEclipsed, nullptr);
  if (gateEclipsed && (nEcl < 0 || nEcl > 2)) {
    return false;
  }
  if (nEcl >= 0 && nEcl <= 2) {
    kind = iceKindFromEclipsed(nEcl);
  }
  int tf[4];
  nearest4(templ, n, 0, tf);
  // Template atom 0 is the centre by construction. tf should be 1..4.
  double best = 1e300;
  double best5 = 1e300;
  int bestPerm[kMaxN];
  double bestR[3][3];
  int order[4] = {0, 1, 2, 3};
  do {
    int permTry[kMaxN];
    for (int i = 0; i < n; i++) {
      permTry[i] = -1;
    }
    permTry[0] = hub;
    bool clash = false;
    for (int k = 0; k < 4; k++) {
      const int ti = tf[k];
      const int ej = first[order[k]];
      if (permTry[ti] != -1) {
        clash = true;
        break;
      }
      permTry[ti] = ej;
    }
    if (clash) {
      continue;
    }
    double fixed5[5][3];
    double mob5[5][3];
    int usedT[5];
    usedT[0] = 0;
    for (int k = 0; k < 4; k++) {
      usedT[k + 1] = tf[k];
    }
    for (int s = 0; s < 5; s++) {
      const int ti = usedT[s];
      const int ej = permTry[ti];
      for (int c = 0; c < 3; c++) {
        fixed5[s][c] = templ[ti][c];
        mob5[s][c] = env[ej][c];
      }
    }
    double R5[3][3];
    double rms5 = 0.0;
    if (!kabschEither(fixed5, mob5, 5, R5, rms5)) {
      continue;
    }
    best5 = std::min(best5, rmsdOf(fixed5, mob5, 5, R5));
    double cf[3] = {0, 0, 0};
    double cm[3] = {0, 0, 0};
    for (int s = 0; s < 5; s++) {
      for (int c = 0; c < 3; c++) {
        cf[c] += fixed5[s][c];
        cm[c] += mob5[s][c];
      }
    }
    for (int c = 0; c < 3; c++) {
      cf[c] /= 5.0;
      cm[c] /= 5.0;
    }
    bool taken[kMaxN] = {};
    for (int s = 0; s < 5; s++) {
      taken[permTry[usedT[s]]] = true;
    }
    for (int ti = 0; ti < n; ti++) {
      if (permTry[ti] != -1) {
        continue;
      }
      double q[3] = {templ[ti][0] - cf[0], templ[ti][1] - cf[1], templ[ti][2] - cf[2]};
      // Invert R: since R may be improper, solve by applying R^T if proper.
      // Compare in the mobile frame: R^T * q + cm.
      double rtq[3] = {R5[0][0] * q[0] + R5[1][0] * q[1] + R5[2][0] * q[2],
                       R5[0][1] * q[0] + R5[1][1] * q[1] + R5[2][1] * q[2],
                       R5[0][2] * q[0] + R5[1][2] * q[1] + R5[2][2] * q[2]};
      rtq[0] += cm[0];
      rtq[1] += cm[1];
      rtq[2] += cm[2];
      double bestD = 1e300;
      int bestE = -1;
      for (int ej = 0; ej < n; ej++) {
        if (taken[ej]) {
          continue;
        }
        const double d = dist2(rtq, env[ej]);
        if (d < bestD) {
          bestD = d;
          bestE = ej;
        }
      }
      if (bestE < 0) {
        clash = true;
        break;
      }
      permTry[ti] = bestE;
      taken[bestE] = true;
    }
    if (clash) {
      continue;
    }
    double mob[kMaxN][3];
    double fix[kMaxN][3];
    for (int ti = 0; ti < n; ti++) {
      for (int c = 0; c < 3; c++) {
        fix[ti][c] = templ[ti][c];
        mob[ti][c] = env[permTry[ti]][c];
      }
    }
    double Rall[3][3];
    double rms = 0.0;
    if (!kabschEither(fix, mob, n, Rall, rms)) {
      continue;
    }
    const double measured = rmsdOf(fix, mob, n, Rall);
    if (measured < best) {
      best = measured;
      for (int ti = 0; ti < n; ti++) {
        bestPerm[ti] = permTry[ti];
      }
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          bestR[i][j] = Rall[i][j];
        }
      }
    }
  } while (std::next_permutation(order, order + 4));
  (void)best5;
  rmsd = best;
  if (!(best < 1.0)) {
    return false;
  }
  for (int i = 0; i < n; i++) {
    perm[i] = bestPerm[i];
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      R[i][j] = bestR[i][j];
    }
  }
  double mob[kMaxN][3];
  for (int i = 0; i < n; i++) {
    mob[i][0] = env[perm[i]][0];
    mob[i][1] = env[perm[i]][1];
    mob[i][2] = env[perm[i]][2];
  }
  translationOf(templ, mob, n, R, t);
  rmsd = best;
  return true;
}

void adjacency(const double p[][3], int n, std::array<int, 3> adj[], int &degFail) {
  degFail = 0;
  for (int i = 0; i < n; i++) {
    int m = 0;
    int nbr[8];
    double nd[8];
    for (int j = 0; j < n; j++) {
      if (i == j) {
        continue;
      }
      const double d = dist2(p[i], p[j]);
      if (d > kBondHi * kBondHi) {
        continue;
      }
      if (m < 8) {
        nbr[m] = j;
        nd[m] = d;
        m++;
      }
    }
    for (int a = 0; a < 3 && a < m; a++) {
      int best = a;
      for (int b = a + 1; b < m; b++) {
        if (nd[b] < nd[best]) {
          best = b;
        }
      }
      std::swap(nd[a], nd[best]);
      std::swap(nbr[a], nbr[best]);
    }
    if (m < 3) {
      degFail++;
      adj[i] = {-1, -1, -1};
      continue;
    }
    adj[i] = {nbr[0], nbr[1], nbr[2]};
  }
}

bool cageRecurse(const std::array<int, 3> adjE[], const std::array<int, 3> adjT[], int n,
                 int mapE[], int mapT[], int depth);

bool cageAssign(const std::array<int, 3> adjE[], const std::array<int, 3> adjT[], int n, int e,
                int t, int mapE[], int mapT[], int depth) {
  if (mapE[e] == t) {
    return cageRecurse(adjE, adjT, n, mapE, mapT, depth);
  }
  if (mapE[e] != -1 || mapT[t] != -1) {
    return false;
  }
  mapE[e] = t;
  mapT[t] = e;
  const bool ok = cageRecurse(adjE, adjT, n, mapE, mapT, depth + 1);
  if (!ok) {
    mapE[e] = -1;
    mapT[t] = -1;
  }
  return ok;
}

bool cageRecurse(const std::array<int, 3> adjE[], const std::array<int, 3> adjT[], int n,
                 int mapE[], int mapT[], int depth) {
  if (depth == n) {
    return true;
  }
  int e = -1;
  int freeN = 4;
  for (int i = 0; i < n; i++) {
    if (mapE[i] == -1) {
      continue;
    }
    int open = 0;
    for (int k = 0; k < 3; k++) {
      if (adjE[i][static_cast<std::size_t>(k)] >= 0 &&
          mapE[adjE[i][static_cast<std::size_t>(k)]] == -1) {
        open++;
      }
    }
    if (open > 0 && open < freeN) {
      freeN = open;
      e = i;
    }
  }
  if (e < 0) {
    for (int i = 0; i < n; i++) {
      if (mapE[i] == -1) {
        e = i;
        break;
      }
    }
    if (e < 0) {
      return true;
    }
    for (int t = 0; t < n; t++) {
      if (mapT[t] != -1) {
        continue;
      }
      if (cageAssign(adjE, adjT, n, e, t, mapE, mapT, depth)) {
        return true;
      }
    }
    return false;
  }
  const int t = mapE[e];
  int openE[3];
  int nOpenE = 0;
  for (int k = 0; k < 3; k++) {
    if (adjE[e][static_cast<std::size_t>(k)] >= 0 &&
        mapE[adjE[e][static_cast<std::size_t>(k)]] == -1) {
      openE[nOpenE++] = adjE[e][static_cast<std::size_t>(k)];
    }
  }
  int openT[3];
  int nOpenT = 0;
  for (int k = 0; k < 3; k++) {
    if (adjT[t][static_cast<std::size_t>(k)] >= 0 &&
        mapT[adjT[t][static_cast<std::size_t>(k)]] == -1) {
      openT[nOpenT++] = adjT[t][static_cast<std::size_t>(k)];
    }
  }
  if (nOpenE != nOpenT) {
    return false;
  }
  int order[3] = {0, 1, 2};
  do {
    bool ok = true;
    for (int i = 0; i < nOpenE; i++) {
      if (!cageAssign(adjE, adjT, n, openE[i], openT[order[i]], mapE, mapT, depth)) {
        ok = false;
        break;
      }
    }
    if (ok && cageRecurse(adjE, adjT, n, mapE, mapT, depth)) {
      return true;
    }
    for (int i = 0; i < nOpenE; i++) {
      const int ee = openE[i];
      const int tt = mapE[ee];
      if (tt >= 0 && mapT[tt] == ee) {
        // undo only the ones this call mapped, if the recurse failed after partial assign
      }
    }
  } while (std::next_permutation(order, order + nOpenE));
  return false;
}

bool matchCage(const double env[][3], int n, const double templ[][3], int perm[kMaxN], double R[3][3],
               double t[3], double &rmsd) {
  std::array<int, 3> adjE[kMaxN];
  std::array<int, 3> adjT[kMaxN];
  int failE = 0;
  int failT = 0;
  adjacency(env, n, adjE, failE);
  adjacency(templ, n, adjT, failT);
  if (failE > 0 || failT > 0) {
    return false;
  }
  double best = 1e300;
  int bestPerm[kMaxN];
  double bestR[3][3];
  for (int e0 = 0; e0 < n; e0++) {
    int neigh[3] = {adjE[e0][0], adjE[e0][1], adjE[e0][2]};
    int order[3] = {0, 1, 2};
    const int tNeigh[3] = {adjT[0][0], adjT[0][1], adjT[0][2]};
    do {
      int mapE[kMaxN];
      int mapT[kMaxN];
      for (int i = 0; i < n; i++) {
        mapE[i] = -1;
        mapT[i] = -1;
      }
      mapE[e0] = 0;
      mapT[0] = e0;
      bool seeded = true;
      for (int k = 0; k < 3; k++) {
        const int ee = neigh[order[k]];
        const int tt = tNeigh[k];
        if (mapE[ee] != -1 || mapT[tt] != -1) {
          seeded = false;
          break;
        }
        mapE[ee] = tt;
        mapT[tt] = ee;
      }
      if (!seeded) {
        continue;
      }
      if (!cageRecurse(adjE, adjT, n, mapE, mapT, 4)) {
        continue;
      }
      int permTry[kMaxN];
      bool bad = false;
      for (int ti = 0; ti < n; ti++) {
        if (mapT[ti] < 0) {
          bad = true;
          break;
        }
        permTry[ti] = mapT[ti];
      }
      if (bad) {
        continue;
      }
      double mob[kMaxN][3];
      for (int ti = 0; ti < n; ti++) {
        mob[ti][0] = env[permTry[ti]][0];
        mob[ti][1] = env[permTry[ti]][1];
        mob[ti][2] = env[permTry[ti]][2];
      }
      double Rtry[3][3];
      double rms = 0.0;
      if (!kabschEither(templ, mob, n, Rtry, rms)) {
        continue;
      }
      const double measured = rmsdOf(templ, mob, n, Rtry);
      if (measured < best) {
        best = measured;
        for (int ti = 0; ti < n; ti++) {
          bestPerm[ti] = permTry[ti];
        }
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            bestR[i][j] = Rtry[i][j];
          }
        }
      }
      // One graph isomorphism is the cage. Further starts repeat it.
      if (best < 0.08) {
        e0 = n;
        break;
      }
    } while (std::next_permutation(order, order + 3));
  }
  if (!(best < 0.5)) {
    return false;
  }
  for (int i = 0; i < n; i++) {
    perm[i] = bestPerm[i];
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      R[i][j] = bestR[i][j];
    }
  }
  double mob[kMaxN][3];
  for (int i = 0; i < n; i++) {
    mob[i][0] = env[perm[i]][0];
    mob[i][1] = env[perm[i]][1];
    mob[i][2] = env[perm[i]][2];
  }
  translationOf(templ, mob, n, R, t);
  rmsd = best;
  return true;
}

void fillOverlay(symm::Overlay &out, symm::Kind kind, const double R[3][3], const double t[3],
                 const int perm[], int n, double rmsd) {
  out.status = 0;
  out.kind = kind;
  out.name = symm::name(kind);
  out.rmsd = rmsd;
  out.hausdorff = rmsd;
  out.backend = "kernel";
  out.rotation.setZero();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      out.rotation(i, j) = R[i][j];
    }
  }
  out.translation = Eigen::Vector3d(t[0], t[1], t[2]);
  out.permutation.assign(perm, perm + n);
}

struct Grid {
  int nx = 1;
  int ny = 1;
  int nz = 1;
  double ax[3]{};
  double ay[3]{};
  double az[3]{};
  double inv[3][3]{};
  std::vector<int> head;
  std::vector<int> next;
  std::vector<double> px;
  std::vector<double> py;
  std::vector<double> pz;
};

} // namespace

namespace symm {

Overlay matchKernel(const Eigen::MatrixXd &environment) {
  Overlay out;
  double env[kMaxN][3];
  int n = 0;
  if (!loadPoints(environment, env, n)) {
    return out;
  }
  const Kind kinds[] = {Kind::iceIh, Kind::iceIc, Kind::iceIsd, Kind::cage512,
                        Kind::cage51262, Kind::cage51264, Kind::ring5, Kind::ring6};
  Kind chosen = Kind::iceIc;
  bool found = false;
  for (Kind kind : kinds) {
    const Eigen::MatrixXd ref = points(kind);
    if (ref.rows() != n) {
      continue;
    }
    found = true;
    chosen = kind;
    break;
  }
  if (!found) {
    return out;
  }
  // Ice shells share one atom count. The eclipsed count selects which.
  if (n == 17) {
    const int hub = hubOf(env, n);
    int first[4];
    nearest4(env, n, hub, first);
    const int nEcl = eclipsedCountPts(env, n, hub, first, kEclipsed, nullptr);
    if (nEcl >= 0 && nEcl <= 2) {
      chosen = iceKindFromEclipsed(nEcl);
    }
  }
  double templ[kMaxN][3];
  loadEigen(points(chosen), templ);
  int perm[kMaxN];
  double R[3][3];
  double t[3];
  double rmsd = 0.0;
  bool ok = false;
  if (n == 5 || n == 6) {
    ok = matchRing(env, n, templ, perm, R, t, rmsd);
  } else if (n == 17) {
    Kind kind = chosen;
    ok = matchIce(env, n, templ, perm, R, t, rmsd, kind);
    chosen = kind;
  } else if (n == 20 || n == 24 || n == 28) {
    ok = matchCage(env, n, templ, perm, R, t, rmsd);
  }
  if (!ok) {
    return out;
  }
  fillOverlay(out, chosen, R, t, perm, n, rmsd);
  double haus = 0.0;
  for (int i = 0; i < n; i++) {
    double y[3] = {env[perm[i]][0], env[perm[i]][1], env[perm[i]][2]};
    double ry[3];
    applyR(R, y, ry);
    const double dx = templ[i][0] - (ry[0] + t[0]);
    const double dy = templ[i][1] - (ry[1] + t[1]);
    const double dz = templ[i][2] - (ry[2] + t[2]);
    haus = std::max(haus, std::sqrt(dx * dx + dy * dy + dz * dz));
  }
  out.hausdorff = haus;
  return out;
}

int bondCosines(const Eigen::MatrixXd &environment, double out[4]) {
  double env[kMaxN][3];
  int n = 0;
  if (!loadPoints(environment, env, n)) {
    return -1;
  }
  const int hub = hubOf(env, n);
  int first[4];
  nearest4(env, n, hub, first);
  for (int k = 0; k < 4; k++) {
    out[k] = -2.0;
  }
  const int nBond = eclipsedCountPts(env, n, hub, first, 2.0, out);
  (void)nBond;
  int nFound = 0;
  for (int k = 0; k < 4; k++) {
    if (out[k] > -1.5) {
      nFound++;
    }
  }
  std::sort(out, out + 4, [](double a, double b) { return a > b; });
  return nFound;
}

double iceTemplateRmsd(const Eigen::MatrixXd &environment, Kind kind) {
  if (kind != Kind::iceIh && kind != Kind::iceIc && kind != Kind::iceIsd) {
    return -1.0;
  }
  double env[kMaxN][3];
  int n = 0;
  if (!loadPoints(environment, env, n) || n != 17) {
    return -1.0;
  }
  double templ[kMaxN][3];
  loadEigen(points(kind), templ);
  int perm[kMaxN];
  double R[3][3];
  double t[3];
  double rmsd = 1e300;
  Kind named = kind;
  (void)matchIce(env, n, templ, perm, R, t, rmsd, named, false);
  if (!(rmsd < 10.0)) {
    return -1.0;
  }
  return rmsd;
}

double shapeMeasure(const Eigen::MatrixXd &environment) {
  const Overlay hit = matchKernel(environment);
  if (hit.status != 0) {
    return 100.0;
  }
  const Eigen::MatrixXd ref = points(hit.kind);
  const int n = static_cast<int>(ref.rows());
  double gp = 0.0;
  double gq = 0.0;
  double dot = 0.0;
  Eigen::Vector3d cf = Eigen::Vector3d::Zero();
  Eigen::Vector3d cm = Eigen::Vector3d::Zero();
  for (int i = 0; i < n; i++) {
    cf += ref.row(i).transpose();
    const int j = hit.permutation[static_cast<size_t>(i)];
    cm += environment.row(j).transpose();
  }
  cf /= static_cast<double>(n);
  cm /= static_cast<double>(n);
  for (int i = 0; i < n; i++) {
    const Eigen::Vector3d p = ref.row(i).transpose() - cf;
    const int j = hit.permutation[static_cast<size_t>(i)];
    const Eigen::Vector3d q = hit.rotation * (environment.row(j).transpose() - cm);
    gp += p.squaredNorm();
    gq += q.squaredNorm();
    dot += p.dot(q);
  }
  if (!(gp > 0.0) || !(gq > 0.0)) {
    return 100.0;
  }
  const double scale = dot / gq;
  const double residual = gp - 2.0 * scale * dot + scale * scale * gq;
  return 100.0 * std::max(0.0, residual) / gp;
}

// Even permutations of four bonds: the rotation group of the tetrahedron.
const int kEven4[12][4] = {
    {0, 1, 2, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 3, 2, 0},
    {2, 0, 1, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 2, 1, 0},
};

void loadOneShell(symm::Kind kind, double first[4][3], double second[12][3], int parent[12]) {
  double p[kMaxN][3];
  loadEigen(points(kind), p);
  int hub[4];
  nearest4(p, 17, 0, hub);
  bool used[17] = {};
  used[0] = true;
  for (int k = 0; k < 4; k++) {
    used[hub[k]] = true;
    for (int c = 0; c < 3; c++) {
      first[k][c] = p[hub[k]][c] - p[0][c];
    }
  }
  int ns = 0;
  for (int i = 1; i < 17 && ns < 12; i++) {
    if (used[i]) {
      continue;
    }
    for (int c = 0; c < 3; c++) {
      second[ns][c] = p[i][c] - p[0][c];
    }
    int owner = 0;
    double od = 1e300;
    for (int k = 0; k < 4; k++) {
      const double dx = second[ns][0] - first[k][0];
      const double dy = second[ns][1] - first[k][1];
      const double dz = second[ns][2] - first[k][2];
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 < od) {
        od = d2;
        owner = k;
      }
    }
    parent[ns] = owner;
    ns++;
  }
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

void prepareCanon(CanonShell &out) {
  const symm::Kind kinds[3] = {symm::Kind::iceIh, symm::Kind::iceIc, symm::Kind::iceIsd};
  for (int k = 0; k < 3; k++) {
    loadOneShell(kinds[k], out.first[k], out.second[k], out.parent[k]);
    double raw[kMaxN][3];
    loadEigen(points(kinds[k]), raw);
    for (int i = 0; i < 17; i++) {
      for (int c = 0; c < 3; c++) {
        out.pts[k][i][c] = raw[i][c];
      }
    }
    for (int a = 0; a < 4; a++) {
      out.nArm[k][a] = 0;
    }
    for (int s = 0; s < 12; s++) {
      const int a = out.parent[k][s];
      if (out.nArm[k][a] < 3) {
        out.arm[k][a][out.nArm[k][a]++] = s;
      }
    }
    for (int a = 0; a < 4; a++) {
      double e1[3];
      double e2[3];
      bondFrame(out.first[k][a], out.first[k][(a + 1) % 4], e1, e2);
      const int nA = out.nArm[k][a];
      for (int i = 0; i < nA; i++) {
        for (int j = i + 1; j < nA; j++) {
          const int si = out.arm[k][a][i];
          const int sj = out.arm[k][a][j];
          const double vi[3] = {out.second[k][si][0] - out.first[k][a][0],
                                out.second[k][si][1] - out.first[k][a][1],
                                out.second[k][si][2] - out.first[k][a][2]};
          const double vj[3] = {out.second[k][sj][0] - out.first[k][a][0],
                                out.second[k][sj][1] - out.first[k][a][1],
                                out.second[k][sj][2] - out.first[k][a][2]};
          const double ai = bondAngle(vi, e1, e2);
          const double aj = bondAngle(vj, e1, e2);
          if (aj < ai) {
            out.arm[k][a][i] = sj;
            out.arm[k][a][j] = si;
          }
        }
      }
    }
  }
}

// The canon shell does not depend on the frame. Each calling thread keeps
// its own copy. The first use on that thread fills it; a second thread does
// not write the first thread's tables.
CanonShell &canonShell() {
  thread_local CanonShell canon = [] {
    CanonShell out{};
    prepareCanon(out);
    return out;
  }();
  return canon;
}

double shellGate(double shellRmsd) {
  return std::isfinite(shellRmsd) && shellRmsd > 0.0 ? shellRmsd : kShellRmsd;
}

// The per-atom search lives in shell_classify.cpp, which does not include
// Eigen and is built with AVX2. This wrapper only copies the coordinates.
void classifyIce(const Eigen::MatrixXd &cartesian, const Eigen::Matrix3d &vectors,
                 std::vector<int> &kind, double shellRmsd, double thermalQ, double bondRatio) {
  const int n = static_cast<int>(cartesian.rows());
  kind.assign(static_cast<size_t>(std::max(n, 0)), -1);
  if (n <= 0 || cartesian.cols() != 3) {
    return;
  }
  const double lat[3][3] = {{vectors(0, 0), vectors(0, 1), vectors(0, 2)},
                            {vectors(1, 0), vectors(1, 1), vectors(1, 2)},
                            {vectors(2, 0), vectors(2, 1), vectors(2, 2)}};
  classifyPrepared(cartesian.col(0).data(), cartesian.col(1).data(), cartesian.col(2).data(), n, lat,
                   canonShell(), shellGate(shellRmsd), kind.data(), thermalQ, bondRatio);
}

void classifyIceCloud(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                      std::vector<int> &kind) {
  const int n = yCloud.nop;
  kind.assign(static_cast<std::size_t>(std::max(n, 0)), -1);
  if (n <= 0) {
    return;
  }
  std::vector<double> x(static_cast<std::size_t>(n));
  std::vector<double> y(x.size());
  std::vector<double> z(x.size());
  for (int i = 0; i < n; ++i) {
    const auto &point = yCloud.pts[static_cast<std::size_t>(i)];
    x[static_cast<std::size_t>(i)] = point.x;
    y[static_cast<std::size_t>(i)] = point.y;
    z[static_cast<std::size_t>(i)] = point.z;
  }
  double lat[3][3] = {};
  const auto &box = yCloud.box;
  lat[0][0] = box.size() > 0 ? box[0] : 1.0;
  lat[1][1] = box.size() > 1 ? box[1] : lat[0][0];
  lat[2][2] = box.size() > 2 ? box[2] : lat[0][0];
  if (box.size() >= 6) {
    lat[1][0] = box[3];
    lat[2][0] = box[4];
    lat[2][1] = box[5];
  }
  classifyIceArrays(x.data(), y.data(), z.data(), n, lat, kind.data());
}

void classifyIceArrays(const double *x, const double *y, const double *z, int n, const double lat[3][3],
                       int *kind, double shellRmsd, double thermalQ, double bondRatio) {
  if (n <= 0 || x == nullptr || y == nullptr || z == nullptr || kind == nullptr) {
    return;
  }
  classifyPrepared(x, y, z, n, lat, canonShell(), shellGate(shellRmsd), kind, thermalQ, bondRatio);
}

void matchAll(const std::vector<Eigen::MatrixXd> &environments, std::vector<Overlay> &out) {
  const int n = static_cast<int>(environments.size());
  out.resize(static_cast<size_t>(std::max(n, 0)));
#ifdef SEAMS_HAS_OPENMP
  int threads = 1;
  if (!omp_in_parallel()) {
    threads = omp_get_max_threads();
  }
#pragma omp parallel for schedule(static) if (n >= 32 && threads > 1)
#endif
  for (int i = 0; i < n; i++) {
    out[static_cast<size_t>(i)] = matchKernel(environments[static_cast<size_t>(i)]);
  }
}

} // namespace symm
