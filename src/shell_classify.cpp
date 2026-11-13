//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Per-atom ice label. The neighbour list is the orthorhombic cutoff CSR:
// cells no wider than the 3.5 A bond, and the 27 cells around each atom.
// One thread walks each pair once. A skewed or tiny box keeps the finer
// 2.8 A grid. The fit is a 4-point Horn quaternion on the ambiguous shell.
// Eigen is not included in this file, so the translation unit can be built
// with AVX2.

#include <shell_hot.hpp>
#include <thread_team.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <immintrin.h>
#include <new>
#include <cstddef>
#include <memory>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace symm {
namespace {

bool symTime() {
  // One initialization. Two classify calls must not race on the flag.
  static const int on = [] {
    const char *e = std::getenv("SEAMS_SYM_TIME");
    return (e != nullptr && e[0] != '\0' && e[0] != '0') ? 1 : 0;
  }();
  return on != 0;
}

double msSince(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

constexpr double kBond = 3.5;
// The common neighbour sits inside 3.35 A. The list still reaches 3.5 A when
// the four nearest are not all inside that, so a stretched crystal bond is kept.
constexpr double kNear = 3.35;
constexpr double kWide = 4.05;
constexpr double kCentre = 5.6;
constexpr double kCentreHex = 6.2;
constexpr double kGrid = 2.8;
constexpr int kMaxBond = 4;
constexpr int kMaxWide = 4;

struct Bond {
  int j;
  double dx;
  double dy;
  double dz;
};

double bondLength2(const Bond &b) {
  return b.dx * b.dx + b.dy * b.dy + b.dz * b.dz;
}

// Keep the closest bonds. Used when one atom records a pair for both ends.
void keepClosest(Bond *dest, int &nb, const Bond *extra, int nExtra) {
  Bond tmp[32];
  int nt = 0;
  for (int a = 0; a < nb && nt < 32; a++) {
    tmp[nt++] = dest[a];
  }
  for (int a = 0; a < nExtra && nt < 32; a++) {
    tmp[nt++] = extra[a];
  }
  for (int a = 1; a < nt; a++) {
    const Bond b = tmp[a];
    const double d2 = bondLength2(b);
    int at = a;
    while (at > 0 && d2 < bondLength2(tmp[at - 1])) {
      tmp[at] = tmp[at - 1];
      at--;
    }
    tmp[at] = b;
  }
  nb = nt < kMaxBond ? nt : kMaxBond;
  for (int a = 0; a < nb; a++) {
    dest[a] = tmp[a];
  }
}

void pushClosest(Bond *row, int &nb, int other, double dx, double dy, double dz, double d2) {
  if (nb < kMaxBond) {
    row[nb++] = Bond{other, dx, dy, dz};
    return;
  }
  int worst = 0;
  double wd = bondLength2(row[0]);
  for (int c = 1; c < kMaxBond; c++) {
    const double d = bondLength2(row[c]);
    if (d > wd) {
      wd = d;
      worst = c;
    }
  }
  if (d2 < wd) {
    row[worst] = Bond{other, dx, dy, dz};
  }
}

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

struct Packed;

// One workspace per calling thread, grown when a larger frame arrives and
// reused after that. Two classify calls may run at once. The OpenMP loops
// inside one call write this thread's buffers; they do not take another
// workspace.
struct Workspace {
  size_t nCap = 0;
  size_t cCap = 0;
  double *px = nullptr;
  double *py = nullptr;
  double *pz = nullptr;
  int *bin = nullptr;
  int *nbon = nullptr;
  int *cellStart = nullptr;
  int *fill = nullptr;
  Packed *packed = nullptr;
  Bond *bonds = nullptr;
  double *cx = nullptr;
  double *cy = nullptr;
  double *cz = nullptr;
  int *cid = nullptr;
  double *bd = nullptr;
  int *js = nullptr;
  double *vx = nullptr;
  double *vy = nullptr;
  double *vz = nullptr;
  int *rawK = nullptr;
  char *crystal = nullptr;
  char *tried = nullptr;

  ~Workspace() {
    ::operator delete(px);
    ::operator delete(py);
    ::operator delete(pz);
    ::operator delete(bin);
    ::operator delete(nbon);
    ::operator delete(cellStart);
    ::operator delete(fill);
    ::operator delete(packed);
    ::operator delete(bonds);
    ::operator delete(cx);
    ::operator delete(cy);
    ::operator delete(cz);
    ::operator delete(cid);
    ::operator delete(bd);
    ::operator delete(js);
    ::operator delete(vx);
    ::operator delete(vy);
    ::operator delete(vz);
    ::operator delete(rawK);
    ::operator delete(crystal);
    ::operator delete(tried);
  }

  template <class T>
  static T *renew(T *p, size_t n) {
    ::operator delete(p);
    return static_cast<T *>(::operator new(n * sizeof(T)));
  }

  void fit(size_t n, size_t ncell) {
    if (n > nCap) {
      px = renew(px, n);
      py = renew(py, n);
      pz = renew(pz, n);
      bin = renew(bin, n);
      nbon = renew(nbon, n);
      packed = renew(packed, n);
      bonds = renew(bonds, n * static_cast<size_t>(kMaxBond));
      cx = renew(cx, n);
      cy = renew(cy, n);
      cz = renew(cz, n);
      cid = renew(cid, n);
      bd = renew(bd, n * static_cast<size_t>(kMaxBond));
      js = renew(js, n * static_cast<size_t>(kMaxBond));
      // The pair walk stores the minimum-image vector beside each kept neighbour.
      vx = renew(vx, n * static_cast<size_t>(kMaxBond));
      vy = renew(vy, n * static_cast<size_t>(kMaxBond));
      vz = renew(vz, n * static_cast<size_t>(kMaxBond));
      rawK = renew(rawK, n);
      crystal = renew(crystal, n);
      tried = renew(tried, n);
      nCap = n;
    }
    if (ncell + 1 > cCap) {
      cellStart = renew(cellStart, ncell + 1);
      fill = renew(fill, ncell);
      cCap = ncell + 1;
    }
  }
};

Workspace &workspace() {
  thread_local Workspace mem;
  return mem;
}

template <class T>
struct View {
  T *p;
  T &operator[](size_t i) const { return p[i]; }
  T *get() const { return p; }
  T *data() const { return p; }
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
  // Two sweeps name the same held-out shells as three. The 4-point fit stays at one.
  jacobi4(K, Q, n <= 4 ? 1 : 2);
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
  // Only the four points are read. Zeroing the other thirteen used to sit on
  // every correspondence.
  Shell17 a;
  Shell17 b;
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

double tetraQ(const Bond *mine) {
  double u[4][3];
  for (int k = 0; k < 4; k++) {
    const double nrm = std::sqrt(mine[k].dx * mine[k].dx + mine[k].dy * mine[k].dy + mine[k].dz * mine[k].dz);
    if (!(nrm > 1e-8)) {
      return -1.0;
    }
    const double inv = 1.0 / nrm;
    u[k][0] = mine[k].dx * inv;
    u[k][1] = mine[k].dy * inv;
    u[k][2] = mine[k].dz * inv;
  }
  double acc = 0.0;
  for (int a = 0; a < 4; a++) {
    for (int b = a + 1; b < 4; b++) {
      const double c = u[a][0] * u[b][0] + u[a][1] * u[b][1] + u[a][2] * u[b][2];
      const double d = c + (1.0 / 3.0);
      acc += d * d;
    }
  }
  return 1.0 - 0.375 * acc;
}

double bondAngle(const double v[3], const double e1[3], const double e2[3]) {
  return std::atan2(v[0] * e2[0] + v[1] * e2[1] + v[2] * e2[2],
                    v[0] * e1[0] + v[1] * e1[1] + v[2] * e1[2]);
}

// True when one of the three grouped seconds on bond `a` is eclipsed.
// The cosine is compared as a square, so the bond and the reference are not normalised.
// Ideal eclipsed is 1 and ideal staggered is 0.5. Thermal staggered bonds stay at or below 0.90.
bool armEclipsed(const double envF[4][3], const double envS[12][3], int a) {
  const double ux = envF[a][0];
  const double uy = envF[a][1];
  const double uz = envF[a][2];
  const double ul2 = ux * ux + uy * uy + uz * uz;
  if (!(ul2 > 1e-16)) {
    return false;
  }
  constexpr double kEcl2 = 0.90 * 0.90;
  for (int b = 0; b < 4; b++) {
    if (b == a) {
      continue;
    }
    const double rd = (envF[b][0] * ux + envF[b][1] * uy + envF[b][2] * uz) / ul2;
    const double rx = envF[b][0] - rd * ux;
    const double ry = envF[b][1] - rd * uy;
    const double rz = envF[b][2] - rd * uz;
    const double rl2 = rx * rx + ry * ry + rz * rz;
    if (rl2 < 1e-12) {
      continue;
    }
    for (int e = 0; e < 3; e++) {
      const int s = 3 * a + e;
      const double vx0 = envS[s][0] - ux;
      const double vy0 = envS[s][1] - uy;
      const double vz0 = envS[s][2] - uz;
      const double vd = (vx0 * ux + vy0 * uy + vz0 * uz) / ul2;
      const double vx = vx0 - vd * ux;
      const double vy = vy0 - vd * uy;
      const double vz = vz0 - vd * uz;
      const double vl2 = vx * vx + vy * vy + vz * vz;
      if (vl2 < 1e-12) {
        continue;
      }
      const double dot = rx * vx + ry * vy + rz * vz;
      if (dot > 0.0 && dot * dot > kEcl2 * rl2 * vl2) {
        return true;
      }
    }
  }
  return false;
}

// Largest azimuthal cosine of the three second neighbours on bond `a`.
// The seconds are grouped, three per first-shell bond. Ideal eclipsed is 1
// and ideal staggered is 0.5. Thermal staggered bonds stay at or below 0.90.
double armMaxCosine(const double envF[4][3], const double envS[12][3], int a) {
  double u[3] = {envF[a][0], envF[a][1], envF[a][2]};
  const double un = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  if (!(un > 1e-8)) {
    return -1.0;
  }
  u[0] /= un;
  u[1] /= un;
  u[2] /= un;
  double mx = -1.0;
  for (int b = 0; b < 4; b++) {
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
    for (int e = 0; e < 3; e++) {
      const int s = 3 * a + e;
      double v[3] = {envS[s][0] - envF[a][0], envS[s][1] - envF[a][1], envS[s][2] - envF[a][2]};
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
  return mx;
}

// The hexagonal template has one eclipsed arm. The index is fixed for the
// canonical shell, so it is computed once. The initializer runs on one
// thread; a later call does not write it.
int templateEclipsedArm(const CanonShell &canon) {
  static const int cached = [&canon] {
    double envF[4][3];
    double envS[12][3];
    for (int a = 0; a < 4; a++) {
      for (int c = 0; c < 3; c++) {
        envF[a][c] = canon.first[0][a][c];
      }
      for (int e = 0; e < 3; e++) {
        const int s = canon.arm[0][a][e];
        for (int c = 0; c < 3; c++) {
          envS[3 * a + e][c] = canon.second[0][s][c];
        }
      }
    }
    int arm = 0;
    for (int a = 0; a < 4; a++) {
      if (armMaxCosine(envF, envS, a) > 0.90) {
        arm = a;
      }
    }
    return arm;
  }();
  return cached;
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

} // namespace

// Second-shell diamond test on the four neighbours, the same split OVITO's
// identifier and pyscal3 use: the twelve second neighbours of a cubic
// diamond are an FCC shell (every bond staggered) and those of a hexagonal
// diamond are an HCP shell (one bond eclipsed). A clean shell is a single
// radius. Anything wider falls through to the full fit.
// Maras, Trushin, Stukowski, Albe and Nordlund, Comput. Phys. Commun. 205, 13 (2016),
// Appendix A. The centre and its twelve second neighbours are tested as an FCC
// cluster (cubic diamond, ice Ic) or an HCP cluster (hexagonal diamond, ice Ih).
// Each pair uses the adaptive cutoff (1 + sqrt(2)) / 2 times that pair's length.
// A shell is not rejected because its radius varies. Twelve (4,2,1) signatures
// are cubic. Six (4,2,1) and six (4,2,2) are hexagonal. Eight (4,2,1), or
// eight with the hexagonal mix, is the partial label. Anything else with at
// least `ambiguous` of those signatures is worth the full fit.

int diamondCna(const double *sx, const double *sy, const double *sz, int ambiguous, int *n421o = nullptr,
               int *n422o = nullptr, int *nGoodo = nullptr) {
  alignas(32) double r2[16];
  double pair[12][12];
  alignas(32) double px[16];
  alignas(32) double py[16];
  alignas(32) double pz[16];
  for (int i = 0; i < 12; i++) {
    px[i] = sx[i];
    py[i] = sy[i];
    pz[i] = sz[i];
    r2[i] = sx[i] * sx[i] + sy[i] * sy[i] + sz[i] * sz[i];
  }
  for (int i = 12; i < 16; i++) {
    px[i] = 1.0e6;
    py[i] = 1.0e6;
    pz[i] = 1.0e6;
    r2[i] = 1.0e300;
  }
  for (int a = 0; a < 12; a++) {
    const __m256d xa = _mm256_set1_pd(px[a]);
    const __m256d ya = _mm256_set1_pd(py[a]);
    const __m256d za = _mm256_set1_pd(pz[a]);
    for (int b = 0; b < 12; b += 4) {
      const __m256d dx = _mm256_sub_pd(xa, _mm256_load_pd(px + b));
      const __m256d dy = _mm256_sub_pd(ya, _mm256_load_pd(py + b));
      const __m256d dz = _mm256_sub_pd(za, _mm256_load_pd(pz + b));
      const __m256d d2 = _mm256_add_pd(_mm256_mul_pd(dx, dx),
                                      _mm256_add_pd(_mm256_mul_pd(dy, dy), _mm256_mul_pd(dz, dz)));
      _mm256_storeu_pd(&pair[a][b], d2);
    }
  }
  // (1 + sqrt(2)) / 2, the adaptive common-neighbour cutoff. It lies between
  // the first and second shells of the FCC/HCP cluster.
  constexpr double kAdapt2 = 1.2071067811865475 * 1.2071067811865475;
  int n421 = 0;
  int n422 = 0;
  int nGood = 0;
  for (int j = 0; j < 12; j++) {
    const __m256d xa = _mm256_set1_pd(px[j]);
    const __m256d ya = _mm256_set1_pd(py[j]);
    const __m256d za = _mm256_set1_pd(pz[j]);
    for (int b = 0; b < 12; b += 4) {
      const __m256d dx = _mm256_sub_pd(xa, _mm256_load_pd(px + b));
      const __m256d dy = _mm256_sub_pd(ya, _mm256_load_pd(py + b));
      const __m256d dz = _mm256_sub_pd(za, _mm256_load_pd(pz + b));
      const __m256d d2 = _mm256_add_pd(_mm256_mul_pd(dx, dx),
                                      _mm256_add_pd(_mm256_mul_pd(dy, dy), _mm256_mul_pd(dz, dz)));
      _mm256_storeu_pd(&pair[j][b], d2);
    }
    const double cut2 = r2[j] * kAdapt2;
    const __m256d cut = _mm256_set1_pd(cut2);
    int mask = 0;
    for (int b = 0; b < 12; b += 4) {
      const __m256d pk = _mm256_loadu_pd(&pair[j][b]);
      const __m256d rk = _mm256_load_pd(r2 + b);
      const __m256d ok = _mm256_and_pd(_mm256_cmp_pd(pk, cut, _CMP_LT_OQ),
                                      _mm256_cmp_pd(rk, cut, _CMP_LT_OQ));
      mask |= _mm256_movemask_pd(ok) << b;
    }
    mask &= ~(1 << j);
    const int nc = __builtin_popcount(static_cast<unsigned>(mask));
    int common[4];
    if (nc == 4) {
      int slot = 0;
      for (int k = 0; k < 12; k++) {
        if (mask & (1 << k)) {
          common[slot++] = k;
        }
      }
    }
    bool good = false;
    if (nc == 4) {
      int nBonds = 0;
      int deg[4] = {};
      for (int a = 0; a < 4; a++) {
        for (int b = a + 1; b < 4; b++) {
          const int ia = common[a];
          const int ib = common[b];
          const double dx = sx[ia] - sx[ib];
          const double dy = sy[ia] - sy[ib];
          const double dz = sz[ia] - sz[ib];
          if (dx * dx + dy * dy + dz * dz < cut2) {
            nBonds++;
            deg[a]++;
            deg[b]++;
          }
        }
      }
      if (nBonds == 2) {
        int chain = 1;
        for (int a = 0; a < 4; a++) {
          if (deg[a] >= 2) {
            chain = 2;
          }
        }
        good = true;
        nGood++;
        if (chain == 1) {
          n421++;
        } else {
          n422++;
        }
      }
    }
    // Eight good signatures can still be a partial label. `ambiguous` is the
    // higher count that only earns a fit. Stop once neither is reachable.
    const int keep = ambiguous < 8 ? ambiguous : 8;
    if (!good && nGood + (11 - j) < keep) {
      return -1;
    }
  }
  if (n421o) {
    *n421o = n421;
  }
  if (n422o) {
    *n422o = n422;
  }
  if (nGoodo) {
    *nGoodo = nGood;
  }
  if (n421 == 12) {
    return 1;
  }
  if (n421 == 6 && n422 == 6) {
    return 0;
  }
  // Mostly one signature. The caller accepts it only for a sharp tetrahedron.
  if (nGood >= 8 && n421 >= 6 && n421 >= n422) {
    return 3;
  }
  if (nGood >= 8 && n422 >= 4 && n422 > n421) {
    return 4;
  }
  if (nGood >= ambiguous) {
    return -2;
  }
  return -1;
}

namespace {

int diamondSecondShell(const double first[4][3], const double *sx, const double *sy, const double *sz,
                      double radiusRatio, double eclipsedCos, double staggeredCos) {
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
  // A thermal crystal still has one shell; liquid does not.
  (void)radiusRatio;
  if (!(rmin > 1.0)) {
    return -1;
  }
  int nEclipsed = 0;
  int nStaggered = 0;
  // Project each arm's three second bonds once, then take the best cosine
  // against the other three bonds. eclipsed: cos >= 0.95. staggered: best
  // cos in (0.20, 0.65].
  const double kEcl2 = eclipsedCos * eclipsedCos;
  const double kStg2 = staggeredCos * staggeredCos;
  constexpr double kLow2 = 0.20 * 0.20;
  for (int a = 0; a < 4; a++) {
    const double u0 = first[a][0];
    const double u1 = first[a][1];
    const double u2 = first[a][2];
    const double un2 = u0 * u0 + u1 * u1 + u2 * u2;
    if (!(un2 > 1e-8)) {
      return -1;
    }
    double px[3], py[3], pz[3], vn2[3];
    for (int e = 0; e < 3; e++) {
      const int s = 3 * a + e;
      const double vx0 = sx[s] - u0;
      const double vy0 = sy[s] - u1;
      const double vz0 = sz[s] - u2;
      const double vd = vx0 * u0 + vy0 * u1 + vz0 * u2;
      px[e] = vx0 * un2 - vd * u0;
      py[e] = vy0 * un2 - vd * u1;
      pz[e] = vz0 * un2 - vd * u2;
      vn2[e] = px[e] * px[e] + py[e] * py[e] + pz[e] * pz[e];
    }
    bool eclipsed = false;
    bool low = true;
    bool any = false;
    for (int b = 0; b < 4 && !eclipsed; b++) {
      if (b == a) {
        continue;
      }
      const double rd = first[b][0] * u0 + first[b][1] * u1 + first[b][2] * u2;
      const double rx = first[b][0] * un2 - rd * u0;
      const double ry = first[b][1] * un2 - rd * u1;
      const double rz = first[b][2] * un2 - rd * u2;
      const double rn2 = rx * rx + ry * ry + rz * rz;
      if (rn2 < 1e-18) {
        continue;
      }
      for (int e = 0; e < 3; e++) {
        if (vn2[e] < 1e-18) {
          continue;
        }
        const double dot = rx * px[e] + ry * py[e] + rz * pz[e];
        if (dot <= 0.0) {
          continue;
        }
        const double dot2 = dot * dot;
        const double scale = rn2 * vn2[e];
        if (dot2 >= kEcl2 * scale) {
          eclipsed = true;
          break;
        }
        if (dot2 > kStg2 * scale) {
          low = false;
        } else if (dot2 >= kLow2 * scale) {
          any = true;
        }
      }
    }
    if (eclipsed) {
      nEclipsed++;
    } else if (low && any) {
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

// Ideal ice has equal bonds and one second-shell radius. The cosine here is
// the raw dot product of a second-shell bond with another first-shell bond.
// The fit uses the azimuthal cosine instead (ideal staggered 1/2, ideal
// eclipsed 1). Staggered raw cosine is 1/3 and eclipsed raw cosine is 7/9,
// so the projection is unnecessary.
__attribute__((always_inline)) inline int idealKind(const double first[4][3], const double *sx,
                                                            const double *sy, const double *sz) {
  // Squares of those raw cosines: eclipsed when the best square is at least
  // 0.36 (|cos| >= 0.6), staggered when it is at most 0.25 (|cos| <= 0.5).
  // 7/9 squared is 49/81 and 1/3 squared is 1/9. No square root on this path.
  constexpr double kEcl2 = 0.36;
  constexpr double kStg2 = 0.25;
  double b2[4];
  for (int k = 0; k < 4; k++) {
    b2[k] = first[k][0] * first[k][0] + first[k][1] * first[k][1] + first[k][2] * first[k][2];
  }
  int nEclipsed = 0;
  int nStaggered = 0;
  for (int a = 0; a < 4; a++) {
    double best = 0.0;
    for (int e = 0; e < 3; e++) {
      const int s = 3 * a + e;
      const double vx = sx[s] - first[a][0];
      const double vy = sy[s] - first[a][1];
      const double vz = sz[s] - first[a][2];
      const double vn2 = vx * vx + vy * vy + vz * vz;
      for (int b = 0; b < 4; b++) {
        if (b == a) {
          continue;
        }
        const double dot = vx * first[b][0] + vy * first[b][1] + vz * first[b][2];
        if (dot <= 0.0) {
          continue;
        }
        const double ratio = (dot * dot) / (vn2 * b2[b]);
        if (ratio > best) {
          best = ratio;
        }
      }
    }
    if (best >= kEcl2) {
      nEclipsed++;
    } else if (best <= kStg2) {
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
             bool recover = false, bool grouped = false, bool shallow = false) {
  const int m = std::min(nNb, 16);
  int order[16];
  double dd[16];
  for (int i = 0; i < m; i++) {
    order[i] = i;
    dd[i] = dist2(center, nb[i]);
  }
  // Grouped input is already the four bonds, then three seconds of each.
  if (!grouped) {
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
  // Bond-angle fingerprint. A shell whose four arms are all eclipsed is the
  // liquid pattern and does not enter the permutation search. A clean cubic
  // or hexagonal count names the one template, and the eclipsed arm fixes
  // which template bond that arm maps onto, so the other templates and the
  // other correspondences are not fitted.
  int onlyKind = -1;
  int anchorEnv = -1;
  int anchorTem = -1;
  if (grouped || recover) {
    int nEcl = 0;
    int eArm = -1;
    for (int a = 0; a < 4; a++) {
      double mx = -1.0;
      if (grouped && armEclipsed(envF, envS, a)) {
        mx = 1.0;
      }
      if (!grouped) {
        double u[3] = {envF[a][0], envF[a][1], envF[a][2]};
        const double un = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        if (!(un > 1e-8)) {
          continue;
        }
        u[0] /= un;
        u[1] /= un;
        u[2] /= un;
        for (int b = 0; b < 4 && mx < 0.90; b++) {
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
      }
      if (mx > 0.90) {
        nEcl++;
        eArm = a;
      }
    }
    if (nEcl == 4 && !recover) {
      return -1;
    }
    if (nEcl == 0) {
      onlyKind = 1;
    } else if (nEcl == 1 || (nEcl == 4 && recover)) {
      onlyKind = 0;
      if (nEcl == 4) {
        eArm = -1;
      }
    }
    if (nEcl == 1) {
      onlyKind = 0;
      anchorEnv = eArm;
      anchorTem = templateEclipsedArm(canon);
    }
  }
  double best = 1e300;
  int bestK = -1;
  const int kinds[3] = {1, 0, 2};
  double closest = 1e300;
  const int limitedKind = onlyKind;
  for (int pass = 0; pass < 2 && !(best <= gate); pass++) {
    if (pass == 1) {
      // A first template that never came close is liquid. A near miss still
      // gets the other templates; that is the hexagonal shell the fingerprint
      // did not name.
      if (!recover || limitedKind < 0 || closest > 40.0) {
        break;
      }
      onlyKind = -1;
      anchorEnv = -1;
      closest = 1e300;
    }
  for (int ko = 0; ko < 3 && !(best <= gate); ko++) {
    const int kindPeek = kinds[ko];
    if (onlyKind >= 0 && kindPeek != onlyKind) {
      continue;
    }
    if (onlyKind < 0 && ko > 0 && closest > 40.0) {
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
      if (anchorEnv >= 0 && pidx[anchorTem] != anchorEnv) {
        continue;
      }
      tried++;
      // The first few correspondences name a crystal. Later ones do not
      // change the held-out labels, and a shell that is still far after the
      // first match is liquid.
      if (tried > 12 || (tried > 1 && closest > 40.0)) {
        break;
      }
      // A deferred shell whose first correspondences stay far from the template
      // does not become crystal on a later permutation.
      if (!recover && tried > 4 && closest > 30.0) {
        break;
      }
      // One crystal neighbour is not enough to justify the later permutations.
      if (shallow && tried > 6 && closest > 28.0) {
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
      Shell17 fixed17;
      Shell17 mobile17;
      for (int c = 0; c < 3; c++) {
        fixed17.p[0][c] = 0.0;
        mobile17.p[0][c] = 0.0;
      }
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
        int nC = 0;
        if (grouped) {
          cand[0] = 3 * envArm;
          cand[1] = 3 * envArm + 1;
          cand[2] = 3 * envArm + 2;
          nC = 3;
        } else {
          double candD[3];
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
        }
        if (nC < 3 || canon.nArm[kind][a] < 3) {
          paired = false;
          break;
        }
        // Azimuth around the bond. The unit frame is only a positive scale on
        // both axes, and the sort reads signs, so the square roots drop out.
        const double *axis = envF[envArm];
        const double *refv = envF[pidx[(a + 1) % 4]];
        const double al2 = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
        const double proj = al2 > 1e-18 ? (refv[0] * axis[0] + refv[1] * axis[1] + refv[2] * axis[2]) / al2 : 0.0;
        const double rx = refv[0] - proj * axis[0];
        const double ry = refv[1] - proj * axis[1];
        const double rz = refv[2] - proj * axis[2];
        const double wx = axis[1] * rz - axis[2] * ry;
        const double wy = axis[2] * rx - axis[0] * rz;
        const double wz = axis[0] * ry - axis[1] * rx;
        double ax[3];
        double ay[3];
        for (int i = 0; i < 3; i++) {
          const double ox = envS[cand[i]][0] - axis[0];
          const double oy = envS[cand[i]][1] - axis[1];
          const double oz = envS[cand[i]][2] - axis[2];
          ax[i] = ox * rx + oy * ry + oz * rz;
          ay[i] = ox * wx + oy * wy + oz * wz;
        }
        for (int i = 0; i < 3; i++) {
          for (int j = i + 1; j < 3; j++) {
            const bool negi = ay[i] < 0.0 || (ay[i] == 0.0 && ax[i] < 0.0);
            const bool negj = ay[j] < 0.0 || (ay[j] == 0.0 && ax[j] < 0.0);
            const bool less = negi != negj ? negi : (ax[i] * ay[j] - ay[i] * ax[j]) > 0.0;
            if (!less) {
              std::swap(ax[i], ax[j]);
              std::swap(ay[i], ay[j]);
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
        // A crystal's first arm already sits on the template. Liquid does not,
        // and the other three arms are not paired.
        if (!recover && a == 0 && bestShiftErr > 30.0) {
          paired = false;
        }
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
  }
  if (bestK < 0 || !(best <= gate)) {
    return -1;
  }
  return bestK;
}

bool gatherTwelve(int self, const Bond *mine, const Bond *bonds, const int *nbon, bool wide, double *fx,
                  double *fy, double *fz) {
  int nf = 0;
  for (int k = 0; k < 4; k++) {
    const int f = mine[k].j;
    const int na = nbon[static_cast<size_t>(f)];
    const Bond *arm = bonds + static_cast<size_t>(f) * static_cast<size_t>(kMaxBond);
    int taken = 0;
    for (int t = 0; t < na && taken < 3; t++) {
      if (wide) {
        const int j = arm[t].j;
        if (j == self || j == mine[0].j || j == mine[1].j || j == mine[2].j || j == mine[3].j) {
          continue;
        }
      }
      const double x = mine[k].dx + arm[t].dx;
      const double y = mine[k].dy + arm[t].dy;
      const double z = mine[k].dz + arm[t].dz;
      if (!wide) {
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
      }
      if (nf >= 12) {
        return false;
      }
      fx[nf] = x;
      fy[nf] = y;
      fz[nf] = z;
      nf++;
      taken++;
    }
    if (taken < 3) {
      return false;
    }
  }
  return nf == 12;
}

} // namespace

int bondsBeyond(int atom, const Packed *packed, const int *cellStart, int nx, int ny, int nz, int reach,
               const int *bin, const double lat[3][3], const double *px, const double *py, const double *pz,
               Bond *out);

void recoverHex(int n, const Bond *bonds, const int *nbon, const Packed *packed, const int *cellStart,
                int nx, int ny, int nz, int reach, const int *bin, const double lat[3][3], const double *px,
                const double *py, const double *pz, int *kind, const CanonShell &canon, double gate, char *crystal,
                char *tried);

void fillHalf(int ncell, int nx, int ny, int nz, const int *cellStart, const Packed *packed,
              const double lat[3][3], Bond *bonds, int *nbon) {
  const double boxX = lat[0][0];
  const double boxY = lat[1][1];
  const double boxZ = lat[2][2];
  const double cut2 = kBond * kBond;
  const double min2 = 2.15 * 2.15;
  for (int cell = 0; cell < ncell; cell++) {
    const int ix0 = cell % nx;
    const int iy0 = (cell / nx) % ny;
    const int iz0 = cell / (nx * ny);
    const int here0 = cellStart[cell];
    const int here1 = cellStart[cell + 1];
    for (int mine = here0; mine < here1; mine++) {
      const int i = packed[mine].id;
      const double ci0 = packed[mine].x;
      const double ci1 = packed[mine].y;
      const double ci2 = packed[mine].z;
      for (int dz = -1; dz <= 1; dz++) {
        int iz = iz0 + dz;
        int sz = 0;
        if (iz < 0) {
          iz += nz;
          sz = -1;
        } else if (iz >= nz) {
          iz -= nz;
          sz = 1;
        }
        for (int dy = -1; dy <= 1; dy++) {
          int iy = iy0 + dy;
          int sy = 0;
          if (iy < 0) {
            iy += ny;
            sy = -1;
          } else if (iy >= ny) {
            iy -= ny;
            sy = 1;
          }
          for (int dx = -1; dx <= 1; dx++) {
            int ix = ix0 + dx;
            int sx = 0;
            if (ix < 0) {
              ix += nx;
              sx = -1;
            } else if (ix >= nx) {
              ix -= nx;
              sx = 1;
            }
            const int nbcell = (iz * ny + iy) * nx + ix;
            const int begin = cellStart[nbcell];
            const int end = cellStart[nbcell + 1];
            const double ox = sx * boxX;
            const double oy = sy * boxY;
            const double oz = sz * boxZ;
            for (int m = begin; m < end; m++) {
              const int j = packed[m].id;
              if (j < i) {
                continue;
              }
              if (j == i && sx == 0 && sy == 0 && sz == 0) {
                continue;
              }
              const double ddx = packed[m].x + ox - ci0;
              const double ddy = packed[m].y + oy - ci1;
              const double ddz = packed[m].z + oz - ci2;
              const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
              if (d2 > cut2 || d2 < min2) {
                continue;
              }
              pushClosest(bonds + static_cast<size_t>(i) * kMaxBond, nbon[i], j, ddx, ddy, ddz, d2);
              if (j != i) {
                pushClosest(bonds + static_cast<size_t>(j) * kMaxBond, nbon[j], i, -ddx, -ddy, -ddz, d2);
              }
            }
          }
        }
      }
    }
  }
}

__attribute__((always_inline)) inline void keepId(Bond *row, double *rd, int &nb, int j, double d2) {
  if (nb < 4) {
    row[nb].j = j;
    rd[nb] = d2;
    nb++;
    return;
  }
  int worst = rd[1] > rd[0] ? 1 : 0;
  if (rd[2] > rd[worst]) {
    worst = 2;
  }
  if (rd[3] > rd[worst]) {
    worst = 3;
  }
  if (d2 < rd[worst]) {
    row[worst].j = j;
    rd[worst] = d2;
  }
}

// Neighbour ids from one half of the cutoff CSR. Interior cells use fixed
// offsets, so no per-cell span is stored. Border cells keep fourteen ranges
// on the stack. Bond vectors and the twelve-neighbour signature are filled
// from these ids and the SoA coordinates, not from a second list.
__attribute__((always_inline)) inline void keepNear(int *js, double *rd, int &nb, int j, double d2) {
  if (nb < 4) {
    js[nb] = j;
    rd[nb] = d2;
    nb++;
    return;
  }
  int worst = rd[1] > rd[0] ? 1 : 0;
  if (rd[2] > rd[worst]) {
    worst = 2;
  }
  if (rd[3] > rd[worst]) {
    worst = 3;
  }
  if (d2 < rd[worst]) {
    js[worst] = j;
    rd[worst] = d2;
  }
}

__attribute__((always_inline)) inline void keepVec(int *js, double *rd, double *vx, double *vy, double *vz, int &nb,
                                                   int j, double d2, double dx, double dy, double dz) {
  if (nb < 4) {
    js[nb] = j;
    rd[nb] = d2;
    vx[nb] = dx;
    vy[nb] = dy;
    vz[nb] = dz;
    nb++;
    return;
  }
  int worst = rd[1] > rd[0] ? 1 : 0;
  if (rd[2] > rd[worst]) {
    worst = 2;
  }
  if (rd[3] > rd[worst]) {
    worst = 3;
  }
  if (d2 < rd[worst]) {
    js[worst] = j;
    rd[worst] = d2;
    vx[worst] = dx;
    vy[worst] = dy;
    vz[worst] = dz;
  }
}

void fillHalfFast(int ncell, int nx, int ny, int nz, double lx, double ly, double lz, const int *cellStart,
                  const double *cx, const double *cy, const double *cz, const int *cid, int *js, int *nbon,
                  double *rd, double *vx, double *vy, double *vz) {
  const double cut2 = kBond * kBond;
  const double min2 = 2.15 * 2.15;
  const int nxy = nx * ny;
  const int off[13] = {1,
                       nx - 1, nx, nx + 1,
                       nxy - nx - 1, nxy - nx, nxy - nx + 1,
                       nxy - 1, nxy, nxy + 1,
                       nxy + nx - 1, nxy + nx, nxy + nx + 1};
  for (int iz0 = 1; iz0 < nz - 1; iz0++) {
    for (int iy0 = 1; iy0 < ny - 1; iy0++) {
      for (int ix0 = 1; ix0 < nx - 1; ix0++) {
        const int cell = (iz0 * ny + iy0) * nx + ix0;
        const int here0 = cellStart[cell];
        const int here1 = cellStart[cell + 1];
        if (here0 == here1) {
          continue;
        }
        for (int mine = here0; mine < here1; mine++) {
          const int i = cid[mine];
          const double xi = cx[mine];
          const double yi = cy[mine];
          const double zi = cz[mine];
          for (int m = mine + 1; m < here1; m++) {
            const double dx = cx[m] - xi;
            const double dy = cy[m] - yi;
            const double dz = cz[m] - zi;
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 > cut2 || d2 < min2) {
              continue;
            }
            const int j = cid[m];
            keepVec(js + static_cast<size_t>(i) * 4, rd + static_cast<size_t>(i) * 4, vx + static_cast<size_t>(i) * 4,
                    vy + static_cast<size_t>(i) * 4, vz + static_cast<size_t>(i) * 4, nbon[i], j, d2, dx, dy, dz);
            keepVec(js + static_cast<size_t>(j) * 4, rd + static_cast<size_t>(j) * 4, vx + static_cast<size_t>(j) * 4,
                    vy + static_cast<size_t>(j) * 4, vz + static_cast<size_t>(j) * 4, nbon[j], i, d2, -dx, -dy, -dz);
          }
          for (int t = 0; t < 13; t++) {
            const int nb = cell + off[t];
            const int b = cellStart[nb];
            const int e = cellStart[nb + 1];
            for (int m = b; m < e; m++) {
              const double dx = cx[m] - xi;
              const double dy = cy[m] - yi;
              const double dz = cz[m] - zi;
              const double d2 = dx * dx + dy * dy + dz * dz;
              if (d2 > cut2 || d2 < min2) {
                continue;
              }
              const int j = cid[m];
              keepVec(js + static_cast<size_t>(i) * 4, rd + static_cast<size_t>(i) * 4, vx + static_cast<size_t>(i) * 4,
                      vy + static_cast<size_t>(i) * 4, vz + static_cast<size_t>(i) * 4, nbon[i], j, d2, dx, dy, dz);
              keepVec(js + static_cast<size_t>(j) * 4, rd + static_cast<size_t>(j) * 4, vx + static_cast<size_t>(j) * 4,
                      vy + static_cast<size_t>(j) * 4, vz + static_cast<size_t>(j) * 4, nbon[j], i, d2, -dx, -dy, -dz);
            }
          }
        }
      }
    }
  }
  struct Nb {
    int b;
    int e;
    double ox;
    double oy;
    double oz;
    int same;
  };
  for (int cell = 0; cell < ncell; cell++) {
    const int ix0 = cell % nx;
    const int iy0 = (cell / nx) % ny;
    const int iz0 = cell / nxy;
    if (ix0 > 0 && iy0 > 0 && iz0 > 0 && ix0 + 1 < nx && iy0 + 1 < ny && iz0 + 1 < nz) {
      continue;
    }
    const int here0 = cellStart[cell];
    const int here1 = cellStart[cell + 1];
    if (here0 == here1) {
      continue;
    }
    Nb sp[14];
    int ns = 0;
    for (int dz = -1; dz <= 1; dz++) {
      int iz = iz0 + dz;
      double oz = 0.0;
      if (iz < 0) {
        iz += nz;
        oz = -lz;
      } else if (iz >= nz) {
        iz -= nz;
        oz = lz;
      }
      for (int dy = -1; dy <= 1; dy++) {
        int iy = iy0 + dy;
        double oy = 0.0;
        if (iy < 0) {
          iy += ny;
          oy = -ly;
        } else if (iy >= ny) {
          iy -= ny;
          oy = ly;
        }
        for (int dx = -1; dx <= 1; dx++) {
          if (!(dz > 0 || (dz == 0 && dy > 0) || (dz == 0 && dy == 0 && dx >= 0))) {
            continue;
          }
          int ix = ix0 + dx;
          double ox = 0.0;
          if (ix < 0) {
            ix += nx;
            ox = -lx;
          } else if (ix >= nx) {
            ix -= nx;
            ox = lx;
          }
          const int nb = (iz * ny + iy) * nx + ix;
          const int b = cellStart[nb];
          const int e = cellStart[nb + 1];
          if (b == e) {
            continue;
          }
          sp[ns++] = Nb{b, e, ox, oy, oz, nb == cell && ox == 0.0 && oy == 0.0 && oz == 0.0};
        }
      }
    }
    for (int mine = here0; mine < here1; mine++) {
      const int i = cid[mine];
      const double xi = cx[mine];
      const double yi = cy[mine];
      const double zi = cz[mine];
      for (int t = 0; t < ns; t++) {
        const double ox = sp[t].ox;
        const double oy = sp[t].oy;
        const double oz = sp[t].oz;
        const int same = sp[t].same;
        for (int m = sp[t].b; m < sp[t].e; m++) {
          if (same && m <= mine) {
            continue;
          }
          const int j = cid[m];
          if (j == i) {
            continue;
          }
          const double dx = cx[m] + ox - xi;
          const double dy = cy[m] + oy - yi;
          const double dz = cz[m] + oz - zi;
          const double d2 = dx * dx + dy * dy + dz * dz;
          if (d2 > cut2 || d2 < min2) {
            continue;
          }
          keepVec(js + static_cast<size_t>(i) * 4, rd + static_cast<size_t>(i) * 4, vx + static_cast<size_t>(i) * 4,
                  vy + static_cast<size_t>(i) * 4, vz + static_cast<size_t>(i) * 4, nbon[i], j, d2, dx, dy, dz);
          keepVec(js + static_cast<size_t>(j) * 4, rd + static_cast<size_t>(j) * 4, vx + static_cast<size_t>(j) * 4,
                  vy + static_cast<size_t>(j) * 4, vz + static_cast<size_t>(j) * 4, nbon[j], i, d2, -dx, -dy, -dz);
        }
      }
    }
  }
}

// The pair walk already stored the minimum-image vector with each kept neighbour.
// Sorting those four rows replaces a second pass over every atom.
void commitBonds(int n, Bond *bonds, const int *nbon, int *js, double *rd, double *vx, double *vy, double *vz) {
  for (int i = 0; i < n; i++) {
    const int nb = nbon[i];
    int *J = js + static_cast<size_t>(i) * 4;
    double *R = rd + static_cast<size_t>(i) * 4;
    double *X = vx + static_cast<size_t>(i) * 4;
    double *Y = vy + static_cast<size_t>(i) * 4;
    double *Z = vz + static_cast<size_t>(i) * 4;
    for (int a = 1; a < nb; a++) {
      const int j = J[a];
      const double d2 = R[a];
      const double dx = X[a];
      const double dy = Y[a];
      const double dz = Z[a];
      int b = a;
      while (b > 0 && R[b - 1] > d2) {
        J[b] = J[b - 1];
        R[b] = R[b - 1];
        X[b] = X[b - 1];
        Y[b] = Y[b - 1];
        Z[b] = Z[b - 1];
        --b;
      }
      J[b] = j;
      R[b] = d2;
      X[b] = dx;
      Y[b] = dy;
      Z[b] = dz;
    }
    Bond *row = bonds + static_cast<size_t>(i) * 4;
    for (int a = 0; a < nb; a++) {
      row[a] = Bond{J[a], X[a], Y[a], Z[a]};
    }
  }
}

// Four nearest bonds inside the cutoff. Interior cells share one 27-offset
// stencil, the same orthorhombic CSR walk, and do not rebuild it per cell.
void fillFastBonds(int ncell, int nx, int ny, int nz, double lx, double ly, double lz, const int *cellStart,
                   const double *cx, const double *cy, const double *cz, const int *cid, Bond *bonds, int *nbon) {
  const double cut2 = kBond * kBond;
  const double min2 = 2.15 * 2.15;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (ncell >= 64)
#endif
  for (int cell = 0; cell < ncell; cell++) {
    const int here0 = cellStart[cell];
    const int here1 = cellStart[cell + 1];
    if (here0 == here1) {
      continue;
    }
    const int ix0 = cell % nx;
    const int iy0 = (cell / nx) % ny;
    const int iz0 = cell / (nx * ny);
    const bool interior = ix0 > 0 && iy0 > 0 && iz0 > 0 && ix0 + 1 < nx && iy0 + 1 < ny && iz0 + 1 < nz;
    struct Span {
      int b, e;
      double ox, oy, oz;
    };
    Span sp[27];
    int ns = 0;
    for (int dz = -1; dz <= 1; dz++) {
      int iz = iz0 + dz;
      double oz = 0.0;
      if (!interior) {
        if (iz < 0) {
          iz += nz;
          oz = -1.0;
        } else if (iz >= nz) {
          iz -= nz;
          oz = 1.0;
        }
      }
      for (int dy = -1; dy <= 1; dy++) {
        int iy = iy0 + dy;
        double oy = 0.0;
        if (!interior) {
          if (iy < 0) {
            iy += ny;
            oy = -1.0;
          } else if (iy >= ny) {
            iy -= ny;
            oy = 1.0;
          }
        }
        for (int dx = -1; dx <= 1; dx++) {
          int ix = ix0 + dx;
          double ox = 0.0;
          if (!interior) {
            if (ix < 0) {
              ix += nx;
              ox = -1.0;
            } else if (ix >= nx) {
              ix -= nx;
              ox = 1.0;
            }
          }
          const int nb = (iz * ny + iy) * nx + ix;
          const int b = cellStart[nb];
          const int e = cellStart[nb + 1];
          if (b != e) {
            sp[ns++] = Span{b, e, ox, oy, oz};
          }
        }
      }
    }
    for (int mine = here0; mine < here1; mine++) {
      const int i = cid[mine];
      const double xi = cx[mine];
      const double yi = cy[mine];
      const double zi = cz[mine];
      int nb = 0;
      int bj[4];
      double bx[4], by[4], bz[4], bd[4];
      for (int t = 0; t < ns; t++) {
        const int b = sp[t].b;
        const int e = sp[t].e;
        const double ox = sp[t].ox;
        const double oy = sp[t].oy;
        const double oz = sp[t].oz;
        for (int m = b; m < e; m++) {
          if (cid[m] == i && ox == 0.0 && oy == 0.0 && oz == 0.0) {
            continue;
          }
          const double dx = cx[m] + ox * lx - xi;
          const double dy = cy[m] + oy * ly - yi;
          const double dz = cz[m] + oz * lz - zi;
          const double d2 = dx * dx + dy * dy + dz * dz;
          if (d2 > cut2 || d2 < min2) {
            continue;
          }
          if (nb < 4) {
            bj[nb] = cid[m];
            bx[nb] = dx;
            by[nb] = dy;
            bz[nb] = dz;
            bd[nb] = d2;
            nb++;
            continue;
          }
          int worst = 0;
          if (bd[1] > bd[0]) {
            worst = 1;
          }
          if (bd[2] > bd[worst]) {
            worst = 2;
          }
          if (bd[3] > bd[worst]) {
            worst = 3;
          }
          if (d2 < bd[worst]) {
            bj[worst] = cid[m];
            bx[worst] = dx;
            by[worst] = dy;
            bz[worst] = dz;
            bd[worst] = d2;
          }
        }
      }
      for (int a = 1; a < nb; a++) {
        const int j = bj[a];
        const double x = bx[a];
        const double y = by[a];
        const double z = bz[a];
        const double d2 = bd[a];
        int at = a;
        while (at > 0 && d2 < bd[at - 1]) {
          bj[at] = bj[at - 1];
          bx[at] = bx[at - 1];
          by[at] = by[at - 1];
          bz[at] = bz[at - 1];
          bd[at] = bd[at - 1];
          at--;
        }
        bj[at] = j;
        bx[at] = x;
        by[at] = y;
        bz[at] = z;
        bd[at] = d2;
      }
      nbon[i] = nb;
      Bond *dest = bonds + static_cast<size_t>(i) * 4;
      for (int a = 0; a < nb; a++) {
        dest[a] = Bond{bj[a], bx[a], by[a], bz[a]};
      }
    }
  }
}

void topUpBonds(int self, double ci0, double ci1, double ci2, int ix0, int iy0, int iz0, int nx, int ny, int nz,
                const double lat[3][3], const int *cellStart, const Packed *packed, Bond *row, double *rd,
                int &nb) {
  const double near2 = kNear * kNear;
  const double cut2 = kBond * kBond;
  for (int dz = -2; dz <= 2; dz++) {
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
    for (int dy = -2; dy <= 2; dy++) {
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
      for (int dx = -2; dx <= 2; dx++) {
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
        const int nbcell = (iz * ny + iy) * nx + ix;
        const int begin = cellStart[nbcell];
        const int end = cellStart[nbcell + 1];
        const double ox = sx * lat[0][0] + sy * lat[1][0] + sz * lat[2][0];
        const double oy = sx * lat[0][1] + sy * lat[1][1] + sz * lat[2][1];
        const double oz = sx * lat[0][2] + sy * lat[1][2] + sz * lat[2][2];
        for (int m = begin; m < end; m++) {
          const Packed &nbr = packed[m];
          if (nbr.id == self && sx == 0 && sy == 0 && sz == 0) {
            continue;
          }
          const double ddx = nbr.x + ox - ci0;
          const double ddy = nbr.y + oy - ci1;
          const double ddz = nbr.z + oz - ci2;
          const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
          if (d2 <= near2 || d2 > cut2 || d2 < 2.15 * 2.15) {
            continue;
          }
          if (nb < kMaxBond) {
            row[nb] = Bond{nbr.id, ddx, ddy, ddz};
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
              row[worst] = Bond{nbr.id, ddx, ddy, ddz};
              rd[worst] = d2;
            }
          }
        }
      }
    }
  }
}

void classifyPrepared(const double *x, const double *y, const double *z, int n,
                      const double lat[3][3], const CanonShell &canon, double gate, int *kind,
                      double thermalQ, double bondRatio) {
  if (n <= 0) {
    return;
  }
  // The shell walk is OpenMP. A few thousand atoms do not feed eight
  // threads, and this call does not start a second pool.
  seams::prepareCompute();
  [[maybe_unused]] const seams::TeamWidth team(n);
  const bool timed = symTime();
  const auto tAll = std::chrono::steady_clock::now();
  // Nine of the twelve second-neighbour signatures is the ambiguous band
  // that still enters the fit. Exact FCC and HCP fingerprints do not.
  constexpr int kAmbiguous = 9;
  constexpr double kQFit = 0.90;
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
    // One cell is at most one 3.5 A bond wide, the same 27-cell cutoff
    // stencil as the shared orthorhombic CSR. The four nearest inside that
    // cutoff are the shell; a second search out to the bond length is not
    // required.
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
  Workspace &mem = workspace();
  mem.fit(static_cast<size_t>(n), static_cast<size_t>(ncell));
  std::memset(mem.cellStart, 0, (static_cast<size_t>(ncell) + 1) * sizeof(int));
  std::memset(mem.nbon, 0, static_cast<size_t>(n) * sizeof(int));
  View<double> px{mem.px};
  View<double> py{mem.py};
  View<double> pz{mem.pz};
  View<int> bin{mem.bin};
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 4096)
#endif
  for (int i = 0; i < n; i++) {
    const double c0 = x[i];
    const double c1 = y[i];
    const double c2 = z[i];
    double f[3] = {0, 0, 0};
    double wrapped[3] = {0, 0, 0};
    if (ortho) {
      f[0] = c0 / lat[0][0];
      f[1] = c1 / lat[1][1];
      f[2] = c2 / lat[2][2];
      for (int a = 0; a < 3; a++) {
        f[a] -= std::floor(f[a]);
      }
      wrapped[0] = f[0] * lat[0][0];
      wrapped[1] = f[1] * lat[1][1];
      wrapped[2] = f[2] * lat[2][2];
    } else {
      f[0] = invT[0][0] * c0 + invT[0][1] * c1 + invT[0][2] * c2;
      f[1] = invT[1][0] * c0 + invT[1][1] * c1 + invT[1][2] * c2;
      f[2] = invT[2][0] * c0 + invT[2][1] * c1 + invT[2][2] * c2;
      for (int a = 0; a < 3; a++) {
        f[a] -= std::floor(f[a]);
      }
      for (int c = 0; c < 3; c++) {
        for (int a = 0; a < 3; a++) {
          wrapped[c] += f[a] * lat[a][c];
        }
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
  View<int> cellStart{mem.cellStart};
  for (int i = 0; i < n; i++) {
    ++cellStart[static_cast<size_t>(bin[static_cast<size_t>(i)]) + 1];
  }
  for (int c = 0; c < ncell; c++) {
    cellStart[static_cast<size_t>(c) + 1] += cellStart[static_cast<size_t>(c)];
  }
  View<int> fill{mem.fill};
  std::memcpy(fill.p, cellStart.p, static_cast<size_t>(ncell) * sizeof(int));
  View<Packed> packed{mem.packed};
  for (int i = 0; i < n; i++) {
    const int slot = fill[static_cast<size_t>(bin[static_cast<size_t>(i)])]++;
    packed[static_cast<size_t>(slot)] = Packed{px[static_cast<size_t>(i)], py[static_cast<size_t>(i)],
                                               pz[static_cast<size_t>(i)], i};
  }
  const int span = 2 * reach + 1;
  const int nImage = span * span * span;
  View<Bond> bonds{mem.bonds};
  View<int> nbon{mem.nbon};
#ifdef SEAMS_HAS_OPENMP
  const int threads = omp_get_max_threads();
#else
  const int threads = 1;
#endif
  const bool halfPairs = false;
  if (fast) {
    double *cx = mem.cx;
    double *cy = mem.cy;
    double *cz = mem.cz;
    int *cid = mem.cid;
    for (int s = 0; s < n; s++) {
      const Packed &p = packed[static_cast<size_t>(s)];
      cx[s] = p.x;
      cy[s] = p.y;
      cz[s] = p.z;
      cid[s] = p.id;
    }
    if (threads == 1) {
      fillHalfFast(ncell, nx, ny, nz, lat[0][0], lat[1][1], lat[2][2], cellStart.data(), cx, cy, cz, cid, mem.js,
                   nbon.data(), mem.bd, mem.vx, mem.vy, mem.vz);
      commitBonds(n, bonds.data(), nbon.data(), mem.js, mem.bd, mem.vx, mem.vy, mem.vz);
    } else {
      fillFastBonds(ncell, nx, ny, nz, lat[0][0], lat[1][1], lat[2][2], cellStart.data(), cx, cy, cz, cid,
                    bonds.data(), nbon.data());
    }
  }
  const int cellWalk = fast ? 0 : ncell;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (cellWalk >= 64 && threads > 1)
#endif
  for (int cell = 0; cell < cellWalk; cell++) {
    const int ix0 = cell % nx;
    const int iy0 = (cell / nx) % ny;
    const int iz0 = cell / (nx * ny);
    const bool interior = fast && ix0 > 0 && iy0 > 0 && iz0 > 0 && ix0 + 1 < nx && iy0 + 1 < ny &&
                          iz0 + 1 < nz;
    if (interior) {
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
        for (int dz = -1; dz <= 1; dz++) {
          const int iz = iz0 + dz;
          for (int dy = -1; dy <= 1; dy++) {
            const int iy = iy0 + dy;
            const int row0 = (iz * ny + iy) * nx + ix0;
            for (int dx = -1; dx <= 1; dx++) {
              const int nbcell = row0 + dx;
              const int begin = cellStart[static_cast<size_t>(nbcell)];
              const int end = cellStart[static_cast<size_t>(nbcell) + 1];
              for (int m = begin; m < end; m++) {
                const Packed &nbr = packed[static_cast<size_t>(m)];
                if (nbr.id == i && dx == 0 && dy == 0 && dz == 0) {
                  continue;
                }
                if (halfPairs && nbr.id < i) {
                  continue;
                }
                const double ddx = nbr.x - ci0;
                const double ddy = nbr.y - ci1;
                const double ddz = nbr.z - ci2;
                const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
                if (d2 > kBond * kBond || d2 < 2.15 * 2.15) {
                  continue;
                }
                if (nb < kMaxBond) {
                  row[nb] = Bond{nbr.id, ddx, ddy, ddz};
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
                    row[worst] = Bond{nbr.id, ddx, ddy, ddz};
                    rd[worst] = d2;
                  }
                }
                if (halfPairs && nbr.id != i) {
                  int &onb = nbon[static_cast<size_t>(nbr.id)];
                  pushClosest(&bonds[static_cast<size_t>(nbr.id) * static_cast<size_t>(kMaxBond)], onb, i,
                              -ddx, -ddy, -ddz, d2);
                }
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
        Bond *dest = &bonds[static_cast<size_t>(i) * static_cast<size_t>(kMaxBond)];
        if (halfPairs) {
          int kept = nbon[static_cast<size_t>(i)];
          keepClosest(dest, kept, row, nb);
          nbon[static_cast<size_t>(i)] = kept;
        } else {
          nbon[static_cast<size_t>(i)] = nb;
          for (int a = 0; a < nb; a++) {
            dest[a] = row[a];
          }
        }
      }
      continue;
    }
    Image local[729];
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
          local[slot++] = Image{(iz * ny + iy) * nx + ix, sx, sy, sz};
        }
      }
    }
    const Image *im = local;
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
        if (halfPairs && nbr.id < i) {
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
        if (halfPairs && nbr.id != i) {
          int &onb = nbon[static_cast<size_t>(nbr.id)];
          pushClosest(&bonds[static_cast<size_t>(nbr.id) * static_cast<size_t>(kMaxBond)], onb, i, -dx, -dy,
                      -dz, d2);
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
    Bond *dest = &bonds[static_cast<size_t>(i) * static_cast<size_t>(kMaxBond)];
    if (halfPairs) {
      int kept = nbon[static_cast<size_t>(i)];
      keepClosest(dest, kept, row, nb);
      nbon[static_cast<size_t>(i)] = kept;
    } else {
      nbon[static_cast<size_t>(i)] = nb;
      for (int a = 0; a < nb; a++) {
        dest[a] = row[a];
      }
    }
    }
  }
  const auto tBonds = std::chrono::steady_clock::now();
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 64 && threads > 1)
#endif
  for (int i = 0; i < n; i++) {
    const int nb = nbon[static_cast<size_t>(i)];
    const Bond *mine = &bonds[static_cast<size_t>(i) * static_cast<size_t>(kMaxBond)];
    const double ci0 = px[static_cast<size_t>(i)];
    const double ci1 = py[static_cast<size_t>(i)];
    const double ci2 = pz[static_cast<size_t>(i)];
    if (nb < 4) {
      kind[i] = -1;
      continue;
    }
    // First shell only. A bond far from the diamond length is not a centre,
    // so the second shell is not read.
    double bMin = 1e300;
    double bMax = 0.0;
    for (int k = 0; k < 4; k++) {
      const double l2 = mine[k].dx * mine[k].dx + mine[k].dy * mine[k].dy + mine[k].dz * mine[k].dz;
      if (l2 < bMin) {
        bMin = l2;
      }
      if (l2 > bMax) {
        bMax = l2;
      }
    }
    if (bMin < 2.30 * 2.30 || bMax > 3.35 * 3.35) {
      kind[i] = -1;
      continue;
    }
    // Liquid is not tetrahedral. The gate only skips an obvious liquid
    // before the second shell is read. Ideal ice has equal bonds and does
    // not pay for q; the staggered test names it. kThermalQ is under the
    // thermal-crystal body and above the liquid median.
    double qEarly = -1.0;
    if (!(bMax <= bMin * 1.00001)) {
      qEarly = tetraQ(mine);
      if (qEarly < thermalQ) {
        kind[i] = -1;
        continue;
      }
    }
    {
      double fx[12];
      double fy[12];
      double fz[12];
      const bool wide = la > 12.0 && lb > 12.0 && lc > 12.0;
      if (!gatherTwelve(i, mine, bonds.get(), nbon.data(), wide, fx, fy, fz)) {
        kind[i] = -1;
        continue;
      }
      double sLo = 1e300;
      double sHi = 0.0;
      for (int s = 0; s < 12; s++) {
        const double r2s = fx[s] * fx[s] + fy[s] * fy[s] + fz[s] * fz[s];
        if (r2s < sLo) {
          sLo = r2s;
        }
        if (r2s > sHi) {
          sHi = r2s;
        }
      }
      if (bMax <= bMin * 1.00001 && sHi <= sLo * 1.00001) {
        double first[4][3];
        for (int k = 0; k < 4; k++) {
          first[k][0] = mine[k].dx;
          first[k][1] = mine[k].dy;
          first[k][2] = mine[k].dz;
        }
        const int even = idealKind(first, fx, fy, fz);
        if (even >= 0) {
          kind[i] = even;
          continue;
        }
      }
      // A wide second shell is not an FCC or HCP fingerprint. The fit still
      // sees a shell whose radius varies by less than 1.6.
      if (sHi > sLo * (1.45 * 1.45)) {
        const double bondCut2 = bondRatio * bondRatio;
        if (qEarly >= thermalQ && bMax <= bMin * bondCut2 &&
            sHi <= sLo * (1.6 * 1.6)) {
          kind[i] = -3;
        } else {
          kind[i] = -1;
        }
        continue;
      }
      int n421 = 0;
      int n422 = 0;
      int nGood = 0;
      const int cna = diamondCna(fx, fy, fz, kAmbiguous, &n421, &n422, &nGood);
      if (cna == 0 || cna == 1) {
        kind[i] = cna;
        continue;
      }
      const double qNow = qEarly >= 0.0 ? qEarly : tetraQ(mine);
      if (qNow < thermalQ) {
        kind[i] = -1;
        continue;
      }
      if ((cna == 3 || cna == 4) && qNow >= 0.97) {
        // Six cubic signatures are the liquid false positives on the training
        // seeds. Those shells still go to the fit; eight or more are labelled here.
        if (cna == 3 && n421 < 8) {
          kind[i] = -3;
        } else {
          kind[i] = cna == 3 ? 1 : 0;
        }
        continue;
      }
      const bool tightBonds = bMax <= bMin * (bondRatio * bondRatio);
      double sMin = 1e300;
      double sMax = 0.0;
      for (int s = 0; s < 12; s++) {
        const double r2s = fx[s] * fx[s] + fy[s] * fy[s] + fz[s] * fz[s];
        if (r2s < sMin) {
          sMin = r2s;
        }
        if (r2s > sMax) {
          sMax = r2s;
        }
      }
      // -3: worth a fit, but only once a crystal neighbour is known.
      // The neighbour test, not a quenched-crystal q, keeps bulk liquid
      // out of Kabsch. Thermal ice sits below 0.95.
      if (qNow >= thermalQ && tightBonds && sMax <= sMin * (1.6 * 1.6)) {
        kind[i] = -3;
        continue;
      }
      kind[i] = -1;
      continue;
    }
  }
  const auto tLabel = std::chrono::steady_clock::now();
  // The fit is only for a shell that already touches a labelled crystal.
  // Bulk liquid has a high q and would otherwise all enter Kabsch.
  const bool wideBox = la > 12.0 && lb > 12.0 && lc > 12.0;
  for (int round = 0; round < 4; round++) {
    for (int i = 0; i < n; i++) {
      if (kind[i] != -3) {
        continue;
      }
      const Bond *mine = &bonds[static_cast<size_t>(i) * static_cast<size_t>(kMaxBond)];
      int nCrystal = 0;
      for (int k = 0; k < 4; k++) {
        if (kind[mine[k].j] >= 0) {
          nCrystal++;
        }
      }
      if (nCrystal < 1) {
        continue;
      }
      double fx[12];
      double fy[12];
      double fz[12];
      if (!gatherTwelve(i, mine, bonds.get(), nbon.data(), wideBox, fx, fy, fz)) {
        kind[i] = -1;
        continue;
      }
      double nbFit[16][3];
      const double c0 = px[static_cast<size_t>(i)];
      const double c1 = py[static_cast<size_t>(i)];
      const double c2 = pz[static_cast<size_t>(i)];
      for (int k = 0; k < 4; k++) {
        nbFit[k][0] = c0 + mine[k].dx;
        nbFit[k][1] = c1 + mine[k].dy;
        nbFit[k][2] = c2 + mine[k].dz;
      }
      for (int s = 0; s < 12; s++) {
        nbFit[4 + s][0] = c0 + fx[s];
        nbFit[4 + s][1] = c1 + fy[s];
        nbFit[4 + s][2] = c2 + fz[s];
      }
      const double center[3] = {c0, c1, c2};
      const int fit = shellFit(center, nbFit, 16, canon, gate, false, true, nCrystal < 2);
      kind[i] = fit >= 0 ? fit : -1;
    }
  }
  const auto tFit = std::chrono::steady_clock::now();
  for (int i = 0; i < n; i++) {
    if (kind[i] == -3) {
      kind[i] = -1;
    }
  }
  // One crystal neighbour among the four bonds. Isolated liquid labels go.
  // The test reads the labels from the fit, not the filtered ones.
  int *raw = mem.rawK;
  std::memcpy(raw, kind, static_cast<size_t>(n) * sizeof(int));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 64 && threads > 1)
#endif
  for (int i = 0; i < n; i++) {
    if (raw[i] < 0) {
      continue;
    }
    const int nb = nbon[static_cast<size_t>(i)];
    const Bond *mine = &bonds[static_cast<size_t>(i) * static_cast<size_t>(kMaxBond)];
    int nCrystal = 0;
    const int lim = nb < 4 ? nb : 4;
    for (int k = 0; k < lim; k++) {
      if (raw[mine[k].j] >= 0) {
        nCrystal++;
      }
    }
    if (nCrystal < 1) {
      kind[i] = -1;
    }
  }
  recoverHex(n, bonds.get(), nbon.data(), packed.get(), cellStart.data(), nx, ny, nz, reach, bin.get(),
             lat, px.get(), py.get(), pz.get(), kind, canon, gate, mem.crystal, mem.tried);
  if (timed) {
    const auto tEnd = std::chrono::steady_clock::now();
    std::fprintf(stderr, "sym_time n=%d bonds=%.3f label=%.3f fit=%.3f recover=%.3f total=%.3f\n", n,
                 msSince(tAll) > 0 ? std::chrono::duration<double, std::milli>(tBonds - tAll).count() : 0,
                 std::chrono::duration<double, std::milli>(tLabel - tBonds).count(),
                 std::chrono::duration<double, std::milli>(tFit - tLabel).count(),
                 std::chrono::duration<double, std::milli>(tEnd - tFit).count(),
                 std::chrono::duration<double, std::milli>(tEnd - tAll).count());
  }
}

int bondsBeyond(int atom, const Packed *packed, const int *cellStart, int nx, int ny, int nz, int reach,
               const int *bin, const double lat[3][3], const double *px, const double *py, const double *pz,
               Bond *out) {
  Bond row[kMaxWide];
  double rd[kMaxWide];
  int nw = 0;
  const int cell = bin[atom];
  const int ix0 = cell % nx;
  const int iy0 = (cell / nx) % ny;
  const int iz0 = cell / (nx * ny);
  const double ci0 = px[atom];
  const double ci1 = py[atom];
  const double ci2 = pz[atom];
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
      for (int dxs = -reach; dxs <= reach; dxs++) {
        int ix = ix0 + dxs;
        int sx = 0;
        while (ix < 0) {
          ix += nx;
          sx--;
        }
        while (ix >= nx) {
          ix -= nx;
          sx++;
        }
        const int nbcell = (iz * ny + iy) * nx + ix;
        const int begin = cellStart[nbcell];
        const int end = cellStart[nbcell + 1];
    for (int m = begin; m < end; m++) {
      const Packed &nbr = packed[m];
      if (nbr.id == atom && sx == 0 && sy == 0 && sz == 0) {
        continue;
      }
      const double dx = nbr.x + sx * lat[0][0] + sy * lat[1][0] + sz * lat[2][0] - ci0;
      const double dy = nbr.y + sx * lat[0][1] + sy * lat[1][1] + sz * lat[2][1] - ci1;
      const double dz = nbr.z + sx * lat[0][2] + sy * lat[1][2] + sz * lat[2][2] - ci2;
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
                int nx, int ny, int nz, int reach, const int *bin, const double lat[3][3], const double *px,
                const double *py, const double *pz, int *kind, const CanonShell &canon, double gate, char *crystal,
                char *tried) {
  int unlabeled = 0;
  for (int i = 0; i < n; i++) {
    if (kind[i] < 0) {
      unlabeled++;
    }
  }
  if (unlabeled == 0) {
    return;
  }
  std::memset(tried, 0, static_cast<size_t>(n));
  for (int pass = 0; pass < 1; pass++) {
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
    int nCrystal = 0;
    const int lim = nb < 4 ? nb : 4;
    for (int k = 0; k < lim; k++) {
      if (crystal[static_cast<size_t>(mine[k].j)]) {
        nCrystal++;
      }
    }
    // The wide search is only for a stretched fourth bond next to a crystal.
    // Liquid has fewer than three crystal neighbours and does not pay for it.
    if (nCrystal < 3) {
      continue;
    }
    Bond mineWbuf[kMaxWide];
    const int nw = bondsBeyond(i, packed, cellStart, nx, ny, nz, reach, bin, lat, px, py, pz, mineWbuf);
    const Bond *mineW = mineWbuf;
    if (nb + nw < 4) {
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
      const int naw = bondsBeyond(f, packed, cellStart, nx, ny, nz, reach, bin, lat, px, py, pz, armWbuf);
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
