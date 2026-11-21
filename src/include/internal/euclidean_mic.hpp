//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#ifndef SEAMS_EUCLIDEAN_MIC_H_
#define SEAMS_EUCLIDEAN_MIC_H_

// Euclidean minimum image for a triclinic cell. The fractional wrap is the
// nearest image only inside the Smith ball. Past that ball the closest
// lattice point is found by the Sommer-Feder-Shalvi slicer over the
// Voronoi-relevant vectors, started from the wrap itself when no cell edge is
// longer than the reduced superbasis, and otherwise from Babai rounding in a
// Selling superbasis.
//
// The relevant vectors are found once per cell. Every Voronoi-relevant vector
// is, up to sign, the shortest member of one of the seven nonzero cosets of
// L / 2L. A superbasis that is obtuse beyond rounding gives those as its
// subset sums; otherwise Fincke-Pohst enumeration finds each coset's
// shortest member. A point that none of them shortens lies in the Voronoi
// cell of the origin, so the slicer's stopping point is the minimum image.
// The reduced and relevant vectors are re-evaluated from their integer
// coordinates in the cell, so reduction rounding does not leave the lattice.
//
// The slicer runs in binary64 and carries the integer coordinates of the
// lattice translate it has applied. Its binary64 result is kept only when
// little cancelled (see closestDisplacement); otherwise the image is
// re-evaluated as wrap - n . (a, b, c) with error-free products and
// compensated sums, and the slicer runs again on that vector with extended
// precision norms. Rounding in Babai or in a long descent, which scales with
// the superbasis and not with the wrap, therefore does not reach the result.
//
// A slicer started from the wrap falls back to Babai after four steps. No
// other loop has an iteration cap. Every accepted Lagrange, Selling, or
// slicer step strictly lowers a non-negative, finite squared norm, which
// bounds the step count. A slicer step whose rounded norm does not fall is
// skipped for the next candidate, and the loop stops when no candidate
// lowers it. Non-finite values, singular bases, and lattice coefficients at
// or above 2^52 throw explicitly without poisoning the per-thread basis
// cache.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

namespace gen {
namespace mic {
namespace {

constexpr double kIntegral = 4503599627370496.0;

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

inline double n2(Vec3 v) { return v.x * v.x + v.y * v.y + v.z * v.z; }

inline long double n2Long(Vec3 v) {
  return static_cast<long double>(v.x) * v.x +
         static_cast<long double>(v.y) * v.y +
         static_cast<long double>(v.z) * v.z;
}

inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

inline Vec3 add(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }

inline Vec3 scale(double s, Vec3 v) { return {s * v.x, s * v.y, s * v.z}; }

inline double maxAbs(Vec3 v) {
  return std::max({std::abs(v.x), std::abs(v.y), std::abs(v.z)});
}

// Half away from zero, matching a truncation after adding just under 1/2.
inline double roundAway(double x) {
  if (!std::isfinite(x) || std::abs(x) >= kIntegral) {
    throw std::runtime_error(
        "minimum-image lattice coefficient exceeds binary64 integer precision");
  }
  const double t = x + std::copysign(0.49999999999999994, x);
  return static_cast<double>(static_cast<long long>(t));
}

// Veltkamp split of a into hi + lo, each with at most 26 significant bits.
struct Split {
  double hi = 0.0;
  double lo = 0.0;
};

inline Split split(double a) {
  constexpr double kSplit = 134217729.0;
  const double c = kSplit * a;
  const double hi = c - (c - a);
  return {hi, a - hi};
}

// The rounding error of product = a * b, exactly. Without a hardware FMA,
// std::fma is a library call, so the product is split instead (Dekker); a
// target without FMA cannot contract that split into one.
inline double productError(double a, Split as, double b, Split bs, double product) {
#if defined(FP_FAST_FMA) || defined(__FMA__) || defined(__aarch64__)
  (void)as;
  (void)bs;
  return std::fma(a, b, -product);
#else
  (void)a;
  (void)b;
  return ((as.hi * bs.hi - product) + as.hi * bs.lo + as.lo * bs.hi) + as.lo * bs.lo;
#endif
}

inline double productError(double a, double b, double product) {
  return productError(a, split(a), b, split(b), product);
}

// sum n[i] * x[i] for integers n[i] below 2^52, with the products and partial
// sums carried error-free (Ogita, Rump and Oishi, Dot2), so the result is
// within about one rounding of itself however large the cancelling terms.
template <std::size_t N>
inline double integerDot(const double (&n)[N], const double (&x)[N]) {
  double sum = 0.0;
  double error = 0.0;
  for (std::size_t i = 0; i < N; i++) {
    const double product = n[i] * x[i];
    const double next = sum + product;
    const double back = next - sum;
    error += productError(n[i], x[i], product) + ((sum - (next - back)) + (product - back));
    sum = next;
  }
  return sum + error;
}

struct Obtuse {
  Vec3 v[4]{};
  // Integer coordinates of v[k] and of each relevant vector in (a, b, c).
  double coef[4][3]{};
  double rel[3][8]{};
  double relCoef[3][8]{};
  double halfN2[8]{};
  double inv[3][3]{};
  // Rows of the upper-triangular cell, a = (lx, 0, 0), b = (xy, ly, 0),
  // c = (xz, yz, lz), by component.
  double cellX[3]{};
  double cellY[3]{};
  double cellZ[3]{};
  Split splitX[3]{};
  Split splitY[3]{};
  Split splitZ[3]{};
  double maxRelN2 = 0.0;
  double maxCoef = 0.0;
  // Largest squared row of inv: |Babai coefficient|^2 <= n2(y) * maxInvN2.
  double maxInvN2 = 0.0;
  // max |component| of each vector of the Babai triple.
  double babaiMax[3]{};
  int idx[3] = {0, 1, 2};
  bool fractional = false;
  bool ok = false;
};

// a + b as sum + error, exactly (Knuth TwoSum).
inline double sumError(double a, double b, double sum) {
  const double back = sum - a;
  return (a - (sum - back)) + (b - back);
}

// w - n . (a, b, c), each component within about one rounding of itself. The
// cell is upper triangular, so x has three products, y two, and z one.
inline Vec3 latticeTranslate(const Obtuse &s, Vec3 w, const double n[3]) {
  const double m0 = -n[0];
  const double m1 = -n[1];
  const double m2 = -n[2];
  const Split s0 = split(m0);
  const Split s1 = split(m1);
  const Split s2 = split(m2);
  const double px0 = m0 * s.cellX[0];
  const double px1 = m1 * s.cellX[1];
  const double px2 = m2 * s.cellX[2];
  const double py1 = m1 * s.cellY[1];
  const double py2 = m2 * s.cellY[2];
  const double pz2 = m2 * s.cellZ[2];
  const double ex = productError(m0, s0, s.cellX[0], s.splitX[0], px0) +
                    productError(m1, s1, s.cellX[1], s.splitX[1], px1) +
                    productError(m2, s2, s.cellX[2], s.splitX[2], px2);
  const double ey = productError(m1, s1, s.cellY[1], s.splitY[1], py1) +
                    productError(m2, s2, s.cellY[2], s.splitY[2], py2);
  const double ez = productError(m2, s2, s.cellZ[2], s.splitZ[2], pz2);
  const double x01 = px0 + px1;
  const double x2w = w.x + px2;
  const double sx = x01 + x2w;
  const double sy = w.y + py1;
  const double ty = sy + py2;
  const double sz = w.z + pz2;
  return {sx + (ex + sumError(px0, px1, x01) + sumError(w.x, px2, x2w) + sumError(x01, x2w, sx)),
          ty + (ey + sumError(w.y, py1, sy) + sumError(sy, py2, ty)),
          sz + (ez + sumError(w.z, pz2, sz))};
}

inline Vec3 latticeVector(const Obtuse &s, const double n[3]) {
  const double negated[3] = {-n[0], -n[1], -n[2]};
  return latticeTranslate(s, Vec3{}, negated);
}

inline Vec3 combination(const Vec3 b[3], long long z0, long long z1, long long z2) {
  const long double c0 = static_cast<long double>(z0);
  const long double c1 = static_cast<long double>(z1);
  const long double c2 = static_cast<long double>(z2);
  return {static_cast<double>(c0 * b[0].x + c1 * b[1].x + c2 * b[2].x),
          static_cast<double>(c0 * b[0].y + c1 * b[1].y + c2 * b[2].y),
          static_cast<double>(c0 * b[0].z + c1 * b[1].z + c2 * b[2].z)};
}

inline long long firstWithParity(long double low, int parity) {
  long long n = static_cast<long long>(std::ceil(low));
  if (((n % 2) + 2) % 2 != parity) {
    n++;
  }
  return n;
}

// Shortest member of each nonzero coset of L / 2L, the parity class of the
// coordinates in `basis`. Fincke-Pohst runs over the norm-sorted basis with
// each coset's shortest {-1, 0, 1}^3 member as its radius. `out` holds the
// coordinates of each member in `basis`. Returns false when the Gram-Schmidt
// norms are not positive and finite.
inline bool cosetMinima(const Vec3 basis[3], long long out[7][3]) {
  constexpr long double kWidth = kIntegral;
  int order[3] = {0, 1, 2};
  std::sort(order, order + 3,
            [&](int p, int q) { return n2(basis[p]) < n2(basis[q]); });
  const Vec3 b[3] = {basis[order[0]], basis[order[1]], basis[order[2]]};
  long double g[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      g[i][j] = static_cast<long double>(b[i].x) * b[j].x +
                static_cast<long double>(b[i].y) * b[j].y +
                static_cast<long double>(b[i].z) * b[j].z;
    }
  }
  const long double q0 = g[0][0];
  const long double m10 = g[1][0] / q0;
  const long double q1 = g[1][1] - m10 * g[1][0];
  const long double m20 = g[2][0] / q0;
  const long double m21 = (g[2][1] - m20 * g[1][0]) / q1;
  const long double q2 = g[2][2] - m20 * g[2][0] - m21 * m21 * q1;
  if (!(std::isfinite(q0) && std::isfinite(q1) && std::isfinite(q2) && q0 > 0.0L &&
        q1 > 0.0L && q2 > 0.0L && std::isfinite(m10) && std::isfinite(m20) &&
        std::isfinite(m21))) {
    return false;
  }

  long long best[8][3]{};
  long double bestN2[8];
  for (long double &n : bestN2) {
    n = std::numeric_limits<long double>::infinity();
  }
  const auto keep = [&](int coset, long long z0, long long z1, long long z2) {
    const long double n = n2Long(combination(b, z0, z1, z2));
    if (n < bestN2[coset]) {
      bestN2[coset] = n;
      best[coset][0] = z0;
      best[coset][1] = z1;
      best[coset][2] = z2;
    }
  };
  for (int z0 = -1; z0 <= 1; z0++) {
    for (int z1 = -1; z1 <= 1; z1++) {
      for (int z2 = -1; z2 <= 1; z2++) {
        const int coset = (z0 & 1) | ((z1 & 1) << 1) | ((z2 & 1) << 2);
        if (coset != 0) {
          keep(coset, z0, z1, z2);
        }
      }
    }
  }

  for (int coset = 1; coset < 8; coset++) {
    const int p0 = coset & 1;
    const int p1 = (coset >> 1) & 1;
    const int p2 = (coset >> 2) & 1;
    const long double bound = bestN2[coset] * (1.0L + 1e-12L);
    const long double w2 = std::sqrt(bound / q2);
    if (!(w2 < kWidth)) {
      return false;
    }
    for (long long z2 = firstWithParity(-w2, p2); z2 <= w2; z2 += 2) {
      const long double t2 = static_cast<long double>(z2);
      const long double r2 = bound - q2 * t2 * t2;
      if (r2 < 0.0L) {
        continue;
      }
      const long double c1 = -m21 * t2;
      const long double w1 = std::sqrt(r2 / q1);
      if (!(std::abs(c1) + w1 < kWidth)) {
        return false;
      }
      for (long long z1 = firstWithParity(c1 - w1, p1); z1 <= c1 + w1; z1 += 2) {
        const long double e1 = static_cast<long double>(z1) - c1;
        const long double r1 = r2 - q1 * e1 * e1;
        if (r1 < 0.0L) {
          continue;
        }
        const long double c0 = -(m10 * static_cast<long double>(z1) + m20 * t2);
        const long double w0 = std::sqrt(r1 / q0);
        if (!(std::abs(c0) + w0 < kWidth)) {
          return false;
        }
        for (long long z0 = firstWithParity(c0 - w0, p0); z0 <= c0 + w0; z0 += 2) {
          keep(coset, z0, z1, z2);
        }
      }
    }
  }

  for (int k = 0; k < 7; k++) {
    for (int i = 0; i < 3; i++) {
      out[k][order[i]] = best[k + 1][i];
    }
  }
  return true;
}

inline double roundingSlack(Vec3 a, Vec3 b) {
  return 16.0 * std::numeric_limits<double>::epsilon() * std::sqrt(n2(a) * n2(b));
}

// A superbasis whose dot products are all negative beyond rounding is
// strictly obtuse, and its seven subset sums are then the Voronoi-relevant
// vectors (Conway and Sloane). Otherwise the coset minima are enumerated.
// Each relevant vector is evaluated from its integer coordinates in the cell.
inline bool relevantOf(Obtuse &o) {
  const Vec3 *v = o.v;
  bool strict = true;
  for (int i = 0; i < 3 && strict; i++) {
    for (int j = i + 1; j < 4; j++) {
      if (!(dot(v[i], v[j]) < -roundingSlack(v[i], v[j]))) {
        strict = false;
        break;
      }
    }
  }
  double coef[7][3];
  if (strict) {
    const int pairs[3] = {1, 2, 3};
    for (int k = 0; k < 4; k++) {
      for (int m = 0; m < 3; m++) {
        coef[k][m] = o.coef[k][m];
      }
    }
    for (int k = 0; k < 3; k++) {
      for (int m = 0; m < 3; m++) {
        coef[4 + k][m] = o.coef[0][m] + o.coef[pairs[k]][m];
      }
    }
  } else {
    const Vec3 basis[3] = {v[o.idx[0]], v[o.idx[1]], v[o.idx[2]]};
    long long z[7][3];
    if (!cosetMinima(basis, z)) {
      return false;
    }
    for (int k = 0; k < 7; k++) {
      for (int m = 0; m < 3; m++) {
        coef[k][m] = static_cast<double>(z[k][0]) * o.coef[o.idx[0]][m] +
                     static_cast<double>(z[k][1]) * o.coef[o.idx[1]][m] +
                     static_cast<double>(z[k][2]) * o.coef[o.idx[2]][m];
      }
    }
  }
  o.maxRelN2 = 0.0;
  for (int k = 0; k < 8; k++) {
    Vec3 r{};
    if (k < 7) {
      for (int m = 0; m < 3; m++) {
        if (!(std::abs(coef[k][m]) < kIntegral)) {
          return false;
        }
      }
      r = latticeVector(o, coef[k]);
    }
    o.rel[0][k] = r.x;
    o.rel[1][k] = r.y;
    o.rel[2][k] = r.z;
    for (int m = 0; m < 3; m++) {
      o.relCoef[m][k] = (k < 7) ? coef[k][m] : 0.0;
    }
    o.halfN2[k] = (k < 7) ? 0.5 * n2(r) : std::numeric_limits<double>::infinity();
    if (k < 7) {
      if (!std::isfinite(o.halfN2[k])) {
        return false;
      }
      o.maxRelN2 = std::max(o.maxRelN2, n2(r));
    }
  }
  return true;
}

inline bool alreadyObtuse(const Vec3 v[4]) {
  for (int i = 0; i < 3; i++) {
    for (int j = i + 1; j < 4; j++) {
      if (dot(v[i], v[j]) > 0.0) {
        return false;
      }
    }
  }
  return true;
}

inline void lagrangeReduce(Vec3 v[3]) {
  for (;;) {
    for (int i = 1; i < 3; i++) {
      int j = i;
      while (j > 0 && n2(v[j]) < n2(v[j - 1])) {
        const Vec3 tmp = v[j];
        v[j] = v[j - 1];
        v[j - 1] = tmp;
        j--;
      }
    }
    bool changed = false;
    for (int i = 1; i < 3; i++) {
      for (int j = i - 1; j >= 0; j--) {
        const double lj = n2(v[j]);
        if (lj < 1e-30) {
          continue;
        }
        const double q = roundAway(dot(v[i], v[j]) / lj);
        if (q != 0.0) {
          const Vec3 reduced = add(v[i], scale(-q, v[j]));
          const double before = n2(v[i]);
          const double after = n2(reduced);
          if (!std::isfinite(after)) {
            throw std::runtime_error(
                "minimum-image Lagrange reduction reached binary64 limits");
          }
          if (after < before) {
            v[i] = reduced;
            changed = true;
          }
        }
      }
    }
    if (!changed) {
      break;
    }
  }
}

inline void sellingSteps(Vec3 v[4]) {
  for (;;) {
    int si = -1;
    int sj = -1;
    double bestDot = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = i + 1; j < 4; j++) {
        const double d = dot(v[i], v[j]);
        if (d > roundingSlack(v[i], v[j]) && d > bestDot) {
          bestDot = d;
          si = i;
          sj = j;
        }
      }
    }
    if (si < 0) {
      return;
    }
    long double before = 0.0L;
    for (int i = 0; i < 4; i++) {
      before += n2Long(v[i]);
    }
    const Vec3 taken = v[si];
    Vec3 saved[4];
    for (int k = 0; k < 4; k++) {
      saved[k] = v[k];
    }
    for (int k = 0; k < 4; k++) {
      if (k != si && k != sj) {
        v[k] = add(v[k], taken);
      }
    }
    v[si] = scale(-1.0, taken);
    long double after = 0.0L;
    for (int i = 0; i < 4; i++) {
      after += n2Long(v[i]);
    }
    if (!std::isfinite(after)) {
      throw std::runtime_error(
          "minimum-image Selling reduction reached binary64 limits");
    }
    if (!(after < before)) {
      for (int k = 0; k < 4; k++) {
        v[k] = saved[k];
      }
      return;
    }
  }
}

// Inverse of the matrix with columns cols[0..2]. inv is stored so that
// coeff k is sum_i inv[i][k] * y[i]. Returns |det|, or 0 when singular.
inline double invertColumns(const Vec3 cols[3], double inv[3][3]) {
  const Vec3 &a = cols[0];
  const Vec3 &b = cols[1];
  const Vec3 &c = cols[2];
  const double det = a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) +
                     a.z * (b.x * c.y - b.y * c.x);
  if (!std::isfinite(det) || std::abs(det) < 1e-18) {
    return 0.0;
  }
  const double invdet = 1.0 / det;
  const double cof[3][3] = {
      {b.y * c.z - b.z * c.y, a.z * c.y - a.y * c.z, a.y * b.z - a.z * b.y},
      {b.z * c.x - b.x * c.z, a.x * c.z - a.z * c.x, a.z * b.x - a.x * b.z},
      {b.x * c.y - b.y * c.x, a.y * c.x - a.x * c.y, a.x * b.y - a.y * b.x}};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv[i][j] = cof[i][j] * invdet;
      if (!std::isfinite(inv[i][j])) {
        return 0.0;
      }
    }
  }
  return std::abs(det);
}

inline Obtuse finish(Vec3 a, Vec3 b, Vec3 c, const Vec3 v[4], const double coef[4][3]) {
  Obtuse o;
  const Vec3 cell[3] = {a, b, c};
  for (int k = 0; k < 3; k++) {
    o.cellX[k] = cell[k].x;
    o.cellY[k] = cell[k].y;
    o.cellZ[k] = cell[k].z;
    o.splitX[k] = split(cell[k].x);
    o.splitY[k] = split(cell[k].y);
    o.splitZ[k] = split(cell[k].z);
  }
  double longestEdge = 0.0;
  double longestReduced = 0.0;
  for (int i = 0; i < 4; i++) {
    o.v[i] = v[i];
    longestReduced = std::max(longestReduced, n2(v[i]));
    for (int m = 0; m < 3; m++) {
      o.coef[i][m] = coef[i][m];
      o.maxCoef = std::max(o.maxCoef, std::abs(coef[i][m]));
    }
    if (i < 3) {
      longestEdge = std::max(longestEdge, n2(cell[i]));
    }
  }
  o.fractional = longestEdge <= longestReduced;
  // Every triple of a superbasis spans the same cell volume, so Babai takes
  // the most orthogonal one, |det| / (|b0| |b1| |b2|).
  double bestAbs = 0.0;
  double bestScore = 0.0;
  bool have = false;
  for (int drop = 0; drop < 4; drop++) {
    int idx[3];
    int n = 0;
    for (int t = 0; t < 4; t++) {
      if (t != drop) {
        idx[n++] = t;
      }
    }
    const Vec3 cols[3] = {v[idx[0]], v[idx[1]], v[idx[2]]};
    double inv[3][3];
    const double absDet = invertColumns(cols, inv);
    const double score =
        absDet / (std::sqrt(n2(cols[0])) * std::sqrt(n2(cols[1])) * std::sqrt(n2(cols[2])));
    if (absDet > 0.0 && score > bestScore) {
      bestScore = score;
      bestAbs = absDet;
      have = true;
      for (int k = 0; k < 3; k++) {
        o.idx[k] = idx[k];
        for (int j = 0; j < 3; j++) {
          o.inv[k][j] = inv[k][j];
        }
      }
    }
  }
  if (!have || !(bestAbs > 1e-18)) {
    return o;
  }
  for (int k = 0; k < 3; k++) {
    o.babaiMax[k] = maxAbs(v[o.idx[k]]);
    o.maxInvN2 = std::max(o.maxInvN2, o.inv[0][k] * o.inv[0][k] + o.inv[1][k] * o.inv[1][k] +
                                          o.inv[2][k] * o.inv[2][k]);
  }
  o.ok = relevantOf(o);
  return o;
}

// Lagrange and Selling form v - q w in binary64, and with q near 1e7 that
// rounding moves a reduced vector off the lattice by far more than one ulp of
// its own length. Each vector is re-evaluated from its integer coordinates in
// the upper-triangular cell (a, b, c), recovered by back-substitution and
// returned in coef. The result must still span the cell volume.
inline void snapToLattice(Vec3 a, Vec3 b, Vec3 c, Vec3 v[4], double coef[4][3]) {
  for (int k = 0; k < 4; k++) {
    const long double n2c = std::round(static_cast<long double>(v[k].z) / c.z);
    const long double n1 =
        std::round((static_cast<long double>(v[k].y) - n2c * c.y) / b.y);
    const long double n0 =
        std::round((static_cast<long double>(v[k].x) - n1 * b.x - n2c * c.x) / a.x);
    if (!(std::abs(n0) < kIntegral && std::abs(n1) < kIntegral &&
          std::abs(n2c) < kIntegral)) {
      throw std::runtime_error(
          "minimum-image lattice coefficient exceeds binary64 integer precision");
    }
    const double n[3] = {static_cast<double>(n0), static_cast<double>(n1),
                         static_cast<double>(n2c)};
    const double xs[3] = {a.x, b.x, c.x};
    const double ys[3] = {0.0, b.y, c.y};
    v[k] = {integerDot(n, xs), integerDot(n, ys), n[2] * c.z};
    for (int m = 0; m < 3; m++) {
      coef[k][m] = n[m];
    }
  }
  const Vec3 cols[3] = {v[0], v[1], v[2]};
  double inv[3][3];
  const double volume = std::abs(a.x * b.y * c.z);
  if (!(std::abs(invertColumns(cols, inv) - volume) <= 1e-6 * volume)) {
    throw std::runtime_error("minimum-image reduction left the lattice in binary64");
  }
}

inline Obtuse obtuseSuperbasis(Vec3 a, Vec3 b, Vec3 c) {
  const Vec3 original[4] = {a, b, c, scale(-1.0, add(add(a, b), c))};
  for (const Vec3 &v : original) {
    if (!std::isfinite(n2(v))) {
      throw std::runtime_error("minimum-image cell exceeds binary64 range");
    }
  }
  if (alreadyObtuse(original)) {
    const double identity[4][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, -1, -1}};
    return finish(a, b, c, original, identity);
  }
  Vec3 reduced[3] = {a, b, c};
  lagrangeReduce(reduced);
  Vec3 v[4] = {reduced[0], reduced[1], reduced[2],
               scale(-1.0, add(add(reduced[0], reduced[1]), reduced[2]))};
  sellingSteps(v);
  double coef[4][3];
  snapToLattice(a, b, c, v, coef);
  return finish(a, b, c, v, coef);
}

// Slicer steps from x over the relevant vectors, adding the integer
// coordinates of each step to n. The step with the largest gain is taken when
// it lowers the norm, binary64 or, when Precise, long double; otherwise the
// next largest gain is tried. Stops when no candidate is left, or after `cap`
// steps, and returns the number of steps taken.
template <bool Precise>
inline int descend(const Obtuse &s, Vec3 &x, double n[3], int cap) {
  // Norms only fall from here, so this bounds every projection and step.
  constexpr double kRange = std::numeric_limits<double>::max() / 16.0;
  if (!(n2(x) < kRange && s.maxRelN2 < kRange)) {
    throw std::runtime_error("minimum-image slicer projection exceeds binary64 range");
  }
  int steps = 0;
  unsigned rejected = 0;
  while (steps < cap) {
    double bestGain = 0.0;
    int at = -1;
    double sign = 0.0;
    for (int k = 0; k < 8; k++) {
      const double projection =
          x.x * s.rel[0][k] + x.y * s.rel[1][k] + x.z * s.rel[2][k];
      const double gain = std::abs(projection) - s.halfN2[k];
      if (gain > bestGain && !(rejected & (1u << k))) {
        bestGain = gain;
        at = k;
        sign = projection > 0.0 ? 1.0 : -1.0;
      }
    }
    if (at < 0) {
      break;
    }
    const Vec3 next = {x.x - sign * s.rel[0][at], x.y - sign * s.rel[1][at],
                       x.z - sign * s.rel[2][at]};
    const bool lower = Precise ? n2Long(next) < n2Long(x) : n2(next) < n2(x);
    if (!lower) {
      rejected |= 1u << at;
      continue;
    }
    rejected = 0;
    x = next;
    n[0] += sign * s.relCoef[0][at];
    n[1] += sign * s.relCoef[1][at];
    n[2] += sign * s.relCoef[2][at];
    steps++;
  }
  return steps;
}

// Babai rounding of y in the chosen superbasis triple: returns y minus the
// rounded lattice vector in binary64, its coordinates in the triple in z, and
// `scale` = sum |z_i| max|b_i|, which bounds the cancellation.
inline Vec3 babai(const Obtuse &s, Vec3 y, double z[3], double &scale) {
  double c[3];
  const double yv[3] = {y.x, y.y, y.z};
  for (int k = 0; k < 3; k++) {
    c[k] = s.inv[0][k] * yv[0] + s.inv[1][k] * yv[1] + s.inv[2][k] * yv[2];
    if (!std::isfinite(c[k])) {
      throw std::runtime_error(
          "minimum-image Babai coefficient exceeds binary64 range");
    }
    z[k] = roundAway(c[k]);
  }
  const Vec3 b0 = s.v[s.idx[0]];
  const Vec3 b1 = s.v[s.idx[1]];
  const Vec3 b2 = s.v[s.idx[2]];
  scale = std::abs(z[0]) * s.babaiMax[0] + std::abs(z[1]) * s.babaiMax[1] +
          std::abs(z[2]) * s.babaiMax[2];
  return {y.x - (z[0] * b0.x + z[1] * b1.x + z[2] * b2.x),
          y.y - (z[0] * b0.y + z[1] * b1.y + z[2] * b2.y),
          y.z - (z[0] * b0.z + z[1] * b1.z + z[2] * b2.z)};
}

// Adds the cell coordinates of the Babai triple combination z to n.
inline void addBabai(const Obtuse &s, const double z[3], double n[3]) {
  const double largest = std::max({std::abs(z[0]), std::abs(z[1]), std::abs(z[2])});
  if (!(largest * s.maxCoef < kIntegral)) {
    throw std::runtime_error(
        "minimum-image lattice coefficient exceeds binary64 integer precision");
  }
  for (int m = 0; m < 3; m++) {
    n[m] += z[0] * s.coef[s.idx[0]][m] + z[1] * s.coef[s.idx[1]][m] +
            z[2] * s.coef[s.idx[2]][m];
  }
}

inline void requireIntegral(const double n[3]) {
  if (!(std::abs(n[0]) < kIntegral && std::abs(n[1]) < kIntegral &&
        std::abs(n[2]) < kIntegral)) {
    throw std::runtime_error(
        "minimum-image lattice coefficient exceeds binary64 integer precision");
  }
}

// The image w - n . (a, b, c) without cancellation error, and the slicer run
// again on that vector until its exactly evaluated norm stops falling.
inline Vec3 refine(const Obtuse &s, Vec3 w, double n[3]) {
  requireIntegral(n);
  Vec3 best = latticeTranslate(s, w, n);
  for (;;) {
    Vec3 x = best;
    double m[3] = {n[0], n[1], n[2]};
    if (descend<true>(s, x, m, std::numeric_limits<int>::max()) == 0) {
      return best;
    }
    requireIntegral(m);
    const Vec3 next = latticeTranslate(s, w, m);
    if (!(n2Long(next) < n2Long(best))) {
      return best;
    }
    best = next;
    n[0] = m[0];
    n[1] = m[1];
    n[2] = m[2];
  }
}

constexpr int kFractionalSteps = 4;

// Binary64 results are kept only when nothing cancelled far below |w|: the
// Babai rounding subtracts lattice vectors with sum |z_i| max|b_i| at most
// 2 max|w|, the slicer takes at most two steps from a start no longer than
// the wrap (its norms only fall, and each relevant vector it applies is
// shorter than twice the current point), and the image is at least |w| / 16.
// Every other image is re-evaluated from its lattice coefficients.
inline Vec3 closestDisplacement(const Obtuse &s, Vec3 w) {
  if (!s.ok) {
    throw std::runtime_error(
        "minimum-image cell is numerically singular in binary64");
  }
  const double w2 = n2(w);
  if (!(w2 * s.maxInvN2 < kIntegral * kIntegral)) {
    throw std::runtime_error(
        "minimum-image lattice coefficient exceeds binary64 integer precision");
  }
  const auto kept = [w2](Vec3 x) { return 256.0 * n2(x) >= w2; };
  double n[3] = {0.0, 0.0, 0.0};
  Vec3 x = w;
  if (s.fractional) {
    const int steps = descend<false>(s, x, n, kFractionalSteps);
    if (steps == 0 || (steps <= 2 && kept(x))) {
      return x;
    }
    if (steps < kFractionalSteps) {
      return refine(s, w, n);
    }
  }
  double z[3];
  double scale = 0.0;
  const Vec3 rounded = babai(s, x, z, scale);
  const bool plain = scale <= 2.0 * maxAbs(w);
  double stepped[3] = {0.0, 0.0, 0.0};
  if (plain) {
    x = rounded;
    const int steps = descend<false>(s, x, stepped, std::numeric_limits<int>::max());
    if (steps <= 2 && n2(rounded) <= w2 && kept(x)) {
      return x;
    }
  }
  addBabai(s, z, n);
  requireIntegral(n);
  if (plain) {
    n[0] += stepped[0];
    n[1] += stepped[1];
    n[2] += stepped[2];
    return refine(s, w, n);
  }
  const Vec3 start = latticeTranslate(s, w, n);
  x = start;
  if (descend<false>(s, x, n, std::numeric_limits<int>::max()) == 0) {
    return start;
  }
  return refine(s, w, n);
}

struct Cached {
  double key[6]{};
  Obtuse s{};
  bool filled = false;
};

inline bool sameKey(const double *a, const double *b) {
  for (int i = 0; i < 6; i++) {
    if (a[i] != b[i]) {
      return false;
    }
  }
  return true;
}

inline const Cached &cachedBasis(double lx, double ly, double lz, double xy, double xz,
                                 double yz) {
  thread_local Cached slots[4];
  thread_local int hand = 0;
  const double key[6] = {lx, ly, lz, xy, xz, yz};
  for (const Cached &slot : slots) {
    if (slot.filled && sameKey(slot.key, key)) {
      return slot;
    }
  }
  const int i = hand;
  hand = (hand + 1) % 4;
  Cached &slot = slots[i];
  const Vec3 a{lx, 0.0, 0.0};
  const Vec3 b{xy, ly, 0.0};
  const Vec3 c{xz, yz, lz};
  slot.filled = false;
  const Obtuse s = obtuseSuperbasis(a, b, c);
  for (int k = 0; k < 6; k++) {
    slot.key[k] = key[k];
  }
  slot.s = s;
  slot.filled = true;
  return slot;
}

} // namespace

// Cartesian displacement from (xj, yj, zj) to (xi, yi, zi) on the
// recovered dump cell. `wrapped` is the fractional image of that
// difference. The fractional vector is kept when it is not strictly
// longer than the Selling image.
inline std::array<double, 3>
euclideanIfShorter(double lx, double ly, double lz, double xy, double xz, double yz,
                   const std::array<double, 3> &raw, const std::array<double, 3> &wrapped) {
  for (double value : wrapped) {
    if (!std::isfinite(value)) {
      throw std::runtime_error("minimum-image candidate is not finite");
    }
  }
  const Cached &slot = cachedBasis(lx, ly, lz, xy, xz, yz);
  (void)raw;
  const Vec3 w{wrapped[0], wrapped[1], wrapped[2]};
  const Vec3 euc = closestDisplacement(slot.s, w);
  const double f2 = n2(w);
  const double e2 = n2(euc);
  if (!std::isfinite(f2) || !std::isfinite(e2)) {
    throw std::runtime_error(
        "minimum-image squared norm exceeds binary64 range");
  }
  if (e2 < f2) {
    return {euc.x, euc.y, euc.z};
  }
  return wrapped;
}

} // namespace mic
} // namespace gen

#endif
