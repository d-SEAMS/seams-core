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
// lattice point is Babai rounding in a Selling superbasis, followed by the
// Sommer-Feder-Shalvi slicer over the Voronoi-relevant vectors.
//
// The relevant vectors are found once per cell. Every Voronoi-relevant vector
// is, up to sign, the shortest member of one of the seven nonzero cosets of
// L / 2L, and Fincke-Pohst enumeration finds each coset's shortest member.
// A point that none of them shortens lies in the Voronoi cell of the origin,
// so the slicer's stopping point is the minimum image.
// The reduced vectors are re-evaluated from their integer coordinates in the
// cell before any of this, so reduction rounding does not leave the lattice.
//
// None of the loops has an iteration cap. Every accepted Lagrange, Selling,
// or slicer step strictly lowers a non-negative, finite squared norm over
// binary64 vectors, which bounds the step count. A step whose rounded norm
// does not fall changes the squared distance by no more than binary64
// rounding of |x|^2, so the loop stops there and keeps the current vector.
// Non-finite values, singular bases, and lattice coefficients at or above
// 2^52 throw explicitly without poisoning the per-thread basis cache.

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

// Half away from zero, matching a truncation after adding just under 1/2.
inline double roundAway(double x) {
  constexpr double kIntegral = 4503599627370496.0;
  if (!std::isfinite(x) || std::abs(x) >= kIntegral) {
    throw std::runtime_error(
        "minimum-image lattice coefficient exceeds binary64 integer precision");
  }
  const double t = x + std::copysign(0.49999999999999994, x);
  return static_cast<double>(static_cast<long long>(t));
}

struct Obtuse {
  Vec3 v[4]{};
  double rel[3][8]{};
  double halfN2[8]{};
  double inv[3][3]{};
  int idx[3] = {0, 1, 2};
  bool ok = false;
};

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
// each coset's shortest {-1, 0, 1}^3 member as its radius. Returns false when
// the Gram-Schmidt norms are not positive and finite.
inline bool cosetMinima(const Vec3 basis[3], Vec3 out[7]) {
  constexpr long double kWidth = 4503599627370496.0L;
  Vec3 b[3] = {basis[0], basis[1], basis[2]};
  std::sort(b, b + 3, [](Vec3 p, Vec3 q) { return n2(p) < n2(q); });
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

  Vec3 best[8]{};
  long double bestN2[8];
  for (long double &n : bestN2) {
    n = std::numeric_limits<long double>::infinity();
  }
  for (int z0 = -1; z0 <= 1; z0++) {
    for (int z1 = -1; z1 <= 1; z1++) {
      for (int z2 = -1; z2 <= 1; z2++) {
        const int coset = (z0 & 1) | ((z1 & 1) << 1) | ((z2 & 1) << 2);
        if (coset == 0) {
          continue;
        }
        const Vec3 u = combination(b, z0, z1, z2);
        const long double n = n2Long(u);
        if (n < bestN2[coset]) {
          bestN2[coset] = n;
          best[coset] = u;
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
          const Vec3 u = combination(b, z0, z1, z2);
          const long double n = n2Long(u);
          if (n < bestN2[coset]) {
            bestN2[coset] = n;
            best[coset] = u;
          }
        }
      }
    }
  }

  for (int k = 0; k < 7; k++) {
    out[k] = best[k + 1];
  }
  return true;
}

inline double roundingSlack(Vec3 a, Vec3 b) {
  return 16.0 * std::numeric_limits<double>::epsilon() * std::sqrt(n2(a) * n2(b));
}

inline bool relevantOf(const Vec3 v[4], const int idx[3], Obtuse &o) {
  const Vec3 basis[3] = {v[idx[0]], v[idx[1]], v[idx[2]]};
  Vec3 u[7];
  if (!cosetMinima(basis, u)) {
    return false;
  }
  for (int k = 0; k < 8; k++) {
    const Vec3 r = (k < 7) ? u[k] : Vec3{};
    o.rel[0][k] = r.x;
    o.rel[1][k] = r.y;
    o.rel[2][k] = r.z;
    o.halfN2[k] = (k < 7) ? 0.5 * n2(r) : std::numeric_limits<double>::infinity();
    if (k < 7 && !std::isfinite(o.halfN2[k])) {
      return false;
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

inline Obtuse finish(const Vec3 v[4]) {
  Obtuse o;
  for (int i = 0; i < 4; i++) {
    o.v[i] = v[i];
  }
  double bestAbs = 0.0;
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
    if (absDet > bestAbs) {
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
  o.ok = relevantOf(v, o.idx, o);
  return o;
}

// sum n[i] * x[i] for integers n[i] below 2^52, with the products and partial
// sums carried error-free (Ogita, Rump and Oishi, Dot2), so the result is
// within about one rounding of itself however large the cancelling terms.
inline double integerDot(const double n[3], const double x[3]) {
  double sum = 0.0;
  double error = 0.0;
  for (int i = 0; i < 3; i++) {
    const double product = n[i] * x[i];
    const double productError = std::fma(n[i], x[i], -product);
    const double next = sum + product;
    const double back = next - sum;
    error += productError + ((sum - (next - back)) + (product - back));
    sum = next;
  }
  return sum + error;
}

// Lagrange and Selling form v - q w in binary64, and with q near 1e7 that
// rounding moves a reduced vector off the lattice by far more than one ulp of
// its own length. Each vector is re-evaluated from its integer coordinates in
// the upper-triangular cell (a, b, c), recovered by back-substitution. The
// result must still span the cell volume.
inline void snapToLattice(Vec3 a, Vec3 b, Vec3 c, Vec3 v[4]) {
  for (int k = 0; k < 4; k++) {
    const long double n2c = std::round(static_cast<long double>(v[k].z) / c.z);
    const long double n1 =
        std::round((static_cast<long double>(v[k].y) - n2c * c.y) / b.y);
    const long double n0 =
        std::round((static_cast<long double>(v[k].x) - n1 * b.x - n2c * c.x) / a.x);
    if (!(std::abs(n0) < 4503599627370496.0L && std::abs(n1) < 4503599627370496.0L &&
          std::abs(n2c) < 4503599627370496.0L)) {
      throw std::runtime_error(
          "minimum-image lattice coefficient exceeds binary64 integer precision");
    }
    const double n[3] = {static_cast<double>(n0), static_cast<double>(n1),
                         static_cast<double>(n2c)};
    const double xs[3] = {a.x, b.x, c.x};
    const double ys[3] = {0.0, b.y, c.y};
    v[k] = {integerDot(n, xs), integerDot(n, ys), n[2] * c.z};
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
    return finish(original);
  }
  Vec3 reduced[3] = {a, b, c};
  lagrangeReduce(reduced);
  Vec3 v[4] = {reduced[0], reduced[1], reduced[2],
               scale(-1.0, add(add(reduced[0], reduced[1]), reduced[2]))};
  sellingSteps(v);
  snapToLattice(a, b, c, v);
  return finish(v);
}

inline Vec3 descend(const Obtuse &s, Vec3 x) {
  for (;;) {
    double bestGain = 0.0;
    int at = -1;
    double sign = 0.0;
    for (int k = 0; k < 8; k++) {
      const double projection =
          x.x * s.rel[0][k] + x.y * s.rel[1][k] + x.z * s.rel[2][k];
      if (!std::isfinite(projection)) {
        throw std::runtime_error(
            "minimum-image slicer projection exceeds binary64 range");
      }
      const double gain = std::abs(projection) - s.halfN2[k];
      if (gain > bestGain) {
        bestGain = gain;
        at = k;
        sign = projection > 0.0 ? 1.0 : -1.0;
      }
    }
    if (at < 0) {
      return x;
    }
    const Vec3 next = {x.x - sign * s.rel[0][at], x.y - sign * s.rel[1][at],
                       x.z - sign * s.rel[2][at]};
    const long double before = n2Long(x);
    const long double after = n2Long(next);
    if (!std::isfinite(n2(next))) {
      throw std::runtime_error("minimum-image slicer reached binary64 limits");
    }
    if (!(after < before)) {
      return x;
    }
    x = next;
  }
}

inline Vec3 babai(const Obtuse &s, Vec3 y) {
  double c[3];
  const double yv[3] = {y.x, y.y, y.z};
  for (int k = 0; k < 3; k++) {
    c[k] = s.inv[0][k] * yv[0] + s.inv[1][k] * yv[1] + s.inv[2][k] * yv[2];
    if (!std::isfinite(c[k])) {
      throw std::runtime_error(
          "minimum-image Babai coefficient exceeds binary64 range");
    }
  }
  const double z0 = roundAway(c[0]);
  const double z1 = roundAway(c[1]);
  const double z2 = roundAway(c[2]);
  const Vec3 b0 = s.v[s.idx[0]];
  const Vec3 b1 = s.v[s.idx[1]];
  const Vec3 b2 = s.v[s.idx[2]];
  return {y.x - (z0 * b0.x + z1 * b1.x + z2 * b2.x),
          y.y - (z0 * b0.y + z1 * b1.y + z2 * b2.y),
          y.z - (z0 * b0.z + z1 * b1.z + z2 * b2.z)};
}

inline Vec3 closestDisplacement(const Obtuse &s, Vec3 y) {
  if (!s.ok) {
    throw std::runtime_error(
        "minimum-image cell is numerically singular in binary64");
  }
  return descend(s, babai(s, y));
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
