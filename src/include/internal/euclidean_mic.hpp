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
// Sommer–Feder–Shalvi slicer.
//
// None of the reduction loops has an iteration cap. Every accepted Lagrange,
// Selling, or slicer step must strictly lower a non-negative, finite squared
// norm, evaluated in long double over binary64 vectors. There are finitely
// many binary64 vectors below the initial norm, which is a (very loose)
// termination bound. A slicer gain below long-double resolution is converged
// at that numeric precision. Other failures to decrease, and bases that cannot
// be represented as full rank, throw instead of returning an unchecked image.

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

inline void relevantOf(const Vec3 v[4], Obtuse &o) {
  const Vec3 vecs[8] = {v[0],
                        v[1],
                        v[2],
                        v[3],
                        add(v[0], v[1]),
                        add(v[0], v[2]),
                        add(v[0], v[3]),
                        {}};
  for (int k = 0; k < 8; k++) {
    o.rel[0][k] = vecs[k].x;
    o.rel[1][k] = vecs[k].y;
    o.rel[2][k] = vecs[k].z;
    o.halfN2[k] =
        (k < 7) ? 0.5 * n2(vecs[k]) : std::numeric_limits<double>::infinity();
  }
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
          if (!std::isfinite(after) || after > before) {
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
    for (int i = 0; i < 3 && si < 0; i++) {
      for (int j = i + 1; j < 4; j++) {
        if (dot(v[i], v[j]) > 0.0) {
          si = i;
          sj = j;
          break;
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
    if (!std::isfinite(after) || !(after < before)) {
      throw std::runtime_error(
          "minimum-image Selling reduction reached binary64 limits");
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
  relevantOf(v, o);
  o.ok = have && bestAbs > 1e-18;
  return o;
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
      const long double resolution =
          16.0L * std::numeric_limits<long double>::epsilon() *
          std::max(before, 1.0L);
      if (static_cast<long double>(bestGain) <= resolution) {
        return x;
      }
      throw std::runtime_error("minimum-image slicer reached binary64 limits");
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

inline bool engineIsShort(Vec3 a, Vec3 b, Vec3 c) {
  const double scaleLen = std::max(std::abs(a.x) + std::abs(b.y) + std::abs(c.z), 1.0);
  const double tol = 1e-10 * scaleLen;
  return std::abs(a.y) <= tol && std::abs(a.z) <= tol && std::abs(b.z) <= tol &&
         std::abs(b.x) <= 0.5 * std::abs(a.x) + tol &&
         std::abs(c.x) <= 0.5 * std::abs(a.x) + tol &&
         std::abs(c.y) <= 0.5 * std::abs(b.y) + tol;
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
  for (int k = 0; k < 6; k++) {
    slot.key[k] = key[k];
  }
  const Vec3 a{lx, 0.0, 0.0};
  const Vec3 b{xy, ly, 0.0};
  const Vec3 c{xz, yz, lz};
  slot.s = obtuseSuperbasis(a, b, c);
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
