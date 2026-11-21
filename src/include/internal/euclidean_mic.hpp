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
// Voronoi-relevant vectors, started from the wrap itself when probes along
// the cell edges and corners settle within the wrap step cap, and otherwise
// from Babai rounding in a Selling superbasis.
//
// The relevant vectors are found once per cell. Every Voronoi-relevant vector
// is, up to sign, the shortest member of one of the seven nonzero cosets of
// L / 2L. A superbasis that is obtuse beyond rounding gives those as its
// subset sums; otherwise LLL-reduced Fincke-Pohst enumeration finds each
// coset's shortest member. A point that none of them shortens lies in the
// Voronoi cell of the origin, so the slicer's stopping point is the minimum
// image. The reduced and relevant vectors are re-evaluated from their integer
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

inline long double dotLong(Vec3 a, Vec3 b) {
  return static_cast<long double>(a.x) * b.x +
         static_cast<long double>(a.y) * b.y +
         static_cast<long double>(a.z) * b.z;
}

// Evaluates |a|^2 < |b|^2 as (a - b) . (a + b) < 0 in extended precision so
// steps along a short basis vector do not cancel against a large orthogonal
// component.
inline bool shorterLong(Vec3 a, Vec3 b) {
  const long double dx = static_cast<long double>(a.x) - b.x;
  const long double dy = static_cast<long double>(a.y) - b.y;
  const long double dz = static_cast<long double>(a.z) - b.z;
  return dx * (static_cast<long double>(a.x) + b.x) +
         dy * (static_cast<long double>(a.y) + b.y) +
         dz * (static_cast<long double>(a.z) + b.z) < 0.0L;
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
  // Hot wrap-path fields first so they stay contiguous in L1 cache.
  double rel[3][8]{};
  double halfN2[8]{};
  // dot(rel_i, rel_j), so a binary64 slicer step can update every
  // projection with one subtraction instead of a fresh dot product.
  double gram[8][8]{};
  // Largest squared row of inv: |Babai coefficient|^2 <= n2(y) * maxInvN2.
  double maxInvN2 = 0.0;
  bool fractional = false;
  bool ok = false;
  Vec3 v[4]{};
  // Integer coordinates of v[k] and of each relevant vector in (a, b, c).
  double coef[4][3]{};
  double relCoef[3][8]{};
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
  // max |component| of each vector of the Babai triple.
  double babaiMax[3]{};
  int idx[3] = {0, 1, 2};
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

inline long long firstWithParity(long double low, int parity) {
  long long n = static_cast<long long>(std::ceil(low));
  if (((n % 2) + 2) % 2 != parity) {
    n++;
  }
  return n;
}

inline long long nearestWithParity(long double c, int parity) {
  const long double shifted = (c - static_cast<long double>(parity)) * 0.5L;
  return 2LL * static_cast<long long>(std::llround(shifted)) + parity;
}

// Shortest member of each nonzero coset of L / 2L, the parity class of the
// coordinates in `basis`. The basis is first 3D LLL-reduced in long double
// with an exact integer unimodular transform (which preserves L / 2L), so
// the Gram-Schmidt norms satisfy q0 <= (4/3) q1 and q1 <= (4/3) q2.
// Enumeration then compares candidates by exact Gram-Schmidt norm
// differences without adding and subtracting q2 * z2^2, and minimizes z0
// directly, keeping the work O(1) independent of aspect ratio.
inline bool cosetMinima(const Vec3 basis[3], const double basisCoef[3][3], long long out[7][3]) {
  constexpr long double kWidth = kIntegral;
  constexpr long double kEps = std::numeric_limits<long double>::epsilon();
  long long U[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  long double B[3][3]{};
  long double Bs[3][3]{};
  long double mu[3][3]{};
  long double q[3]{};
  const auto dot3 = [](const long double a[3], const long double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  };
  const auto evalB = [&]() {
    for (int i = 0; i < 3; i++) {
      const long double u0 = static_cast<long double>(U[i][0]);
      const long double u1 = static_cast<long double>(U[i][1]);
      const long double u2 = static_cast<long double>(U[i][2]);
      B[i][0] = u0 * basis[0].x + u1 * basis[1].x + u2 * basis[2].x;
      B[i][1] = u0 * basis[0].y + u1 * basis[1].y + u2 * basis[2].y;
      B[i][2] = u0 * basis[0].z + u1 * basis[1].z + u2 * basis[2].z;
    }
  };
  const auto gso = [&]() -> bool {
    for (int i = 0; i < 3; i++) {
      for (int m = 0; m < 3; m++) {
        Bs[i][m] = B[i][m];
      }
      for (int j = 0; j < i; j++) {
        mu[i][j] = dot3(B[i], Bs[j]) / q[j];
        if (!std::isfinite(mu[i][j])) {
          return false;
        }
        for (int m = 0; m < 3; m++) {
          Bs[i][m] -= mu[i][j] * Bs[j][m];
        }
      }
      for (int j = 0; j < i; j++) {
        const long double corr = dot3(Bs[i], Bs[j]) / q[j];
        if (!std::isfinite(corr)) {
          return false;
        }
        mu[i][j] += corr;
        for (int m = 0; m < 3; m++) {
          Bs[i][m] -= corr * Bs[j][m];
        }
      }
      q[i] = dot3(Bs[i], Bs[i]);
      if (!(std::isfinite(q[i]) && q[i] > 0.0L)) {
        return false;
      }
    }
    return true;
  };

  evalB();
  if (!gso()) {
    return false;
  }
  int k = 1;
  int guard = 0;
  while (k < 3 && guard++ < 1000) {
    for (int j = k - 1; j >= 0; j--) {
      if (std::abs(mu[k][j]) > 0.5001L) {
        const long double rq = std::round(mu[k][j]);
        if (!(std::abs(rq) < kWidth)) {
          return false;
        }
        const long long qi = static_cast<long long>(rq);
        if (qi != 0) {
          for (int m = 0; m < 3; m++) {
            const long double next =
                static_cast<long double>(U[k][m]) - rq * static_cast<long double>(U[j][m]);
            if (!(std::abs(next) < kWidth)) {
              return false;
            }
            U[k][m] -= qi * U[j][m];
          }
          evalB();
          if (!gso()) {
            return false;
          }
        }
      }
    }
    if (q[k] >= (0.99L - mu[k][k - 1] * mu[k][k - 1]) * q[k - 1]) {
      k++;
    } else {
      for (int m = 0; m < 3; m++) {
        std::swap(U[k][m], U[k - 1][m]);
      }
      evalB();
      if (!gso()) {
        return false;
      }
      k = k > 1 ? k - 1 : 1;
    }
  }
  if (guard >= 1000) {
    return false;
  }

  long double C[3][3];
  for (int i = 0; i < 3; i++) {
    for (int m = 0; m < 3; m++) {
      C[i][m] = static_cast<long double>(U[i][0]) * basisCoef[0][m] +
                static_cast<long double>(U[i][1]) * basisCoef[1][m] +
                static_cast<long double>(U[i][2]) * basisCoef[2][m];
    }
  }
  const auto cellL1 = [&](long long z0, long long z1, long long z2) -> long double {
    long double s = 0.0L;
    for (int m = 0; m < 3; m++) {
      s += std::abs(static_cast<long double>(z0) * C[0][m] +
                    static_cast<long double>(z1) * C[1][m] +
                    static_cast<long double>(z2) * C[2][m]);
    }
    return s;
  };

  for (int coset = 1; coset < 8; coset++) {
    const int p0 = coset & 1;
    const int p1 = (coset >> 1) & 1;
    const int p2 = (coset >> 2) & 1;

    long long z2Best = p2;
    const long double c1Seed = -mu[2][1] * static_cast<long double>(z2Best);
    if (!(std::abs(c1Seed) < kWidth)) {
      return false;
    }
    long long z1Best = nearestWithParity(c1Seed, p1);
    long double e1Best = static_cast<long double>(z1Best) - c1Seed;
    const long double c0Seed = -(mu[1][0] * static_cast<long double>(z1Best) +
                                 mu[2][0] * static_cast<long double>(z2Best));
    if (!(std::abs(c0Seed) < kWidth)) {
      return false;
    }
    long long z0Best = nearestWithParity(c0Seed, p0);
    long double e0Best = static_cast<long double>(z0Best) - c0Seed;
    long double l1Best = cellL1(z0Best, z1Best, z2Best);

    const auto remNorm2 = [&]() -> long double {
      return q[1] * e1Best * e1Best + q[0] * e0Best * e0Best;
    };
    const auto z2Limit = [&]() -> long double {
      const long double rem = remNorm2();
      return std::sqrt(static_cast<long double>(z2Best * z2Best) +
                       (rem / q[2]) * (1.0L + 1e-9L));
    };
    long double w2 = z2Limit();
    if (!(std::isfinite(w2) && w2 < kWidth)) {
      return false;
    }
    // Sign symmetry (-z0, -z1, -z2) preserves both parity and norm, so z2 >= 0 suffices.
    for (long long z2 = p2; static_cast<long double>(z2) <= w2; z2 += 2) {
      const long long dz2 = z2 * z2 - z2Best * z2Best;
      const long double rem = remNorm2();
      const long double budget1 =
          (dz2 == 0 ? 0.0L : -q[2] * static_cast<long double>(dz2)) +
          rem * (1.0L + 64.0L * kEps);
      if (budget1 < 0.0L) {
        continue;
      }
      const long double c1 = -mu[2][1] * static_cast<long double>(z2);
      const long double w1 = std::sqrt((budget1 / q[1]) * (1.0L + 1e-9L));
      if (!(std::isfinite(w1) && std::abs(c1) + w1 < kWidth)) {
        return false;
      }
      for (long long z1 = firstWithParity(c1 - w1, p1);
           static_cast<long double>(z1) <= c1 + w1; z1 += 2) {
        const long double e1 = static_cast<long double>(z1) - c1;
        const long double c0 = -(mu[1][0] * static_cast<long double>(z1) +
                                 mu[2][0] * static_cast<long double>(z2));
        if (!(std::abs(c0) < kWidth)) {
          return false;
        }
        const long long z0Near = nearestWithParity(c0, p0);
        for (long long off = -2; off <= 2; off += 2) {
          const long long z0 = z0Near + off;
          const long double e0 = static_cast<long double>(z0) - c0;
          const long double diff =
              (z2 == z2Best ? 0.0L
                            : q[2] * static_cast<long double>(z2 * z2 - z2Best * z2Best)) +
              q[1] * (e1 * e1 - e1Best * e1Best) +
              q[0] * (e0 * e0 - e0Best * e0Best);
          const long double diffScale =
              (z2 == z2Best ? 0.0L
                            : q[2] * static_cast<long double>(z2 * z2 + z2Best * z2Best)) +
              q[1] * (e1 * e1 + e1Best * e1Best) +
              q[0] * (e0 * e0 + e0Best * e0Best);
          const long double curTol = 64.0L * kEps * diffScale;
          bool better = diff < -curTol;
          long double l1Cand = 0.0L;
          if (better) {
            l1Cand = cellL1(z0, z1, z2);
          } else if (diff <= curTol) {
            l1Cand = cellL1(z0, z1, z2);
            better = l1Cand < l1Best - 0.5L;
          }
          if (better) {
            z0Best = z0;
            z1Best = z1;
            z2Best = z2;
            e1Best = e1;
            e0Best = e0;
            l1Best = l1Cand;
            w2 = z2Limit();
          }
        }
      }
    }

    for (int m = 0; m < 3; m++) {
      const long double cm = static_cast<long double>(z0Best) * U[0][m] +
                             static_cast<long double>(z1Best) * U[1][m] +
                             static_cast<long double>(z2Best) * U[2][m];
      if (!(std::abs(cm) < kWidth)) {
        return false;
      }
      out[coset - 1][m] = z0Best * U[0][m] + z1Best * U[1][m] + z2Best * U[2][m];
    }
  }
  return true;
}

inline double roundingSlack(Vec3 a, Vec3 b) {
  return 16.0 * std::numeric_limits<double>::epsilon() * std::sqrt(n2(a) * n2(b));
}

// A superbasis whose dot products are all non-positive beyond rounding is
// obtuse (including exact-zero orthogonal pairs), and its seven subset sums
// are coset minima of L / 2L containing every Voronoi-relevant vector
// (Conway and Sloane). Otherwise the coset minima are enumerated.
// Each relevant vector is evaluated from its integer coordinates in the cell.
inline bool relevantOf(Obtuse &o) {
  const Vec3 *v = o.v;
  bool obtuse = true;
  for (int m = 0; m < 3; m++) {
    if (o.coef[0][m] + o.coef[1][m] + o.coef[2][m] + o.coef[3][m] != 0.0) {
      obtuse = false;
      break;
    }
  }
  for (int i = 0; i < 3 && obtuse; i++) {
    for (int j = i + 1; j < 4; j++) {
      const long double px = static_cast<long double>(v[i].x) * v[j].x;
      const long double py = static_cast<long double>(v[i].y) * v[j].y;
      const long double pz = static_cast<long double>(v[i].z) * v[j].z;
      const long double d = px + py + pz;
      const long double absSum = std::abs(px) + std::abs(py) + std::abs(pz);
      if (!(d <= -16.0L * std::numeric_limits<long double>::epsilon() * absSum)) {
        obtuse = false;
        break;
      }
    }
  }
  double coef[7][3];
  if (obtuse) {
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
    double basisCoef[3][3];
    for (int i = 0; i < 3; i++) {
      for (int m = 0; m < 3; m++) {
        basisCoef[i][m] = o.coef[o.idx[i]][m];
      }
    }
    long long z[7][3];
    if (!cosetMinima(basis, basisCoef, z)) {
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
      if (!(std::isfinite(o.halfN2[k]) &&
            o.halfN2[k] >= std::numeric_limits<double>::min())) {
        return false;
      }
      o.maxRelN2 = std::max(o.maxRelN2, n2(r));
    }
  }
  if (!(o.maxRelN2 < std::numeric_limits<double>::max() / 64.0)) {
    return false;
  }
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      o.gram[i][j] =
          o.rel[0][i] * o.rel[0][j] + o.rel[1][i] * o.rel[1][j] + o.rel[2][i] * o.rel[2][j];
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
        if (!(lj > 0.0)) {
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

inline long double detLong(Vec3 a, Vec3 b, Vec3 c) {
  const long double ax = a.x, ay = a.y, az = a.z;
  const long double bx = b.x, by = b.y, bz = b.z;
  const long double cx = c.x, cy = c.y, cz = c.z;
  return ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) +
         az * (bx * cy - by * cx);
}

// Inverse of the matrix with columns cols[0..2]. inv is stored so that
// coeff k is sum_i inv[i][k] * y[i]. Returns the scale-invariant
// orthogonality |det| / (|a| |b| |c|), or 0 when singular.
inline double invertColumns(const Vec3 cols[3], double inv[3][3]) {
  const Vec3 &a = cols[0];
  const Vec3 &b = cols[1];
  const Vec3 &c = cols[2];
  const long double ax = a.x, ay = a.y, az = a.z;
  const long double bx = b.x, by = b.y, bz = b.z;
  const long double cx = c.x, cy = c.y, cz = c.z;
  const long double det = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) +
                          az * (bx * cy - by * cx);
  const long double normProd =
      std::sqrt(n2Long(a)) * std::sqrt(n2Long(b)) * std::sqrt(n2Long(c));
  if (!(std::isfinite(det) && std::isfinite(normProd) && normProd > 0.0L)) {
    return 0.0;
  }
  const long double score = std::abs(det) / normProd;
  if (!(std::isfinite(score) && score > 1e-18L)) {
    return 0.0;
  }
  const long double cof[3][3] = {
      {by * cz - bz * cy, az * cy - ay * cz, ay * bz - az * by},
      {bz * cx - bx * cz, ax * cz - az * cx, az * bx - ax * bz},
      {bx * cy - by * cx, ay * cx - ax * cy, ax * by - ay * bx}};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv[i][j] = static_cast<double>(cof[i][j] / det);
      if (!std::isfinite(inv[i][j])) {
        return 0.0;
      }
    }
  }
  return static_cast<double>(score);
}

constexpr int kFractionalSteps = 4;

inline int descendFast(const Obtuse &s, Vec3 &x, int cap);

inline Obtuse finish(Vec3 a, Vec3 b, Vec3 c, const Vec3 v[4], const double coef[4][3]) {
  Obtuse o;
  const Vec3 cell[3] = {a, b, c};
  double maxEdgeN2 = 0.0;
  for (int k = 0; k < 3; k++) {
    const double edge = n2(cell[k]);
    if (!(std::isfinite(edge) && edge >= std::numeric_limits<double>::min())) {
      return o;
    }
    maxEdgeN2 = std::max(maxEdgeN2, edge);
    o.cellX[k] = cell[k].x;
    o.cellY[k] = cell[k].y;
    o.cellZ[k] = cell[k].z;
    o.splitX[k] = split(cell[k].x);
    o.splitY[k] = split(cell[k].y);
    o.splitZ[k] = split(cell[k].z);
  }
  for (int i = 0; i < 4; i++) {
    o.v[i] = v[i];
    for (int m = 0; m < 3; m++) {
      o.coef[i][m] = coef[i][m];
    }
  }
  // Every triple of a superbasis spans the same cell volume, so Babai takes
  // the most orthogonal one, |det| / (|b0| |b1| |b2|).
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
    const double score = invertColumns(cols, inv);
    if (score > bestScore) {
      bestScore = score;
      have = true;
      for (int k = 0; k < 3; k++) {
        o.idx[k] = idx[k];
        for (int j = 0; j < 3; j++) {
          o.inv[k][j] = inv[k][j];
        }
      }
    }
  }
  if (!have || !(bestScore > 1e-18)) {
    return o;
  }
  for (int k = 0; k < 3; k++) {
    o.babaiMax[k] = maxAbs(v[o.idx[k]]);
    o.maxInvN2 = std::max(o.maxInvN2, o.inv[0][k] * o.inv[0][k] + o.inv[1][k] * o.inv[1][k] +
                                          o.inv[2][k] * o.inv[2][k]);
  }
  if (!(std::isfinite(o.maxInvN2) && o.maxInvN2 > 0.0)) {
    return o;
  }
  o.ok = relevantOf(o);
  // Scale-invariant wrap-path gate: directly probe 0.375 * each cell edge and
  // 0.375 * the four antipodal cell corners with descendFast. Compact cells and
  // arbitrarily long near-orthogonal slabs settle in 0-3 steps; skewed cells
  // where wrap points routinely hit the step cap go straight to Babai.
  constexpr double kRange = std::numeric_limits<double>::max() / 16.0;
  if (o.ok && maxEdgeN2 < kRange) {
    bool fits = true;
    for (int j = 0; j < 3; j++) {
      Vec3 probe = scale(0.375, cell[j]);
      if (descendFast(o, probe, kFractionalSteps) >= kFractionalSteps) {
        fits = false;
        break;
      }
    }
    if (fits) {
      const Vec3 corners[4] = {
          {0.375 * (a.x + b.x + c.x), 0.375 * (a.y + b.y + c.y), 0.375 * (a.z + b.z + c.z)},
          {0.375 * (a.x + b.x - c.x), 0.375 * (a.y + b.y - c.y), 0.375 * (a.z + b.z - c.z)},
          {0.375 * (a.x - b.x + c.x), 0.375 * (a.y - b.y + c.y), 0.375 * (a.z - b.z + c.z)},
          {0.375 * (a.x - b.x - c.x), 0.375 * (a.y - b.y - c.y), 0.375 * (a.z - b.z - c.z)},
      };
      int cornerHits = 0;
      for (int j = 0; j < 4; j++) {
        Vec3 probe = corners[j];
        if (descendFast(o, probe, kFractionalSteps) >= kFractionalSteps) {
          cornerHits++;
        }
      }
      fits = cornerHits <= 1;
    }
    o.fractional = fits;
  }
  return o;
}

// Lagrange and Selling form v - q w in binary64, and with q near 1e7 that
// rounding moves a reduced vector off the lattice by far more than one ulp of
// its own length. Each vector is re-evaluated from its integer coordinates in
// the upper-triangular cell (a, b, c), recovered by back-substitution and
// returned in coef. The result must still span the cell volume.
inline void snapToLattice(Vec3 a, Vec3 b, Vec3 c, Vec3 v[4], double coef[4][3]) {
  const double xs[3] = {a.x, b.x, c.x};
  const double ys[3] = {0.0, b.y, c.y};
  const auto evalFromCoef = [&](const double (&n)[3]) -> Vec3 {
    if (!(std::abs(n[0]) < kIntegral && std::abs(n[1]) < kIntegral &&
          std::abs(n[2]) < kIntegral)) {
      throw std::runtime_error(
          "minimum-image lattice coefficient exceeds binary64 integer precision");
    }
    return {integerDot(n, xs), integerDot(n, ys), n[2] * c.z};
  };
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
    v[k] = evalFromCoef(n);
    for (int m = 0; m < 3; m++) {
      coef[k][m] = n[m];
    }
  }
  const long double volume =
      std::abs(static_cast<long double>(a.x) * b.y * c.z);
  const long double detV = std::abs(detLong(v[0], v[1], v[2]));
  if (!(volume > 0.0L && std::isfinite(volume) && std::isfinite(detV) &&
        std::abs(detV - volume) <= 1e-6L * volume)) {
    throw std::runtime_error("minimum-image reduction left the lattice in binary64");
  }
  // Exact-integer Lagrange size-reduction in extended precision: if binary64
  // rounding in Lagrange/Selling shifted a long reduced vector by many short
  // vectors during back-substitution, pull it back to mutual orthogonality
  // without ever leaving the lattice.
  for (int pass = 0; pass < 32; pass++) {
    for (int i = 1; i < 3; i++) {
      int j = i;
      while (j > 0 && n2Long(v[j]) < n2Long(v[j - 1])) {
        std::swap(v[j], v[j - 1]);
        for (int m = 0; m < 3; m++) {
          std::swap(coef[j][m], coef[j - 1][m]);
        }
        j--;
      }
    }
    bool changed = false;
    for (int i = 1; i < 3; i++) {
      for (int j = i - 1; j >= 0; j--) {
        const long double lj = n2Long(v[j]);
        if (!(lj > 0.0L)) {
          continue;
        }
        const long double qd = std::round(dotLong(v[i], v[j]) / lj);
        if (!std::isfinite(qd) || qd == 0.0L) {
          continue;
        }
        double nextCoef[3];
        for (int m = 0; m < 3; m++) {
          const long double nc =
              static_cast<long double>(coef[i][m]) - qd * static_cast<long double>(coef[j][m]);
          if (!(std::abs(nc) < kIntegral)) {
            throw std::runtime_error(
                "minimum-image lattice coefficient exceeds binary64 integer precision");
          }
          nextCoef[m] = static_cast<double>(nc);
        }
        const Vec3 next = evalFromCoef(nextCoef);
        if (shorterLong(next, v[i])) {
          v[i] = next;
          for (int m = 0; m < 3; m++) {
            coef[i][m] = nextCoef[m];
          }
          changed = true;
        }
      }
    }
    if (!changed) {
      break;
    }
  }
  for (int m = 0; m < 3; m++) {
    coef[3][m] = -(coef[0][m] + coef[1][m] + coef[2][m]);
  }
  v[3] = evalFromCoef(coef[3]);
  // Keep the superbasis obtuse after any integer Lagrange cleanup.
  for (int iter = 0; iter < 64; iter++) {
    int pi = -1;
    int pj = -1;
    long double best = 0.0L;
    for (int i = 0; i < 3; i++) {
      for (int j = i + 1; j < 4; j++) {
        const long double d = dotLong(v[i], v[j]);
        if (d > best) {
          best = d;
          pi = i;
          pj = j;
        }
      }
    }
    if (pi < 0) {
      break;
    }
    double nextCoef[4][3];
    Vec3 nextV[4];
    for (int k = 0; k < 4; k++) {
      for (int m = 0; m < 3; m++) {
        if (k == pi) {
          nextCoef[k][m] = -coef[pi][m];
        } else if (k == pj) {
          nextCoef[k][m] = coef[pj][m];
        } else {
          nextCoef[k][m] = coef[k][m] + coef[pi][m];
        }
      }
      nextV[k] = evalFromCoef(nextCoef[k]);
    }
    long double before = 0.0L;
    long double after = 0.0L;
    for (int k = 0; k < 4; k++) {
      before += n2Long(v[k]);
      after += n2Long(nextV[k]);
    }
    if (!(after < before)) {
      break;
    }
    for (int k = 0; k < 4; k++) {
      v[k] = nextV[k];
      for (int m = 0; m < 3; m++) {
        coef[k][m] = nextCoef[k][m];
      }
    }
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
  double proj[8];
  const auto project = [&]() {
    for (int k = 0; k < 8; k++) {
      proj[k] = x.x * s.rel[0][k] + x.y * s.rel[1][k] + x.z * s.rel[2][k];
    }
  };
  project();
  int steps = 0;
  unsigned rejected = 0;
  while (steps < cap) {
    double bestGain = 0.0;
    int at = -1;
    double sign = 0.0;
    for (int k = 0; k < 8; k++) {
      const double gain = std::abs(proj[k]) - s.halfN2[k];
      if (gain > bestGain && !(rejected & (1u << k))) {
        bestGain = gain;
        at = k;
        sign = proj[k] > 0.0 ? 1.0 : -1.0;
      }
    }
    if (at < 0) {
      break;
    }
    const Vec3 next = {x.x - sign * s.rel[0][at], x.y - sign * s.rel[1][at],
                       x.z - sign * s.rel[2][at]};
    const bool lower = shorterLong(next, x);
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
    // A binary64 step moves x by exactly sign * rel[at], so the stored
    // Gram row updates every projection. The extended-precision norm does
    // not, and a rejected step leaves x where it was.
    if constexpr (Precise) {
      project();
    } else {
      for (int k = 0; k < 8; k++) {
        proj[k] -= sign * s.gram[at][k];
      }
    }
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
// Each stored coefficient is an integer, and z[k] multiplies only row
// coef[idx[k]]. The product that has to be an exact binary64 integer is
// that term, not max|z| times the largest coefficient of a different
// vector. A term at or above 2^52 is rejected; requireIntegral then
// rejects a sum that reaches the same line.
inline void addBabai(const Obtuse &s, const double z[3], double n[3]) {
  for (int m = 0; m < 3; m++) {
    double sum = 0.0;
    for (int k = 0; k < 3; k++) {
      const double term = z[k] * s.coef[s.idx[k]][m];
      if (!(std::abs(term) < kIntegral)) {
        throw std::runtime_error(
            "minimum-image lattice coefficient exceeds binary64 integer precision");
      }
      sum += term;
    }
    n[m] += sum;
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
    if (!shorterLong(next, best)) {
      return best;
    }
    best = next;
    n[0] = m[0];
    n[1] = m[1];
    n[2] = m[2];
  }
}

// Wrap-path slicer. Uses 8-element arrays so the compiler vectorizes the
// initial projections, gain differences, and Gram updates, selects the best
// gain with a branchless 7-way reduction tree, and exits immediately when
// no projection exceeds its half-norm. Hitting `cap` falls back to Babai.
inline int descendFast(const Obtuse &s, Vec3 &x, int cap) {
  double t[8];
  for (int k = 0; k < 8; k++) {
    t[k] = x.x * s.rel[0][k] + x.y * s.rel[1][k] + x.z * s.rel[2][k];
  }
  const auto larger = [](int i, double gi, int j, double gj) {
    return gj > gi ? std::pair<int, double>{j, gj} : std::pair<int, double>{i, gi};
  };
  for (int iter = 0; iter < cap; iter++) {
    double g[8];
    for (int k = 0; k < 8; k++) {
      g[k] = std::abs(t[k]) - s.halfN2[k];
    }
    const auto p01 = larger(0, g[0], 1, g[1]);
    const auto p23 = larger(2, g[2], 3, g[3]);
    const auto p45 = larger(4, g[4], 5, g[5]);
    const auto p03 = larger(p01.first, p01.second, p23.first, p23.second);
    const auto p46 = larger(p45.first, p45.second, 6, g[6]);
    const auto best = larger(p03.first, p03.second, p46.first, p46.second);
    if (!(best.second > 0.0)) {
      return iter;
    }
    const int at = best.first;
    const double f = t[at] > 0.0 ? 1.0 : -1.0;
    x.x -= f * s.rel[0][at];
    x.y -= f * s.rel[1][at];
    x.z -= f * s.rel[2][at];
    for (int k = 0; k < 8; k++) {
      t[k] -= f * s.gram[at][k];
    }
    unsigned any = 0;
    for (int k = 0; k < 8; k++) {
      any |= static_cast<unsigned>(std::abs(t[k]) > s.halfN2[k]);
    }
    if (!any) {
      return iter + 1;
    }
  }
  return cap;
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
inline Vec3 closestBabai(const Obtuse &s, Vec3 w, double w2) {
  constexpr double kRange = std::numeric_limits<double>::max() / 16.0;
  const auto kept = [w2](double x2) { return x2 >= (1.0 / 256.0) * w2; };
  double n[3] = {0.0, 0.0, 0.0};
  double z[3];
  double scale = 0.0;
  const Vec3 rounded = babai(s, w, z, scale);
  if (scale <= 2.0 * maxAbs(w) && w2 < kRange) {
    const double r2 = n2(rounded);
    if (r2 <= w2) {
      Vec3 y = rounded;
      const int steps = descendFast(s, y, 3);
      if (steps < 3) {
        const double y2 = n2(y);
        if (kept(y2)) {
          return y2 < w2 ? y : w;
        }
      }
    }
  }
  addBabai(s, z, n);
  requireIntegral(n);
  Vec3 x = latticeTranslate(s, w, n);
  descend<false>(s, x, n, std::numeric_limits<int>::max());
  const Vec3 ref = refine(s, w, n);
  return n2(ref) < w2 ? ref : w;
}

// Binary64 results are kept only when nothing cancelled far below |w|: the
// Babai rounding subtracts lattice vectors with sum |z_i| max|b_i| at most
// 2 max|w|, a Babai descent takes at most two steps, a wrap descent takes
// at most three steps, and the image is at least |w| / 16. Every other
// image is re-evaluated from its lattice coefficients.
inline Vec3 closestDisplacement(const Obtuse &s, Vec3 w, double w2) {
  constexpr double kRange = std::numeric_limits<double>::max() / 16.0;
  if (!s.ok) {
    throw std::runtime_error(
        "minimum-image cell is numerically singular in binary64");
  }
  if (!(w2 * s.maxInvN2 < kIntegral * kIntegral)) {
    throw std::runtime_error(
        "minimum-image lattice coefficient exceeds binary64 integer precision");
  }
  if (s.fractional && w2 < kRange) {
    Vec3 y = w;
    const int steps = descendFast(s, y, kFractionalSteps);
    // Stopped inside the cap (at most three steps): keep y when it has not
    // cancelled below |w|/16. Heavy cancellation or hitting the cap discards
    // the partial move and starts Babai from the wrap.
    if (steps < kFractionalSteps) {
      const double y2 = n2(y);
      if (y2 >= (1.0 / 256.0) * w2) {
        return y2 < w2 ? y : w;
      }
    }
  }
  return closestBabai(s, w, w2);
}

inline Vec3 closestDisplacement(const Obtuse &s, Vec3 w) {
  return closestDisplacement(s, w, n2(w));
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
  (void)raw;
  const Vec3 w{wrapped[0], wrapped[1], wrapped[2]};
  const double f2 = n2(w);
  constexpr double kRange = std::numeric_limits<double>::max() / 16.0;
  if (!(f2 < kRange)) {
    for (double value : wrapped) {
      if (!std::isfinite(value)) {
        throw std::runtime_error("minimum-image candidate is not finite");
      }
    }
    throw std::runtime_error(
        "minimum-image squared norm exceeds binary64 range");
  }
  const Cached &slot = cachedBasis(lx, ly, lz, xy, xz, yz);
  const Vec3 euc = closestDisplacement(slot.s, w, f2);
  return {euc.x, euc.y, euc.z};
}

} // namespace mic
} // namespace gen

#endif
