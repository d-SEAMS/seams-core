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
// A slicer started from the wrap falls back to Babai after four steps. A
// cell-coordinate residual already inside the certified Smith ball (open
// ball of radius half the shortest relevant vector, shrunk by a rounding
// margin) is that ball's unique closest image and is accepted from an
// error-bounded evaluation, with no second descent and no certificate
// replay. The non-cancelling wrap descent does not reference that test or
// the Babai and recertify fallbacks; both sit in a separate call reached
// after the keep check. No
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
// Fractional residuals this close to an integer are the near-lattice case
// (a lattice vector plus a small Cartesian offset). Random fractional wraps
// fail the first axis. The length check, not this tolerance, is the certificate.
constexpr double kSmithTol = 0.05;

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

// v as rounded on its own. With FMA hardware the compiler may fuse the
// product that produced v into the sum that consumes it, and the error
// terms computed for the separately rounded product and sum would then no
// longer add up to the exact result.
inline double rounded(double v) {
#if (defined(__GNUC__) || defined(__clang__)) && defined(__SSE2__)
  __asm__("" : "+x"(v));
#elif (defined(__GNUC__) || defined(__clang__)) && defined(__aarch64__)
  __asm__("" : "+w"(v));
#elif defined(__GNUC__) || defined(__clang__)
  __asm__("" : "+m"(v));
#else
  volatile double kept = v;
  v = kept;
#endif
  return v;
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
    const double product = rounded(n[i] * x[i]);
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
  // Precomputed per-box Voronoi certificate bounds for the final slicer
  // projections t[0..6] after at most three wrap steps.
  double tCut[8]{};
  double tSlope[8]{};
  double maxRelLo1 = 0.0;
  double tRes = 0.0;
  // rel[m][k] + relLo[m][k] is the relevant vector to about eps^2.
  double relLo[3][8]{};
  // Voronoi test with rounding margins (see voronoiRoom).
  double cut[8]{};
  double slope[8]{};
  double invNorm[8]{};
  // Largest squared row of inv: |Babai coefficient|^2 <= n2(y) * maxInvN2.
  double maxInvN2 = 0.0;
  bool fractional = false;
  // Wrap images often cancel below |w| / 16, so the wrap descent records
  // its steps on the way down.
  bool recordWrap = false;
  // Compact fractional cells try the Smith-ball accept before descending.
  // Needle cells leave it on the cold side so the hot descent stays a leaf.
  bool smithHot = false;
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
  // Certified Smith ball. smithLim2 is a strict upper bound on the squared
  // length of a lattice image that is still inside the open ball of radius
  // λ1/2. inv* = 1/edge, for the fractional residual.
  double smithLim2 = 0.0;
  double smithZ0 = -1.0;
  double smithZ1 = std::numeric_limits<double>::infinity();
  double invX = 0.0;
  double invY = 0.0;
  double invZ = 0.0;
};

// a + b as sum + error, exactly (Knuth TwoSum).
inline double sumError(double a, double b, double sum) {
  const double back = sum - a;
  return (a - (sum - back)) + (b - back);
}

// w - n . (a, b, c), each component within about one rounding of itself, and
// in lo the rounding of that component, to about eps^2 of the summed terms.
// The cell is upper triangular, so x has three products, y two, and z one.
inline Vec3 latticeTranslate(const Obtuse &s, Vec3 w, const double n[3], Vec3 &lo) {
  const double m0 = -n[0];
  const double m1 = -n[1];
  const double m2 = -n[2];
  const Split s0 = split(m0);
  const Split s1 = split(m1);
  const Split s2 = split(m2);
  const double px0 = rounded(m0 * s.cellX[0]);
  const double px1 = rounded(m1 * s.cellX[1]);
  const double px2 = rounded(m2 * s.cellX[2]);
  const double py1 = rounded(m1 * s.cellY[1]);
  const double py2 = rounded(m2 * s.cellY[2]);
  const double pz2 = rounded(m2 * s.cellZ[2]);
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
  const double cx = ex + sumError(px0, px1, x01) + sumError(w.x, px2, x2w) + sumError(x01, x2w, sx);
  const double cy = ey + sumError(w.y, py1, sy) + sumError(sy, py2, ty);
  const double cz = ez + sumError(w.z, pz2, sz);
  const Vec3 hi = {sx + cx, ty + cy, sz + cz};
  lo = {sumError(sx, cx, hi.x), sumError(ty, cy, hi.y), sumError(sz, cz, hi.z)};
  return hi;
}

inline Vec3 latticeTranslate(const Obtuse &s, Vec3 w, const double n[3]) {
  Vec3 lo;
  return latticeTranslate(s, w, n, lo);
}

inline Vec3 latticeVector(const Obtuse &s, const double n[3], Vec3 &lo) {
  const double negated[3] = {-n[0], -n[1], -n[2]};
  return latticeTranslate(s, Vec3{}, negated, lo);
}

inline Vec3 latticeVector(const Obtuse &s, const double n[3]) {
  Vec3 lo;
  return latticeVector(s, n, lo);
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
      const long double z2b = static_cast<long double>(z2Best);
      return std::sqrt(z2b * z2b + (rem / q[2]) * (1.0L + 1e-9L));
    };
    long double w2 = z2Limit();
    if (!(std::isfinite(w2) && w2 < kWidth)) {
      return false;
    }
    // z2 is below 2^52, so its square is formed in long double, not long long.
    const auto square = [](long long z) {
      const long double l = static_cast<long double>(z);
      return l * l;
    };
    // Sign symmetry (-z0, -z1, -z2) preserves both parity and norm, so z2 >= 0 suffices.
    for (long long z2 = p2; static_cast<long double>(z2) <= w2; z2 += 2) {
      const long double rem = remNorm2();
      const long double budget1 =
          (z2 == z2Best ? 0.0L : -q[2] * (square(z2) - square(z2Best))) +
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
              (z2 == z2Best ? 0.0L : q[2] * (square(z2) - square(z2Best))) +
              q[1] * (e1 * e1 - e1Best * e1Best) +
              q[0] * (e0 * e0 - e0Best * e0Best);
          const long double diffScale =
              (z2 == z2Best ? 0.0L : q[2] * (square(z2) + square(z2Best))) +
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

    // Each product and the sum stay far inside long long before the exact
    // integer evaluation.
    constexpr long double kTerm = 0x1p62L;
    for (int m = 0; m < 3; m++) {
      const long double t0 = static_cast<long double>(z0Best) * U[0][m];
      const long double t1 = static_cast<long double>(z1Best) * U[1][m];
      const long double t2 = static_cast<long double>(z2Best) * U[2][m];
      if (!(std::abs(t0) < kTerm && std::abs(t1) < kTerm && std::abs(t2) < kTerm &&
            std::abs(t0 + t1 + t2) < kWidth)) {
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
    Vec3 lo{};
    if (k < 7) {
      for (int m = 0; m < 3; m++) {
        if (!(std::abs(coef[k][m]) < kIntegral)) {
          return false;
        }
      }
      r = latticeVector(o, coef[k], lo);
    }
    o.rel[0][k] = r.x;
    o.rel[1][k] = r.y;
    o.rel[2][k] = r.z;
    o.relLo[0][k] = lo.x;
    o.relLo[1][k] = lo.y;
    o.relLo[2][k] = lo.z;
    for (int m = 0; m < 3; m++) {
      o.relCoef[m][k] = (k < 7) ? coef[k][m] : 0.0;
    }
    o.halfN2[k] = (k < 7) ? 0.5 * n2(r) : std::numeric_limits<double>::infinity();
    o.cut[k] = std::numeric_limits<double>::infinity();
    o.slope[k] = 0.0;
    o.invNorm[k] = 1.0;
    if (k < 7) {
      if (!(std::isfinite(o.halfN2[k]) &&
            o.halfN2[k] >= std::numeric_limits<double>::min())) {
        return false;
      }
      o.maxRelN2 = std::max(o.maxRelN2, n2(r));
      constexpr double kMargin = 16.0 * std::numeric_limits<double>::epsilon();
      o.cut[k] = o.halfN2[k] - kMargin * n2(r) - 0x1p-1000;
      o.slope[k] = kMargin * std::sqrt(n2(r));
      o.invNorm[k] = 1.0 / std::sqrt(n2(r));
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
  // Precompute rigorous per-box margins so a short wrap descent (<= 3 steps)
  // can certify Voronoi interior directly from its final projections t[0..6].
  {
    constexpr long double kEpsD = std::numeric_limits<double>::epsilon();
    constexpr long double kEpsL = std::numeric_limits<long double>::epsilon();
    long double relTail[7]{};
    long double relErrInf[7]{};
    long double relL1[7]{};
    long double relInf[7]{};
    long double maxLo1 = 0.0L;
    long double maxTail = 0.0L;
    for (int k = 0; k < 7; k++) {
      const long double termSum =
          std::abs(static_cast<long double>(coef[k][0])) *
              std::abs(static_cast<long double>(o.cellX[0])) +
          std::abs(static_cast<long double>(coef[k][1])) *
              (std::abs(static_cast<long double>(o.cellX[1])) +
               std::abs(static_cast<long double>(o.cellY[1]))) +
          std::abs(static_cast<long double>(coef[k][2])) *
              (std::abs(static_cast<long double>(o.cellX[2])) +
               std::abs(static_cast<long double>(o.cellY[2])) +
               std::abs(static_cast<long double>(o.cellZ[2])));
      relTail[k] = 64.0L * kEpsD * kEpsD * termSum + 0x1p-1000L;
      maxTail = std::max(maxTail, relTail[k]);
      long double loInf = 0.0L;
      long double loL1 = 0.0L;
      long double rInf = 0.0L;
      long double l1 = 0.0L;
      for (int m = 0; m < 3; m++) {
        const long double rm = std::abs(static_cast<long double>(o.rel[m][k]));
        const long double lm = std::abs(static_cast<long double>(o.relLo[m][k]));
        rInf = std::max(rInf, rm);
        loInf = std::max(loInf, lm);
        loL1 += lm;
        l1 += rm + lm;
      }
      relInf[k] = rInf;
      relL1[k] = l1;
      relErrInf[k] = loInf + relTail[k];
      maxLo1 = std::max(maxLo1, loL1 + 3.0L * relTail[k]);
    }
    o.maxRelLo1 = std::isfinite(static_cast<double>(maxLo1))
                      ? static_cast<double>(maxLo1)
                      : std::numeric_limits<double>::infinity();
    const long double resScaled = 32.0L * maxLo1 + 0x1p55L * maxTail;
    o.tRes = std::isfinite(static_cast<double>(resScaled))
                 ? static_cast<double>(resScaled)
                 : std::numeric_limits<double>::infinity();
    for (int k = 0; k < 7; k++) {
      long double maxG = 0.0L;
      long double maxGErr = 0.0L;
      long double exactSelf = 0.0L;
      long double absSelf = 0.0L;
      for (int a = 0; a < 7; a++) {
        long double exactG = 0.0L;
        long double absG = 0.0L;
        for (int m = 0; m < 3; m++) {
          const long double ra =
              static_cast<long double>(o.rel[m][a]) + static_cast<long double>(o.relLo[m][a]);
          const long double rk =
              static_cast<long double>(o.rel[m][k]) + static_cast<long double>(o.relLo[m][k]);
          exactG += ra * rk;
          absG += std::abs(ra * rk);
        }
        if (a == k) {
          exactSelf = exactG;
          absSelf = absG;
        }
        const long double gErr =
            std::abs(static_cast<long double>(o.gram[a][k]) - exactG) +
            16.0L * kEpsL * absG + relL1[a] * relTail[k] + relL1[k] * relTail[a] +
            relTail[a] * relTail[k];
        maxGErr = std::max(maxGErr, gErr);
        maxG = std::max(maxG, std::abs(static_cast<long double>(o.gram[a][k])));
      }
      const long double halfExact =
          0.5L * exactSelf -
          (relL1[k] * relTail[k] + 0.5L * relTail[k] * relTail[k] + 8.0L * kEpsL * absSelf);
      const long double baseCut = std::min(static_cast<long double>(o.cut[k]), halfExact);
      const long double slope = 2.0L * (8.0L * kEpsD * relInf[k] + relErrInf[k]);
      const long double cst = 2.0L * (3.0L * maxGErr + 8.0L * kEpsD * maxG) + 0x1p-1000L;
      const double tc = static_cast<double>(baseCut - cst);
      const double ts = static_cast<double>(slope);
      if (std::isfinite(tc) && std::isfinite(ts)) {
        o.tCut[k] = tc;
        o.tSlope[k] = ts;
      } else {
        o.tCut[k] = -1.0;
        o.tSlope[k] = 0.0;
      }
    }
    o.tCut[7] = std::numeric_limits<double>::infinity();
    o.tSlope[7] = 0.0;
  }
  o.invX = 1.0 / o.cellX[0];
  o.invY = 1.0 / o.cellY[1];
  o.invZ = 1.0 / o.cellZ[2];
  if (!(std::isfinite(o.invX) && std::isfinite(o.invY) && std::isfinite(o.invZ))) {
    return true;
  }
  // λ1 is the shortest relevant vector: every nearest lattice vector is
  // Voronoi-relevant, and those are the seven coset minima stored in rel.
  // rel + relLo is that vector within about eps^2; the factors below sit
  // far above that error and far below a geometric margin.
  double minLower = std::numeric_limits<double>::infinity();
  for (int k = 0; k < 7; k++) {
    const double x = o.rel[0][k];
    const double y = o.rel[1][k];
    const double z = o.rel[2][k];
    const double v2 = x * x + y * y + z * z;
    if (!(v2 > 0.0) || !std::isfinite(v2)) {
      minLower = 0.0;
      break;
    }
    const double lo1 =
        std::abs(o.relLo[0][k]) + std::abs(o.relLo[1][k]) + std::abs(o.relLo[2][k]);
    const double scale = std::abs(x) + std::abs(y) + std::abs(z) + lo1;
    const double lower = std::sqrt(v2) * (1.0 - 0x1p-40) - lo1 - scale * 0x1p-80 - 0x1p-900;
    if (!(lower > 0.0)) {
      minLower = 0.0;
      break;
    }
    minLower = std::min(minLower, lower);
  }
  if (minLower > 0.0 && std::isfinite(minLower)) {
    const double radius = 0.5 * minLower;
    o.smithLim2 = radius * radius * (1.0 - 0x1p-40);
    const double lzAbs = std::abs(o.cellZ[2]);
    o.smithZ0 = kSmithTol * lzAbs;
    o.smithZ1 = (1.0 - kSmithTol) * lzAbs;
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
// A binary64 descent from the Babai point that has not settled after this
// many steps continues in extended precision.
constexpr int kSlicerSteps = 16;

template <bool Track = false, bool Record = false, bool OutT = false>
inline int descendFast(const Obtuse &s, Vec3 &x, int cap, double *n = nullptr,
                       unsigned long long *path = nullptr, double *tOut = nullptr);

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
    // A fractional wrap is at most half the summed edge lengths long. When
    // the shortest superbasis vector is below 1/16 of that, wrap images
    // cancel below |w| / 16 often enough that the descent records its steps
    // on the way down rather than running again.
    double shortest = n2(v[0]);
    for (int i = 1; i < 4; i++) {
      shortest = std::min(shortest, n2(v[i]));
    }
    const double la = std::sqrt(n2(a));
    const double lb = std::sqrt(n2(b));
    const double lc = std::sqrt(n2(c));
    const double edges = la + lb + lc;
    // Only a fractional cell records. A cell that already failed the wrap
    // probes (skew-both) stays on Babai; recordWrap must not steal that path.
    o.recordWrap = o.fractional && 1024.0 * shortest < edges * edges;
    const double mx = std::max(la, std::max(lb, lc));
    const double mn = std::min(la, std::min(lb, lc));
    // Aspect under 10: the inline residual test is paid back by the compact
    // near-lattice pool. Needle cells (aspect-5 tilt, long skew) skip it.
    // Aspect under 40. The residual test rejects a random wrap on the first
    // axis, so the Smith call is not on that path, and the wider band covers
    // the moderate skew cells whose near-lattice images otherwise descend first.
    o.smithHot = o.fractional && !o.recordWrap && o.smithLim2 > 0.0 && mn > 0.0 &&
                 mx < 40.0 * mn;
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
// With Track, the integer cell coordinates of each step are added to n.
// With Record (cap at most 16), path gets 4 bits per step, the first step
// highest: 2 * k + (sign < 0) for the relevant vector k taken.
template <bool Track, bool Record, bool OutT>
#if defined(__GNUC__) || defined(__clang__)
__attribute__((optimize("O3,tree-vectorize")))
#endif
inline int descendFast(const Obtuse &s, Vec3 &x, int cap, double *n, unsigned long long *path,
                       double *tOut) {
  unsigned long long code = 0;
  double xx = x.x;
  double xy = x.y;
  double xz = x.z;
  double t[8];
  for (int k = 0; k < 8; k++) {
    t[k] = xx * s.rel[0][k] + xy * s.rel[1][k] + xz * s.rel[2][k];
  }
  for (int iter = 0; iter < cap; iter++) {
    double g[8];
    for (int k = 0; k < 8; k++) {
      g[k] = std::abs(t[k]) - s.halfN2[k];
    }
    const bool c01 = g[1] > g[0];
    const double m01 = c01 ? g[1] : g[0];
    const int i01 = c01 ? 1 : 0;
    const bool c23 = g[3] > g[2];
    const double m23 = c23 ? g[3] : g[2];
    const int i23 = c23 ? 3 : 2;
    const bool c45 = g[5] > g[4];
    const double m45 = c45 ? g[5] : g[4];
    const int i45 = c45 ? 5 : 4;
    const bool c03 = m23 > m01;
    const double m03 = c03 ? m23 : m01;
    const int i03 = c03 ? i23 : i01;
    const bool c46 = g[6] > m45;
    const double m46 = c46 ? g[6] : m45;
    const int i46 = c46 ? 6 : i45;
    const bool c06 = m46 > m03;
    const double bestG = c06 ? m46 : m03;
    const int at = c06 ? i46 : i03;
    if (!(bestG > 0.0)) {
      x = {xx, xy, xz};
      if constexpr (Record) {
        *path = code;
      }
      if constexpr (OutT) {
        for (int k = 0; k < 8; k++) {
          tOut[k] = t[k];
        }
      }
      return iter;
    }
    const bool pos = t[at] > 0.0;
    const double f = pos ? 1.0 : -1.0;
    if constexpr (Record) {
      code = (code << 4) | static_cast<unsigned long long>(2 * at + (pos ? 0 : 1));
    }
    xx -= f * s.rel[0][at];
    xy -= f * s.rel[1][at];
    xz -= f * s.rel[2][at];
    if constexpr (Track) {
      for (int m = 0; m < 3; m++) {
        n[m] += f * s.relCoef[m][at];
      }
    }
    if (iter + 1 < cap) {
      for (int k = 0; k < 8; k++) {
        t[k] -= f * s.gram[at][k];
      }
      unsigned any = 0;
      for (int k = 0; k < 8; k++) {
        any |= static_cast<unsigned>(std::abs(t[k]) > s.halfN2[k]);
      }
      if (!any) {
        x = {xx, xy, xz};
        if constexpr (Record) {
          *path = code;
        }
        if constexpr (OutT) {
          for (int k = 0; k < 8; k++) {
            tOut[k] = t[k];
          }
        }
        return iter + 1;
      }
    }
  }
  x = {xx, xy, xz};
  if constexpr (Record) {
    *path = code;
  }
  return cap;
}

// Distance by which x clears the Voronoi cell of the origin with room for
// rounding, positive when it lies strictly inside: x is within about one
// rounding of its lattice image, each stored relevant vector within about
// one rounding of its lattice vector, and the binary64 dot product within
// 3 eps |x| |r|. Those move x . r by under 8 eps |x| |r| and |r|^2 / 2 by
// under 4 eps |r|^2, half the margins in cut and slope, so the exact image
// is then the unique minimum image, and so is any lattice image within the
// returned distance of x. |x|_1 bounds |x|. The 2^-1000 in cut covers
// absolute rounding where products are subnormal.
inline double voronoiRoom(const Obtuse &s, Vec3 x) {
  const double xn = std::abs(x.x) + std::abs(x.y) + std::abs(x.z);
  double room[8];
  for (int k = 0; k < 8; k++) {
    const double p = x.x * s.rel[0][k] + x.y * s.rel[1][k] + x.z * s.rel[2][k];
    room[k] = (s.cut[k] - (std::abs(p) + s.slope[k] * xn)) * s.invNorm[k];
  }
  return std::min(std::min(std::min(room[0], room[1]), std::min(room[2], room[3])),
                  std::min(std::min(room[4], room[5]), room[6]));
}

inline bool insideVoronoi(const Obtuse &s, Vec3 x) { return voronoiRoom(s, x) > 0.0; }

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
inline Vec3 closestBabai(const Obtuse &s, const Vec3 &w, double w2) {
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
  Vec3 y = x;
  const int steps = descendFast<true>(s, y, kSlicerSteps, n);
  if (steps < kSlicerSteps) {
    requireIntegral(n);
    if (steps > 0) {
      x = latticeTranslate(s, w, n);
    }
    if (insideVoronoi(s, x)) {
      return n2(x) < w2 ? x : w;
    }
  }
  const Vec3 ref = refine(s, w, n);
  return n2(ref) < w2 ? ref : w;
}

// Replay helper for the rare near-tie fallback: recover integer cell
// coordinates n from the recorded path and run extended-precision refine.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
inline Vec3 refineFromPath(const Obtuse &s, Vec3 w, double w2, unsigned long long path,
                           int steps) {
  double n[3] = {0.0, 0.0, 0.0};
  for (int i = steps - 1; i >= 0; i--) {
    const unsigned code = static_cast<unsigned>(path >> (4 * i)) & 15u;
    const int at = static_cast<int>(code >> 1);
    const double f = (code & 1u) != 0 ? -1.0 : 1.0;
    for (int m = 0; m < 3; m++) {
      n[m] += f * s.relCoef[m][at];
    }
  }
  const Vec3 ref = refine(s, w, n);
  return n2(ref) < w2 ? ref : w;
}

// Slow fallback when the fast tFinal certificate does not clear its margin:
// check voronoiRoom(y) > 2 |lo|_1 first, then recover integer coordinates n
// from the recorded path and run extended-precision refine.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
inline Vec3 certifySlow(const Obtuse &s, const Vec3 &w, double w2, const Vec3 &y,
                        const Vec3 &x, const Vec3 &lo, double moved,
                        unsigned long long path, int steps) {
  if (x.x == y.x && x.y == y.y && x.z == y.z && voronoiRoom(s, y) > 2.0 * moved) {
    const Vec3 image = add(y, lo);
    return n2(image) < w2 ? image : w;
  }
  return refineFromPath(s, w, w2, path, steps);
}

// A short wrap descent from w whose image y cancelled below |w| / 16, with
// its `steps` steps recorded in `path` and final projections in `tFinal`.
// The steps subtract relevant vectors with sign +-1, so replaying them
// repeats every rounding exactly; the rounding of each subtraction and the
// low part of each relevant vector give the image as y + lo to about eps^2
// per step. When tFinal clears the precomputed per-box Voronoi margins,
// y + lo is the unique minimum image without calling voronoiRoom or refine.
// For a single step with exact/tiny relLo, IEEE-754 round-to-nearest already
// guarantees |y - y_exact| <= 1.5 eps |y|_1, avoiding the TwoSum replay too.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, optimize("O3,tree-vectorize")))
#endif
inline Vec3 certifyWrapT(const Obtuse &s, double wx, double wy, double wz, double w2,
                         double yx, double yy, double yz, unsigned long long path, int steps,
                         const double tFinal[8]) {
  const double wn = std::abs(wx) + std::abs(wy) + std::abs(wz);
  unsigned tOk = 1u;
  for (int k = 0; k < 8; k++) {
    tOk &= static_cast<unsigned>(std::abs(tFinal[k]) + s.tSlope[k] * wn < s.tCut[k]);
  }
  const Vec3 y{yx, yy, yz};
  const Vec3 w{wx, wy, wz};
  if (steps == 1 && tOk != 0u) {
    const double y1 = std::abs(yx) + std::abs(yy) + std::abs(yz);
    if (y1 > 0x1p-900 && s.maxRelLo1 <= 0x1p-52 * y1) {
      return n2(y) < w2 ? y : w;
    }
  }
  static constexpr double kSign[2] = {1.0, -1.0};
  double xx = wx;
  double xy = wy;
  double xz = wz;
  double lox = 0.0;
  double loy = 0.0;
  double loz = 0.0;
  const auto replayStep = [&](unsigned code) {
    const int at = static_cast<int>(code >> 1);
    const double f = kSign[code & 1u];
    const double rx = -f * s.rel[0][at];
    const double ry = -f * s.rel[1][at];
    const double rz = -f * s.rel[2][at];
    const double nx = xx + rx;
    const double ny = xy + ry;
    const double nz = xz + rz;
    lox += sumError(xx, rx, nx) - f * s.relLo[0][at];
    loy += sumError(xy, ry, ny) - f * s.relLo[1][at];
    loz += sumError(xz, rz, nz) - f * s.relLo[2][at];
    xx = nx;
    xy = ny;
    xz = nz;
  };
  if (steps == 2) {
    replayStep(static_cast<unsigned>(path >> 4) & 15u);
    replayStep(static_cast<unsigned>(path) & 15u);
  } else if (steps == 3) {
    replayStep(static_cast<unsigned>(path >> 8) & 15u);
    replayStep(static_cast<unsigned>(path >> 4) & 15u);
    replayStep(static_cast<unsigned>(path) & 15u);
  } else {
    for (int i = steps - 1; i >= 0; i--) {
      replayStep(static_cast<unsigned>(path >> (4 * i)) & 15u);
    }
  }
  const double moved = std::abs(lox) + std::abs(loy) + std::abs(loz);
  const Vec3 image{yx + lox, yy + loy, yz + loz};
  const double img1 = std::abs(image.x) + std::abs(image.y) + std::abs(image.z);
#if defined(__GNUC__) || defined(__clang__)
  if (__builtin_expect(xx == yx && xy == yy && xz == yz && tOk != 0u &&
                           moved <= 0x1p-20 * img1 && wn * 0x1p-47 + s.tRes <= img1,
                       1)) {
#else
  if (xx == yx && xy == yy && xz == yz && tOk != 0u && moved <= 0x1p-20 * img1 &&
      wn * 0x1p-47 + s.tRes <= img1) {
#endif
    return n2(image) < w2 ? image : w;
  }
  return certifySlow(s, w, w2, y, Vec3{xx, xy, xz}, Vec3{lox, loy, loz}, moved, path, steps);
}

// certifyWrapT for a cell that does not record its wrap steps: the descent
// from w runs again, recording.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, optimize("O3,tree-vectorize")))
#endif
inline Vec3 recertifyWrap(const Obtuse &s, const Vec3 &w, double w2) {
  Vec3 y = w;
  unsigned long long path = 0;
  double tFinal[8];
  const int steps = descendFast<false, true, true>(s, y, kFractionalSteps, nullptr, &path, tFinal);
  if (steps < kFractionalSteps) {
    return certifyWrapT(s, w.x, w.y, w.z, w2, y.x, y.y, y.z, path, steps, tFinal);
  }
  return closestBabai(s, w, w2);
}

// Nearest integer, half away from zero, or false when x is not a safe
// binary64 integer coefficient.
inline bool nearestInt(double x, double &n) {
  if (!(std::abs(x) < kIntegral)) {
    return false;
  }
  const double t = x + std::copysign(0.49999999999999994, x);
  n = static_cast<double>(static_cast<long long>(t));
  return true;
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline))
#endif
inline bool smithCandidate(const Obtuse &s, Vec3 w, double &nx, double &ny, double &nz) {
  const double z = w.z * s.invZ;
  if (!nearestInt(z, nz) || !(std::abs(z - nz) < kSmithTol)) {
    return false;
  }
  const double y = (w.y - z * s.cellY[2]) * s.invY;
  if (!nearestInt(y, ny) || !(std::abs(y - ny) < kSmithTol)) {
    return false;
  }
  const double x = (w.x - y * s.cellX[1] - z * s.cellX[2]) * s.invX;
  return nearestInt(x, nx) && std::abs(x - nx) < kSmithTol;
}

// Image of w - n.(a,b,c) for coefficients in {-1,0,1}. Those products are
// exact. TwoSum on each sum leaves a residual of order eps^2 |w|, returned
// in radBound as a strict upper bound on the absolute error of every
// component (far above the residual, far below a geometric margin).
inline bool unitCoeffImage(const Obtuse &s, Vec3 w, double nx, double ny, double nz, Vec3 &out,
                           double &radBound) {
  const double p0 = rounded(nx * s.cellX[0]);
  const double p1 = rounded(ny * s.cellX[1]);
  const double p2 = rounded(nz * s.cellX[2]);
  const double q1 = rounded(ny * s.cellY[1]);
  const double q2 = rounded(nz * s.cellY[2]);
  const double r2 = rounded(nz * s.cellZ[2]);
  const double s01 = p0 + p1;
  const double e01 = sumError(p0, p1, s01);
  const double px = s01 + p2;
  const double epx = sumError(s01, p2, px);
  const double ix = w.x - px;
  const double eix = sumError(w.x, -px, ix);
  const double xx = rounded(ix + ((eix - epx) - e01));
  const double py = q1 + q2;
  const double epy = sumError(q1, q2, py);
  const double iy = w.y - py;
  const double eiy = sumError(w.y, -py, iy);
  const double yy = rounded(iy + (eiy - epy));
  const double iz = w.z - r2;
  const double eiz = sumError(w.z, -r2, iz);
  const double zz = rounded(iz + eiz);
  if (!(std::isfinite(xx) && std::isfinite(yy) && std::isfinite(zz))) {
    return false;
  }
  out = {xx, yy, zz};
  const double scale = std::abs(w.x) + std::abs(w.y) + std::abs(w.z) + std::abs(p0) +
                       std::abs(p1) + std::abs(p2) + std::abs(q1) + std::abs(q2) + std::abs(r2);
  radBound = scale * 0x1p-90 + 0x1p-900;
  return true;
}

// Accept the image of the coefficients from smithCandidate when a rigorous
// upper bound on its length sits under smithLim2. The ball of radius λ1/2
// is inside the Voronoi cell, so that image is the unique minimum image.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline))
#endif
inline bool smithFromN(const Obtuse &s, Vec3 w, double w2, double nx, double ny, double nz,
                       Vec3 &out) {
  if (!(s.smithLim2 > 0.0)) {
    return false;
  }
  Vec3 image;
  double radBound = 0.0;
  if (std::abs(nx) <= 1.0 && std::abs(ny) <= 1.0 && std::abs(nz) <= 1.0) {
    if (!unitCoeffImage(s, w, nx, ny, nz, image, radBound)) {
      return false;
    }
  } else {
    const double n[3] = {nx, ny, nz};
    Vec3 lo;
    const Vec3 hi = latticeTranslate(s, w, n, lo);
    image = {hi.x + lo.x, hi.y + lo.y, hi.z + lo.z};
    if (!(std::isfinite(image.x) && std::isfinite(image.y) && std::isfinite(image.z))) {
      return false;
    }
    const double scale = std::abs(w.x) + std::abs(w.y) + std::abs(w.z) +
                         std::abs(nx) * (std::abs(s.cellX[0]) + std::abs(s.cellX[1]) +
                                         std::abs(s.cellX[2])) +
                         std::abs(ny) * (std::abs(s.cellY[1]) + std::abs(s.cellY[2])) +
                         std::abs(nz) * std::abs(s.cellZ[2]);
    radBound = scale * 0x1p-80 + 0x1p-900;
  }
  const double bx = std::abs(image.x) + radBound;
  const double by = std::abs(image.y) + radBound;
  const double bz = std::abs(image.z) + radBound;
  const double upper = (bx * bx + by * by) + bz * bz;
  if (!(upper < s.smithLim2 && upper * (1.0 + 0x1p-40) < s.smithLim2)) {
    return false;
  }
  if (!(n2(image) < w2)) {
    return false;
  }
  out = image;
  return true;
}

// Non-cancelling plain wrap. No recertify or Babai symbol.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline))
#endif
inline bool hotPlain(const Obtuse &s, Vec3 w, double w2, Vec3 &out, int &steps) {
  Vec3 y = w;
  steps = descendFast(s, y, kFractionalSteps);
  out = y;
  if (steps < kFractionalSteps) {
    const double y2 = n2(y);
    if (y2 >= (1.0 / 256.0) * w2) {
      if (!(y2 < w2)) {
        out = w;
      }
      return true;
    }
  }
  return false;
}

// Speed-optimized, not cold: long-skew cells take the failure edge on a
// sizable fraction of wraps, and a size-optimized clone shows up in the median.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
inline Vec3 coldResolve(const Obtuse &s, const Vec3 &w, double w2, const Vec3 &partial, int steps) {
  constexpr double kRange = std::numeric_limits<double>::max() / 16.0;
  if (s.smithLim2 > 0.0 && w2 > s.smithLim2 && w2 < kRange) {
    double nx = 0.0;
    double ny = 0.0;
    double nz = 0.0;
    Vec3 image;
    if (smithCandidate(s, w, nx, ny, nz) && smithFromN(s, w, w2, nx, ny, nz, image)) {
      return image;
    }
  }
  // The hot budget is four steps. Every measured wrap that hits it is the
  // minimum one or two slicer steps later, so continue from that image
  // instead of starting Babai. A cancelled image is not kept.
  if (s.fractional && w2 < kRange && steps >= kFractionalSteps) {
    Vec3 y = partial;
    const int st = descendFast(s, y, 2);
    if (st < 2) {
      const double y2 = n2(y);
      if (y2 >= (1.0 / 256.0) * w2) {
        return y2 < w2 ? y : w;
      }
    }
    return closestBabai(s, w, w2);
  }
  if (s.fractional && w2 < kRange) {
    return recertifyWrap(s, w, w2);
  }
  return closestBabai(s, w, w2);
}

// Smith accept as its own speed-optimized leaf, so the residual test and the
// compensated image are not in the straight-line wrap descent.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
inline bool trySmith(const Obtuse &s, const Vec3 &w, double w2, Vec3 &out) {
  double nx = 0.0;
  double ny = 0.0;
  double nz = 0.0;
  if (!smithCandidate(s, w, nx, ny, nz)) {
    return false;
  }
  return smithFromN(s, w, w2, nx, ny, nz, out);
}

// Recording slicer. Runs at O3+vectorize with the descent loop directly
// inlined so there is no call overhead or stack copy of t[8] on kept wraps.
// Near-lattice residuals inside the Smith ball are filtered in 3 instructions
// on |w.z| before descending; cancelled descents certify from t[0..6].
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, optimize("O3,tree-vectorize")))
#endif
inline Vec3 recordDisplacement(const Obtuse &s, double wx, double wy, double wz, double w2) {
  const double az = std::abs(wz);
#if defined(__GNUC__) || defined(__clang__)
  if (__builtin_expect(az < s.smithZ0 || az > s.smithZ1, 0)) {
#else
  if (az < s.smithZ0 || az > s.smithZ1) {
#endif
    const double z = wz * s.invZ;
    double nz = 0.0;
    if (nearestInt(z, nz) && std::abs(z - nz) < kSmithTol) {
      const double y = (wy - z * s.cellY[2]) * s.invY;
      double ny = 0.0;
      if (nearestInt(y, ny) && std::abs(y - ny) < kSmithTol) {
        Vec3 image;
        if (trySmith(s, Vec3{wx, wy, wz}, w2, image)) {
          return image;
        }
      }
    }
  }
  unsigned long long code = 0;
  double xx = wx;
  double xy = wy;
  double xz = wz;
  double t[8];
  for (int k = 0; k < 8; k++) {
    t[k] = xx * s.rel[0][k] + xy * s.rel[1][k] + xz * s.rel[2][k];
  }
  int steps = kFractionalSteps;
  for (int iter = 0; iter < kFractionalSteps - 1; iter++) {
    double g[8];
    for (int k = 0; k < 8; k++) {
      g[k] = std::abs(t[k]) - s.halfN2[k];
    }
    const bool c01 = g[1] > g[0];
    const double m01 = c01 ? g[1] : g[0];
    const int i01 = c01 ? 1 : 0;
    const bool c23 = g[3] > g[2];
    const double m23 = c23 ? g[3] : g[2];
    const int i23 = c23 ? 3 : 2;
    const bool c45 = g[5] > g[4];
    const double m45 = c45 ? g[5] : g[4];
    const int i45 = c45 ? 5 : 4;
    const bool c03 = m23 > m01;
    const double m03 = c03 ? m23 : m01;
    const int i03 = c03 ? i23 : i01;
    const bool c46 = g[6] > m45;
    const double m46 = c46 ? g[6] : m45;
    const int i46 = c46 ? 6 : i45;
    const bool c06 = m46 > m03;
    const double bestG = c06 ? m46 : m03;
    const int at = c06 ? i46 : i03;
    if (!(bestG > 0.0)) {
      steps = iter;
      break;
    }
    const bool pos = t[at] > 0.0;
    const double f = pos ? 1.0 : -1.0;
    code = (code << 4) | static_cast<unsigned long long>(2 * at + (pos ? 0 : 1));
    xx -= f * s.rel[0][at];
    xy -= f * s.rel[1][at];
    xz -= f * s.rel[2][at];
    for (int k = 0; k < 8; k++) {
      t[k] -= f * s.gram[at][k];
    }
  }
  if (steps == kFractionalSteps) {
    unsigned any = 0;
    for (int k = 0; k < 8; k++) {
      any |= static_cast<unsigned>(std::abs(t[k]) > s.halfN2[k]);
    }
    if (!any) {
      steps = kFractionalSteps - 1;
    }
  }
  if (steps < kFractionalSteps) {
    const double y2 = xx * xx + xy * xy + xz * xz;
#if defined(__GNUC__) || defined(__clang__)
    if (__builtin_expect(y2 >= (1.0 / 256.0) * w2, 1)) {
#else
    if (y2 >= (1.0 / 256.0) * w2) {
#endif
      return y2 < w2 ? Vec3{xx, xy, xz} : Vec3{wx, wy, wz};
    }
    return certifyWrapT(s, wx, wy, wz, w2, xx, xy, xz, code, steps, t);
  }
  return closestBabai(s, Vec3{wx, wy, wz}, w2);
}

// Non-cancelling wrap. Names neither recertify nor Babai; the failure
// edge is one separate call.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline))
#endif
inline Vec3 plainDisplacement(const Obtuse &s, const Vec3 &w, double w2) {
  Vec3 out;
  int steps = -1;
#if defined(__GNUC__) || defined(__clang__)
  if (__builtin_expect(hotPlain(s, w, w2, out, steps), 1)) {
#else
  if (hotPlain(s, w, w2, out, steps)) {
#endif
    return out;
  }
  return coldResolve(s, w, w2, out, steps);
}

// Binary64 results are kept only when nothing cancelled far below |w|: the
// Babai rounding subtracts lattice vectors with sum |z_i| max|b_i| at most
// 2 max|w|, a Babai descent takes at most two steps, a wrap descent takes
// at most three steps, and the image is at least |w| / 16. A short wrap
// that cancelled is certified from its recorded steps and final projections,
// or recorded on a second descent when the cell does not record. A wrap that
// uses the whole four-step budget is continued for two slicer steps and kept
// when that image has not cancelled; otherwise it goes to Babai. Any vector
// already inside the certified Smith ball (|w|^2 < smithLim2) is its own
// unique minimum image and returns immediately.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline))
#endif
inline Vec3 closestDisplacement(const Obtuse &s, const Vec3 &w, double w2) {
  if (!s.ok) {
    throw std::runtime_error(
        "minimum-image cell is numerically singular in binary64");
  }
  if (!(w2 * s.maxInvN2 < kIntegral * kIntegral)) {
    throw std::runtime_error(
        "minimum-image lattice coefficient exceeds binary64 integer precision");
  }
  if (w2 < s.smithLim2) {
    return w;
  }
  constexpr double kRange = std::numeric_limits<double>::max() / 16.0;
#if defined(__GNUC__) || defined(__clang__)
  if (__builtin_expect(s.recordWrap && w2 < kRange, 0)) {
#else
  if (s.recordWrap && w2 < kRange) {
#endif
    return recordDisplacement(s, w.x, w.y, w.z, w2);
  }
#if defined(__GNUC__) || defined(__clang__)
  if (__builtin_expect(!(s.fractional && w2 < kRange), 0)) {
#else
  if (!(s.fractional && w2 < kRange)) {
#endif
    return closestBabai(s, w, w2);
  }
  // Compact cells only. The 3-instruction |w.z| filter rejects 90% of random
  // wraps without division or integer conversion; a near-lattice residual
  // falls through into trySmith.
  if (s.smithHot) {
    const double az = std::abs(w.z);
#if defined(__GNUC__) || defined(__clang__)
    if (__builtin_expect(az < s.smithZ0 || az > s.smithZ1, 0)) {
#else
    if (az < s.smithZ0 || az > s.smithZ1) {
#endif
      Vec3 image;
      if (trySmith(s, w, w2, image)) {
        return image;
      }
    }
  }
  return plainDisplacement(s, w, w2);
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

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
inline const Cached &fillCachedBasis(Cached (&slots)[4], int &hand, const Cached *&lastSlot,
                                     double lx, double ly, double lz,
                                     double xy, double xz, double yz) {
  const double key[6] = {lx, ly, lz, xy, xz, yz};
  for (int idx = 0; idx < 4; idx++) {
    if (slots[idx].filled && sameKey(slots[idx].key, key)) {
      lastSlot = &slots[idx];
      return slots[idx];
    }
  }
  const int i = hand;
  hand = (hand + 1) % 4;
  Cached &slot = slots[i];
  const Vec3 a{lx, 0.0, 0.0};
  const Vec3 b{xy, ly, 0.0};
  const Vec3 c{xz, yz, lz};
  if (lastSlot == &slot) {
    lastSlot = nullptr;
  }
  slot.filled = false;
  const Obtuse s = obtuseSuperbasis(a, b, c);
  for (int k = 0; k < 6; k++) {
    slot.key[k] = key[k];
  }
  slot.s = s;
  slot.filled = true;
  lastSlot = &slot;
  return slot;
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline))
#endif
inline const Cached &cachedBasis(double lx, double ly, double lz, double xy, double xz,
                                 double yz) {
  thread_local Cached slots[4];
  thread_local int hand = 0;
  thread_local const Cached *lastSlot = nullptr;
  const Cached *mru = lastSlot;
  if (mru != nullptr && mru->key[0] == lx && mru->key[1] == ly && mru->key[2] == lz &&
      mru->key[3] == xy && mru->key[4] == xz && mru->key[5] == yz) {
    return *mru;
  }
  return fillCachedBasis(slots, hand, lastSlot, lx, ly, lz, xy, xz, yz);
}

} // namespace

// Cartesian displacement from (xj, yj, zj) to (xi, yi, zi) on the
// recovered dump cell. `wrapped` is the fractional image of that
// difference. The fractional vector is kept when it is not strictly
// longer than the Selling image.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline))
#endif
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
