// Copyright (c) d-SEAMS developers
//
// Binary128 reference for the Euclidean minimum image on a triclinic cell
// a = (lx, 0, 0), b = (xy, ly, 0), c = (xz, yz, lz). It shares no code with
// src/include/internal/euclidean_mic.hpp: the cell is LLL-reduced over
// integer coefficients in __float128, each reduced vector is evaluated from
// its coefficients, and a Fincke-Pohst enumeration over the reduced basis
// visits every lattice translate of the displacement inside the Babai bound,
// so the shortest image is exact to binary128 rounding.
//
// A size-reduction step whose integer would not fit __int128 (a needle such
// as lx 1e160 against unit edges) updates that vector in binary128 and marks
// it inexact. An image that uses an inexact vector is refused (ok false);
// for a displacement much shorter than that vector the enumeration never
// uses it, and the image stays exact.
#pragma once

#include <quadmath.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace micoracle {

using Q = __float128;
using I = __int128;

struct Cell {
  double lx, ly, lz, xy, xz, yz;
};

struct V {
  Q x = 0, y = 0, z = 0;
};

inline V operator-(V a, V b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V operator+(V a, V b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline V operator*(Q s, V a) { return {s * a.x, s * a.y, s * a.z}; }
inline Q dot(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Q qabs(Q v) { return v < 0 ? -v : v; }
inline Q qround(Q v) { return roundq(v); }

// n . (a, b, c).
inline V combo(const Cell &c, const I n[3]) {
  const Q n0 = (Q)n[0], n1 = (Q)n[1], n2 = (Q)n[2];
  return {n0 * (Q)c.lx + n1 * (Q)c.xy + n2 * (Q)c.xz, n1 * (Q)c.ly + n2 * (Q)c.yz,
          n2 * (Q)c.lz};
}

struct Basis {
  I coef[3][3]{};  // row i: reduced vector i in (a, b, c), when exact[i]
  bool exact[3] = {true, true, true};
  V v[3];
  V gs[3];
  Q gn[3]{};
  Q mu[3][3]{};
  bool ok = false;
};

inline void gramSchmidt(Basis &r) {
  for (int i = 0; i < 3; i++) {
    V g = r.v[i];
    for (int j = 0; j < i; j++) {
      r.mu[i][j] = dot(r.v[i], r.gs[j]) / r.gn[j];
      g = g - r.mu[i][j] * r.gs[j];
    }
    r.gs[i] = g;
    r.gn[i] = dot(g, g);
    r.mu[i][i] = 1;
  }
}

// LLL with delta 0.99, size-reducing against every earlier vector.
inline Basis reduce(const Cell &cell) {
  Basis r;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      r.coef[i][j] = i == j;
    }
    r.v[i] = combo(cell, r.coef[i]);
  }
  gramSchmidt(r);
  int k = 1;
  for (int guard = 0; k < 3; guard++) {
    if (guard > 100000 || !(r.gn[0] > 0 && r.gn[1] > 0 && r.gn[2] > 0)) {
      return r;
    }
    for (int j = k - 1; j >= 0; j--) {
      const Q q = qround(r.mu[k][j]);
      if (q != 0) {
        const Q small = (Q)0x1p60;
        bool fits = r.exact[k] && r.exact[j] && qabs(q) < small;
        for (int t = 0; t < 3 && fits; t++) {
          fits = qabs((Q)r.coef[k][t]) < small && qabs((Q)r.coef[j][t]) < small;
        }
        if (fits) {
          const I qi = (I)q;
          for (int t = 0; t < 3; t++) {
            r.coef[k][t] -= qi * r.coef[j][t];
          }
          r.v[k] = combo(cell, r.coef[k]);
        } else {
          r.v[k] = r.v[k] - q * r.v[j];
          r.exact[k] = false;
        }
        gramSchmidt(r);
      }
    }
    if (r.gn[k] >= ((Q)0.99 - r.mu[k][k - 1] * r.mu[k][k - 1]) * r.gn[k - 1]) {
      k++;
    } else {
      for (int t = 0; t < 3; t++) {
        std::swap(r.coef[k][t], r.coef[k - 1][t]);
      }
      std::swap(r.exact[k], r.exact[k - 1]);
      std::swap(r.v[k], r.v[k - 1]);
      gramSchmidt(r);
      k = std::max(k - 1, 1);
    }
  }
  r.ok = r.gn[0] > 0 && r.gn[1] > 0 && r.gn[2] > 0;
  return r;
}

struct Image {
  V d;           // w - n . (a, b, c)
  Q n2 = 0;      // |d|^2
  I n[3]{};
  Q x[3]{};      // coefficients over the reduced basis
  bool ok = false;
};

// Calls visit(x0, x1, x2, |w - x . v|^2) for every translate of w by the
// reduced basis with squared length at most bound (Fincke-Pohst).
template <class Visit>
inline void enumerate(const Basis &r, V w, Q bound, Visit &&visit) {
  Q t[3];
  for (int k = 0; k < 3; k++) {
    t[k] = dot(w, r.gs[k]) / r.gn[k];
  }
  Q cur[3] = {0, 0, 0};
  const auto centre = [&](int k) {
    Q c = t[k];
    for (int i = k + 1; i < 3; i++) {
      c -= cur[i] * r.mu[i][k];
    }
    return c;
  };
  const auto span = [&](int k, Q partial) {
    const Q room = (bound - partial) / r.gn[k];
    return room > 0 ? sqrtq(room) : (Q)0;
  };
  const Q c2 = centre(2);
  const Q s2 = span(2, 0);
  for (Q x2 = ceilq(c2 - s2); x2 <= floorq(c2 + s2); x2 += 1) {
    cur[2] = x2;
    const Q p2 = r.gn[2] * (c2 - x2) * (c2 - x2);
    const Q c1 = centre(1);
    const Q s1 = span(1, p2);
    for (Q x1 = ceilq(c1 - s1); x1 <= floorq(c1 + s1); x1 += 1) {
      cur[1] = x1;
      const Q p1 = p2 + r.gn[1] * (c1 - x1) * (c1 - x1);
      const Q c0 = centre(0);
      const Q s0 = span(0, p1);
      for (Q x0 = ceilq(c0 - s0); x0 <= floorq(c0 + s0); x0 += 1) {
        V d = w;
        d = d - x0 * r.v[0];
        d = d - x1 * r.v[1];
        d = d - x2 * r.v[2];
        visit(x0, x1, x2, dot(d, d));
      }
    }
  }
}

// The shortest translate of w by the lattice. The search radius is the
// Babai nearest-plane residual with a binary128 rounding allowance only: a
// relative slack would scale with the unavoidable component along a long
// basis vector and enumerate the whole short sublattice under it.
inline Image closest(const Cell &cell, const Basis &r, V w) {
  Image out;
  if (!r.ok) {
    return out;
  }
  Q x[3] = {0, 0, 0};
  V res = w;
  for (int i = 2; i >= 0; i--) {
    x[i] = qround(dot(res, r.gs[i]) / r.gn[i]);
    res = res - x[i] * r.v[i];
  }
  const Q bound = dot(res, res) * (Q)(1 + 1e-30) + (Q)1e-300 * (Q)1e-300;
  Q best = bound;
  Q bestX[3] = {x[0], x[1], x[2]};
  enumerate(r, w, bound, [&](Q x0, Q x1, Q x2, Q n2) {
    if (n2 < best) {
      best = n2;
      bestX[0] = x0;
      bestX[1] = x1;
      bestX[2] = x2;
    }
  });
  I n[3] = {0, 0, 0};
  for (int i = 0; i < 3; i++) {
    if (bestX[i] == 0) {
      continue;
    }
    if (!r.exact[i] || !(qabs(bestX[i]) < (Q)0x1p60)) {
      return out;
    }
    for (int j = 0; j < 3; j++) {
      n[j] += (I)bestX[i] * r.coef[i][j];
    }
  }
  out.d = w - combo(cell, n);
  out.n2 = dot(out.d, out.d);
  for (int j = 0; j < 3; j++) {
    out.n[j] = n[j];
    out.x[j] = bestX[j];
  }
  out.ok = true;
  return out;
}

// True when another translate of w is within margin of best in squared
// length. Moving best by the shortest reduced vector v or -v changes the
// squared length by at most |v|^2, so a margin that large is always a tie.
inline bool nearTie(const Basis &r, V w, const Image &best, Q margin) {
  Q shortest = dot(r.v[0], r.v[0]);
  for (int i = 1; i < 3; i++) {
    shortest = std::min(shortest, dot(r.v[i], r.v[i]));
  }
  if (margin >= shortest) {
    return true;
  }
  bool tie = false;
  enumerate(r, w, best.n2 + margin, [&](Q x0, Q x1, Q x2, Q) {
    tie = tie || x0 != best.x[0] || x1 != best.x[1] || x2 != best.x[2];
  });
  return tie;
}

inline V fromDouble(const double v[3]) { return {(Q)v[0], (Q)v[1], (Q)v[2]}; }

// Per component k, the sum over i of |n_i| times the k-th component of
// lattice vector i: the size of the numbers a binary64 evaluation of
// n . (a, b, c) adds up.
inline std::array<double, 3> translationScale(const Cell &c, const I n[3]) {
  const double n0 = std::fabs((double)n[0]), n1 = std::fabs((double)n[1]),
               n2 = std::fabs((double)n[2]);
  return {n0 * std::fabs(c.lx) + n1 * std::fabs(c.xy) + n2 * std::fabs(c.xz),
          n1 * std::fabs(c.ly) + n2 * std::fabs(c.yz), n2 * std::fabs(c.lz)};
}

// Per component k, the sum over i of |t_i| times the k-th component of
// lattice vector i, for t the fractional coordinates of v: the size of the
// numbers a wrap through fractional coordinates adds up to get v.
inline std::array<double, 3> fractionalScale(const Cell &c, const double v[3]) {
  const Q tz = (Q)v[2] / (Q)c.lz;
  const Q ty = ((Q)v[1] - tz * (Q)c.yz) / (Q)c.ly;
  const Q tx = ((Q)v[0] - ty * (Q)c.xy - tz * (Q)c.xz) / (Q)c.lx;
  const double t0 = (double)qabs(tx), t1 = (double)qabs(ty), t2 = (double)qabs(tz);
  return {t0 * std::fabs(c.lx) + t1 * std::fabs(c.xy) + t2 * std::fabs(c.xz),
          t1 * std::fabs(c.ly) + t2 * std::fabs(c.yz), t2 * std::fabs(c.lz)};
}

// The largest component of raw - candidate - n . (a, b, c) for the nearest
// integer n, in units of 2^-52 times the largest |raw_k|, |candidate_k|,
// translation scale and fractional scale of raw and of the candidate. A
// lattice image computed in binary64 stays within a few of these units.
inline double offLattice(const Cell &c, const double raw[3], const double candidate[3]) {
  const V diff = fromDouble(raw) - fromDouble(candidate);
  I n[3];
  const Q nz = qround(diff.z / (Q)c.lz);
  const Q ny = qround((diff.y - nz * (Q)c.yz) / (Q)c.ly);
  const Q nx = qround((diff.x - ny * (Q)c.xy - nz * (Q)c.xz) / (Q)c.lx);
  const Q limit = (Q)1e30;
  if (!(qabs(nx) < limit && qabs(ny) < limit && qabs(nz) < limit)) {
    return INFINITY;
  }
  n[0] = (I)nx;
  n[1] = (I)ny;
  n[2] = (I)nz;
  const V r = diff - combo(c, n);
  const auto t = translationScale(c, n);
  const auto fr = fractionalScale(c, raw);
  const auto fc = fractionalScale(c, candidate);
  double scale = 0;
  for (int k = 0; k < 3; k++) {
    scale = std::max({scale, std::fabs(raw[k]), std::fabs(candidate[k]), t[k], fr[k], fc[k]});
  }
  const Q unit = (Q)scale * (Q)0x1p-52;
  if (!(unit > 0)) {
    return qabs(r.x) + qabs(r.y) + qabs(r.z) == 0 ? 0.0 : INFINITY;
  }
  return (double)(std::max({qabs(r.x), qabs(r.y), qabs(r.z)}) / unit);
}

}  // namespace micoracle
