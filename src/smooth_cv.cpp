//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <smooth_cv.hpp>

#include <franzblau.hpp>
#include <generic.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace {

// Forward-mode value with derivatives in 38 slots: the 36 bonds of a
// hexagon pair and its two ring weights.
struct Dual {
  static constexpr int kSlots = 38;
  double v = 0.0;
  std::array<double, kSlots> g{};
  Dual() = default;
  Dual(double value) : v(value) {}
  static Dual seed(double value, int slot) {
    Dual d(value);
    d.g[static_cast<std::size_t>(slot)] = 1.0;
    return d;
  }
};

Dual operator*(const Dual &a, const Dual &b) {
  Dual r(a.v * b.v);
  for (std::size_t k = 0; k < Dual::kSlots; ++k) {
    r.g[k] = a.g[k] * b.v + a.v * b.g[k];
  }
  return r;
}

Dual operator-(double a, const Dual &b) {
  Dual r(a - b.v);
  for (std::size_t k = 0; k < Dual::kSlots; ++k) {
    r.g[k] = -b.g[k];
  }
  return r;
}

// The soft form of the test in ring::basalConditions. bond[i][j] is the
// weight of the bond between member i of the first ring and member j of the
// second. For q in {0, 1} the members q + 1, q + 3, q + 5 of the second ring
// must have no bond to the first ring, and some alternate triplet of the
// first ring must be bonded one to one to the members q, q + 2, q + 4. Each
// "and" is a product and each "some" is 1 - prod (1 - term), so 0/1 weights
// give the boolean test exactly.
template <class T> T basalSoft(const T (&bond)[6][6]) {
  static constexpr int kPerm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                                      {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  T noneQ(1.0);
  for (int q = 0; q < 2; ++q) {
    T quiet(1.0);
    for (int t = 0; t < 3; ++t) {
      const int j = (q + 1 + 2 * t) % 6;
      for (int i = 0; i < 6; ++i) {
        quiet = quiet * (1.0 - bond[i][j]);
      }
    }
    T noneMatch(1.0);
    for (int p = 0; p < 2; ++p) {
      for (const auto &s : kPerm) {
        noneMatch = noneMatch * (1.0 - bond[p][q + 2 * s[0]] * bond[p + 2][q + 2 * s[1]] *
                                           bond[p + 4][q + 2 * s[2]]);
      }
    }
    noneQ = noneQ * (1.0 - quiet * (1.0 - noneMatch));
  }
  return 1.0 - noneQ;
}

template <class T>
T hexPairTerm(const smooth::Topology::HexPair &pair, const std::vector<double> &s, const T &wa,
              const T &wb, bool seedBonds) {
  T fwd[6][6];
  T rev[6][6];
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      const int e = pair.bond[static_cast<std::size_t>(6 * i + j)];
      const double weight = e < 0 ? 0.0 : s[static_cast<std::size_t>(e)];
      T b(weight);
      if constexpr (std::is_same_v<T, Dual>) {
        if (seedBonds) {
          b = Dual::seed(weight, 6 * i + j);
        }
      }
      fwd[i][j] = b;
      rev[j][i] = b;
    }
  }
  const T either = 1.0 - (1.0 - basalSoft(fwd)) * (1.0 - basalSoft(rev));
  return wa * wb * either;
}

// The orientation primitive::ringNetwork reports: lowest member first,
// toward the lesser of its two ring neighbours.
std::vector<int> canonicalCycle(std::vector<int> ring) {
  const auto low = std::min_element(ring.begin(), ring.end());
  std::rotate(ring.begin(), low, ring.end());
  if (ring.size() > 2 && ring.back() < ring[1]) {
    std::reverse(ring.begin() + 1, ring.end());
  }
  return ring;
}

double logistic(double z) {
  if (z >= 0.0) {
    return 1.0 / (1.0 + std::exp(-z));
  }
  const double e = std::exp(z);
  return e / (1.0 + e);
}

bool sizeMatches(int ringSize, int size) { return size == 0 || ringSize == size; }

} // namespace

namespace smooth {

void RationalSwitch::validate() const {
  const bool finite = std::isfinite(r0) && std::isfinite(d0) && std::isfinite(dOn) &&
                      std::isfinite(dMax);
  if (!finite) {
    throw std::invalid_argument("switching parameters must be finite");
  }
  if (!(r0 > 0.0)) {
    throw std::invalid_argument("r0 must be positive");
  }
  if (n < 1 || m <= n) {
    throw std::invalid_argument("the exponents need n >= 1 and m > n");
  }
  if (!(d0 >= 0.0 && d0 <= dOn && dOn < dMax)) {
    throw std::invalid_argument("the radii need 0 <= d0 <= d_on < d_max");
  }
}

RationalSwitch makeSwitch(double r0, double d0, int n, int m, double dOn, double dMax) {
  RationalSwitch sw;
  sw.r0 = r0;
  sw.d0 = d0;
  sw.n = n;
  sw.m = m == 0 ? 2 * n : m;
  sw.dOn = dOn;
  sw.dMax = dMax;
  sw.validate();
  return sw;
}

double RationalSwitch::operator()(double r, double &dsdr) const {
  dsdr = 0.0;
  if (r >= dMax) {
    return 0.0;
  }
  double rational = 1.0;
  double dRational = 0.0;
  if (r > d0) {
    const double x = (r - d0) / r0;
    // Horner for 1 + x + ... + x^(k-1) and its derivative.
    auto geometric = [x](int k, double &deriv) {
      double sum = 0.0;
      deriv = 0.0;
      for (int i = 0; i < k; ++i) {
        deriv = deriv * x + sum;
        sum = sum * x + 1.0;
      }
      return sum;
    };
    double dn = 0.0;
    double dm = 0.0;
    const double sn = geometric(n, dn);
    const double sm = geometric(m, dm);
    rational = sn / sm;
    dRational = (dn * sm - sn * dm) / (sm * sm) / r0;
  }
  double taper = 1.0;
  double dTaper = 0.0;
  if (r > dOn) {
    const double width = dMax - dOn;
    const double u = (r - dOn) / width;
    const double u2 = u * u;
    taper = 1.0 - u2 * u * (10.0 - 15.0 * u + 6.0 * u2);
    const double v = 1.0 - u;
    dTaper = -30.0 * u2 * v * v / width;
  }
  dsdr = dRational * taper + rational * dTaper;
  return rational * taper;
}

int Topology::edgeOf(int i, int j) const {
  if (i < 0 || i >= nAtoms) {
    return -1;
  }
  const auto begin = adjAtom.begin() + adjOff[static_cast<std::size_t>(i)];
  const auto end = adjAtom.begin() + adjOff[static_cast<std::size_t>(i) + 1];
  const auto it = std::lower_bound(begin, end, j);
  if (it == end || *it != j) {
    return -1;
  }
  return adjEdge[static_cast<std::size_t>(it - adjAtom.begin())];
}

Topology buildTopology(const std::vector<std::vector<int>> &rows, int maxRing) {
  if (maxRing < 3 || maxRing > 10) {
    throw std::invalid_argument("max_ring must be from 3 to 10");
  }
  Topology topo;
  topo.nAtoms = static_cast<int>(rows.size());
  topo.maxRing = maxRing;
  const int n = topo.nAtoms;

  std::vector<std::vector<int>> nbrs(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const auto &row = rows[static_cast<std::size_t>(i)];
    for (std::size_t k = 1; k < row.size(); ++k) {
      const int j = row[k];
      if (j < 0 || j >= n || j == i) {
        throw std::invalid_argument("a neighbour row names an atom outside the graph");
      }
      nbrs[static_cast<std::size_t>(i)].push_back(j);
    }
  }
  for (int i = 0; i < n; ++i) {
    auto &row = nbrs[static_cast<std::size_t>(i)];
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
  }
  for (int i = 0; i < n; ++i) {
    for (const int j : nbrs[static_cast<std::size_t>(i)]) {
      const auto &back = nbrs[static_cast<std::size_t>(j)];
      if (!std::binary_search(back.begin(), back.end(), i)) {
        throw std::invalid_argument("neighbour rows are not symmetric");
      }
      if (i < j) {
        topo.edges.push_back({i, j});
      }
    }
  }
  topo.adjOff.assign(static_cast<std::size_t>(n) + 1, 0);
  for (int i = 0; i < n; ++i) {
    topo.adjOff[static_cast<std::size_t>(i) + 1] =
        topo.adjOff[static_cast<std::size_t>(i)] +
        static_cast<int>(nbrs[static_cast<std::size_t>(i)].size());
  }
  topo.adjAtom.resize(static_cast<std::size_t>(topo.adjOff.back()));
  topo.adjEdge.assign(topo.adjAtom.size(), -1);
  for (int i = 0; i < n; ++i) {
    std::copy(nbrs[static_cast<std::size_t>(i)].begin(), nbrs[static_cast<std::size_t>(i)].end(),
              topo.adjAtom.begin() + topo.adjOff[static_cast<std::size_t>(i)]);
  }
  for (int e = 0; e < static_cast<int>(topo.edges.size()); ++e) {
    const auto [i, j] = topo.edges[static_cast<std::size_t>(e)];
    for (const auto [a, b] : {std::pair<int, int>{i, j}, std::pair<int, int>{j, i}}) {
      const auto begin = topo.adjAtom.begin() + topo.adjOff[static_cast<std::size_t>(a)];
      const auto end = topo.adjAtom.begin() + topo.adjOff[static_cast<std::size_t>(a) + 1];
      const auto it = std::lower_bound(begin, end, b);
      topo.adjEdge[static_cast<std::size_t>(it - topo.adjAtom.begin())] = e;
    }
  }

  // Every simple cycle, from the backtracking search behind the Franzblau
  // rings. The shortest-path filter is replaced by the soft factors below.
  std::vector<std::vector<int>> graphRows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    auto &row = graphRows[static_cast<std::size_t>(i)];
    row.push_back(i);
    row.insert(row.end(), nbrs[static_cast<std::size_t>(i)].begin(),
               nbrs[static_cast<std::size_t>(i)].end());
  }
  auto graph = primitive::countAllRingsFromIndex(graphRows, maxRing);
  std::vector<std::vector<int>> cycles;
  cycles.reserve(graph.rings.size());
  for (auto &ring : graph.rings) {
    if (ring.size() >= 3) {
      cycles.push_back(canonicalCycle(std::move(ring)));
    }
  }
  std::sort(cycles.begin(), cycles.end(), [](const auto &a, const auto &b) {
    return a.size() != b.size() ? a.size() < b.size() : a < b;
  });
  cycles.erase(std::unique(cycles.begin(), cycles.end()), cycles.end());

  std::vector<int> path;
  std::vector<char> onPath(static_cast<std::size_t>(n), 0);
  for (const auto &ring : cycles) {
    const int k = static_cast<int>(ring.size());
    for (int t = 0; t < k; ++t) {
      topo.ringAtoms.push_back(ring[static_cast<std::size_t>(t)]);
      topo.ringEdges.push_back(topo.edgeOf(ring[static_cast<std::size_t>(t)],
                                           ring[static_cast<std::size_t>((t + 1) % k)]));
    }
    topo.ringOff.push_back(static_cast<int>(topo.ringAtoms.size()));
    // Shortcuts: simple paths between members a and b shorter than their
    // separation d around the ring. Paths may pass other members; such a
    // path holds a shorter shortcut of its own, so the product is unchanged
    // in the 0/1 limit and stays a polynomial in the weights.
    for (int pa = 0; pa < k; ++pa) {
      for (int pb = pa + 1; pb < k; ++pb) {
        const int d = std::min(pb - pa, k - (pb - pa));
        if (d < 2) {
          continue;
        }
        const int from = ring[static_cast<std::size_t>(pa)];
        const int goal = ring[static_cast<std::size_t>(pb)];
        const int limit = d - 1;
        path.clear();
        onPath[static_cast<std::size_t>(from)] = 1;
        auto walk = [&](auto &&self, int u) -> void {
          for (int a = topo.adjOff[static_cast<std::size_t>(u)];
               a < topo.adjOff[static_cast<std::size_t>(u) + 1]; ++a) {
            const int v = topo.adjAtom[static_cast<std::size_t>(a)];
            const int e = topo.adjEdge[static_cast<std::size_t>(a)];
            if (onPath[static_cast<std::size_t>(v)] != 0) {
              continue;
            }
            if (v == goal) {
              topo.pathEdges.insert(topo.pathEdges.end(), path.begin(), path.end());
              topo.pathEdges.push_back(e);
              topo.pathEdgeOff.push_back(static_cast<int>(topo.pathEdges.size()));
              continue;
            }
            if (static_cast<int>(path.size()) + 1 < limit) {
              onPath[static_cast<std::size_t>(v)] = 1;
              path.push_back(e);
              self(self, v);
              path.pop_back();
              onPath[static_cast<std::size_t>(v)] = 0;
            }
          }
        };
        walk(walk, from);
        onPath[static_cast<std::size_t>(from)] = 0;
      }
    }
    topo.pathOff.push_back(static_cast<int>(topo.pathEdgeOff.size()) - 1);
  }

  // Hexagon pairs that can pass the basal test: vertex-disjoint, with every
  // bond of at least one matched triplet in the graph. A pair without one
  // has weight zero until the next rebuild.
  if (maxRing >= 6) {
    std::vector<std::vector<int>> hexOf(static_cast<std::size_t>(n));
    std::vector<int> hexes;
    for (int r = 0; r < topo.ringCount(); ++r) {
      if (topo.ringSize(r) == 6) {
        hexes.push_back(r);
        for (int t = 0; t < 6; ++t) {
          hexOf[static_cast<std::size_t>(topo.ringAtoms[static_cast<std::size_t>(
                    topo.ringOff[static_cast<std::size_t>(r)] + t)])]
              .push_back(r);
        }
      }
    }
    static constexpr int kPerm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                                        {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    auto reachable = [](const std::array<int, 36> &bond, bool transpose) {
      auto at = [&](int i, int j) {
        return transpose ? bond[static_cast<std::size_t>(6 * j + i)]
                         : bond[static_cast<std::size_t>(6 * i + j)];
      };
      for (int q = 0; q < 2; ++q) {
        for (int p = 0; p < 2; ++p) {
          for (const auto &s : kPerm) {
            if (at(p, q + 2 * s[0]) >= 0 && at(p + 2, q + 2 * s[1]) >= 0 &&
                at(p + 4, q + 2 * s[2]) >= 0) {
              return true;
            }
          }
        }
      }
      return false;
    };
    for (const int a : hexes) {
      const int *ra = topo.ringAtoms.data() + topo.ringOff[static_cast<std::size_t>(a)];
      std::set<int> partners;
      for (int t = 0; t < 6; ++t) {
        const int atom = ra[t];
        for (int x = topo.adjOff[static_cast<std::size_t>(atom)];
             x < topo.adjOff[static_cast<std::size_t>(atom) + 1]; ++x) {
          const int nbr = topo.adjAtom[static_cast<std::size_t>(x)];
          for (const int b : hexOf[static_cast<std::size_t>(nbr)]) {
            if (b > a) {
              partners.insert(b);
            }
          }
        }
      }
      for (const int b : partners) {
        const int *rb = topo.ringAtoms.data() + topo.ringOff[static_cast<std::size_t>(b)];
        bool disjoint = true;
        for (int i = 0; i < 6 && disjoint; ++i) {
          for (int j = 0; j < 6; ++j) {
            if (ra[i] == rb[j]) {
              disjoint = false;
              break;
            }
          }
        }
        if (!disjoint) {
          continue;
        }
        Topology::HexPair pair;
        pair.a = a;
        pair.b = b;
        for (int i = 0; i < 6; ++i) {
          for (int j = 0; j < 6; ++j) {
            pair.bond[static_cast<std::size_t>(6 * i + j)] = topo.edgeOf(ra[i], rb[j]);
          }
        }
        if (reachable(pair.bond, false) || reachable(pair.bond, true)) {
          topo.hexPairs.push_back(pair);
        }
      }
    }
  }
  return topo;
}

std::vector<double> ringWeights(const Topology &topo, const std::vector<double> &s) {
  std::vector<double> w(static_cast<std::size_t>(topo.ringCount()), 0.0);
  for (int r = 0; r < topo.ringCount(); ++r) {
    double bonds = 1.0;
    for (int t = topo.ringOff[static_cast<std::size_t>(r)];
         t < topo.ringOff[static_cast<std::size_t>(r) + 1]; ++t) {
      bonds *= s[static_cast<std::size_t>(topo.ringEdges[static_cast<std::size_t>(t)])];
    }
    if (bonds == 0.0) {
      continue;
    }
    double open = 1.0;
    for (int p = topo.pathOff[static_cast<std::size_t>(r)];
         p < topo.pathOff[static_cast<std::size_t>(r) + 1] && open != 0.0; ++p) {
      double through = 1.0;
      for (int q = topo.pathEdgeOff[static_cast<std::size_t>(p)];
           q < topo.pathEdgeOff[static_cast<std::size_t>(p) + 1]; ++q) {
        through *= s[static_cast<std::size_t>(topo.pathEdges[static_cast<std::size_t>(q)])];
      }
      open *= 1.0 - through;
    }
    w[static_cast<std::size_t>(r)] = bonds * open;
  }
  return w;
}

namespace {

// prefix[k] * suffix[k] is the product of every factor but factor k.
void leaveOneOut(const std::vector<double> &factors, std::vector<double> &out) {
  const std::size_t n = factors.size();
  out.assign(n, 1.0);
  double run = 1.0;
  for (std::size_t k = 0; k < n; ++k) {
    out[k] = run;
    run *= factors[k];
  }
  run = 1.0;
  for (std::size_t k = n; k-- > 0;) {
    out[k] *= run;
    run *= factors[k];
  }
}

} // namespace

void ringWeightsAdjoint(const Topology &topo, const std::vector<double> &s,
                        const std::vector<double> &wbar, std::vector<double> &sbar) {
  std::vector<double> factors;
  std::vector<double> others;
  std::vector<double> pathOthers;
  std::vector<double> edgeOthers;
  std::vector<double> through;
  for (int r = 0; r < topo.ringCount(); ++r) {
    const double gbar = wbar[static_cast<std::size_t>(r)];
    if (gbar == 0.0) {
      continue;
    }
    const int r0 = topo.ringOff[static_cast<std::size_t>(r)];
    const int r1 = topo.ringOff[static_cast<std::size_t>(r) + 1];
    factors.clear();
    for (int t = r0; t < r1; ++t) {
      factors.push_back(s[static_cast<std::size_t>(topo.ringEdges[static_cast<std::size_t>(t)])]);
    }
    leaveOneOut(factors, others);
    const double bonds = others.front() * factors.front();
    const int p0 = topo.pathOff[static_cast<std::size_t>(r)];
    const int p1 = topo.pathOff[static_cast<std::size_t>(r) + 1];
    through.clear();
    factors.clear();
    for (int p = p0; p < p1; ++p) {
      double prod = 1.0;
      for (int q = topo.pathEdgeOff[static_cast<std::size_t>(p)];
           q < topo.pathEdgeOff[static_cast<std::size_t>(p) + 1]; ++q) {
        prod *= s[static_cast<std::size_t>(topo.pathEdges[static_cast<std::size_t>(q)])];
      }
      through.push_back(prod);
      factors.push_back(1.0 - prod);
    }
    leaveOneOut(factors, pathOthers);
    const double open = factors.empty() ? 1.0 : pathOthers.front() * factors.front();
    for (int t = r0; t < r1; ++t) {
      sbar[static_cast<std::size_t>(topo.ringEdges[static_cast<std::size_t>(t)])] +=
          gbar * open * others[static_cast<std::size_t>(t - r0)];
    }
    if (bonds == 0.0) {
      continue;
    }
    for (int p = p0; p < p1; ++p) {
      const double dPath = -gbar * bonds * pathOthers[static_cast<std::size_t>(p - p0)];
      if (dPath == 0.0) {
        continue;
      }
      const int q0 = topo.pathEdgeOff[static_cast<std::size_t>(p)];
      const int q1 = topo.pathEdgeOff[static_cast<std::size_t>(p) + 1];
      factors.clear();
      for (int q = q0; q < q1; ++q) {
        factors.push_back(s[static_cast<std::size_t>(topo.pathEdges[static_cast<std::size_t>(q)])]);
      }
      leaveOneOut(factors, edgeOthers);
      for (int q = q0; q < q1; ++q) {
        sbar[static_cast<std::size_t>(topo.pathEdges[static_cast<std::size_t>(q)])] +=
            dPath * edgeOthers[static_cast<std::size_t>(q - q0)];
      }
    }
  }
}

double hexBasalPairs(const Topology &topo, const std::vector<double> &s,
                     const std::vector<double> &w, std::vector<double> *sbar,
                     std::vector<double> *wbar) {
  double total = 0.0;
  for (const auto &pair : topo.hexPairs) {
    const double wa = w[static_cast<std::size_t>(pair.a)];
    const double wb = w[static_cast<std::size_t>(pair.b)];
    if (sbar == nullptr || wbar == nullptr) {
      if (wa != 0.0 && wb != 0.0) {
        total += hexPairTerm<double>(pair, s, wa, wb, false);
      }
      continue;
    }
    const Dual term =
        hexPairTerm<Dual>(pair, s, Dual::seed(wa, 36), Dual::seed(wb, 37), true);
    total += term.v;
    for (int k = 0; k < 36; ++k) {
      const int e = pair.bond[static_cast<std::size_t>(k)];
      if (e >= 0) {
        (*sbar)[static_cast<std::size_t>(e)] += term.g[static_cast<std::size_t>(k)];
      }
    }
    (*wbar)[static_cast<std::size_t>(pair.a)] += term.g[36];
    (*wbar)[static_cast<std::size_t>(pair.b)] += term.g[37];
  }
  return total;
}

std::vector<double> ringParticipation(const Topology &topo, const std::vector<double> &w,
                                      int size) {
  std::vector<double> none(static_cast<std::size_t>(topo.nAtoms), 1.0);
  for (int r = 0; r < topo.ringCount(); ++r) {
    if (!sizeMatches(topo.ringSize(r), size)) {
      continue;
    }
    for (int t = topo.ringOff[static_cast<std::size_t>(r)];
         t < topo.ringOff[static_cast<std::size_t>(r) + 1]; ++t) {
      none[static_cast<std::size_t>(topo.ringAtoms[static_cast<std::size_t>(t)])] *=
          1.0 - w[static_cast<std::size_t>(r)];
    }
  }
  for (double &p : none) {
    p = 1.0 - p;
  }
  return none;
}

void ringParticipationAdjoint(const Topology &topo, const std::vector<double> &w, int size,
                              const std::vector<double> &pbar, std::vector<double> &wbar) {
  std::vector<std::vector<int>> through(static_cast<std::size_t>(topo.nAtoms));
  for (int r = 0; r < topo.ringCount(); ++r) {
    if (!sizeMatches(topo.ringSize(r), size)) {
      continue;
    }
    for (int t = topo.ringOff[static_cast<std::size_t>(r)];
         t < topo.ringOff[static_cast<std::size_t>(r) + 1]; ++t) {
      through[static_cast<std::size_t>(topo.ringAtoms[static_cast<std::size_t>(t)])].push_back(r);
    }
  }
  std::vector<double> factors;
  std::vector<double> others;
  for (int i = 0; i < topo.nAtoms; ++i) {
    const double g = pbar[static_cast<std::size_t>(i)];
    const auto &rings = through[static_cast<std::size_t>(i)];
    if (g == 0.0 || rings.empty()) {
      continue;
    }
    factors.clear();
    for (const int r : rings) {
      factors.push_back(1.0 - w[static_cast<std::size_t>(r)]);
    }
    leaveOneOut(factors, others);
    for (std::size_t k = 0; k < rings.size(); ++k) {
      wbar[static_cast<std::size_t>(rings[k])] += g * others[k];
    }
  }
}

Workspace::Workspace(const RationalSwitch &sw, int maxRing) : sw_(sw), maxRing_(maxRing) {
  sw_.validate();
  if (maxRing < 3 || maxRing > 10) {
    throw std::invalid_argument("max_ring must be from 3 to 10");
  }
}

void Workspace::rebuild(const Cloud &cloud, std::vector<int> selected, double skin) {
  if (!std::isfinite(skin) || skin < 0.0) {
    throw std::invalid_argument("skin must be finite and not negative");
  }
  const gen::FracBox box = gen::makeFracBox(cloud);
  if (!box.ok) {
    throw std::invalid_argument("the cell has no volume");
  }
  const double listCutoff = sw_.dMax + skin;
  if (!(listCutoff < box.halfMin)) {
    throw std::invalid_argument("d_max + skin must be below half the narrowest cell width");
  }
  const int n = cloud.nop;
  std::vector<int> mask(static_cast<std::size_t>(n), 0);
  for (const int i : selected) {
    if (i < 0 || i >= n) {
      throw std::invalid_argument("a selected atom is outside the frame");
    }
    mask[static_cast<std::size_t>(i)] = 1;
  }
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    rows[static_cast<std::size_t>(i)].push_back(i);
  }
  if (!selected.empty()) {
    linkcell::PairCsr csr;
    nneigh::linkcellCutoffCsr(cloud, listCutoff, false, {true, true, true}, mask, 1, csr);
    for (int i = 0; i < n; ++i) {
      for (auto k = csr.offsets[static_cast<std::size_t>(i)];
           k < csr.offsets[static_cast<std::size_t>(i) + 1]; ++k) {
        rows[static_cast<std::size_t>(i)].push_back(
            static_cast<int>(csr.cols[static_cast<std::size_t>(k)]));
      }
    }
  }
  Topology topo = buildTopology(rows, maxRing_);
  reference_.resize(3 * static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    reference_[3 * static_cast<std::size_t>(i)] = p.x;
    reference_[3 * static_cast<std::size_t>(i) + 1] = p.y;
    reference_[3 * static_cast<std::size_t>(i) + 2] = p.z;
  }
  topo_ = std::move(topo);
  selected_ = std::move(selected);
  box_ = cloud.box;
  boxLow_ = cloud.boxLow;
  skin_ = skin;
  built_ = true;
}

bool Workspace::stale(const Cloud &cloud) const {
  if (!built_ || cloud.nop != topo_.nAtoms || cloud.box != box_ || cloud.boxLow != boxLow_) {
    return true;
  }
  const gen::FracBox box = gen::makeFracBox(cloud);
  const double limit = 0.25 * skin_ * skin_;
  for (const int i : selected_) {
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    const double *r = reference_.data() + 3 * static_cast<std::size_t>(i);
    const auto d = gen::fracDelta(box, p.x, p.y, p.z, r[0], r[1], r[2]);
    if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > limit) {
      return true;
    }
  }
  return false;
}

void Workspace::checkKind(Kind kind, int arg) const {
  const bool size = arg == 0 || (arg >= 3 && arg <= maxRing_);
  switch (kind) {
  case Kind::RingCount:
  case Kind::RingAtoms:
    if (!size) {
      throw std::invalid_argument("the ring size must be 0 or from 3 to max_ring");
    }
    return;
  case Kind::HexBasalPairs:
    if (arg != 0 || maxRing_ < 6) {
      throw std::invalid_argument("hexagon pairs take arg 0 and need max_ring >= 6");
    }
    return;
  case Kind::GyrationRg2:
  case Kind::GyrationVolume:
    if (!(arg == -1 || size)) {
      throw std::invalid_argument("gyration weights are -1, 0 or a ring size up to max_ring");
    }
    return;
  case Kind::State:
    if (arg < 0 || arg >= stateCount()) {
      throw std::invalid_argument("no state with that index");
    }
    return;
  }
  throw std::invalid_argument("unknown collective variable kind");
}

void Workspace::edgeWeights(const Cloud &cloud) {
  const gen::FracBox box = gen::makeFracBox(cloud);
  const std::size_t ne = topo_.edges.size();
  s_.assign(ne, 0.0);
  ds_.assign(ne, 0.0);
  unit_.assign(ne, {0.0, 0.0, 0.0});
  for (std::size_t e = 0; e < ne; ++e) {
    const auto &pi = cloud.pts[static_cast<std::size_t>(topo_.edges[e][0])];
    const auto &pj = cloud.pts[static_cast<std::size_t>(topo_.edges[e][1])];
    const auto d = gen::fracDelta(box, pj.x, pj.y, pj.z, pi.x, pi.y, pi.z);
    const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    double dsdr = 0.0;
    s_[e] = sw_(r, dsdr);
    ds_[e] = dsdr;
    if (r > 0.0) {
      unit_[e] = {d[0] / r, d[1] / r, d[2] / r};
    }
  }
}

void Workspace::edgesToAtoms(const std::vector<double> &sbar, std::vector<double> &grad) const {
  for (std::size_t e = 0; e < topo_.edges.size(); ++e) {
    const double g = sbar[e] * ds_[e];
    if (g == 0.0) {
      continue;
    }
    const std::size_t i = 3 * static_cast<std::size_t>(topo_.edges[e][0]);
    const std::size_t j = 3 * static_cast<std::size_t>(topo_.edges[e][1]);
    for (std::size_t k = 0; k < 3; ++k) {
      grad[j + k] += g * unit_[e][k];
      grad[i + k] -= g * unit_[e][k];
    }
  }
}

double Workspace::gyration(const Cloud &cloud, bool volume, int arg, std::vector<double> *grad) {
  const int n = topo_.nAtoms;
  std::vector<double> weight(static_cast<std::size_t>(n), 0.0);
  std::vector<double> w;
  std::vector<double> participation;
  if (arg >= 0) {
    w = ringWeights(topo_, s_);
    participation = ringParticipation(topo_, w, arg);
  }
  for (const int i : selected_) {
    weight[static_cast<std::size_t>(i)] =
        arg >= 0 ? participation[static_cast<std::size_t>(i)] : 1.0;
  }
  if (selected_.empty()) {
    return 0.0;
  }
  // Unwrap along graph bonds, which are shorter than half the cell; a
  // component the graph does not reach is placed by the minimum image of
  // its first atom from the first selected atom.
  const gen::FracBox box = gen::makeFracBox(cloud);
  std::vector<std::array<double, 3>> u(static_cast<std::size_t>(n));
  std::vector<char> seen(static_cast<std::size_t>(n), 0);
  const auto &origin = cloud.pts[static_cast<std::size_t>(selected_.front())];
  std::vector<int> queue;
  for (const int start : selected_) {
    if (seen[static_cast<std::size_t>(start)] != 0) {
      continue;
    }
    const auto &ps = cloud.pts[static_cast<std::size_t>(start)];
    const auto d0 = gen::fracDelta(box, ps.x, ps.y, ps.z, origin.x, origin.y, origin.z);
    u[static_cast<std::size_t>(start)] = {origin.x + d0[0], origin.y + d0[1], origin.z + d0[2]};
    seen[static_cast<std::size_t>(start)] = 1;
    queue.assign(1, start);
    for (std::size_t head = 0; head < queue.size(); ++head) {
      const int a = queue[head];
      const auto &pa = cloud.pts[static_cast<std::size_t>(a)];
      for (int x = topo_.adjOff[static_cast<std::size_t>(a)];
           x < topo_.adjOff[static_cast<std::size_t>(a) + 1]; ++x) {
        const int b = topo_.adjAtom[static_cast<std::size_t>(x)];
        if (seen[static_cast<std::size_t>(b)] != 0) {
          continue;
        }
        const auto &pb = cloud.pts[static_cast<std::size_t>(b)];
        const auto d = gen::fracDelta(box, pb.x, pb.y, pb.z, pa.x, pa.y, pa.z);
        const auto &ua = u[static_cast<std::size_t>(a)];
        u[static_cast<std::size_t>(b)] = {ua[0] + d[0], ua[1] + d[1], ua[2] + d[2]};
        seen[static_cast<std::size_t>(b)] = 1;
        queue.push_back(b);
      }
    }
  }
  double total = 0.0;
  std::array<double, 3> c{0.0, 0.0, 0.0};
  for (const int i : selected_) {
    const double wi = weight[static_cast<std::size_t>(i)];
    total += wi;
    for (std::size_t k = 0; k < 3; ++k) {
      c[k] += wi * u[static_cast<std::size_t>(i)][k];
    }
  }
  if (!(total > 0.0)) {
    return 0.0;
  }
  for (double &ck : c) {
    ck /= total;
  }
  double g[3][3] = {};
  for (const int i : selected_) {
    const double wi = weight[static_cast<std::size_t>(i)];
    const auto &ui = u[static_cast<std::size_t>(i)];
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        g[a][b] += wi * (ui[static_cast<std::size_t>(a)] - c[static_cast<std::size_t>(a)]) *
                   (ui[static_cast<std::size_t>(b)] - c[static_cast<std::size_t>(b)]);
      }
    }
  }
  for (auto &row : g) {
    for (double &x : row) {
      x /= total;
    }
  }
  const double rg2 = g[0][0] + g[1][1] + g[2][2];
  const double cof[3][3] = {{g[1][1] * g[2][2] - g[1][2] * g[2][1],
                             g[0][2] * g[2][1] - g[0][1] * g[2][2],
                             g[0][1] * g[1][2] - g[0][2] * g[1][1]},
                            {g[1][2] * g[2][0] - g[1][0] * g[2][2],
                             g[0][0] * g[2][2] - g[0][2] * g[2][0],
                             g[0][2] * g[1][0] - g[0][0] * g[1][2]},
                            {g[1][0] * g[2][1] - g[1][1] * g[2][0],
                             g[0][1] * g[2][0] - g[0][0] * g[2][1],
                             g[0][0] * g[1][1] - g[0][1] * g[1][0]}};
  const double det = g[0][0] * cof[0][0] + g[0][1] * cof[1][0] + g[0][2] * cof[2][0];
  // A uniform solid ellipsoid with semi-axes a, b, c has gyration
  // eigenvalues a^2/5, b^2/5, c^2/5, so its volume is this multiple of
  // sqrt(det G).
  const double kEllipsoid = 4.0 * 3.14159265358979323846 / 3.0 * std::pow(5.0, 1.5);
  const double vol = det > 0.0 ? kEllipsoid * std::sqrt(det) : 0.0;
  const double result = volume ? vol : rg2;
  if (grad == nullptr || (volume && !(det > 0.0))) {
    return result;
  }
  std::vector<double> wbar(static_cast<std::size_t>(n), 0.0);
  for (const int i : selected_) {
    const double wi = weight[static_cast<std::size_t>(i)];
    const auto &ui = u[static_cast<std::size_t>(i)];
    const double x[3] = {ui[0] - c[0], ui[1] - c[1], ui[2] - c[2]};
    double *gi = grad->data() + 3 * static_cast<std::size_t>(i);
    if (!volume) {
      for (int k = 0; k < 3; ++k) {
        gi[k] += 2.0 * wi * x[k] / total;
      }
      wbar[static_cast<std::size_t>(i)] = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - rg2) / total;
    } else {
      double gx[3];
      for (int a = 0; a < 3; ++a) {
        gx[a] = (cof[a][0] * x[0] + cof[a][1] * x[1] + cof[a][2] * x[2]) / det;
      }
      for (int k = 0; k < 3; ++k) {
        gi[k] += vol * wi * gx[k] / total;
      }
      wbar[static_cast<std::size_t>(i)] =
          0.5 * vol * (x[0] * gx[0] + x[1] * gx[1] + x[2] * gx[2] - 3.0) / total;
    }
  }
  if (arg >= 0) {
    std::vector<double> ringBar(w.size(), 0.0);
    ringParticipationAdjoint(topo_, w, arg, wbar, ringBar);
    std::vector<double> sbar(s_.size(), 0.0);
    ringWeightsAdjoint(topo_, s_, ringBar, sbar);
    edgesToAtoms(sbar, *grad);
  }
  return result;
}

double Workspace::channel(const Cloud &cloud, Kind kind, int arg, std::vector<double> *grad) {
  if (grad != nullptr) {
    grad->assign(3 * static_cast<std::size_t>(cloud.nop), 0.0);
  }
  if (kind == Kind::GyrationRg2 || kind == Kind::GyrationVolume) {
    return gyration(cloud, kind == Kind::GyrationVolume, arg, grad);
  }
  const std::vector<double> w = ringWeights(topo_, s_);
  std::vector<double> wbar(w.size(), 0.0);
  std::vector<double> sbar(s_.size(), 0.0);
  double total = 0.0;
  if (kind == Kind::RingCount) {
    for (int r = 0; r < topo_.ringCount(); ++r) {
      if (sizeMatches(topo_.ringSize(r), arg)) {
        total += w[static_cast<std::size_t>(r)];
        wbar[static_cast<std::size_t>(r)] = 1.0;
      }
    }
  } else if (kind == Kind::RingAtoms) {
    const auto p = ringParticipation(topo_, w, arg);
    std::vector<double> pbar(p.size(), 0.0);
    for (const int i : selected_) {
      total += p[static_cast<std::size_t>(i)];
      pbar[static_cast<std::size_t>(i)] = 1.0;
    }
    if (grad != nullptr) {
      ringParticipationAdjoint(topo_, w, arg, pbar, wbar);
    }
  } else if (kind == Kind::HexBasalPairs) {
    total = hexBasalPairs(topo_, s_, w, grad != nullptr ? &sbar : nullptr,
                          grad != nullptr ? &wbar : nullptr);
  }
  if (grad != nullptr) {
    ringWeightsAdjoint(topo_, s_, wbar, sbar);
    edgesToAtoms(sbar, *grad);
  }
  return total;
}

void Workspace::addClause(int state, Kind kind, int arg, double threshold, double width,
                          bool above) {
  if (state < 0 || state > stateCount()) {
    throw std::invalid_argument("a clause opens the next state or extends an existing one");
  }
  if (!std::isfinite(threshold) || !std::isfinite(width) || width < 0.0) {
    throw std::invalid_argument("threshold and width must be finite, width not negative");
  }
  if (kind == Kind::State) {
    if (arg < 0 || arg >= state) {
      throw std::invalid_argument("a state clause refers to an earlier state");
    }
  } else {
    checkKind(kind, arg);
    if (!(width > 0.0)) {
      throw std::invalid_argument("width must be positive");
    }
  }
  if (state == stateCount()) {
    states_.emplace_back();
  }
  states_[static_cast<std::size_t>(state)].push_back({kind, arg, threshold, width, above});
}

void Workspace::stateValues(const Cloud &cloud, std::vector<double> &chi,
                            std::vector<std::vector<double>> *grads) {
  const std::size_t dim = 3 * static_cast<std::size_t>(cloud.nop);
  std::map<std::pair<int, int>, std::pair<double, std::vector<double>>> cache;
  chi.assign(states_.size(), 0.0);
  if (grads != nullptr) {
    grads->assign(states_.size(), std::vector<double>(dim, 0.0));
  }
  std::vector<double> factors;
  std::vector<double> slopes;
  std::vector<double> others;
  for (std::size_t k = 0; k < states_.size(); ++k) {
    const auto &clauses = states_[k];
    factors.clear();
    slopes.clear();
    std::vector<const std::vector<double> *> sources;
    for (const auto &clause : clauses) {
      double v = 0.0;
      const std::vector<double> *source = nullptr;
      if (clause.kind == Kind::State) {
        v = chi[static_cast<std::size_t>(clause.arg)];
        if (grads != nullptr) {
          source = &(*grads)[static_cast<std::size_t>(clause.arg)];
        }
      } else {
        const auto key = std::make_pair(static_cast<int>(clause.kind), clause.arg);
        auto it = cache.find(key);
        if (it == cache.end()) {
          std::vector<double> g;
          const double val =
              channel(cloud, clause.kind, clause.arg, grads != nullptr ? &g : nullptr);
          it = cache.emplace(key, std::make_pair(val, std::move(g))).first;
        }
        v = it->second.first;
        source = grads != nullptr ? &it->second.second : nullptr;
      }
      double f = 0.0;
      double slope = 0.0;
      if (clause.width == 0.0) {
        f = clause.above ? v : 1.0 - v;
        slope = clause.above ? 1.0 : -1.0;
      } else {
        const double sig = logistic((v - clause.threshold) / clause.width);
        const double dsig = sig * (1.0 - sig) / clause.width;
        f = clause.above ? sig : 1.0 - sig;
        slope = clause.above ? dsig : -dsig;
      }
      factors.push_back(f);
      slopes.push_back(slope);
      sources.push_back(source);
    }
    leaveOneOut(factors, others);
    chi[k] = factors.empty() ? 1.0 : others.front() * factors.front();
    if (grads != nullptr) {
      auto &out = (*grads)[k];
      for (std::size_t c = 0; c < factors.size(); ++c) {
        const double scale = others[c] * slopes[c];
        if (scale == 0.0 || sources[c] == nullptr) {
          continue;
        }
        const auto &src = *sources[c];
        for (std::size_t d = 0; d < dim; ++d) {
          out[d] += scale * src[d];
        }
      }
    }
  }
}

double Workspace::value(const Cloud &cloud, Kind kind, int arg, double *grad) {
  checkKind(kind, arg);
  edgeWeights(cloud);
  const std::size_t dim = 3 * static_cast<std::size_t>(cloud.nop);
  if (kind == Kind::State) {
    std::vector<double> chi;
    std::vector<std::vector<double>> grads;
    stateValues(cloud, chi, grad != nullptr ? &grads : nullptr);
    if (grad != nullptr) {
      std::copy(grads[static_cast<std::size_t>(arg)].begin(),
                grads[static_cast<std::size_t>(arg)].end(), grad);
    }
    return chi[static_cast<std::size_t>(arg)];
  }
  std::vector<double> g;
  const double v = channel(cloud, kind, arg, grad != nullptr ? &g : nullptr);
  if (grad != nullptr) {
    std::copy(g.begin(), g.begin() + static_cast<std::ptrdiff_t>(dim), grad);
  }
  return v;
}

void Workspace::states(const Cloud &cloud, double *chi, double *grad) {
  edgeWeights(cloud);
  std::vector<double> values;
  std::vector<std::vector<double>> grads;
  stateValues(cloud, values, grad != nullptr ? &grads : nullptr);
  std::copy(values.begin(), values.end(), chi);
  if (grad != nullptr) {
    const std::size_t dim = 3 * static_cast<std::size_t>(cloud.nop);
    for (std::size_t k = 0; k < grads.size(); ++k) {
      std::copy(grads[k].begin(), grads[k].end(), grad + k * dim);
    }
  }
}

} // namespace smooth
