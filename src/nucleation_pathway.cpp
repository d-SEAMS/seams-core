//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <nucleation_pathway.hpp>

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <local_symmetry.hpp>
#include <franzblau.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace nucleation {
namespace {

struct Tag {
  bool cubic = false;
  bool hex = false;
};

struct Component {
  int size = 0;
  int cubic = 0;
  int hex = 0;
};

int findRoot(std::vector<int> &parent, int i) {
  int root = i;
  while (parent[static_cast<std::size_t>(root)] != root) {
    root = parent[static_cast<std::size_t>(root)];
  }
  while (parent[static_cast<std::size_t>(i)] != root) {
    const int next = parent[static_cast<std::size_t>(i)];
    parent[static_cast<std::size_t>(i)] = root;
    i = next;
  }
  return root;
}

void unite(std::vector<int> &parent, std::vector<int> &sz, int a, int b) {
  a = findRoot(parent, a);
  b = findRoot(parent, b);
  if (a == b) {
    return;
  }
  if (sz[static_cast<std::size_t>(a)] < sz[static_cast<std::size_t>(b)]) {
    std::swap(a, b);
  }
  parent[static_cast<std::size_t>(b)] = a;
  sz[static_cast<std::size_t>(a)] += sz[static_cast<std::size_t>(b)];
}

Component largestComponent(int n, const std::vector<int> &offsets,
                           const std::vector<int> &cols,
                           const std::vector<char> &mask,
                           const std::vector<Tag> &tag) {
  Component out;
  if (n <= 0) {
    return out;
  }
  std::vector<int> parent(static_cast<std::size_t>(n));
  std::vector<int> sz(static_cast<std::size_t>(n), 0);
  std::iota(parent.begin(), parent.end(), 0);
  for (int i = 0; i < n; ++i) {
    if (mask[static_cast<std::size_t>(i)]) {
      sz[static_cast<std::size_t>(i)] = 1;
    }
  }
  for (int i = 0; i < n; ++i) {
    if (!mask[static_cast<std::size_t>(i)]) {
      continue;
    }
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    for (int k = a; k < b; ++k) {
      const int j = cols[static_cast<std::size_t>(k)];
      if (j > i && j < n && mask[static_cast<std::size_t>(j)]) {
        unite(parent, sz, i, j);
      }
    }
  }
  int bestRoot = -1;
  int bestSize = 0;
  for (int i = 0; i < n; ++i) {
    if (!mask[static_cast<std::size_t>(i)]) {
      continue;
    }
    if (parent[static_cast<std::size_t>(i)] != i) {
      continue;
    }
    const int s = sz[static_cast<std::size_t>(i)];
    if (s > bestSize || (s == bestSize && (bestRoot < 0 || i < bestRoot))) {
      bestSize = s;
      bestRoot = i;
    }
  }
  if (bestRoot < 0) {
    return out;
  }
  out.size = bestSize;
  for (int i = 0; i < n; ++i) {
    if (!mask[static_cast<std::size_t>(i)]) {
      continue;
    }
    if (findRoot(parent, i) != bestRoot) {
      continue;
    }
    out.cubic += tag[static_cast<std::size_t>(i)].cubic ? 1 : 0;
    out.hex += tag[static_cast<std::size_t>(i)].hex ? 1 : 0;
  }
  return out;
}

void tallyChill(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
                int &cubic, int &hex, int &inter, int &clath, int &interClath,
                int &water) {
  cubic = hex = inter = clath = interClath = water = 0;
  for (int i = 0; i < cloud.nop; ++i) {
    switch (cloud.pts[static_cast<std::size_t>(i)].iceType) {
    case molSys::atom_state_type::cubic:
      ++cubic;
      break;
    case molSys::atom_state_type::hexagonal:
      ++hex;
      break;
    case molSys::atom_state_type::interfacial:
      ++inter;
      break;
    case molSys::atom_state_type::clathrate:
      ++clath;
      break;
    case molSys::atom_state_type::interClathrate:
      ++interClath;
      break;
    case molSys::atom_state_type::water:
    case molSys::atom_state_type::unclassified:
    case molSys::atom_state_type::reCubic:
    case molSys::atom_state_type::reHex:
      ++water;
      break;
    }
  }
}

std::vector<Tag> tagsFromChill(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  std::vector<Tag> tag(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    const auto kind = cloud.pts[static_cast<std::size_t>(i)].iceType;
    tag[static_cast<std::size_t>(i)].cubic =
        kind == molSys::atom_state_type::cubic;
    tag[static_cast<std::size_t>(i)].hex =
        kind == molSys::atom_state_type::hexagonal;
  }
  return tag;
}

std::vector<char> maskCore(const std::vector<Tag> &tag) {
  std::vector<char> mask(tag.size(), 0);
  for (std::size_t i = 0; i < tag.size(); ++i) {
    mask[i] = (tag[i].cubic || tag[i].hex) ? 1 : 0;
  }
  return mask;
}

std::vector<char> maskExtended(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  std::vector<char> mask(static_cast<std::size_t>(cloud.nop), 0);
  for (int i = 0; i < cloud.nop; ++i) {
    const auto kind = cloud.pts[static_cast<std::size_t>(i)].iceType;
    mask[static_cast<std::size_t>(i)] =
        (kind == molSys::atom_state_type::cubic ||
         kind == molSys::atom_state_type::hexagonal ||
         kind == molSys::atom_state_type::interfacial)
            ? 1
            : 0;
  }
  return mask;
}

void storeComponent(int &n, int &c, int &h, const Component &part) {
  n = part.size;
  c = part.cubic;
  h = part.hex;
}

double medianOf(std::vector<double> values) {
  if (values.empty()) {
    return 0;
  }
  const std::size_t mid = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid),
                   values.end());
  if (values.size() % 2 == 1) {
    return values[mid];
  }
  const double upper = values[mid];
  std::nth_element(values.begin(),
                   values.begin() + static_cast<std::ptrdiff_t>(mid - 1),
                   values.end());
  return 0.5 * (values[mid - 1] + upper);
}

double pearson(const std::vector<double> &x, const std::vector<double> &y) {
  const int n = static_cast<int>(x.size());
  if (n < 5 || static_cast<int>(y.size()) != n) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  double mx = 0;
  double my = 0;
  for (int i = 0; i < n; ++i) {
    mx += x[static_cast<std::size_t>(i)];
    my += y[static_cast<std::size_t>(i)];
  }
  mx /= n;
  my /= n;
  double num = 0;
  double dx = 0;
  double dy = 0;
  for (int i = 0; i < n; ++i) {
    const double a = x[static_cast<std::size_t>(i)] - mx;
    const double b = y[static_cast<std::size_t>(i)] - my;
    num += a * b;
    dx += a * a;
    dy += b * b;
  }
  if (!(dx > 0) || !(dy > 0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return num / std::sqrt(dx * dy);
}

std::vector<double> firstDifference(const std::vector<double> &series) {
  std::vector<double> out;
  if (series.size() < 2) {
    return out;
  }
  out.resize(series.size() - 1);
  for (std::size_t i = 0; i + 1 < series.size(); ++i) {
    out[i] = series[i + 1] - series[i];
  }
  return out;
}

} // namespace

double stackFraction(int cubic, int hexagonal) {
  const int den = cubic + hexagonal;
  if (den <= 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return static_cast<double>(cubic) / static_cast<double>(den);
}

HalfRise halfRise(const std::vector<double> &series) {
  HalfRise out;
  if (series.empty()) {
    return out;
  }
  int imax = 0;
  for (int i = 1; i < static_cast<int>(series.size()); ++i) {
    if (series[static_cast<std::size_t>(i)] >
        series[static_cast<std::size_t>(imax)]) {
      imax = i;
    }
  }
  out.start = series[0];
  out.maximum = series[static_cast<std::size_t>(imax)];
  if (!(out.maximum > out.start)) {
    return out;
  }
  out.level = out.start + 0.5 * (out.maximum - out.start);
  for (int i = 0; i < static_cast<int>(series.size()); ++i) {
    if (series[static_cast<std::size_t>(i)] >= out.level) {
      out.ok = true;
      out.index = i;
      return out;
    }
  }
  return out;
}

Ascent finalAscent(const std::vector<double> &series) {
  Ascent out;
  if (series.empty()) {
    out.index = -1;
    return out;
  }
  int imax = 0;
  for (int i = 1; i < static_cast<int>(series.size()); ++i) {
    if (series[static_cast<std::size_t>(i)] >
        series[static_cast<std::size_t>(imax)]) {
      imax = i;
    }
  }
  int i = imax;
  while (i > 0 && series[static_cast<std::size_t>(i - 1)] <
                      series[static_cast<std::size_t>(i)]) {
    --i;
  }
  out.index = i;
  out.maxIndex = imax;
  out.baseline = series[static_cast<std::size_t>(i)];
  out.maximum = series[static_cast<std::size_t>(imax)];
  out.climbed = out.maximum > out.baseline;
  return out;
}

GrowBasin growBasin(const std::vector<double> &series) {
  GrowBasin out;
  const Ascent climb = finalAscent(series);
  if (!climb.climbed) {
    return out;
  }
  out.occupied = true;
  out.barrier = climb.baseline;
  out.index = climb.index;
  return out;
}

CommittorCall committorOf(const std::vector<double> &series, int molecules) {
  CommittorCall out;
  out.basin = growBasin(series);
  if (series.empty() || molecules <= 0) {
    return out;
  }
  const double cap = 2.0 * static_cast<double>(molecules);
  double maximum = series[0];
  for (double value : series) {
    if (value > maximum) {
      maximum = value;
    }
  }
  // A series that stops on 2N, and never goes past it, was cut by the
  // old wall. It is not evidence that the nucleus committed.
  out.cutAtTwiceN = series.back() == cap && maximum == cap;
  if (out.cutAtTwiceN) {
    out.fate = CommittorFate::truncated;
    return out;
  }
  // Ending back on the foot has given the climb back. The replica is
  // committed only while the final size is still strictly above it.
  if (out.basin.occupied && series.back() > out.basin.barrier) {
    out.fate = CommittorFate::committed;
  }
  return out;
}

Lead ringLead(const std::vector<double> &ringFraction,
              const std::vector<double> &iceFraction) {
  Lead out;
  if (ringFraction.size() < 8 || ringFraction.size() != iceFraction.size()) {
    return out;
  }
  const Ascent ring = finalAscent(ringFraction);
  const Ascent ice = finalAscent(iceFraction);
  const HalfRise ringHalf = halfRise(ringFraction);
  const HalfRise iceHalf = halfRise(iceFraction);
  if (!ringHalf.ok || !iceHalf.ok) {
    return out;
  }
  out.ok = true;
  out.halfFrames = iceHalf.index - ringHalf.index;
  if (ring.climbed && ice.climbed) {
    out.ascentFrames = ice.index - ring.index;
  }
  const auto dRing = firstDifference(ringFraction);
  const auto dIce = firstDifference(iceFraction);
  const int n = static_cast<int>(dRing.size());
  // A lag is kept only when at least five differences overlap, and the
  // search stops at 20 frames.
  const int maxLag = std::max(0, std::min(20, n - 5));
  double best = -2;
  int bestLag = 0;
  for (int lag = -maxLag; lag <= maxLag; ++lag) {
    const int t0 = std::max(0, -lag);
    const int t1 = std::min(n, n - lag);
    if (t1 - t0 < 5) {
      continue;
    }
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(static_cast<std::size_t>(t1 - t0));
    y.reserve(static_cast<std::size_t>(t1 - t0));
    for (int t = t0; t < t1; ++t) {
      x.push_back(dRing[static_cast<std::size_t>(t)]);
      y.push_back(dIce[static_cast<std::size_t>(t + lag)]);
    }
    const double c = pearson(x, y);
    if (!std::isfinite(c)) {
      continue;
    }
    const bool better = c > best + 1e-12;
    const bool tie = std::abs(c - best) <= 1e-12 &&
                     (std::abs(lag) < std::abs(bestLag) ||
                      (std::abs(lag) == std::abs(bestLag) && lag > bestLag));
    if (better || (std::isfinite(best) && tie) || best < -1) {
      best = c;
      bestLag = lag;
    }
  }
  if (best < -1) {
    out.corrLag = 0;
    out.corr = std::numeric_limits<double>::quiet_NaN();
    return out;
  }
  out.corrLag = bestLag;
  out.corr = best;
  return out;
}

Shift definitionShift(int sizeA, int frameA, int sizeB, int frameB) {
  Shift out;
  out.dSize = sizeB - sizeA;
  out.dFrame = frameB - frameA;
  return out;
}

std::vector<double>
q6DotScore(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
           const std::vector<int> &offsets, const std::vector<int> &cols) {
  constexpr int kL = 6;
  const int n = cloud.nop;
  // m = 0..6 only. Y_6^{-m} = (-1)^m conj(Y_6^m), so the negative orders
  // double the real part of each positive-order product.
  std::vector<std::complex<double>> qlm(static_cast<std::size_t>(n) * (kL + 1));
  std::vector<double> norm(static_cast<std::size_t>(n), 0);
  const double lx = cloud.box[0];
  const double ly = cloud.box[1];
  const double lz = cloud.box[2];
  for (int i = 0; i < n; ++i) {
    const auto &pi = cloud.pts[static_cast<std::size_t>(i)];
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    std::complex<double> *q = &qlm[static_cast<std::size_t>(i) * (kL + 1)];
    for (int k = a; k < b; ++k) {
      const auto &pj = cloud.pts[static_cast<std::size_t>(cols[static_cast<std::size_t>(k)])];
      double dx = pj.x - pi.x;
      double dy = pj.y - pi.y;
      double dz = pj.z - pi.z;
      dx -= lx * std::round(dx / lx);
      dy -= ly * std::round(dy / ly);
      dz -= lz * std::round(dz / lz);
      const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (r <= 0) {
        continue;
      }
      const double theta = std::acos(std::clamp(dz / r, -1.0, 1.0));
      const double phi = std::atan2(dy, dx);
      for (int m = 0; m <= kL; ++m) {
        const double y = std::sph_legendre(kL, m, theta);
        q[m] += std::complex<double>(y * std::cos(m * phi), y * std::sin(m * phi));
      }
    }
    double s = std::norm(q[0]);
    for (int m = 1; m <= kL; ++m) {
      s += 2.0 * std::norm(q[m]);
    }
    norm[static_cast<std::size_t>(i)] = std::sqrt(s);
  }
  std::vector<double> score(static_cast<std::size_t>(n), 0);
  for (int i = 0; i < n; ++i) {
    const double ni = norm[static_cast<std::size_t>(i)];
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    if (ni <= 0 || b <= a) {
      continue;
    }
    const std::complex<double> *qi = &qlm[static_cast<std::size_t>(i) * (kL + 1)];
    double sum = 0;
    for (int k = a; k < b; ++k) {
      const int j = cols[static_cast<std::size_t>(k)];
      const double nj = norm[static_cast<std::size_t>(j)];
      if (nj <= 0) {
        continue;
      }
      const std::complex<double> *qj = &qlm[static_cast<std::size_t>(j) * (kL + 1)];
      double dot = (qi[0] * std::conj(qj[0])).real();
      for (int m = 1; m <= kL; ++m) {
        dot += 2.0 * (qi[m] * std::conj(qj[m])).real();
      }
      sum += dot / (ni * nj);
    }
    score[static_cast<std::size_t>(i)] = sum / static_cast<double>(b - a);
  }
  return score;
}

FrameRecord analyzeFrame(molSys::PointCloud<molSys::Point<double>, double> &cloud,
                         double cutoff, int ringDepth, int typeI, double q6Cut,
                         int frame) {
  FrameRecord rec;
  rec.frame = frame;
  rec.nop = cloud.nop;
  if (cloud.nop <= 0) {
    return rec;
  }
  if (ringDepth < 3 || ringDepth > 7) {
    throw std::invalid_argument("ring depth must be from 3 to 7");
  }
  std::vector<int> offsets;
  std::vector<int> cols;
  if (!nneigh::orthoFullCsr(cloud, nullptr, cutoff, offsets, cols, false)) {
    throw std::runtime_error("cutoff neighbour list failed");
  }
  const auto rings = primitive::ringNetwork(offsets, cols, ringDepth);
  for (const auto &ring : rings) {
    const int sz = static_cast<int>(ring.size());
    if (sz >= 0 && sz <= ringDepth) {
      rec.ringCount[sz] += 1;
    }
  }

  std::vector<int> kind;
  symm::classifyIceCloud(cloud, kind);
  std::vector<Tag> symTag(static_cast<std::size_t>(cloud.nop));
  std::vector<char> symMask(static_cast<std::size_t>(cloud.nop), 0);
  for (int i = 0; i < cloud.nop; ++i) {
    const int k = kind[static_cast<std::size_t>(i)];
    if (k == 1) {
      ++rec.symCubic;
      symTag[static_cast<std::size_t>(i)].cubic = true;
      symMask[static_cast<std::size_t>(i)] = 1;
    } else if (k == 0) {
      ++rec.symHex;
      symTag[static_cast<std::size_t>(i)].hex = true;
      symMask[static_cast<std::size_t>(i)] = 1;
    }
  }
  storeComponent(rec.nSym, rec.cSym, rec.hSym,
                 largestComponent(cloud.nop, offsets, cols, symMask, symTag));

  const auto q6 = chill::steinhardtQlCutoff(cloud, cutoff, 6, typeI, true);
  std::vector<double> finite;
  finite.reserve(static_cast<std::size_t>(cloud.nop));
  double sum = 0;
  int nFinite = 0;
  std::vector<char> qMask(static_cast<std::size_t>(cloud.nop), 0);
  for (int i = 0; i < cloud.nop; ++i) {
    const double q = q6.qlBar[static_cast<std::size_t>(i)];
    if (!std::isfinite(q)) {
      continue;
    }
    sum += q;
    ++nFinite;
    finite.push_back(q);
    if (q > q6Cut) {
      ++rec.q6barIce;
      qMask[static_cast<std::size_t>(i)] = 1;
    }
  }
  rec.q6barMean = nFinite > 0 ? sum / static_cast<double>(nFinite) : 0;
  rec.q6barP50 = medianOf(std::move(finite));
  storeComponent(rec.nQ6, rec.cQ6, rec.hQ6,
                 largestComponent(cloud.nop, offsets, cols, qMask, symTag));

  {
    std::vector<int> solidOffsets;
    std::vector<int> solidCols;
    if (!nneigh::orthoFullCsr(cloud, nullptr, kQ6SolidCutoffA, solidOffsets, solidCols,
                              false)) {
      throw std::runtime_error("q6 solid neighbour list failed");
    }
    const auto score = q6DotScore(cloud, solidOffsets, solidCols);
    std::vector<char> solid(static_cast<std::size_t>(cloud.nop), 0);
    for (int i = 0; i < cloud.nop; ++i) {
      solid[static_cast<std::size_t>(i)] = score[static_cast<std::size_t>(i)] > kQ6SolidCut;
    }
    storeComponent(rec.nQ6Solid, rec.cQ6Solid, rec.hQ6Solid,
                   largestComponent(cloud.nop, solidOffsets, solidCols, solid, symTag));
  }

  std::vector<int> indexToID(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    indexToID[static_cast<std::size_t>(i)] =
        cloud.pts[static_cast<std::size_t>(i)].atomID;
  }
  const auto idRows = nneigh::idRowsFromCsr(cloud, indexToID, offsets, cols);

  chill::applyChillPlus(cloud, idRows, chill::ChillPlusRule::Ovito, false);
  tallyChill(cloud, rec.ovitoCubic, rec.ovitoHex, rec.ovitoInter, rec.ovitoClath,
             rec.ovitoInterClath, rec.ovitoWater);
  {
    const auto tag = tagsFromChill(cloud);
    storeComponent(rec.nOvitoCore, rec.cOvitoCore, rec.hOvitoCore,
                   largestComponent(cloud.nop, offsets, cols, maskCore(tag), tag));
    storeComponent(rec.nOvitoExt, rec.cOvitoExt, rec.hOvitoExt,
                   largestComponent(cloud.nop, offsets, cols, maskExtended(cloud),
                                    tag));
  }

  chill::applyChillPlus(cloud, idRows, chill::ChillPlusRule::Published, false);
  tallyChill(cloud, rec.chillCubic, rec.chillHex, rec.chillInter, rec.chillClath,
             rec.chillInterClath, rec.chillWater);
  {
    const auto tag = tagsFromChill(cloud);
    storeComponent(rec.nChillCore, rec.cChillCore, rec.hChillCore,
                   largestComponent(cloud.nop, offsets, cols, maskCore(tag), tag));
    storeComponent(rec.nChillExt, rec.cChillExt, rec.hChillExt,
                   largestComponent(cloud.nop, offsets, cols, maskExtended(cloud),
                                    tag));
  }
  return rec;
}

} // namespace nucleation
