//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//-----------------------------------------------------------------------------------

#include <science_front.hpp>

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <interface_order.hpp>
#include <local_symmetry.hpp>
#include <neighbours.hpp>
#include <nucleation_pathway.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>

namespace science {
namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

int chillCode(molSys::atom_state_type type) { return static_cast<int>(type); }

bool isCubic(int code) { return code == chillCode(molSys::atom_state_type::cubic); }

bool isHex(int code) { return code == chillCode(molSys::atom_state_type::hexagonal); }

bool isInter(int code) { return code == chillCode(molSys::atom_state_type::interfacial); }

bool isCrystal(int code) { return isCubic(code) || isHex(code); }

bool isExtended(int code) { return isCrystal(code) || isInter(code); }

double axisOf(const molSys::Point<double> &p, int axis) {
  if (axis == 0) {
    return p.x;
  }
  if (axis == 1) {
    return p.y;
  }
  return p.z;
}

double lengthOf(const Cloud &cloud, int axis) {
  if (axis < 0 || static_cast<std::size_t>(axis) >= cloud.box.size()) {
    return 0.0;
  }
  return cloud.box[static_cast<std::size_t>(axis)];
}

double originOf(const Cloud &cloud, int axis) {
  if (axis < 0 || static_cast<std::size_t>(axis) >= cloud.boxLow.size()) {
    return 0.0;
  }
  return cloud.boxLow[static_cast<std::size_t>(axis)];
}

int longAxis(const Cloud &cloud) {
  int axis = 0;
  for (int i = 1; i < 3; ++i) {
    if (lengthOf(cloud, i) > lengthOf(cloud, axis)) {
      axis = i;
    }
  }
  return axis;
}

std::uint64_t xorshift(std::uint64_t &state) {
  state ^= state << 7;
  state ^= state >> 9;
  state ^= state << 8;
  return state;
}

double quantile7(const std::vector<double> &sorted, double q) {
  if (sorted.empty()) {
    return 0.0;
  }
  if (sorted.size() == 1) {
    return sorted.front();
  }
  const double pos = q * static_cast<double>(sorted.size() - 1);
  const int lo = static_cast<int>(pos);
  const double frac = pos - static_cast<double>(lo);
  if (lo + 1 >= static_cast<int>(sorted.size())) {
    return sorted.back();
  }
  return sorted[static_cast<std::size_t>(lo)] * (1.0 - frac) +
         sorted[static_cast<std::size_t>(lo + 1)] * frac;
}

double crossing(double s0, double p0, double s1, double p1, double level) {
  const double t = (level - p0) / (p1 - p0);
  return s0 + t * (s1 - s0);
}

struct Molecule {
  int chill = 2;
  int shell = -1;
  double q6 = 0;
  double q6bar = 0;
};

struct Prepared {
  std::vector<Molecule> mol;
  std::vector<int> offsets;
  std::vector<int> cols;
};

Prepared prepare(Cloud &cloud, double q6Cut) {
  Prepared out;
  out.mol.resize(static_cast<std::size_t>(std::max(cloud.nop, 0)));
  if (cloud.nop <= 0) {
    return out;
  }
  if (!nneigh::orthoFullCsr(cloud, nullptr, kCutoffA, out.offsets, out.cols, false)) {
    out.offsets.assign(static_cast<std::size_t>(cloud.nop + 1), 0);
  }
  std::vector<int> shell;
  symm::classifyIceCloud(cloud, shell);
  const auto q6 = chill::steinhardtQlCutoff(cloud, kCutoffA, 6, 1, true);
  std::vector<int> indexToID(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    indexToID[static_cast<std::size_t>(i)] = cloud.pts[static_cast<std::size_t>(i)].atomID;
  }
  const auto idRows = nneigh::idRowsFromCsr(cloud, indexToID, out.offsets, out.cols);
  chill::applyChillPlus(cloud, idRows, chill::ChillPlusRule::Published, false);
  (void)q6Cut;
  for (int i = 0; i < cloud.nop; ++i) {
    Molecule &m = out.mol[static_cast<std::size_t>(i)];
    m.chill = chillCode(cloud.pts[static_cast<std::size_t>(i)].iceType);
    if (i < static_cast<int>(shell.size())) {
      m.shell = shell[static_cast<std::size_t>(i)];
    }
    if (i < static_cast<int>(q6.ql.size()) && std::isfinite(q6.ql[static_cast<std::size_t>(i)])) {
      m.q6 = q6.ql[static_cast<std::size_t>(i)];
    }
    if (i < static_cast<int>(q6.qlBar.size()) &&
        std::isfinite(q6.qlBar[static_cast<std::size_t>(i)])) {
      m.q6bar = q6.qlBar[static_cast<std::size_t>(i)];
    }
  }
  return out;
}

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

ClusterPart largestPart(int n, const std::vector<int> &offsets, const std::vector<int> &cols,
                        const std::vector<char> &mask, const std::vector<Molecule> &mol) {
  ClusterPart out;
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
  auto unite = [&](int a, int b) {
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
  };
  for (int i = 0; i < n; ++i) {
    if (!mask[static_cast<std::size_t>(i)]) {
      continue;
    }
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    for (int k = a; k < b; ++k) {
      const int j = cols[static_cast<std::size_t>(k)];
      if (j > i && j < n && mask[static_cast<std::size_t>(j)]) {
        unite(i, j);
      }
    }
  }
  int bestRoot = -1;
  int bestSize = 0;
  for (int i = 0; i < n; ++i) {
    if (!mask[static_cast<std::size_t>(i)] || parent[static_cast<std::size_t>(i)] != i) {
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
    if (!mask[static_cast<std::size_t>(i)] || findRoot(parent, i) != bestRoot) {
      continue;
    }
    const Molecule &m = mol[static_cast<std::size_t>(i)];
    out.cubic += isCubic(m.chill) ? 1 : 0;
    out.hexagonal += isHex(m.chill) ? 1 : 0;
    out.interfacial += isInter(m.chill) ? 1 : 0;
    if (m.shell == 1) {
      ++out.shellCubic;
    } else if (m.shell == 0) {
      ++out.shellHex;
    } else if (m.shell == 2) {
      ++out.shellStack;
    } else {
      ++out.shellOther;
    }
  }
  return out;
}

std::vector<char> maskOf(const std::vector<Molecule> &mol, auto keep) {
  std::vector<char> mask(mol.size(), 0);
  for (std::size_t i = 0; i < mol.size(); ++i) {
    mask[i] = keep(mol[i]) ? 1 : 0;
  }
  return mask;
}

void addShell(LayerRow &row, int shell) {
  if (shell == 1) {
    ++row.shellCubic;
  } else if (shell == 0) {
    ++row.shellHex;
  } else if (shell == 2) {
    ++row.shellStack;
  } else {
    ++row.shellOther;
  }
}

} // namespace

Layers fitSolidFraction(const std::vector<double> &phi, const std::vector<double> &mid,
                        double span) {
  Layers fit;
  const int n = static_cast<int>(std::min(phi.size(), mid.size()));
  if (n == 0 || !(span > 0.0)) {
    return fit;
  }
  int imax = 0;
  for (int i = 1; i < n; ++i) {
    if (phi[static_cast<std::size_t>(i)] > phi[static_cast<std::size_t>(imax)]) {
      imax = i;
    }
  }
  if (phi[static_cast<std::size_t>(imax)] < 0.5) {
    return fit;
  }
  double left50 = 0;
  double right50 = 0;
  bool haveL = false;
  bool haveR = false;
  for (int dir = -1; dir <= 1; dir += 2) {
    double prevP = phi[static_cast<std::size_t>(imax)];
    bool have90 = false;
    bool have50 = false;
    bool have10 = false;
    double s90 = 0;
    double s50 = 0;
    double s10 = 0;
    double prevS = 0;
    for (int step = 1; step < n; ++step) {
      const int b = (imax + dir * step + n * 8) % n;
      const double s = step * (span / static_cast<double>(n));
      const double p = phi[static_cast<std::size_t>(b)];
      if (!have90 && prevP >= 0.9 && p < 0.9 && p != prevP) {
        s90 = crossing(prevS, prevP, s, p, 0.9);
        have90 = true;
      }
      if (!have50 && prevP >= 0.5 && p < 0.5 && p != prevP) {
        s50 = crossing(prevS, prevP, s, p, 0.5);
        have50 = true;
      }
      if (!have10 && prevP >= 0.1 && p < 0.1 && p != prevP) {
        s10 = crossing(prevS, prevP, s, p, 0.1);
        have10 = true;
      }
      prevS = s;
      prevP = p;
      if (have50 && have10) {
        break;
      }
    }
    if (!have50) {
      continue;
    }
    Front front;
    double z = mid[static_cast<std::size_t>(imax)] + static_cast<double>(dir) * s50;
    z = std::fmod(z, span);
    if (z < 0.0) {
      z += span;
    }
    front.position = z;
    front.defined = true;
    if (have90 && have10) {
      front.width = std::abs(s10 - s90);
    }
    fit.fronts.push_back(front);
    if (dir < 0) {
      left50 = s50;
      haveL = true;
    } else {
      right50 = s50;
      haveR = true;
    }
  }
  if (haveL && haveR) {
    fit.thickness = left50 + right50;
    fit.defined = true;
  }
  return fit;
}

Layers layersOf(const Cloud &cloud, int axis, double binWidth) {
  Layers out;
  Cloud work = cloud;
  const Prepared prep = prepare(work, kQ6BarIceCut);
  out.n = cloud.nop;
  if (cloud.nop <= 0) {
    return out;
  }
  if (axis < 0 || axis > 2) {
    axis = longAxis(cloud);
  }
  out.axis = axis;
  const double span = lengthOf(cloud, axis);
  out.length = span;
  if (!(span > 0.0)) {
    return out;
  }
  const double dr = binWidth > 0.0 ? binWidth : 1.0;
  const int nb = std::max(8, static_cast<int>(std::lround(span / dr)));
  const double width = span / static_cast<double>(nb);
  out.rows.assign(static_cast<std::size_t>(nb), LayerRow{});
  std::vector<double> q6Sum(static_cast<std::size_t>(nb), 0.0);
  std::vector<double> barSum(static_cast<std::size_t>(nb), 0.0);
  std::vector<int> q6N(static_cast<std::size_t>(nb), 0);
  const double origin = originOf(cloud, axis);
  for (int b = 0; b < nb; ++b) {
    out.rows[static_cast<std::size_t>(b)].mid = (static_cast<double>(b) + 0.5) * width;
  }
  for (int i = 0; i < cloud.nop; ++i) {
    double q = axisOf(cloud.pts[static_cast<std::size_t>(i)], axis) - origin;
    q = std::fmod(q, span);
    if (q < 0.0) {
      q += span;
    }
    int bin = static_cast<int>(q / width);
    if (bin >= nb) {
      bin = nb - 1;
    }
    if (bin < 0) {
      bin = 0;
    }
    LayerRow &row = out.rows[static_cast<std::size_t>(bin)];
    const Molecule &m = prep.mol[static_cast<std::size_t>(i)];
    ++row.n;
    if (isCubic(m.chill)) {
      ++row.cubic;
    } else if (isHex(m.chill)) {
      ++row.hexagonal;
    } else if (isInter(m.chill)) {
      ++row.interfacial;
    } else {
      ++row.liquid;
    }
    addShell(row, m.shell);
    q6Sum[static_cast<std::size_t>(bin)] += m.q6;
    barSum[static_cast<std::size_t>(bin)] += m.q6bar;
    ++q6N[static_cast<std::size_t>(bin)];
  }
  std::vector<double> phi(static_cast<std::size_t>(nb), 0.0);
  std::vector<double> mid(static_cast<std::size_t>(nb), 0.0);
  for (int b = 0; b < nb; ++b) {
    LayerRow &row = out.rows[static_cast<std::size_t>(b)];
    if (q6N[static_cast<std::size_t>(b)] > 0) {
      const double inv = 1.0 / static_cast<double>(q6N[static_cast<std::size_t>(b)]);
      row.q6 = q6Sum[static_cast<std::size_t>(b)] * inv;
      row.q6bar = barSum[static_cast<std::size_t>(b)] * inv;
    }
    const int solid = row.cubic + row.hexagonal + row.interfacial;
    phi[static_cast<std::size_t>(b)] =
        row.n > 0 ? static_cast<double>(solid) / static_cast<double>(row.n) : 0.0;
    mid[static_cast<std::size_t>(b)] = row.mid;
  }
  // A one-angstrom bin between ice planes can be empty. That gap is not
  // liquid. Fill it from the occupied bins on either side so the front
  // is the crystal-liquid crossing.
  for (int b = 0; b < nb; ++b) {
    if (out.rows[static_cast<std::size_t>(b)].n > 0) {
      continue;
    }
    int left = -1;
    int right = -1;
    for (int step = 1; step < nb; ++step) {
      const int j = b - step;
      if (left < 0 && j >= 0 && out.rows[static_cast<std::size_t>(j)].n > 0) {
        left = j;
      }
      const int k = b + step;
      if (right < 0 && k < nb && out.rows[static_cast<std::size_t>(k)].n > 0) {
        right = k;
      }
      if (left >= 0 && right >= 0) {
        break;
      }
    }
    if (left >= 0 && right >= 0) {
      phi[static_cast<std::size_t>(b)] =
          0.5 * (phi[static_cast<std::size_t>(left)] + phi[static_cast<std::size_t>(right)]);
    } else if (left >= 0) {
      phi[static_cast<std::size_t>(b)] = phi[static_cast<std::size_t>(left)];
    } else if (right >= 0) {
      phi[static_cast<std::size_t>(b)] = phi[static_cast<std::size_t>(right)];
    }
  }
  const Layers fitted = fitSolidFraction(phi, mid, span);
  out.fronts = fitted.fronts;
  out.thickness = fitted.thickness;
  out.defined = fitted.defined;
  return out;
}

SlopeInterval bootstrapSlope(const std::vector<double> &x, const std::vector<double> &y,
                             int nDraw, std::uint64_t seed) {
  SlopeInterval out;
  const int n = static_cast<int>(std::min(x.size(), y.size()));
  out.n = n;
  if (n < 2 || nDraw <= 0) {
    return out;
  }
  auto slopeOf = [](const std::vector<double> &xs, const std::vector<double> &ys) {
    const int m = static_cast<int>(xs.size());
    double mx = 0.0;
    double my = 0.0;
    for (int i = 0; i < m; ++i) {
      mx += xs[static_cast<std::size_t>(i)];
      my += ys[static_cast<std::size_t>(i)];
    }
    mx /= static_cast<double>(m);
    my /= static_cast<double>(m);
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < m; ++i) {
      const double dx = xs[static_cast<std::size_t>(i)] - mx;
      num += dx * (ys[static_cast<std::size_t>(i)] - my);
      den += dx * dx;
    }
    if (!(den > 0.0)) {
      return std::pair<double, bool>{0.0, false};
    }
    return std::pair<double, bool>{num / den, true};
  };
  std::vector<double> xs(x.begin(), x.begin() + n);
  std::vector<double> ys(y.begin(), y.begin() + n);
  const auto point = slopeOf(xs, ys);
  if (!point.second) {
    return out;
  }
  out.slope = point.first;
  std::uint64_t state = seed == 0 ? 1 : seed;
  std::vector<double> draws;
  draws.reserve(static_cast<std::size_t>(nDraw));
  std::vector<double> rx(static_cast<std::size_t>(n));
  std::vector<double> ry(static_cast<std::size_t>(n));
  for (int draw = 0; draw < nDraw; ++draw) {
    for (int i = 0; i < n; ++i) {
      const auto pick = static_cast<std::size_t>(xorshift(state) % static_cast<std::uint64_t>(n));
      rx[static_cast<std::size_t>(i)] = xs[pick];
      ry[static_cast<std::size_t>(i)] = ys[pick];
    }
    const auto one = slopeOf(rx, ry);
    if (one.second) {
      draws.push_back(one.first);
    }
  }
  if (draws.empty()) {
    out.lo = out.slope;
    out.hi = out.slope;
    return out;
  }
  std::sort(draws.begin(), draws.end());
  out.lo = quantile7(draws, 0.025);
  out.hi = quantile7(draws, 0.975);
  return out;
}

CriticalSize criticalSize(const std::vector<double> &series, int nDraw, std::uint64_t seed) {
  CriticalSize out;
  if (series.empty()) {
    return out;
  }
  // One spontaneous dump does not measure a committor. The estimate is
  // the mean size in the contiguous band around the smallest largest
  // cluster. The band stops where the size exceeds that minimum by
  // the larger of 5 molecules and 15 percent. A minimum in the first
  // two frames, followed by a larger cluster, is an upper bound: the
  // dump starts after the nucleus is already at least this large.
  int imin = 0;
  for (int i = 1; i < static_cast<int>(series.size()); ++i) {
    if (series[static_cast<std::size_t>(i)] < series[static_cast<std::size_t>(imin)]) {
      imin = i;
    }
  }
  const double floor = series[static_cast<std::size_t>(imin)];
  const double cap = floor + std::max(5.0, 0.15 * std::abs(floor));
  int lo = imin;
  while (lo > 0 && series[static_cast<std::size_t>(lo - 1)] <= cap) {
    --lo;
  }
  int hi = imin;
  while (hi + 1 < static_cast<int>(series.size()) &&
         series[static_cast<std::size_t>(hi + 1)] <= cap) {
    ++hi;
  }
  std::vector<double> window;
  window.reserve(static_cast<std::size_t>(hi - lo + 1));
  for (int i = lo; i <= hi; ++i) {
    window.push_back(series[static_cast<std::size_t>(i)]);
  }
  out.defined = true;
  out.onset = imin;
  out.window = static_cast<int>(window.size());
  double sum = 0.0;
  for (double v : window) {
    sum += v;
  }
  out.nStar = sum / static_cast<double>(window.size());
  double var = 0.0;
  for (double v : window) {
    const double d = v - out.nStar;
    var += d * d;
  }
  if (window.size() > 1) {
    out.sigma = std::sqrt(var / static_cast<double>(window.size() - 1));
  }
  const double later = *std::max_element(series.begin(), series.end());
  out.upperBound = imin <= 1 && later > out.nStar + out.sigma;
  const iface::BootstrapInterval band = iface::bootstrapMean(window, nDraw, seed);
  out.lo = band.lo;
  out.hi = band.hi;
  if (window.size() == 1) {
    out.lo = out.nStar;
    out.hi = out.nStar;
  }
  return out;
}

NucleusFrame nucleusOf(Cloud &cloud, int frame, double q6Cut) {
  NucleusFrame rec;
  rec.frame = frame;
  rec.nop = cloud.nop;
  const Prepared prep = prepare(cloud, q6Cut);
  if (cloud.nop <= 0) {
    return rec;
  }
  double sum = 0.0;
  int nFinite = 0;
  for (const Molecule &m : prep.mol) {
    if (isCubic(m.chill)) {
      ++rec.chillCubic;
    } else if (isHex(m.chill)) {
      ++rec.chillHex;
    } else if (isInter(m.chill)) {
      ++rec.chillInter;
    } else {
      ++rec.chillLiquid;
    }
    if (m.shell == 1) {
      ++rec.shellCubic;
    } else if (m.shell == 0) {
      ++rec.shellHex;
    } else if (m.shell == 2) {
      ++rec.shellStack;
      if (isCubic(m.chill)) {
        ++rec.stackAsCubic;
      } else if (isHex(m.chill)) {
        ++rec.stackAsHex;
      } else if (isInter(m.chill)) {
        ++rec.stackAsInter;
      } else {
        ++rec.stackAsLiquid;
      }
    } else {
      ++rec.shellOther;
    }
    if (m.q6bar > q6Cut) {
      ++rec.q6Ice;
    }
    if (std::isfinite(m.q6bar)) {
      sum += m.q6bar;
      ++nFinite;
    }
    if ((m.shell == 0 || m.shell == 1) && isCrystal(m.chill)) {
      ++rec.polytypeBoth;
      const bool agree = (m.shell == 0 && isHex(m.chill)) || (m.shell == 1 && isCubic(m.chill));
      rec.polytypeAgree += agree ? 1 : 0;
    }
  }
  rec.q6barMean = nFinite > 0 ? sum / static_cast<double>(nFinite) : 0.0;
  const int n = cloud.nop;
  rec.core = largestPart(n, prep.offsets, prep.cols,
                         maskOf(prep.mol, [](const Molecule &m) { return isCrystal(m.chill); }),
                         prep.mol);
  rec.extended =
      largestPart(n, prep.offsets, prep.cols,
                  maskOf(prep.mol, [](const Molecule &m) { return isExtended(m.chill); }), prep.mol);
  rec.shell = largestPart(
      n, prep.offsets, prep.cols,
      maskOf(prep.mol, [](const Molecule &m) { return m.shell == 0 || m.shell == 1 || m.shell == 2; }),
      prep.mol);
  rec.q6 = largestPart(n, prep.offsets, prep.cols,
                       maskOf(prep.mol, [q6Cut](const Molecule &m) { return m.q6bar > q6Cut; }),
                       prep.mol);
  return rec;
}

std::vector<AtomRow> labelAtoms(Cloud cloud, double q6Cut) {
  std::vector<AtomRow> rows;
  if (cloud.nop <= 0) {
    return rows;
  }
  Cloud ovitoCloud = cloud;
  const Prepared published = prepare(cloud, q6Cut);
  std::vector<int> indexToID(static_cast<std::size_t>(ovitoCloud.nop));
  for (int i = 0; i < ovitoCloud.nop; ++i) {
    indexToID[static_cast<std::size_t>(i)] = ovitoCloud.pts[static_cast<std::size_t>(i)].atomID;
  }
  std::vector<int> offsets;
  std::vector<int> cols;
  nneigh::orthoFullCsr(ovitoCloud, nullptr, kCutoffA, offsets, cols, false);
  const auto idRows = nneigh::idRowsFromCsr(ovitoCloud, indexToID, offsets, cols);
  chill::applyChillPlus(ovitoCloud, idRows, chill::ChillPlusRule::Ovito, false);

  const auto extendedMask =
      maskOf(published.mol, [](const Molecule &m) { return isExtended(m.chill); });
  const ClusterPart part =
      largestPart(cloud.nop, published.offsets, published.cols, extendedMask, published.mol);
  std::vector<char> in(static_cast<std::size_t>(cloud.nop), 0);
  if (part.size > 0) {
    // Recompute membership. largestPart does not return the root, so
    // mark every extended molecule and keep those in a component of
    // that size whose root is the smallest such root.
    std::vector<int> parent(static_cast<std::size_t>(cloud.nop));
    std::vector<int> sz(static_cast<std::size_t>(cloud.nop), 0);
    std::iota(parent.begin(), parent.end(), 0);
    for (int i = 0; i < cloud.nop; ++i) {
      if (extendedMask[static_cast<std::size_t>(i)]) {
        sz[static_cast<std::size_t>(i)] = 1;
      }
    }
    auto unite = [&](int a, int b) {
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
    };
    for (int i = 0; i < cloud.nop; ++i) {
      if (!extendedMask[static_cast<std::size_t>(i)]) {
        continue;
      }
      const int a = published.offsets[static_cast<std::size_t>(i)];
      const int b = published.offsets[static_cast<std::size_t>(i) + 1];
      for (int k = a; k < b; ++k) {
        const int j = published.cols[static_cast<std::size_t>(k)];
        if (j > i && j < cloud.nop && extendedMask[static_cast<std::size_t>(j)]) {
          unite(i, j);
        }
      }
    }
    int bestRoot = -1;
    for (int i = 0; i < cloud.nop; ++i) {
      if (!extendedMask[static_cast<std::size_t>(i)] || parent[static_cast<std::size_t>(i)] != i) {
        continue;
      }
      if (sz[static_cast<std::size_t>(i)] == part.size && (bestRoot < 0 || i < bestRoot)) {
        bestRoot = i;
      }
    }
    if (bestRoot >= 0) {
      for (int i = 0; i < cloud.nop; ++i) {
        if (extendedMask[static_cast<std::size_t>(i)] && findRoot(parent, i) == bestRoot) {
          in[static_cast<std::size_t>(i)] = 1;
        }
      }
    }
  }
  rows.resize(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    AtomRow &row = rows[static_cast<std::size_t>(i)];
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    const Molecule &m = published.mol[static_cast<std::size_t>(i)];
    row.id = p.atomID;
    row.x = p.x;
    row.y = p.y;
    row.z = p.z;
    row.chill = m.chill;
    row.chillOvito = chillCode(ovitoCloud.pts[static_cast<std::size_t>(i)].iceType);
    row.shell = m.shell;
    row.q6 = m.q6;
    row.q6bar = m.q6bar;
    row.inExtended = in[static_cast<std::size_t>(i)];
  }
  return rows;
}

} // namespace science
