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

#include <brine_front.hpp>

#include <franzblau.hpp>
#include <generic.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <unordered_set>
#include <utility>

namespace site {
namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

double axisOf(const molSys::Point<double> &p, int axis) {
  if (axis == 0) {
    return p.x;
  }
  if (axis == 1) {
    return p.y;
  }
  return p.z;
}

double boxLength(const Cloud &cloud, int axis) {
  if (axis < 0 || static_cast<std::size_t>(axis) >= cloud.box.size()) {
    return 0;
  }
  return cloud.box[static_cast<std::size_t>(axis)];
}

double boxOrigin(const Cloud &cloud, int axis) {
  if (axis < 0 || static_cast<std::size_t>(axis) >= cloud.boxLow.size()) {
    return 0;
  }
  return cloud.boxLow[static_cast<std::size_t>(axis)];
}

double fold(double z, double origin, double length) {
  if (length <= 0) {
    return z;
  }
  double u = (z - origin) / length;
  u -= std::floor(u);
  if (u < 0) {
    u += 1;
  }
  if (u >= 1) {
    u = 0;
  }
  return origin + u * length;
}

double molalityOf(int nIon, int nWater, int ionsPerFormula) {
  if (nWater <= 0 || ionsPerFormula <= 0) {
    return 0;
  }
  const double formulas = static_cast<double>(nIon) / static_cast<double>(ionsPerFormula);
  const double kg = static_cast<double>(nWater) * kWaterMolarMassKg;
  return formulas / kg;
}

struct Cells {
  double cell = 1;
  int nx = 1;
  int ny = 1;
  int nz = 1;
  std::vector<std::vector<int>> buckets;
};

Cells makeCells(const Cloud &cloud, const std::vector<int> &members, double cutoff) {
  Cells cells;
  const double lx = cloud.box.size() > 0 ? cloud.box[0] : 1;
  const double ly = cloud.box.size() > 1 ? cloud.box[1] : 1;
  const double lz = cloud.box.size() > 2 ? cloud.box[2] : 1;
  cells.cell = cutoff > 0 ? cutoff : 1;
  cells.nx = std::max(1, static_cast<int>(std::floor(lx / cells.cell)));
  cells.ny = std::max(1, static_cast<int>(std::floor(ly / cells.cell)));
  cells.nz = std::max(1, static_cast<int>(std::floor(lz / cells.cell)));
  cells.cell = cutoff;
  cells.buckets.resize(static_cast<std::size_t>(cells.nx * cells.ny * cells.nz));
  auto bin = [&](double value, double origin, double length, int n) {
    if (length <= 0 || n <= 1) {
      return 0;
    }
    int b = static_cast<int>(std::floor((fold(value, origin, length) - origin) / length * n));
    if (b < 0) {
      b = 0;
    }
    if (b >= n) {
      b = n - 1;
    }
    return b;
  };
  const double x0 = boxOrigin(cloud, 0);
  const double y0 = boxOrigin(cloud, 1);
  const double z0 = boxOrigin(cloud, 2);
  for (int index : members) {
    const auto &p = cloud.pts[static_cast<std::size_t>(index)];
    const int ix = bin(p.x, x0, lx, cells.nx);
    const int iy = bin(p.y, y0, ly, cells.ny);
    const int iz = bin(p.z, z0, lz, cells.nz);
    const int key = (ix * cells.ny + iy) * cells.nz + iz;
    cells.buckets[static_cast<std::size_t>(key)].push_back(index);
  }
  return cells;
}

template <typename Visit>
void visitNearby(const Cells &cells, int index, const Cloud &cloud, Visit visit) {
  const double lx = cloud.box.size() > 0 ? cloud.box[0] : 1;
  const double ly = cloud.box.size() > 1 ? cloud.box[1] : 1;
  const double lz = cloud.box.size() > 2 ? cloud.box[2] : 1;
  const double x0 = boxOrigin(cloud, 0);
  const double y0 = boxOrigin(cloud, 1);
  const double z0 = boxOrigin(cloud, 2);
  const auto &p = cloud.pts[static_cast<std::size_t>(index)];
  auto bin = [](double value, double origin, double length, int n) {
    if (length <= 0 || n <= 1) {
      return 0;
    }
    int b = static_cast<int>(std::floor((fold(value, origin, length) - origin) / length * n));
    if (b < 0) {
      b = 0;
    }
    if (b >= n) {
      b = n - 1;
    }
    return b;
  };
  const int cx = bin(p.x, x0, lx, cells.nx);
  const int cy = bin(p.y, y0, ly, cells.ny);
  const int cz = bin(p.z, z0, lz, cells.nz);
  // A repeated cell only happens when an axis has fewer than three bins.
  // Visiting it twice is safe: the caller rejects a self pair.
  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dz = -1; dz <= 1; ++dz) {
        const int ix = (cx + dx + cells.nx) % cells.nx;
        const int iy = (cy + dy + cells.ny) % cells.ny;
        const int iz = (cz + dz + cells.nz) % cells.nz;
        const int key = (ix * cells.ny + iy) * cells.nz + iz;
        for (int other : cells.buckets[static_cast<std::size_t>(key)]) {
          visit(other);
        }
      }
    }
  }
}

} // namespace

ConcentrationProfile concentrationProfile(const Cloud &cloud, const std::vector<bool> &iceFlag,
                                          const std::vector<int> &waterIndices,
                                          const std::vector<int> &ionIndices, int nBins, int axis,
                                          int ionsPerFormula) {
  ConcentrationProfile profile;
  if (nBins < 1) {
    nBins = 1;
  }
  const double origin = boxOrigin(cloud, axis);
  const double length = boxLength(cloud, axis);
  profile.zOrigin = origin;
  profile.length = length;
  profile.binWidth = length > 0 ? length / static_cast<double>(nBins) : 0;
  profile.bins.resize(static_cast<std::size_t>(nBins));
  for (int i = 0; i < nBins; ++i) {
    auto &bin = profile.bins[static_cast<std::size_t>(i)];
    bin.zLo = origin + profile.binWidth * i;
    bin.zHi = bin.zLo + profile.binWidth;
    bin.zMid = 0.5 * (bin.zLo + bin.zHi);
  }
  auto place = [&](int index) {
    if (index < 0 || static_cast<std::size_t>(index) >= cloud.pts.size()) {
      return -1;
    }
    if (profile.binWidth <= 0) {
      return 0;
    }
    const double z = fold(axisOf(cloud.pts[static_cast<std::size_t>(index)], axis), origin, length);
    int bin = static_cast<int>(std::floor((z - origin) / profile.binWidth));
    if (bin < 0) {
      bin = 0;
    }
    if (bin >= nBins) {
      bin = nBins - 1;
    }
    return bin;
  };
  for (int water : waterIndices) {
    const int bin = place(water);
    if (bin < 0) {
      continue;
    }
    auto &slot = profile.bins[static_cast<std::size_t>(bin)];
    slot.nWater += 1;
    if (static_cast<std::size_t>(water) < iceFlag.size() && iceFlag[static_cast<std::size_t>(water)]) {
      slot.nIceWater += 1;
    }
  }
  for (int ion : ionIndices) {
    const int bin = place(ion);
    if (bin < 0) {
      continue;
    }
    profile.bins[static_cast<std::size_t>(bin)].nIon += 1;
  }
  for (auto &bin : profile.bins) {
    bin.iceFraction = bin.nWater > 0 ? static_cast<double>(bin.nIceWater) / static_cast<double>(bin.nWater)
                                     : 0;
    bin.molality = molalityOf(bin.nIon, bin.nWater, ionsPerFormula);
  }
  return profile;
}

FrontPair locateFronts(const ConcentrationProfile &profile, double cut) {
  FrontPair fronts;
  const int n = static_cast<int>(profile.bins.size());
  if (n == 0 || profile.length <= 0) {
    return fronts;
  }
  // Bins with no water are gaps between crystal planes, not liquid.
  // The front is located on the occupied bins only.
  std::vector<int> occupied;
  occupied.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    if (profile.bins[static_cast<std::size_t>(i)].nWater > 0) {
      occupied.push_back(i);
    }
  }
  const int m = static_cast<int>(occupied.size());
  if (m == 0) {
    return fronts;
  }
  std::vector<char> solid(static_cast<std::size_t>(m), 0);
  int nSolid = 0;
  for (int i = 0; i < m; ++i) {
    const auto &bin = profile.bins[static_cast<std::size_t>(occupied[static_cast<std::size_t>(i)])];
    solid[static_cast<std::size_t>(i)] = bin.iceFraction >= cut;
    nSolid += solid[static_cast<std::size_t>(i)] ? 1 : 0;
  }
  if (nSolid == 0) {
    return fronts;
  }
  const double origin = profile.zOrigin;
  const double length = profile.length;
  if (nSolid == m) {
    fronts.found = true;
    fronts.wrapped = true;
    fronts.thickness = length;
    fronts.zEnter = origin;
    fronts.zExit = origin;
    return fronts;
  }
  int bestStart = 0;
  int bestLen = -1;
  for (int i = 0; i < m; ++i) {
    const int prev = (i + m - 1) % m;
    if (!solid[static_cast<std::size_t>(i)] || solid[static_cast<std::size_t>(prev)]) {
      continue;
    }
    int len = 0;
    while (len < m && solid[static_cast<std::size_t>((i + len) % m)]) {
      ++len;
    }
    if (len > bestLen) {
      bestLen = len;
      bestStart = i;
    }
  }
  auto crossing = [&](int left, int right) {
    const int leftBin = occupied[static_cast<std::size_t>(left)];
    const int rightBin = occupied[static_cast<std::size_t>(right)];
    double zLeft = profile.bins[static_cast<std::size_t>(leftBin)].zMid;
    double zRight = profile.bins[static_cast<std::size_t>(rightBin)].zMid;
    if (zRight < zLeft) {
      zRight += length;
    }
    const double fLeft = profile.bins[static_cast<std::size_t>(leftBin)].iceFraction;
    const double fRight = profile.bins[static_cast<std::size_t>(rightBin)].iceFraction;
    const double denom = fRight - fLeft;
    double frac = std::abs(denom) < 1e-12 ? 0.5 : (cut - fLeft) / denom;
    if (frac < 0) {
      frac = 0;
    }
    if (frac > 1) {
      frac = 1;
    }
    return fold(zLeft + frac * (zRight - zLeft), origin, length);
  };
  const int enterSolid = bestStart;
  const int enterLiquid = (bestStart + m - 1) % m;
  const int exitSolid = (bestStart + bestLen - 1) % m;
  const int exitLiquid = (exitSolid + 1) % m;
  fronts.zEnter = crossing(enterLiquid, enterSolid);
  fronts.zExit = crossing(exitSolid, exitLiquid);
  fronts.wrapped = bestStart + bestLen > m;
  fronts.found = true;
  if (fronts.wrapped) {
    fronts.thickness = (origin + length - fronts.zEnter) + (fronts.zExit - origin);
  } else {
    fronts.thickness = fronts.zExit - fronts.zEnter;
    if (fronts.thickness < 0) {
      fronts.thickness += length;
    }
  }
  return fronts;
}

RejectionCoefficient rejectionCoefficient(const ConcentrationProfile &profile, double bulkIce,
                                          double bulkLiquid) {
  RejectionCoefficient out;
  for (const auto &bin : profile.bins) {
    if (bin.nWater <= 0) {
      continue;
    }
    if (bin.iceFraction >= bulkIce) {
      out.nIceWater += bin.nWater;
      out.nIceIon += bin.nIon;
    } else if (bin.iceFraction <= bulkLiquid) {
      out.nLiquidWater += bin.nWater;
      out.nLiquidIon += bin.nIon;
    }
  }
  constexpr int kFormula = 2;
  out.cIce = molalityOf(out.nIceIon, out.nIceWater, kFormula);
  out.cLiquid = molalityOf(out.nLiquidIon, out.nLiquidWater, kFormula);
  if (out.nLiquidWater > 0 && out.cLiquid > 0) {
    out.k = out.cIce / out.cLiquid;
    out.defined = true;
  } else if (out.nIceIon == 0) {
    out.k = 0;
    out.defined = out.nLiquidWater > 0 || out.nIceWater > 0;
  }
  return out;
}

LinearFit linearFit(const std::vector<double> &time, const std::vector<double> &value) {
  LinearFit fit;
  const int n = static_cast<int>(time.size());
  fit.n = n;
  if (n < 2 || value.size() != time.size()) {
    fit.n = 0;
    return fit;
  }
  double sumT = 0;
  double sumY = 0;
  for (int i = 0; i < n; ++i) {
    sumT += time[static_cast<std::size_t>(i)];
    sumY += value[static_cast<std::size_t>(i)];
  }
  const double meanT = sumT / static_cast<double>(n);
  const double meanY = sumY / static_cast<double>(n);
  double ssxx = 0;
  double ssxy = 0;
  for (int i = 0; i < n; ++i) {
    const double dt = time[static_cast<std::size_t>(i)] - meanT;
    ssxx += dt * dt;
    ssxy += dt * (value[static_cast<std::size_t>(i)] - meanY);
  }
  if (ssxx <= 0) {
    fit.intercept = meanY;
    return fit;
  }
  fit.slope = ssxy / ssxx;
  fit.intercept = meanY - fit.slope * meanT;
  if (n < 3) {
    return fit;
  }
  double sse = 0;
  for (int i = 0; i < n; ++i) {
    const double residual =
        value[static_cast<std::size_t>(i)] - (fit.intercept + fit.slope * time[static_cast<std::size_t>(i)]);
    sse += residual * residual;
  }
  const double variance = sse / static_cast<double>(n - 2);
  fit.slopeStderr = std::sqrt(variance / ssxx);
  return fit;
}

PocketCensus brinePockets(const Cloud &cloud, const std::vector<bool> &iceFlag,
                          const std::vector<int> &waterIndices, const std::vector<int> &ionIndices,
                          double cutoff) {
  PocketCensus census;
  std::vector<int> liquid;
  for (int water : waterIndices) {
    if (water < 0 || static_cast<std::size_t>(water) >= iceFlag.size()) {
      continue;
    }
    if (!iceFlag[static_cast<std::size_t>(water)]) {
      liquid.push_back(water);
    }
  }
  const int n = static_cast<int>(liquid.size());
  std::vector<int> parent(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    parent[static_cast<std::size_t>(i)] = i;
  }
  auto find = [&](int x) {
    int root = x;
    while (parent[static_cast<std::size_t>(root)] != root) {
      root = parent[static_cast<std::size_t>(root)];
    }
    while (parent[static_cast<std::size_t>(x)] != root) {
      const int next = parent[static_cast<std::size_t>(x)];
      parent[static_cast<std::size_t>(x)] = root;
      x = next;
    }
    return root;
  };
  auto unite = [&](int a, int b) {
    a = find(a);
    b = find(b);
    if (a != b) {
      parent[static_cast<std::size_t>(b)] = a;
    }
  };
  const Cells cells = makeCells(cloud, liquid, cutoff);
  const double cut2 = cutoff * cutoff;
  const gen::FracBox frame = gen::makeFracBox(cloud);
  // Cloud index to the liquid slot. The cell list already returns cloud
  // indices, so the old linear scan of `liquid` dominated the frame.
  std::vector<int> slot(cloud.pts.size(), -1);
  for (int i = 0; i < n; ++i) {
    const int index = liquid[static_cast<std::size_t>(i)];
    if (index >= 0 && static_cast<std::size_t>(index) < slot.size()) {
      slot[static_cast<std::size_t>(index)] = i;
    }
  }
  for (int i = 0; i < n; ++i) {
    visitNearby(cells, liquid[static_cast<std::size_t>(i)], cloud, [&](int other) {
      if (other < 0 || static_cast<std::size_t>(other) >= slot.size()) {
        return;
      }
      const int j = slot[static_cast<std::size_t>(other)];
      if (j < 0 || j == i) {
        return;
      }
      if (gen::periodicDistSq(frame, cloud, liquid[static_cast<std::size_t>(i)], other) <= cut2) {
        unite(i, j);
      }
    });
  }
  std::vector<std::vector<int>> groups(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    groups[static_cast<std::size_t>(find(i))].push_back(liquid[static_cast<std::size_t>(i)]);
  }
  int reservoir = -1;
  int reservoirSize = -1;
  for (int i = 0; i < n; ++i) {
    const int size = static_cast<int>(groups[static_cast<std::size_t>(i)].size());
    if (size > reservoirSize) {
      reservoirSize = size;
      reservoir = i;
    }
  }
  if (reservoir >= 0) {
    census.reservoirWaters = reservoirSize;
  }
  std::vector<int> componentOf(cloud.pts.size(), -1);
  for (int i = 0; i < n; ++i) {
    if (groups[static_cast<std::size_t>(i)].empty()) {
      continue;
    }
    const int id = i == reservoir ? -2 : i;
    for (int water : groups[static_cast<std::size_t>(i)]) {
      componentOf[static_cast<std::size_t>(water)] = id;
    }
    if (i == reservoir) {
      continue;
    }
    BrinePocket pocket;
    pocket.waters = groups[static_cast<std::size_t>(i)];
    census.pockets.push_back(std::move(pocket));
  }
  for (int ion : ionIndices) {
    if (ion < 0 || static_cast<std::size_t>(ion) >= cloud.pts.size()) {
      continue;
    }
    int nearest = -1;
    double best = cutoff * cutoff;
    visitNearby(cells, ion, cloud, [&](int other) {
      if (other < 0 || static_cast<std::size_t>(other) >= slot.size()) {
        return;
      }
      if (slot[static_cast<std::size_t>(other)] < 0) {
        return;
      }
      const double r2 = gen::periodicDistSq(frame, cloud, ion, other);
      if (r2 <= best) {
        best = r2;
        nearest = other;
      }
    });
    if (nearest < 0) {
      continue;
    }
    const int comp = componentOf[static_cast<std::size_t>(nearest)];
    if (comp < 0) {
      continue;
    }
    for (auto &pocket : census.pockets) {
      if (!pocket.waters.empty() &&
          componentOf[static_cast<std::size_t>(pocket.waters.front())] == comp) {
        pocket.nIons += 1;
        break;
      }
    }
  }
  for (const auto &pocket : census.pockets) {
    census.sizes.push_back(static_cast<int>(pocket.waters.size()));
  }
  std::sort(census.sizes.begin(), census.sizes.end());
  census.nPockets = static_cast<int>(census.pockets.size());
  return census;
}

FrontLayerReport frontLayer(const Cloud &cloud, const std::vector<int> &waterIndices,
                            const FrontPair &fronts, double halfWidth, double cutoff, int maxRing,
                            int axis) {
  FrontLayerReport report;
  report.coordination.assign(static_cast<std::size_t>(std::max(maxRing, 4)) + 1, 0);
  report.ringCensus.assign(static_cast<std::size_t>(std::max(maxRing, 0)) + 1, 0);
  if (!fronts.found || halfWidth <= 0) {
    return report;
  }
  const double origin = boxOrigin(cloud, axis);
  const double length = boxLength(cloud, axis);
  auto nearFront = [&](int index) {
    const double z = fold(axisOf(cloud.pts[static_cast<std::size_t>(index)], axis), origin, length);
    auto gap = [&](double plane) {
      double d = z - plane;
      if (length > 0) {
        d -= length * std::round(d / length);
      }
      return std::abs(d);
    };
    return gap(fronts.zEnter) <= halfWidth || gap(fronts.zExit) <= halfWidth;
  };
  std::vector<char> inLayer(cloud.pts.size(), 0);
  for (int water : waterIndices) {
    if (water < 0 || static_cast<std::size_t>(water) >= cloud.pts.size()) {
      continue;
    }
    if (nearFront(water)) {
      inLayer[static_cast<std::size_t>(water)] = 1;
      report.waters.push_back(water);
    }
  }
  const Cells cells = makeCells(cloud, waterIndices, cutoff);
  const double cut2 = cutoff * cutoff;
  const gen::FracBox frame = gen::makeFracBox(cloud);
  // Rings are searched on the front layer only. A row per atom in the
  // whole cell made the primitive-ring walk pay for the bulk ice.
  std::vector<int> localOf(cloud.pts.size(), -1);
  for (int i = 0; i < static_cast<int>(report.waters.size()); ++i) {
    localOf[static_cast<std::size_t>(report.waters[static_cast<std::size_t>(i)])] = i;
  }
  std::vector<std::vector<int>> rows(report.waters.size());
  for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
    rows[static_cast<std::size_t>(i)].push_back(i);
  }
  auto addCoord = [&](int count) {
    const int cap = static_cast<int>(report.coordination.size()) - 1;
    const int bin = std::min(std::max(count, 0), cap);
    report.coordination[static_cast<std::size_t>(bin)] += 1;
    if (count == 4) {
      report.nFour += 1;
    } else {
      report.nDefect += 1;
    }
  };
  for (int water : report.waters) {
    int neighbours = 0;
    std::unordered_set<int> seen;
    visitNearby(cells, water, cloud, [&](int other) {
      if (other == water || !seen.insert(other).second) {
        return;
      }
      if (gen::periodicDistSq(frame, cloud, water, other) > cut2) {
        return;
      }
      neighbours += 1;
      const int local = localOf[static_cast<std::size_t>(other)];
      if (local >= 0) {
        rows[static_cast<std::size_t>(localOf[static_cast<std::size_t>(water)])].push_back(local);
      }
    });
    addCoord(neighbours);
    if (static_cast<std::size_t>(water) < cloud.pts.size()) {
      const auto label = cloud.pts[static_cast<std::size_t>(water)].iceType;
      if (label == molSys::atom_state_type::hexagonal || label == molSys::atom_state_type::reHex) {
        report.nHexagonal += 1;
      } else if (label == molSys::atom_state_type::cubic || label == molSys::atom_state_type::reCubic) {
        report.nCubic += 1;
      } else if (label == molSys::atom_state_type::interfacial ||
                 label == molSys::atom_state_type::interClathrate) {
        report.nInterfacial += 1;
      } else if (label == molSys::atom_state_type::clathrate) {
        report.nClathrate += 1;
      } else if (label == molSys::atom_state_type::water) {
        report.nLiquidLabel += 1;
      }
    }
  }
  if (maxRing >= 3 && !report.waters.empty()) {
    const auto rings = primitive::ringNetwork(rows, maxRing + 1);
    for (const auto &ring : rings) {
      if (ring.empty()) {
        continue;
      }
      const int size = static_cast<int>(ring.size());
      if (size >= 0 && size < static_cast<int>(report.ringCensus.size())) {
        report.ringCensus[static_cast<std::size_t>(size)] += 1;
      }
    }
  }
  return report;
}

std::vector<TrapEvent> trappingTransitions(const IonEnvironment &before, const IonEnvironment &after) {
  std::vector<TrapEvent> events;
  const std::size_t n = std::min(before.state.size(), after.state.size());
  events.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (before.state[i] == after.state[i]) {
      continue;
    }
    TrapEvent event;
    event.ion = static_cast<int>(i);
    event.from = before.state[i];
    event.to = after.state[i];
    events.push_back(event);
  }
  return events;
}

FrontSnapshot analyseFront(const Cloud &cloud, const std::vector<bool> &iceFlag,
                           const std::vector<int> &waterIndices, const std::vector<int> &ionIndices,
                           const FrontOptions &options) {
  FrontSnapshot snap;
  snap.profile = concentrationProfile(cloud, iceFlag, waterIndices, ionIndices, options.nBins,
                                      options.axis, options.ionsPerFormula);
  snap.fronts = locateFronts(snap.profile, options.frontCut);
  snap.rejection = rejectionCoefficient(snap.profile, options.bulkIce, options.bulkLiquid);
  snap.ions = ionEnvironment(cloud, iceFlag, ionIndices, 0, options.shellCutoff);
  // ionEnvironment's waterType 0 counts every non-ion. Restrict the flag
  // array so hydrogens are not a second species in the shell: pass a cloud
  // view by clearing ice flags of non-waters and letting waterType 0 skip
  // ions only. Callers that include hydrogens should pass oxygen indices
  // and a cloud of oxygens plus ions. The shell cutoff still uses whatever
  // atoms are not listed as ions.
  snap.pockets = brinePockets(cloud, iceFlag, waterIndices, ionIndices, options.shellCutoff);
  snap.layer = frontLayer(cloud, waterIndices, snap.fronts, options.layerHalfWidth,
                          options.shellCutoff, options.maxRing, options.axis);
  return snap;
}

} // namespace site
