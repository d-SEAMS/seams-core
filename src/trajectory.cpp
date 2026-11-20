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

#include <trajectory.hpp>

#include <generic.hpp>
#ifdef SEAMS_HAS_MINIMAGE
#include <minimage.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <utility>

namespace traj {

namespace {

bool finiteArray(const double *ptr, std::size_t n) noexcept {
  if (n > 0 && ptr == nullptr) {
    return false;
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(ptr[i])) {
      return false;
    }
  }
  return true;
}

std::uint64_t packEdge(int u, int v, bool directed) noexcept {
  if (!directed && u > v) {
    std::swap(u, v);
  }
  const auto uu = static_cast<std::uint32_t>(u);
  const auto vv = static_cast<std::uint32_t>(v);
  return (static_cast<std::uint64_t>(uu) << 32) | static_cast<std::uint64_t>(vv);
}

std::size_t sortedIntersectionCount(const std::vector<std::uint64_t> &a,
                                    const std::vector<std::uint64_t> &b) noexcept {
  std::size_t ia = 0;
  std::size_t ib = 0;
  std::size_t count = 0;
  while (ia < a.size() && ib < b.size()) {
    if (a[ia] == b[ib]) {
      ++count;
      ++ia;
      ++ib;
    } else if (a[ia] < b[ib]) {
      ++ia;
    } else {
      ++ib;
    }
  }
  return count;
}

void inplaceSortedIntersection(std::vector<std::uint64_t> &alive,
                               const std::vector<std::uint64_t> &cur) noexcept {
  std::size_t ia = 0;
  std::size_t ic = 0;
  std::size_t write = 0;
  while (ia < alive.size() && ic < cur.size()) {
    if (alive[ia] == cur[ic]) {
      alive[write++] = alive[ia];
      ++ia;
      ++ic;
    } else if (alive[ia] < cur[ic]) {
      ++ia;
    } else {
      ++ic;
    }
  }
  alive.resize(write);
}

struct LinearRegression {
  double intercept = 0.0;
  double slope = 0.0;
  double r2 = 0.0;
  bool ok = false;
};

LinearRegression fitWeightedLine(const std::vector<double> &x,
                                 const std::vector<double> &y,
                                 const std::vector<double> &w) noexcept {
  LinearRegression res;
  if (x.size() < 2 || x.size() != y.size() || x.size() != w.size()) {
    return res;
  }
  double sw = 0.0;
  double swx = 0.0;
  double swy = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (!(w[i] > 0.0) || !std::isfinite(x[i]) || !std::isfinite(y[i])) {
      continue;
    }
    sw += w[i];
    swx += w[i] * x[i];
    swy += w[i] * y[i];
  }
  if (!(sw > 0.0)) {
    return res;
  }
  const double mx = swx / sw;
  const double my = swy / sw;
  double sxx = 0.0;
  double sxy = 0.0;
  double syy = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (!(w[i] > 0.0) || !std::isfinite(x[i]) || !std::isfinite(y[i])) {
      continue;
    }
    const double dx = x[i] - mx;
    const double dy = y[i] - my;
    sxx += w[i] * dx * dx;
    sxy += w[i] * dx * dy;
    syy += w[i] * dy * dy;
  }
  if (!(sxx > 0.0) || !std::isfinite(sxx)) {
    return res;
  }
  res.slope = sxy / sxx;
  res.intercept = my - res.slope * mx;
  if (!std::isfinite(res.slope) || !std::isfinite(res.intercept)) {
    return res;
  }
  if (syy > 0.0) {
    double ssRes = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      if (!(w[i] > 0.0) || !std::isfinite(x[i]) || !std::isfinite(y[i])) {
        continue;
      }
      const double err = y[i] - (res.intercept + res.slope * x[i]);
      ssRes += w[i] * err * err;
    }
    res.r2 = std::max(0.0, 1.0 - ssRes / syy);
  } else {
    res.r2 = 1.0;
  }
  res.ok = true;
  return res;
}

bool resolveTimeWindow(const std::vector<double> &lagTime, double tStart,
                       double tEnd, double &loOut, double &hiOut) noexcept {
  if (lagTime.size() < 2) {
    return false;
  }
  const double tMax = lagTime.back();
  if (!(tMax > 0.0) || !std::isfinite(tMax)) {
    return false;
  }
  if (!std::isfinite(tStart) || !std::isfinite(tEnd)) {
    return false;
  }
  if (tStart == 0.0 && tEnd == 0.0) {
    loOut = 0.1 * tMax;
    hiOut = 0.5 * tMax;
    std::size_t count = 0;
    for (double t : lagTime) {
      if (t >= loOut && t <= hiOut) {
        ++count;
      }
    }
    if (count < 2) {
      loOut = lagTime[1];
      hiOut = tMax;
    }
    return true;
  }
  if (tStart < 0.0 || tEnd <= tStart) {
    return false;
  }
  loOut = tStart;
  hiOut = tEnd;
  return true;
}

} // namespace

BoxCell BoxCell::fromBox9(const double box9[9]) noexcept {
  BoxCell c;
  if (box9 == nullptr || !finiteArray(box9, 9)) {
    return c;
  }
  c.ox = box9[0];
  c.oy = box9[1];
  c.oz = box9[2];
  c.lx = box9[3];
  c.ly = box9[4];
  c.lz = box9[5];
  c.xy = box9[6];
  c.xz = box9[7];
  c.yz = box9[8];
  const double vol = c.lx * c.ly * c.lz;
  if (!(c.lx > 0.0 && c.ly > 0.0 && c.lz > 0.0) || !(vol > 0.0) ||
      !std::isfinite(vol)) {
    return c;
  }
  c.triclinic = (std::abs(c.xy) > 0.0 || std::abs(c.xz) > 0.0 ||
                 std::abs(c.yz) > 0.0);
  const double bcz = c.xy * c.yz - c.ly * c.xz;
  const double crossBc = std::sqrt((c.ly * c.lz) * (c.ly * c.lz) +
                                   (c.xy * c.lz) * (c.xy * c.lz) + bcz * bcz);
  const double crossAc = std::sqrt(c.lz * c.lz + c.yz * c.yz);
  if (!(crossBc > 0.0 && crossAc > 0.0) || !std::isfinite(crossBc) ||
      !std::isfinite(crossAc)) {
    return c;
  }
  const double wx = vol / crossBc;
  const double wy = (c.ly * c.lz) / crossAc;
  const double wz = c.lz;
  c.halfMin = 0.5 * std::min(wx, std::min(wy, wz));
  if (!(c.halfMin > 0.0) || !std::isfinite(c.halfMin)) {
    return c;
  }
  c.ok = true;
  return c;
}

BoxCell BoxCell::fromCloud(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud) noexcept {
  const gen::FracBox fb = gen::makeFracBox(cloud);
  if (!fb.ok) {
    return BoxCell{};
  }
  const double b9[9] = {fb.ox, fb.oy, fb.oz, fb.lx, fb.ly,
                        fb.lz, fb.xy, fb.xz, fb.yz};
  return fromBox9(b9);
}

void BoxCell::toBox9(double out[9]) const noexcept {
  out[0] = ox;
  out[1] = oy;
  out[2] = oz;
  out[3] = lx;
  out[4] = ly;
  out[5] = lz;
  out[6] = xy;
  out[7] = xz;
  out[8] = yz;
}

std::array<double, 3> BoxCell::micDelta(double dx, double dy,
                                        double dz) const noexcept {
  if (!triclinic) {
    return {dx - lx * std::nearbyint(dx / lx),
            dy - ly * std::nearbyint(dy / ly),
            dz - lz * std::nearbyint(dz / lz)};
  }
  const double sz = dz / lz;
  const double sy = (dy - yz * sz) / ly;
  const double sx = (dx - xy * sy - xz * sz) / lx;
  const double dsx = sx - std::nearbyint(sx);
  const double dsy = sy - std::nearbyint(sy);
  const double dsz = sz - std::nearbyint(sz);
  const double rx = lx * dsx + xy * dsy + xz * dsz;
  const double ry = ly * dsy + yz * dsz;
  const double rz = lz * dsz;
  if (std::sqrt(rx * rx + ry * ry + rz * rz) + 1e-12 < halfMin) {
    return {rx, ry, rz};
  }
#ifdef SEAMS_HAS_MINIMAGE
  // minimage takes LAMMPS dump bounds, which widen lx and ly by the tilts.
  const double xlo = std::min({0.0, xy, xz, xy + xz});
  const double xhi = std::max({0.0, xy, xz, xy + xz});
  const double ylo = std::min(0.0, yz);
  const double yhi = std::max(0.0, yz);
  mi_cell cell;
  if (mi_cell_from_lammps_bounds(lx + xhi - xlo, ly + yhi - ylo, lz, xy, xz,
                                 yz, ox + xlo, oy + ylo, oz, &cell) == 0) {
    const double zero[3] = {0.0, 0.0, 0.0};
    const double target[3] = {dx, dy, dz};
    double dr[3] = {0.0, 0.0, 0.0};
    if (mi_displacement_euclidean(&cell, zero, target, dr) == 0) {
      return {dr[0], dr[1], dr[2]};
    }
  }
#endif
  return micDeltaSearch(dx, dy, dz);
}

std::array<double, 3> BoxCell::micDeltaSearch(double dx, double dy,
                                              double dz) const noexcept {
  if (!triclinic) {
    return {dx - lx * std::nearbyint(dx / lx),
            dy - ly * std::nearbyint(dy / ly),
            dz - lz * std::nearbyint(dz / lz)};
  }
  const double sz = dz / lz;
  const double sy = (dy - yz * sz) / ly;
  const double sx = (dx - xy * sy - xz * sz) / lx;
  const double dsx = sx - std::nearbyint(sx);
  const double dsy = sy - std::nearbyint(sy);
  const double dsz = sz - std::nearbyint(sz);
  double rx = lx * dsx + xy * dsy + xz * dsz;
  double ry = ly * dsy + yz * dsz;
  double rz = lz * dsz;
  double r2 = rx * rx + ry * ry + rz * rz;
  // The minimum image v = H t is no longer than the fractional image r, and
  // |t_i| <= |v| * |row i of H^-1|, which bounds every lattice shift to try.
  const double g0 = std::sqrt(1.0 / (lx * lx) + (xy * xy) / (lx * lx * ly * ly) +
                              std::pow((xy * yz - xz * ly) / (lx * ly * lz), 2));
  const double g1 = std::sqrt(1.0 / (ly * ly) + (yz * yz) / (ly * ly * lz * lz));
  const double g2 = 1.0 / lz;
  const double rmax = std::sqrt(r2) * (1.0 + 1e-12);
  const auto lo = [](double ds, double span) {
    return static_cast<long>(std::ceil(ds - span));
  };
  const auto hi = [](double ds, double span) {
    return static_cast<long>(std::floor(ds + span));
  };
  const long ix0 = lo(dsx, rmax * g0), ix1 = hi(dsx, rmax * g0);
  const long iy0 = lo(dsy, rmax * g1), iy1 = hi(dsy, rmax * g1);
  const long iz0 = lo(dsz, rmax * g2), iz1 = hi(dsz, rmax * g2);
  for (long ix = ix0; ix <= ix1; ++ix) {
    for (long iy = iy0; iy <= iy1; ++iy) {
      for (long iz = iz0; iz <= iz1; ++iz) {
        const double tx = dsx - static_cast<double>(ix);
        const double ty = dsy - static_cast<double>(iy);
        const double tz = dsz - static_cast<double>(iz);
        const double cx = lx * tx + xy * ty + xz * tz;
        const double cy = ly * ty + yz * tz;
        const double cz = lz * tz;
        const double c2 = cx * cx + cy * cy + cz * cz;
        if (c2 < r2) {
          r2 = c2;
          rx = cx;
          ry = cy;
          rz = cz;
        }
      }
    }
  }
  return {rx, ry, rz};
}

TrajectoryUnwrapper::TrajectoryUnwrapper(std::size_t nParticles) {
  reset(nParticles);
}

void TrajectoryUnwrapper::reset(std::size_t nParticles) {
  nParticles_ = nParticles;
  framesSeen_ = 0;
  prevWrapped_.assign(3 * nParticles, 0.0);
  currentUnwrapped_.assign(3 * nParticles, 0.0);
}

bool TrajectoryUnwrapper::pushFrame(const double *xyz, std::size_t nParticles,
                                    const BoxCell &cell, double *unwrappedOut) {
  if (!cell.ok || nParticles == 0 || xyz == nullptr) {
    return false;
  }
  if (!finiteArray(xyz, 3 * nParticles)) {
    return false;
  }
  if (framesSeen_ == 0 && nParticles_ == 0) {
    reset(nParticles);
  }
  if (nParticles != nParticles_) {
    return false;
  }
  if (framesSeen_ == 0) {
    std::copy_n(xyz, 3 * nParticles_, prevWrapped_.data());
    std::copy_n(xyz, 3 * nParticles_, currentUnwrapped_.data());
  } else {
    for (std::size_t i = 0; i < nParticles_; ++i) {
      const std::size_t base = 3 * i;
      const double dx = xyz[base] - currentUnwrapped_[base];
      const double dy = xyz[base + 1] - currentUnwrapped_[base + 1];
      const double dz = xyz[base + 2] - currentUnwrapped_[base + 2];
      const auto step = cell.micDelta(dx, dy, dz);
      currentUnwrapped_[base] += step[0];
      currentUnwrapped_[base + 1] += step[1];
      currentUnwrapped_[base + 2] += step[2];
      prevWrapped_[base] = xyz[base];
      prevWrapped_[base + 1] = xyz[base + 1];
      prevWrapped_[base + 2] = xyz[base + 2];
    }
  }
  if (unwrappedOut != nullptr) {
    std::copy_n(currentUnwrapped_.data(), 3 * nParticles_, unwrappedOut);
  }
  ++framesSeen_;
  return true;
}

bool unwrapPositions(const double *xyzFrames, std::size_t nFrames,
                     std::size_t nAtoms, const double *boxes,
                     std::size_t nBoxes, double *unwrappedOut) {
  if (xyzFrames == nullptr || boxes == nullptr || unwrappedOut == nullptr) {
    return false;
  }
  if (nFrames == 0 || nAtoms == 0 || (nBoxes != 1 && nBoxes != nFrames)) {
    return false;
  }
  TrajectoryUnwrapper unwrapper(nAtoms);
  const std::size_t stride = 3 * nAtoms;
  for (std::size_t f = 0; f < nFrames; ++f) {
    const double *b9 = boxes + 9 * (nBoxes == 1 ? 0 : f);
    const BoxCell cell = BoxCell::fromBox9(b9);
    if (!unwrapper.pushFrame(xyzFrames + f * stride, nAtoms, cell,
                             unwrappedOut + f * stride)) {
      return false;
    }
  }
  return true;
}

bool moleculeCentres(const double *xyz, std::size_t nAtoms, const BoxCell &cell,
                     const int *molIds, const double *masses,
                     std::vector<int> &uniqueMolIds,
                     std::vector<double> &centresOut) {
  uniqueMolIds.clear();
  centresOut.clear();
  if (!cell.ok || nAtoms == 0 || xyz == nullptr || molIds == nullptr) {
    return false;
  }
  if (!finiteArray(xyz, 3 * nAtoms)) {
    return false;
  }
  if (masses != nullptr) {
    for (std::size_t i = 0; i < nAtoms; ++i) {
      if (!std::isfinite(masses[i]) || !(masses[i] > 0.0)) {
        return false;
      }
    }
  }
  std::unordered_map<int, std::size_t> molToIdx;
  molToIdx.reserve(nAtoms);
  std::vector<std::size_t> anchorAtom;
  std::vector<double> weightSums;

  for (std::size_t i = 0; i < nAtoms; ++i) {
    const int mid = molIds[i];
    const double w = (masses != nullptr) ? masses[i] : 1.0;
    auto [it, inserted] = molToIdx.emplace(mid, uniqueMolIds.size());
    if (inserted) {
      uniqueMolIds.push_back(mid);
      anchorAtom.push_back(i);
      weightSums.push_back(w);
      centresOut.push_back(w * xyz[3 * i]);
      centresOut.push_back(w * xyz[3 * i + 1]);
      centresOut.push_back(w * xyz[3 * i + 2]);
    } else {
      const std::size_t m = it->second;
      const std::size_t a0 = anchorAtom[m];
      const double rx0 = xyz[3 * a0];
      const double ry0 = xyz[3 * a0 + 1];
      const double rz0 = xyz[3 * a0 + 2];
      const auto dr = cell.micDelta(xyz[3 * i] - rx0, xyz[3 * i + 1] - ry0,
                                    xyz[3 * i + 2] - rz0);
      weightSums[m] += w;
      centresOut[3 * m] += w * (rx0 + dr[0]);
      centresOut[3 * m + 1] += w * (ry0 + dr[1]);
      centresOut[3 * m + 2] += w * (rz0 + dr[2]);
    }
  }

  for (std::size_t m = 0; m < uniqueMolIds.size(); ++m) {
    const double invW = 1.0 / weightSums[m];
    centresOut[3 * m] *= invW;
    centresOut[3 * m + 1] *= invW;
    centresOut[3 * m + 2] *= invW;
  }
  return true;
}

bool moleculeCentres(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const double *masses, std::vector<int> &uniqueMolIds,
    std::vector<double> &centresOut) {
  const BoxCell cell = BoxCell::fromCloud(cloud);
  if (!cell.ok || cloud.nop <= 0 ||
      cloud.pts.size() < static_cast<std::size_t>(cloud.nop)) {
    uniqueMolIds.clear();
    centresOut.clear();
    return false;
  }
  const auto n = static_cast<std::size_t>(cloud.nop);
  std::vector<double> xyz(3 * n);
  std::vector<int> mols(n);
  for (std::size_t i = 0; i < n; ++i) {
    xyz[3 * i] = cloud.pts[i].x;
    xyz[3 * i + 1] = cloud.pts[i].y;
    xyz[3 * i + 2] = cloud.pts[i].z;
    mols[i] = cloud.pts[i].molID;
  }
  return moleculeCentres(xyz.data(), n, cell, mols.data(), masses, uniqueMolIds,
                         centresOut);
}

bool unwrapMoleculeCentres(const double *xyzFrames, std::size_t nFrames,
                           std::size_t nAtoms, const double *boxes,
                           std::size_t nBoxes, const int *molIds,
                           const double *masses,
                           std::vector<int> &uniqueMolIds,
                           std::vector<double> &unwrappedCentresOut) {
  uniqueMolIds.clear();
  unwrappedCentresOut.clear();
  if (xyzFrames == nullptr || boxes == nullptr || molIds == nullptr) {
    return false;
  }
  if (nFrames == 0 || nAtoms == 0 || (nBoxes != 1 && nBoxes != nFrames)) {
    return false;
  }

  std::vector<int> frameMols;
  std::vector<double> frameCentres;
  TrajectoryUnwrapper unwrapper;
  const std::size_t atomStride = 3 * nAtoms;

  for (std::size_t f = 0; f < nFrames; ++f) {
    const double *b9 = boxes + 9 * (nBoxes == 1 ? 0 : f);
    const BoxCell cell = BoxCell::fromBox9(b9);
    if (!moleculeCentres(xyzFrames + f * atomStride, nAtoms, cell, molIds,
                         masses, frameMols, frameCentres)) {
      uniqueMolIds.clear();
      unwrappedCentresOut.clear();
      return false;
    }
    const std::size_t nMols = frameMols.size();
    if (f == 0) {
      uniqueMolIds = frameMols;
      unwrappedCentresOut.resize(nFrames * 3 * nMols);
      unwrapper.reset(nMols);
    }
    if (!unwrapper.pushFrame(frameCentres.data(), nMols, cell,
                             unwrappedCentresOut.data() + f * 3 * nMols)) {
      uniqueMolIds.clear();
      unwrappedCentresOut.clear();
      return false;
    }
  }
  return true;
}

bool waterOrientationVectors(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    int oxygenType, int hydrogenType, WaterVectorsResult &out,
    double covalentCutoff) {
  out = WaterVectorsResult{};
  if (oxygenType == hydrogenType || !(covalentCutoff > 0.0) ||
      !std::isfinite(covalentCutoff)) {
    return false;
  }
  const BoxCell cell = BoxCell::fromCloud(cloud);
  if (!cell.ok || cloud.nop <= 0 ||
      cloud.pts.size() < static_cast<std::size_t>(cloud.nop)) {
    return false;
  }
  const auto n = static_cast<std::size_t>(cloud.nop);
  std::vector<std::size_t> oxygens;
  std::vector<std::size_t> hydrogens;
  std::unordered_map<int, std::vector<std::size_t>> hByMol;
  std::unordered_map<int, std::size_t> oCountByMol;

  for (std::size_t i = 0; i < n; ++i) {
    const auto &pt = cloud.pts[i];
    if (!std::isfinite(pt.x) || !std::isfinite(pt.y) || !std::isfinite(pt.z)) {
      return false;
    }
    if (pt.type == oxygenType) {
      oxygens.push_back(i);
      oCountByMol[pt.molID] += 1;
    } else if (pt.type == hydrogenType) {
      hydrogens.push_back(i);
      hByMol[pt.molID].push_back(i);
    }
  }
  if (oxygens.empty()) {
    return true;
  }

  bool hasMolGrouping = false;
  for (std::size_t oi : oxygens) {
    const int mid = cloud.pts[oi].molID;
    auto hit = hByMol.find(mid);
    if (oCountByMol[mid] == 1 && hit != hByMol.end() &&
        hit->second.size() == 2) {
      hasMolGrouping = true;
      break;
    }
  }

  auto appendWater = [&](std::size_t oi, std::size_t h1, std::size_t h2) {
    const auto &po = cloud.pts[oi];
    const auto &p1 = cloud.pts[h1];
    const auto &p2 = cloud.pts[h2];
    const auto d1 = cell.micDelta(p1.x - po.x, p1.y - po.y, p1.z - po.z);
    const auto d2 = cell.micDelta(p2.x - po.x, p2.y - po.y, p2.z - po.z);

    const double mx = 0.5 * (d1[0] + d2[0]);
    const double my = 0.5 * (d1[1] + d2[1]);
    const double mz = 0.5 * (d1[2] + d2[2]);
    const double mNorm = std::sqrt(mx * mx + my * my + mz * mz);

    const double hx = d2[0] - d1[0];
    const double hy = d2[1] - d1[1];
    const double hz = d2[2] - d1[2];
    const double hNorm = std::sqrt(hx * hx + hy * hy + hz * hz);

    out.molIds.push_back(po.molID);
    out.oxygenIndices.push_back(static_cast<int>(oi));
    out.ohMidVectors.insert(out.ohMidVectors.end(), {mx, my, mz});
    if (mNorm > 1e-15) {
      out.bisectorUnit.insert(out.bisectorUnit.end(),
                              {mx / mNorm, my / mNorm, mz / mNorm});
    } else {
      out.bisectorUnit.insert(out.bisectorUnit.end(), {0.0, 0.0, 0.0});
    }
    out.hhVectors.insert(out.hhVectors.end(), {hx, hy, hz});
    if (hNorm > 1e-15) {
      out.hhUnit.insert(out.hhUnit.end(),
                        {hx / hNorm, hy / hNorm, hz / hNorm});
    } else {
      out.hhUnit.insert(out.hhUnit.end(), {0.0, 0.0, 0.0});
    }
  };

  if (hasMolGrouping) {
    for (std::size_t oi : oxygens) {
      const int mid = cloud.pts[oi].molID;
      auto hit = hByMol.find(mid);
      if (oCountByMol[mid] != 1 || hit == hByMol.end() ||
          hit->second.size() != 2) {
        continue;
      }
      std::size_t h1 = hit->second[0];
      std::size_t h2 = hit->second[1];
      if (cloud.pts[h1].atomID > cloud.pts[h2].atomID) {
        std::swap(h1, h2);
      }
      appendWater(oi, h1, h2);
    }
    return true;
  }

  const double cut2 = covalentCutoff * covalentCutoff;
  std::vector<std::uint8_t> hUsed(hydrogens.size(), 0);
  for (std::size_t oi : oxygens) {
    const auto &po = cloud.pts[oi];
    double best1 = cut2;
    double best2 = cut2;
    std::size_t idx1 = hydrogens.size();
    std::size_t idx2 = hydrogens.size();
    for (std::size_t k = 0; k < hydrogens.size(); ++k) {
      if (hUsed[k] != 0) {
        continue;
      }
      const auto &ph = cloud.pts[hydrogens[k]];
      const auto dr = cell.micDelta(ph.x - po.x, ph.y - po.y, ph.z - po.z);
      const double r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
      if (r2 < best1) {
        best2 = best1;
        idx2 = idx1;
        best1 = r2;
        idx1 = k;
      } else if (r2 < best2) {
        best2 = r2;
        idx2 = k;
      }
    }
    if (idx1 < hydrogens.size() && idx2 < hydrogens.size()) {
      hUsed[idx1] = 1;
      hUsed[idx2] = 1;
      std::size_t h1 = hydrogens[idx1];
      std::size_t h2 = hydrogens[idx2];
      if (cloud.pts[h1].atomID > cloud.pts[h2].atomID) {
        std::swap(h1, h2);
      }
      appendWater(oi, h1, h2);
    }
  }
  return true;
}

bool moleculeDipoleVectors(const double *xyz, std::size_t nAtoms,
                           const BoxCell &cell, const int *molIds,
                           const double *charges, DipoleVectorsResult &out) {
  out = DipoleVectorsResult{};
  if (!cell.ok || nAtoms == 0 || xyz == nullptr || molIds == nullptr) {
    return false;
  }
  if (!finiteArray(xyz, 3 * nAtoms)) {
    return false;
  }
  if (charges != nullptr && !finiteArray(charges, nAtoms)) {
    return false;
  }

  std::unordered_map<int, std::size_t> molToIdx;
  std::vector<std::vector<std::size_t>> molAtoms;
  for (std::size_t i = 0; i < nAtoms; ++i) {
    const int mid = molIds[i];
    auto [it, inserted] = molToIdx.emplace(mid, out.molIds.size());
    if (inserted) {
      out.molIds.push_back(mid);
      molAtoms.push_back({i});
    } else {
      molAtoms[it->second].push_back(i);
    }
  }

  const std::size_t nMols = out.molIds.size();
  out.dipoles.resize(3 * nMols, 0.0);
  out.unitVectors.resize(3 * nMols, 0.0);
  out.magnitudes.resize(nMols, 0.0);

  std::vector<std::array<double, 3>> relCoords;
  for (std::size_t m = 0; m < nMols; ++m) {
    const auto &atoms = molAtoms[m];
    const std::size_t k = atoms.size();
    const std::size_t a0 = atoms[0];
    const double x0 = xyz[3 * a0];
    const double y0 = xyz[3 * a0 + 1];
    const double z0 = xyz[3 * a0 + 2];

    relCoords.resize(k);
    relCoords[0] = {0.0, 0.0, 0.0};
    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;
    for (std::size_t j = 1; j < k; ++j) {
      const std::size_t aj = atoms[j];
      relCoords[j] = cell.micDelta(xyz[3 * aj] - x0, xyz[3 * aj + 1] - y0,
                                   xyz[3 * aj + 2] - z0);
      cx += relCoords[j][0];
      cy += relCoords[j][1];
      cz += relCoords[j][2];
    }

    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    if (charges != nullptr) {
      const double invK = 1.0 / static_cast<double>(k);
      cx *= invK;
      cy *= invK;
      cz *= invK;
      for (std::size_t j = 0; j < k; ++j) {
        const double q = charges[atoms[j]];
        dx += q * (relCoords[j][0] - cx);
        dy += q * (relCoords[j][1] - cy);
        dz += q * (relCoords[j][2] - cz);
      }
    } else if (k > 1) {
      const double invOther = 1.0 / static_cast<double>(k - 1);
      dx = cx * invOther;
      dy = cy * invOther;
      dz = cz * invOther;
    }

    const double norm = std::sqrt(dx * dx + dy * dy + dz * dz);
    out.dipoles[3 * m] = dx;
    out.dipoles[3 * m + 1] = dy;
    out.dipoles[3 * m + 2] = dz;
    out.magnitudes[m] = norm;
    if (norm > 1e-15) {
      out.unitVectors[3 * m] = dx / norm;
      out.unitVectors[3 * m + 1] = dy / norm;
      out.unitVectors[3 * m + 2] = dz / norm;
    }
  }
  return true;
}

bool moleculeDipoleVectors(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const int *molIdsOverride, const double *charges,
    DipoleVectorsResult &out) {
  const BoxCell cell = BoxCell::fromCloud(cloud);
  if (!cell.ok || cloud.nop <= 0 ||
      cloud.pts.size() < static_cast<std::size_t>(cloud.nop)) {
    out = DipoleVectorsResult{};
    return false;
  }
  const auto n = static_cast<std::size_t>(cloud.nop);
  std::vector<double> xyz(3 * n);
  std::vector<int> storedMols;
  if (molIdsOverride == nullptr) {
    storedMols.resize(n);
  }
  for (std::size_t i = 0; i < n; ++i) {
    xyz[3 * i] = cloud.pts[i].x;
    xyz[3 * i + 1] = cloud.pts[i].y;
    xyz[3 * i + 2] = cloud.pts[i].z;
    if (molIdsOverride == nullptr) {
      storedMols[i] = cloud.pts[i].molID;
    }
  }
  const int *mols =
      (molIdsOverride != nullptr) ? molIdsOverride : storedMols.data();
  return moleculeDipoleVectors(xyz.data(), n, cell, mols, charges, out);
}

bool axisOrderParameters(const double *vectors, std::size_t nVectors,
                         const double axis[3], OrientationOrderResult &out) {
  out = OrientationOrderResult{};
  if (axis == nullptr || !finiteArray(axis, 3)) {
    return false;
  }
  const double aNorm =
      std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  if (!(aNorm > 0.0) || !std::isfinite(aNorm)) {
    return false;
  }
  if (nVectors == 0) {
    return true;
  }
  if (vectors == nullptr || !finiteArray(vectors, 3 * nVectors)) {
    return false;
  }

  const double ex = axis[0] / aNorm;
  const double ey = axis[1] / aNorm;
  const double ez = axis[2] / aNorm;
  const double nan = std::numeric_limits<double>::quiet_NaN();

  out.cosTheta.assign(nVectors, nan);
  out.p2.assign(nVectors, nan);
  double sumP1 = 0.0;
  double sumP2 = 0.0;
  std::size_t validCount = 0;

  for (std::size_t i = 0; i < nVectors; ++i) {
    const double vx = vectors[3 * i];
    const double vy = vectors[3 * i + 1];
    const double vz = vectors[3 * i + 2];
    const double vNorm = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (!(vNorm > 1e-15)) {
      continue;
    }
    const double c =
        std::clamp((vx * ex + vy * ey + vz * ez) / vNorm, -1.0, 1.0);
    const double p2Val = 0.5 * (3.0 * c * c - 1.0);
    out.cosTheta[i] = c;
    out.p2[i] = p2Val;
    sumP1 += c;
    sumP2 += p2Val;
    ++validCount;
  }

  out.nValid = validCount;
  if (validCount > 0) {
    const double inv = 1.0 / static_cast<double>(validCount);
    out.meanP1 = sumP1 * inv;
    out.meanP2 = sumP2 * inv;
  } else {
    out.meanP1 = nan;
    out.meanP2 = nan;
  }
  return true;
}

bool fitDiffusionCurve(const std::vector<double> &lagTime,
                       const std::vector<double> &msdTotal,
                       const std::vector<double> &msdAxial,
                       const std::vector<double> &msdLateral, double tStart,
                       double tEnd, DiffusionFitResult &out) {
  out = DiffusionFitResult{};
  const std::size_t n = lagTime.size();
  if (n < 2 || msdTotal.size() != n || msdAxial.size() != n ||
      msdLateral.size() != n) {
    return false;
  }
  double lo = 0.0;
  double hi = 0.0;
  if (!resolveTimeWindow(lagTime, tStart, tEnd, lo, hi)) {
    return false;
  }

  std::vector<double> tx;
  std::vector<double> yTot;
  std::vector<double> yAx;
  std::vector<double> yLat;
  std::vector<double> ones;
  std::vector<double> logT;
  std::vector<double> logMsd;
  std::vector<double> logW;

  for (std::size_t i = 0; i < n; ++i) {
    const double t = lagTime[i];
    if (t < lo - 1e-12 * hi || t > hi + 1e-12 * hi) {
      continue;
    }
    if (!std::isfinite(t) || !std::isfinite(msdTotal[i]) ||
        !std::isfinite(msdAxial[i]) || !std::isfinite(msdLateral[i])) {
      continue;
    }
    tx.push_back(t);
    yTot.push_back(msdTotal[i]);
    yAx.push_back(msdAxial[i]);
    yLat.push_back(msdLateral[i]);
    ones.push_back(1.0);
    if (t > 0.0 && msdTotal[i] > 1e-20) {
      logT.push_back(std::log(t));
      logMsd.push_back(std::log(msdTotal[i]));
      logW.push_back(1.0);
    }
  }

  if (tx.size() < 2) {
    return false;
  }
  const auto regTot = fitWeightedLine(tx, yTot, ones);
  const auto regAx = fitWeightedLine(tx, yAx, ones);
  const auto regLat = fitWeightedLine(tx, yLat, ones);
  if (!regTot.ok || !regAx.ok || !regLat.ok) {
    return false;
  }
  out.dTotal = regTot.slope / 6.0;
  out.dAxial = regAx.slope / 2.0;
  out.dLateral = regLat.slope / 4.0;
  out.r2Total = regTot.r2;

  if (logT.size() >= 2) {
    const auto regAlpha = fitWeightedLine(logT, logMsd, logW);
    if (regAlpha.ok) {
      out.alpha = regAlpha.slope;
    }
  }
  return true;
}

bool fitCorrelationDecay(const std::vector<double> &lagTime,
                         const std::vector<double> &corr, double tStart,
                         double tEnd, DecayFitResult &out) {
  out = DecayFitResult{};
  const std::size_t n = lagTime.size();
  if (n < 2 || corr.size() != n) {
    return false;
  }
  double lo = 0.0;
  double hi = 0.0;
  if (tStart == 0.0 && tEnd == 0.0) {
    lo = 0.0;
    hi = lagTime.back();
    if (!(hi > 0.0) || !std::isfinite(hi)) {
      return false;
    }
  } else if (!resolveTimeWindow(lagTime, tStart, tEnd, lo, hi)) {
    return false;
  }

  // Trapezoidal integral up to first non-positive crossing.
  double integral = 0.0;
  for (std::size_t i = 1; i < n; ++i) {
    if (!std::isfinite(lagTime[i - 1]) || !std::isfinite(lagTime[i]) ||
        !std::isfinite(corr[i - 1]) || !std::isfinite(corr[i])) {
      break;
    }
    const double dtStep = lagTime[i] - lagTime[i - 1];
    if (!(dtStep > 0.0)) {
      break;
    }
    if (corr[i - 1] <= 0.0) {
      break;
    }
    if (corr[i] <= 0.0) {
      const double frac = corr[i - 1] / (corr[i - 1] - corr[i]);
      integral += 0.5 * corr[i - 1] * (dtStep * frac);
      break;
    }
    integral += 0.5 * (corr[i - 1] + corr[i]) * dtStep;
  }
  out.tauIntegral = integral;

  std::vector<double> tx;
  std::vector<double> logC;
  std::vector<double> weights;
  for (std::size_t i = 0; i < n; ++i) {
    const double t = lagTime[i];
    const double c = corr[i];
    if (t < lo - 1e-12 * std::max(1.0, hi) ||
        t > hi + 1e-12 * std::max(1.0, hi)) {
      continue;
    }
    if (!std::isfinite(t) || !std::isfinite(c)) {
      continue;
    }
    if (c <= 1e-6) {
      if (t > lo) {
        break;
      }
      continue;
    }
    tx.push_back(t);
    logC.push_back(std::log(c));
    weights.push_back(c * c);
  }

  const double nan = std::numeric_limits<double>::quiet_NaN();
  if (tx.size() >= 2) {
    const auto reg = fitWeightedLine(tx, logC, weights);
    if (reg.ok && reg.slope < -1e-15) {
      out.tauExp = -1.0 / reg.slope;
      out.amplitude = std::exp(reg.intercept);
      out.r2 = reg.r2;
    } else {
      out.tauExp = nan;
      out.amplitude = nan;
      out.r2 = reg.r2;
    }
  } else {
    out.tauExp = nan;
    out.amplitude = nan;
    out.r2 = 0.0;
  }
  return true;
}

MsdAccumulator::MsdAccumulator(std::size_t nParticles, std::size_t maxLag,
                               std::size_t originStride, double dt,
                               const double axis[3])
    : nParticles_(nParticles), maxLag_(maxLag), originStride_(originStride),
      dt_(dt), unwrapper_(nParticles) {
  if (nParticles == 0 || originStride == 0 || !(dt > 0.0) ||
      !std::isfinite(dt) || axis == nullptr || !finiteArray(axis, 3)) {
    return;
  }
  const double aNorm =
      std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  if (!(aNorm > 0.0) || !std::isfinite(aNorm)) {
    return;
  }
  axisUnit_ = {axis[0] / aNorm, axis[1] / aNorm, axis[2] / aNorm};
  scratchUnwrapped_.resize(3 * nParticles, 0.0);
  if (maxLag_ > 0) {
    sumTotal_.assign(maxLag_ + 1, 0.0);
    sumAxial_.assign(maxLag_ + 1, 0.0);
    sumLateral_.assign(maxLag_ + 1, 0.0);
    originCounts_.assign(maxLag_ + 1, 0);
  }
  valid_ = true;
}

std::size_t MsdAccumulator::numLags() const noexcept {
  if (!valid_ || framesSeen_ == 0) {
    return 0;
  }
  if (maxLag_ > 0) {
    return std::min(framesSeen_, maxLag_ + 1);
  }
  return framesSeen_;
}

bool MsdAccumulator::push(const double *xyz, std::size_t nParticles,
                          const double *box9OrNull) {
  if (!valid_ || xyz == nullptr || nParticles != nParticles_) {
    return false;
  }
  if (!finiteArray(xyz, 3 * nParticles_)) {
    return false;
  }

  const double *coords = xyz;
  if (box9OrNull != nullptr) {
    if (framesSeen_ > 0 && !usedBox_) {
      return false;
    }
    const BoxCell cell = BoxCell::fromBox9(box9OrNull);
    if (!unwrapper_.pushFrame(xyz, nParticles_, cell,
                              scratchUnwrapped_.data())) {
      return false;
    }
    usedBox_ = true;
    coords = scratchUnwrapped_.data();
  } else if (framesSeen_ > 0 && usedBox_) {
    return false;
  }

  const std::size_t curFrame = framesSeen_;
  if (maxLag_ == 0) {
    sumTotal_.resize(curFrame + 1, 0.0);
    sumAxial_.resize(curFrame + 1, 0.0);
    sumLateral_.resize(curFrame + 1, 0.0);
    originCounts_.resize(curFrame + 1, 0);
  } else {
    while (!activeOrigins_.empty() &&
           curFrame - activeOrigins_.front().frameIndex > maxLag_) {
      activeOrigins_.pop_front();
    }
  }

  const double ex = axisUnit_[0];
  const double ey = axisUnit_[1];
  const double ez = axisUnit_[2];

  for (const auto &origin : activeOrigins_) {
    const std::size_t lag = curFrame - origin.frameIndex;
    if (maxLag_ > 0 && lag > maxLag_) {
      continue;
    }
    double tot = 0.0;
    double ax = 0.0;
    double lat = 0.0;
    const double *orig = origin.coords.data();
    for (std::size_t i = 0; i < nParticles_; ++i) {
      const std::size_t base = 3 * i;
      const double dx = coords[base] - orig[base];
      const double dy = coords[base + 1] - orig[base + 1];
      const double dz = coords[base + 2] - orig[base + 2];
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double proj = dx * ex + dy * ey + dz * ez;
      const double ax2 = proj * proj;
      const double lat2 = std::max(0.0, r2 - ax2);
      tot += r2;
      ax += ax2;
      lat += lat2;
    }
    sumTotal_[lag] += tot;
    sumAxial_[lag] += ax;
    sumLateral_[lag] += lat;
    originCounts_[lag] += 1;
  }

  if (curFrame % originStride_ == 0) {
    OriginSlot slot;
    slot.frameIndex = curFrame;
    slot.coords.assign(coords, coords + 3 * nParticles_);
    activeOrigins_.push_back(std::move(slot));
    originCounts_[0] += 1;
  }

  ++framesSeen_;
  return true;
}

void MsdAccumulator::curve(std::vector<double> &lagTime,
                           std::vector<double> &total,
                           std::vector<double> &axial,
                           std::vector<double> &lateral,
                           std::vector<std::size_t> &nOrigins) const {
  const std::size_t m = numLags();
  lagTime.resize(m);
  total.resize(m);
  axial.resize(m);
  lateral.resize(m);
  nOrigins.resize(m);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (std::size_t k = 0; k < m; ++k) {
    lagTime[k] = static_cast<double>(k) * dt_;
    nOrigins[k] = originCounts_[k];
    if (k == 0) {
      total[k] = originCounts_[0] > 0 ? 0.0 : nan;
      axial[k] = originCounts_[0] > 0 ? 0.0 : nan;
      lateral[k] = originCounts_[0] > 0 ? 0.0 : nan;
    } else if (originCounts_[k] > 0) {
      const double inv =
          1.0 / (static_cast<double>(originCounts_[k]) *
                 static_cast<double>(nParticles_));
      total[k] = sumTotal_[k] * inv;
      axial[k] = sumAxial_[k] * inv;
      lateral[k] = sumLateral_[k] * inv;
    } else {
      total[k] = nan;
      axial[k] = nan;
      lateral[k] = nan;
    }
  }
}

bool MsdAccumulator::fit(double tStart, double tEnd,
                         DiffusionFitResult &out) const {
  std::vector<double> lag;
  std::vector<double> tot;
  std::vector<double> ax;
  std::vector<double> lat;
  std::vector<std::size_t> counts;
  curve(lag, tot, ax, lat, counts);
  return fitDiffusionCurve(lag, tot, ax, lat, tStart, tEnd, out);
}

HbondLifetimeAccumulator::HbondLifetimeAccumulator(std::size_t maxLag,
                                                   std::size_t originStride,
                                                   double dt, bool directed)
    : maxLag_(maxLag), originStride_(originStride), dt_(dt),
      directed_(directed) {
  if (originStride == 0 || !(dt > 0.0) || !std::isfinite(dt)) {
    return;
  }
  if (maxLag_ > 0) {
    intermittentNum_.assign(maxLag_ + 1, 0.0);
    continuousNum_.assign(maxLag_ + 1, 0.0);
    denom_.assign(maxLag_ + 1, 0.0);
    originCounts_.assign(maxLag_ + 1, 0);
  }
  valid_ = true;
}

std::size_t HbondLifetimeAccumulator::numLags() const noexcept {
  if (!valid_ || framesSeen_ == 0) {
    return 0;
  }
  if (maxLag_ > 0) {
    return std::min(framesSeen_, maxLag_ + 1);
  }
  return framesSeen_;
}

bool HbondLifetimeAccumulator::push(const int *donors, const int *acceptors,
                                    std::size_t nEdges) {
  if (!valid_) {
    return false;
  }
  if (nEdges > 0 && (donors == nullptr || acceptors == nullptr)) {
    return false;
  }

  std::vector<std::uint64_t> curEdges;
  curEdges.reserve(nEdges);
  for (std::size_t i = 0; i < nEdges; ++i) {
    curEdges.push_back(packEdge(donors[i], acceptors[i], directed_));
  }
  std::sort(curEdges.begin(), curEdges.end());
  curEdges.erase(std::unique(curEdges.begin(), curEdges.end()), curEdges.end());

  const std::size_t curFrame = framesSeen_;
  if (maxLag_ == 0) {
    intermittentNum_.resize(curFrame + 1, 0.0);
    continuousNum_.resize(curFrame + 1, 0.0);
    denom_.resize(curFrame + 1, 0.0);
    originCounts_.resize(curFrame + 1, 0);
  } else {
    while (!activeOrigins_.empty() &&
           curFrame - activeOrigins_.front().frameIndex > maxLag_) {
      activeOrigins_.pop_front();
    }
  }

  for (auto &origin : activeOrigins_) {
    const std::size_t lag = curFrame - origin.frameIndex;
    if (maxLag_ > 0 && lag > maxLag_) {
      continue;
    }
    const std::size_t nInter =
        sortedIntersectionCount(origin.initialEdges, curEdges);
    inplaceSortedIntersection(origin.survivingEdges, curEdges);
    const std::size_t nCont = origin.survivingEdges.size();

    intermittentNum_[lag] += static_cast<double>(nInter);
    continuousNum_[lag] += static_cast<double>(nCont);
    denom_[lag] += static_cast<double>(origin.initialEdges.size());
    originCounts_[lag] += 1;
  }

  if (curFrame % originStride_ == 0) {
    OriginSlot slot;
    slot.frameIndex = curFrame;
    slot.initialEdges = curEdges;
    slot.survivingEdges = std::move(curEdges);
    const auto nInit = static_cast<double>(slot.initialEdges.size());
    intermittentNum_[0] += nInit;
    continuousNum_[0] += nInit;
    denom_[0] += nInit;
    originCounts_[0] += 1;
    activeOrigins_.push_back(std::move(slot));
  }

  ++framesSeen_;
  return true;
}

void HbondLifetimeAccumulator::curve(std::vector<double> &lagTime,
                                     std::vector<double> &cIntermittent,
                                     std::vector<double> &sContinuous,
                                     std::vector<std::size_t> &nOrigins) const {
  const std::size_t m = numLags();
  lagTime.resize(m);
  cIntermittent.resize(m);
  sContinuous.resize(m);
  nOrigins.resize(m);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (std::size_t k = 0; k < m; ++k) {
    lagTime[k] = static_cast<double>(k) * dt_;
    nOrigins[k] = originCounts_[k];
    if (denom_[k] > 0.0) {
      cIntermittent[k] = intermittentNum_[k] / denom_[k];
      sContinuous[k] = continuousNum_[k] / denom_[k];
    } else {
      cIntermittent[k] = nan;
      sContinuous[k] = nan;
    }
  }
}

bool HbondLifetimeAccumulator::fit(double tStart, double tEnd,
                                   DecayFitResult &intermittentFit,
                                   DecayFitResult &continuousFit) const {
  std::vector<double> lag;
  std::vector<double> cHb;
  std::vector<double> sHb;
  std::vector<std::size_t> counts;
  curve(lag, cHb, sHb, counts);
  if (!fitCorrelationDecay(lag, cHb, tStart, tEnd, intermittentFit)) {
    return false;
  }
  return fitCorrelationDecay(lag, sHb, tStart, tEnd, continuousFit);
}

ReorientationAccumulator::ReorientationAccumulator(std::size_t nVectors,
                                                   std::size_t maxLag,
                                                   std::size_t originStride,
                                                   double dt)
    : nVectors_(nVectors), maxLag_(maxLag), originStride_(originStride),
      dt_(dt) {
  if (nVectors == 0 || originStride == 0 || !(dt > 0.0) || !std::isfinite(dt)) {
    return;
  }
  if (maxLag_ > 0) {
    sumC1_.assign(maxLag_ + 1, 0.0);
    sumC2_.assign(maxLag_ + 1, 0.0);
    pairCounts_.assign(maxLag_ + 1, 0);
    originCounts_.assign(maxLag_ + 1, 0);
  }
  valid_ = true;
}

std::size_t ReorientationAccumulator::numLags() const noexcept {
  if (!valid_ || framesSeen_ == 0) {
    return 0;
  }
  if (maxLag_ > 0) {
    return std::min(framesSeen_, maxLag_ + 1);
  }
  return framesSeen_;
}

bool ReorientationAccumulator::push(const double *vectors,
                                    std::size_t nVectors) {
  if (!valid_ || vectors == nullptr || nVectors != nVectors_) {
    return false;
  }
  if (!finiteArray(vectors, 3 * nVectors_)) {
    return false;
  }

  std::vector<double> units(3 * nVectors_, 0.0);
  std::vector<std::uint8_t> validMask(nVectors_, 0);
  std::size_t nValidNow = 0;
  for (std::size_t i = 0; i < nVectors_; ++i) {
    const double vx = vectors[3 * i];
    const double vy = vectors[3 * i + 1];
    const double vz = vectors[3 * i + 2];
    const double norm = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (norm > 1e-15) {
      units[3 * i] = vx / norm;
      units[3 * i + 1] = vy / norm;
      units[3 * i + 2] = vz / norm;
      validMask[i] = 1;
      ++nValidNow;
    }
  }

  const std::size_t curFrame = framesSeen_;
  if (maxLag_ == 0) {
    sumC1_.resize(curFrame + 1, 0.0);
    sumC2_.resize(curFrame + 1, 0.0);
    pairCounts_.resize(curFrame + 1, 0);
    originCounts_.resize(curFrame + 1, 0);
  } else {
    while (!activeOrigins_.empty() &&
           curFrame - activeOrigins_.front().frameIndex > maxLag_) {
      activeOrigins_.pop_front();
    }
  }

  for (const auto &origin : activeOrigins_) {
    const std::size_t lag = curFrame - origin.frameIndex;
    if (maxLag_ > 0 && lag > maxLag_) {
      continue;
    }
    double c1Acc = 0.0;
    double c2Acc = 0.0;
    std::size_t nPairs = 0;
    for (std::size_t i = 0; i < nVectors_; ++i) {
      if (origin.mask[i] == 0 || validMask[i] == 0) {
        continue;
      }
      const double dot = std::clamp(
          origin.unitVectors[3 * i] * units[3 * i] +
              origin.unitVectors[3 * i + 1] * units[3 * i + 1] +
              origin.unitVectors[3 * i + 2] * units[3 * i + 2],
          -1.0, 1.0);
      c1Acc += dot;
      c2Acc += 0.5 * (3.0 * dot * dot - 1.0);
      ++nPairs;
    }
    sumC1_[lag] += c1Acc;
    sumC2_[lag] += c2Acc;
    pairCounts_[lag] += nPairs;
    originCounts_[lag] += 1;
  }

  if (curFrame % originStride_ == 0) {
    sumC1_[0] += static_cast<double>(nValidNow);
    sumC2_[0] += static_cast<double>(nValidNow);
    pairCounts_[0] += nValidNow;
    originCounts_[0] += 1;
    OriginSlot slot;
    slot.frameIndex = curFrame;
    slot.unitVectors = std::move(units);
    slot.mask = std::move(validMask);
    activeOrigins_.push_back(std::move(slot));
  }

  ++framesSeen_;
  return true;
}

void ReorientationAccumulator::curve(std::vector<double> &lagTime,
                                     std::vector<double> &c1,
                                     std::vector<double> &c2,
                                     std::vector<std::size_t> &nOrigins) const {
  const std::size_t m = numLags();
  lagTime.resize(m);
  c1.resize(m);
  c2.resize(m);
  nOrigins.resize(m);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (std::size_t k = 0; k < m; ++k) {
    lagTime[k] = static_cast<double>(k) * dt_;
    nOrigins[k] = originCounts_[k];
    if (pairCounts_[k] > 0) {
      const double inv = 1.0 / static_cast<double>(pairCounts_[k]);
      c1[k] = sumC1_[k] * inv;
      c2[k] = sumC2_[k] * inv;
    } else {
      c1[k] = nan;
      c2[k] = nan;
    }
  }
}

bool ReorientationAccumulator::fit(double tStart, double tEnd,
                                   DecayFitResult &c1Fit,
                                   DecayFitResult &c2Fit) const {
  std::vector<double> lag;
  std::vector<double> c1Vec;
  std::vector<double> c2Vec;
  std::vector<std::size_t> counts;
  curve(lag, c1Vec, c2Vec, counts);
  if (!fitCorrelationDecay(lag, c1Vec, tStart, tEnd, c1Fit)) {
    return false;
  }
  return fitCorrelationDecay(lag, c2Vec, tStart, tEnd, c2Fit);
}

} // namespace traj
