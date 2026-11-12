//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//-----------------------------------------------------------------------------------

#include <interface_order.hpp>

#include <bop.hpp>
#include <generic.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <vector>

namespace iface {
namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

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

std::uint64_t xorshift(std::uint64_t &state) {
  state ^= state << 13;
  state ^= state >> 7;
  state ^= state << 17;
  return state;
}

double tenNinety(const std::vector<double> &coord, const std::vector<double> &value) {
  if (coord.size() < 4 || coord.size() != value.size()) {
    return 0.0;
  }
  int imax = 0;
  for (int i = 1; i < static_cast<int>(value.size()); ++i) {
    if (value[static_cast<std::size_t>(i)] > value[static_cast<std::size_t>(imax)]) {
      imax = i;
    }
  }
  const double hi = value[static_cast<std::size_t>(imax)];
  double tail = 0.0;
  int nTail = 0;
  const int tail0 = static_cast<int>(value.size() * 3 / 4);
  for (int i = tail0; i < static_cast<int>(value.size()); ++i) {
    tail += value[static_cast<std::size_t>(i)];
    ++nTail;
  }
  const double lo = nTail > 0 ? tail / static_cast<double>(nTail) : value.back();
  if (!(hi > lo)) {
    return 0.0;
  }
  const double v90 = lo + 0.9 * (hi - lo);
  const double v10 = lo + 0.1 * (hi - lo);
  bool have90 = false;
  bool have10 = false;
  double s90 = 0.0;
  double s10 = 0.0;
  double prevS = coord[static_cast<std::size_t>(imax)];
  double prevV = hi;
  for (int i = imax + 1; i < static_cast<int>(value.size()); ++i) {
    const double s = coord[static_cast<std::size_t>(i)];
    const double v = value[static_cast<std::size_t>(i)];
    if (!have90 && prevV >= v90 && v < v90 && v != prevV) {
      const double t = (v90 - prevV) / (v - prevV);
      s90 = prevS + t * (s - prevS);
      have90 = true;
    }
    if (!have10 && prevV >= v10 && v < v10 && v != prevV) {
      const double t = (v10 - prevV) / (v - prevV);
      s10 = prevS + t * (s - prevS);
      have10 = true;
    }
    prevS = s;
    prevV = v;
  }
  if (have90 && have10) {
    return std::abs(s10 - s90);
  }
  return 0.0;
}

double gauss(double d2, double xi) {
  const double s2 = xi * xi;
  const double norm = std::pow(2.0 * 3.14159265358979323846 * s2, -1.5);
  return norm * std::exp(-0.5 * d2 / s2);
}

bool solidLabel(int label) {
  return label == static_cast<int>(molSys::atom_state_type::cubic) ||
         label == static_cast<int>(molSys::atom_state_type::hexagonal) ||
         label == static_cast<int>(molSys::atom_state_type::interfacial);
}

} // namespace

double gaussianStepWidth(double xi) { return 2.0 * 1.2815515655446004 * xi; }

BootstrapInterval bootstrapMean(const std::vector<double> &samples, int nDraw, std::uint64_t seed) {
  BootstrapInterval out;
  if (samples.empty() || nDraw <= 0) {
    return out;
  }
  double sum = 0.0;
  for (double v : samples) {
    sum += v;
  }
  out.mean = sum / static_cast<double>(samples.size());
  std::uint64_t state = seed == 0 ? 1 : seed;
  std::vector<double> means;
  means.reserve(static_cast<std::size_t>(nDraw));
  const int n = static_cast<int>(samples.size());
  for (int draw = 0; draw < nDraw; ++draw) {
    double acc = 0.0;
    for (int i = 0; i < n; ++i) {
      const auto pick = static_cast<std::size_t>(xorshift(state) % static_cast<std::uint64_t>(n));
      acc += samples[pick];
    }
    means.push_back(acc / static_cast<double>(n));
  }
  std::sort(means.begin(), means.end());
  out.lo = quantile7(means, 0.025);
  out.hi = quantile7(means, 0.975);
  out.excludesZero = out.lo > 0.0 || out.hi < 0.0;
  return out;
}

BootstrapInterval pairedDifference(const std::vector<double> &a, const std::vector<double> &b,
                                   int nDraw, std::uint64_t seed) {
  const int n = static_cast<int>(std::min(a.size(), b.size()));
  std::vector<double> diff(static_cast<std::size_t>(std::max(n, 0)));
  for (int i = 0; i < n; ++i) {
    diff[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i)] - b[static_cast<std::size_t>(i)];
  }
  return bootstrapMean(diff, nDraw, seed);
}

WillardChandler willardChandler(Cloud cloud, double xi) {
  WillardChandler surface;
  surface.xi = xi;
  if (cloud.nop < 16 || cloud.box.size() < 3 || !(xi > 0.0)) {
    return surface;
  }
  const auto rows = nneigh::neighListO(3.5, cloud, 1);
  Cloud published = cloud;
  chill::getCorrelPlus(published, rows, false, 4, true);
  chill::getIceTypePlusNoPrint(published, rows, false, true);
  Cloud ovito = cloud;
  chill::getCorrelPlus(ovito, rows, false, 0, false);
  chill::getIceTypePlusNoPrint(ovito, rows, false, false);

  int anchor = -1;
  for (int i = 0; i < cloud.nop; ++i) {
    const int label = static_cast<int>(published.pts[static_cast<std::size_t>(i)].iceType);
    if (label == static_cast<int>(molSys::atom_state_type::cubic) ||
        label == static_cast<int>(molSys::atom_state_type::hexagonal)) {
      anchor = i;
      break;
    }
  }
  if (anchor < 0) {
    return surface;
  }
  std::array<double, 3> acc{0.0, 0.0, 0.0};
  int nCrystal = 0;
  for (int i = 0; i < cloud.nop; ++i) {
    const int label = static_cast<int>(published.pts[static_cast<std::size_t>(i)].iceType);
    if (label != static_cast<int>(molSys::atom_state_type::cubic) &&
        label != static_cast<int>(molSys::atom_state_type::hexagonal)) {
      continue;
    }
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    const auto dr = gen::relDistFromPoint(cloud, anchor, p.x, p.y, p.z);
    acc[0] += dr[0];
    acc[1] += dr[1];
    acc[2] += dr[2];
    ++nCrystal;
  }
  if (nCrystal < 8) {
    return surface;
  }
  const double cx = cloud.pts[static_cast<std::size_t>(anchor)].x - acc[0] / nCrystal;
  const double cy = cloud.pts[static_cast<std::size_t>(anchor)].y - acc[1] / nCrystal;
  const double cz = cloud.pts[static_cast<std::size_t>(anchor)].z - acc[2] / nCrystal;
  const double rMax =
      0.5 * std::min(cloud.box[0], std::min(cloud.box[1], cloud.box[2]));
  const double cut2 = (3.0 * xi) * (3.0 * xi);
  const int nTheta = 12;
  const int nPhi = 24;
  const double pi = 3.14159265358979323846;
  const double dr = 0.5;
  const int nR = std::max(4, static_cast<int>(rMax / dr));

  auto densityAt = [&](double x, double y, double z) {
    double rho = 0.0;
    for (int i = 0; i < cloud.nop; ++i) {
      const auto drv = gen::relDistFromPoint(cloud, i, x, y, z);
      const double d2 = drv[0] * drv[0] + drv[1] * drv[1] + drv[2] * drv[2];
      if (d2 < cut2) {
        rho += gauss(d2, xi);
      }
    }
    return rho;
  };

  double rhoIn = 0.0;
  double rhoOut = 0.0;
  int nIn = 0;
  int nOut = 0;
  const int nRay = nTheta * nPhi;
  std::vector<double> radius(static_cast<std::size_t>(nRay), -1.0);
  std::vector<std::array<double, 3>> direction(radius.size());
  std::vector<std::vector<double>> along(static_cast<std::size_t>(nRay),
                                         std::vector<double>(static_cast<std::size_t>(nR), 0.0));
  int ray = 0;
  for (int it = 0; it < nTheta; ++it) {
    const double theta = (static_cast<double>(it) + 0.5) * pi / nTheta;
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    for (int ip = 0; ip < nPhi; ++ip, ++ray) {
      const double phi = (static_cast<double>(ip) + 0.5) * 2.0 * pi / nPhi;
      const double ux = st * std::cos(phi);
      const double uy = st * std::sin(phi);
      const double uz = ct;
      direction[static_cast<std::size_t>(ray)] = {ux, uy, uz};
      for (int ir = 1; ir < nR; ++ir) {
        const double r = ir * dr;
        const double rho = densityAt(cx + ux * r, cy + uy * r, cz + uz * r);
        along[static_cast<std::size_t>(ray)][static_cast<std::size_t>(ir)] = rho;
        if (r < 4.0) {
          rhoIn += rho;
          ++nIn;
        }
        if (r > 0.4 * rMax) {
          rhoOut += rho;
          ++nOut;
        }
      }
    }
  }
  if (nIn == 0 || nOut == 0) {
    return surface;
  }
  rhoIn /= static_cast<double>(nIn);
  rhoOut /= static_cast<double>(nOut);
  const double threshold = 0.5 * (rhoIn + rhoOut);
  if (!(rhoIn > rhoOut)) {
    return surface;
  }
  ray = 0;
  std::vector<double> kept;
  for (int it = 0; it < nTheta; ++it) {
    for (int ip = 0; ip < nPhi; ++ip, ++ray) {
      double prev = along[static_cast<std::size_t>(ray)][1];
      double prevR = dr;
      for (int ir = 2; ir < nR; ++ir) {
        const double r = ir * dr;
        const double rho = along[static_cast<std::size_t>(ray)][static_cast<std::size_t>(ir)];
        if (prev >= threshold && rho < threshold && rho != prev) {
          const double t = (threshold - prev) / (rho - prev);
          const double rc = prevR + t * (r - prevR);
          radius[static_cast<std::size_t>(ray)] = rc;
          kept.push_back(rc);
          break;
        }
        prev = rho;
        prevR = r;
      }
    }
  }
  if (kept.size() < 8) {
    return surface;
  }
  double mean = 0.0;
  for (double r : kept) {
    mean += r;
  }
  mean /= static_cast<double>(kept.size());
  double var = 0.0;
  for (double r : kept) {
    const double d = r - mean;
    var += d * d;
  }
  surface.meanRadius = mean;
  surface.surfaceRms = std::sqrt(var / static_cast<double>(kept.size()));
  surface.nRays = static_cast<int>(kept.size());

  const int nBin = 48;
  const double d0 = -8.0;
  const double bw = 0.5;
  std::vector<int> count(static_cast<std::size_t>(nBin), 0);
  std::vector<int> solidP(count.size(), 0);
  std::vector<int> solidO(count.size(), 0);
  for (int i = 0; i < cloud.nop; ++i) {
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    const auto drv = gen::relDistFromPoint(cloud, i, cx, cy, cz);
    const double vx = drv[0];
    const double vy = drv[1];
    const double vz = drv[2];
    const double r = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (!(r > 1.0e-8)) {
      continue;
    }
    int best = 0;
    double bestDot = -2.0;
    for (int k = 0; k < nTheta * nPhi; ++k) {
      if (radius[static_cast<std::size_t>(k)] < 0.0) {
        continue;
      }
      const auto &u = direction[static_cast<std::size_t>(k)];
      const double dot = (vx * u[0] + vy * u[1] + vz * u[2]) / r;
      if (dot > bestDot) {
        bestDot = dot;
        best = k;
      }
    }
    const double delta = r - radius[static_cast<std::size_t>(best)];
    int bin = static_cast<int>((delta - d0) / bw);
    if (bin < 0 || bin >= nBin) {
      continue;
    }
    ++count[static_cast<std::size_t>(bin)];
    const int lp = static_cast<int>(published.pts[static_cast<std::size_t>(i)].iceType);
    const int lo = static_cast<int>(ovito.pts[static_cast<std::size_t>(i)].iceType);
    solidP[static_cast<std::size_t>(bin)] += solidLabel(lp) ? 1 : 0;
    solidO[static_cast<std::size_t>(bin)] += solidLabel(lo) ? 1 : 0;
  }
  std::vector<double> coord(static_cast<std::size_t>(nBin));
  std::vector<double> dens(static_cast<std::size_t>(nBin), 0.0);
  std::vector<double> fracP(static_cast<std::size_t>(nBin), 0.0);
  std::vector<double> fracO(static_cast<std::size_t>(nBin), 0.0);
  for (int b = 0; b < nBin; ++b) {
    coord[static_cast<std::size_t>(b)] = d0 + (static_cast<double>(b) + 0.5) * bw;
    dens[static_cast<std::size_t>(b)] = static_cast<double>(count[static_cast<std::size_t>(b)]);
    if (count[static_cast<std::size_t>(b)] > 0) {
      fracP[static_cast<std::size_t>(b)] =
          static_cast<double>(solidP[static_cast<std::size_t>(b)]) / count[static_cast<std::size_t>(b)];
      fracO[static_cast<std::size_t>(b)] =
          static_cast<double>(solidO[static_cast<std::size_t>(b)]) / count[static_cast<std::size_t>(b)];
    }
  }
  std::vector<double> fieldX;
  std::vector<double> fieldY;
  for (double delta = -8.0; delta <= 12.0; delta += 0.5) {
    double acc = 0.0;
    int nField = 0;
    for (int k = 0; k < nRay; ++k) {
      if (radius[static_cast<std::size_t>(k)] < 0.0) {
        continue;
      }
      const double r = radius[static_cast<std::size_t>(k)] + delta;
      const double idx = r / dr;
      const int lo = static_cast<int>(idx);
      if (lo < 1 || lo + 1 >= nR) {
        continue;
      }
      const double frac = idx - static_cast<double>(lo);
      const auto &row = along[static_cast<std::size_t>(k)];
      acc += row[static_cast<std::size_t>(lo)] * (1.0 - frac) +
             row[static_cast<std::size_t>(lo + 1)] * frac;
      ++nField;
    }
    fieldX.push_back(delta);
    fieldY.push_back(nField > 0 ? acc / static_cast<double>(nField) : 0.0);
  }
  surface.densityWidth = tenNinety(fieldX, fieldY);
  (void)dens;
  surface.solidWidthPub = tenNinety(coord, fracP);
  surface.solidWidthOvi = tenNinety(coord, fracO);
  surface.defined = surface.densityWidth > 0.0;
  return surface;
}

std::string formatWillardChandler(const WillardChandler &surface) {
  std::ostringstream os;
  os.setf(std::ios::fixed);
  os.precision(6);
  os << "# willard-chandler\n";
  os << "xi " << surface.xi << "\n";
  os << "defined " << (surface.defined ? 1 : 0) << "\n";
  os << "nrays " << surface.nRays << "\n";
  os << "mean_radius " << surface.meanRadius << "\n";
  os << "surface_rms " << surface.surfaceRms << "\n";
  os << "density_width " << surface.densityWidth << "\n";
  os << "solid_width_pub " << surface.solidWidthPub << "\n";
  os << "solid_width_ovi " << surface.solidWidthOvi << "\n";
  return os.str();
}

} // namespace iface
