//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <order_parameter.hpp>
#include <generic.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_map>

/**
 * @details The average height of prism blocks remains relatively constant. We
 * have observed a average prism heights of 2.7-2.85 Angstrom for prisms
 * irrespective of the number of nodes. The equation is given by:
 *
 *  @f[
 *  Height_{n}% = \frac{N_n}{N_{max}} \times 100
 *  @f]
 *
 * Here, @f$N_{max} = H_{SWCT}/h_{avg}f$ and @f$N_{n}$ is the
 * number of prism blocks for n-sided prismatic phase.
 *
 * This means that the normalization factor, is the same for
 * every node number @f$n@f$.
 */
double topoparam::normHeightPercent(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int nPrisms,
    double avgPrismHeight) {
  //
  double hPercent;       // Normalized height percent
  double nanoTubeHeight; // Height of the SWCT
  double numberMax; // Maximum number possible, given the average prism height

  // ---------------------------------------
  // Calculate the height of the SWCT
  // This is the longest recovered cell length
  double lengths[3];
  nneigh::dumpCellLengths(yCloud.box, yCloud.boxLow, lengths);
  nanoTubeHeight = std::max(lengths[0], std::max(lengths[1], lengths[2]));
  // ---------------------------------------
  // Calculate the maximum possible height, given the average prism height
  // and the height of the nanotube
  numberMax = nanoTubeHeight / avgPrismHeight;
  // ---------------------------------------
  // Calculate the normalized height percentage
  hPercent = nPrisms / numberMax * 100.0;

  return hPercent;
}

/**
 * @details The average height of prism blocks remains relatively constant. We
 * have observed a average prism heights of 2.7-2.85 Angstrom for prisms
 * irrespective of the number of nodes. The equation is given by:
 *
 *  @f[
 *  Height_{n}% = \frac{N_n}{N_{max}} \times 100
 *  @f]
 *
 * Here, \f$N_{max} = H_{SWCT}/h_{avg}f$ and \f$N_{n}$ is the
 * number of prism blocks for n-sided prismatic phase.
 *
 * This means that the normalization factor, is the same for
 * every node number \f$n\f$.
 */
std::vector<double> topoparam::calcCoverageArea(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &rings, double sheetArea) {
  //
  double areaXY, areaXZ, areaYZ;   // Total coverage area
  std::vector<double> singleAreas; // Area of single rings

  // ---------------------------------------
  // Initialization
  areaXY = 0.0;
  areaXZ = 0.0;
  areaYZ = 0.0;
  // ---------------------------------------
  // Loop through all the rings
  for (int iring = 0; iring < rings.size(); iring++) {
    // Get the coverage area for the current ring
    singleAreas = topoparam::projAreaSingleRing(yCloud, rings[iring]);
    // Add these to the total coverage area
    areaXY += singleAreas[0];
    areaXZ += singleAreas[1];
    areaYZ += singleAreas[2];
  } // end of loop through all the rings
  // ---------------------------------------
  // Normalize the coverage area by the sheet area
  areaXY = areaXY / sheetArea * 100.0;
  areaXZ = areaXZ / sheetArea * 100.0;
  areaYZ = areaYZ / sheetArea * 100.0;

  return {areaXY, areaXZ, areaYZ};
}

/**
 * @details Calculates the coverage area/ projected area of a single ring
 *  given the ring and the PointCloud.
 */
std::vector<double> topoparam::projAreaSingleRing(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &ring) {
  //
  int iatomIndex, jatomIndex; // Atom indices of the i^th and j^th atoms
  int ringSize = ring.size(); // Number of nodes in the ring
  double areaXY, areaXZ, areaYZ;
  double x_iatom, y_iatom, z_iatom; // Coordinates of iatom
  double x_jatom, y_jatom, z_jatom; // Coordinates of jatom
  // ----------------------------------------
  // Calculate projected area onto the XY, YZ and XZ planes for basal1

  // Init the projected area
  areaXY = 0.0;
  areaXZ = 0.0;
  areaYZ = 0.0;

  jatomIndex = ring[0];

  // All points except the first pair
  for (int k = 1; k < ringSize; k++) {
    iatomIndex = ring[k]; // Current vertex

    // --------------------------------------------------------------------
    // SHIFT PARTICLES TEMPORARILY (IN CASE OF UNWRAPPED COORDINATES)
    gen::unwrappedCoordShift(yCloud, iatomIndex, jatomIndex, &x_iatom, &y_iatom,
                             &z_iatom, &x_jatom, &y_jatom, &z_jatom);
    // --------------------------------------------------------------------

    // Add to the polygon area
    // ------
    // XY plane
    areaXY += (x_jatom + x_iatom) * (y_jatom - y_iatom);
    // ------
    // XZ plane
    areaXZ += (x_jatom + x_iatom) * (z_jatom - z_iatom);
    // ------
    // YZ plane
    areaYZ += (y_jatom + y_iatom) * (z_jatom - z_iatom);
    // ------
    jatomIndex = iatomIndex;
  }

  // Closure point
  iatomIndex = ring[0];
  // Unwrapped coordinates needed
  gen::unwrappedCoordShift(yCloud, iatomIndex, jatomIndex, &x_iatom, &y_iatom,
                           &z_iatom, &x_jatom, &y_jatom, &z_jatom);
  // ------
  // XY plane
  areaXY += (x_jatom + x_iatom) * (y_jatom - y_iatom);
  // ------
  // XZ plane
  areaXZ += (x_jatom + x_iatom) * (z_jatom - z_iatom);
  // ------
  // YZ plane
  areaYZ += (y_jatom + y_iatom) * (z_jatom - z_iatom);
  // ------
  // The actual projected area is half of this
  areaXY *= 0.5;
  areaXZ *= 0.5;
  areaYZ *= 0.5;

  // Absolute area
  areaXY = std::abs(areaXY);
  areaXZ = std::abs(areaXZ);
  areaYZ = std::abs(areaYZ);

  return {areaXY, areaXZ, areaYZ};
}

namespace {

std::array<double, 3> cross3(const std::array<double, 3> &a,
                             const std::array<double, 3> &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

double dot3(const std::array<double, 3> &a, const std::array<double, 3> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm2(const std::array<double, 3> &a) { return dot3(a, a); }

// Dihedral of H1-O1-O2-H2 from MIC displacements O1-H1, O2-O1, H2-O2.
bool dihedralCos3(const std::array<double, 3> &b1, const std::array<double, 3> &b2,
                  const std::array<double, 3> &b3, double &out) {
  const auto n1 = cross3(b1, b2);
  const auto n2 = cross3(b2, b3);
  const double n1s = norm2(n1);
  const double n2s = norm2(n2);
  const double b2s = norm2(b2);
  if (n1s == 0.0 || n2s == 0.0 || b2s == 0.0) {
    return false;
  }
  const auto n1xn2 = cross3(n1, n2);
  const double phi = std::atan2(dot3(n1xn2, b2) / std::sqrt(b2s), dot3(n1, n2));
  out = std::cos(3.0 * phi);
  return true;
}

int outerHydrogen(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                  const std::vector<int> &hs, int otherO) {
  int best = -1;
  double bestD2 = -1.0;
  for (int h : hs) {
    const auto dr = gen::relDist(yCloud, h, otherO);
    const double d2 = norm2(dr);
    if (d2 > bestD2) {
      bestD2 = d2;
      best = h;
    }
  }
  return best;
}

} // namespace

std::vector<double>
topoparam::rodgerF4(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    const std::vector<std::vector<int>> &nList, int oxygenType,
                    int hydrogenType) {
  std::vector<double> out(static_cast<std::size_t>(yCloud.nop),
                          std::numeric_limits<double>::quiet_NaN());
  if (yCloud.nop <= 0) {
    return out;
  }
  std::unordered_map<int, std::vector<int>> hydrogens;
  for (int i = 0; i < yCloud.nop; i++) {
    if (yCloud.pts[static_cast<std::size_t>(i)].type == hydrogenType) {
      hydrogens[yCloud.pts[static_cast<std::size_t>(i)].molID].push_back(i);
    }
  }
  for (int i = 0; i < yCloud.nop; i++) {
    const auto &pi = yCloud.pts[static_cast<std::size_t>(i)];
    if (pi.type != oxygenType) {
      continue;
    }
    const auto hit = hydrogens.find(pi.molID);
    if (hit == hydrogens.end() || hit->second.empty()) {
      continue;
    }
    if (static_cast<std::size_t>(i) >= nList.size() ||
        nList[static_cast<std::size_t>(i)].size() < 2) {
      continue;
    }
    double acc = 0.0;
    int nPair = 0;
    for (std::size_t m = 1; m < nList[static_cast<std::size_t>(i)].size(); m++) {
      const int jid = nList[static_cast<std::size_t>(i)][m];
      const auto it = yCloud.idIndexMap.find(jid);
      if (it == yCloud.idIndexMap.end()) {
        continue;
      }
      const int j = it->second;
      if (j < 0 || j >= yCloud.nop || j == i) {
        continue;
      }
      const auto &pj = yCloud.pts[static_cast<std::size_t>(j)];
      if (pj.type != oxygenType) {
        continue;
      }
      const auto hjt = hydrogens.find(pj.molID);
      if (hjt == hydrogens.end() || hjt->second.empty()) {
        continue;
      }
      const int h1 = outerHydrogen(yCloud, hit->second, j);
      const int h2 = outerHydrogen(yCloud, hjt->second, i);
      if (h1 < 0 || h2 < 0) {
        continue;
      }
      const auto o1h1 = gen::relDist(yCloud, i, h1);
      const auto o2o1 = gen::relDist(yCloud, j, i);
      const auto h2o2 = gen::relDist(yCloud, h2, j);
      double c3 = 0.0;
      if (dihedralCos3(o1h1, o2o1, h2o2, c3)) {
        acc += c3;
        ++nPair;
      }
    }
    if (nPair > 0) {
      out[static_cast<std::size_t>(i)] = acc / static_cast<double>(nPair);
    }
  }
  return out;
}

double topoparam::meanFinite(const std::vector<double> &values) {
  double acc = 0.0;
  int n = 0;
  for (double v : values) {
    if (std::isfinite(v)) {
      acc += v;
      ++n;
    }
  }
  if (n == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return acc / static_cast<double>(n);
}

topoparam::LayerStack
topoparam::layerCubicity(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int axis,
    double layerWidth) {
  LayerStack out;
  if (yCloud.nop <= 0 || axis < 0 || axis > 2 || yCloud.box.size() < 3) {
    return out;
  }
  const double L = yCloud.box[static_cast<std::size_t>(axis)];
  if (!(L > 0.0) || !(layerWidth > 0.0)) {
    return out;
  }
  const int nLayers = std::max(1, static_cast<int>(std::lround(L / layerWidth)));
  const double w = L / static_cast<double>(nLayers);
  out.cubicPerLayer.assign(static_cast<std::size_t>(nLayers), 0);
  out.hexPerLayer.assign(static_cast<std::size_t>(nLayers), 0);
  int nC = 0;
  int nH = 0;
  const double lo = yCloud.boxLow.size() > static_cast<std::size_t>(axis)
                        ? yCloud.boxLow[static_cast<std::size_t>(axis)]
                        : 0.0;
  for (int i = 0; i < yCloud.nop; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    const bool cubic = p.iceType == molSys::atom_state_type::cubic ||
                       p.iceType == molSys::atom_state_type::reCubic;
    const bool hex = p.iceType == molSys::atom_state_type::hexagonal ||
                     p.iceType == molSys::atom_state_type::reHex;
    if (!cubic && !hex) {
      continue;
    }
    const double coord = axis == 0 ? p.x : (axis == 1 ? p.y : p.z);
    double u = coord - lo;
    u -= L * std::floor(u / L);
    if (u < 0.0) {
      u += L;
    }
    if (u >= L) {
      u = 0.0;
    }
    int layer = static_cast<int>(std::floor(u / w));
    if (layer < 0) {
      layer = 0;
    }
    if (layer >= nLayers) {
      layer = nLayers - 1;
    }
    if (cubic) {
      ++out.cubicPerLayer[static_cast<std::size_t>(layer)];
      ++nC;
    } else {
      ++out.hexPerLayer[static_cast<std::size_t>(layer)];
      ++nH;
    }
  }
  out.phiC = (nC + nH) > 0
                 ? static_cast<double>(nC) / static_cast<double>(nC + nH)
                 : 0.0;
  out.sequence.resize(static_cast<std::size_t>(nLayers), '.');
  for (int k = 0; k < nLayers; k++) {
    const int c = out.cubicPerLayer[static_cast<std::size_t>(k)];
    const int h = out.hexPerLayer[static_cast<std::size_t>(k)];
    if (c == 0 && h == 0) {
      out.sequence[static_cast<std::size_t>(k)] = '.';
    } else if (c > h) {
      out.sequence[static_cast<std::size_t>(k)] = 'C';
    } else if (h > c) {
      out.sequence[static_cast<std::size_t>(k)] = 'H';
    } else {
      out.sequence[static_cast<std::size_t>(k)] = 'M';
    }
  }
  return out;
}

topoparam::LayerStack
topoparam::tumLayerStack(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &rings, const std::vector<bool> &basal,
    const std::vector<bool> &equatorial, int axis, double layerWidth) {
  LayerStack out;
  if (yCloud.nop <= 0 || axis < 0 || axis > 2 || yCloud.box.size() < 3) {
    return out;
  }
  const double L = yCloud.box[static_cast<std::size_t>(axis)];
  if (!(L > 0.0) || !(layerWidth > 0.0)) {
    return out;
  }
  const int nLayers = std::max(1, static_cast<int>(std::lround(L / layerWidth)));
  const double w = L / static_cast<double>(nLayers);
  out.cubicPerLayer.assign(static_cast<std::size_t>(nLayers), 0);
  out.hexPerLayer.assign(static_cast<std::size_t>(nLayers), 0);
  const double lo = yCloud.boxLow.size() > static_cast<std::size_t>(axis)
                        ? yCloud.boxLow[static_cast<std::size_t>(axis)]
                        : 0.0;
  int nBasal = 0;
  int nEq = 0;
  const std::size_t nR = rings.size();
  for (std::size_t r = 0; r < nR; r++) {
    const bool h = r < basal.size() && basal[r];
    const bool c = r < equatorial.size() && equatorial[r];
    if (!h && !c) {
      continue;
    }
    if (rings[r].empty()) {
      continue;
    }
    const int a0 = rings[r][0];
    if (a0 < 0 || a0 >= yCloud.nop) {
      continue;
    }
    double acc = axis == 0 ? yCloud.pts[static_cast<std::size_t>(a0)].x
                           : (axis == 1 ? yCloud.pts[static_cast<std::size_t>(a0)].y
                                        : yCloud.pts[static_cast<std::size_t>(a0)].z);
    for (std::size_t k = 1; k < rings[r].size(); k++) {
      const int a = rings[r][k];
      if (a < 0 || a >= yCloud.nop) {
        continue;
      }
      const auto dr = gen::relDist(yCloud, a, a0);
      acc += (axis == 0 ? yCloud.pts[static_cast<std::size_t>(a0)].x + dr[0]
                        : (axis == 1 ? yCloud.pts[static_cast<std::size_t>(a0)].y + dr[1]
                                     : yCloud.pts[static_cast<std::size_t>(a0)].z + dr[2]));
    }
    acc /= static_cast<double>(rings[r].size());
    double u = acc - lo;
    u -= L * std::floor(u / L);
    if (u < 0.0) {
      u += L;
    }
    if (u >= L) {
      u = 0.0;
    }
    int layer = static_cast<int>(std::floor(u / w));
    if (layer < 0) {
      layer = 0;
    }
    if (layer >= nLayers) {
      layer = nLayers - 1;
    }
    if (c) {
      ++out.cubicPerLayer[static_cast<std::size_t>(layer)];
      ++nEq;
    } else {
      ++out.hexPerLayer[static_cast<std::size_t>(layer)];
      ++nBasal;
    }
  }
  out.phiC = (nBasal + nEq) > 0
                 ? static_cast<double>(nEq) / static_cast<double>(nBasal + nEq)
                 : 0.0;
  out.sequence.resize(static_cast<std::size_t>(nLayers), '.');
  for (int k = 0; k < nLayers; k++) {
    const int c = out.cubicPerLayer[static_cast<std::size_t>(k)];
    const int h = out.hexPerLayer[static_cast<std::size_t>(k)];
    if (c == 0 && h == 0) {
      out.sequence[static_cast<std::size_t>(k)] = '.';
    } else if (c > h) {
      out.sequence[static_cast<std::size_t>(k)] = 'C';
    } else if (h > c) {
      out.sequence[static_cast<std::size_t>(k)] = 'H';
    } else {
      out.sequence[static_cast<std::size_t>(k)] = 'M';
    }
  }
  return out;
}
