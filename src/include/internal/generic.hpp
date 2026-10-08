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

#ifndef SEAMS_GENERIC_H_
#define SEAMS_GENERIC_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <mol_sys.hpp>

#ifdef SEAMS_HAS_MINIMAGE
#include <minimage.h>
#endif

// C++20
#include <numbers>
// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>

/** @file generic.hpp
 *   @brief File for containing generic or common functions.
 */

/**
 *  @addtogroup gen
 *  @{
 */

/** \brief Small generic functions that are shared by all namespaces.
 *
 * These are general functions (eg. for finding the periodic distance) which are
 * required by several namespaces.
 *
 *  ### Changelog ###
 *
 *  - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 14, 2019
 *  - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

namespace gen {

/**
 *  Uses Boost to get the value of pi.
 */
constexpr double pi = std::numbers::pi;

/**
 *  Inline function for converting radians->degrees.
 *  @param[in] angle The input angle, in radians
 *  @return The input angle, in degrees
 */
inline double radDeg(double angle) { return (angle * 180) / gen::pi; }

//! Eigen function for getting the angle (in radians) between the O--O and O-H vectors
double eigenVecAngle(std::vector<double> OO, std::vector<double> OH);

//! Get the average, after excluding the outliers, using quartiles
double getAverageWithoutOutliers(std::vector<double> inpVec);

/**
 *  @brief Inline generic function for calculating the median given a vector of
 * the values
 *  @param[in] yCloud The input PointCloud, which contains the particle
 * coordinates, simulation box lengths etc.
 *  @param[in] input The input vector with the values
 *  @return The median value
 */
inline double calcMedian(std::vector<double> *input) {
  int n = (*input).size(); // Number of elements
  double median;           // Output median value

  // Sort a copy (avoid mutating input)
  std::vector<double> sorted = *input;
  std::sort(sorted.begin(), sorted.end());

  // Calculate the median
  // For even values, the median is the average of the two middle values
  if (n % 2 == 0) {
    median = 0.5 * (sorted[n / 2] + sorted[n / 2 - 1]);
  } else {
    median = sorted[(n + 1) / 2 - 1];
  }

  return median;
}

#ifdef SEAMS_HAS_MINIMAGE
// Dump bound spans plus optional tilt onto an mi_cell.
inline bool pointCloudCell(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    mi_cell *out) {
  const auto &box = yCloud.box;
  const auto &boxLow = yCloud.boxLow;
  const double xlo_b = boxLow.size() > 0 ? boxLow[0] : 0.0;
  const double ylo_b = boxLow.size() > 1 ? boxLow[1] : 0.0;
  const double zlo_b = boxLow.size() > 2 ? boxLow[2] : 0.0;
  if (box.size() >= 6) {
    return mi_cell_from_lammps_bounds(box[0], box[1], box[2], box[3], box[4],
                                      box[5], xlo_b, ylo_b, zlo_b, out) == 0;
  }
  if (box.size() >= 3 && box[0] > 0.0 && box[1] > 0.0 && box[2] > 0.0) {
    *out = mi_cell_ortho(box[0], box[1], box[2]);
    out->ox = xlo_b;
    out->oy = ylo_b;
    out->oz = zlo_b;
    return true;
  }
  return false;
}
#endif

// One dump box, recovered once per frame. Orthorhombic lengths are the
// bound spans. A tilt dump stores those spans in box[0..2] and xy, xz, yz
// in box[3..5]; lx, ly, lz are the recovered edge lengths, and wx, wy, wz
// the separations of opposite faces, which tilt makes shorter.
struct FracBox {
  double lx = 0.0;
  double ly = 0.0;
  double lz = 0.0;
  double xy = 0.0;
  double xz = 0.0;
  double yz = 0.0;
  double ox = 0.0;
  double oy = 0.0;
  double oz = 0.0;
  double wx = 0.0;
  double wy = 0.0;
  double wz = 0.0;
  double halfMin = 0.0;
  bool triclinic = false;
  bool ok = false;
};

inline FracBox
makeFracBox(const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  FracBox b;
  const auto &box = yCloud.box;
  const auto &boxLow = yCloud.boxLow;
  b.ox = boxLow.size() > 0 ? boxLow[0] : 0.0;
  b.oy = boxLow.size() > 1 ? boxLow[1] : 0.0;
  b.oz = boxLow.size() > 2 ? boxLow[2] : 0.0;
  if (box.size() >= 6) {
    const double xy = box[3];
    const double xz = box[4];
    const double yz = box[5];
    const double xmin = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
    const double xmax = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
    const double ymin = std::min(0.0, yz);
    const double ymax = std::max(0.0, yz);
    b.lx = box[0] - xmax + xmin;
    b.ly = box[1] - ymax + ymin;
    b.lz = box[2];
    b.xy = xy;
    b.xz = xz;
    b.yz = yz;
    b.ox -= xmin;
    b.oy -= ymin;
    b.triclinic = true;
  } else if (box.size() >= 3) {
    b.lx = box[0];
    b.ly = box[1];
    b.lz = box[2];
  } else {
    return b;
  }
  if (!(b.lx > 0.0 && b.ly > 0.0 && b.lz > 0.0)) {
    return b;
  }
  const double bcz = b.xy * b.yz - b.ly * b.xz;
  b.wx = b.lx * b.ly * b.lz /
         std::sqrt(b.ly * b.lz * b.ly * b.lz + b.xy * b.lz * b.xy * b.lz +
                   bcz * bcz);
  b.wy = b.ly * b.lz / std::sqrt(b.lz * b.lz + b.yz * b.yz);
  b.wz = b.lz;
  b.halfMin = 0.5 * std::min(b.wx, std::min(b.wy, b.wz));
  b.ok = true;
  return b;
}

// Displacement from (xj, yj, zj) to (xi, yi, zi). Orthorhombic axes wrap
// independently. A tilt box wraps the fractional difference and maps it
// back through H. That vector is the Euclidean minimum image while its
// length stays below half the shortest edge.
inline std::array<double, 3> fracDelta(const FracBox &b, double xi, double yi,
                                       double zi, double xj, double yj,
                                       double zj) {
  if (!b.triclinic) {
    std::array<double, 3> dr = {xi - xj, yi - yj, zi - zj};
    const double len[3] = {b.lx, b.ly, b.lz};
    for (int k = 0; k < 3; k++) {
      dr[static_cast<std::size_t>(k)] -=
          len[k] * std::floor(dr[static_cast<std::size_t>(k)] / len[k] + 0.5);
    }
    return dr;
  }
  const double sz_i = (zi - b.oz) / b.lz;
  const double sy_i = (yi - b.oy - b.yz * sz_i) / b.ly;
  const double sx_i = (xi - b.ox - b.xy * sy_i - b.xz * sz_i) / b.lx;
  const double sz_j = (zj - b.oz) / b.lz;
  const double sy_j = (yj - b.oy - b.yz * sz_j) / b.ly;
  const double sx_j = (xj - b.ox - b.xy * sy_j - b.xz * sz_j) / b.lx;
  double dsx = sx_i - sx_j;
  double dsy = sy_i - sy_j;
  double dsz = sz_i - sz_j;
  dsx -= std::round(dsx);
  dsy -= std::round(dsy);
  dsz -= std::round(dsz);
  return {b.lx * dsx + b.xy * dsy + b.xz * dsz, b.ly * dsy + b.yz * dsz,
          b.lz * dsz};
}

inline double fracDistSq(const FracBox &b, double xi, double yi, double zi,
                         double xj, double yj, double zj) {
  const auto dr = fracDelta(b, xi, yi, zi, xj, yj, zj);
  return dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
}

// Smith, CCP5 1989: below half the narrowest face separation the fractional
// wrap is the only lattice image in the ball, and every image inside that
// ball is a fractional wrap. An orthorhombic wrap is that image on every
// axis, including past half an edge.
inline bool smithInside(const FracBox &b, double r2) {
  return b.ok && (!b.triclinic || std::sqrt(r2) + 1e-12 < b.halfMin);
}

#ifdef SEAMS_HAS_MINIMAGE
inline bool euclideanDelta(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double xi,
    double yi, double zi, double xj, double yj, double zj, double dr[3]) {
  mi_cell cell;
  if (!pointCloudCell(yCloud, &cell)) {
    return false;
  }
  const double p[3] = {xj, yj, zj};
  const double q[3] = {xi, yi, zi};
  return mi_displacement_euclidean(&cell, p, q, dr) == 0;
}

// One frame on minimage's cell, built once, with the positions packed for
// its batch kernels.
struct CellFrame {
  mi_cell cell{};
  std::vector<double> xyz;
  bool ok = false;
};

inline CellFrame
makeCellFrame(const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  CellFrame f;
  if (!pointCloudCell(yCloud, &f.cell)) {
    return f;
  }
  const std::size_t n = yCloud.pts.size();
  f.xyz.resize(3 * n);
  for (std::size_t i = 0; i < n; i++) {
    f.xyz[3 * i] = yCloud.pts[i].x;
    f.xyz[3 * i + 1] = yCloud.pts[i].y;
    f.xyz[3 * i + 2] = yCloud.pts[i].z;
  }
  f.ok = true;
  return f;
}

// Squared Euclidean minimum image from j to i on the frame's cell, or -1.
inline double euclideanDistSq(const CellFrame &f, int i, int j) {
  double dr[3];
  if (mi_displacement_euclidean(&f.cell, &f.xyz[3 * static_cast<std::size_t>(j)],
                                &f.xyz[3 * static_cast<std::size_t>(i)],
                                dr) != 0) {
    return -1.0;
  }
  return dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
}
#endif

// Recover H (columns a, b, c) and origin from a LAMMPS dump box with the
// same mapping as nneigh::lammpsBoxToLcCell: box[0..3] are bound spans,
// box[3..6] are tilt factors xy, xz, yz, boxLow is the bound lo.
inline std::array<double, 3> triclinicMinImage(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double xi,
    double yi, double zi, double xj, double yj, double zj) {
  const FracBox b = makeFracBox(yCloud);
  if (b.ok && b.triclinic) {
    const auto dr = fracDelta(b, xi, yi, zi, xj, yj, zj);
    const double r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    if (smithInside(b, r2)) {
      return dr;
    }
#ifdef SEAMS_HAS_MINIMAGE
    double ed[3] = {0.0, 0.0, 0.0};
    if (euclideanDelta(yCloud, xi, yi, zi, xj, yj, zj, ed)) {
      return {ed[0], ed[1], ed[2]};
    }
#endif
    return dr;
  }
  const auto &box = yCloud.box;
  const auto &boxLow = yCloud.boxLow;
  const double xspan = box[0];
  const double yspan = box[1];
  const double zspan = box[2];
  const double xlo_b = boxLow.size() > 0 ? boxLow[0] : 0.0;
  const double ylo_b = boxLow.size() > 1 ? boxLow[1] : 0.0;
  const double zlo_b = boxLow.size() > 2 ? boxLow[2] : 0.0;
  const double xy = box[3];
  const double xz = box[4];
  const double yz = box[5];
  const double xmin = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
  const double xmax = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
  const double ymin = std::min(0.0, yz);
  const double ymax = std::max(0.0, yz);
  const double lx = xspan - xmax + xmin;
  const double ly = yspan - ymax + ymin;
  const double lz = zspan;
  const double ox = xlo_b - xmin;
  const double oy = ylo_b - ymin;
  const double oz = zlo_b;
  const double sz_i = (zi - oz) / lz;
  const double sy_i = (yi - oy - yz * sz_i) / ly;
  const double sx_i = (xi - ox - xy * sy_i - xz * sz_i) / lx;
  const double sz_j = (zj - oz) / lz;
  const double sy_j = (yj - oy - yz * sz_j) / ly;
  const double sx_j = (xj - ox - xy * sy_j - xz * sz_j) / lx;
  double dsx = sx_i - sx_j;
  double dsy = sy_i - sy_j;
  double dsz = sz_i - sz_j;
  dsx -= std::round(dsx);
  dsy -= std::round(dsy);
  dsz -= std::round(dsz);
  return {lx * dsx + xy * dsy + xz * dsz, ly * dsy + yz * dsz, lz * dsz};
}

// Generic function for getting the unwrapped distance
/**
 *  @brief Inline generic function for obtaining the unwrapped periodic distance
 *  between two particles, whose indices (not IDs) have been given.
 *  @param[in] yCloud The input PointCloud, which contains the particle
 * coordinates, simulation box lengths etc.
 *  @param[in] iatom The index of the @f$ i^{th} @f$ atom.
 *  @param[in] jatom The index of the @f$ j^{th} @f$ atom.
 *  @return The unwrapped periodic distance.
 */
inline double
periodicDistSq(const FracBox &b,
               const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               int iatom, int jatom) {
  const auto &pi = yCloud.pts[static_cast<std::size_t>(iatom)];
  const auto &pj = yCloud.pts[static_cast<std::size_t>(jatom)];
  if (b.ok) {
    const double r2 = fracDistSq(b, pi.x, pi.y, pi.z, pj.x, pj.y, pj.z);
    if (smithInside(b, r2)) {
      return r2;
    }
#ifdef SEAMS_HAS_MINIMAGE
    double ed[3] = {0.0, 0.0, 0.0};
    if (euclideanDelta(yCloud, pi.x, pi.y, pi.z, pj.x, pj.y, pj.z, ed)) {
      return ed[0] * ed[0] + ed[1] * ed[1] + ed[2] * ed[2];
    }
#endif
    return r2;
  }
  if (yCloud.box.size() >= 6) {
    const auto dr =
        triclinicMinImage(yCloud, pi.x, pi.y, pi.z, pj.x, pj.y, pj.z);
    return dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
  }
  std::array<double, 3> dr;
  double r2 = 0.0;
  dr[0] = std::fabs(pi.x - pj.x);
  dr[1] = std::fabs(pi.y - pj.y);
  dr[2] = std::fabs(pi.z - pj.z);
  for (int k = 0; k < 3; k++) {
    dr[static_cast<std::size_t>(k)] -=
        yCloud.box[static_cast<std::size_t>(k)] *
        std::round(dr[static_cast<std::size_t>(k)] /
                   yCloud.box[static_cast<std::size_t>(k)]);
    r2 += dr[static_cast<std::size_t>(k)] * dr[static_cast<std::size_t>(k)];
  }
  return r2;
}

inline double
periodicDistSq(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               int iatom, int jatom) {
  return periodicDistSq(makeFracBox(yCloud), yCloud, iatom, jatom);
}

/**
 *  @brief Inline generic function for obtaining the unwrapped periodic distance
 *  between two particles, whose indices (not IDs) have been given.
 *  @param[in] yCloud The input PointCloud, which contains the particle
 * coordinates, simulation box lengths etc.
 *  @param[in] iatom The index of the @f$ i^{th} @f$ atom.
 *  @param[in] jatom The index of the @f$ j^{th} @f$ atom.
 *  @return The unwrapped periodic distance.
 */
inline double
periodicDist(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
             int iatom, int jatom) {
  return std::sqrt(periodicDistSq(yCloud, iatom, jatom));
}

// Tilt batches share one mi_cell. mi_dist2_many is the fractional
// image for every candidate. A result outside the Smith ball is
// replaced by the Euclidean image. Orthorhombic batches stay on
// periodicDistSq; the Highway kernel is the ortho difference path.
inline void batchPeriodicDistSq(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int iatom,
    const int *jatom, std::size_t n, double *distSq) {
  const FracBox b = makeFracBox(yCloud);
#ifdef SEAMS_HAS_MINIMAGE
  if (n > 0 && b.ok && b.triclinic) {
    mi_cell cell;
    if (pointCloudCell(yCloud, &cell)) {
      const auto &pi = yCloud.pts[static_cast<std::size_t>(iatom)];
      const double p[3] = {pi.x, pi.y, pi.z};
      std::vector<double> qs(n * 3);
      for (std::size_t k = 0; k < n; ++k) {
        const auto &pj = yCloud.pts[static_cast<std::size_t>(jatom[k])];
        qs[3 * k] = pj.x;
        qs[3 * k + 1] = pj.y;
        qs[3 * k + 2] = pj.z;
      }
      if (mi_dist2_many(&cell, p, qs.data(), n, distSq) == 0) {
        for (std::size_t k = 0; k < n; ++k) {
          if (smithInside(b, distSq[k])) {
            continue;
          }
          const double q[3] = {qs[3 * k], qs[3 * k + 1], qs[3 * k + 2]};
          double ed[3] = {0.0, 0.0, 0.0};
          if (mi_displacement_euclidean(&cell, q, p, ed) == 0) {
            distSq[k] = ed[0] * ed[0] + ed[1] * ed[1] + ed[2] * ed[2];
          }
        }
        return;
      }
    }
  }
#endif
  for (std::size_t k = 0; k < n; k++) {
    distSq[k] = periodicDistSq(b, yCloud, iatom, jatom[k]);
  }
}

// Bound spans, then tilt when box.size() >= 6.
inline std::string formatDumpBox(const std::vector<double> &box) {
  std::ostringstream oss;
  if (box.size() >= 3) {
    oss << box[0] << ' ' << box[1] << ' ' << box[2];
  }
  if (box.size() >= 6) {
    oss << " xy " << box[3] << " xz " << box[4] << " yz " << box[5];
  }
  return oss.str();
}

// ITEM line plus three bound lines. Tilt is a third field per line.
inline void writeDumpBoxBounds(
    std::ostream &os,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  const bool tilt = yCloud.box.size() >= 6;
  if (tilt) {
    os << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n";
  } else {
    os << "ITEM: BOX BOUNDS pp pp pp\n";
  }
  for (int k = 0; k < 3; k++) {
    const double lo = (static_cast<std::size_t>(k) < yCloud.boxLow.size())
                          ? yCloud.boxLow[k]
                          : 0.0;
    const double len = (static_cast<std::size_t>(k) < yCloud.box.size())
                           ? yCloud.box[k]
                           : 0.0;
    os << lo << ' ' << lo + len;
    if (tilt) {
      os << ' ' << yCloud.box[static_cast<std::size_t>(k + 3)];
    }
    os << '\n';
  }
}

/**
 *  Inline generic function for obtaining
 *  the unwrapped periodic distance between one particle and another point,
 *  whose index has been given.
 *  @param[in] yCloud The input PointCloud, which contains the particle
 *  coordinates, simulation box lengths etc.
 *  @param[in] iatom The index of the \f$ i^{th} \f$ atom.
 *  @param[in] singlePoint Vector containing coordinate values
 *  \return The unwrapped periodic distance.
 */
inline std::array<double, 3> relDistFromPoint(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int iatom,
    double xj, double yj, double zj) {
  const FracBox b = makeFracBox(yCloud);
  if (b.ok) {
    const auto &pi = yCloud.pts[static_cast<std::size_t>(iatom)];
    const auto dr = fracDelta(b, pi.x, pi.y, pi.z, xj, yj, zj);
    const double r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    if (smithInside(b, r2)) {
      return dr;
    }
#ifdef SEAMS_HAS_MINIMAGE
    double ed[3] = {0.0, 0.0, 0.0};
    if (euclideanDelta(yCloud, pi.x, pi.y, pi.z, xj, yj, zj, ed)) {
      return {ed[0], ed[1], ed[2]};
    }
#endif
    return dr;
  }
  if (yCloud.box.size() >= 6) {
    return triclinicMinImage(yCloud, yCloud.pts[iatom].x, yCloud.pts[iatom].y,
                             yCloud.pts[iatom].z, xj, yj, zj);
  }

  std::array<double, 3> dr = {yCloud.pts[iatom].x - xj, yCloud.pts[iatom].y - yj,
                              yCloud.pts[iatom].z - zj};
  for (int k = 0; k < 3; k++) {
    if (dr[k] < -yCloud.box[k] * 0.5) {
      dr[k] += yCloud.box[k];
    }
    if (dr[k] >= yCloud.box[k] * 0.5) {
      dr[k] -= yCloud.box[k];
    }
  }
  return dr;
}

inline double unWrappedDistFromPoint(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int iatom,
    const std::vector<double> &singlePoint) {
  const auto dr =
      relDistFromPoint(yCloud, iatom, singlePoint[0], singlePoint[1],
                       singlePoint[2]);
  return std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
}

// Generic function for getting the distance (no PBCs applied)
/**
 * @brief Inline generic function for obtaining the wrapped distance between two
 * particles WITHOUT applying PBCs, whose indices (not IDs) have been given.
 *  @param[in] yCloud The input PointCloud, which contains the particle
 coordinates, simulation box lengths etc.
 *  @param[in] iatom The index of the \f$ i^{th} \f$ atom.
 *  @param[in] jatom The index of the \f$ j^{th} \f$ atom.
 *  @return The wrapped distance.
 */
inline double
distance(const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int iatom,
         int jatom) {
  std::array<double, 3> dr;
  double r2 = 0.0; // Squared absolute distance

  // Get x1-x2 etc
  dr[0] = fabs(yCloud.pts[iatom].x - yCloud.pts[jatom].x);
  dr[1] = fabs(yCloud.pts[iatom].y - yCloud.pts[jatom].y);
  dr[2] = fabs(yCloud.pts[iatom].z - yCloud.pts[jatom].z);

  // Get the squared absolute distance
  for (int k = 0; k < 3; k++) {
    r2 += pow(dr[k], 2.0);
  }

  return sqrt(r2);
}

// Generic function for getting the relative coordinates
/**
 *  Inline generic function for getting the relative unwrapped distance between
 *  two particles for each dimension. The indices (not IDs) of the particles
 * have been given.
 *  @param[in] yCloud The input PointCloud, which contains the particle
 *  coordinates, simulation box lengths etc.
 *  @param[in] iatom The index of the \f$ i^{th} \f$ atom.
 *  @param[in] jatom The index of the \f$ j^{th} \f$ atom.
 *  @return The unwrapped relative distances for each dimension.
 */
inline std::array<double, 3>
relDist(const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int iatom,
        int jatom) {
  return relDistFromPoint(yCloud, iatom, yCloud.pts[jatom].x,
                          yCloud.pts[jatom].y, yCloud.pts[jatom].z);
}

// Function for sorting according to atom ID
// Comparator for std::sort
/**
 *  Inline generic function for sorting or comparing two particles, according to
 *  the atom ID when the entire Point objects have been passed.
 *  @param[in] a The input Point for A.
 *  @param[in] b The input Point for B.
 *  @return True if the atom ID of A is less than the atom ID of B
 */
inline bool compareByAtomID(const molSys::Point<double> &a,
                            const molSys::Point<double> &b) {
  return a.atomID < b.atomID;
}

//! Generic function for printing all the struct information
[[nodiscard]] int prettyPrintYoda(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    std::string outFile);

//! Shift particles (unwrapped coordinates)
[[nodiscard]] int unwrappedCoordShift(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int iatomIndex,
    int jatomIndex, double *x_i, double *y_i, double *z_i, double *x_j,
    double *y_j, double *z_j);

//! Function for getting the angular distance between two quaternions. Returns
//! the result in degrees
double angDistDegQuaternions(std::vector<double> quat1,
                             std::vector<double> quat2);

/**
 * @brief Function for tokenizing line strings into words (strings) delimited by
 * whitespace. This returns a vector with the words in it.
 * @param[in] line The string containing the line to be tokenized
 */
inline std::vector<std::string> tokenizer(std::string line) {
  std::istringstream iss(line);
  std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                  std::istream_iterator<std::string>{}};
  return tokens;
}

/**
 *  @brief Function for tokenizing line strings into a vector of doubles.
 *  @param[in] line The string containing the line to be tokenized
 */
inline std::vector<double> tokenizerDouble(std::string line) {
  std::istringstream iss(line);
  std::vector<double> tokens;
  double number; // Each number being read in from the line
  while (iss >> number) {
    tokens.push_back(number);
  }
  return tokens;
}

/**
 * @brief Function for tokenizing line strings into a vector of ints.
 * @param[in] line The string containing the line to be tokenized
 */
inline std::vector<int> tokenizerInt(std::string line) {
  std::istringstream iss(line);
  std::vector<int> tokens;
  int number; // Each number being read in from the line
  while (iss >> number) {
    tokens.push_back(number);
  }
  return tokens;
}

/**
 *  @brief Function for checking if a file exists or not.
 *  @param[in] name The name of the file
 */
inline bool file_exists(const std::string &name) {
  return std::filesystem::exists(name);
}

/**
 *   Calculates the complex vector, normalized by the number of nearest
 * neighbours, of length @f$2l+1@f$.
 *   @param[in] v The complex vector to be normalized, of length @f$2l+1@f$
 *   @param[in] l A free integer parameter
 *   @param[in] neigh The number of nearest neighbours
 *   @return length @f$2l+1@f$, normalized by the number of nearest neighbours
 */
inline std::vector<std::complex<double>>
avgVector(std::vector<std::complex<double>> v, int l, int neigh) {
  if (neigh == 0) {
    return v;
  }
  for (int m = 0; m < 2 * l + 1; m++) {
    v[m] = (1.0 / static_cast<double>(neigh)) * v[m];
  }

  return v;
}

} // namespace gen

#endif // SEAMS_GENERIC_H_
