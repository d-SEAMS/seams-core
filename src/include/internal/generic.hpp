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
#include <filesystem>
#include <iostream>
#include <mol_sys.hpp>

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

// Recover H (columns a, b, c) and origin from a LAMMPS dump box with the
// same mapping as nneigh::lammpsBoxToLcCell: box[0..3] are bound spans,
// box[3..6] are tilt factors xy, xz, yz, boxLow is the bound lo.
// Convert each cartesian point to fractional coordinates relative to
// origin, wrap the fractional difference, and map back through H.
inline std::array<double, 3> triclinicMinImage(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double xi,
    double yi, double zi, double xj, double yj, double zj) {
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

  auto toFrac = [&](double x, double y, double z) {
    const double sz = (z - oz) / lz;
    const double sy = (y - oy - yz * sz) / ly;
    const double sx = (x - ox - xy * sy - xz * sz) / lx;
    return std::array<double, 3>{sx, sy, sz};
  };
  const auto si = toFrac(xi, yi, zi);
  const auto sj = toFrac(xj, yj, zj);
  double dsx = si[0] - sj[0];
  double dsy = si[1] - sj[1];
  double dsz = si[2] - sj[2];
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
periodicDistSq(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               int iatom, int jatom) {
  if (yCloud.box.size() >= 6) {
    const auto dr = triclinicMinImage(
        yCloud, yCloud.pts[iatom].x, yCloud.pts[iatom].y, yCloud.pts[iatom].z,
        yCloud.pts[jatom].x, yCloud.pts[jatom].y, yCloud.pts[jatom].z);
    return dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
  }

  std::array<double, 3> dr;
  double r2 = 0.0; // Squared absolute distance

  // Get x1-x2 etc
  dr[0] = std::fabs(yCloud.pts[iatom].x - yCloud.pts[jatom].x);
  dr[1] = std::fabs(yCloud.pts[iatom].y - yCloud.pts[jatom].y);
  dr[2] = std::fabs(yCloud.pts[iatom].z - yCloud.pts[jatom].z);

  // Three-axis wrap for an orthorhombic length-3 box
  for (int k = 0; k < 3; k++) {
    dr[k] -= yCloud.box[k] * std::round(dr[k] / yCloud.box[k]);
    r2 += dr[k] * dr[k];
  }

  return r2;
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
inline double unWrappedDistFromPoint(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int iatom,
    std::vector<double> singlePoint) {
  if (yCloud.box.size() >= 6) {
    const auto dr = triclinicMinImage(yCloud, yCloud.pts[iatom].x,
                                      yCloud.pts[iatom].y, yCloud.pts[iatom].z,
                                      singlePoint[0], singlePoint[1],
                                      singlePoint[2]);
    return std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
  }

  std::array<double, 3> dr;
  double r2 = 0.0; // Squared absolute distance

  // Get x1-x2 etc
  dr[0] = fabs(yCloud.pts[iatom].x - singlePoint[0]);
  dr[1] = fabs(yCloud.pts[iatom].y - singlePoint[1]);
  dr[2] = fabs(yCloud.pts[iatom].z - singlePoint[2]);

  // Three-axis wrap for an orthorhombic length-3 box
  for (int k = 0; k < 3; k++) {
    dr[k] -= yCloud.box[k] * round(dr[k] / yCloud.box[k]);
    r2 += pow(dr[k], 2.0);
  }

  return sqrt(r2);
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
  if (yCloud.box.size() >= 6) {
    return triclinicMinImage(
        yCloud, yCloud.pts[iatom].x, yCloud.pts[iatom].y, yCloud.pts[iatom].z,
        yCloud.pts[jatom].x, yCloud.pts[jatom].y, yCloud.pts[jatom].z);
  }

  std::array<double, 3> dr;
  std::array<double, 3> box = {yCloud.box[0], yCloud.box[1], yCloud.box[2]};

  // Get x1-x2 etc
  dr[0] = yCloud.pts[iatom].x - yCloud.pts[jatom].x;
  dr[1] = yCloud.pts[iatom].y - yCloud.pts[jatom].y;
  dr[2] = yCloud.pts[iatom].z - yCloud.pts[jatom].z;

  // Three-axis wrap for an orthorhombic length-3 box
  for (int k = 0; k < 3; k++) {
    if (dr[k] < -box[k] * 0.5) {
      dr[k] = dr[k] + box[k];
    }
    if (dr[k] >= box[k] * 0.5) {
      dr[k] = dr[k] - box[k];
    }
  }

  return dr;
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
