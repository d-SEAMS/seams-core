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

#ifndef __GENERIC_H_
#define __GENERIC_H_

#include <array>
#include <iostream>
#include <math.h>
#include <mol_sys.hpp>

// Boost
#include <boost/math/constants/constants.hpp>
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
const double pi = boost::math::constants::pi<double>();

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

  // Sort the vector
  std::sort((*input).begin(), (*input).end());

  // Calculate the median
  // For even values, the median is the average of the two middle values
  if (n % 2 == 0) {
    median = 0.5 * ((*input)[n / 2] + (*input)[n / 2 - 1]); // n/2+n/2-1
  } // median is average of middle values
  else {
    median = (*input)[(n + 1) / 2 -
                      1]; // middle value of 7 elements is the 4th element
  }                       // if odd, it is the middle value

  return median;
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
periodicDist(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
             int iatom, int jatom) {
  std::array<double, 3> dr;
  double r2 = 0.0; // Squared absolute distance

  // Get x1-x2 etc
  dr[0] = fabs(yCloud->pts[iatom].x - yCloud->pts[jatom].x);
  dr[1] = fabs(yCloud->pts[iatom].y - yCloud->pts[jatom].y);
  dr[2] = fabs(yCloud->pts[iatom].z - yCloud->pts[jatom].z);

  // Get the squared absolute distance
  for (int k = 0; k < 3; k++) {
    // Correct for periodicity
    dr[k] -= yCloud->box[k] * round(dr[k] / yCloud->box[k]);
    r2 += pow(dr[k], 2.0);
  }

  return sqrt(r2);
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
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatom,
    std::vector<double> singlePoint) {
  std::array<double, 3> dr;
  double r2 = 0.0; // Squared absolute distance

  // Get x1-x2 etc
  dr[0] = fabs(yCloud->pts[iatom].x - singlePoint[0]);
  dr[1] = fabs(yCloud->pts[iatom].y - singlePoint[1]);
  dr[2] = fabs(yCloud->pts[iatom].z - singlePoint[2]);

  // Get the squared absolute distance
  for (int k = 0; k < 3; k++) {
    // Correct for periodicity
    dr[k] -= yCloud->box[k] * round(dr[k] / yCloud->box[k]);
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
distance(molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatom,
         int jatom) {
  std::array<double, 3> dr;
  double r2 = 0.0; // Squared absolute distance

  // Get x1-x2 etc
  dr[0] = fabs(yCloud->pts[iatom].x - yCloud->pts[jatom].x);
  dr[1] = fabs(yCloud->pts[iatom].y - yCloud->pts[jatom].y);
  dr[2] = fabs(yCloud->pts[iatom].z - yCloud->pts[jatom].z);

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
relDist(molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatom,
        int jatom) {
  std::array<double, 3> dr;
  std::array<double, 3> box = {yCloud->box[0], yCloud->box[1], yCloud->box[2]};
  double r2 = 0.0; // Squared absolute distance

  // Get x1-x2 etc
  dr[0] = yCloud->pts[iatom].x - yCloud->pts[jatom].x;
  dr[1] = yCloud->pts[iatom].y - yCloud->pts[jatom].y;
  dr[2] = yCloud->pts[iatom].z - yCloud->pts[jatom].z;

  // Get the relative distance
  for (int k = 0; k < 3; k++) {
    //
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
int prettyPrintYoda(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::string outFile);

//! Shift particles (unwrapped coordinates)
int unwrappedCoordShift(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatomIndex,
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
  // Replace by boost function later
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
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
    v[m] = (1.0 / (double)neigh) * v[m];
  }

  return v;
}

} // namespace gen

#endif // __NEIGHBOURS_H_
