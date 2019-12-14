#ifndef __GENERIC_H_
#define __GENERIC_H_

#include <math.h>
#include <array>
#include <iostream>
#include <mol_sys.hpp>

// Boost
#include <gsl/gsl_blas.h>
#include <boost/math/constants/constants.hpp>

/*! \file generic.hpp
    \brief File for containing generic or common functions.

    Details.
*/

/*!
 *  \addtogroup gen
 *  @{
 */

/*! \brief Small generic functions that are shared by all namespaces.
 *
These are general functions (eg. for finding the periodic distance) which are
required by several namespaces.

  ### Changelog ###

  - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 14, 2019
 */

namespace gen {

/********************************************/ /**
 *  Uses Boost to
 get the value of pi.
 ***********************************************/
const double pi = boost::math::constants::pi<double>();

/********************************************/ /**
 *  Inline function for converting
 radians->degrees.
 *  @param[in] angle The input angle, in radians
 *  \return The input angle, in degrees
 ***********************************************/
inline double radDeg(double angle) { return (angle * 180) / gen::pi; }

// GSL for getting the angle (in radians) between the O--O and O-H vectors
double gslVecAngle(std::vector<double> OO, std::vector<double> OH);

// Generic function for getting the unwrapped distance
/********************************************/ /**
 *  Inline generic function for obtaining
 the unwrapped periodic distance between two particles, whose indices (not IDs)
 have been given.
 *  @param[in] yCloud The input PointCloud, which contains the particle
 coordinates, simulation box lengths etc.
 *  @param[in] iatom The index of the \f$ i^{th} \f$ atom.
 *  @param[in] jatom The index of the \f$ j^{th} \f$ atom.
 *  \return The unwrapped periodic distance.
 ***********************************************/
inline double periodicDist(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatom,
    int jatom) {
  std::array<double, 3> dr;
  double r2 = 0.0;  // Squared absolute distance

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

// Generic function for getting the distance (no PBCs applied)
/********************************************/ /**
 *  Inline generic function for obtaining
 the wrapped distance between two particles WITHOUT applying PBCs, whose indices
 (not IDs) have been given.
 *  @param[in] yCloud The input PointCloud, which contains the particle
 coordinates, simulation box lengths etc.
 *  @param[in] iatom The index of the \f$ i^{th} \f$ atom.
 *  @param[in] jatom The index of the \f$ j^{th} \f$ atom.
 *  \return The wrapped distance.
 ***********************************************/
inline double distance(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatom,
    int jatom) {
  std::array<double, 3> dr;
  double r2 = 0.0;  // Squared absolute distance

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
/********************************************/ /**
 *  Inline generic function for getting
 the relative unwrapped distance between two particles for each dimension. The
 indices (not IDs) of the particles have been given.
 *  @param[in] yCloud The input PointCloud, which contains the particle
 coordinates, simulation box lengths etc.
 *  @param[in] iatom The index of the \f$ i^{th} \f$ atom.
 *  @param[in] jatom The index of the \f$ j^{th} \f$ atom.
 *  \return The unwrapped relative distances for each dimension.
 ***********************************************/
inline std::array<double, 3> relDist(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatom,
    int jatom) {
  std::array<double, 3> dr;
  std::array<double, 3> box = {yCloud->box[0], yCloud->box[1], yCloud->box[2]};
  double r2 = 0.0;  // Squared absolute distance

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
/********************************************/ /**
 *  Inline generic function for sorting or comparing two particles, according to
 the atom ID when the entire Point objects have been passed.
 *  @param[in] a The input Point for A.
 *  @param[in] b The input Point for B.
 *  \return A bool which is true if the atom ID of A is less than the atom ID of
 B.
 ***********************************************/
inline bool compareByAtomID(const molSys::Point<double> &a,
                            const molSys::Point<double> &b) {
  return a.atomID < b.atomID;
}

// Generic function for printing all the struct information
int prettyPrintYoda(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::string outFile);

// Shift particles (unwrapped coordinates)
int unwrappedCoordShift(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatomIndex,
    int jatomIndex, double *x_i, double *y_i, double *z_i, double *x_j,
    double *y_j, double *z_j);

/********************************************/ /**
                                                *  Function for tokenizing line
                                                *strings into words (strings)
                                                *delimited by whitespace. This
                                                *returns a vector with the words
                                                *in it.
                                                *  @param[in] line The string
                                                *containing the line to be
                                                *tokenized
                                                ***********************************************/
inline std::vector<std::string> tokenizer(std::string line) {
  std::istringstream iss(line);
  std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                  std::istream_iterator<std::string>{}};
  return tokens;
}

/********************************************/ /**
                                                *  Function for tokenizing line
                                                *strings into a vector of
                                                *doubles.
                                                *  @param[in] line The string
                                                *containing the line to be
                                                *tokenized
                                                ***********************************************/
inline std::vector<double> tokenizerDouble(std::string line) {
  std::istringstream iss(line);
  std::vector<double> tokens;
  double number;  // Each number being read in from the line
  while (iss >> number) {
    tokens.push_back(number);
  }
  return tokens;
}

/********************************************/ /**
                                                *  Function for tokenizing line
                                                *strings into a vector of
                                                *ints.
                                                *  @param[in] line The string
                                                *containing the line to be
                                                *tokenized
                                                ***********************************************/
inline std::vector<int> tokenizerInt(std::string line) {
  std::istringstream iss(line);
  std::vector<int> tokens;
  int number;  // Each number being read in from the line
  while (iss >> number) {
    tokens.push_back(number);
  }
  return tokens;
}

/********************************************/ /**
                                                *  Function for checking if a
                                                *file exists or not.
                                                *  @param[in] name The name of
                                                *the file
                                                ***********************************************/
inline bool file_exists(const std::string &name) {
  // Replace by boost function later
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
}

/********************************************/ /**
                                                *  Calculates the complex
                                                *vector, normalized by the
                                                *number of nearest neighbours,
                                                *of length \f$2l+1\f$.
                                                *
                                                *  @param[in] v The complex
                                                *vector to be normalized, of
                                                *length \f$2l+1\f$
                                                *  @param[in] l A free integer
                                                *parameter
                                                *  @param[in] neigh The number
                                                *of nearest neighbours \return a
                                                *complex vector, of length
                                                *\f$2l+1\f$, normalized by the
                                                *number of nearest neighbours
                                                ***********************************************/
inline std::vector<std::complex<double>> avgVector(
    std::vector<std::complex<double>> v, int l, int neigh) {
  if (neigh == 0) {
    return v;
  }
  for (int m = 0; m < 2 * l + 1; m++) {
    v[m] = (1.0 / (double)neigh) * v[m];
  }

  return v;
}

}  // namespace gen

#endif  // __NEIGHBOURS_H_
