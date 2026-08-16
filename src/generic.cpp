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

#include <generic.hpp>
#include <algorithm>
#include <iostream>

/**
 * @details Function for printing out
 *  info in a PointCloud object.
 * @param[in] yCloud The input PointCloud to be printed.
 * @param[in] outFile The name of the output file to which the information will
 *  be printed.
 */
int gen::prettyPrintYoda(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    std::string outFile) {
  std::ofstream outputFile;
  // Create a new file in the output directory
  outputFile.open(outFile);

  if (outputFile.is_open()) {
    // First line
    outputFile << "# Frame\tAtomID\tx\ty\tz\tcij\ticeType\n";
    // Write out all the information out line by line
    for (int i = 0; i < yCloud.nop; i++) {
      outputFile << yCloud.currentFrame << "\t" << yCloud.pts[i].atomID
                 << "\t" << yCloud.pts[i].x << "\t" << yCloud.pts[i].y << "\t"
                 << yCloud.pts[i].z << "\t";
      // Print out cij
      // for(int c=0; c<yCloud.pts[i].c_ij.size(); c++){outputFile <<
      // yCloud.pts[i].c_ij[c]<<"\t";} Print out the classifier
      // TODO: Should print string representation
      outputFile << static_cast<int>(yCloud.pts[i].iceType) << "\n";
    }
  }
  // Close the file
  outputFile.close();
  return 0;
}

/**
 *  @details Function for getting the unwrapped coordinates
 *   of a pair of atoms.
 *  @param[in] yCloud The input PointCloud to be printed.
 *  @param[in] iatomIndex Index of the \f$ i^{th} \f$ atom.
 *  @param[in] jatomIndex Index of the \f$ j^{th} \f$ atom.
 *  @param[in, out] x_i X Coordinate of the \f$ i^{th} \f$ atom corresponding to
 *   the unwrapped distance.
 *  @param[in, out] y_i Y Coordinate of the \f$ i^{th} \f$ atom corresponding to
 *   the unwrapped distance.
 *  @param[in, out] z_i Z Coordinate of the \f$ i^{th} \f$ atom corresponding to
 *   the unwrapped distance.
 *  @param[in, out] x_j X Coordinate of the \f$ j^{th} \f$ atom corresponding to
 *   the unwrapped distance.
 *  @param[in, out] y_j Y Coordinate of the \f$ j^{th} \f$ atom corresponding to
 *   the unwrapped distance.
 *  @param[in, out] z_j Z Coordinate of the \f$ j^{th} \f$ atom corresponding to
 *   the unwrapped distance.
 */
int gen::unwrappedCoordShift(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int iatomIndex,
    int jatomIndex, double *x_i, double *y_i, double *z_i, double *x_j,
    double *y_j, double *z_j) {
  //
  const double x_iatom = yCloud.pts[iatomIndex].x;
  const double y_iatom = yCloud.pts[iatomIndex].y;
  const double z_iatom = yCloud.pts[iatomIndex].z;
  const auto dr = gen::relDist(yCloud, iatomIndex, jatomIndex);
  *x_i = x_iatom;
  *y_i = y_iatom;
  *z_i = z_iatom;
  *x_j = x_iatom - dr[0];
  *y_j = y_iatom - dr[1];
  *z_j = z_iatom - dr[2];

  return 0;
}

/**
 * @details Function for obtaining the angle between two input vectors
 * (std::vector). Internally, the vectors are converted to Eigen vectors. The
 * dot product between the input vectors is used to calculate the angle between
 * them.
 * @param[in] OO The O--O vector (but can be any vector, in general).
 * @param[in] OH The O-H vector (but can be any vector, in general).
 * @return The output angle between the input vectors, in radians
 */
double gen::eigenVecAngle(std::vector<double> OO, std::vector<double> OH) {
  Eigen::Vector3d eigOO = Eigen::Map<Eigen::Vector3d>(OO.data(), OO.size());
  Eigen::Vector3d eigOH = Eigen::Map<Eigen::Vector3d>(OH.data(), OH.size());
  double cosAngle = eigOO.dot(eigOH) / (eigOO.norm() * eigOH.norm());
  // Clamp to [-1, 1] to guard against floating-point imprecision in acos
  cosAngle = std::clamp(cosAngle, -1.0, 1.0);
  return acos(cosAngle);
}

/**
 * @details Get the average, after excluding
 *  the outliers, using quartiles
 * @param[in] inpVec The vector containing values whose average is desired
 * @return The desired average value
 */
double gen::getAverageWithoutOutliers(std::vector<double> inpVec) {
  double avgVal = 0.0;
  const int n = static_cast<int>(inpVec.size());
  if (n == 0) {
    return 0.0;
  }
  std::vector<double> sorted = inpVec;
  std::sort(sorted.begin(), sorted.end());
  if (n < 4) {
    double sumVal = 0.0;
    for (double v : sorted) {
      sumVal += v;
    }
    return sumVal / n;
  }
  std::vector<double> lowerRange;
  std::vector<double> upperRange;
  double firstQuartile, thirdQuartile;
  double iqr;
  double outlierLimLow, outlierLimHigh;
  int numOfObservations = 0;
  if (n % 2 == 0) {
    for (int i = 0; i < n / 2; i++) {
      lowerRange.push_back(sorted[i]);
      upperRange.push_back(sorted[n / 2 + i]);
    }
  } else {
    const int mid = n / 2;
    for (int i = 0; i < mid; i++) {
      lowerRange.push_back(sorted[i]);
      upperRange.push_back(sorted[mid + 1 + i]);
    }
  }
  // ----------------------
  // Calculate the first and third quartiles, and interquartile range
  //
  // First quartile
  firstQuartile = calcMedian(&lowerRange);
  // Third quartile
  thirdQuartile = calcMedian(&upperRange);
  // Interquartile range
  iqr = thirdQuartile - firstQuartile;
  // ----------------------
  // Calculate the average without outliers
  // Outliers are defined as values which
  // are less than Q1 - 1.5IQR
  // or greater than Q3 + 1.5IQR
  //
  // Get the limits for the outliers
  outlierLimLow = firstQuartile - 1.5 * iqr;
  outlierLimHigh = thirdQuartile + 1.5 * iqr;
  //
  // Loop through the values in inpVec to get the average, excluding outliers
  for (int i = 0; i < n; i++) {
    if (sorted[i] < outlierLimLow || sorted[i] > outlierLimHigh) {
      continue;
    }
    numOfObservations++;
    avgVal += sorted[i];
  }
  // ----------------------
  // This fails if there are not enough observations (ring size = 3)
  if (numOfObservations == 0) {
    double sumVal = 0.0;
    // Loop through all the values and sum
    for (int i = 0; i < n; i++) {
      sumVal += inpVec[i];
    } // end of sum
    // Normalize
    avgVal = sumVal / n;
    return avgVal;
  } // for triangular prism blocks
  // ----------------------
  // Divide by the number of observations used
  avgVal /= numOfObservations;

  return avgVal;
}

/**
 * @details Function for getting the angular distance between two quaternions.
 * Returns the result in degrees.
 * @param[in] quat1 The first quaternion
 * @param[in] quat2 The second quaternion
 * @return The output angle between the input quaternions, in degrees
 */
double gen::angDistDegQuaternions(std::vector<double> quat1,
                                  std::vector<double> quat2) {
  //
  double prod; // Product of quat1 and conjugate of quat2
  // The angular distance is
  // angularDistance = 2*cosInverse(quat1*conj(quat2))
  prod = quat1[0] * quat2[0] - quat1[1] * quat2[1] - quat1[2] * quat2[2] -
         quat1[3] * quat2[3];
  // Clamp to [-1, 1] to guard against floating-point imprecision in acos
  prod = std::clamp(prod, -1.0, 1.0);
  // The angular distance is:
  double angDist = 2 * acos(prod) * 180.0 / (gen::pi);
  // Return the angular distance
  return angDist;
}
