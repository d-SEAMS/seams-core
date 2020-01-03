#include <generic.hpp>
#include <iostream>

/********************************************/ /**
*  Function for printing out
 info in a PointCloud object.
 *  @param[in] yCloud The input PointCloud to be printed.
 *  @param[in] outFile The name of the output file to which the information will
be printed.
***********************************************/
int gen::prettyPrintYoda(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::string outFile) {
  std::ofstream outputFile;
  // Create a new file in the output directory
  outputFile.open(outFile);

  if (outputFile.is_open()) {
    // First line
    outputFile << "# Frame\tAtomID\tx\ty\tz\tcij\ticeType\n";
    // Write out all the information out line by line
    for (int i = 0; i < yCloud->nop; i++) {
      outputFile << yCloud->currentFrame << "\t" << yCloud->pts[i].atomID
                 << "\t" << yCloud->pts[i].x << "\t" << yCloud->pts[i].y << "\t"
                 << yCloud->pts[i].z << "\t";
      // Print out cij
      // for(int c=0; c<yCloud->pts[i].c_ij.size(); c++){outputFile <<
      // yCloud->pts[i].c_ij[c]<<"\t";} Print out the classifier
      outputFile << yCloud->pts[i].iceType << "\n";
    }
  }
  // Close the file
  outputFile.close();
  return 0;
}

/********************************************/ /**
*  Function for getting the unwrapped coordinates
 of a pair of atoms.
 *  @param[in] yCloud The input PointCloud to be printed.
 *  @param[in] iatomIndex Index of the \f$ i^{th} \f$ atom.
 *  @param[in] jatomIndex Index of the \f$ j^{th} \f$ atom.
 *  @param[in, out] x_i X Coordinate of the \f$ i^{th} \f$ atom corresponding to
the unwrapped distance.
 *  @param[in, out] y_i Y Coordinate of the \f$ i^{th} \f$ atom corresponding to
the unwrapped distance.
 *  @param[in, out] z_i Z Coordinate of the \f$ i^{th} \f$ atom corresponding to
the unwrapped distance.
 *  @param[in, out] x_j X Coordinate of the \f$ j^{th} \f$ atom corresponding to
the unwrapped distance.
 *  @param[in, out] y_j Y Coordinate of the \f$ j^{th} \f$ atom corresponding to
the unwrapped distance.
 *  @param[in, out] z_j Z Coordinate of the \f$ j^{th} \f$ atom corresponding to
the unwrapped distance.
***********************************************/
int gen::unwrappedCoordShift(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatomIndex,
    int jatomIndex, double *x_i, double *y_i, double *z_i, double *x_j,
    double *y_j, double *z_j) {
  //
  double x_iatom, y_iatom, z_iatom;
  double x_jatom, y_jatom, z_jatom;
  double x_ij, y_ij, z_ij;  // Relative distance
  std::vector<double> box = yCloud->box;
  double xPBC, yPBC, zPBC;  // Actual unwrapped distance

  // ----------------------------------------------------------------------
  // INIT
  // iatom
  x_iatom = yCloud->pts[iatomIndex].x;
  y_iatom = yCloud->pts[iatomIndex].y;
  z_iatom = yCloud->pts[iatomIndex].z;
  // jatom
  x_jatom = yCloud->pts[jatomIndex].x;
  y_jatom = yCloud->pts[jatomIndex].y;
  z_jatom = yCloud->pts[jatomIndex].z;
  // ----------------------------------------------------------------------
  // GET RELATIVE DISTANCE
  x_ij = x_iatom - x_jatom;
  y_ij = y_iatom - y_jatom;
  z_ij = z_iatom - z_jatom;
  // ----------------------------------------------------------------------
  // SHIFT COORDINATES IF REQUIRED
  // Shift x
  if (fabs(x_ij) > 0.5 * box[0]) {
    // Get the actual distance
    xPBC = box[0] - fabs(x_ij);
    if (x_ij < 0) {
      x_jatom = x_iatom - xPBC;
    }  // To the -x side of currentIndex
    else {
      x_jatom = x_iatom + xPBC;
    }  // Add to the + side
  }    // Shift nextElement
  //
  // Shift y
  if (fabs(y_ij) > 0.5 * box[1]) {
    // Get the actual distance
    yPBC = box[1] - fabs(y_ij);
    if (y_ij < 0) {
      y_jatom = y_iatom - yPBC;
    }  // To the -y side of currentIndex
    else {
      y_jatom = y_iatom + yPBC;
    }  // Add to the + side
  }    // Shift nextElement
  //
  // Shift z
  if (fabs(z_ij) > 0.5 * box[2]) {
    // Get the actual distance
    zPBC = box[2] - fabs(z_ij);
    if (z_ij < 0) {
      z_jatom = z_iatom - zPBC;
    }  // To the -z side of currentIndex
    else {
      z_jatom = z_iatom + zPBC;
    }  // Add to the + side
  }    // Shift nextElement
  // ----------------------------------------------------------------------
  // Assign values
  *x_i = x_iatom;
  *y_i = y_iatom;
  *z_i = z_iatom;
  *x_j = x_jatom;
  *y_j = y_jatom;
  *z_j = z_jatom;

  return 0;
}

/********************************************/ /**
 *  Function for obtaining the angle between two input vectors (std::vector).
 Internally, the vectors are converted to GSL vectors. The dot product between
 the input vectors is used to calculate the angle between them.
 *  @param[in] OO The O--O vector (but can be any vector, in general).
 *  @param[in] OH The O-H vector (but can be any vector, in general).
 *  \return The output angle between the input vectors, in radians
 ***********************************************/
double gen::gslVecAngle(std::vector<double> OO, std::vector<double> OH) {
  gsl_vector *gOO = gsl_vector_alloc(3);
  gsl_vector *gOH = gsl_vector_alloc(3);
  double norm_gOO, norm_gOH, xDummy, angle;
  for (int i = 0; i < 3; i++) {
    gsl_vector_set(gOO, i, OO[i]);
    gsl_vector_set(gOH, i, OH[i]);
  }
  norm_gOO = gsl_blas_dnrm2(gOO);
  norm_gOH = gsl_blas_dnrm2(gOH);
  gsl_blas_ddot(gOO, gOH, &xDummy);
  angle = acos(xDummy / (norm_gOO * norm_gOH));
  gsl_vector_free(gOO);
  gsl_vector_free(gOH);
  return angle;
}

/********************************************/ /**
 *  Get the average, after excluding
 the outliers, using quartiles
 *  @param[in] inpVec The vector containing values whose average is desired
 *  \return The desired average value
 ***********************************************/
double gen::getAverageWithoutOutliers(std::vector<double> inpVec) {
  //
  double avgVal = 0.0;    // Average value, excluding the outliers
  double median;          // Median value
  int n = inpVec.size();  // Number of values
  std::vector<double> lowerRange,
      upperRange;                        // n/2 smallest and n/2 largest numbers
  double firstQuartile, thirdQuartile;   // First and third quartiles
  double iqr;                            // Interquartile range
  double outlierLimLow, outlierLimHigh;  // Outliers limit
  int numOfObservations = 0;  // Number of observations used for the average
  // ----------------------
  // Calculate the median (the vector is sorted inside the function)
  median = calcMedian(&inpVec);
  // ----------------------
  // Get the n/2 smallest and largest numbers
  //
  if (n % 2 == 0) {
    for (int i = 0; i < n / 2; i++) {
      // n/2 smallest numbers
      lowerRange.push_back(inpVec[i]);
      // n/2 largest numbers
      upperRange.push_back(inpVec[n / 2 + i]);
    }  // end of loop to fill up the n/2 smallest and n/2 largest
  }    // even
  else {
    //
    int halfN = (n + 1) / 2;
    // Exclude the median
    for (int i = 0; i < halfN; i++) {
      // (n+1)/2 smallest numbers
      lowerRange.push_back(inpVec[i]);
      // (n+1)/2 largest numbers
      upperRange.push_back(inpVec[halfN + i]);
    }  // end of filling up the smallest and largest half-ranges
  }    // for odd numbers
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
    //
    if (inpVec[i] < outlierLimLow) {
      continue;
    }  // lower limit outlier
    else if (inpVec[i] > outlierLimHigh) {
      continue;
    }  // higher limit outlier
    else {
      // Number of observations added
      numOfObservations++;
      // Add to the average
      avgVal += inpVec[i];
    }  // take the average
  }    // end of loop for getting the average
  //
  // Divide by the number of observations used
  avgVal /= numOfObservations;
  // ----------------------

  return avgVal;
}

/********************************************/ /**
 *  Function for getting the angular distance between two quaternions. Returns
 the result in degrees.
 *  @param[in] quat1 The first quaternion
 *  @param[in] quat2 The second quaternion
 *  \return The output angle between the input quaternions, in degrees
 ***********************************************/
double gen::angDistDegQuaternions(std::vector<double> quat1,
                                  std::vector<double> quat2) {
  //
  double prod;  // Product of quat1 and conjugate of quat2
  // The angular distance is
  // angularDistance = 2*cosInverse(quat1*conj(quat2))
  prod = quat1[0] * quat2[0] - quat1[1] * quat2[1] - quat1[2] * quat2[2] -
         quat1[3] * quat2[3];
  // The angular distance is:
  double angDist = 2 * acos(prod) * 180.0 / (gen::pi);
  // Return the angular distance
  return angDist;
}