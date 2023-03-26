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

#include <absOrientation.hpp>

/**
 * @details Get the absolute orientation using Horn's algorithm (with
 *  quaternions). The algorithm is described in the [linked
 *  paper](http://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf).
 *  @param[in] refPoints The point set of the reference system (or right
 *   system). This is a @f$ (n \times 3) @f$ Eigen matrix. Here, @f$ n @f$ is
 *   the number of particles.
 *  @param[in] targetPoints @f$ (n \times 3) @f$ Eigen matrix of the
 *   candidate/test system (or left system).
 *  @param[in, out] quat The quaternion for the optimum rotation of the left
 *   system into the right system.
 *  @param[in, out] scale The scale factor obtained from Horn's algorithm.
 */
int absor::hornAbsOrientation(const Eigen::MatrixXd &refPoints,
                              const Eigen::MatrixXd &targetPoints,
                              std::vector<double> *quat, double *rmsd,
                              std::vector<double> *rmsdList, double *scale) {
  int nop =
      refPoints.rows(); // Number of particles (equal to the number of rows)
  int dim =
      refPoints.cols(); // Number of dimensions (equal to the number of columns)
  Eigen::MatrixXd centeredRefPnts(
      nop, dim); // Reference point set after centering wrt the centroid
  Eigen::MatrixXd centeredTargetPnts(
      nop, dim); // Target point set after centering wrt the centroid
  Eigen::MatrixXd S(dim,
                    dim); // Matrix containing sums of products of coordinates
  Eigen::MatrixXd N(
      4, 4); // 4x4 Matrix, whose largest eigenvector must be calculated
  Eigen::VectorXd calcEigenVec(
      4); // This should have 4 components (eigen vector calculated from N)
  // -----
  // Check that the sizes of the reference point set (right point system) and
  // the target point set (left point system) are the same
  if (refPoints.rows() != targetPoints.rows() ||
      refPoints.cols() != targetPoints.cols()) {
    // Throw error
    std::cerr
        << "The reference and target point sets are not of the same size.\n";
    return 1;
  } // unequal size; error!
  // -----
  //
  // ---------------------------------------------------
  // FINDING THE CENTROIDS AND THE NEW COORDINATES WRT THE CENTROIDS
  centeredRefPnts = absor::centerWRTcentroid(refPoints);
  centeredTargetPnts = absor::centerWRTcentroid(targetPoints);
  // ---------------------------------------------------
  // FINDING THE ROTATION MATRIX
  //
  // Find the 3x3 matrix S
  S = absor::calcMatrixS(centeredRefPnts, centeredTargetPnts, nop, dim);
  // Calculate the 4x4 symmetric matrix N, whose largest eigenvector yields the
  // quaternion in the same direction
  N = absor::calcMatrixN(S);
  // --------
  // Calculate the eigenvector corresponding to the largest eigenvalue
  //
  // Construct matrix operation object (op) using the wrapper class
  // DenseSymMatProd for the matrix N
  Spectra::DenseSymMatProd<double> op(N);
  //
  // Construct eigen solver object, requesting the largest 1 eigenvalue and
  // eigenvector
  Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE,
                         Spectra::DenseSymMatProd<double>>
      eigs(&op, 1, 4);
  //
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  // Get the eigenvalue and eigenvector
  if (eigs.info() == Spectra::SUCCESSFUL) {
    Eigen::VectorXd calcEigenValue = eigs.eigenvalues(); // Eigenvalue
    calcEigenVec = eigs.eigenvectors();
  } // end of eigenvector calculation
  //
  // --------
  // Normalize the eigenvector calculated
  double qNorm = sqrt(calcEigenVec.dot(calcEigenVec));
  calcEigenVec /= qNorm; // Divide by the square root of the sum
  // Update the quaternion with the normalized eigenvector
  (*quat).resize(4); // Output quaternion update
  for (int i = 0; i < 4; i++) {
    (*quat)[i] = calcEigenVec(i);
  } // end of quaternion update
  // --------
  // ---------------------------------------------------
  // COMPUTE THE OPTIMUM SCALE
  (*scale) = absor::calcScaleFactor(centeredRefPnts, centeredTargetPnts, nop);
  // ---------------------------------------------------
  // GETTING THE ERROR
  (*rmsd) = absor::getRMSD(centeredRefPnts, centeredTargetPnts, calcEigenVec,
                           rmsdList, nop, (*scale));
  // ---------------------------------------------------
  return 0;
} // end of function

/**
 * @details Compute the matrix S, or M, whose elements are the sums of products
 *  of coordinates measured in the left and right systems. The matrix S is :
 *  S=[Sxx  Sxy      Sxz;...
 *     Syx  Syy      Syz;...
 *     Szx  Szy      Szz];
 *  @param[in] centeredRefPnts The point set of the reference system (or right
 *   system), centered with respect to the centroid. This is a @f$ (n \times 3)
 *   @f$ Eigen matrix. Here, @f$ n @f$ is the number of particles.
 *  @param[in] centeredTargetPnts @f$ (n \times 3) @f$ Eigen matrix of the
 *   candidate/test system (or left system), centered with respect to the
 *   centroid.
 *  @param[in] nop The number of particles.
 *  @param[in] dim The number of dimensions = 3.
 *  @return a @f$ (3 \times 3) @f$ Eigen matrix, or S.
 */
Eigen::MatrixXd absor::calcMatrixS(const Eigen::MatrixXd &centeredRefPnts,
                                   const Eigen::MatrixXd &centeredTargetPnts,
                                   int nop, int dim) {
  Eigen::MatrixXd S(nop, dim);      // Output matrix S
  Eigen::VectorXd targetCoord(nop); // Column of the target point set
  Eigen::VectorXd refCoord(nop);    // Column of the reference point set
  double Svalue; // Current value being filled (Sxx, Sxy etc.)

  // Calculate Sxx, Sxy, Sxz etc
  for (int iCol = 0; iCol < dim; iCol++) {
    //
    for (int jCol = 0; jCol < dim; jCol++) {
      targetCoord =
          centeredTargetPnts.col(iCol); // iCol^th column of target point set
      refCoord = centeredRefPnts.col(jCol); // jCol^th of reference point set
      Svalue = targetCoord.dot(refCoord);
      S(iCol, jCol) = Svalue;
    } // end column wise filling
  }   // end of filling

  // Output matrix
  return S;
} // end of function

/**
 * @details Compute the matrix @f$(4 \times 4)@f$ N, whose largest eigenvector
 *  corresponds to the optimal rotation. It is calculated from the elements of
 *  the @f$(3 \times 3)@f$ matrix S. The matrix N is:
 *
 * N=[(Sxx+Syy+Szz)  (Syz-Szy) (Szx-Sxz) (Sxy-Syx);...
 *
 *    (Syz-Szy)      (Sxx-Syy-Szz)  (Sxy+Syx)      (Szx+Sxz);...
 *
 *    (Szx-Sxz)      (Sxy+Syx)     (-Sxx+Syy-Szz)  (Syz+Szy);...
 *
 *    (Sxy-Syx)      (Szx+Sxz)      (Syz+Szy)      (-Sxx-Syy+Szz)];
 *
 *  @param[in] S (3 \times 3) Eigen matrix whose elements are the sums of
 *   products of the coordinates measured in the left and right systems.
 *  @return a @f$ (4 \times 4) @f$ Eigen matrix, or N.
 */
Eigen::MatrixXd absor::calcMatrixN(const Eigen::MatrixXd &S) {
  //
  Eigen::MatrixXd N(4, 4); // Output matrix N
  // Components of S
  double Sxx = S(0, 0);
  double Sxy = S(0, 1);
  double Sxz = S(0, 2);
  double Syx = S(1, 0);
  double Syy = S(1, 1);
  double Syz = S(1, 2);
  double Szx = S(2, 0);
  double Szy = S(2, 1);
  double Szz = S(2, 2);

  // N=[(Sxx+Syy+Szz)  (Syz-Szy)      (Szx-Sxz)      (Sxy-Syx);...
  //          (Syz-Szy)      (Sxx-Syy-Szz)  (Sxy+Syx)      (Szx+Sxz);...
  //          (Szx-Sxz)      (Sxy+Syx)     (-Sxx+Syy-Szz)  (Syz+Szy);...
  //          (Sxy-Syx)      (Szx+Sxz)      (Syz+Szy)      (-Sxx-Syy+Szz)];

  // ------------------
  // Calculating N:
  // Diagonal elements
  N(0, 0) = Sxx + Syy + Szz;
  N(1, 1) = Sxx - Syy - Szz;
  N(2, 2) = -Sxx + Syy - Szz;
  N(3, 3) = -Sxx - Syy + Szz;
  // Other elements
  // First row
  N(0, 1) = Syz - Szy;
  N(0, 2) = Szx - Sxz;
  N(0, 3) = Sxy - Syx;
  // Second row
  N(1, 0) = Syz - Szy;
  N(1, 2) = Sxy + Syx;
  N(1, 3) = Szx + Sxz;
  // Third row
  N(2, 0) = Szx - Sxz;
  N(2, 1) = Sxy + Syx;
  N(2, 3) = Syz + Szy;
  // Fourth row
  N(3, 0) = Sxy - Syx;
  N(3, 1) = Szx + Sxz;
  N(3, 2) = Syz + Szy;
  // ------------------
  // Output matrix
  return N;
} // end of function

/**
 *  @details Centers a point set (which is an Eigen matrix),
 *   with respect to the centroid.
 *  @param[in] pointSet @f$ (n \times 3) @f$ Eigen matrix for the point set.
 *   Here @f$ n @f$ is the number of points.
 *  @return a @f$ (n \times 3) @f$ Eigen matrix of the same size as the input
 *   point set.
 */
Eigen::MatrixXd absor::centerWRTcentroid(const Eigen::MatrixXd &pointSet) {
  int nop = pointSet.rows();                  // Number of particles
  int dim = pointSet.cols();                  // Number of dimensions
  Eigen::MatrixXd centeredPointSet(nop, dim); // Output point set
  Eigen::VectorXd vecOfOnes(nop);             // vector of ones
  std::vector<double> centroid;
  double coordValue;
  double centeredVal;
  //
  centroid.resize(dim); // Init to zero
  vecOfOnes = Eigen::VectorXd::Ones(nop);
  // --------------------------------

  for (int i = 0; i < nop; i++) {
    for (int k = 0; k < dim; k++) {
      coordValue = pointSet(i, k);
      centroid[k] += coordValue;
    } // loop through columns
  }   // end of loop through rows
  // Divide by the total number of particles
  centroid[0] /= nop; // x
  centroid[1] /= nop; // y
  centroid[2] /= nop; // z
  // --------------------------------
  // Subtract the centroid from the coordinates to get the centered point set
  for (int i = 0; i < nop; i++) {
    for (int k = 0; k < dim; k++) {
      coordValue = pointSet(i, k);
      centeredVal = coordValue - centroid[k];
      centeredPointSet(i, k) = centeredVal;
    } // end of loop through columns (dimensions)
  }   // end of loop through the rows
  // --------------------------------
  return centeredPointSet;
} // end of function

/**
 *  @details Calculate the scale factor from the centered right
 *   and left point sets.
 *  @param[in] centeredRefPnts The point set of the reference system (or right
 *   system), centered with respect to the centroid. This is a @f$ (n \times 3)
 *   @f$ Eigen matrix. Here, @f$ n @f$ is the number of particles.
 *  @param[in] centeredTargetPnts @f$ (n \times 3) @f$ Eigen matrix of the
 *   candidate/test system (or left system), centered with respect to the
 *   centroid.
 *  @param[in] n The number of points.
 *  @return the scale factor.
 */
double absor::calcScaleFactor(const Eigen::MatrixXd &rightSys,
                              const Eigen::MatrixXd &leftSys, int n) {
  double scale;  // Output scale
  double v1, v2; // Sum of the length of the vector
  Eigen::VectorXd rightVec(
      3);                     // Vector of the i^th particle in the right system
  Eigen::VectorXd leftVec(3); // Vector of the i^th particle in the right system

  // scale = (sigma_to_n ||r_r||^2 / ||r_l||^2)^0.5
  // ref: http://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf

  // Loop through all the points, and get the dot product of the vector for each
  // point
  for (int i = 0; i < n; i++) {
    //
    rightVec = rightSys.row(i); // i^th row of the right system
    leftVec = leftSys.row(i);   // i^th row of the left system
    v1 += rightVec.dot(rightVec);
    v2 += leftVec.dot(leftVec);
  } // end of loop through all points

  // The optimum scale is the ratio of v1 and v2
  scale = std::sqrt(v1 / v2);

  return scale;
} // end of function

/**
 *  @details Get a \f$ (3 \times 3) \f$ rotation matrix from
 *   a unit quaternion.
 *  @param[in] quat The Eigen vector of length 4, for the input quaternion.
 *  @return the rotation matrix.
 */
Eigen::MatrixXd absor::quat2RotMatrix(const Eigen::VectorXd &quat) {
  //
  Eigen::MatrixXd R(3, 3); // Rotation matrix
  // Components of the quaternion
  double q0 = quat(0);
  double qx = quat(1);
  double qy = quat(2);
  double qz = quat(3);

  // Quaternion derived rotation matrix, when q (q=q0+qx*i+qy*j+qz*k)
  // is a unit quaternion:
  // R=[(q0^2+qx^2+qy^2+qz^2)      2(qx*qy-q0*qz)        2(qx*qz+q0*qy);...
  //    2(qy*qx-q0*qz)            (q0^2+qx^2+qy^2+qz^2)  2(qy*qz-q0*qx);...
  //    2(qz*qx-q0*qy)            2(qz*qy+q0*qz)         (q0^2+qx^2+qy^2+qz^2)];

  // Fill up the rotation matrix R according to the above formula
  //
  // First row
  R(0, 0) = q0 * q0 + qx * qx + qy * qy + qz * qz;
  R(0, 1) = 2 * (qx * qy - q0 * qz);
  R(0, 2) = 2 * (qx * qz + q0 * qy);
  // Second row
  R(1, 0) = 2 * (qy * qx + q0 * qz);
  R(1, 1) = q0 * q0 + qx * qx + qy * qy + qz * qz;
  R(1, 2) = 2 * (qy * qz - q0 * qy);
  // Third row
  R(2, 0) = 2 * (qz * qx - q0 * qy);
  R(2, 1) = 2 * (qz * qy + q0 * qz);
  R(2, 2) = q0 * q0 + qx * qx + qy * qy + qz * qz;

  // return the rotation matrix
  return R;
} // end of function

/**
 *  @details Get the root mean square of the errors from Horn's absolute
 *   orientation algorithm.
 *  @param[in] quat The Eigen vector of length 4, for the input quaternion.
 *  @param[in] centeredRefPnts The point set of the reference system (or right
 *   system), centered with respect to the centroid. This is a @f$ (n \times 3)
 *   @f$ Eigen matrix. Here, @f$ n @f$ is the number of particles.
 *  @param[in] centeredTargetPnts @f$ (n \times 3) @f$ Eigen matrix of the
 *   candidate/test system (or left system), centered with respect to the
 *   centroid.
 *  @param[in] quat The Eigen vector of length 4, for the quaternion denoting
 *   the rotation.
 *  @param[in] nop The number of particles.
 *  @param[in] scale The scale factor
 *  @return the least RMSD.
 */
double absor::getRMSD(const Eigen::MatrixXd &centeredRefPnts,
                      const Eigen::MatrixXd &centeredTargetPnts,
                      const Eigen::VectorXd &quat,
                      std::vector<double> *rmsdList, int nop, double scale) {
  Eigen::MatrixXd R(3, 3); // The (3x3) rotation vector
  // The RMSD per atom is filled in this vector
  (*rmsdList).resize(nop);
  //
  R = absor::quat2RotMatrix(
      quat); // orthonormal rotation matrix from the quaternion
  // The total error is:
  // sum_over_all_n (r'_r - s*rotated_r_l')
  double rmsd = 0.0;           // Error
  Eigen::VectorXd errorVec(3); // The vector which is the r_r -s*R
  Eigen::VectorXd rotatedLeft(3);
  Eigen::VectorXd targetCol(3);
  Eigen::VectorXd refCol(3);
  //
  for (int i = 0; i < nop; i++) {
    //
    // Rotate the left system coordinate using the rotation matrix
    targetCol = (centeredTargetPnts.row(i)).transpose();
    refCol = (centeredRefPnts.row(i)).transpose();
    //
    rotatedLeft = R * targetCol;
    errorVec = refCol - scale * rotatedLeft;
    rmsd += errorVec.dot(errorVec);                // Total error
    (*rmsdList)[i] = sqrt(errorVec.dot(errorVec)); // Error per atom
  } // end of loop through every row
  //
  return sqrt(rmsd / nop);
} // end of function
