#include <absOrientation.hpp>

// Get the absolute orientation using Horn's algorithm (with quaternions)
int absor::hornAbsOrientation(const Eigen::MatrixXd& refPoints,
                              const Eigen::MatrixXd& targetPoints,
                              std::vector<double>* quat) {
  //
  int nop =
      refPoints.rows();  // Number of particles (equal to the number of rows)
  int dim =
      refPoints
          .cols();  // Number of dimensions (equal to the number of columns)
  Eigen::MatrixXd centeredRefPnts(
      nop, dim);  // Reference point set after centering wrt the centroid
  Eigen::MatrixXd centeredTargetPnts(
      nop, dim);  // Target point set after centering wrt the centroid
  Eigen::MatrixXd S(dim,
                    dim);  // Matrix containing sums of products of coordinates
  Eigen::MatrixXd N(
      4, 4);  // 4x4 Matrix, whose largest eigenvector must be calculated
  Eigen::VectorXd calcEigenVec(
      4);  // This should have 4 components (eigen vector calculated from N)
  // -----
  // Check that the sizes of the reference point set (right point system) and
  // the target point set (left point system) are the same
  if (refPoints.rows() != targetPoints.rows() ||
      refPoints.cols() != targetPoints.cols()) {
    // Throw error
    std::cerr
        << "The reference and target point sets are not of the same size.\n";
    return 1;
  }  // unequal size; error!
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
                         Spectra::DenseSymMatProd<double> >
      eigs(&op, 1, 4);
  //
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  // Get the eigenvalue and eigenvector
  if (eigs.info() == Spectra::SUCCESSFUL) {
    Eigen::VectorXd calcEigenValue = eigs.eigenvalues();  // Eigenvalue
    calcEigenVec = eigs.eigenvectors();
  }  // end of eigenvector calculation
  //
  // --------
  // Normalize the eigenvector calculated
  double qNorm = sqrt(calcEigenVec.dot(calcEigenVec));
  calcEigenVec /= qNorm;  // Divide by the square root of the sum
  // Update the quaternion with the normalized eigenvector
  (*quat).resize(4);  // Output quaternion update
  for (int i = 0; i < 4; i++) {
    (*quat)[i] = calcEigenVec(i);
  }  // end of quaternion update
  // --------
  // ---------------------------------------------------
  return 0;
}  // end of function

// Compute the matrix S, or M, whose elements are the sums of products of
// coordinates measured in the left and right systems
Eigen::MatrixXd absor::calcMatrixS(const Eigen::MatrixXd& centeredRefPnts,
                                   const Eigen::MatrixXd& centeredTargetPnts,
                                   int nop, int dim) {
  //
  Eigen::MatrixXd S(nop, dim);       // Output matrix S
  Eigen::VectorXd targetCoord(nop);  // Column of the target point set
  Eigen::VectorXd refCoord(nop);     // Column of the reference point set
  double Svalue;  // Current value being filled (Sxx, Sxy etc.)

  // Calculate Sxx, Sxy, Sxz etc
  for (int iCol = 0; iCol < dim; iCol++) {
    //
    for (int jCol = 0; jCol < dim; jCol++) {
      targetCoord =
          centeredTargetPnts.col(iCol);  // iCol^th column of target point set
      refCoord = centeredRefPnts.col(jCol);  // jCol^th of reference point set
      Svalue = targetCoord.dot(refCoord);
      S(iCol, jCol) = Svalue;
    }  // end column wise filling
  }    // end of filling

  // Output matrix
  return S;
}  // end of function

// Compute the matrix S, or M, whose elements are the sums of products of
// coordinates measured in the left and right systems
Eigen::MatrixXd absor::calcMatrixN(const Eigen::MatrixXd& S) {
  //
  Eigen::MatrixXd N(4, 4);  // Output matrix N
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
}  // end of function

// Center a point set wrt the centroid
Eigen::MatrixXd absor::centerWRTcentroid(const Eigen::MatrixXd& pointSet) {
  //
  int nop = pointSet.rows();                   // Number of particles
  int dim = pointSet.cols();                   // Number of dimensions
  Eigen::MatrixXd centeredPointSet(nop, dim);  // Output point set
  Eigen::VectorXd vecOfOnes(nop);              // vector of ones
  std::vector<double> centroid;
  double coordValue;
  double centeredVal;
  //
  centroid.resize(dim);  // Init to zero
  vecOfOnes = Eigen::VectorXd::Ones(nop);
  // --------------------------------

  for (int i = 0; i < nop; i++) {
    for (int k = 0; k < dim; k++) {
      coordValue = pointSet(i, k);
      centroid[k] += coordValue;
    }  // loop through columns
  }    // end of loop through rows
  // Divide by the total number of particles
  centroid[0] /= nop;  // x
  centroid[1] /= nop;  // y
  centroid[2] /= nop;  // z
  // --------------------------------
  // Subtract the centroid from the coordinates to get the centered point set
  for (int i = 0; i < nop; i++) {
    for (int k = 0; k < dim; k++) {
      coordValue = pointSet(i, k);
      centeredVal = coordValue - centroid[k];
      centeredPointSet(i, k) = centeredVal;
    }  // end of loop through columns (dimensions)
  }    // end of loop through the rows
  // --------------------------------
  return centeredPointSet;
}  // end of function

// Get a rotation matrix from a unit quaternion
Eigen::MatrixXd absor::quat2RotMatrix(const Eigen::VectorXd& quat) {
  //
  Eigen::MatrixXd R(3, 3);  // Rotation matrix
  // Components of the quaternion
  double q_r = quat(0);
  double q_i = quat(1);
  double q_j = quat(2);
  double q_k = quat(3);

  // Quaternion derived rotation matrix, when q (q=qr+qi*i+qj*j+qk*k)
  // is a unit quaternion:
  // R=[1-2(qj^2+qk^2)      2(qi*qj-qk*qr)      2(qi*qk+qj*qr);...
  //    2(qi*qj-qk*qr)      1-2(qi^2-qk^2)      2(qj*qk-qi*qr);...
  //    2(qi*qk-qj*qr)      2(qj*qk+qi*qr)      1-2(qi^2+qj^2)];

  // Fill up the rotation matrix R according to the above formula
  //
  // First row
  R(0, 0) = 1 - 2 * (q_j * q_j + q_k + q_k);
  R(1, 0) = 2 * (q_i * q_j - q_k * q_r);
  R(2, 0) = 2 * (q_i * q_k + q_j * q_r);
  // Second row
  R(1, 0) = 2 * (q_i * q_j + q_k * q_r);
  R(1, 1) = 1 - 2 * (q_i * q_i + q_k * q_k);
  R(1, 2) = 2 * (q_j * q_k - q_i * q_r);
  // Third row
  R(2, 0) = 2 * (q_i * q_k - q_j * q_r);
  R(2, 1) = 2 * (q_j * q_k + q_i * q_r);
  R(2, 2) = 1 - 2 * (q_i * q_i + q_j * q_j);

  // return the rotation matrix
  return R;
}  // end of function