#include <absOrientation.hpp>

// Get the absolute orientation using Horn's algorithm (with quaternions)
int absor::hornAbsOrientation(const Eigen::MatrixXd& refPoints,
                              const Eigen::MatrixXd& targetPoints) {
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