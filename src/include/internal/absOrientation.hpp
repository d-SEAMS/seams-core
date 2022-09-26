//-----------------------------------------------------------------------------------
// d-SEAMS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef __ABSORIENTATION_H_
#define __ABSORIENTATION_H_

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <math.h>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

// External
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <mol_sys.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

// Inspired by AStar_Dual_Tree_HandPose by jsupancic

namespace absor {

/**
 * @brief Get the absolute orientation using Horn's algorithm (with quaternions)
 */
int hornAbsOrientation(const Eigen::MatrixXd &refPoints,
                       const Eigen::MatrixXd &targetPoints,
                       std::vector<double> *quat, double *rmsd,
                       std::vector<double> *rmsdList, double *scale);

/**
 * @brief Get the absolute orientation using Horn's algorithm (with quaternions)
 * Row major Eigen matrices taken as input
 */
int hornAbsOrientationRowMajor(const Eigen::MatrixXdRowMajor &refPoints,
                       const Eigen::MatrixXdRowMajor &targetPoints,
                       std::vector<double> *quat, double *rmsd,
                       std::vector<double> *rmsdList, double *scale);

//! Compute the matrix S, or M, whose elements are the sums of products of
//! coordinates measured in the left and right systems
Eigen::MatrixXd calcMatrixS(const Eigen::MatrixXd &centeredRefPnts,
                            const Eigen::MatrixXd &centeredTargetPnts, int nop,
                            int dim);

//! Compute the matrix N, a 4x4 symmetric matrix, by combining sums (saved as
//! elements in the matrix S)
Eigen::MatrixXd calcMatrixN(const Eigen::MatrixXd &S);

//! Center a point set wrt the centroid
Eigen::MatrixXd centerWRTcentroid(const Eigen::MatrixXd &pointSet);

//! Center a point set (row-major) wrt the centroid
Eigen::MatrixXd centerWRTcentroid(const Eigen::MatrixXdRowMajor &pointSet);

//! Calculate the scale factor from the centered left and right point sets
double calcScaleFactor(const Eigen::MatrixXd &rightSys,
                       const Eigen::MatrixXd &leftSys, int n);

//! Get a rotation matrix from a unit quaternion
Eigen::MatrixXd quat2RotMatrix(const Eigen::VectorXd &quat);

//! Calculate the RMSD
double getRMSD(const Eigen::MatrixXd &centeredRefPnts,
               const Eigen::MatrixXd &centeredTargetPnts,
               const Eigen::VectorXd &quat, std::vector<double> *rmsdList,
               int nop, double scale);

} // namespace absor

#endif // __ABSORIENTATION_H_
