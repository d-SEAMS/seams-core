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
#include <Eigen/Core>
#include <Eigen/Dense>

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
