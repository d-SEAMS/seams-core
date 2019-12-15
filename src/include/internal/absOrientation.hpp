#ifndef __ABSORIENTATION_H_
#define __ABSORIENTATION_H_

#include <math.h>
#include <sys/stat.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// External
#include <Spectra/SymEigsShiftSolver.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <mol_sys.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

// Inspired by AStar_Dual_Tree_HandPose by jsupancic

namespace absor {

// Get the absolute orientation using Horn's algorithm (with quaternions)
int hornAbsOrientation(const Eigen::MatrixXd& refPoints,
                       const Eigen::MatrixXd& targetPoints);

// Compute the matrix S, or M, whose elements are the sums of products of
// coordinates measured in the left and right systems
Eigen::MatrixXd calcMatrixS(const Eigen::MatrixXd& centeredRefPnts,
                            const Eigen::MatrixXd& centeredTargetPnts, int nop,
                            int dim);

// Compute the matrix N, a 4x4 symmetric matrix, by combining sums (saved as
// elements in the matrix S)
Eigen::MatrixXd calcMatrixN(const Eigen::MatrixXd& S);

// Center a point set wrt the centroid
Eigen::MatrixXd centerWRTcentroid(const Eigen::MatrixXd& pointSet);

}  // namespace absor

#endif  // __ABSORIENTATION_H_
