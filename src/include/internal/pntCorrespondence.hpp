#ifndef __PNTCORRESPONDENCE_H_
#define __PNTCORRESPONDENCE_H_

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
#include <eigen3/Eigen/Core>

#include <cage.hpp>
#include <mol_sys.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

namespace pntToPnt {

// Fills up an eigen matrix point set for an HC, according to an input
// pointCloud, the relative order given by the basal rings, beginning from the
// startingIndex
Eigen::MatrixXd fillPointSetHC(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> relOrder, int startingIndex);

// Fills up an eigen matrix point set a reference ring, which is a regular
// n-gonal polygon; where n is the number of nodes in the ring
Eigen::MatrixXd getPointSetRefRing(int n);

}  // namespace pntToPnt

#endif  // __PNTCORRESPONDENCE_H_
