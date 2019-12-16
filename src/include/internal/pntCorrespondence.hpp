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

// Get the relative ordering of a pair of basal rings for a deformed
// prism/perfect prism. Outputs a vector of vectors of indices, such that the
// first vector is for the first basal ring, and the second vector is for the
// second basal ring. The input neighbour list is with respect to indices, not
// IDs
int relOrderPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2,
    std::vector<std::vector<int>> nList, std::vector<int> *outBasal1,
    std::vector<int> *outBasal2);

}  // namespace pntToPnt

#endif  // __PNTCORRESPONDENCE_H_
