//-----------------------------------------------------------------------------------
// d-SEAMS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef __SHAPEMATCH_H_
#define __SHAPEMATCH_H_

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

#include <absOrientation.hpp>
#include <mol_sys.hpp>
#include <pntCorrespondence.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

namespace match {

//! Shape-matching for a pair of polygon basal rings. Returns true if the pair
//! of basal rings form a prism block.
bool matchPrism(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                std::vector<std::vector<int>> nList,
                const Eigen::MatrixXd &refPoints, std::vector<int> *basal1,
                std::vector<int> *basal2, std::vector<double> *rmsdPerAtom,
                bool isPerfect = true);

//! Shape-matching for a pair of polygon basal rings. Returns true if the pair
//! of basal rings form a prism block.
bool matchUntetheredPrism(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, const Eigen::MatrixXd &refPoints,
    std::vector<int> *basal1, std::vector<int> *basal2,
    std::vector<double> *rmsdPerAtom);

//! Shape-matching for a pair of polygon basal rings, comparing with a complete
//! prism block. Returns true if the pair of basal rings form a prism block.
bool matchPrismBlock(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                     std::vector<std::vector<int>> nList,
                     const Eigen::MatrixXd &refPoints, std::vector<int> *basal1,
                     std::vector<int> *basal2, int *beginIndex);

//! Update the per-particle RMSD for a prism block basal ring.
int updatePerAtomRMSDRing(std::vector<int> basalRing, int startingIndex,
                          std::vector<double> rmsdFromMatch,
                          std::vector<double> *rmsdPerAtom);

//! Update the RMSD of each particle in a prism block basal ring with the RMSD
//! of the ring.
int updateRMSDRing(std::vector<int> basalRing, int startingIndex,
                   double rmsdVal, std::vector<double> *rmsdPerAtom);

} // namespace match

#endif // __SHAPEMATCH_H_
