//-----------------------------------------------------------------------------------
// d-SEAMS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef __SELECTION_H_
#define __SELECTION_H_

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

#include <mol_sys.hpp>
#include <order_parameter.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

/** @file selection.hpp
 *  @brief File containing functions used to define 'selections'
 *   in a given range, using ring information.
 */

/**
 *  @addtogroup gen
 *  @{
 */

namespace gen {

//! Given a pointCloud containing certain atom types,
//! this returns a pointCloud containing atoms of only the desired type
molSys::PointCloud<molSys::Point<double>, double>
getPointCloudOneAtomType(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    molSys::PointCloud<molSys::Point<double>, double> *outCloud,
    int atomTypeI, bool isSlice = false,
    std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
    std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

//! Given a pointCloud set the inSlice bool for every atom,
//! if the molecules are inside the specified (single) region. 
//! If even one atom of a molecule is inside the region, then all
//! atoms of that molecule will be inside the region (irrespective of type)
molSys::PointCloud<molSys::Point<double>, double>
moleculesInSingleSlice(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    bool clearPreviousSliceSelection=true,
    std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
    std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

}  // namespace gen

/**
 *  @addtogroup ring
 *  @{
 */

// namespace ring {

// //! Find out which rings are prisms.
// //! Returns a vector containing all the ring IDs which are prisms
// std::vector<int> findPrisms(
//     std::vector<std::vector<int>> rings, std::vector<strucType> *ringType,
//     int *nPerfectPrisms, int *nImperfectPrisms,
//     std::vector<std::vector<int>> nList,
//     molSys::PointCloud<molSys::Point<double>, double> *yCloud,
//     std::vector<double> *rmsdPerAtom, bool doShapeMatching = false);

// }  // namespace ring

#endif  // __SELECTION_H_
