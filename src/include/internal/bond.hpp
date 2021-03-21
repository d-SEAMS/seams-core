//-----------------------------------------------------------------------------------
// d-SEAMS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef __BONDING_H_
#define __BONDING_H_

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

// Internal
#include <cage.hpp>
#include <mol_sys.hpp>
#include <seams_input.hpp>

/** @file bond.hpp
 *  @brief File for bond-related analyses (hydrogen bonds, bonded atoms for data
 * file write-outs etc.).
 */

/**
 *  @addtogroup bond
 *  @{
 */

/** @brief Functions for bond-related analyses
 *
 *  @details This namespace contains functions that are used for determining
 the hydrogen bonds, or lists of bonded atoms for write-outs to data files and
 visualization.
 *
 * These functions are distinct from the neighbour-list building functions.
 Bonds can be built from cages, ring lists etc.
 *
 * Hydrogen bonds are determined using a strict geometric criterion.
 *
 * A hydrogen bond between two water molecules exists when:
 *
 * 1. The distance between the donor oxygen atom and the acceptor hydrogen atom
 is less than 2.42 Angstrom
 * 2. The angle between the O--O vector and the O-H vector should be less than
 30 degrees
 *
 * ### Changelog ###
 *
 * - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 13, 2019
 * - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

namespace bond {

//! Create a vector of vectors (similar to the neighbour list conventions)
//! containing the hydrogen bond connectivity information. Decides the
//! existence of the hydrogen bond depending on the O--O and O--H vectors from
//! the neighbour list already constructed
std::vector<std::vector<int>>
populateHbonds(std::string filename,
               molSys::PointCloud<molSys::Point<double>, double> *yCloud,
               std::vector<std::vector<int>> nList, int targetFrame, int Htype);

//! Calculates the distance of the hydrogen bond between O and H (of different
//! atoms), given the respective pointClouds and the indices to each atom
double
getHbondDistanceOH(molSys::PointCloud<molSys::Point<double>, double> *oCloud,
                   molSys::PointCloud<molSys::Point<double>, double> *hCloud,
                   int oAtomIndex, int hAtomIndex);

//! Create a vector of vectors containing bond connectivity information. May
//! contain duplicates! Gets the bond information from the vector of vectors
//! containing the neighbour list by index
std::vector<std::vector<int>>
populateBonds(std::vector<std::vector<int>> nList,
              molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Create a vector of vectors containing bond connectivity information
//! Gets the bond information from the vector of vectors
//! containing the neighbour list by index. Bonds between dummy atoms are not
//! filled.
std::vector<std::vector<int>>
populateBonds(std::vector<std::vector<int>> nList,
              molSys::PointCloud<molSys::Point<double>, double> *yCloud,
              std::vector<cage::iceType> atomTypes);

//! Creates a vector of vectors containing bond connectivity information from
//! the rings vector of vectors and cage information
std::vector<std::vector<int>>
createBondsFromCages(std::vector<std::vector<int>> rings,
                     std::vector<cage::Cage> *cageList, cage::cageType type,
                     int *nRings);

//! Remove duplicate bonds
std::vector<std::vector<int>> trimBonds(std::vector<std::vector<int>> bonds);

} // namespace bond

#endif // __BONDING_H_
