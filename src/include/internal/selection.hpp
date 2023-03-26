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
void moleculesInSingleSlice(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    bool clearPreviousSliceSelection=true,
    std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
    std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

//! Given a pointCloud set the inSlice bool for every atom,
//! if the atoms are inside the specified (single) region. 
//! Does not handle atoms in molecules straddling the boundary
void atomsInSingleSlice(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    bool clearPreviousSliceSelection=true,
    std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
    std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

//! Given a particular molecule ID and a pointCloud set the inSlice bool for all atoms,
//! with that molecule ID 
void setAtomsWithSameMolID(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::unordered_multimap<int,int> molIDAtomIDmap,
    int molID, bool inSliceValue=true);

}  // namespace gen

/**
 *  @addtogroup ring
 *  @{
 */

namespace ring {

//! Select edge molecules and atoms which are part of rings, such that rings
//! formed with even one atom in the slice will be included in the selection
//! Modifies the inSlice bool of a given PointCloud (this may be the same)
//! as the given oxygen atom PointCloud which was used to construct the 
//! neighbour list used to construct the rings vector of vectors.
//! We assume that the PointCloud structs have the inSlice bool values set according
//! to the presence of the atom in the slice 
//! (this can be done using the gen::moleculesInSingleSlice function. 
void getEdgeMoleculesInRings(
    std::vector<std::vector<int>> rings, molSys::PointCloud<molSys::Point<double>, double> *oCloud,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, 
    std::array<double, 3> coordLow, std::array<double, 3> coordHigh,
    bool identicalCloud=false);

//! Master function for selecting edge molecules and atoms which are part of rings, such that rings
//! formed with even one atom in the slice will be included in the selection
//! Modifies the inSlice bool of a given PointCloud (this may be the same)
//! as the given oxygen atom PointCloud which was used to construct the 
//! neighbour list used to construct the rings vector of vectors (calls ring::getEdgeMoleculesInRings)
//! Prints out molecule IDs individually of molecules in the slice, and also prints out a LAMMPS
//! data file of just the molecules and atoms in the slice  
void printSliceGetEdgeMoleculesInRings(
    std::string path, std::vector<std::vector<int>> rings, 
    molSys::PointCloud<molSys::Point<double>, double> *oCloud,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, 
    std::array<double, 3> coordLow, std::array<double, 3> coordHigh,
    bool identicalCloud=false);

}  // namespace ring

#endif  // __SELECTION_H_
