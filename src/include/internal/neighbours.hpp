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

#ifndef __NEIGHBOURS_H_
#define __NEIGHBOURS_H_

#include <generic.hpp>
#include <mol_sys.hpp>

/** @file neighbours.hpp
 *  @brief Header file for neighbour list generation.
 */

/**
 *  @addtogroup nneigh
 *  @{
 */

/** @brief Functions for building neighbour lists.
 * This namespace contains functions that build neighbour lists (using
 * brute-force), saving either the atom IDs or atom indices (according to a
 * PointCloud) in a row-ordered vector of vectors.
 * Whether the atom IDs or atom indices (i.e. the indices of the elements in the
 * vector pts inside the PointCloud) are saved, the neighbour lists are
 * constructed such that the first element is the 'central atom', whose
 * neighbours are being saved on that particular line. The central atom is
 * followed by the atom IDs or indices of the nearest neighbours.
 *
 * In a 'full' neighbour list, if 1 is a neighbour of 2 then 1 is saved in the
 * neighbour list of 2 AND 2 is also saved in the neighbour list of 1.
 *
 * In a 'half' neighbour list, if 1 is a neighbour of 2, then 2 is saved in the
 * neighbour list of 1 but not vice versa.
 *
 *   ### Changelog ###
 *
 * - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 14, 2019
 * - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

namespace nneigh {
//! All these functions use atom IDs and not indices

//! Inefficient @f$O(n^2)@f$ implementation of neighbour lists when there are
//! two different types of atoms The neighbour list does not differentiate
//! between the types of atoms
std::vector<std::vector<int>> neighList(
    double rcutoff, molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    int typeI, int typeJ);

//! Inefficient @f$O(n^2)@f$ implementation of neighbour lists
//! You can only use this for neighbour lists with one type
std::vector<std::vector<int>> neighListO(
    double rcutoff, molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    int typeI);

//! Inefficient @f$O(n^2)@f$ implementation of neighbour lists
//! You can only use this for neighbour lists with one type
std::vector<std::vector<int>> halfNeighList(
    double rcutoff, molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    int typeI = 1);

//! The following function outputs a neighbour list using indices and NOT atom
//! IDs

//! Converts the neighbour list build with atom IDs into a neighbour list of
//! atom indices, according to the pointCloud
std::vector<std::vector<int>> neighbourListByIndex(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList);

//! Gets a neighbour list by index, according to a pointCloud given as the
//! input. Assume no slices or other skullduggery
std::vector<std::vector<int>> getNewNeighbourListByIndex(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, double cutoff);

//! Erases memory for a vector of vectors for the neighbour list
int clearNeighbourList(std::vector<std::vector<int>> &nList);

}  // namespace nneigh

#endif  // __NEIGHBOURS_H_
