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

#ifndef __TOPO_ONE_DIM_H_
#define __TOPO_ONE_DIM_H_

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
#include <shapeMatch.hpp>

/** @file topo_one_dim.hpp
 *  @brief File containing functions used specific to quasi-one-dimensional
 *   topological network critera (the prism identification scheme).
 */

/**
 *  @addtogroup ring
 *  @{
 */

namespace ring {

//! Find out which rings are prisms.
//! Returns a vector containing all the ring IDs which are prisms
std::vector<int> findPrisms(
    std::vector<std::vector<int>> rings, std::vector<strucType> *ringType,
    int *nPerfectPrisms, int *nImperfectPrisms,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<double> *rmsdPerAtom, bool doShapeMatching = false);

//! Tests whether two rings are basal rings (true) or not (false) for a prism
//! (strict criterion)
bool basalPrismConditions(std::vector<std::vector<int>> nList,
                          std::vector<int> *basal1, std::vector<int> *basal2);

//! Reduced criterion: Two candidate basal rings of a prism block should have at
//! least one bond between them
bool relaxedPrismConditions(std::vector<std::vector<int>> nList,
                            std::vector<int> *basal1, std::vector<int> *basal2);

//! Checks whether two 4-membered rings are parallel in one dimension or not to
//! prevent overcounting
bool discardExtraTetragonBlocks(
    std::vector<int> *basal1, std::vector<int> *basal2,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Saves only axial rings out of all possible rings
std::vector<std::vector<int>> keepAxialRingsOnly(
    std::vector<std::vector<int>> rings,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Find out which rings are prisms, looping through all ring sizes upto the
//! maxDepth The input ringsAllSizes array has rings of every size.
int prismAnalysis(std::string path, std::vector<std::vector<int>> rings,
                  std::vector<std::vector<int>> nList,
                  molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                  int maxDepth, int *atomID, int firstFrame, int currentFrame,
                  bool doShapeMatching = false);

//! Assign an atomType (equal to the number of nodes in the ring)
//! given a vector with a list of indices of rings comprising the prisms
int assignPrismType(std::vector<std::vector<int>> rings,
                    std::vector<int> listPrism, int ringSize,
                    std::vector<ring::strucType> ringType,
                    std::vector<int> *atomTypes,
                    std::vector<ring::strucType> *atomState);

//! Get the atom type values for deformed prisms
int deformedPrismTypes(std::vector<ring::strucType> atomState,
                       std::vector<int> *atomTypes, int maxDepth);

//! Shift the entire ice nanotube and remove axial translations
int rmAxialTranslations(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int *atomID,
    int firstFrame, int currentFrame);

}  // namespace ring

#endif  // __TOPOCONFINED_H_
