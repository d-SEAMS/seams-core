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

#ifndef __TOPO_TWO_DIM_H_
#define __TOPO_TWO_DIM_H_

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

#include <mol_sys.hpp>
#include <order_parameter.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

/** @file topo_one_dim.hpp
 *  @brief File containing functions used specific to quasi-one-dimensional
 *   topological network critera (the prism identification scheme).
 */

/**
 *  @addtogroup ring
 *  @{
 */

namespace ring {

//! Find out which rings are prisms, looping through all ring sizes upto the
//! maxDepth The input ringsAllSizes array has rings of every size.
int polygonRingAnalysis(
    std::string path, std::vector<std::vector<int>> rings,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int maxDepth,
    double sheetArea, int firstFrame);

} // namespace ring

#endif // __TOPOCONFINED_H_
