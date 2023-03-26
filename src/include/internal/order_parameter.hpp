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

#ifndef __ORDER_PARAMETER_H_
#define __ORDER_PARAMETER_H_

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

#include <cage.hpp>
#include <mol_sys.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

namespace topoparam {

//! Calculates the height%, an average measure of filled volume. The average
//! height of a prism can be taken to be 2.75-2.85 Angstrom. (Koga et. al.,
//! 2001)
double
normHeightPercent(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                  int nPrisms, double avgPrismHeight);

//! Calculates the coverage area%, an area-based measure of relative proportion
//! of monolayer ices.
std::vector<double>
calcCoverageArea(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 std::vector<std::vector<int>> rings, double sheetArea);

//! Calculates the projected area on the XY, YZ and XZ planes
std::vector<double>
projAreaSingleRing(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                   std::vector<int> ring);

} // namespace topoparam

#endif // __ORDER_PARAMETER_H_
