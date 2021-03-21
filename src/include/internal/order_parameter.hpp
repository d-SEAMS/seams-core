//-----------------------------------------------------------------------------------
// d-SEAMS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
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
