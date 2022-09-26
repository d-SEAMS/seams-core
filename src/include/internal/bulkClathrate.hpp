//-----------------------------------------------------------------------------------
// d-SEAMS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef __BULKCLATHRATE_H_
#define __BULKCLATHRATE_H_

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
#include <utility> 

#include <franzblau.hpp>
#include <neighbours.hpp>
#include <pntCorrespondence.hpp>
#include <topo_bulk.hpp>

/** @file bulkClathrate.hpp
 *   @brief File containing functions used specific to criteria for clathrates
 (specifically SII clathrates for now)
 */

/*!
 *  @addtogroup clath
 *  @{
 */

namespace clath {

  //! Shape-matching algorithms for S2 clathrates
  void shapeMatchS2ClathrateSystem(std::string templateFileName, std::string templateFileO, int oxygenAtomType);

  //! Build a reference SII large cage (5^12 6^4) reading in from a template trajectory file
  std::pair<Eigen::MatrixXdRowMajor, Eigen::MatrixXdRowMajor>
  buildRefS2CageLammpsTrj(std::string filename, std::string filenameO, int oxygenAtomType);

} // namespace clath

#endif // __BULKCLATHRATE_H_
