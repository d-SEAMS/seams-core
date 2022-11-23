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
  void shapeMatchS2ClathrateSystem(std::string path,
    std::vector<std::vector<int>> nList, molSys::PointCloud<molSys::Point<double>, double> yCloud, 
    std::string filename, int targetFrame,
    int atomTypeI, bool isSlice, std::array<double, 3> coordLow, std::array<double, 3> coordHigh,
    std::string templateFileO, int oxygenAtomType);

  //! Build a reference SII large cage (5^12 6^4) reading in from a template trajectory file
  std::pair<Eigen::MatrixXdRowMajor, Eigen::MatrixXdRowMajor>
  buildRefS2CageLammpsTrj(std::string filename, std::string filenameO, int oxygenAtomType);

  //! Build a reference SII large cage (5^12 6^4) reading in from a template trajectory file
  std::tuple<molSys::PointCloud<molSys::Point<double>, double>,std::vector<std::vector<int>> , Eigen::MatrixXdRowMajor>
  buildRefS2Cage(std::string filename, int oxygenAtomType);

  //! Add a given vector to the current rings vector of vectors for a 
  //! target clathrate structure. Match the reference and target structures.
  void matchClathrateLastRing(std::vector<std::vector<int>> targetRings, std::vector<int> lastTargetVec,
  std::vector<std::vector<int>> refRings, molSys::PointCloud<molSys::Point<double>, double> targetCloud,
  molSys::PointCloud<molSys::Point<double>, double> refCloud,
  std::vector<double> *quat, double *rmsd,
  std::vector<double> *rmsdList, double *scale);

} // namespace clath

namespace misc {
  
  //! Function for getting the COM of molecules given a particular atom type 
std::vector<std::vector<double>> 
getCentroidMolecules(std::string filename, int targetFrame,
  int atomTypeI, bool isSlice = false, std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
              std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

} // namespace misc 

#endif // __BULKCLATHRATE_H_
