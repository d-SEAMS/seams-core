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
#include <cage.hpp>

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
    molSys::PointCloud<molSys::Point<double>, double> yCloud, 
    std::string filename, int targetFrame,
    int atomTypeI, bool isSlice, std::array<double, 3> coordLow, std::array<double, 3> coordHigh,
    std::string templateFileO, int oxygenAtomType, double rcutoff=3.5);

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

  //! Generate a vector of vectors containing ordered 6-membered rings and a last ring with
  //! 4 water molecules corresponding to the 4 water molecules not present in the 6-membered rings 
  std::vector<std::vector<int>>
  getShapeMatchedRingsInCage(molSys::PointCloud<molSys::Point<double>, double> targetCloud, molSys::PointCloud<molSys::Point<double>, double> refCloud, 
  std::vector<std::vector<int>> ringsRef, double rcutoff, int oxygenAtomType);

  //! Update the atom types of Clathrate S2 cage atoms 
  void assignAtomTypes(std::vector<int> atomIndices,
                               std::vector<cage::iceType> *atomTypes,
                               cage::iceType clathType = cage::clathS2);

} // namespace clath

namespace misc {
  
  //! Function for getting the centroid of molecules given a particular atom type 
  std::vector<std::vector<double>> 
  getCentroidMolecules(std::string filename, int targetFrame,
    int atomTypeI, bool isSlice = false, std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
                std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

  //! Function for finding the k closest points (of type atomType) in a pointCloud, from a given target point (x y z coordinates).
  //! Returns a vector of k atom indices in yCloud corresponding to the k closest points 
  std::pair<std::vector<int>, molSys::PointCloud<molSys::Point<double>, double>>
  kClosestPoints(molSys::PointCloud<molSys::Point<double>, double> yCloud, int atomType, 
    std::vector<double> targetPointCoord, int k, double maxCutoff);

  //! Function for reordering a vector of atomIndices given a
  //! vector of vectors corresponding to the rings 
  std::vector<int> reorderAtomIndices(std::vector<int> atomIndices, 
    std::vector<std::vector<int>> rings);

  //! Function for updating the RMSD per atom for every O atom in a
  //! clathrate cage. 
  void updateRMSDatom(std::vector<int> atomIndices, 
  std::vector<std::vector<int>> rings,
  std::vector<double> rmsdList, std::vector<double> *rmsdPerAtom,
  std::vector<int> *noOfCommonAtoms);

} // namespace misc 

#endif // __BULKCLATHRATE_H_
