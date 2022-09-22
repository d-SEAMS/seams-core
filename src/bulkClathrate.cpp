//-----------------------------------------------------------------------------------
// d-SEAMS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <bulkClathrate.hpp>

// -----------------------------------------------------------------------------------------------------
// ALGORITHMS FOR CLATHRATES 
// -----------------------------------------------------------------------------------------------------

/**
 * @details Build a reference SII cage, consisting of 28 water molecules, reading it in from a template
 * file saved in the templates directory 
 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> 
clath::buildRefS2CageLammpsTrj(std::string filename, std::string filenameO, int oxygenAtomType) {
  //
  Eigen::MatrixXd refPntsO(28, 3); // Reference point set of just O atoms (Eigen matrix)
  Eigen::MatrixXd refPntsWat(84, 3); // Reference point set of O H H water atoms (Eigen matrix)
  // Get the reference HC point set
  molSys::PointCloud<molSys::Point<double>, double>
      waterCloud; // PointCloud for holding the reference point values for all the water molecules
      molSys::PointCloud<molSys::Point<double>, double>
      oCloud; // PointCloud for holding the reference point values for just the O atoms
  //
  // read in all the water molecules into the pointCloud setCloud
  //
  // All water molecules 
  waterCloud = sinp::readLammpsTrj(filename, 1, &waterCloud);
  oCloud = sinp::readLammpsTrjO(filenameO, 1, &oCloud, oxygenAtomType);

  return std::make_pair (refPntsO, refPntsWat);
}
