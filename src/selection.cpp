//-----------------------------------------------------------------------------------
// d-SEAMS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <selection.hpp>

// -----------------------------------------------------------------------------------------------------
// FUNCTIONS FOR SELECTIONS
// -----------------------------------------------------------------------------------------------------

/**
 * @details Function that loops through a given input pointCloud and 
 * returns a new pointCloud only containing atoms of a given atom type ID.  
 * This is registered as a Lua function, and is exposed to the user directly. 
 * The function calls the following functions internally:
 * - ring::clearRingList (Clears the vector of vectors for rings of a single
 * type, to prevent excessive memory being blocked).
 * - ring::getSingleRingSize (Fill a vector of vectors for rings of a particular
 * ring size).
 * - ring::findPrisms (Now that rings of a particular size have been obtained,
 * prism blocks are found and saved).
 * - topoparam::normHeightPercent (Gets the height% for the prism blocks).
 * - ring::assignPrismType (Assigns a prism type to each atom type).
 * - sout::writePrismNum (Write out the prism information for the current
 * frame).
 * - sout::writeLAMMPSdataAllPrisms (Writes out a LAMMPS data file for the
 * current frame, which can be visualized in OVITO).
 * @param[in] yCloud The given input PointCloud
 * @param[in] atomTypeI The type ID of the atoms to save in the output PointCloud
 * @param[in] isSlice This decides whether a slice will be used or not
 * @param[in] coordLow Contains the lower limits of the slice, if a slice is to
 *  be created
 * @param[in] coordHigh Contains the upper limits of the slice, if a slice is
 *  to be created
 */
molSys::PointCloud<molSys::Point<double>, double>
gen::getPointCloudOneAtomType(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    int atomTypeI, bool isSlice, std::array<double, 3> coordLow,
    std::array<double, 3> coordHigh) {
  //
  molSys::PointCloud<molSys::Point<double>, double> outCloud; // Output PointCloud
  int natoms = 0;             // Number of atoms of the desired type 
  bool pointInSlice = true; // If the current point is in the slice, this is true (otherwise this is false)

  // Loop through every iatom and check the type
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    // Skip if the atom is not of the desired type 
    if (yCloud->pts[iatom].type != atomTypeI) {
      continue;
    } // end of checking for the type 
    // ----------------
    // only if a slice has been requested
    if (isSlice) {
      pointInSlice = sinp::atomInSlice(yCloud->pts[iatom].x, yCloud->pts[iatom].y, yCloud->pts[iatom].z,
                                               coordLow, coordHigh);
      // Skip if the atom is not part of the slice
      if (!pointInSlice) {
        continue;
      } // skipped for atoms not in the slice 
    } // end of slice handling 
    //
    // Actually add the atoms to the outCloud
    natoms++; // Update the number of atoms in outCloud
    outCloud.pts.push_back(yCloud->pts[iatom]); // Update the pts vector
    outCloud.idIndexMap[yCloud->pts[iatom].atomID] =
                outCloud.pts.size() - 1; // array index
  }

  // Update the number of particles in the PointCloud
  outCloud.nop = outCloud.pts.size();

  // Update the frame number
  outCloud.currentFrame = yCloud->currentFrame; 

  // nop should be the same as natoms

  return outCloud;
}