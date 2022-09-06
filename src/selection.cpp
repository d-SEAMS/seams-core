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
#include <icecream.hpp>

// -----------------------------------------------------------------------------------------------------
// FUNCTIONS FOR SELECTIONS
// -----------------------------------------------------------------------------------------------------

/**
 * @details Function that loops through a given input pointCloud and 
 * returns a new pointCloud only containing atoms of a given atom type ID.  
 * This is registered as a Lua function, and is exposed to the user directly. 
 * @param[in] yCloud The given input PointCloud
 * @param[out] outCloud The output PointCloud
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
    molSys::PointCloud<molSys::Point<double>, double> *outCloud,
    int atomTypeI, bool isSlice, std::array<double, 3> coordLow,
    std::array<double, 3> coordHigh) {
  //
  int natoms = 0;             // Number of atoms of the desired type 
  bool pointInSlice = true; // If the current point is in the slice, this is true (otherwise this is false)

  // --------
  // Before filling up the PointCloud, if the vectors are filled
  // empty them
  *outCloud = molSys::clearPointCloud(outCloud);
  // --------

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
    outCloud->pts.push_back(yCloud->pts[iatom]); // Update the pts vector
    outCloud->idIndexMap[yCloud->pts[iatom].atomID] =
                outCloud->pts.size() - 1; // array index
  }

  // Update the number of particles in the PointCloud
  outCloud->nop = outCloud->pts.size();

  // Box and box lengths 
  outCloud->box = yCloud->box;
  outCloud->boxLow = yCloud->boxLow;

  // Update the frame number
  outCloud->currentFrame = yCloud->currentFrame; 

  return *outCloud;
}

/**
 * @details Function that loops through a given input PointCloud and 
 * sets the inSlice bool for every Point according to whether the molecule  
 * is in the specified (single) slice or not. If even one atom of a molecule 
 * is inside the region, then all atoms belonging to that molecule should be
 * inside the slice as well (therefore, inSlice would be set to true)
 * @param[in] yCloud The given input PointCloud
 * @param[in] clearPreviousSliceSelection sets all inSlice bool values to false before 
 * adding Points to the slice
 * @param[in] coordLow Contains the lower limits of the slice, if a slice is to
 *  be created
 * @param[in] coordHigh Contains the upper limits of the slice, if a slice is
 *  to be created
 */
molSys::PointCloud<molSys::Point<double>, double>
gen::moleculesInSingleSlice(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    bool clearPreviousSliceSelection,
    std::array<double, 3> coordLow,
    std::array<double, 3> coordHigh) {
  //
  bool pointInSlice = true; // If the current point is in the slice, this is true (otherwise this is false)
  std::unordered_multimap<int, int>
      molIDAtomIDmap; // Unordered multimap with molecule IDs of the atoms as the keys and the
                  // atom IDs as the values. More than one atom can have the same molecule ID
  int iatomMolID; // molID of the current iatom
  int jatomID;    // atom ID of the current jatom
  int jatomIndex; // index of jatom

  // --------------------
  // Get the unordered map of the atom IDs (keys) and the molecular IDs
  // (values)
  molIDAtomIDmap = molSys::createMolIDAtomIDMultiMap(yCloud);
  // --------------------
  // Set inSlice to false for every atom first
  if (clearPreviousSliceSelection)
  {
    for (int iatom = 0; iatom < yCloud->nop; iatom++) {
      yCloud->pts[iatom].inSlice = false;
    }
  } // end of init
  // --------------------

  // Loop through every iatom and check if the atom is in the slice or not 
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    // Check if the current atom is in the slice 
    pointInSlice = sinp::atomInSlice(yCloud->pts[iatom].x, yCloud->pts[iatom].y, yCloud->pts[iatom].z,
                                               coordLow, coordHigh);
    //
    // Set the inSlice bool for the particular Point 
    if (pointInSlice)
    {
      yCloud->pts[iatom].inSlice = true; // iatom is inside the slice
      // Find mol ID and atoms with that particular molecular ID 
      iatomMolID = yCloud->pts[iatom].molID; // molecule ID 
      // -----------
      // Find all atoms with iatomMolID and set the inSlice bool
      // to true
      auto range = molIDAtomIDmap.equal_range(iatomMolID); 
      // Loop through all atoms with iatomMolID
      for (auto it = range.first; it != range.second; it++)
      {
        // it->second gives the value (in this case, the atom ID)
        jatomID = it->second; // Atom ID with molecule ID equal to iatomMolID
        auto gotJ = yCloud->idIndexMap.find(jatomID);
        jatomIndex = gotJ->second;
        // Set the jatom inSlice bool to true
        yCloud->pts[jatomIndex].inSlice = true; // jatomIndex is inside the slice 
      }
      // -----------
    } // the atom is in the slice 
    else{
      yCloud->pts[iatom].inSlice = false; // iatom is not in the slice  
    } // atom is not in the slice
  }

  return *yCloud;
}