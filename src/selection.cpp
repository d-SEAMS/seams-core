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
 * sets the inSlice bool for every Point according to whether the atom  
 * is in the specified (single) slice or not. Not inclusive of atoms in molecules 
 * @param[in] yCloud The given input PointCloud
 * @param[in] clearPreviousSliceSelection sets all inSlice bool values to false before 
 * adding Points to the slice
 * @param[in] coordLow Contains the lower limits of the slice, if a slice is to
 *  be created
 * @param[in] coordHigh Contains the upper limits of the slice, if a slice is
 *  to be created
 */
void gen::atomsInSingleSlice(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    bool clearPreviousSliceSelection,
    std::array<double, 3> coordLow,
    std::array<double, 3> coordHigh) {
  //
  bool pointInSlice = true; // If the current point is in the slice, this is true (otherwise this is false)

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

    yCloud->pts[iatom].inSlice = pointInSlice; // iatom is inside the slice
    //
  }

  return;
}

/**
 * @details Function that loops through a given input PointCloud and 
 * sets the inSlice bool for every Point according to whether the molecule  
 * is in the specified (single) slice or not. If even one atom of a molecule 
 * is inside the region, then all atoms belonging to that molecule should be
 * inside the slice as well (therefore, inSlice would be set to true)
 * NOTE: THIS DOES NOT WORK. ERROR
 * @param[in] yCloud The given input PointCloud
 * @param[in] clearPreviousSliceSelection sets all inSlice bool values to false before 
 * adding Points to the slice
 * @param[in] coordLow Contains the lower limits of the slice, if a slice is to
 *  be created
 * @param[in] coordHigh Contains the upper limits of the slice, if a slice is
 *  to be created
 */
void gen::moleculesInSingleSlice(
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

  return;
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
void gen::setAtomsWithSameMolID(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::unordered_multimap<int,int> molIDAtomIDmap,
    int molID, bool inSliceValue) {
  //
  int jatomID;    // atom ID of the current jatom
  int jatomIndex; // index of jatom

  // Find all atoms with molID and set the inSlice bool to inSliceValue
  auto range = molIDAtomIDmap.equal_range(molID); 

  // Loop through all atoms with molID
  for (auto it = range.first; it != range.second; it++){
    // it->second gives the value (in this case, the atom ID)
    jatomID = it->second; // Atom ID with molecule ID equal to iatomMolID
    auto gotJ = yCloud->idIndexMap.find(jatomID);
    jatomIndex = gotJ->second;
    // Set the jatom inSlice bool to true
    yCloud->pts[jatomIndex].inSlice = inSliceValue; // jatomIndex is assigned inSliceValue
  } // end of loop through all atoms with molID

  // IC(yCloud->pts[jatomIndex]);

  return;
}

/**
 * @details Function that loops through the PointCloud used to construct the
 * neighbour list (used to generate primitive rings) and sets the inSlice bool values
 * of edge atoms which belong to rings that are formed by atoms in the slice. 
 * The output PointCloud may not be the same as the PointCloud used to 
 * construct the nList, and the inSlice bool value can be set for both.
 * @param[in] rings Vector of vectors of the primitive rings (by index) according to oCloud
 * @param[in] oCloud PointCloud of O atoms, used to construct the rings vector of vectors
 * @param[in] yCloud The output PointCloud (may contain more than just the O atoms)
 * @param[in] identicalCloud bool value; if this is true then oCloud and yCloud are the same
 * @param[in] coordLow Contains the lower limits of the slice (needed to find all the molecules of oCloud
 * in the slice)
 * @param[in] coordHigh Contains the upper limits of the slice
 */
void ring::getEdgeMoleculesInRings(
    std::vector<std::vector<int>> rings, molSys::PointCloud<molSys::Point<double>, double> *oCloud,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::array<double, 3> coordLow, std::array<double, 3> coordHigh, bool identicalCloud) {
  //
  // A vector of bool values, such that every ring has a value of true (in the slice) or false (not in the slice) 
  std::vector<bool> ringInSlice(rings.size(), false); // all set to false initially.
  int jatomIndex, jatomID; // Index and ID in oCloud 
  int jatomIndex1; // Index in yCloud  
  std::unordered_multimap<int, int>
      molIDAtomIDmap; // Unordered multimap with molecule IDs of the atoms as the keys and the
                  // atom IDs as the values. More than one atom can have the same molecule ID

  // --------------------
  // Sets all O atoms within the slice to an inSlice bool value of true. If a single atom of a molecule is in the
  // slice, then all atoms in the molecule will also be inside the slice 
  gen::moleculesInSingleSlice(oCloud, true, coordLow, coordHigh);
  // --------------------
  // Get the unordered map of the atom IDs (keys) and the molecular IDs
  // (values)
  molIDAtomIDmap = molSys::createMolIDAtomIDMultiMap(yCloud);
  // --------------------

  // Loop through every iatom present in the slice in oCloud
  // Check to see if iatom is in any ring. If it is present in a ring, set that ring to true 
  // Change the inSlice bool of every iatom in the ring to true if false. Set also the inSlice bool
  // of the output cloud if it is not identical. 

  // Loop through every iatom present in the slice in oCloud to determine which rings are in the slice
  // The indices of oCloud and rings should match (or this will produce unexpected results)
  for (int iatom = 0; iatom < oCloud->nop; iatom++) {
    // Skip if iatom is not in the slice
    if (!oCloud->pts[iatom].inSlice)
    {
      continue;
    } // skip for iatom not in slice
    //
    // For iatom in the slice, 
    // loop through all rings to find all the rings it is a part of 
    for (int iring = 0; iring < rings.size(); iring++)
    {
      // Skip if iring is in the slice already
      if (ringInSlice[iring])
      {
        continue;
      } // skip for iring in slice 
      // Check and see if iatom is in iring 
      if(std::find(rings[iring].begin(), rings[iring].end(), iatom)!=rings[iring].end()){
        // Found iatom; ring is part of the slice 
        ringInSlice[iring] = true; // update the vector of bool values
        // --------------------------
        // Change the inSlice bool of every iatom in oCloud 
        // (and optionally, yCloud) in the ring to true
        // Loop through the elements of the ring 
        for (int j = 0; j < rings[iring].size(); j++)
        {
          jatomIndex = rings[iring][j]; // Index of the atom in oCloud 
          // Set this to true in oCloud 
          oCloud->pts[jatomIndex].inSlice = true; // part of slice 
          jatomID = oCloud->pts[jatomIndex].atomID; // Atom ID 
          // Now if oCloud and yCloud are not the same, use the
          // atom ID to set the inSlice bool value in yCloud 
          if (!identicalCloud)
          {
            // Find the index corresponding to the same atom in yCloud 
            auto gotJ = yCloud->idIndexMap.find(jatomID);
            jatomIndex1 = gotJ->second;
            // throw if not found ?
            // Set the jatom inSlice bool to true
            yCloud->pts[jatomIndex1].inSlice = true; // jatomIndex is inside the slice 
            // set the inSlice value of all atoms in yCloud with the current molecule ID 
            gen::setAtomsWithSameMolID(yCloud, molIDAtomIDmap, 
              yCloud->pts[jatomIndex1].molID, true);
          } // end of setting values in yCloud 
        } // end of loop through the elements of the current ring 
        // --------------------------
      } // found iatom in the ring
      //
    } // end of loop through all rings searching for iatom 
    //
  } // end of loop through atoms in oCloud in the slice 


  // Optionally add support for molecule slice update if yCloud and oCloud are not identical?
  // Some other bool?  
  // IC(ringInSlice);

  return;
}

// -----------------------------------------------------------------------------------------------------
// FUNCTIONS FOR SELECTIONS (RINGS)
// -----------------------------------------------------------------------------------------------------

/**
 * @details Function that loops through the PointCloud used to construct the
 * neighbour list (used to generate primitive rings) and sets the inSlice bool values
 * of edge atoms which belong to rings that are formed by atoms in the slice. 
 * The selected molecule IDs are printed to a separate file, and the molecules and 
 * atoms in the slice will be output to a LAMMPS data file (with the original box volume)
 * The output PointCloud may not be the same as the PointCloud used to 
 * construct the nList, and the inSlice bool value can be set for both.
 * @param[in] rings Vector of vectors of the primitive rings (by index) according to oCloud
 * @param[in] oCloud PointCloud of O atoms, used to construct the rings vector of vectors
 * @param[in] yCloud The output PointCloud (may contain more than just the O atoms)
 * @param[in] identicalCloud bool value; if this is true then oCloud and yCloud are the same
 * @param[in] coordLow Contains the lower limits of the slice (needed to find all the molecules of oCloud
 * in the slice)
 * @param[in] coordHigh Contains the upper limits of the slice
 */
void ring::printSliceGetEdgeMoleculesInRings(
    std::string path, std::vector<std::vector<int>> rings, 
    molSys::PointCloud<molSys::Point<double>, double> *oCloud,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::array<double, 3> coordLow, std::array<double, 3> coordHigh, bool identicalCloud) {
  //

  //Given the full yCloud PointCloud, set the inSlice bool for every atom,
  // if the molecules are inside the specified (single) region. 
  // gen::atomsInSingleSlice(yCloud, true, coordLow, coordHigh);
  gen::moleculesInSingleSlice(yCloud, true, coordLow, coordHigh);

  // Make sure that molecules which participate in the rings inside the slice are also 
  // in the selection 
  ring::getEdgeMoleculesInRings(rings, oCloud, yCloud, coordLow, coordHigh, identicalCloud);

  // Print out the molecule IDs of all the atoms in the slice
  sout::writeMoleculeIDsInSlice(path, yCloud);
  // Print out a command that could be used for an expression select command in OVITO
  sout::writeMoleculeIDsExpressionSelectOVITO(path, yCloud);

  // Print out the dump of all atoms and molecules, with an inSlice value printed in a separate column
  // H atoms not included in the slice (TODO: fix)
  sout::writeLAMMPSdumpSlice(yCloud, path); 
  
  IC(rings[0]);

  return;
}
