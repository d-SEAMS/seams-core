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

#include <topo_bulk.hpp>

// -----------------------------------------------------------------------------------------------------
// BULK RING SEARCH ONLY
// -----------------------------------------------------------------------------------------------------

/**
 * @details Function that loops through the primitive rings (which is a vector
 * of vectors) of all sizes, upto maxDepth (the largest ring size). 
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
 * @param[in] path The string to the output directory, in which files will be
 *  written out.
 * @param[in] rings Row-ordered vector of vectors for rings of a single type.
 * @param[in] nList Row-ordered neighbour list by index.
 * @param[in] yCloud The input PointCloud.
 * @param[in] maxDepth The maximum possible size of the primitive rings.
 * @param[in] firstFrame The first frame to be analyzed
 */
int ring::bulkPolygonRingAnalysis(
    std::string path, std::vector<std::vector<int>> rings,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int maxDepth,
    int firstFrame) {
  //
  std::vector<std::vector<int>>
      ringsOneType;           // Vector of vectors of rings of a single size
  int nRings;                 // Number of rings of the current type
  std::vector<int> nRingList; // Vector of the values of the number of rings
                              // for a particular frame
  std::vector<int>
      atomTypes; // contains int values for each ring type considered
  // -------------------------------------------------------------------------------
  // Init
  nRingList.resize(
      maxDepth -
      2); // Has a value for every value of ringSize from 3, upto maxDepth
  // The atomTypes vector is the same size as the pointCloud atoms
  atomTypes.resize(yCloud->nop, 1); // The dummy or unclassified value is 1
  // -------------------------------------------------------------------------------
  // Run this loop for rings of sizes upto maxDepth
  // The smallest possible ring is of size 3
  for (int ringSize = 3; ringSize <= maxDepth; ringSize++) {
    // Clear ringsOneType
    ring::clearRingList(ringsOneType);
    // Get rings of the current ring size
    ringsOneType = ring::getSingleRingSize(rings, ringSize);
    //
    // Continue if there are zero rings of ringSize
    if (ringsOneType.size() == 0) {
      nRingList[ringSize - 3] = 0; // Update the number of prisms
      continue;
    } // skip if there are no rings
    //
    // -------------
    // Number of rings with n nodes
    nRings = ringsOneType.size();
    // -------------
    // Now that you have rings of a certain size:
    nRingList[ringSize - 3] = nRings; // Update the number of n-membered rings
    // -------------
  } // end of loop through every possible ringSize

  // Get the atom types for all the ring types
  ring::assignPolygonType(rings, &atomTypes, nRingList);

  // Write out the ring information
  sout::writeRingNumBulk(path, yCloud->currentFrame, nRingList, maxDepth, firstFrame);
  // Write out the lammps data file for the particular frame
  sout::writeLAMMPSdataAllRings(yCloud, nList, atomTypes, maxDepth, path, false);

  return 0;
}

// -----------------------------------------------------------------------------------------------------
// DDC / HC ALGORITHMS
// -----------------------------------------------------------------------------------------------------

/**
 * @details Finds out if rings constitute double-diamond cages or hexagonal
 * cages. Requires a neighbour list (by index) and a vector of vectors of
 * primitive rings which should also be by index. This is registered as a Lua
 * function and is accessible to the user. Internally, this function calls the
 * following functions:
 * - ring::getSingleRingSize (Saves rings of a single ring size into a new
 * vector of vectors, which is subsequently used for finding DDCs, HCs etc).
 * - ring::findDDC (Finds the DDCs).
 * - ring::findHC (Finds the HCs).
 * - ring::findMixedRings (Finds the mixed rings, which are shared by DDCs and
 * HCs both).
 * - ring::getStrucNumbers (Gets the number of structures (DDCs, HCs, mixed
 * rings, basal rings, prismatic rings, to be used for write-outs).
 * - sout::writeTopoBulkData (Writes out the numbers and data obtained for the
 * current frame).
 * - ring::getAtomTypesTopoBulk (Gets the atom type for every atom, to be used
 * for printing out the ice types found).
 * - sout::writeLAMMPSdataTopoBulk (Writes out the atoms, with the classified
 * types, into a LAMMPS data file, which can be visualized in OVITO).
 *  @param[in] path The file path of the output directory to which output files
 * will be written.
 *  @param[in] rings Vector of vectors containing the primitive rings. This
 * contains rings of all sizes.
 *  @param[in] nList Row-ordered neighbour list, by index.
 *  @param[in] yCloud The input PointCloud, with respect to which the indices in
 * the rings and nList vector of vectors have been saved.
 *  @param[in] firstFrame First frame to be analyzed
 *  @param[in] onlyTetrahedral Flag for only finding DDCs and HCs (true) or also
 * finding PNCs (false)
 */
int ring::topoBulkAnalysis(
    std::string path, std::vector<std::vector<int>> rings,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int firstFrame,
    bool onlyTetrahedral) {
  //
  // Ring IDs of each type will be saved in these vectors
  std::vector<int> listDDC; // Vector for ring indices of DDC
  std::vector<int> listHC;  // Vector for ring indices of HC
  std::vector<int>
      listMixed; // Vector for ring indices of rings that are both DDC and HC
  std::vector<ring::strucType>
      ringType; // This vector will have a value for each ring inside
                // ringList, of type enum strucType in gen.hpp
  // Make a list of all the DDCs and HCs
  std::vector<cage::Cage> cageList;
  std::vector<std::vector<int>>
      ringsOneType;    // Vector of vectors of rings of a single size
  int initRingSize;    // Todo or not: calculate the PNCs or not
  int maxRingSize = 6; // DDCs and HCs are for 6-membered rings
  std::vector<cage::iceType>
      atomTypes; // This vector will have a value for every atom
  // Number of types
  int numHC, numDDC, mixedRings, prismaticRings, basalRings;

  // TODO: handle shapeMatching
  bool doShapeMatching = true;
  // Qualifier for the RMSD per atom:
  std::vector<double> rmsdPerAtom;
  //

  // test
  // std::cout << "rings"
  //           << "\n";
  // for (int i = 0; i < rings.size(); i++) {
  //   for (int j = 0; j < rings[i].size(); j++) {
  //     std::cout << rings[i][j] << " ";
  //   }
  //   std::cout << "\n";
  // }
  // std::cout << "\n";

  if (onlyTetrahedral) {
    initRingSize = 6;
  } else {
    initRingSize = 5;
  }

  // Init the atom type vector
  atomTypes.resize(yCloud->nop); // Has a value for each atom
  // Init the rmsd per atom (not used yet)
  rmsdPerAtom.resize(yCloud->nop); // Has a value for each atom

  // ----------------------------------------------
  // Init
  // Get rings of size 5 or 6.
  for (int ringSize = initRingSize; ringSize <= maxRingSize; ringSize++) {
    // Clear ringsOneType
    ring::clearRingList(ringsOneType);
    // Get rings of the current ring size
    ringsOneType = ring::getSingleRingSize(rings, ringSize);
    // Skip for zero rings
    if (ringsOneType.size() == 0) {
      continue;
    } // skip for no rings of ringSize
    //
    // Init the ringType vector
    ringType.resize(
        ringsOneType.size()); // Has a value for each ring. init to zero.
    // ----------------------------------------------
    if (ringSize == 6) {
      // Get the cages

      // Find HC rings, saving the ring IDs (starting from 0) to listHC
      listHC = ring::findHC(ringsOneType, &ringType, nList, &cageList);

      // Find DDC rings, saving the IDs to listDDC
      listDDC = ring::findDDC(ringsOneType, &ringType, listHC, &cageList);

      // Find rings which are both DDCs and HCs (mixed)
      // A dummy value of -10 in the listDDC and listHC vectors for mixed rings
      listMixed =
          ring::findMixedRings(ringsOneType, &ringType, &listDDC, &listHC);

      // Get the number of structures (DDCs, HCs, mixed rings, basal rings,
      // prismatic rings)
      ring::getStrucNumbers(ringType, cageList, &numHC, &numDDC, &mixedRings,
                            &prismaticRings, &basalRings);

      // Write out to a file
      sout::writeTopoBulkData(path, yCloud->currentFrame, numHC, numDDC,
                              mixedRings, basalRings, prismaticRings,
                              firstFrame);
      // Gets the atom type for every atom, to be used for printing out the ice
      // types found
      ring::getAtomTypesTopoBulk(ringsOneType, ringType, &atomTypes);
    }
    // Finding prismatic blocks
    else {
      // Get the prism block classifications
      prism3::findBulkPrisms(ringsOneType, &ringType, nList, yCloud,
                             &rmsdPerAtom);
      // Gets the atom type for every atom, to be used for printing out the ice
      // types found
      ring::getAtomTypesTopoBulk(ringsOneType, ringType, &atomTypes);
    }
  }

  // Print out the lammps data file with the bonds
  sout::writeLAMMPSdataTopoBulk(yCloud, nList, atomTypes, path);
  // To output the bonds between dummy atoms, uncomment the following line
  // sout::writeLAMMPSdataTopoBulk(yCloud, nList, atomTypes, path, true);

  return 0;
}

/**
 * @details Determines which hexagonal rings are DDC rings. This function
 * returns a vector which contains the ring IDs of all the rings which are DDC
 * rings. The ring IDs correspond to the index of the rings inside the vector of
 * vector rings, starting from 0. DDC rings can be found using a three-step
 * procedure, in which equatorial rings and their corresponding rings can be
 * found. Peripheral rings can be shared, so first all the equatorial rings must
 * be found. Whenever an equatorial ring and the corresponding peripheral rings
 * are found, their values must be updated to DDC enum type inside the ringType
 * vector, which has been passed as an input to this function. At the end, the
 * output vector can be updated to avoid adding the same ring index more than
 * once (this may happen for peripheral rings, which can be shared). Internally,
 * the function calls the following:
 * - ring::conditionOneDDC (Step one: Finds all rings which contain each index
 * (m_k) of the equatorial ring, iring, in at least three other rings).
 * - ring::conditionTwoDDC (Step two: For every triplet in iring, there is at
 * least one hexagonal ring other than iring that passes through the triplet.
 * The peripheral rings are stored in order of the starting element of each
 * triplet.)
 * - ring::conditionThreeDDC (Step three: For every triplet in the equatorial
 * ring, there is at least one hexagonal ring other than iring that passes
 * through the triplet. Rings corresponding to triplets need not be searched
 * again since peripheralRings are stored in that order. Rings corresponding to
 * T1, T3, T5 must have a common element. Similarly rings corresponding to T2,
 * T4, T6 must have at least one common element. Alternating rings corresponding
 * to triplets must have at least three common elements).
 *
 * @param[in] rings Vector of vectors containing 6-membered primitive rings
 *  (wherein each ring contains the atom indices of the particles which
 * constitute the ring).
 * @param[in] ringType Vector containing a ring::strucType (structure type) for
 *  each ring.
 * @param[in] cageList Vector in which every cage is saved.
 * @return A vector of all the ring indices which constitute DDCs.
 */
std::vector<int> ring::findDDC(std::vector<std::vector<int>> rings,
                               std::vector<ring::strucType> *ringType,
                               std::vector<int> listHC,
                               std::vector<cage::Cage> *cageList) {
  std::vector<int> listDDC;
  int totalRingNum = rings.size();  // Total number of hexagonal rings
  std::vector<int> peripheralRings; // Indices which may be peripheral rings
  bool cond1, cond2, cond3;         // Conditions for DDC
  int jring;                        // Index for peripheral ring being added
  std::vector<int> DDCRings; // Indices of rings which constitute a single DDC,
                             // with the equatorial ring first
  // Vector for checking if a ring is basal, prismatic or peripheral
  std::vector<bool>
      notEquatorial; // true if the ring is prismatic, basal or peripheral.
  notEquatorial.resize(totalRingNum); // Initialized to false

  // --------
  // Set all basal and prismatic rings to true (which are part of an HC)
  for (int i = 0; i < listHC.size(); i++) {
    int currentRingIndex = listHC[i];
    // Set this to true
    notEquatorial[currentRingIndex] = true;
  } // end of update of notEquatorial
  // --------

  // To search for equatorial rings, loop through all
  // the hexagonal rings
  for (int iring = 0; iring < totalRingNum; iring++) {
    // ------------
    // Step zero: If the ring has been classified as a basal or prismatic ring
    // in an HC or is a peripheral ring, then it cannot be the equatiorial ring
    // in a DDC
    //
    if (notEquatorial[iring]) {
      continue;
    } // skip for rings which are not equatorial
    //
    // ------------
    // Init
    peripheralRings.clear();
    // ------------
    // Step one: Find all rings which contain each index (m_k) of the equatorial
    // ring, iring, in at least three other rings
    cond1 = ring::conditionOneDDC(rings, &peripheralRings, iring);
    if (cond1 == false) {
      continue;
    }
    // ------------
    // Step two: For every triplet in iring, there is at least one
    // hexagonal ring other than iring that passes through the triplet.
    // The peripheral rings are stored in order of the starting element
    // of each triplet.
    cond2 = ring::conditionTwoDDC(rings, &peripheralRings, iring);
    if (cond2 == false) {
      continue;
    }
    // ------------
    // Step three: For every triplet in the equatorial ring, there is at least
    // one hexagonal ring other than iring that passes through the triplet.
    // Rings corresponding to triplets need not be searched again since
    // peripheralRings are stored in that order. Rings corresponding to T1, T3,
    // T5 must have a common element. Similarly rings corresponding to T2, T4,
    // T6 must have at least one common element. Alternating rings corresponding
    // to triplets must have at least three common elements
    cond3 = ring::conditionThreeDDC(rings, &peripheralRings);
    if (cond3 == false) {
      continue;
    }
    // ------------
    // If the peripheral rings are duplicates, skip everything
    sort(peripheralRings.begin(), peripheralRings.end());
    peripheralRings.erase(
        unique(peripheralRings.begin(), peripheralRings.end()),
        peripheralRings.end());
    // There should be 6 unique peripheral rings
    if (peripheralRings.size() != 6) {
      continue;
    }
    // ------------
    // If iring is an equatorial ring, add it to the listDDC vector
    if ((*ringType)[iring] == ring::strucType::unclassified) {
      (*ringType)[iring] = ring::strucType::DDC;
      listDDC.push_back(iring);
    }
    // Add the peripheral ring IDs too
    for (int j = 0; j < peripheralRings.size(); j++) {
      jring = peripheralRings[j];
      if ((*ringType)[jring] == ring::strucType::unclassified) {
        (*ringType)[jring] = ring::strucType::DDC;
        listDDC.push_back(jring);
      } else if ((*ringType)[jring] == ring::strucType::HCbasal) {
        (*ringType)[jring] = ring::strucType::bothBasal;
        listDDC.push_back(jring);
      } // end of update
      // never true
      else if ((*ringType)[jring] == ring::strucType::HCprismatic) {
        (*ringType)[jring] = ring::strucType::bothPrismatic;
        listDDC.push_back(jring);
      }
      //
      // Update the notEquatorial vector
      notEquatorial[jring] = true;
    } // end of update for peripheral rings
    // Add rings to the cageList vector of struct Cages.
    DDCRings.clear();          // init
    DDCRings.push_back(iring); // Add the equatorial ring first
    DDCRings.insert(std::end(DDCRings), std::begin(peripheralRings),
                    std::end(peripheralRings));
    (*cageList).push_back({cage::cageType::DoubleDiaC, DDCRings});
    // ------------
  } // end of loop through all hexagonal rings

  return listDDC;
}

/**
 * @details For a given ring, which is being tested as the equatorial ring,
 * this function tests if each element of the ring \f$ (m_k) \f$ is present
 * in at least three other rings or not. Returns false if this is not true.
 * The ring indices of rings that have the common element are saved inside
 * the periperal ring vector (peripheralRings) as potential peripheral rings,
 * which is then passed as an input to the subsequent functions testing
 * conditions.
 * @param[in] rings Vector of vectors containing the 6-membered primitive
 *  rings.
 * @param[in] peripheralRings Vector containing the indices of rings which are
 *  potential peripheral rings.
 * @param[in] iring Index of the ring being tested as equatorial.
 * @return A bool; true if ring being tested as equatorial, iring, satisfies
 *  the above condition for being an equatorial ring, and false otherwise.
 */
bool ring::conditionOneDDC(std::vector<std::vector<int>> rings,
                           std::vector<int> *peripheralRings, int iring) {
  int index; // Atom ID to be compared
  int noOfCommonRings =
      0;        // No of rings in which the element to be matched has been found
  int jElement; // Atom ID being compared to index

  // Loop through each element of iring for finding matches
  for (int m = 0; m < 6; m++) {
    index = rings[iring][m]; // Atom Index to be compared and matched with
    noOfCommonRings = 0;     // init to zero.
    // Loop through every ring except iring
    for (int jring = 0; jring < rings.size(); jring++) {
      if (iring == jring) {
        continue;
      } // Skip for iring
      // -------
      // Search every element of jring
      for (int k = 0; k < 6; k++) {
        jElement = rings[jring][k];
        if (jElement == index) {
          noOfCommonRings++;
          (*peripheralRings).push_back(jring);
          break;
        } // if index is found inside jring
        else {
          continue;
        }
      } // end of loop through every element of jring
      // -------
    } // end of loop through all rings except iring
    // If less than 3 rings have been found for each element, then this is
    // not an equatorial ring
    if (noOfCommonRings < 3) {
      return false;
    } // end of check for common ring number per element in iring
  }   // end of loop through each element of iring

  // iring is an equatorial ring. The duplicate ring IDs inside
  // peripheralRings should be removed
  std::vector<int>::iterator ip; // Iterator to find the last element upto
                                 // which unique elements are present
  // Duplicate IDs must be removed
  int numElements =
      (*peripheralRings).size(); // number of elements in peripheralRings
  sort((*peripheralRings).begin(), (*peripheralRings).end());
  ip = std::unique((*peripheralRings).begin(),
                   (*peripheralRings).begin() + numElements);
  // Resize peripheral rings to remove undefined terms
  (*peripheralRings).resize(std::distance((*peripheralRings).begin(), ip));

  return true;
}

/**
 * @details For a given ring, which is being tested as the equatorial ring,
 * this function tests if each triplet that can be formed from the ring is
 * common to at least one other ring or not. Returns false if this is not true.
 * The ring indices of rings that have the common triplet are ultimately saved
 * inside the periperal ring vector as potential peripheral rings, which is
 * passed as an input to the subsequent condition-testing functions. The
 * function calls the following function internally:
 * - ring::findTripletInRing (Searches inside another ring with index jring for
 * the current triplet, which was obtained from iring).
 * @param[in] rings Vector of vectors containing the 6-membered primitive
 *  rings.
 * @param[in] peripheralRings Vector containing the indices of rings which are
 *  potential peripheral rings.
 * @param[in] iring Index of the ring being tested as equatorial.
 * @return A bool; true if ring being tested as equatorial, iring, satisfies
 *  the above condition for being an equatorial ring, and false otherwise.
 */
bool ring::conditionTwoDDC(std::vector<std::vector<int>> rings,
                           std::vector<int> *peripheralRings, int iring) {
  std::vector<int> triplet; //  Triplet formed from iring
  int ringSize = 6;         // Here, all the rings are hexagons
  int j;                    // Used for making the triplet
  int jring;                // Peripheral ring ID to be searched
  int count;                // Number of  rings found that match the triplet
  std::vector<int>
      newPeripherals; // Vector in which the new peripheral ring IDs are saved.
                      // This will be swapped with peripheralRings later

  for (int k = 0; k < ringSize; k++) {
    triplet.clear(); // Clear the triplet
    // Get a triplet
    for (int i = k; i < k + 3; i++) {
      j = i;
      if (i >= ringSize) {
        j = i - ringSize;
      }
      triplet.push_back(rings[iring][j]);
    } // end of getting a triplet from k
    // -------------
    // Compare the triplet with every possible peripheral
    // ring inside peripheralRings.
    count = 0; // init to zero
    // Loop through all possible peripheral rings
    for (int m = 0; m < (*peripheralRings).size(); m++) {
      jring = (*peripheralRings)[m]; // Ring ID of ring to be searched
      // Search inside the ring with index jring for the triplet
      bool foundTriplet = ring::findTripletInRing(rings[jring], triplet);

      // If the triplet has been found inside jring
      if (foundTriplet) {
        newPeripherals.push_back(jring); // Update new peripheral vector
        count++;
        break;
      } // end of ring found
    }   // end of loop through all possible peripheral rings
    // If count is 0, then the triplet was not found in any peripheral ring
    if (count == 0) {
      return false;
    } // Return false since the triplet was not found
    // -------------
  } // end of looping through 0-6 to get triplets

  // Swap the old peripheral rings vector with the new one
  (*peripheralRings).swap(newPeripherals);

  // If there are more than 6 peripheral rings, the code will fail
  // Comment this out if you want
  if ((*peripheralRings).size() > 6) {
    std::cerr
        << "There are more than 6 peripheral rings. The code will fail. \n";
    return false;
  } // end of check for more than 6 peripherals

  return true;
}

/**
 * @details For a given ring, which is being tested as the equatorial ring,
 * this function tests the following, given peripheralRings stored in increasing
 * order of the triplet starting element:
 * 1. Rings corresponding to T1, T3, T5 should have at least one common element.
 * 2. Rings corresponding to T2, T4, T6 should have at least one common element.
 * 3. The following rings should have at least three common elements- {T1, T3},
 * {T2, T4}, {T3, T5}, {T4, T6}.
 * Internally calls the following functions:
 * - ring::commonElementsInThreeRings (CONDITION 1: Rings corresponding to T1,
 * T3, T5 should have at least one common element).
 * - ring::commonElementsInThreeRings (CONDITION 2: Rings corresponding to T1,
 * T3, T5 should have at least one common element).
 * - ring::findsCommonElements (Gets the common elements between a pair of
 * rings).
 * @param[in] rings Vector of vectors containing the 6-membered primitive
 *  rings.
 * @param[in] peripheralRings Vector containing the indices of rings which are
 *  potential peripheral rings.
 * @return A bool; true if ring being tested as equatorial satisfies
 *  the above condition for being an equatorial ring, and false otherwise.
 */
bool ring::conditionThreeDDC(std::vector<std::vector<int>> rings,
                             std::vector<int> *peripheralRings) {
  // New
  std::vector<int> common; // Vector containing common elements
  bool hasCommon;          // true if the rings have a common element
  int iring, jring;        // Pairs of peripheral rings
  // ----------------------------------------------------------------------------
  // CONDITION 1: Rings corresponding to T1, T3, T5 should have at least one
  // common element.
  hasCommon = ring::commonElementsInThreeRings(rings[(*peripheralRings)[0]],
                                               rings[(*peripheralRings)[2]],
                                               rings[(*peripheralRings)[4]]);

  // If T1, T3, T5 don't have a common element, return false
  if (!hasCommon) {
    return false;
  } // not a DDC
  // ----------------------------------------------------------------------------
  // CONDITION 2: Rings corresponding to T1, T3, T5 should have at least one
  // common element.
  hasCommon = ring::commonElementsInThreeRings(rings[(*peripheralRings)[1]],
                                               rings[(*peripheralRings)[3]],
                                               rings[(*peripheralRings)[5]]);

  // If T1, T3, T5 don't have a common element, return false
  if (!hasCommon) {
    return false;
  } // not a DDC
  // ----------------------------------------------------------------------------
  // CONDITION 3: Rings corresponding to {T1, T3}, {T2, T4}, {T3, T5}, {T4,
  // T6}
  // must have three elements in common amongst them

  // Loops to get pairs of rings corresponding to the right triplets
  for (int i = 0; i <= 3; i++) {
    common.clear(); // init
    // Pairs of rings corresponding to triplets.
    iring = (*peripheralRings)[i];
    jring = (*peripheralRings)[i + 2];
    // Get the common elements
    common = ring::findsCommonElements(rings[iring], rings[jring]);
    // There should be at least three elements
    if (common.size() < 3) {
      return false;
    } // not a DDC
  }   // end of getting iring and jring
  // ----------------------------------------------------------------------------

  // iring is an equatorial ring and peripheralRings has the 6 peripheral rings
  return true;
}

/**
 * @details Determines which hexagonal rings are HCs. This function
 * returns a vector which contains the ring IDs of all the rings which are HC
 * rings. The ring IDs correspond to the index of the rings inside the vector of
 * vector rings, starting from 0. HC rings can be found using a three-step
 * procedure, in which first two basal rings are found. Prismatic rings are
 * simply rings which share every face made by upper and lower triplets of the
 * basal rings The neighbour list is also required as an input, which is a
 * vector of vectors, containing atom IDs. The first element of the neighbour
 * list is the atomID of the atom for which the other elements are nearest
 * neighbours. The function calls the following functions internally:
 * - ring::hasCommonElements (Step one: Checks to see if basal1 and basal2, two
 * candidate basal rings, have common elements or not. If they don't, then they
 * cannot be basal rings).
 * - ring::basalConditions (Step two and three: One of the elements of basal2
 * must be the nearest neighbour of either the first (index0; \f$ l_1 \f$) or
 * second (index1; \f$ l_2 \f$) element of basal1. If \f$ m_k \f$ is the nearest
 * neighbour of \f$ l_1 \f$, \f$ m_{k+2} \f$ and \f$ m_{k+4} \f$ must be
 * neighbours of \f$ l_3 \f$ and \f$ l_5 \f$(\f$ l_5 \f$ or \f$ l_3 \f$).
 * Similarly modified for \f$ l_2 \f$).
 * @param[in] rings Vector of vectors containing 6-membered primitive rings
 *  (wherein each ring contains the atom indices of the particles which
 * constitute the ring).
 * @param[in] ringType Vector containing a ring::strucType (structure type) for
 *  each ring.
 * @param[in] nList Row-ordered vector of vectors of the neighbour list (by
 *  index).
 * @param[in] cageList Vector in which every cage is saved.
 * @return A vector of all the ring indices which constitute HCs.
 */
std::vector<int> ring::findHC(std::vector<std::vector<int>> rings,
                              std::vector<ring::strucType> *ringType,
                              std::vector<std::vector<int>> nList,
                              std::vector<cage::Cage> *cageList) {
  std::vector<int> listHC;
  int totalRingNum = rings.size(); // Total number of hexagonal rings
  std::vector<int> basal1;         // First basal ring
  std::vector<int> basal2;         // Second basal ring
  bool cond1, cond2; // Conditions for rings to be basal (true) or not (false)
  std::vector<int>
      HCRings; // Indices of rings which constitute a single HC, with the basal
               // rings first, followed by prismatic rings
  std::vector<int> prismaticRings; // Ring indices of prismatic rings
  int kring;                       // Ring index of the prismatic rings
  std::vector<bool>
      isPrismatic; // Flag for checking if the ring is prismatic (true) or not
                   // (false), since the basal rings are checked
  isPrismatic.resize(totalRingNum); // Initialized to false

  // Two loops through all the rings are required to find pairs of basal rings
  for (int iring = 0; iring < totalRingNum - 1; iring++) {
    // -----------------------
    // Skip if iring is prismatic
    if (isPrismatic[iring]) {
      continue;
    } // Skip if prismatic
    // -----------------------
    cond1 = false;
    cond2 = false;
    basal1 = rings[iring]; // Assign iring to basal1
    // Loop through the other rings to get a pair
    for (int jring = iring + 1; jring < totalRingNum; jring++) {
      // -----------------------
      // Skip if iring is prismatic
      if (isPrismatic[jring]) {
        continue;
      }                      // Skip if prismatic
                             // -----------------------
      basal2 = rings[jring]; // Assign jring to basal2
      // ------------
      // Step one: Check to see if basal1 and basal2 have common
      // elements or not. If they don't, then they cannot be basal rings
      cond1 = ring::hasCommonElements(basal1, basal2);
      if (cond1 == true) {
        continue;
      }
      // -----------
      // Step two and three: One of the elements of basal2 must be the nearest
      // neighbour of either the first (index0; l1) or second (index1; l2)
      // element of basal1. If m_k is the nearest neighbour of l1, m_{k+2} and
      // m_{k+4} must be neighbours of l3 and l5(l5 or l3). Modify for l2.
      cond2 = ring::basalConditions(nList, &basal1, &basal2);
      if (cond2 == false) {
        continue;
      }
      // -----------
      // iring and jring are basal rings!
      // Update iring
      if ((*ringType)[iring] == ring::strucType::unclassified) {
        (*ringType)[iring] = ring::strucType::HCbasal;
        listHC.push_back(iring);
      } else if ((*ringType)[iring] == ring::strucType::DDC) {
        (*ringType)[iring] = ring::strucType::bothBasal;
        listHC.push_back(iring);
      }
      // Update jring
      if ((*ringType)[jring] == ring::strucType::unclassified) {
        (*ringType)[jring] = ring::strucType::HCbasal;
        listHC.push_back(jring);
      } else if ((*ringType)[jring] == ring::strucType::DDC) {
        (*ringType)[jring] = ring::strucType::bothBasal;
        listHC.push_back(jring);
      }
      // Find the prismatic rings
      prismaticRings.clear(); // Clear the prismatic ring vector first
      ring::findPrismatic(rings, &listHC, ringType, iring, jring,
                          &prismaticRings);
      // Update the prismatic rings
      for (int k = 0; k < prismaticRings.size(); k++) {
        kring =
            prismaticRings[k]; // Current ring index of the (3) prismatic rings
        // Update kring
        if ((*ringType)[kring] == ring::strucType::unclassified) {
          (*ringType)[kring] = ring::strucType::HCprismatic;
          listHC.push_back(kring);
        } else if ((*ringType)[kring] == ring::strucType::DDC) {
          (*ringType)[kring] = ring::strucType::bothPrismatic;
          listHC.push_back(kring);
        }
        //
        // Update the isPrismatic vector
        isPrismatic[kring] = true;
        //
      } // End of update of prismatic rings in listHC
      // -----------
      // Update the cageList vector of Cages
      // Update the basal rings
      HCRings.clear();
      HCRings.push_back(iring);
      HCRings.push_back(jring);
      // Add the prismaticRings
      HCRings.insert(std::end(HCRings), std::begin(prismaticRings),
                     std::end(prismaticRings));
      (*cageList).push_back({cage::cageType::HexC, HCRings});
      // -----------
    } // end of loop through rest of the rings to get the second basal ring
  }   // end of loop through all rings for first basal ring

  sort(listHC.begin(), listHC.end());
  auto ip = std::unique(listHC.begin(), listHC.end());
  // Resize peripheral rings to remove undefined terms
  listHC.resize(std::distance(listHC.begin(), ip));

  return listHC;
}

/**
 * @details Check to see if two basal rings are HCs or not, using the neighbour
 * list information. The neighbour list nList is a vector of vectors, containing
 * atom indices (not atom IDs!). The first element of each subvector in nList is
 * the atom index of the particle for which the other elements are the nearest
 * neighbours.
 * Internally, the following functions are called:
 * - ring::basalNeighbours (CONDITION1: @f$ m_{k+2} @f$ and @f$ m_{k+4} @f$ must
 * be bonded to @f$ l_3 @f$ and @f$ l_5 @f$ (if @f$ l_1 @f$ is a neighbour) or
 * @f$ m_{k+2} @f$ and @f$ m_{k+4} @f$ must be bonded to @f$ l_4 @f$ and @f$ l_6
 * @f$ (if @f$ l_2 @f$ is a neighbour)).
 * - ring::notNeighboursOfRing (CONDITION2: @f$ m_{k+1} @f$, @f$ m_{k+3} @f$ and
 * @f$ m_{k+5}@f$ must NOT be bonded to any element in basal1).
 * @param[in] nList Vector of vectors of the neighbour list (by index).
 * @param[in] basal1 Vector containing the first candidate basal ring.
 * @param[in] basal2 Vector containing the second candidate basal ring.
 * @return A bool; true if the basal rings being tested fulfill this condition
 *  for being the basal rings of an HC.
 */
bool ring::basalConditions(std::vector<std::vector<int>> nList,
                           std::vector<int> *basal1, std::vector<int> *basal2) {
  int l1 = (*basal1)[0]; // first element of basal1 ring
  int l2 = (*basal1)[1]; // second element of basal1 ring
  int ringSize = 6;      // Size of the ring; each ring contains 6 elements
  int m_k;               // Atom Index (in pointCloud) of element in basal2
  int kIndex;            // Index of m_k in basal2, corresponding to m_k
  int currentKindex;     // Current k index when finding alternating elements of
                         // basal2
  std::vector<int> evenTriplet; // contains m_k, m_{k+2}, m_{k+4}
  std::vector<int> oddTriplet;  // contains m_{k+1}, m_{k+3}, m_{k+5}
  int compare1, compare2;       // l3 and l5 OR l4 and l6
  int index;
  bool l1_neighbour, l2_neighbour; // m_k is a neighbour of l1(true) or not
                                   // (false); m_k is a neighbour of l2(true)
  bool isNeigh, notNeigh; // Used to check if the rings are basal or not

  // ---------------------------------------------
  // SEARCH FOR L1_NEIGHBOUR OR L2_NEIGHBOUR
  // Search for whether an element of basal2 is a neighbour of l1 or l2
  for (int k = 0; k < ringSize; k++) {
    // init
    l1_neighbour = false;
    l2_neighbour = false;
    m_k = (*basal2)[k];

    // ---------------
    // CHECK IF M_K MATCHES L1 NEIGHBOURS
    auto it1 = std::find(nList[l1].begin() + 1, nList[l1].end(), m_k);
    // If m_k was found in l1's nList
    if (it1 != nList[l1].end()) {
      compare1 = (*basal1)[2]; // l3
      compare2 = (*basal1)[4]; // l5
      kIndex = k;              // Saving the array index of m_k
      l1_neighbour = true;
      break;
    } // m_k found in l1's nList
    // ---------------
    // CHECK IF M_K MATCHES L2 NEIGHBOURS
    auto it2 = std::find(nList[l2].begin() + 1, nList[l2].end(), m_k);
    // If m_k was found in l1's nList
    if (it2 != nList[l2].end()) {
      compare1 = (*basal1)[3]; // l4
      compare2 = (*basal1)[5]; // l6
      kIndex = k;              // Saving the array index of m_k
      l2_neighbour = true;
      break;
    } // m_k found in l1's nList
    // ---------------
  } // End of search for basal2 elements in l1 or l2's nList
  // ---------------------------------------------

  // Return false if neither l1 nor l2 have any neighbours
  // in basal2

  if (l1_neighbour == false && l2_neighbour == false) {
    return false;
  } // basal conditions not fulfilled

  // Get the alternating elements starting with kIndex.
  // 'evenTriplet': m_k, m_{k+2}, m_{k+4} - neighbours of compare1 and compare2.
  // 'oddTriplet': m_{k+1}, m_{k+3}, m_{k+5}- cannot be neighbours of basal1
  for (int k = 0; k <= 5; k++) {
    currentKindex = kIndex + k; // k
    // Wrap-around
    if (currentKindex >= ringSize) {
      currentKindex -= ringSize;
    } // end of wrap-around of k
    //
    // Update 'evenTriplet'
    if (k % 2 == 0) {
      evenTriplet.push_back((*basal2)[currentKindex]);
    } // end of update of evenTriplet
    // Update 'oddTriplet'
    else {
      oddTriplet.push_back((*basal2)[currentKindex]);
    } // end of update of oddTriplet
  }   // End of getting alternating triplets

  // ---------------------------------------------
  // CONDITION1: m_{k+2} and m_{k+4} must be bonded to l3 and l5 (if l1 is a
  // neighbour) or m_{k+2} and m_{k+4} must be bonded to l4 and l6 (if l2 is a
  // neighbour) Basically, this boils down to checking whether compare1 and
  // compare2 are in the neighbour lists of the last two elements of evenTriplet

  isNeigh = ring::basalNeighbours(nList, &evenTriplet, compare1, compare2);

  // If condition1 is not true, then the candidate
  // rings are not part of an HC
  if (!isNeigh) {
    return false;
  } // Not an HC

  // ---------------------------------------------
  // CONDITION2: m_{k+1}, m_{k+3} and m_{k+5} must NOT be bonded to any element
  // in basal1.
  // Basically, this boils down to checking whether the elements of oddTriplet
  // are in the neighbour lists of all the elements of basal1.

  // condition 2. This must be true for an HC
  notNeigh = ring::notNeighboursOfRing(nList, &oddTriplet, basal1);

  // If condition2 is not true, the the candidate rings
  // are not part of an HC
  if (!notNeigh) {
    return false;
  } // Not an HC

  // Otherwise, all the conditions are true and this is an HC
  return true;

} // end of function

/**
 * @details Tests whether the last two elements of a triplet are neighbours of
 *two atom indices which have been passed in as inputs.
 * @param[in] nList Vector of vectors of the neighbour list (by index).
 * @param[in] triplet Vector containing the current triplet being tested.
 * @param[in] atomOne Index of the first atom.
 * @param[in] atomTwo Index of the second atom.
 * @return A bool; true if the condition is met and false otherwise.
 */
bool ring::basalNeighbours(std::vector<std::vector<int>> nList,
                           std::vector<int> *triplet, int atomOne,
                           int atomTwo) {
  // Search for needles in a haystack :)
  int needle1 = (*triplet)[1];
  int needle2 = (*triplet)[2];
  bool neighbourFound = false;
  bool neighOne = false,
       neighTwo = false; // Neighbour for atomOne, or neighbour for atomTwo
  // ----------------------------
  // For first element needle1, which must belong to either atomOne's or
  // atomTwo's neighbour list Search atomOne's neighbours
  auto it =
      std::find(nList[atomOne].begin() + 1, nList[atomOne].end(), needle1);
  if (it != nList[atomOne].end()) {
    neighbourFound = true;
    neighOne = true;
  } // atomOne's neighbour
  // If it is not atomOne's neighbour, it might be atomTwo's neighbour
  if (!neighOne) {
    it = std::find(nList[atomTwo].begin() + 1, nList[atomTwo].end(), needle1);
    if (it != nList[atomTwo].end()) {
      neighbourFound = true;
      neighTwo = true;
    } // end of check to see if neighbour was found
  }   // End of check to see if needle1 is atomTwo's neighbour
  // ----------------------------
  // If needle1 is not a neighbour of atomOne or atomTwo, return false
  if (neighbourFound == false) {
    return false;
  }

  // Check to see if needle2 is a neighbour of either atomOne or atomTwo
  // ===============
  // if atomOne was a neighbour of needle1, needle2 must be a neighbour of
  // atomTwo
  if (neighOne) {
    it = std::find(nList[atomTwo].begin() + 1, nList[atomTwo].end(), needle2);
    // It is a neighbour
    if (it != nList[atomTwo].end()) {
      return true;
    }
    // It is not a neighbour
    else {
      return false;
    }
  } // End of check for neighbour of atomTwo
  // ===============
  // if atomTwo was a neighbour of needle1, needle2 must be a neighbour of
  // atomOne
  else {
    it = std::find(nList[atomOne].begin() + 1, nList[atomOne].end(), needle2);
    // It is a neighbour
    if (it != nList[atomOne].end()) {
      return true;
    }
    // It is not a neighbour
    else {
      return false;
    }
  }
  // ===============
}

/**
 * @details Checks to make sure that the elements of the triplet are NOT
 * neighbours of any elements inside a vector (ring) passed in (false)
 * If any of them are neighbours, this function returns false.
 * @param[in] nList Vector of vectors of the neighbour list (by index).
 * @param[in] triplet Vector containing the current triplet being tested.
 * @param[in] ring Ring passed in.
 * @return A bool; true if the condition is met and false otherwise.
 */
bool ring::notNeighboursOfRing(std::vector<std::vector<int>> nList,
                               std::vector<int> *triplet,
                               std::vector<int> *ring) {
  int iatom; // AtomID of the atom to be searched for inside the neighbour
             // lists
  int jatom; // AtomID of in whose neighbour list iatom will be searched for
  std::vector<int>::iterator it;

  for (int i = 0; i < (*triplet).size(); i++) {
    iatom = (*triplet)[i]; // AtomID to be searched for
    // iatom must be searched for in the neighbour lists of all elements
    // of the ring vector
    for (int j = 0; j < (*ring).size(); j++) {
      jatom = (*ring)[j];
      // ------------------
      // Search for iatom in the neighbour list of jatom
      it = std::find(nList[jatom].begin() + 1, nList[jatom].end(), iatom);
      // It is a neighbour!
      if (it != nList[jatom].end()) {
        return false;
      }
      // ------------------
    } // end of loop through every element of ring
  }   // end of loop through all triplet elements

  return true;
}

/**
 * @details Finding prismatic rings when passed the information in the
 *ringType input vector.
 * @param[in] rings The 6-membered primitive rings.
 * @param[in] listHC Vector containing the ring indices of rings which are part
 *  of HCs.
 * @param[in] ringType Contains a structure type for each ring.
 * @param[in] iring Index of the \f$ i^{th} \f$ ring.
 * @param[in] jring Index of the \f$ j^{th} \f$ ring.
 * @param[in] prismaticRings Vector containing the indices of rings which are
 *  prismatic.
 */
int ring::findPrismatic(std::vector<std::vector<int>> rings,
                        std::vector<int> *listHC,
                        std::vector<ring::strucType> *ringType, int iring,
                        int jring, std::vector<int> *prismaticRings) {
  int iIndex, jIndex;             // Used for making rings to be searched
  int ringSize = rings[0].size(); // This is 6 for hexagons
  std::vector<int> iTriplet;      // triplet formed from iring
  std::vector<int> jTriplet;      // triplet formed from jring
  std::vector<int> common;        // Common elements

  // Make all possible triplets out of iring
  for (int i = 0; i < ringSize; i++) {
    // init
    iTriplet.clear();
    // ------
    // Get a triplet
    for (int m = 0; m < 3; m++) {
      iIndex = i + m;
      if (iIndex >= ringSize) {
        iIndex -= ringSize;
      }
      iTriplet.push_back(rings[iring][iIndex]);
    } // end of getting a triplet from iring

    // -------------------------------------------
    // Now that a triplet has been found, find all rings with that triplet in
    // it!
    for (int kring = 0; kring < rings.size(); kring++) {
      // If this is the same as iring or jring, skip
      if (kring == iring || kring == jring) {
        continue;
      } // is not prismatic
      //
      // Now find out whether kring has the triplet or not!
      common = ring::findsCommonElements(iTriplet, rings[kring]);

      // If this triplet is not shared by  kring
      // skip
      if (common.size() != 3) {
        continue;
      } // kring does not have iTriplet

      // -----------------
      // If kring does have the triplet, check to see if at least three other
      // elements of kring are shared by jring
      jTriplet.clear(); // init
      // Make jTriplet
      for (int j = 0; j < ringSize; j++) {
        int katom = rings[kring][j];
        auto it = std::find(iTriplet.begin(), iTriplet.end(), katom);

        // If not found, add it to jTriplet
        if (it == iTriplet.end()) {
          jTriplet.push_back(katom);
        } // update jTriplet
      }   // end of making jTriplet out of kring
      // -----------------

      // Now search for jTriplet inside jring
      common = ring::findsCommonElements(jTriplet, rings[jring]);

      // Update the prismatic rings
      if (common.size() == 3) {
        //
        (*listHC).push_back(kring);         // Update listHC vector
        (*prismaticRings).push_back(kring); // Update prismatic rings
        // Update the type inside ringType
        // If the ring is already a DDC ring, it is a mixed ring
        if ((*ringType)[kring] == ring::strucType::DDC) {
          (*ringType)[kring] = ring::strucType::bothPrismatic;
        }
        // If it is unclassified, it is just a prismatic ring
        if ((*ringType)[kring] == ring::strucType::unclassified) {
          (*ringType)[kring] = ring::strucType::HCprismatic;
        } // end ring update
      }   // add kring to the list of prismatic rings
    }     // end of searching through rings for kring
    // -------------------------------------------
  } // end of getting pairs of triplets

  return 0;
}

/**
 * @details Determines which hexagonal rings are both DDCs and HCs. This
 * function returns a vector which contains the ring IDs of all the rings which
 * are both. The ring IDs correspond to the index of the rings inside the vector
 * of vector rings, starting from 0. Rings which are both are of enum type
 * bothBasal OR bothPrismatic. Reassign rings which are mixed in listDDC and
 * listHC as the dummy value -10.
 * @param[in] rings The 6-membered primitive rings.
 * @param[in] ringType Structure type for every ring.
 * @param[in] listDDC Contains a vector of ring indices of DDCs.
 * @param[in] listHC Contains a vector of ring indices of HCs.
 * @return A vector containing the indices of rings which are mixed.
 */
std::vector<int> ring::findMixedRings(std::vector<std::vector<int>> rings,
                                      std::vector<strucType> *ringType,
                                      std::vector<int> *listDDC,
                                      std::vector<int> *listHC) {
  std::vector<int> listMixed;
  int dummyValue = -10;

  // Loop through all rings in the ringType and
  // adds the ring Indices of all rings which are both DDCs and HCs
  for (int iring = 0; iring < (*ringType).size(); iring++) {
    // If iring is of mixed type, add it to the listMixed vector
    if ((*ringType)[iring] == ring::strucType::bothBasal ||
        (*ringType)[iring] == ring::strucType::bothPrismatic) {
      listMixed.push_back(iring);

      //-----------------
      // Search for iring in listDDC
      std::sort((*listDDC).begin(), (*listDDC).end());
      auto iter = std::find((*listDDC).begin(), (*listDDC).end(), iring);
      if (iter != (*listDDC).end()) {
        *iter = dummyValue; // Assign dummy value to the mixed ring
      }                     // found in listDDC
      //-----------------
      //-----------------
      // Search for iring in listHC
      std::sort((*listHC).begin(), (*listHC).end());
      auto itr = std::find((*listHC).begin(), (*listHC).end(), iring);
      if (itr != (*listHC).end()) {
        *itr = dummyValue; // Assign dummy value to the mixed ring
      }                    // found in listHC
      //-----------------

    } // end of check for type
  }   // end of loop through all every ring

  return listMixed;
}

/**
 * @details Assigns a type (cage::iceType) to each atom, according to the
 * previous classification of every ring in ringType. Each subring or vector
 * inside the vector of vector rings, is by index, meaning that the atoms are
 * saved by their indices starting from 0 in the PointCloud.
 * @param[in] rings Vector of vectors of 6-membered rings.
 * @param[in] ringType Vector containing the structural type for each ring.
 * @param[in] atomTypes Structural type for each atom.
 */
int ring::getAtomTypesTopoBulk(std::vector<std::vector<int>> rings,
                               std::vector<ring::strucType> ringType,
                               std::vector<cage::iceType> *atomTypes) {
  cage::iceType iRingType;        // Type of the current ring
  int iatom;                      // Current ring
  int ringSize = rings[0].size(); // Size of the ring

  // Loop through every ring in ringType
  for (int iring = 0; iring < ringType.size(); iring++) {
    //
    // Skip if the ring is unclassified
    if (ringType[iring] == ring::strucType::unclassified) {
      continue;
    } // skip for unclassified rings

    // ------------
    // Get the current ring type
    // DDC
    if (ringType[iring] == ring::strucType::DDC) {
      iRingType = cage::iceType::ddc;
    } // DDC atoms
    //
    // HC
    else if (ringType[iring] == ring::strucType::HCbasal ||
             ringType[iring] == ring::strucType::HCprismatic) {
      iRingType = cage::iceType::hc;
    } // HC atoms
    //
    // Mixed
    else if (ringType[iring] == ring::strucType::bothBasal ||
             ringType[iring] == ring::strucType::bothPrismatic) {
      iRingType = cage::iceType::mixed;
    } // HC atoms
    // Prism
    else if (ringType[iring] == ring::strucType::Prism ||
             ringType[iring] == ring::strucType::deformedPrism ||
             ringType[iring] == ring::strucType::mixedPrismRing) {
      iRingType = cage::iceType::pnc; // 5 membered pnc
    }                        // prism
    // Should never go here
    else {
      continue;
    } //
    // ------------
    // Otherwise, loop through every inside the ring and assign atomTypes the
    // iRingType
    for (int i = 0; i < ringSize; i++) {
      iatom = rings[iring][i]; // Atom index in ring
      if ((*atomTypes)[iatom] == cage::iceType::mixed ||
          (*atomTypes)[iatom] == cage::iceType::mixed2) {
        continue;
      } // Don't reassign
      // For atoms shared by PNCs and DDCs/HCs
      if (ringSize == 6) {
        if ((*atomTypes)[iatom] == cage::iceType::pnc) {
          (*atomTypes)[iatom] = cage::iceType::mixed2;
        } else {
          (*atomTypes)[iatom] = iRingType;
        }
      } else {
        (*atomTypes)[iatom] = iRingType;
      }
    } // end of loop thorugh the current ring
  }   // end of loop through every ring

  return 0;
}

/**
 * @details Determines the number of HCs, DDCs from the cageList vector,
 * containing a list of cages.
 * The number of mixed rings, prismatic rings and basal rings are obtained
 * from the ringType vector.
 * @param[in] ringType Vector containing the structural type for each ring.
 * @param[in] cageList Contains all the cages (DDCs and HCs).
 * @param[in] numHC The number of HCs.
 * @param[in] numDDC The number of DDCs.
 * @param[in] mixedRings The number of mixed rings.
 * @param[in] prismaticRings The number of prismatic rings (of HCs/mixed
 *  rings).
 * @param[in] basalRings TThe number of basal rings (of HCs/mixed
 *  rings).
 */
int ring::getStrucNumbers(std::vector<ring::strucType> ringType,
                          std::vector<cage::Cage> cageList, int *numHC,
                          int *numDDC, int *mixedRings, int *prismaticRings,
                          int *basalRings) {
  // Determines the number of HCs, DDCs, Mixed rings, prismatic and basal rings
  // Init
  *numHC = 0;
  *numDDC = 0;
  *mixedRings = 0;
  *prismaticRings = 0;
  *basalRings = 0;
  // ------------------------------------
  // GETTING THE CAGES (DDCs and HCs)
  // Loop through cages
  for (int icage = 0; icage < cageList.size(); icage++) {
    //
    // HC
    if (cageList[icage].type == cage::cageType::HexC) {
      *numHC += 1;
    } // end of updating HC number
    //
    // DDC
    if (cageList[icage].type == cage::cageType::DoubleDiaC) {
      *numDDC += 1;
    } // end of updating DDC number
  }   // end of loop through cages
  // ------------------------------------
  // GETTING THE RINGSS (Mixed, Prismatic and Basal rings)
  // Loop through the rings
  for (int iring = 0; iring < ringType.size(); iring++) {
    // Mixed
    if (ringType[iring] == ring::strucType::bothBasal ||
        ringType[iring] == ring::strucType::bothPrismatic) {
      *mixedRings += 1;
      // Also update basal rings
      if (ringType[iring] == ring::strucType::bothBasal) {
        *basalRings += 1;
      } // mixed basal rings
      // Also update prismatic rings
      if (ringType[iring] == ring::strucType::bothPrismatic) {
        *prismaticRings += 1;
      } // mixed prismatic rings
    }   // end of updating mixed
    //
    // HCs
    if (ringType[iring] == ring::strucType::HCprismatic) {
      *prismaticRings += 1;
    } // HC prismatic
    // basal HCs
    if (ringType[iring] == ring::strucType::HCbasal) {
      *basalRings += 1;
    } // HC basal
  }   // end of loop through every ring
  // ------------------------------------

  return 0;
} // end of function

/**
 * @details Determines which rings are n-sided prisms. This function
 * returns a vector which contains the ring indices of all the rings which are
 * prisms. The ring indices correspond to the index of the rings inside the
 * vector of vector rings, starting from 0. Prism rings can be found using a
 * three-step procedure, in which first two basal rings are found. Prismatic
 * rings are simply rings which share every face made by upper and lower
 * triplets of the basal rings The neighbour list is also required as an input,
 * which is a vector of vectors, containing atom IDs. The first element of the
 * neighbour list is the atom index of
 * the atom for which the other elements are nearest neighbours.\
 * @param[in] rings The input vector of vectors containing the primitive rings
 *  of a single ring size (number of nodes).
 * @param[in] ringType A vector containing a ring::strucType value (a
 *  classification type) for each ring.
 * @param[in] nPrisms The number of prism blocks identified for the number of
 *  nodes.
 * @param[in] nList The row-ordered neighbour list (by atom index).
 * @param[in] yCloud The input PointCloud.
 * @return A vector containing the ring indices of all the rings which have
 *  been classified as prisms. The indices are with respect to the input rings
 *  vector of vectors.
 */
int prism3::findBulkPrisms(
    std::vector<std::vector<int>> rings, std::vector<ring::strucType> *ringType,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<double> *rmsdPerAtom, double heightCutoff) {
  int totalRingNum = rings.size(); // Total number of rings
  std::vector<int> basal1;         // First basal ring
  std::vector<int> basal2;         // Second basal ring
  bool cond1, cond2; // Conditions for rings to be basal (true) or not (false)
  int ringSize = rings[0].size(); // Number of nodes in each ring
  // Matrix for the reference ring for a given ringSize.
  Eigen::MatrixXd refPointSet(ringSize, 3);

  int axialDim = 2; // Default=z
  refPointSet = pntToPnt::getPointSetRefRing(ringSize, axialDim);
  //

  // Two loops through all the rings are required to find pairs of basal rings
  for (int iring = 0; iring < totalRingNum - 1; iring++) {
    cond1 = false;
    cond2 = false;
    basal1 = rings[iring]; // Assign iring to basal1
    // Loop through the other rings to get a pair
    for (int jring = iring + 1; jring < totalRingNum; jring++) {
      basal2 = rings[jring]; // Assign jring to basal2
      // ------------
      // Put extra check for axial basal rings if shapeMatching is being done
      // ------------
      // Step one: Check to see if basal1 and basal2 have common
      // elements or not. If they don't, then they cannot be basal rings
      cond1 = ring::hasCommonElements(basal1, basal2);
      if (cond1 == true) {
        continue;
      }

      // ------------
      bool smallDist =
          prism3::basalRingsSeparation(yCloud, basal1, basal2, heightCutoff);
      if (!smallDist) {
        continue;
      } // the basal rings are too far apart

      // Otherwise
      // Do shape matching here
      bool isPrism = match::matchUntetheredPrism(yCloud, nList, refPointSet,
                                                 &basal1, &basal2, rmsdPerAtom);

      // Success! The rings are basal rings of a prism!
      if (isPrism) {
        //
        // Update iring
        if ((*ringType)[iring] == ring::strucType::unclassified) {
          (*ringType)[iring] = ring::strucType::Prism;
        }
        // Update jring
        if ((*ringType)[jring] == ring::strucType::unclassified) {
          (*ringType)[jring] = ring::strucType::Prism;
        }
      } // end of reduced criteria
      // Strict criteria
      else {
        cond2 = prism3::basalPrismConditions(nList, &basal1, &basal2);
        // If the condition is false then the strict criterion has not been met
        if (!cond2) {
          continue;
        }
        // Update iring
        if ((*ringType)[iring] == ring::strucType::unclassified) {
          (*ringType)[iring] = ring::strucType::Prism;
        }
        // Update jring
        if ((*ringType)[jring] == ring::strucType::unclassified) {
          (*ringType)[jring] = ring::strucType::Prism;
        }
        //
        // Shape-matching to get the RMSD (if shape-matching is desired)

        // bool isKnownPrism = match::matchPrism(
        //     yCloud, nList, refPointSet, &basal1, &basal2, rmsdPerAtom, true);

        // -----------
      } // end of strict criteria

    } // end of loop through rest of the rings to get the second basal ring
  }   // end of loop through all rings for first basal ring

  return 0;
}

/**
 * @details A function that checks to see if two basal rings are basal rings of
 * a prism block or not, using the neighbour list information. The neighbour
 * list nList is a row-ordered vector of vectors, containing atom indices (not
 * atom IDs!). The first element of each subvector in nList is the atom index of
 * the particle for which the other elements are the nearest neighbours.
 * @param[in] nList Row-ordered neighbour list by atom index.
 * @param[in] basal1 The vector for one of the basal rings.
 * @param[in] basal2 The vector for the other basal ring.
 * @return A value that is true if the basal rings constitute a prism block,
 *  and false if they do not make up a prism block.
 */
bool prism3::basalPrismConditions(std::vector<std::vector<int>> nList,
                                  std::vector<int> *basal1,
                                  std::vector<int> *basal2) {
  int l1 = (*basal1)[0]; // first element of basal1 ring
  int ringSize =
      (*basal1).size(); // Size of the ring; each ring contains n elements
  int m_k;              // Atom ID of element in basal2
  bool l1_neighbour;    // m_k is a neighbour of l1(true) or not (false)

  // isNeighbour is initialized to false for all basal2 elements; indication if
  // basal2 elements are neighbours of basal1
  std::vector<bool> isNeighbour(ringSize, false);
  int kIndex;  // m_k index
  int lAtomID; // atomID of the current element of basal1
  int kAtomID; // atomID of the current element of basal2

  // ---------------------------------------------
  // COMPARISON OF basal2 ELEMENTS WITH l1
  for (int k = 0; k < ringSize; k++) {
    l1_neighbour = false;
    m_k = (*basal2)[k];
    // =================================
    // Checking to seee if m_k is be a neighbour of l1
    // Find m_k inside l1 neighbour list
    auto it = std::find(nList[l1].begin() + 1, nList[l1].end(), m_k);

    // If the element has been found, for l1
    if (it != nList[l1].end()) {
      l1_neighbour = true;
      kIndex = k;
      break;
    }
  } // l1 is a neighbour of m_k

  // If there is no nearest neighbour, then the two rings are not part of the
  // prism
  if (!l1_neighbour) {
    return false;
  }

  // ---------------------------------------------
  // NEIGHBOURS of basal1 in basal2
  isNeighbour[kIndex] = true;

  // All elements of basal1 must be neighbours of basal2
  for (int i = 1; i < ringSize; i++) {
    lAtomID = (*basal1)[i]; // element of basal1 ring
    for (int k = 0; k < ringSize; k++) {
      // Skip if already a neighbour
      if (isNeighbour[k]) {
        continue;
      }
      // Get the comparison basal2 element
      kAtomID = (*basal2)[k]; // element of basal2 ring;

      // Checking to see if kAtomID is a neighbour of lAtomID
      // Find kAtomID inside lAtomID neighbour list
      auto it1 =
          std::find(nList[lAtomID].begin() + 1, nList[lAtomID].end(), kAtomID);

      // If the element has been found, for l1
      if (it1 != nList[lAtomID].end()) {
        isNeighbour[k] = true;
      }
    } // Loop through basal2
  }   // end of check for neighbours of basal1

  // ---------------------------------------------

  // They should all be neighbours
  for (int k = 0; k < ringSize; k++) {
    // Check to see if any element is false
    if (!isNeighbour[k]) {
      return false;
    }
  }

  // Everything works out!
  return true;
}

/**
 * @details Relaxed criteria for deformed
 * prism blocks: at least one bond
 * should exist between the basal
 * rings.
 */
bool prism3::relaxedPrismConditions(std::vector<std::vector<int>> nList,
                                    std::vector<int> *basal1,
                                    std::vector<int> *basal2) {
  int ringSize =
      (*basal1).size();     // Size of the ring; each ring contains n elements
  int m_k;                  // Atom ID of element in basal2
  bool isNeighbour = false; // This is true if there is at least one bond
                            // between the basal rings
  int l_k;                  // Atom ID of element in basal1

  // ---------------------------------------------
  // COMPARISON OF basal2 ELEMENTS (m_k) WITH basal1 ELEMENTS (l_k)
  // Loop through all the elements of basal1
  for (int l = 0; l < ringSize; l++) {
    l_k = (*basal1)[l];
    // Search for the nearest neighbour of l_k in basal2
    // Loop through basal2 elements
    for (int m = 0; m < ringSize; m++) {
      m_k = (*basal2)[m];
      // Find m_k inside l_k neighbour list
      auto it = std::find(nList[l_k].begin() + 1, nList[l_k].end(), m_k);

      // If the element has been found, for l1
      if (it != nList[l_k].end()) {
        isNeighbour = true;
        break;
      } // found element
    }   // end of loop through all the elements of basal2

    // If a neighbour has been found then
    if (isNeighbour) {
      return true;
    }
  } // end of loop through all the elements of basal1

  // If a neighbour has not been found, return false
  return false;
}

//! Check to see that candidate basal prisms are not really far from each other
//! Return true if the basal rings are within the heightCutoff
bool prism3::basalRingsSeparation(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2, double heightCutoff) {
  //
  int ringSize = basal1.size();
  int l_k, m_k; // Atom indices
  double infHugeNumber = 100000;
  double leastDist = infHugeNumber;
  int index = -1; // starting index
  // For the first element of basal1:

  l_k = basal1[0]; // This is the atom particle C++ index

  // Search for the nearest neighbour of l_k in basal2
  // Loop through basal2 elements
  for (int m = 0; m < ringSize; m++) {
    m_k = basal2[m]; // Atom index to find in the neighbour list of iatom

    // Calculate the distance
    double dist = gen::periodicDist(yCloud, l_k, m_k);

    // Update the least distance
    if (leastDist > dist) {
      leastDist = dist; // This is the new least distance
      index = m;
    } // end of update of the least distance

  } // found element

  // If the element has been found, for l1
  if (leastDist < heightCutoff) {
    return true;
  } // end of check
  else {
    return false;
  }
}
