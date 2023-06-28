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

#include <topo_one_dim.hpp>

// -----------------------------------------------------------------------------------------------------
// PRISM ALGORITHMS
// -----------------------------------------------------------------------------------------------------

/**
 * @details Function that loops through the primitive rings (which is a vector
 * of vectors) of all sizes, upto maxDepth (the largest ring size). The function
 * returns a vector which contains the ring indices of all the rings which
 * constitute prism blocks. The ring IDs correspond to the index of the rings
 * inside the vector of vector rings, starting from 0. This is registered as a
 * Lua function, and is exposed to the user directly. The function calls the
 * following functions internally:
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
 * @param[in] maxDepth The first frame.
 */
int ring::prismAnalysis(
    std::string path, std::vector<std::vector<int>> rings,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int maxDepth,
    int *atomID, int firstFrame, int currentFrame, bool doShapeMatching) {
  //
  std::vector<std::vector<int>>
      ringsOneType;           // Vector of vectors of rings of a single size
  std::vector<int> listPrism; // Vector for ring indices of n-sided prism
  std::vector<ring::strucType>
      ringType; // This vector will have a value for each ring inside
  std::vector<ring::strucType>
      atomState; // This vector will have a value for each atom, depending on
                 // the ring type
  int nPerfectPrisms;          // Number of perfect prisms of each type
  int nImperfectPrisms;        // Number of deformed prisms of each type
  std::vector<int> nPrismList; // Vector of the values of the number of perfect
                               // prisms for a particular frame
  std::vector<int> nDefPrismList; // Vector of the values of the number of
                                  // deformed prisms for a particular frame
  std::vector<double>
      heightPercent; // Height percent for a particular n and frame
  std::vector<int>
      atomTypes; // contains int values for each prism type considered
  double avgPrismHeight = 2.845; // A value of 2.7-2.85 Angstrom is reasonable
  // Qualifier for the RMSD per atom:
  std::vector<double> rmsdPerAtom;
  // -------------------------------------------------------------------------------
  // Init
  nPrismList.resize(
      maxDepth -
      2); // Has a value for every value of ringSize from 3, upto maxDepth
  nDefPrismList.resize(maxDepth - 2);
  heightPercent.resize(maxDepth - 2);
  // The atomTypes vector is the same size as the pointCloud atoms
  atomTypes.resize(yCloud->nop, 1); // The dummy or unclassified value is 1
  // The rmsdPerAtom vector is the same size as the pointCloud atoms and has an
  // RMSD value for every atom
  rmsdPerAtom.resize(yCloud->nop, -1);
  // Resize the atom state vector
  atomState.resize(yCloud->nop); // Dummy or unclassified
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
      nPrismList[ringSize - 3] = 0;      // Update the number of prisms
      nDefPrismList[ringSize - 3] = 0;   // Update the number of deformed prisms
      heightPercent[ringSize - 3] = 0.0; // Update the height%
      continue;
    } // skip if there are no rings
    //
    // -------------
    // Init of variables specific to ringSize prisms
    listPrism.resize(0);
    ringType.resize(0);
    nPerfectPrisms = 0;
    nImperfectPrisms = 0;
    ringType.resize(
        ringsOneType.size()); // Has a value for each ring. init to zero.
    // -------------
    // Now that you have rings of a certain size:
    // Find prisms, saving the ring indices to listPrism
    listPrism = ring::findPrisms(ringsOneType, &ringType, &nPerfectPrisms,
                                 &nImperfectPrisms, nList, yCloud, &rmsdPerAtom,
                                 doShapeMatching);
    // -------------
    nPrismList[ringSize - 3] =
        nPerfectPrisms; // Update the number of perfect prisms
    nDefPrismList[ringSize - 3] =
        nImperfectPrisms; // Update the number of defromed prisms
    int totalPrisms = nPerfectPrisms + nImperfectPrisms;
    // Update the height% for the phase
    heightPercent[ringSize - 3] =
        topoparam::normHeightPercent(yCloud, totalPrisms, avgPrismHeight);
    // Continue if there are no prism units
    if (nPerfectPrisms + nImperfectPrisms == 0) {
      continue;
    } // skip for no prisms
    // Do a bunch of write-outs and calculations
    // TODO: Write out each individual prism as data files (maybe with an
    // option)
    // Get the atom types for a particular prism type
    ring::assignPrismType(ringsOneType, listPrism, ringSize, ringType,
                          &atomTypes, &atomState);
    // -------------
  } // end of loop through every possible ringSize

  // Calculate the height%

  // Write out the prism information
  sout::writePrismNum(path, nPrismList, nDefPrismList, heightPercent, maxDepth,
                      yCloud->currentFrame, firstFrame);
  // Reassign the prism block types for the deformed prisms
  if (doShapeMatching) {
    ring::deformedPrismTypes(atomState, &atomTypes, maxDepth);
  } // reassign prism block types for deformed prisms

  // Get rid of translations along the axial dimension
  ring::rmAxialTranslations(yCloud, atomID, firstFrame, currentFrame);

  // Write out dump files with the RMSD per atom for the shape-matching
  // criterion
  if (doShapeMatching) {
    sout::writeLAMMPSdumpINT(yCloud, rmsdPerAtom, atomTypes, maxDepth, path);
  } // reassign prism block types for deformed prisms

  // Write out the lammps data file for the particular frame
  sout::writeLAMMPSdataAllPrisms(yCloud, nList, atomTypes, maxDepth, path,
                                 doShapeMatching);

  return 0;
}

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
std::vector<int>
ring::findPrisms(std::vector<std::vector<int>> rings,
                 std::vector<ring::strucType> *ringType, int *nPerfectPrisms,
                 int *nImperfectPrisms, std::vector<std::vector<int>> nList,
                 molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 std::vector<double> *rmsdPerAtom, bool doShapeMatching) {
  std::vector<int> listPrism;
  int totalRingNum = rings.size(); // Total number of rings
  std::vector<int> basal1;         // First basal ring
  std::vector<int> basal2;         // Second basal ring
  bool cond1, cond2; // Conditions for rings to be basal (true) or not (false)
  bool relaxedCond;  // Condition so that at least one bond exists between the
                     // two basal rings
  bool isAxialPair;  // Basal rings should be parallel in one dimension to
                     // prevent overcounting
  int ringSize = rings[0].size(); // Number of nodes in each ring
  *nImperfectPrisms = 0;          // Number of undeformed prisms
  *nPerfectPrisms = 0;            // Number of undeformed prisms
  // Matrix for the reference ring for a given ringSize.
  Eigen::MatrixXd refPointSet(ringSize, 3);

  // Get the reference ring point set for a given ring size.
  // Get the axial dimension
  int axialDim = std::max_element(yCloud->box.begin(), yCloud->box.end()) -
                 yCloud->box.begin();
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
      if (doShapeMatching == true || ringSize == 4) {
        isAxialPair = false; // init
        isAxialPair =
            ring::discardExtraTetragonBlocks(&basal1, &basal2, yCloud);
        if (isAxialPair == false) {
          continue;
        }
      } // end of check for tetragonal prism blocks
      // ------------
      // Step one: Check to see if basal1 and basal2 have common
      // elements or not. If they don't, then they cannot be basal rings
      cond1 = ring::hasCommonElements(basal1, basal2);
      if (cond1 == true) {
        continue;
      }
      // -----------
      // Step two and three: One of the elements of basal2 must be the nearest
      // neighbour of the first (index0; l1) If m_k is the nearest neighbour of
      // l1, m_{k+1} ... m_{k+(n-1)} must be neighbours of l_i+1 etc or l_i-1
      cond2 = ring::basalPrismConditions(nList, &basal1, &basal2);
      // If cond2 is false, the strict criteria for prisms has not been met
      if (cond2 == false) {
        // Skip if shape-matching is not desired
        if (!doShapeMatching) {
          continue;
        } // shape-matching not desired
        //
        // If shape-matching is to be done:
        // Check for the reduced criteria fulfilment
        relaxedCond = ring::relaxedPrismConditions(nList, &basal1, &basal2);
        // Skip if relaxed criteria are not met
        if (relaxedCond == false) {
          continue;
        } // end of skipping if the prisms do not fulfil relaxed criteria

        // Do shape matching here
        bool isDeformedPrism = match::matchPrism(
            yCloud, nList, refPointSet, &basal1, &basal2, rmsdPerAtom, false);

        // Success! The rings are basal rings of a deformed prism!
        if (isDeformedPrism) {
          //
          // // Update the number of prism blocks
          *nImperfectPrisms += 1;
          // Update iring
          if ((*ringType)[iring] == ring::strucType::unclassified) {
            (*ringType)[iring] = ring::strucType::deformedPrism;
            listPrism.push_back(iring);
          } else if ((*ringType)[iring] == ring::strucType::Prism) {
            (*ringType)[iring] = ring::strucType::mixedPrismRing;
          } // if it is deformed
          // Update jring
          if ((*ringType)[jring] == ring::strucType::unclassified) {
            (*ringType)[jring] = ring::strucType::deformedPrism;
            listPrism.push_back(jring);
          } else if ((*ringType)[jring] == ring::strucType::Prism) {
            (*ringType)[jring] = ring::strucType::mixedPrismRing;
          } // if it is deformed
        }   // end of update of ring types

        // // Write outs
        // // Now write out axial basal rings
        // sout::writeBasalRingsPrism(&basal1, &basal2, nDeformedPrisms, nList,
        //                            yCloud, true);
      } // end of reduced criteria
      // Strict criteria
      else {
        // Update the number of prism blocks
        *nPerfectPrisms += 1;
        // Update iring
        if ((*ringType)[iring] == ring::strucType::unclassified) {
          (*ringType)[iring] = ring::strucType::Prism;
          listPrism.push_back(iring);
        } else if ((*ringType)[iring] == ring::strucType::deformedPrism) {
          (*ringType)[iring] = ring::strucType::mixedPrismRing;
        } // if it is deformed
        // Update jring
        if ((*ringType)[jring] == ring::strucType::unclassified) {
          (*ringType)[jring] = ring::strucType::Prism;
          listPrism.push_back(jring);
        } else if ((*ringType)[jring] == ring::strucType::deformedPrism) {
          (*ringType)[jring] = ring::strucType::mixedPrismRing;
        } // if it is deformed
        //
        // Shape-matching to get the RMSD (if shape-matching is desired)
        if (doShapeMatching) {
          bool isKnownPrism = match::matchPrism(
              yCloud, nList, refPointSet, &basal1, &basal2, rmsdPerAtom, true);
        } // end of shape-matching to get rmsd
        //
        // // Now write out axial basal rings for convex hull calculations
        // sout::writePrisms(&basal1, &basal2, *nPrisms, yCloud);
        // // Write out prisms for shape-matching
        // sout::writeBasalRingsPrism(&basal1, &basal2, *nPrisms, nList, yCloud,
        //                            false);
        // -----------
      } // end of strict criteria

    } // end of loop through rest of the rings to get the second basal ring
  }   // end of loop through all rings for first basal ring

  sort(listPrism.begin(), listPrism.end());
  auto ip = std::unique(listPrism.begin(), listPrism.end());
  // Resize peripheral rings to remove undefined terms
  listPrism.resize(std::distance(listPrism.begin(), ip));

  return listPrism;
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
bool ring::basalPrismConditions(std::vector<std::vector<int>> nList,
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
 * @details Relaxed criteria for deformed prism blocks: at least one bond should
 * exist between the basal rings.
 */
bool ring::relaxedPrismConditions(std::vector<std::vector<int>> nList,
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

/**
 * @details A function that discards basal tetragonal rings which are not
 * oriented perpendicular to the axial direction. This is will only work for the
 * XY, YZ or XZ planes. If true is returned, then the rings are axial. This
 * function is only needed for tetragonal prism blocks because in the other
 * cases, the basal rings and lateral planes have a different number of nodes.
 * @param[in] basal1 The first basal ring.
 * @param[in] basal2 The other candidate basal ring.
 * @param[in] yCloud The input PointCloud.
 * @return A bool value which is true if the candidate basal rings are axial,
 *  and is false if the planes are lateral planes.
 */
bool ring::discardExtraTetragonBlocks(
    std::vector<int> *basal1, std::vector<int> *basal2,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  int ringSize =
      (*basal1).size(); // Size of the ring; each ring contains n elements
  int iatomIndex,
      jatomIndex;  // Indices of the elements in basal1 and basal2 respectively
  double r_i, r_j; // Coordinates in the axial dimension of iatom and jatom of
                   // basal1 and basal2 respectively
  int axialDim;    // 0 for x, 1 for y and 2 for z dimensions respectively
  // Variables for getting the projected area
  bool axialBasal1, axialBasal2; // bools for checking if basal1 and basal2 are
                                 // axial (true) respectively
  double areaXY, areaXZ,
      areaYZ; // Projected area on the XY, XZ and YZ planes respectively
  // ----------------------------------------
  // Find the axial dimension for a quasi-one-dimensional ice nanotube
  // The axial dimension will have the largest box length
  // Index -> axial dimension
  // 0 -> x dim
  // 1 -> y dim
  // 2 -> z dim
  axialDim = std::max_element(yCloud->box.begin(), yCloud->box.end()) -
             yCloud->box.begin();
  // ----------------------------------------
  // Calculate projected area onto the XY, YZ and XZ planes for basal1
  axialBasal1 = false; // Init to false
  axialBasal2 = false; // Init

  // Init the projected area
  areaXY = 0.0;
  areaXZ = 0.0;
  areaYZ = 0.0;

  jatomIndex = (*basal1)[0];

  // All points except the first pair
  for (int k = 1; k < ringSize; k++) {
    iatomIndex = (*basal1)[k]; // Current vertex

    // Add to the polygon area
    // ------
    // XY plane
    areaXY += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
              (yCloud->pts[jatomIndex].y - yCloud->pts[iatomIndex].y);
    // ------
    // XZ plane
    areaXZ += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
              (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
    // ------
    // YZ plane
    areaYZ += (yCloud->pts[jatomIndex].y + yCloud->pts[iatomIndex].y) *
              (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
    // ------
    jatomIndex = iatomIndex;
  }

  // Closure point
  iatomIndex = (*basal1)[0];
  // ------
  // XY plane
  areaXY += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
            (yCloud->pts[jatomIndex].y - yCloud->pts[iatomIndex].y);
  // ------
  // XZ plane
  areaXZ += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
            (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
  // ------
  // YZ plane
  areaYZ += (yCloud->pts[jatomIndex].y + yCloud->pts[iatomIndex].y) *
            (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
  // ------
  // The actual projected area is half of this
  areaXY *= 0.5;
  areaXZ *= 0.5;
  areaYZ *= 0.5;
  // Get the absolute value
  areaXY = fabs(areaXY);
  areaXZ = fabs(areaXZ);
  areaYZ = fabs(areaYZ);

  // If the axial dimension is x, y, or z:
  // then the maximum basal area should be in the YZ, XZ and XY dimensions
  // respectively
  // x dim
  if (axialDim == 0) {
    if (areaYZ > areaXY && areaYZ > areaXZ) {
      axialBasal1 = true;
    } // end of check for axial ring for basal1
  }   // x dim
  // y dim
  else if (axialDim == 1) {
    if (areaXZ > areaXY && areaXZ > areaYZ) {
      axialBasal1 = true;
    } // end of check for axial ring for basal1
  }   // x dim
  // z dim
  else if (axialDim == 2) {
    if (areaXY > areaXZ && areaXY > areaYZ) {
      axialBasal1 = true;
    } // end of check for axial ring for basal1
  }   // x dim
  else {
    std::cerr << "Could not find the axial dimension.\n";
    return false;
  }
  // ----------------------------------------
  // Calculate projected area onto the XY, YZ and XZ planes for basal2

  // Init the projected area
  areaXY = 0.0;
  areaXZ = 0.0;
  areaYZ = 0.0;

  jatomIndex = (*basal2)[0];

  // All points except the first pair
  for (int k = 1; k < ringSize; k++) {
    iatomIndex = (*basal2)[k]; // Current vertex

    // Add to the polygon area
    // ------
    // XY plane
    areaXY += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
              (yCloud->pts[jatomIndex].y - yCloud->pts[iatomIndex].y);
    // ------
    // XZ plane
    areaXZ += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
              (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
    // ------
    // YZ plane
    areaYZ += (yCloud->pts[jatomIndex].y + yCloud->pts[iatomIndex].y) *
              (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
    // ------
    jatomIndex = iatomIndex;
  }

  // Closure point
  iatomIndex = (*basal2)[0];
  // ------
  // XY plane
  areaXY += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
            (yCloud->pts[jatomIndex].y - yCloud->pts[iatomIndex].y);
  // ------
  // XZ plane
  areaXZ += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
            (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
  // ------
  // YZ plane
  areaYZ += (yCloud->pts[jatomIndex].y + yCloud->pts[iatomIndex].y) *
            (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
  // ------
  // The actual projected area is half of this
  areaXY *= 0.5;
  areaXZ *= 0.5;
  areaYZ *= 0.5;
  // Get the absolute value
  areaXY = fabs(areaXY);
  areaXZ = fabs(areaXZ);
  areaYZ = fabs(areaYZ);

  // Check if xy projected area is the greatest
  // If the axial dimension is x, y, or z:
  // then the maximum basal area should be in the YZ, XZ and XY dimensions
  // respectively
  // x dim
  if (axialDim == 0) {
    if (areaYZ > areaXY && areaYZ > areaXZ) {
      axialBasal2 = true;
    } // end of check for axial ring for basal1
  }   // x dim
  // y dim
  else if (axialDim == 1) {
    if (areaXZ > areaXY && areaXZ > areaYZ) {
      axialBasal2 = true;
    } // end of check for axial ring for basal1
  }   // x dim
  // z dim
  else if (axialDim == 2) {
    if (areaXY > areaXZ && areaXY > areaYZ) {
      axialBasal2 = true;
    } // end of check for axial ring for basal1
  }   // x dim
  else {
    std::cerr << "Could not find the axial dimension.\n";
    return false;
  }
  // ----------------------------------------

  // Now check if basal1 and basal2 are axial or not
  if (axialBasal1 == true && axialBasal2 == true) {
    return true;
  } else {
    return false;
  } // Check for basal1 and basal2
}

/**
 * @details Assign an atomType (equal to the number of nodes in the ring)
 * given a vector with a list of indices of rings comprising the prisms.
 * Note that the ring indices in the listPrism vector correspond to the indices
 * in the input rings vector of vectors.
 * @param[in] rings The vector of vectors containing the primitive rings, of a
 *  particular ring size.
 * @param[in] listPrism The list of prism blocks found.
 * @param[in] ringSize The current ring size or number of nodes in each ring.
 * @param[in, out] atomTypes A vector which contains a type for each atom,
 *  depending on it's type as classified by the prism identification scheme.
 */
int ring::assignPrismType(std::vector<std::vector<int>> rings,
                          std::vector<int> listPrism, int ringSize,
                          std::vector<ring::strucType> ringType,
                          std::vector<int> *atomTypes,
                          std::vector<ring::strucType> *atomState) {
  // Every value in listPrism corresponds to an index in rings.
  // Every ring contains atom indices, corresponding to the indices (not atom
  // IDs) in rings
  int iring;                    // Index of current ring
  int iatom;                    // Index of current atom
  ring::strucType currentState; // Current state of the ring being filled

  // Dummy value corresponds to a value of 1.
  // Each value is initialized to the value of 1.

  // Loop through every ring in rings
  for (int i = 0; i < listPrism.size(); i++) {
    iring = listPrism[i];
    // Get the state of all atoms in the ring
    currentState = ringType[iring];
    // Loop through every element in iring
    for (int j = 0; j < ringSize; j++) {
      iatom = rings[iring][j]; // Atom index
      // Update the atom type
      (*atomTypes)[iatom] = ringSize;
      //
      // Update the state of the atom
      if ((*atomState)[iatom] == ring::strucType::unclassified) {
        (*atomState)[iatom] = currentState;
      } // Update the unclassified atom
      else {
        if ((*atomState)[iatom] != currentState) {
          // For mixed, there is a preference
          if (currentState == ring::strucType::mixedPrismRing) {
            (*atomState)[iatom] = currentState;
          } // fill
          else if ((*atomState)[iatom] == ring::strucType::deformedPrism &&
                   currentState == ring::strucType::Prism) {
            (*atomState)[iatom] = ring::strucType::mixedPrismRing;
          } else if ((*atomState)[iatom] == ring::strucType::Prism &&
                     currentState == ring::strucType::deformedPrism) {
            (*atomState)[iatom] = ring::strucType::mixedPrismRing;
          }
        } //
      }   // already filled?

    } // end of loop through every atom in iring
  }   // end of loop through every ring

  return 0;
} // end of function

// Get the atom type values for deformed prisms
/**
 * @details Assign an atomType value
 * for atoms belonging to deformed prisms.
 * @param[in] rings The vector of vectors containing the primitive rings, of a
 *  particular ring size.
 * @param[in] listPrism The list of prism blocks found.
 * @param[in] ringSize The current ring size or number of nodes in each ring.
 * @param[in, out] atomTypes A vector which contains a type for each atom,
 * depending on it's type as classified by the prism identification scheme.
 */
int ring::deformedPrismTypes(std::vector<ring::strucType> atomState,
                             std::vector<int> *atomTypes, int maxDepth) {
  //
  int nop = atomState.size(); // Number of particles

  // Loop through all particles
  for (int iatom = 0; iatom < nop; iatom++) {
    // Check the atom state
    // Deformed
    if (atomState[iatom] == ring::strucType::deformedPrism) {
      (*atomTypes)[iatom] += maxDepth - 2;
    } // type for a deformed prism atom
    else if (atomState[iatom] == ring::strucType::mixedPrismRing) {
      (*atomTypes)[iatom] = 2;
    } // type for a mixed prism ring
  }   // end of reassignation

  //
  return 0;
}

// Shift the entire ice nanotube and remove axial translations
/**
 * @details Assign an atomType value
 * for atoms belonging to deformed prisms.
 * @param[in] rings The vector of vectors containing the primitive rings, of a
 *  particular ring size.
 * @param[in] listPrism The list of prism blocks found.
 * @param[in] ringSize The current ring size or number of nodes in each ring.
 * @param[in, out] atomTypes A vector which contains a type for each atom,
 *  depending on it's type as classified by the prism identification scheme.
 */
int ring::rmAxialTranslations(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int *atomID,
    int firstFrame, int currentFrame) {
  //
  int atomIndex;          // Index of the atom to be shifted first
  double shiftDistance;   // Value by which all atoms will be shifted
  double distFrmBoundary; // Distance from either the upper or lower boundary
                          // along the axial dimension
  // Get the axial dimension
  int axialDim = std::max_element(yCloud->box.begin(), yCloud->box.end()) -
                 yCloud->box.begin();
  // Check: (default z)
  if (axialDim < 0 || axialDim > 2) {
    axialDim = 2;
  } // end of check to set the axial dimension
  // Lower and higher limits of the box in the axial dimension
  double boxLowAxial = yCloud->boxLow[axialDim];
  double boxHighAxial = boxLowAxial + yCloud->box[axialDim];
  //
  // -----------------------------------
  // Update the atomID of the 'first' atom in yCloud if the current frame is the
  // first frame.
  // Get the atom index of the first atom to be shifted down or up
  if (currentFrame == firstFrame) {
    *atomID = yCloud->pts[0].atomID;
    atomIndex = 0;
  } // end of update of the atom ID to be shifted for the first frame
  else {
    //
    int iatomID = *atomID; // Atom ID of the first particle to be moved down
    // Find the index given the atom ID
    auto it = yCloud->idIndexMap.find(iatomID);

    if (it != yCloud->idIndexMap.end()) {
      atomIndex = it->second;
    } // found jatom
    else {
      std::cerr << "Lost atoms.\n";
      return 1;
    } // error handling
    //
  } // Update of atomIndex for all other frames
  // -----------------------------------
  // Calculate the distance by which the atoms must be shifted (negative value)
  if (axialDim == 0) {
    shiftDistance = boxLowAxial - yCloud->pts[atomIndex].x;
  } // x coordinate
  else if (axialDim == 1) {
    shiftDistance = boxLowAxial - yCloud->pts[atomIndex].y;
  } // y coordinate
  else {
    shiftDistance = boxLowAxial - yCloud->pts[atomIndex].z;
  } // z coordinate
  // -----------------------------------
  // Shift all the particles
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    // Shift the particles along the axial dimension only
    if (axialDim == 0) {
      yCloud->pts[iatom].x += shiftDistance;
      // Wrap-around
      if (yCloud->pts[iatom].x < boxLowAxial) {
        distFrmBoundary = boxLowAxial - yCloud->pts[iatom].x; // positive value
        yCloud->pts[iatom].x = boxHighAxial - distFrmBoundary;
      } // end of wrap-around
    }   // x coordinate
    else if (axialDim == 1) {
      yCloud->pts[iatom].y += shiftDistance;
      // Wrap-around
      if (yCloud->pts[iatom].y < boxLowAxial) {
        distFrmBoundary = boxLowAxial - yCloud->pts[iatom].y; // positive value
        yCloud->pts[iatom].y = boxHighAxial - distFrmBoundary;
      } // end of wrap-around
    }   // y coordinate
    else {
      yCloud->pts[iatom].z += shiftDistance;
      // Wrap-around
      if (yCloud->pts[iatom].z < boxLowAxial) {
        distFrmBoundary = boxLowAxial - yCloud->pts[iatom].z; // positive value
        yCloud->pts[iatom].z = boxHighAxial - distFrmBoundary;
      } // end of wrap-around
    }   // z coordinate
  }     // end of shifting all paritcles downward
  // -----------------------------------
  return 0;
} // end of function
