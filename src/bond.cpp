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

// Internal Libraries
#include <bond.hpp>
#include <generic.hpp>

/**
 * @details Create a vector of vectors containing bond information (outputs
 * bonded atom IDs, not indices!) from the neighbour list vector of vectors
 * (which contains atom INDICES). Moreover, the first element corresponds to the
 * atom whose neighbours have been found.
 */
std::vector<std::vector<int>>
bond::populateBonds(std::vector<std::vector<int>> nList,
                    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  //
  std::vector<std::vector<int>> bonds; // Output vector of vectors
  std::vector<int> currentBond;        // Vector for the current bond
  int first, second;                   // Indices for the bonds
  int iatom, jatom;                    // Elements of the bond
  int iatomID, jatomID; // Atom IDs of the elements which are bonded

  // Error handling
  if (nList.size() == 0) {
    // There is some problem!
    std::cerr << "There are no bonds in the system!\n";
    return bonds;
  }

  // Form of the bonds vector of vectors:
  // 214    272
  // 1       2
  // Meaning that 272 and 214 are bonded; similarly 1 and 2 are bonded
  // To avoid repitition, discard bonds whose first value is greater than the
  // second value

  // Traverse the neighbour list

  // Loop through every atom in the neighbour list by index
  for (int i = 0; i < nList.size(); i++) {
    iatom = nList[i][0]; // Index of the i^th atom
    // Get the neighbours of iatom
    for (int j = 1; j < nList[i].size(); j++) {
      //
      jatom = nList[iatom][j]; // Index of the neighbour
      // To avoid duplicates, skip all bonds such
      // that iatom>jatom
      if (iatom > jatom) {
        continue;
      } // Skip to avoid duplicates

      // Clear the current bond vector
      currentBond.clear();
      // Fill the current bond vector with ATOM IDs
      iatomID = yCloud->pts[iatom].atomID;
      jatomID = yCloud->pts[jatom].atomID;
      currentBond.push_back(iatomID);
      currentBond.push_back(jatomID);
      // Add to the bonds vector of vectors
      bonds.push_back(currentBond);
    } // end of loop the neighbour list
    // The last pair is with the last element and the first element
    // Fill currentBond and update bonds

  } // end of loop through rings

  return bonds;
}

/**
 *  @details Create a vector of vectors containing bond information (outputs
 * bonded atom IDs, not indices!) from the neighbour list vector of vectors
 * (which contains atom INDICES). Bonds between dummy atoms, and between dummy
 * and ice atoms are not added. Moreover, the first element corresponds to the
 * atom whose neighbours have been found.
 *  @param[in] nList Row-ordered neighbour list by ID
 *  @param[in] yCloud The input molSys::PointCloud
 *  @param[in] atomTypes Contains an atom type for each particle in yCloud
 */
std::vector<std::vector<int>>
bond::populateBonds(std::vector<std::vector<int>> nList,
                    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::vector<cage::iceType> atomTypes) {
  //
  std::vector<std::vector<int>> bonds; // Output vector of vectors
  std::vector<int> currentBond;        // Vector for the current bond
  int first, second;                   // Indices for the bonds
  int iatom, jatom;                    // Elements of the bond
  int iatomID, jatomID; // Atom IDs of the elements which are bonded

  // Error handling
  if (nList.size() == 0) {
    // There is some problem!
    std::cerr << "There are no bonds in the system!\n";
    return bonds;
  }

  // Form of the bonds vector of vectors:
  // 214    272
  // 1       2
  // Meaning that 272 and 214 are bonded; similarly 1 and 2 are bonded
  // To avoid repitition, discard bonds whose first value is greater than the
  // second value

  // Traverse the neighbour list

  // Loop through every atom in the neighbour list by index
  for (int i = 0; i < nList.size(); i++) {
    iatom = nList[i][0]; // Index of the i^th atom
    // Skip for dummy atoms
    if (atomTypes[iatom] == cage::iceType::dummy) {
      continue;
    } // Skip for dummy atoms
    // Get the neighbours of iatom
    for (int j = 1; j < nList[i].size(); j++) {
      //
      jatom = nList[iatom][j]; // Index of the neighbour
      // Skip for dummy atoms
      if (atomTypes[jatom] == cage::iceType::dummy) {
        continue;
      } // Skip for dummy atoms
      // To avoid duplicates, skip all bonds such
      // that iatom>jatom
      if (iatom > jatom) {
        continue;
      } // Skip to avoid duplicates

      // Clear the current bond vector
      currentBond.clear();
      // Fill the current bond vector with ATOM IDs
      iatomID = yCloud->pts[iatom].atomID;
      jatomID = yCloud->pts[jatom].atomID;
      currentBond.push_back(iatomID);
      currentBond.push_back(jatomID);
      // Add to the bonds vector of vectors
      bonds.push_back(currentBond);
    } // end of loop the neighbour list
    // The last pair is with the last element and the first element
    // Fill currentBond and update bonds

  } // end of loop through rings

  return bonds;
}

/********************************************/ /**
 *  Create a vector of vectors (similar to the neighbour list conventions). The
 output vector of vectors is row-ordered. The first element is the atom ID of
 the particle for which the neighbours are enumerated (the central atom),
 followed by the central atom's neighbour's IDs (not indices). Decides the
 existence of the hydrogen bond depending on the O--O and O--H vectors from the
 neighbour list (by ID) already constructed.
 *  @param[in] filename Filename of the trajectory, with the hydrogen and oxygen
 coordinates
 *  @param[in] yCloud The input molSys::PointCloud for the oxygen atoms only
 *  @param[in] nList Row-ordered neighbour list by atom ID
 *  @param[in] targetFrame The target or current frame number (starts from 1)
 and is not the timestep value
 *  @param[in] Htype The type ID of the hydrogen atoms
 ***********************************************/
std::vector<std::vector<int>>
bond::populateHbonds(std::string filename,
                     molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                     std::vector<std::vector<int>> nList, int targetFrame,
                     int Htype) {
  //
  std::vector<std::vector<int>>
      hBondNet; // Output vector of vectors containing the HBN
  std::vector<std::vector<int>>
      molIDlist; // Vector of vectors; first element is the molID, and the next
                 // two elements are the hydrogen atom indices
  std::unordered_map<int, int>
      idMolIDmap; // Unordered map with atom IDs of oxygens as the keys and the
                  // molecular IDs as the values
  std::vector<int> currentBondList; // Current bond list for atom
  int nnumNeighbours;   // Number of nearest neighbours for the current atom
  int iatomID, jatomID; // Atom IDs
  int iatomIndex, jatomIndex; // Atomic indices of oxygen atoms
  int hAtomIndex;             // Atom index of hydrogen
  int listIndex;   // Index in molIDlist corresponding to a particular molecular
                   // ID
  int jOxyMolID;   // Molecular ID of the jatom oxygen atom
  double hBondLen; // Length of O-H (between the donor O and acceptor H)
  double distCutoff = 2.42;  // Distance cutoff of O-H, hard-coded
  double angleCutoff = 30;   // Angle cutoff in degrees
  std::vector<double> ooVec; // Array for the O--O vector
  std::vector<double> ohVec; // Array for the O-H vector

  // --------------------
  // Get all the hydrogen atoms in the frame (no slice)
  auto hCloud = sinp::readLammpsTrjreduced(filename, targetFrame, Htype); // point Cloud for the hydrogen atoms

  // Get the unordered map of the oxygen atom IDs (keys) and the molecular IDs
  // (values)
  idMolIDmap = molSys::createIDMolIDmap(yCloud);

  // Get a vector of vectors with the molID in the first column, and the
  // hydrogen atom indices (not ID) in each row
  molIDlist = molSys::hAtomMolList(&hCloud, yCloud);

  // Initialize the vector of vectors hBondNet
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    hBondNet.push_back(std::vector<int>()); // Empty vector for the index iatom
    // Fill the first element with the atom ID of iatom itself
    hBondNet[iatom].push_back(yCloud->pts[iatom].atomID);
  } // end of init of hBondNet

  // Loop through the neighbour list
  for (int iatom = 0; iatom < nList.size(); iatom++) {
    currentBondList.clear();                  // Clear the current bond vector
    iatomID = nList[iatom][0];                // atom ID corresponding to iatom
    nnumNeighbours = nList[iatom].size() - 1; // No. of nearest neighbours
    iatomIndex = iatom;                       // Atomic index
    //
    // Loop through the nearest neighbours
    for (int j = 1; j <= nnumNeighbours; j++) {
      jatomID = nList[iatom][j]; // Atom ID of the nearest neighbour
      // Get the hydrogen atom indices corresponding to the molID of jatomID
      // Find jOxyMolID
      auto it = idMolIDmap.find(jatomID);
      if (it != idMolIDmap.end()) {
        jOxyMolID = it->second;
      } // found molecular ID of jatom oxygen atom
      else {
        continue;
      } // not found

      // Find the index inside the molIDlist corresponding to the molecular ID
      // to look for
      listIndex = molSys::searchMolList(molIDlist, jOxyMolID);

      // Get the atom index of the oxygen atom jatom corresponding jatomID
      auto gotJ = yCloud->idIndexMap.find(jatomID);
      if (gotJ != yCloud->idIndexMap.end()) {
        jatomIndex = gotJ->second;
      } // found atom index of jatomID
      else {
        std::cerr << "Something is wrong with the map.\n";
        continue;
      } // index not found

      // Loop through the hydrogen atoms connected to jatom oxygen atom
      for (int k = 1; k <= 2; k++) {
        hAtomIndex = molIDlist[listIndex][k];
        // --------
        // Condition One: The O-H length (between the donor hydrogen atom and
        // the acceptor oxygen atom) should be less than 2.42 Angstrom
        // (hard-coded)
        hBondLen =
            bond::getHbondDistanceOH(yCloud, &hCloud, iatomIndex, hAtomIndex);

        // If O-H distance is greater than or equal to 2.42 then it is not a
        // hydrogen bond
        if (hBondLen >= distCutoff) {
          continue;
        } // not a hydrogen bond
        // --------
        // Condition Two: The angle between the O-H and O--O vectors is less
        // than 30 degrees (hard-coded)
        //
        ooVec.clear();
        ohVec.clear();
        // Get the O--O and O-H vectors
        // O--O
        ooVec.push_back(yCloud->pts[iatomIndex].x -
                        yCloud->pts[jatomIndex].x); // dx
        ooVec.push_back(yCloud->pts[iatomIndex].y -
                        yCloud->pts[jatomIndex].y); // dy
        ooVec.push_back(yCloud->pts[iatomIndex].z -
                        yCloud->pts[jatomIndex].z); // dz
        // O-H
        ohVec.push_back(yCloud->pts[iatomIndex].x -
                        hCloud.pts[hAtomIndex].x); // dx
        ohVec.push_back(yCloud->pts[iatomIndex].y -
                        hCloud.pts[hAtomIndex].y); // dy
        ohVec.push_back(yCloud->pts[iatomIndex].z -
                        hCloud.pts[hAtomIndex].z); // dz
        // Apply PBCs
        for (int l = 0; l < 3; l++) {
          ooVec[l] -= yCloud->box[l] * round(ooVec[l] / yCloud->box[l]);
          ohVec[l] -= yCloud->box[l] * round(ohVec[l] / yCloud->box[l]);
        } // end of applying PBCs to the O-H and O--O vectors
        //
        // Get the angle between the O--O and O-H vectors
        double eigenAngle = gen::eigenVecAngle(ooVec, ohVec);
        double eigenAngleDeg = gen::radDeg(eigenAngle);

        //
        // A hydrogen bond is formed if the angle is less than 30 degrees
        if (eigenAngleDeg > angleCutoff) {
          continue;
        } // not a hydrogen bond

        // If you have reached this point, then O and H and indeed
        // hydrogen-bonded. This means that jatomID should be saved in the new
        // currentBond
        hBondNet[iatomIndex].push_back(jatomID);
        hBondNet[jatomIndex].push_back(iatomID);
        break; // No need to test the other hydrogen atom if the first has
        // been tested
      } // end of loop through hydrogen atoms

    } // end of loop through the nearest neighbours

  } // end of loop through the neighbour list

  // --------------------

  // Erase all temporary stuff
  hCloud = molSys::clearPointCloud(&hCloud);

  return hBondNet;
}

/********************************************/ /**
 *  Create a vector of vectors (similar to the neighbour list conventions). The
 output vector of vectors is row-ordered. The first element is the atom ID of
 the particle for which the neighbours are enumerated (the central atom),
 followed by the central atom's neighbour's IDs (not indices). Decides the
 existence of the hydrogen bond depending on the O--O and O--H vectors from the
 neighbour list (by ID) already constructed. The only difference between this function and the
 similar populateHbonds function is that this takes in the H atom PointCloud as input
 *  @param[in] yCloud The input molSys::PointCloud for the hydrogen atoms only,
 for the entire system (regardless of whether there is a slice or not)
 *  @param[in] hCloud The input molSys::PointCloud for the oxygen atoms only
 *  @param[in] nList Row-ordered neighbour list by atom ID
 ***********************************************/
std::vector<std::vector<int>>
bond::populateHbondsWithInputClouds(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                     molSys::PointCloud<molSys::Point<double>, double> *hCloud,
                     std::vector<std::vector<int>> nList) {
  //
  std::vector<std::vector<int>>
      hBondNet; // Output vector of vectors containing the HBN
  std::vector<std::vector<int>>
      molIDlist; // Vector of vectors; first element is the molID, and the next
                 // two elements are the hydrogen atom indices
  std::unordered_map<int, int>
      idMolIDmap; // Unordered map with atom IDs of oxygens as the keys and the
                  // molecular IDs as the values
  std::vector<int> currentBondList; // Current bond list for atom
  int nnumNeighbours;   // Number of nearest neighbours for the current atom
  int iatomID, jatomID; // Atom IDs
  int iatomIndex, jatomIndex; // Atomic indices of oxygen atoms
  int hAtomIndex;             // Atom index of hydrogen
  int listIndex;   // Index in molIDlist corresponding to a particular molecular
                   // ID
  int jOxyMolID;   // Molecular ID of the jatom oxygen atom
  double hBondLen; // Length of O-H (between the donor O and acceptor H)
  double distCutoff = 2.42;  // Distance cutoff of O-H, hard-coded
  double angleCutoff = 30;   // Angle cutoff in degrees
  std::vector<double> ooVec; // Array for the O--O vector
  std::vector<double> ohVec; // Array for the O-H vector

  // --------------------
  // Get the unordered map of the oxygen atom IDs (keys) and the molecular IDs
  // (values)
  idMolIDmap = molSys::createIDMolIDmap(yCloud);

  // Get a vector of vectors with the molID in the first column, and the
  // hydrogen atom indices (not ID) in each row
  molIDlist = molSys::hAtomMolList(hCloud, yCloud);

  // Initialize the vector of vectors hBondNet
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    hBondNet.push_back(std::vector<int>()); // Empty vector for the index iatom
    // Fill the first element with the atom ID of iatom itself
    hBondNet[iatom].push_back(yCloud->pts[iatom].atomID);
  } // end of init of hBondNet

  // Loop through the neighbour list
  for (int iatom = 0; iatom < nList.size(); iatom++) {
    currentBondList.clear();                  // Clear the current bond vector
    iatomID = nList[iatom][0];                // atom ID corresponding to iatom
    nnumNeighbours = nList[iatom].size() - 1; // No. of nearest neighbours
    iatomIndex = iatom;                       // Atomic index
    //
    // Loop through the nearest neighbours
    for (int j = 1; j <= nnumNeighbours; j++) {
      jatomID = nList[iatom][j]; // Atom ID of the nearest neighbour
      // Get the hydrogen atom indices corresponding to the molID of jatomID
      // Find jOxyMolID
      auto it = idMolIDmap.find(jatomID);
      if (it != idMolIDmap.end()) {
        jOxyMolID = it->second;
      } // found molecular ID of jatom oxygen atom
      else {
        continue;
      } // not found

      // Find the index inside the molIDlist corresponding to the molecular ID
      // to look for
      listIndex = molSys::searchMolList(molIDlist, jOxyMolID);

      // Get the atom index of the oxygen atom jatom corresponding jatomID
      auto gotJ = yCloud->idIndexMap.find(jatomID);
      if (gotJ != yCloud->idIndexMap.end()) {
        jatomIndex = gotJ->second;
      } // found atom index of jatomID
      else {
        std::cerr << "Something is wrong with the map.\n";
        continue;
      } // index not found

      // Loop through the hydrogen atoms connected to jatom oxygen atom
      for (int k = 1; k <= 2; k++) {
        hAtomIndex = molIDlist[listIndex][k];
        // --------
        // Condition One: The O-H length (between the donor hydrogen atom and
        // the acceptor oxygen atom) should be less than 2.42 Angstrom
        // (hard-coded)
        hBondLen =
            bond::getHbondDistanceOH(yCloud, hCloud, iatomIndex, hAtomIndex);

        // If O-H distance is greater than or equal to 2.42 then it is not a
        // hydrogen bond
        if (hBondLen >= distCutoff) {
          continue;
        } // not a hydrogen bond
        // --------
        // Condition Two: The angle between the O-H and O--O vectors is less
        // than 30 degrees (hard-coded)
        //
        ooVec.clear();
        ohVec.clear();
        // Get the O--O and O-H vectors
        // O--O
        ooVec.push_back(yCloud->pts[iatomIndex].x -
                        yCloud->pts[jatomIndex].x); // dx
        ooVec.push_back(yCloud->pts[iatomIndex].y -
                        yCloud->pts[jatomIndex].y); // dy
        ooVec.push_back(yCloud->pts[iatomIndex].z -
                        yCloud->pts[jatomIndex].z); // dz
        // O-H
        ohVec.push_back(yCloud->pts[iatomIndex].x -
                        hCloud->pts[hAtomIndex].x); // dx
        ohVec.push_back(yCloud->pts[iatomIndex].y -
                        hCloud->pts[hAtomIndex].y); // dy
        ohVec.push_back(yCloud->pts[iatomIndex].z -
                        hCloud->pts[hAtomIndex].z); // dz
        // Apply PBCs
        for (int l = 0; l < 3; l++) {
          ooVec[l] -= yCloud->box[l] * round(ooVec[l] / yCloud->box[l]);
          ohVec[l] -= yCloud->box[l] * round(ohVec[l] / yCloud->box[l]);
        } // end of applying PBCs to the O-H and O--O vectors
        //
        // Get the angle between the O--O and O-H vectors
        double eigenAngle = gen::eigenVecAngle(ooVec, ohVec);
        double eigenAngleDeg = gen::radDeg(eigenAngle);

        //
        // A hydrogen bond is formed if the angle is less than 30 degrees
        if (eigenAngleDeg > angleCutoff) {
          continue;
        } // not a hydrogen bond

        // If you have reached this point, then O and H and indeed
        // hydrogen-bonded. This means that jatomID should be saved in the new
        // currentBond
        hBondNet[iatomIndex].push_back(jatomID);
        hBondNet[jatomIndex].push_back(iatomID);
        break; // No need to test the other hydrogen atom if the first has
        // been tested
      } // end of loop through hydrogen atoms

    } // end of loop through the nearest neighbours

  } // end of loop through the neighbour list

  // --------------------

  return hBondNet;
}

/**
*  Calculates the bond length between a Hydrogen and Oxygen
 atom of two different atoms, given their respective pointClouds and the indices
to each atom.
 *  @param[in] oCloud The molSys::PointCloud for the oxygen atoms only
 *  @param[in] hCloud The molSys::PointCloud for the hydrogen atoms only
 *  @param[in] oAtomIndex The index (in the oCloud) of the oxygen atom
 *  @param[in] hAtomIndex The index (in the hCloud) of the hydrogen atom
*/
double bond::getHbondDistanceOH(
    molSys::PointCloud<molSys::Point<double>, double> *oCloud,
    molSys::PointCloud<molSys::Point<double>, double> *hCloud, int oAtomIndex,
    int hAtomIndex) {
  std::array<double, 3> dr; // relative distance in the X, Y, Z dimensions
  double r2 = 0.0;          // Bond length

  dr[0] = oCloud->pts[oAtomIndex].x - hCloud->pts[hAtomIndex].x;
  dr[1] = oCloud->pts[oAtomIndex].y - hCloud->pts[hAtomIndex].y;
  dr[2] = oCloud->pts[oAtomIndex].z - hCloud->pts[hAtomIndex].z;

  // Apply the PBCs and get the squared area
  for (int k = 0; k < 3; k++) {
    dr[k] -= oCloud->box[k] * round(dr[k] / oCloud->box[k]);
    r2 += pow(dr[k], 2.0);
  } // end of applying the PBCs and getting the squared area
  return sqrt(r2);
}

/**
 *  Create a vector of vectors containing bond information (bonded atom IDs, not
 vector or array indices!) from the ring vector of vectors and cageList
 *  @param[in] rings Row-ordered vector of vectors atom indices of ring
 information. Each row is a ring, containing the indices of the particles in
 that ring
 *  @param[in] cageList A vector of cage::Cage containing a list of HCs or DDCs
 *  @param[in] type The type of cage to get bonds for
 *  @param[in, out] nRings The total number of rings for all the cages, for the
 particular cage type
*/
std::vector<std::vector<int>>
bond::createBondsFromCages(std::vector<std::vector<int>> rings,
                           std::vector<cage::Cage> *cageList,
                           cage::cageType type, int *nRings) {
  std::vector<std::vector<int>> bonds; // Output vector of vectors
  std::vector<int> currentBond;        // Vector for the current bond
  int ringSize = rings[0].size();
  int currentRing; // (vector) index of the current ring in a particular cage

  // Error handling
  if (rings.size() == 0) {
    // There is some problem!
    std::cerr << "There are no rings in the system!\n";
    return bonds;
  }

  // Form of the bonds vector of vectors:
  // 272    214
  // 1       2
  // Meaning that 272 and 214 are bonded; similarly 1 and 2 are bonded

  // Traverse the cageList vector of Cages

  *nRings = 0; // init

  // Loop through all the cages
  for (int icage = 0; icage < (*cageList).size(); icage++) {
    // Skip if the cage is of a different type
    if ((*cageList)[icage].type != type) {
      continue;
    }
    *nRings += (*cageList)[icage].rings.size(); // Add to the number of rings
    //
    // Now loop through a particular ring inside the i^th cage
    for (int iring = 0; iring < (*cageList)[icage].rings.size(); iring++) {
      currentRing = (*cageList)[icage].rings[iring]; // Current ring index
      // Get the first atom of each pair inside currentRing
      for (int k = 0; k < rings[currentRing].size() - 1; k++) {
        // Clear the current bond vector
        currentBond.clear();
        // Fill the current bond vector
        currentBond.push_back(rings[currentRing][k]);
        currentBond.push_back(rings[currentRing][k + 1]);
        std::sort(currentBond.begin(), currentBond.end());
        // Add to the bonds vector of vectors
        bonds.push_back(currentBond);
      } // end of loop through ring elements, except the last one
      // The last pair is with the last element and the first element
      // Fill currentBond and update bonds
      currentBond.clear();
      currentBond.push_back(rings[currentRing][ringSize - 1]);
      currentBond.push_back(rings[currentRing][0]);
      std::sort(currentBond.begin(), currentBond.end());
      bonds.push_back(currentBond);
    } // end of loop through a particular ring
  }   // end of loop through cages

  // This may have duplicates, so the duplicates should be removed
  std::sort(bonds.begin(), bonds.end());
  bonds.erase(std::unique(bonds.begin(), bonds.end()), bonds.end());

  return bonds;
}

/**
 *  The bond 1 2 and 2 1 are the same. To prevent multiple bonds between the
 same atoms, remove all bonds which are duplicates of the reversed vectors
 (denoting individual bonds) within the bonds vector of vectors
 *  @param[in, out] bonds Row-ordered vector of vectors of the bond matrix
 Each row is a ring, containing the indices of the particles in that ring
*/
std::vector<std::vector<int>>
bond::trimBonds(std::vector<std::vector<int>> bonds) {
  std::vector<int>
      reversedBond; // Vector for the current bond in reversed order
  std::vector<bool> isBondFlag;
  int temp = 0;

  std::sort(bonds.begin(), bonds.end());
  bonds.erase(std::unique(bonds.begin(), bonds.end()), bonds.end());

  // // Temp uncommented
  // // Resize the flag vector and init to true
  // isBondFlag.resize(bonds.size(), true);
  // // Form of the bonds vector of vectors:
  // // 272    214
  // // 1       2
  // // Meaning that 272 and 214 are bonded; similarly 1 and 2 are bonded

  // // Traverse the bonds vector of vectors

  // // Loop through all the bonds
  // for(int ibond=0; ibond<bonds.size(); ibond++){

  //   // Skip if the bond has been flagged as false
  //   if(isBondFlag[ibond] == false){continue;}

  //   reversedBond.clear();
  //   reversedBond.push_back( bonds[ibond][1] );
  //   reversedBond.push_back( bonds[ibond][0] );

  //   // Compare with other bonds by looping through
  //   // the rest
  //   for(int jbond=0; jbond<bonds.size(); jbond++){
  //     // Skip if jbond has been flagged
  //     if(isBondFlag[jbond] == false){continue;}
  //     // Compare reversed bond and jbond if jbond has not been flagged yet
  //     if(reversedBond==bonds[jbond]){
  //       isBondFlag[jbond] = false;
  //       temp++;
  //     } // end of check of comparison
  //   } // end of looping through all possible bonds
  // } // end of loop through all possible bonds
  // // temp

  return bonds;
}
