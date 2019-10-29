// Internal Libraries
#include <bond.hpp>
#include <generic.hpp>

/********************************************/ /**
 *  Create a vector of vectors containing bond information (bonded atom IDs, not
 indices!) from the ring vector of vectors
 ***********************************************/
std::vector<std::vector<int>>
bond::populateBonds(std::vector<std::vector<int>> rings) {
  std::vector<std::vector<int>> bonds; // Output vector of vectors
  std::vector<int> currentBond;        // Vector for the current bond
  int ringSize = rings[0].size();

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

  // Traverse the rings vector of vectors

  // Loop through all the rings
  for (int iring = 0; iring < rings.size(); iring++) {
    // Get the first atom of each pair inside iring
    for (int k = 0; k < rings[iring].size() - 1; k++) {
      // Clear the current bond vector
      currentBond.clear();
      // Fill the current bond vector
      currentBond.push_back(rings[iring][k]);
      currentBond.push_back(rings[iring][k + 1]);
      std::sort(currentBond.begin(), currentBond.end());
      // Add to the bonds vector of vectors
      bonds.push_back(currentBond);
    } // end of loop through ring elements, except the last one
    // The last pair is with the last element and the first element
    // Fill currentBond and update bonds
    currentBond.clear();
    currentBond.push_back(rings[iring][ringSize - 1]);
    currentBond.push_back(rings[iring][0]);
    std::sort(currentBond.begin(), currentBond.end());
    bonds.push_back(currentBond);
  } // end of loop through rings

  return bonds;
}

/********************************************/ /**
 *  Create a vector of vectors containing bond information (bonded atom IDs, not
 indices!) from the ring vector of vectors
 ***********************************************/
std::vector<std::vector<int>>
bond::populateHbonds(std::string filename,
                     molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                     std::vector<std::vector<int>> nList, int targetFrame,
                     int Htype) {
  //
  std::vector<std::vector<int>>
      hBondNet; // Output vector of vectors containing the HBN
  molSys::PointCloud<molSys::Point<double>, double>
      hCloud; // point Cloud for the hydrogen atoms
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
  double distCutoff = 2.42; // Distance cutoff of O-H, hard-coded
  double angleCutoff = 30;  // Angle cutoff in degrees
  std::vector<int> ooVec;   // Array for the O--O vector
  std::vector<int> ohVec;   // Array for the O-H vector

  // pi value
  auto pi = std::atan(1) * 4;

  // --------------------
  // Get all the hydrogen atoms in the frame (no slice)
  hCloud = sinp::readLammpsTrjreduced(filename, targetFrame, &hCloud, Htype);

  // Get the unordered map of the oxygen atom IDs (keys) and the molecular IDs
  // (values)
  idMolIDmap = molSys::createIDMolIDmap(yCloud);

  // Get a vector of vectors with the molID in the first column, and the
  // hydrogen atom indices (not ID) in each row
  molIDlist = molSys::hAtomMolList(&hCloud, yCloud);

  // Loop through the neighbour list
  for (int iatom = 0; iatom < nList.size(); iatom++) {
    currentBondList.clear();   // Clear the current bond vector
    iatomID = nList[iatom][0]; // atom ID corresponding to iatom
    currentBondList.push_back(
        iatomID); // The first element should be the atomID of iatom
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
        double dot_product =
            std::inner_product(ooVec.begin(), ooVec.end(), ohVec.begin(), 0);
        double normFactor = sqrt(std::inner_product(ooVec.begin(), ooVec.end(),
                                                    ooVec.begin(), 0)) *
                            sqrt(std::inner_product(ohVec.begin(), ohVec.end(),
                                                    ohVec.begin(), 0));
        double calcAngle = acos(dot_product / normFactor); // in radians
        calcAngle *= 180 / pi;                             // Convert to degrees
        //
        // A hydrogen bond is formed if the angle is less than 30 degrees
        if (calcAngle > angleCutoff) {
          continue;
        } // not a hydrogen bond

        // If you have reached this point, then O and H and indeed
        // hydrogen-bonded. This means that jatomID should be saved in the new
        // currentBond
        currentBondList.push_back(jatomID);
        break; // No need to test the other hydrogen atom if the first has been
               // tested
      }        // end of loop through hydrogen atoms

    } // end of loop through the nearest neighbours

    // Update HBN vector of vectors with currentBondList
    // hBondNet.push_back(std::vector<int>());  // Empty vector init
    hBondNet.push_back(currentBondList);
  } // end of loop through the neighbour list

  // --------------------

  // Erase all temporary stuff
  hCloud = molSys::clearPointCloud(&hCloud);

  return hBondNet;
}

/********************************************/ /**
                                                *  Calculates the bond length
                                                *between a Hydrogen and Oxygen
                                                *atom of two different atoms,
                                                *given their respective
                                                *pointClouds and the indices to
                                                *each atom
                                                ***********************************************/
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

/********************************************/ /**
 *  Create a vector of vectors containing bond information (bonded atom IDs, not
 vector or array indices!) from the ring vector of vectors and cageList
 ***********************************************/
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

/********************************************/ /**
 *  The bond 1 2 and 2 1 are the same. To prevent multiple bonds between the
 same atoms, remove all bonds which are duplicates of the reversed vectors
 (denoting individual bonds) within the bonds vector of vectors
 ***********************************************/
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

/********************************************/ /**
 *  Remove bonds which are diagonal, cutting across rings
 Set the flag to false if the bond is diagonal
 ***********************************************/
int bond::rmDiagBonds(std::vector<std::vector<int>> rings,
                      std::vector<std::vector<int>> bonds,
                      std::vector<bool> *flag) {
  std::vector<int> currentBond;   // Bond to be compared. A diagonal bond
  std::vector<int> reversedBond;  // Reversed order wrt currentBond
  int ringSize = rings[0].size(); // Number of elements in each ring
  int nNonBondPairs =
      (ringSize - 1) - 2; // Number of non-bonding pairs per atom in each ring
  int iatom, jatom;       // Actual atom IDs within a bond
  int j, k;               // Counter

  if (nNonBondPairs == 0) {
    return 0;
  }

  for (int iring = 0; iring < rings.size(); iring++) {
    // Select all but the last atom. Use reversedBond for reversing the order
    for (int i = 0; i < ringSize - 1; i++) {
      iatom = rings[iring][i];
      // Loop over other elements of the ring wrt i
      j = i + 2;
      for (int nbonds = 0; nbonds < nNonBondPairs; nbonds++) {
        k = j + nbonds; // Index
        // Wrap the index around the ring
        if (k >= ringSize) {
          k -= ringSize;
        }
        jatom = rings[iring][k];

        // Update currentBond and reversedBond
        currentBond.clear();
        reversedBond.clear();
        currentBond.push_back(iatom);
        currentBond.push_back(jatom);
        reversedBond.push_back(jatom);
        reversedBond.push_back(iatom);

        // Match with all other bonds
        bond::searchBondMatch(currentBond, bonds, flag);
        bond::searchBondMatch(reversedBond, bonds, flag);
      } // end of selection of the seconf atom in the bond
    }   // Selection of the first atom in the bond
  }     // end of loop through every ring

  // Bond selection requires the selection of a pair of atoms

  return 0;
}

/********************************************/ /**
 *  Search through vector of vector (bonds). If the input vector (matchBond)
 is the same as a bond inside the vector of vectors, set its flag to false
 ***********************************************/
int bond::searchBondMatch(std::vector<int> matchBond,
                          std::vector<std::vector<int>> bonds,
                          std::vector<bool> *flag) {
  for (int ibond; ibond < bonds.size(); ibond++) {
    // Skip a bond if the flag is already false
    if ((*flag)[ibond] == false) {
      continue;
    }

    // Check if the bond matches
    if (matchBond == bonds[ibond]) {
      (*flag)[ibond] = false;
    } // end of check for matching
  }   // end of loop through all possible bonds in the vector of vector

  return 0;
}

/********************************************/ /**
                                                *  Remove bonds which are longer
                                                *than the cutoff
                                                ***********************************************/
int bond::rmLongBonds(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                      std::vector<std::vector<int>> bonds,
                      std::vector<bool> *flag, double cutoff) {
  double bondLength; // Bond distance
  int iatom, jatom;  // Indices are 1 less than the atom ID
  for (int ibond; ibond < bonds.size(); ibond++) {
    // Skip a bond if the flag is already false
    if ((*flag)[ibond] == false) {
      continue;
    }

    iatom = bonds[ibond][0] - 1;
    jatom = bonds[ibond][1] - 1;
    // Calculate the bond length
    bondLength = gen::distance(yCloud, iatom, jatom);
    // Check if the bond matches
    if (bondLength > cutoff) {
      (*flag)[ibond] = false;
    } // end of check for matching
  }   // end of loop through all possible bonds in the vector of vector

  return 0;
}
