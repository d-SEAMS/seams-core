// Internal Libraries
#include <bond.hpp>
#include <generic.hpp>

/********************************************/ /**
 *  Create a vector of vectors containing bond information (bonded atom IDs, not
 indices!) from the ring vector of vectors
 ***********************************************/
std::vector<std::vector<int>> bond::populateBonds(
    std::vector<std::vector<int>> rings) {
  std::vector<std::vector<int>> bonds;  // Output vector of vectors
  std::vector<int> currentBond;         // Vector for the current bond
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
    }  // end of loop through ring elements, except the last one
    // The last pair is with the last element and the first element
    // Fill currentBond and update bonds
    currentBond.clear();
    currentBond.push_back(rings[iring][ringSize - 1]);
    currentBond.push_back(rings[iring][0]);
    std::sort(currentBond.begin(), currentBond.end());
    bonds.push_back(currentBond);
  }  // end of loop through rings

  return bonds;
}

/********************************************/ /**
 *  Create a vector of vectors containing bond information (bonded atom IDs, not
 indices!) from the ring vector of vectors
 ***********************************************/
std::vector<std::vector<int>> bond::populateHbonds(
    std::string filename,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, int targetFrame, int Htype) {
  std::vector<std::vector<int>>
      hBondNet;  // Output vector of vectors containing the HBN
  molSys::PointCloud<molSys::Point<double>, double>
      hCloud;  // point Cloud for the hydrogen atoms
  std::vector<std::vector<int>>
      molIDlist;  // Vector of vectors; first element is the molID, and the next
                  // two elements are the hydrogen atom indices
  std::unordered_map<int, int>
      idMolIDmap;  // Unordered map with atom IDs of oxygens as the keys and the
                   // molecular IDs as the values

  // --------------------
  // Get all the hydrogen atoms in the frame (no slice)
  hCloud = sinp::readLammpsTrjreduced(filename, targetFrame, &hCloud, Htype);

  // Get the unordered map of the oxygen atom IDs (keys) and the molecular IDs
  // (values)
  idMolIDmap = molSys::createIDMolIDmap(yCloud);

  // --------------------

  // Erase all temporary stuff
  hCloud = molSys::clearPointCloud(&hCloud);

  return hBondNet;
}

/********************************************/ /**
 *  Create a vector of vectors containing bond information (bonded atom IDs, not
 vector or array indices!) from the ring vector of vectors and cageList
 ***********************************************/
std::vector<std::vector<int>> bond::createBondsFromCages(
    std::vector<std::vector<int>> rings, std::vector<cage::Cage> *cageList,
    cage::cageType type, int *nRings) {
  std::vector<std::vector<int>> bonds;  // Output vector of vectors
  std::vector<int> currentBond;         // Vector for the current bond
  int ringSize = rings[0].size();
  int currentRing;  // (vector) index of the current ring in a particular cage

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

  *nRings = 0;  // init

  // Loop through all the cages
  for (int icage = 0; icage < (*cageList).size(); icage++) {
    // Skip if the cage is of a different type
    if ((*cageList)[icage].type != type) {
      continue;
    }
    *nRings += (*cageList)[icage].rings.size();  // Add to the number of rings
    //
    // Now loop through a particular ring inside the i^th cage
    for (int iring = 0; iring < (*cageList)[icage].rings.size(); iring++) {
      currentRing = (*cageList)[icage].rings[iring];  // Current ring index
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
      }  // end of loop through ring elements, except the last one
      // The last pair is with the last element and the first element
      // Fill currentBond and update bonds
      currentBond.clear();
      currentBond.push_back(rings[currentRing][ringSize - 1]);
      currentBond.push_back(rings[currentRing][0]);
      std::sort(currentBond.begin(), currentBond.end());
      bonds.push_back(currentBond);
    }  // end of loop through a particular ring
  }    // end of loop through cages

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
std::vector<std::vector<int>> bond::trimBonds(
    std::vector<std::vector<int>> bonds) {
  std::vector<int>
      reversedBond;  // Vector for the current bond in reversed order
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
  std::vector<int> currentBond;    // Bond to be compared. A diagonal bond
  std::vector<int> reversedBond;   // Reversed order wrt currentBond
  int ringSize = rings[0].size();  // Number of elements in each ring
  int nNonBondPairs =
      (ringSize - 1) - 2;  // Number of non-bonding pairs per atom in each ring
  int iatom, jatom;        // Actual atom IDs within a bond
  int j, k;                // Counter

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
        k = j + nbonds;  // Index
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
      }  // end of selection of the seconf atom in the bond
    }    // Selection of the first atom in the bond
  }      // end of loop through every ring

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
    }  // end of check for matching
  }    // end of loop through all possible bonds in the vector of vector

  return 0;
}

/********************************************/ /**
                                                *  Remove bonds which are longer
                                                *than the cutoff
                                                ***********************************************/
int bond::rmLongBonds(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                      std::vector<std::vector<int>> bonds,
                      std::vector<bool> *flag, double cutoff) {
  double bondLength;  // Bond distance
  int iatom, jatom;   // Indices are 1 less than the atom ID
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
    }  // end of check for matching
  }    // end of loop through all possible bonds in the vector of vector

  return 0;
}
