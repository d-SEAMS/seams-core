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

#include <shapeMatch.hpp>

/**
 * @details Shape-matching for a pair of polygon basal rings. Returns true if
 * the pair of basal rings form a prism block.
 */
bool match::matchPrism(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, const Eigen::MatrixXd &refPoints,
    std::vector<int> *basal1, std::vector<int> *basal2,
    std::vector<double> *rmsdPerAtom, bool isPerfect) {
  //
  int ringSize = (*basal1).size(); // Number of nodes in each basal ring
  std::vector<int> matchedBasal1,
      matchedBasal2; // Re-ordered basal rings 1 and 2
  int dim = 3;       // Number of dimensions
  // Matrices for the point sets of the target basal rings
  Eigen::MatrixXd basal1Set(ringSize,
                            dim); // Point set for the first basal ring
  Eigen::MatrixXd basal2Set(ringSize,
                            dim); // Point set for the second basal ring
  int startingIndex; // Index in the basal rings from which the reference point
                     // set is matched
  double cutoffAngle = 15; // Tolerance for the angular distance between the
                           // quaternions of the basal rings
  double angDist; // Calculated angular distance between the two basal rings
  // Variables for the absolute orientation
  std::vector<double> quat1, quat2; // quaternion rotation
  double rmsd1, rmsd2;              // least total RMSD
  std::vector<double> rmsdList1,
      rmsdList2;         // List of RMSD per atom in the order fed in
  double scale1, scale2; // Scale factor
  // Deformed block classification
  bool doAngleCriterion = false; // For deformed blocks
  // -----------------------
  // Getting the target Eigen vectors
  // Get the re-ordered matched basal rings, ordered with respect to each other
  pntToPnt::relOrderPrismBlock(yCloud, *basal1, *basal2, nList, &matchedBasal1,
                               &matchedBasal2);
  // -----------------------
  // Match the basal rings with a complete prism block, given the relatively
  // ordered basal rings This actually only needs to be done for deformed prism
  // blocks
  bool blockMatch = match::matchPrismBlock(
      yCloud, nList, refPoints, &matchedBasal1, &matchedBasal2, &startingIndex);
  // -----------------------
  // Check for deformed prisms
  if (!isPerfect) {
    if (!blockMatch) {
      return blockMatch;
    }
  }
  // If the deformed prism does not fulfil the shape matching criterion, then do
  // not calculate the RMSD per atom
  // -----------------------
  // Section for RMSD per atom, based on basal ring matching

  basal1Set =
      pntToPnt::fillPointSetPrismRing(yCloud, matchedBasal1, startingIndex);
  // Fill up the point set for basal2
  basal2Set =
      pntToPnt::fillPointSetPrismRing(yCloud, matchedBasal2, startingIndex);
  // Use Horn's algorithm to calculate the absolute orientation and RMSD etc.
  absor::hornAbsOrientation(refPoints, basal1Set, &quat1, &rmsd1, &rmsdList1,
                            &scale1); // basal1
  absor::hornAbsOrientation(refPoints, basal2Set, &quat2, &rmsd2, &rmsdList2,
                            &scale2); // basal2
  // // ------------
  // // Update the per-atom RMSD for each particle
  // // Basal1
  // match::updatePerAtomRMSDRing(matchedBasal1, startingIndex, rmsdList1,
  //                              rmsdPerAtom);
  // // Basal2
  // match::updatePerAtomRMSDRing(matchedBasal2, startingIndex, rmsdList2,
  //                              rmsdPerAtom);
  // // ------------
  // ------------
  // Update the RMSD (obtained for each ring) for each particle
  // Basal1
  match::updateRMSDRing(matchedBasal1, startingIndex, rmsd1, rmsdPerAtom);
  // Basal2
  match::updateRMSDRing(matchedBasal2, startingIndex, rmsd2, rmsdPerAtom);
  // ------------

  return true;
} // end of function

/**
 * @details For the pentagonal nanochannels in amorphous ice.
 *
 * These must be aligned.
 *
 * Shape-matching for a pair of polygon basal rings.
 * Returns true if the pair of  basal rings form a prism block.
 */
bool match::matchUntetheredPrism(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, const Eigen::MatrixXd &refPoints,
    std::vector<int> *basal1, std::vector<int> *basal2,
    std::vector<double> *rmsdPerAtom) {
  //
  int ringSize = (*basal1).size(); // Number of nodes in each basal ring
  std::vector<int> matchedBasal1,
      matchedBasal2; // Re-ordered basal rings 1 and 2
  int dim = 3;       // Number of dimensions
  // Matrices for the point sets of the target basal rings
  Eigen::MatrixXd basal1Set(ringSize,
                            dim); // Point set for the first basal ring
  Eigen::MatrixXd basal2Set(ringSize,
                            dim); // Point set for the second basal ring
  int startingIndex; // Index in the basal rings from which the reference point
                     // set is matched
  double angDist;    // Calculated angular distance between the two basal rings
  // Variables for the absolute orientation
  std::vector<double> quat1, quat2; // quaternion rotation
  double rmsd1, rmsd2;              // least total RMSD
  std::vector<double> rmsdList1,
      rmsdList2;         // List of RMSD per atom in the order fed in
  double scale1, scale2; // Scale factor
  // Deformed block classification
  bool doAngleCriterion = false; // For deformed blocks

  // Getting the target Eigen vectors
  // Get the re-ordered matched basal rings, ordered with respect to each other
  pntToPnt::relOrderPrismBlock(yCloud, *basal1, *basal2, &matchedBasal1,
                               &matchedBasal2);
  // -----------------------
  // Match the basal rings with a complete prism block, given the relatively
  // ordered basal rings This actually only needs to be done for deformed prism
  // blocks
  bool blockMatch = match::matchPrismBlock(
      yCloud, nList, refPoints, &matchedBasal1, &matchedBasal2, &startingIndex);
  // -----------------------
  // Check to see if the prism matches the reference prism block
  if (!blockMatch) {
    return blockMatch;
  }

  // If the deformed prism does not fulfil the shape matching criterion, then do
  // not calculate the RMSD per atom
  // -----------------------
  // Section for RMSD per atom, based on basal ring matching

  basal1Set =
      pntToPnt::fillPointSetPrismRing(yCloud, matchedBasal1, startingIndex);
  // Fill up the point set for basal2
  basal2Set =
      pntToPnt::fillPointSetPrismRing(yCloud, matchedBasal2, startingIndex);
  // Use Horn's algorithm to calculate the absolute orientation and RMSD etc.
  absor::hornAbsOrientation(refPoints, basal1Set, &quat1, &rmsd1, &rmsdList1,
                            &scale1); // basal1
  absor::hornAbsOrientation(refPoints, basal2Set, &quat2, &rmsd2, &rmsdList2,
                            &scale2); // basal2
  // // Calculate the angular distance between basal1 and basal2
  // angDist = gen::angDistDegQuaternions(quat1, quat2);
  // // Check if the shapes are aligned
  // if (angDist > cutoffAngle) {
  //   return false;
  // } // If not aligned, it is not a prism block
  // ------------
  // Update the RMSD (obtained for each ring) for each particle
  // Basal1
  match::updateRMSDRing(matchedBasal1, startingIndex, rmsd1, rmsdPerAtom);
  // Basal2
  match::updateRMSDRing(matchedBasal2, startingIndex, rmsd2, rmsdPerAtom);
  // ------------

  return true;
} // end of function

//! Update the per-particle RMSD for a prism block basal ring.
int match::updatePerAtomRMSDRing(std::vector<int> basalRing, int startingIndex,
                                 std::vector<double> rmsdFromMatch,
                                 std::vector<double> *rmsdPerAtom) {
  //
  int atomIndex;                   // Atom particle index, in rmsdPerAtom
  int ringSize = basalRing.size(); // Size of the basal ring
  int index;                       // Corresponding index in the basal ring
  double iRMSD;                    // Per-particle RMSD

  // -----------------
  // Update the RMSD per particle
  for (int i = 0; i < ringSize; i++) {
    index = i + startingIndex; // Corresponding index in the basal ring
    // Wrap-around
    if (index >= ringSize) {
      index -= ringSize;
    } // end of wrap-around
    //
    // Get the atom index
    atomIndex = basalRing[index];
    // Get and update the RMSD per particle
    iRMSD = rmsdFromMatch[i];
    if ((*rmsdPerAtom)[atomIndex] == -1) {
      (*rmsdPerAtom)[atomIndex] = iRMSD;
    } // end of update of RMSD

  } // end of updating the RMSD per particle
  // -----------------
  // finito
  return 0;
} // end of function

//! Update the RMSD for each particle with the RMSD of each ring for a prism
//! block basal ring.
int match::updateRMSDRing(std::vector<int> basalRing, int startingIndex,
                          double rmsdVal, std::vector<double> *rmsdPerAtom) {
  //
  int atomIndex;                   // Atom particle index, in rmsdPerAtom
  int ringSize = basalRing.size(); // Size of the basal ring
  int index;                       // Corresponding index in the basal ring

  // -----------------
  // Update the RMSD per particle
  for (int i = 0; i < ringSize; i++) {
    index = i + startingIndex; // Corresponding index in the basal ring
    // Wrap-around
    if (index >= ringSize) {
      index -= ringSize;
    } // end of wrap-around
    //
    // Get the atom index
    atomIndex = basalRing[index];
    // The RMSD value to update with is rmsdVal
    if ((*rmsdPerAtom)[atomIndex] == -1) {
      (*rmsdPerAtom)[atomIndex] = rmsdVal;
    } // end of update of RMSD

  } // end of updating the RMSD for each particle
  // -----------------
  // finito
  return 0;
} // end of function

//! Shape-matching for a pair of polygon basal rings, comparing with a complete
//! prism block. Returns true if the pair of basal rings form a prism block.
bool match::matchPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, const Eigen::MatrixXd &refPoints,
    std::vector<int> *basal1, std::vector<int> *basal2, int *beginIndex) {
  //
  int ringSize = (*basal1).size(); // Number of nodes in the basal rings
  bool isMatch; // Qualifier for whether the prism block matches or not
  int dim = 3;  // Number of dimensions
  int nop = 2 * ringSize; // Number of particles in the prism block
  Eigen::MatrixXd refPrismBlock(
      nop, dim); // eigen matrix for the reference prism block
  Eigen::MatrixXd targetPrismBlock(
      nop, dim);     // eigen matrix for the target prism block
  int startingIndex; // Index from which the basal rings will be ordered
                     // (permutations)
  // variables for shape-matching
  std::vector<double> quat;     // Quaternion for the rotation
  double rmsd;                  // RMSD value for the entire shape
  std::vector<double> rmsdList; // RMSD value for each point
  double scale;                 // Scale factor
  // ----------------------------------------------------
  // Get the reference prism block
  refPrismBlock =
      pntToPnt::createPrismBlock(yCloud, refPoints, ringSize, *basal1, *basal2);
  // ----------------------------------------------------

  // Loop through possible point-to-point correspondences
  if (ringSize % 2 == 0 || ringSize == 3) {
    startingIndex = 0;
    // Fill up the point set for the target prism block
    targetPrismBlock =
        pntToPnt::fillPointSetPrismBlock(yCloud, *basal1, *basal2, 0);
    // Shape-matching
    absor::hornAbsOrientation(refPrismBlock, targetPrismBlock, &quat, &rmsd,
                              &rmsdList,
                              &scale); // basal2
  }                                    // even or if there are 3 nodes
  else {
    // Define the vector, RMSD etc:
    std::vector<double> currentQuat; // quaternion rotation
    double currentRmsd;              // least total RMSD
    std::vector<double>
        currentRmsdList; // List of RMSD per atom in the order fed in
    double currentScale;
    // Loop through all possible startingIndex
    for (int i = 0; i < ringSize; i++) {
      //
      // Fill up the point set for basal1
      targetPrismBlock =
          pntToPnt::fillPointSetPrismBlock(yCloud, *basal1, *basal2, i);
      // Use Horn's algorithm to calculate the absolute orientation and RMSD
      // etc. for basal1
      absor::hornAbsOrientation(refPrismBlock, targetPrismBlock, &currentQuat,
                                &currentRmsd, &currentRmsdList, &currentScale);
      // Comparison to get the least RMSD for the correct mapping
      if (i == 0) {
        // Init
        quat = currentQuat;
        rmsd = currentRmsd;
        rmsdList = currentRmsdList;
        scale = currentScale;
        startingIndex = i;
      } // init
      else {
        // Check to see if the calculated RMSD is less than the RMSD already
        // saved
        if (currentRmsd < rmsd) {
          quat = currentQuat;
          rmsd = currentRmsd;
          rmsdList = currentRmsdList;
          scale = currentScale;
          startingIndex = i;
        } // end of check to see if the current RMSD is smaller
      }   // end of comparison and filling
    }     // end of loop through startingIndex
  }       // ringSize is odd, so every point must be tried

  // Update the index from which matching is started
  *beginIndex = startingIndex;

  // // TEST DELETE BLOCK LATER
  // if (ringSize == 6) {
  //   std::fstream rmsdFile;
  //   rmsdFile.open("../runOne/topoINT/rmsd6.dat",
  //                 std::fstream::in | std::fstream::out | std::fstream::app);
  //   rmsdFile << rmsd << "\n";
  //   rmsdFile.close();
  // } else if (ringSize == 4) {
  //   std::fstream rmsdFile;
  //   rmsdFile.open("../runOne/topoINT/rmsd4.dat",
  //                 std::fstream::in | std::fstream::out | std::fstream::app);
  //   rmsdFile << rmsd << "\n";
  //   rmsdFile.close();
  // }
  // // END OF BLOCK TO BE DELETED

  // Condition for shape-matching
  if (rmsd <= 6) {
    isMatch = true; // Is a prism block
  }                 // within RMSD range
  else {
    isMatch = false;
  }

  // Return
  return isMatch;
}
