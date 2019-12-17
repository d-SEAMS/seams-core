#include <shapeMatch.hpp>

// Shape-matching for a pair of polygon basal rings. Returns true if the pair of
// basal rings form a prism block.
bool match::matchPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, const Eigen::MatrixXd &refPoints,
    std::vector<int> *basal1, std::vector<int> *basal2,
    std::vector<double> *rmsdPerAtom, bool isPerfect) {
  //
  int ringSize = (*basal1).size();  // Number of nodes in each basal ring
  std::vector<int> matchedBasal1,
      matchedBasal2;  // Re-ordered basal rings 1 and 2
  int dim = 3;        // Number of dimensions
  // Matrices for the point sets of the target basal rings
  Eigen::MatrixXd basal1Set(ringSize,
                            dim);  // Point set for the first basal ring
  Eigen::MatrixXd basal2Set(ringSize,
                            dim);  // Point set for the second basal ring
  int startingIndex;  // Index in the basal rings from which the reference point
                      // set is matched
  double cutoffAngle = 15;  // Tolerance for the angular distance between the
                            // quaternions of the basal rings
  double angDist;  // Calculated angular distance between the two basal rings
  // Variables for the absolute orientation
  std::vector<double> quat1, quat2;  // quaternion rotation
  double rmsd1, rmsd2;               // least total RMSD
  std::vector<double> rmsdList1,
      rmsdList2;          // List of RMSD per atom in the order fed in
  double scale1, scale2;  // Scale factor
  // -----------------------
  // Getting the target Eigen vectors
  // Get the re-ordered matched basal rings, ordered with respect to each other
  pntToPnt::relOrderPrismBlock(yCloud, *basal1, *basal2, nList, &matchedBasal1,
                               &matchedBasal2);
  // -----------------------
  // Loop through possible point-to-point correspondences
  if (ringSize % 2 == 0 || ringSize == 3) {
    startingIndex = 0;
    // Fill up the point set for basal1
    basal1Set =
        pntToPnt::fillPointSetPrismRing(yCloud, matchedBasal1, startingIndex);
    // Fill up the point set for basal2
    basal2Set =
        pntToPnt::fillPointSetPrismRing(yCloud, matchedBasal2, startingIndex);
    // Use Horn's algorithm to calculate the absolute orientation and RMSD etc.
    absor::hornAbsOrientation(refPoints, basal1Set, &quat1, &rmsd1, &rmsdList1,
                              &scale1);  // basal1
    absor::hornAbsOrientation(refPoints, basal2Set, &quat2, &rmsd2, &rmsdList2,
                              &scale2);  // basal2
  }                                      // even or if there are 3 nodes
  else {
    // Define the vector, RMSD etc:
    std::vector<double> currentQuat1;  // quaternion rotation
    double currentRmsd1;               // least total RMSD
    std::vector<double>
        currentRmsdList1;  // List of RMSD per atom in the order fed in
    double currentScale;
    // Loop through all possible startingIndex
    for (int i = 0; i < ringSize; i++) {
      //
      // Fill up the point set for basal1
      basal1Set = pntToPnt::fillPointSetPrismRing(yCloud, matchedBasal1, i);
      // Use Horn's algorithm to calculate the absolute orientation and RMSD
      // etc. for basal1
      absor::hornAbsOrientation(refPoints, basal1Set, &currentQuat1,
                                &currentRmsd1, &currentRmsdList1,
                                &currentScale);
      // Comparison to get the least RMSD for the correct mapping
      if (i == 0) {
        // Init
        quat1 = currentQuat1;
        rmsd1 = currentRmsd1;
        rmsdList1 = currentRmsdList1;
        scale1 = currentScale;
        startingIndex = i;
      }  // init
      else {
        // Check to see if the calculated RMSD is less than the RMSD already
        // saved
        if (currentRmsd1 < rmsd1) {
          quat1 = currentQuat1;
          rmsd1 = currentRmsd1;
          rmsdList1 = currentRmsdList1;
          scale1 = currentScale;
          startingIndex = i;
        }  // end of check to see if the current RMSD is smaller
      }    // end of comparison and filling
    }      // end of loop through startingIndex
           //
    // Now that the startingIndex has been obtained, get the second basal ring
    // absolute orientation too
    // Fill up the point set for basal2
    basal2Set =
        pntToPnt::fillPointSetPrismRing(yCloud, matchedBasal2, startingIndex);
    // Use Horn's algorithm to calculate the absolute orientation and RMSD
    // etc. for basal1
    absor::hornAbsOrientation(refPoints, basal2Set, &quat2, &rmsd2, &rmsdList2,
                              &scale2);
  }  // ringSize is odd, so every point must be tried

  // ------------
  // Update the RMSD for each particle
  // Basal1
  match::updatePerAtomRMSDRing(matchedBasal1, startingIndex, rmsdList1,
                               rmsdPerAtom);
  // Basal2
  match::updatePerAtomRMSDRing(matchedBasal2, startingIndex, rmsdList2,
                               rmsdPerAtom);
  // ------------
  // The basal rings of a deformed block should have an angular distance within
  // a cutoff angle.
  if (!isPerfect) {
    angDist = gen::angDistDegQuaternions(quat1, quat2);
    // within angle cutoff:
    // within angle tolerance of 0 degrees
    if (angDist > -1 * cutoffAngle && angDist < cutoffAngle) {
      return true;
    }  // angle tolerance satisfied
    else if (angDist > 180 - 1 * cutoffAngle && angDist < 180 + cutoffAngle) {
      return true;
    }  // angle tolerance satisfied
    else if (angDist > 360 - 1 * cutoffAngle && angDist < 360 + cutoffAngle) {
      return true;
    }  // angle tolerance satisfied
    else {
      return false;
    }  // angle tolerance not satisfied
  }    // Classification criterion for deformed blocks

  // Delete this later
  angDist = gen::angDistDegQuaternions(quat1, quat2);
  //
  if (angDist > -1 * cutoffAngle && angDist < cutoffAngle) {
    //
  }  // angle tolerance satisfied
  else if (angDist > 180 - 1 * cutoffAngle && angDist < 180 + cutoffAngle) {
    //
  }  // angle tolerance satisfied
  else {
    // std::cout << "ringSize is " << ringSize << " for angle= " << angDist
    //           << "\n";
  }  // angle tolerance not satisfied
  //
  // Write out to file
  std::fstream temp;
  temp.open("../runOne/topoINT/angles.dat",
            std::fstream::in | std::fstream::out | std::fstream::app);
  temp << angDist << "\n";
  temp.close();
  //

  return true;
}  // end of function

// Update the per-particle RMSD for a prism block basal ring.
int match::updatePerAtomRMSDRing(std::vector<int> basalRing, int startingIndex,
                                 std::vector<double> rmsdFromMatch,
                                 std::vector<double> *rmsdPerAtom) {
  //
  int atomIndex;                    // Atom particle index, in rmsdPerAtom
  int ringSize = basalRing.size();  // Size of the basal ring
  int index;                        // Corresponding index in the basal ring
  double iRMSD;                     // Per-particle RMSD

  // -----------------
  // Reorder the basal ring elements
  for (int i = 0; i < ringSize; i++) {
    index = i + startingIndex;  // Corresponding index in the basal ring
    // Wrap-around
    if (index >= ringSize) {
      index -= ringSize;
    }  // end of wrap-around
    //
    // Get the atom index
    atomIndex = basalRing[index];
    // Get and update the RMSD per particle
    iRMSD = rmsdFromMatch[i];
    if ((*rmsdPerAtom)[atomIndex] == -1) {
      (*rmsdPerAtom)[atomIndex] = iRMSD;
    }  // end of update of RMSD

  }  // end of reordering and update
  // -----------------
  // finito
  return 0;
}  // end of function