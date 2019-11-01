#include <topo_one_dim.hpp>

// -----------------------------------------------------------------------------------------------------
// PRISM ALGORITHMS
// -----------------------------------------------------------------------------------------------------

/********************************************/ /**
 *  Determines which rings are n-sided prisms. This function
 returns a vector which contains the ring IDs of all the rings which are prisms.
 The ring IDs correspond to the index of the rings inside the vector of vector
 rings, starting from 0. Prism rings can be found using a three-step procedure,
 in which first two basal rings are found. Prismatic rings are simply rings
 which share every face made by upper and lower triplets of the basal rings The
 neighbour list is also required as an input, which is a vector of vectors,
 containing atom IDs. The first element of the neighbour list is the atomID of
 the atom for which the other elements are nearest neighbours.
 ***********************************************/
std::vector<int> ring::findPrisms(
    std::vector<std::vector<int>> rings, std::vector<ring::strucType> *ringType,
    int *nPrisms, std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  std::vector<int> listPrism;
  int totalRingNum = rings.size();  // Total number of rings
  std::vector<int> basal1;          // First basal ring
  std::vector<int> basal2;          // Second basal ring
  bool cond1, cond2;  // Conditions for rings to be basal (true) or not (false)
  bool relaxedCond;   // Condition so that at least one bond exists between the
                      // two basal rings
  bool isAxialPair;   // Basal rings should be parallel in one dimension to
                      // prevent overcounting
  int ringSize = rings[0].size();  // Number of nodes in each ring
  int nDeformedPrisms = 0;         // Number of undeformed prisms
  *nPrisms = 0;                    // Number of undeformed prisms

  // Two loops through all the rings are required to find pairs of basal rings
  for (int iring = 0; iring < totalRingNum - 1; iring++) {
    cond1 = false;
    cond2 = false;
    basal1 = rings[iring];  // Assign iring to basal1
    // Loop through the other rings to get a pair
    for (int jring = iring + 1; jring < totalRingNum; jring++) {
      basal2 = rings[jring];  // Assign jring to basal2
      // ------------
      // Put extra check for tetragonal prism blocks to prevent overcounting
      if (ringSize == 4) {
        isAxialPair = false;  // init
        isAxialPair =
            ring::discardExtraTetragonBlocks(&basal1, &basal2, yCloud);
        if (isAxialPair == false) {
          continue;
        }
      }  // end of check for tetragonal prism blocks
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
        // Check for the reduced criteria fulfilment
        relaxedCond = ring::relaxedPrismConditions(nList, &basal1, &basal2);
        // Skip if relaxed criteria are not met
        if (relaxedCond == false) {
          continue;
        }  // end of skipping if the prisms do not fulfil relaxed criteria

        // Write outs
        nDeformedPrisms += 1;
        // Now write out axial basal rings
        sout::writeBasalRingsPrism(&basal1, &basal2, nDeformedPrisms, nList,
                                   yCloud, true);
      }  // end of reduced criteria
      // Strict criteria
      else {
        // Update the number of prism blocks
        *nPrisms += 1;
        // Update iring
        if ((*ringType)[iring] == ring::unclassified) {
          (*ringType)[iring] = ring::Prism;
          listPrism.push_back(iring);
        }
        // Update jring
        if ((*ringType)[jring] == ring::unclassified) {
          (*ringType)[jring] = ring::Prism;
          listPrism.push_back(jring);
        }
        // Now write out axial basal rings for convex hull calculations
        sout::writePrisms(&basal1, &basal2, *nPrisms, yCloud);
        // Write out prisms for shape-matching
        sout::writeBasalRingsPrism(&basal1, &basal2, *nPrisms, nList, yCloud,
                                   false);
        // -----------
      }  // end of strict criteria

    }  // end of loop through rest of the rings to get the second basal ring
  }    // end of loop through all rings for first basal ring

  sort(listPrism.begin(), listPrism.end());
  auto ip = std::unique(listPrism.begin(), listPrism.end());
  // Resize peripheral rings to remove undefined terms
  listPrism.resize(std::distance(listPrism.begin(), ip));

  return listPrism;
}

/********************************************/ /**
 *  Check to see if two basal rings are basal rings of a prism block or not,
 using the neighbour list information. The neighbour list nlist is a vector of
 vectors, containing atom IDs (not vector indices!). The first element of each
 subvector in nlist is the atom ID of the particle for which the other elements
 are the nearest neighbours
 ***********************************************/
bool ring::basalPrismConditions(std::vector<std::vector<int>> nList,
                                std::vector<int> *basal1,
                                std::vector<int> *basal2) {
  int l1 = (*basal1)[0];  // first element of basal1 ring
  int ringSize =
      (*basal1).size();  // Size of the ring; each ring contains n elements
  int m_k;               // Atom ID of element in basal2
  bool l1_neighbour;     // m_k is a neighbour of l1(true) or not (false)

  // isNeighbour is initialized to false for all basal2 elements; indication if
  // basal2 elements are neighbours of basal1
  std::vector<bool> isNeighbour(ringSize, false);
  int kIndex;   // m_k index
  int lAtomID;  // atomID of the current element of basal1
  int kAtomID;  // atomID of the current element of basal2

  // ---------------------------------------------
  // COMPARISON OF basal2 ELEMENTS WITH l1
  for (int k = 0; k < ringSize; k++) {
    l1_neighbour = false;
    m_k = (*basal2)[k];
    // =================================
    // Checking to seee if m_k is be a neighbour of l1
    // Find m_k inside l1 neighbour list
    auto it = std::find(nList[l1 - 1].begin() + 1, nList[l1 - 1].end(), m_k);

    // If the element has been found, for l1
    if (it != nList[l1 - 1].end()) {
      l1_neighbour = true;
      kIndex = k;
      break;
    }
  }  // l1 is a neighbour of m_k

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
    lAtomID = (*basal1)[i];  // element of basal1 ring
    for (int k = 0; k < ringSize; k++) {
      // Skip if already a neighbour
      if (isNeighbour[k]) {
        continue;
      }
      // Get the comparison basal2 element
      kAtomID = (*basal2)[k];  // element of basal2 ring;

      // Checking to see if kAtomID is a neighbour of lAtomID
      // Find kAtomID inside lAtomID neighbour list
      auto it1 = std::find(nList[lAtomID - 1].begin() + 1,
                           nList[lAtomID - 1].end(), kAtomID);

      // If the element has been found, for l1
      if (it1 != nList[lAtomID - 1].end()) {
        isNeighbour[k] = true;
      }
    }  // Loop through basal2
  }    // end of check for neighbours of basal1

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

/********************************************/ /**
                                                *  Relaxed criteria for deformed
                                                *prism blocks: at least one bond
                                                *should exist between the basal
                                                *rings
                                                ***********************************************/
bool ring::relaxedPrismConditions(std::vector<std::vector<int>> nList,
                                  std::vector<int> *basal1,
                                  std::vector<int> *basal2) {
  int ringSize =
      (*basal1).size();      // Size of the ring; each ring contains n elements
  int m_k;                   // Atom ID of element in basal2
  bool isNeighbour = false;  // This is true if there is at least one bond
                             // between the basal rings
  int l_k;                   // Atom ID of element in basal1

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
      auto it =
          std::find(nList[l_k - 1].begin() + 1, nList[l_k - 1].end(), m_k);

      // If the element has been found, for l1
      if (it != nList[l_k - 1].end()) {
        isNeighbour = true;
        break;
      }  // found element
    }    // end of loop through all the elements of basal2

    // If a neighbour has been found then
    if (isNeighbour) {
      return true;
    }
  }  // end of loop through all the elements of basal1

  // If a neighbour has not been found, return false
  return false;
}

/********************************************/ /**
 *  Discards basal tetragonal rings which are not oriented perpendicular to the
 z dimension. This is hard-coded for the z-dimension but any dimension is
 equivalent. If true, then the rings are axial.
 ***********************************************/
bool ring::discardExtraTetragonBlocks(
    std::vector<int> *basal1, std::vector<int> *basal2,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  int ringSize =
      (*basal1).size();    // Size of the ring; each ring contains n elements
  bool isHigher, isLower;  // For figuring out if basal1 is in above or beneath
                           // basal2, in the +X direction
  int iatomIndex,
      jatomIndex;   // Indices of the elements in basal1 and basal2 respectively
  double z_i, z_j;  // Z coordinates of iatom and jatom of basal1 and basal2
                    // respectively
  // Variables for getting the projected area
  bool axialBasal1, axialBasal2;  // bools for checking if basal1 and basal2 are
                                  // axial (true) respectively
  double areaXY, areaXZ,
      areaYZ;  // Projected area on the XY, XZ and YZ planes respectively

  // ----------------------------------------
  // Calculate projected area onto the XY, YZ and XZ planes for basal1
  axialBasal1 = false;  // Init to false
  axialBasal2 = false;  // Init

  // Init the projected area
  areaXY = 0.0;
  areaXZ = 0.0;
  areaYZ = 0.0;

  jatomIndex = (*basal1)[0] - 1;

  // All points except the first pair
  for (int k = 1; k < ringSize; k++) {
    iatomIndex = (*basal1)[k] - 1;  // Current vertex

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
  iatomIndex = (*basal1)[0] - 1;
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

  // Hard-coded for the z-axis. Check if xy projected area is the greatest
  if (areaXY > areaXZ && areaXY > areaYZ) {
    axialBasal1 = true;
  }  // end of check for axial ring for basal1
  // ----------------------------------------
  // Calculate projected area onto the XY, YZ and XZ planes for basal2

  // Init the projected area
  areaXY = 0.0;
  areaXZ = 0.0;
  areaYZ = 0.0;

  jatomIndex = (*basal2)[0] - 1;

  // All points except the first pair
  for (int k = 1; k < ringSize; k++) {
    iatomIndex = (*basal2)[k] - 1;  // Current vertex

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
  iatomIndex = (*basal2)[0] - 1;
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

  // Hard-coded for the z-axis. Check if xy projected area is the greatest
  if (areaXY > areaXZ && areaXY > areaYZ) {
    axialBasal2 = true;
  }  // end of check for axial ring for basal1
  // ----------------------------------------

  // Now check if basal1 and basal2 are axial or not
  if (axialBasal1 == true && axialBasal2 == true) {
    return true;
  } else {
    return false;
  }

  // ----------------------------------------

  // Get the one element of basal1
  iatomIndex = (*basal1)[0] - 1;    // C++ indices are one less
  z_i = yCloud->pts[iatomIndex].z;  // z coordinate of element 0 of basal1
  // First element of basal2
  jatomIndex = (*basal2)[0] - 1;    // C++ indices are one less
  z_j = yCloud->pts[jatomIndex].z;  // z coordinate of element 0 of basal2

  // Init
  isHigher = false;
  isLower = false;

  if (z_i > z_j) {
    isHigher = true;
  }  // basal1 is 'above' basal2
  else if (z_i < z_j) {
    isLower = true;
  }  // basal1 is 'below' basal2
  else {
    return false;
  }  // Basal1 and basal2 are 'lateral rings'

  // Check to see that every element of basal1 is above or below basal2
  // Loop through the elements of basal1 and basal2
  for (int k = 0; k < ringSize; k++) {
    iatomIndex = (*basal1)[k] - 1;    // C++ indices are one less
    z_i = yCloud->pts[iatomIndex].z;  // z coordinate of iatom of basal1
    // Now loop through every element of basal2
    for (int l = 0; l < ringSize; l++) {
      jatomIndex = (*basal2)[l] - 1;    // C++ indices are one less
      z_j = yCloud->pts[jatomIndex].z;  // z coordinate of basal2

      // when basal1 is higher than basal2
      if (isHigher) {
        if (z_i <= z_j) {
          return false;
        }
      }
      //
      // when basal1 is lower
      if (isLower) {
        if (z_i >= z_j) {
          return false;
        }
      }
      //
    }  // end of loop through basal2
  }    // end of loop through elements of basal1

  return true;
}
