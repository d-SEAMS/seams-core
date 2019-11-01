#include <ring.hpp>

/********************************************/ /**
                                                *  Deletes the memory of a
                                                *vector of vectors
                                                ***********************************************/
int ring::clearRingList(std::vector<std::vector<int>> &rings) {
  //
  std::vector<std::vector<int>> tempEmpty;

  rings.swap(tempEmpty);

  return 0;
}

/********************************************/ /**
                                                *  Gets rings of a single ring
                                                *size from all primitive rings
                                                *and returns that vector of
                                                *vectors
                                                ***********************************************/
std::vector<std::vector<int>> ring::getSingleRingSize(
    std::vector<std::vector<int>> rings, int ringSize) {
  //
  std::vector<std::vector<int>> ringSingleSize;  // rings of one size

  // rings contains primitive rings of all sizes
  // Only save rings of a given size (ringSize) to the new
  // vector of vectors, ringSingleSize
  for (int iring = 0; iring < rings.size(); iring++) {
    // Check the size of the current ring
    // If it is the correct size, save it in ringSingleSize
    if (rings[iring].size() == ringSize) {
      ringSingleSize.push_back(rings[iring]);
    }  // End of check of the size of iring
  }    // end of loop through all rings in rings

  return ringSingleSize;
}

/********************************************/ /**
 *  Checks if the ring has more than three consecutive
 water molecules or not. Returns false if a quadruplet is within
 the neighbour list.
 ***********************************************/
bool ring::checkRing(std::vector<int> ring,
                     std::vector<std::vector<int>> nList) {
  int ringSize = ring.size();  // Size of the ring (no. of atoms in each ring)
  std::vector<int> subRing;    //  Quadruplet formed from a ring
  int j;

  for (int k = 0; k < ringSize; k++) {
    subRing.clear();  // Clear the quadruplet
    // Get 4 elements
    for (int i = k; i < k + 4; i++) {
      j = i;
      if (i >= ringSize) {
        j = i - ringSize;
      }
      subRing.push_back(ring[j]);
    }  // end of getting a quadruplet from k
    // Check to see that the quadruplet isn't part
    // of the neighbour list
    if (ring::compareQuad(subRing, nList) == true) {
      return false;
    }  // end of comparison
  }    // end of looping through all possible quadruplets in ring

  return true;
}

/********************************************/ /**
 *  Checks if a quadruplet has more than three consecutive
 water molecules or not. Returns true if the quadruplet is within
 the neighbour list. The neighbour list is sorted already
 ***********************************************/
bool ring::compareQuad(std::vector<int> quad,
                       std::vector<std::vector<int>> nList) {
  int iatom = quad[0] - 1;  // Index of the first atom ID in the quadruplet
  int atomID;
  int jatom;
  int flag = 0;
  int nNumNeighbours =
      nList[iatom].size();  // Number of nearest neighbours for iatom+1

  for (int k = 1; k < 4; k++) {
    atomID = quad[k];
    flag = 0;  // init to zero
    for (int j = 1; j < nNumNeighbours; j++) {
      jatom = nList[iatom][j];  // Get j^th element of the iatom row in vector
                                // of vectors nList
      if (jatom == atomID) {
        flag++;
      }
      if (jatom > atomID) {
        continue;
      }  // sorted list
    }    // end of loop through neighbour list
    int r = 3;
    // If the flag is three, then the quadruplet is within the nearest
    // neighbours
    if (flag == 3) {
      return true;
    }
  }  // end of looping through elements in quad

  return false;
}

/********************************************/ /**
 *  Determines which hexagonal rings are DDC rings. This function
 returns a vector which contains the ring IDs of all the rings which are DDC
 rings. The ring IDs correspond to the index of the rings inside the vector of
 vector rings, starting from 0. DDC rings can be found using a three-step
 procedure, in which equatorial rings and their corresponding rings can be
 found. Peripheral rings can be shared, so first all the equatorial rings must
 be found. Whenever an equatorial ring and the corresponding peripheral rings
 are found, their values must be updated to DDC enum type inside the ringType
 vector, which has been passed as an input to this function. At the end, the
 output vector can be updated to avoid adding the same ring index more than once
 (this may happen for peripheral rings, which can be shared)
 ***********************************************/
std::vector<int> ring::findDDC(std::vector<std::vector<int>> rings,
                               std::vector<ring::strucType> *ringType,
                               std::vector<cage::Cage> *cageList) {
  std::vector<int> listDDC;
  int totalRingNum = rings.size();   // Total number of hexagonal rings
  std::vector<int> peripheralRings;  // Indices which may be peripheral rings
  bool cond1, cond2, cond3;          // Conditions for DDC
  int jring;                         // Index for peripheral ring being added
  std::vector<int> DDCRings;  // Indices of rings which constitute a single DDC,
                              // with the equatorial ring first

  // To search for equatorial rings, loop through all
  // the hexagonal rings
  for (int iring = 0; iring < totalRingNum; iring++) {
    peripheralRings.clear();
    // ------------
    // Step one: Find all rings which contain each index (m_k) of the equatorial
    // ring iring in at least three other rings
    cond1 = ring::conditionOneDDC(rings, &peripheralRings, iring);
    if (cond1 == false) {
      continue;
    }
    // ------------
    // Step two: For every triplet in iring, there is at least one
    // hexagonal ring other than iring that passes through the triplet
    cond2 = ring::conditionTwoDDC(rings, &peripheralRings, iring);
    if (cond2 == false) {
      continue;
    }
    // ------------
    // Step three: For every triplet in iring, there is at least one
    // hexagonal ring other than iring that passes through the triplet
    cond3 = ring::conditionThreeDDC(rings, &peripheralRings, iring);
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
    if ((*ringType)[iring] == ring::unclassified) {
      (*ringType)[iring] = ring::DDC;
      listDDC.push_back(iring);
    }
    // Add the peripheral ring IDs too
    for (int j = 0; j < peripheralRings.size(); j++) {
      jring = peripheralRings[j];
      if ((*ringType)[jring] == ring::unclassified) {
        (*ringType)[jring] = ring::DDC;
        listDDC.push_back(jring);
      }
    }  // end of update for peripheral rings
    // Add rings to the cageList vector of struct Cages.
    DDCRings.clear();           // init
    DDCRings.push_back(iring);  // Add the equatorial ring first
    DDCRings.insert(std::end(DDCRings), std::begin(peripheralRings),
                    std::end(peripheralRings));
    (*cageList).push_back({cage::DoubleDiaC, DDCRings});
    // ------------
  }  // end of loop through all hexagonal rings

  return listDDC;
}

/********************************************/ /**
 *  For a given ring, which is being tested as the equatorial ring,
 this function tests if each element of the ring (m_k) is present
 in at least three other rings or not. Returns false if this is not true.
 The ring IDs of rings that have the common element are saved inside
 the periperal ring vector as potential peripheral rings, which is passed
 as an input to the function
 ***********************************************/
bool ring::conditionOneDDC(std::vector<std::vector<int>> rings,
                           std::vector<int> *peripheralRings, int iring) {
  int index;  // Atom ID to be compared
  int noOfCommonRings =
      0;  // No of rings in which the element to be matched has been found
  int jElement;  // Atom ID being compared to index

  // Loop through each element of iring for finding matches
  for (int m = 0; m < 6; m++) {
    index = rings[iring][m];  // Atom ID to be compared and matched with
    noOfCommonRings = 0;      // init to zero.
    // Loop through every ring except iring
    for (int jring = 0; jring < rings.size(); jring++) {
      if (iring == jring) {
        continue;
      }  // Skip for iring
      // -------
      // Search every element of jring
      for (int k = 0; k < 6; k++) {
        jElement = rings[jring][k];
        if (jElement == index) {
          noOfCommonRings++;
          (*peripheralRings).push_back(jring);
          break;
        }  // if index is found inside jring
        else {
          continue;
        }
      }  // end of loop through every element of jring
      // -------
    }  // end of loop through all rings except iring
    // If less than 3 rings have been found for each element, then this is
    // not an equatorial ring
    if (noOfCommonRings < 3) {
      return false;
    }  // end of check for common ring number per element in iring
  }    // end of loop through each element of iring

  // iring is an equatorial ring. The duplicate ring IDs inside
  // peripheralRings should be removed
  std::vector<int>::iterator ip;  // Iterator to find the last element upto
                                  // which unique elements are present
  // Duplicate IDs must be removed
  int numElements =
      (*peripheralRings).size();  // number of elements in peripheralRings
  sort((*peripheralRings).begin(), (*peripheralRings).end());
  ip = std::unique((*peripheralRings).begin(),
                   (*peripheralRings).begin() + numElements);
  // Resize peripheral rings to remove undefined terms
  (*peripheralRings).resize(std::distance((*peripheralRings).begin(), ip));

  return true;
}

/********************************************/ /**
 *  For a given ring, which is being tested as the equatorial ring,
 this function tests if each triplet that can be formed from the ring is
 common to at least one other ring or not. Returns false if this is not true.
 The ring IDs of rings that have the common triplet are ultimately saved inside
 the periperal ring vector as potential peripheral rings, which is passed
 as an input to the function
 ***********************************************/
bool ring::conditionTwoDDC(std::vector<std::vector<int>> rings,
                           std::vector<int> *peripheralRings, int iring) {
  std::vector<int> triplet;  //  Triplet formed from iring
  int ringSize = 6;          // Here, all the rings are hexagons
  int j;                     // Used for making the triplet
  int jring;                 // Peripheral ring ID to be searched
  int count;                 // Number of  rings found that match the triplet
  std::vector<int>
      newPeripherals;  // Vector in which the new peripheral ring IDs are saved.
                       // This will be swapped with peripheralRings later

  for (int k = 0; k < ringSize; k++) {
    triplet.clear();  // Clear the triplet
    // Get a triplet
    for (int i = k; i < k + 3; i++) {
      j = i;
      if (i >= ringSize) {
        j = i - ringSize;
      }
      triplet.push_back(rings[iring][j]);
    }  // end of getting a triplet from k
    // -------------
    // Compare the triplet with every possible peripheral
    // ring inside peripheralRings.
    count = 0;  // init to zero
    // Loop through all possible peripheral rings
    for (int m = 0; m < (*peripheralRings).size(); m++) {
      jring = (*peripheralRings)[m];  // Ring ID of ring to be searched
      // Search inside the ring with index jring for the triplet
      const auto it = std::search(rings[jring].begin(), rings[jring].end(),
                                  triplet.begin(), triplet.end());
      // If the ring has been found inside jring
      if (it != rings[jring].end()) {
        newPeripherals.push_back(jring);  // Update new peripheral vector
        count++;
      }  // end of ring found
    }    // end of loop through all possible peripheral rings
    // If count is 0, then the triplet was not found in any peripheral ring
    if (count == 0) {
      return false;
    }  // Return false since the triplet was not found
    // -------------
  }  // end of looping through 0-6 to get triplets

  // Swap the old peripheral rings vector with the new one
  (*peripheralRings).swap(newPeripherals);

  // If there are more than 6 peripheral rings, the code will fail
  // Comment this out if you want
  if ((*peripheralRings).size() > 6) {
    std::cerr
        << "There are more than 6 peripheral rings. The code will fail. \n";
    return false;
  }  // end of check for more than 6 peripherals

  return true;
}

/********************************************/ /**
 *  For a given ring, which is being tested as the equatorial ring,
 this function tests if a set of even vector-index triplets and odd vector-index
 triplets have at least one element in common or not. Returns false if this is
 not true. If the condition is true, all six peripheral rings and the equatorial
 ring are DDC cages
 ***********************************************/
bool ring::conditionThreeDDC(std::vector<std::vector<int>> rings,
                             std::vector<int> *peripheralRings, int iring) {
  std::vector<int> triplet;  //  Triplet formed from iring
  int ringSize = 6;          // Here, all the rings are hexagons
  // Three peripheral rings, each intersecting with a triplet of iring
  std::vector<int> ring1;
  std::vector<int> ring2;
  std::vector<int> ring3;
  std::vector<int> ringIndex;  // Contains ring IDs of the peripheral indices
                               // matching with the triplets
  int j;                       // Used for making the triplet
  int jring;                   // Peripheral ring ID to be searched

  // ----------------------------------------------------------------------------
  // EVEN-INDEXED TRIPLETS
  // Search within even-indexed triplets and find ring1, ring2, ring3
  for (int k = 0; k < ringSize; k += 2) {
    triplet.clear();  // Clear the triplet
    // Get a triplet
    for (int i = k; i < k + 3; i++) {
      j = i;
      if (i >= ringSize) {
        j = i - ringSize;
      }
      triplet.push_back(rings[iring][j]);
    }  // end of getting a triplet from k
    // -------------
    // Compare the triplet with every possible peripheral
    // ring inside peripheralRings.
    // Loop through all possible peripheral rings
    for (int m = 0; m < (*peripheralRings).size(); m++) {
      jring = (*peripheralRings)[m];  // Ring ID of ring to be searched
      // Search inside the ring with index jring for the triplet
      const auto it = std::search(rings[jring].begin(), rings[jring].end(),
                                  triplet.begin(), triplet.end());
      // If the ring has been found inside jring
      if (it != rings[jring].end()) {
        // Put jring found into ringIndex
        ringIndex.push_back(jring);
        break;
      }  // end of ring found
    }    // end of loop through all possible peripheral rings
    // -------------
  }  // end of looping through 0-6 to get triplets

  // Save the even-indexed triplet-matching rings into ring1, ring2, ring3
  ring1 = rings[ringIndex[0]];
  ring2 = rings[ringIndex[1]];
  ring3 = rings[ringIndex[2]];

  // Get the intersection of 2 rings first
  sort(ring1.begin(), ring1.end());  // Sort ring1 by ID
  sort(ring2.begin(), ring2.end());  // Sort ring2 by ID
  std::vector<int>
      common2;  // Vector for holding common elements from ring1 and ring2

  // Get the intersection of sorted ring1 and ring2
  auto it1 = std::set_intersection(ring1.begin(), ring1.end(), ring2.begin(),
                                   ring2.end(), std::back_inserter(common2));
  // If the common2 vector is empty, there are no elements in common
  if (common2.size() == 0) {
    return false;
  }

  // Get intersection of the common elements from ring1, ring2 and ring3
  sort(ring3.begin(), ring3.end());  // Sort ring3 by ID
  std::vector<int>
      common3;  // Vector for holding common elements of ring1, ring2, ring3
  // Get the intersection of sorted ring1 and ring2 and ring3 and save in
  // common3
  auto it2 =
      std::set_intersection(common2.begin(), common2.end(), ring3.begin(),
                            ring3.end(), std::back_inserter(common3));
  // If there are no elements in common, then return false
  if (common3.size() == 0) {
    return false;
  }
  // ----------------------------------------------------------------------------
  // ODD-INDEXED TRIPLETS
  // ----------------------------------------------------------------------------
  // Clear all vectors used before
  ring1.clear();
  ring2.clear();
  ring3.clear();
  common2.clear();
  common3.clear();
  ringIndex.clear();

  // Search within odd-indexed triplets and find ring1, ring2, ring3
  for (int k = 1; k < ringSize; k += 2) {
    triplet.clear();  // Clear the triplet
    // Get a triplet
    for (int i = k; i < k + 3; i++) {
      j = i;
      if (i >= ringSize) {
        j = i - ringSize;
      }
      triplet.push_back(rings[iring][j]);
    }  // end of getting a triplet from k
    // -------------
    // Compare the triplet with every possible peripheral
    // ring inside peripheralRings.
    // Loop through all possible peripheral rings
    for (int m = 0; m < (*peripheralRings).size(); m++) {
      jring = (*peripheralRings)[m];  // Ring ID of ring to be searched
      // Search inside the ring with index jring for the triplet
      const auto iter = std::search(rings[jring].begin(), rings[jring].end(),
                                    triplet.begin(), triplet.end());
      // If the ring has been found inside jring
      if (iter != rings[jring].end()) {
        // Put jring found into ringIndex
        ringIndex.push_back(jring);
        break;
      }  // end of ring found
    }    // end of loop through all possible peripheral rings
    // -------------
  }  // end of looping through 0-6 to get triplets

  // Save the even-indexed triplet-matching rings into ring1, ring2, ring3
  ring1 = rings[ringIndex[0]];
  ring2 = rings[ringIndex[1]];
  ring3 = rings[ringIndex[2]];

  // Get the intersection of 2 rings first
  sort(ring1.begin(), ring1.end());  // Sort ring1 by ID
  sort(ring2.begin(), ring2.end());  // Sort ring2 by ID

  // Get the intersection of sorted ring1 and ring2
  it1 = std::set_intersection(ring1.begin(), ring1.end(), ring2.begin(),
                              ring2.end(), std::back_inserter(common2));
  // If the common2 vector is empty, there are no elements in common
  if (common2.size() == 0) {
    return false;
  }

  // Get intersection of the common elements from ring1, ring2 and ring3
  sort(ring3.begin(), ring3.end());  // Sort ring3 by ID

  // Get the intersection of sorted ring1 and ring2 and ring3 and save in
  // common3
  it2 = std::set_intersection(common2.begin(), common2.end(), ring3.begin(),
                              ring3.end(), std::back_inserter(common3));
  // If there are no elements in common, then return false
  if (common3.size() == 0) {
    return false;
  }
  // ----------------------------------------------------------------------------

  // iring is an equatorial ring and peripheralRings has the 6 peripheral rings
  return true;
}

/********************************************/ /**
 *  Determines which hexagonal rings are HCs. This function
 returns a vector which contains the ring IDs of all the rings which are HC
 rings. The ring IDs correspond to the index of the rings inside the vector of
 vector rings, starting from 0. HC rings can be found using a three-step
 procedure, in which first two basal rings are found. Prismatic rings are simply
 rings which share every face made by upper and lower triplets of the basal
 rings The neighbour list is also required as an input, which is a vector of
 vectors, containing atom IDs. The first element of the neighbour list is the
 atomID of the atom for which the other elements are nearest neighbours.
 ***********************************************/
std::vector<int> ring::findHC(std::vector<std::vector<int>> rings,
                              std::vector<ring::strucType> *ringType,
                              std::vector<std::vector<int>> nList,
                              std::vector<cage::Cage> *cageList) {
  std::vector<int> listHC;
  int totalRingNum = rings.size();  // Total number of hexagonal rings
  std::vector<int> basal1;          // First basal ring
  std::vector<int> basal2;          // Second basal ring
  bool cond1, cond2;  // Conditions for rings to be basal (true) or not (false)
  std::vector<int>
      HCRings;  // Indices of rings which constitute a single HC, with the basal
                // rings first, followed by prismatic rings
  std::vector<int> prismaticRings;  // Ring indices of prismatic rings
  int kring;                        // Ring index of the prismatic rings

  // Two loops through all the rings are required to find pairs of basal rings
  for (int iring = 0; iring < totalRingNum - 1; iring++) {
    cond1 = false;
    cond2 = false;
    basal1 = rings[iring];  // Assign iring to basal1
    // Loop through the other rings to get a pair
    for (int jring = iring + 1; jring < totalRingNum; jring++) {
      basal2 = rings[jring];  // Assign jring to basal2
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
      if ((*ringType)[iring] == ring::unclassified) {
        (*ringType)[iring] = ring::HCbasal;
        listHC.push_back(iring);
      } else if ((*ringType)[iring] == ring::DDC) {
        (*ringType)[iring] = ring::bothBasal;
        listHC.push_back(iring);
      }
      // Update jring
      if ((*ringType)[jring] == ring::unclassified) {
        (*ringType)[jring] = ring::HCbasal;
        listHC.push_back(jring);
      } else if ((*ringType)[jring] == ring::DDC) {
        (*ringType)[jring] = ring::bothBasal;
        listHC.push_back(jring);
      }
      // Find the prismatic rings
      prismaticRings.clear();  // Clear the prismatic ring vector first
      ring::findPrismatic(rings, &listHC, ringType, iring, jring,
                          &prismaticRings);
      // Update the prismatic rings
      for (int k = 0; k < prismaticRings.size(); k++) {
        kring =
            prismaticRings[k];  // Current ring index of the (3) prismatic rings
        // Update kring
        if ((*ringType)[kring] == ring::unclassified) {
          (*ringType)[kring] = ring::HCprismatic;
          listHC.push_back(kring);
        } else if ((*ringType)[kring] == ring::DDC) {
          (*ringType)[kring] = ring::bothPrismatic;
          listHC.push_back(kring);
        }
      }
      // -----------
      // Update the cageList vector of Cages
      // Update the basal rings
      HCRings.clear();
      HCRings.push_back(iring);
      HCRings.push_back(jring);
      // Add the prismaticRings
      HCRings.insert(std::end(HCRings), std::begin(prismaticRings),
                     std::end(prismaticRings));
      (*cageList).push_back({cage::HexC, HCRings});
      // -----------
    }  // end of loop through rest of the rings to get the second basal ring
  }    // end of loop through all rings for first basal ring

  sort(listHC.begin(), listHC.end());
  auto ip = std::unique(listHC.begin(), listHC.end());
  // Resize peripheral rings to remove undefined terms
  listHC.resize(std::distance(listHC.begin(), ip));

  return listHC;
}

/********************************************/ /**
 *  For two vectors, checks to see if there are common elements (true)
 or not (false)
 ***********************************************/
bool ring::hasCommonElements(std::vector<int> ring1, std::vector<int> ring2) {
  std::vector<int> commonElements;  // Vector containing common elements

  // Sort the vectors before finding common elements
  sort(ring1.begin(), ring1.end());
  sort(ring2.begin(), ring2.end());

  // Find intersection of sorted vectors
  auto it =
      std::set_intersection(ring1.begin(), ring1.end(), ring2.begin(),
                            ring2.end(), std::back_inserter(commonElements));

  // If there are no elements in common, then return false
  if (commonElements.size() == 0) {
    return false;
  }
  // If there are common elements, return true
  else {
    return true;
  }
}

/********************************************/ /**
 *  Check to see if two basal rings are HCs or not, using the neighbour list
 information. The neighbour list nlist is a vector of vectors, containing atom
 IDs (not vector indices!). The first element of each subvector in nlist is the
 atom ID of the particle for which the other elements are the nearest neighbours
 ***********************************************/
bool ring::basalConditions(std::vector<std::vector<int>> nList,
                           std::vector<int> *basal1, std::vector<int> *basal2) {
  int l1 = (*basal1)[0];  // first element of basal1 ring
  int l2 = (*basal1)[1];  // second element of basal1 ring
  int ringSize = 6;       // Size of the ring; each ring contains 6 elements
  int m_k;                // Atom ID of element in basal2
  std::vector<int> evenElements;  // contains m_k, m_{k+2}, m_{k+4}
  std::vector<int> oddElements;   // contains m_{k+1}, m_{k+3}, m_{k+5}
  int compare1, compare2;         // l3 and l5 OR l4 and l6
  int index;
  bool l1_neighbour, l2_neighbour;  // m_k is a neighbour of l1(true) or not
                                    // (false); m_k is a neighbour of l2(true)
  bool isNeigh, notNeigh;  // Used to check if the rings are basal or not

  // ---------------------------------------------
  // COMPARISON OF basal2 ELEMENTS WITH l1
  for (int k = 0; k < ringSize; k++) {
    l1_neighbour = false;
    l2_neighbour = false;
    m_k = (*basal2)[k];
    // =================================
    // Checking to seee if m_k is be a neighbour of l1
    // Find m_k inside l1 neighbour list
    auto it = std::find(nList[l1 - 1].begin() + 1, nList[l1 - 1].end(), m_k);

    // If the element has been found, for l1
    if (it != nList[l1 - 1].end()) {
      l1_neighbour = true;
      compare1 = (*basal1)[2];  // l3
      compare2 = (*basal1)[4];  // l5
    }                           // l1 is a neighbour of m_k

    // =================================
    // Checking to seee if m_k is be a neighbour of l2
    // Do this ONLY if m_k is not a neighbour of l1
    if (l1_neighbour == false) {
      auto iter =
          std::find(nList[l2 - 1].begin() + 1, nList[l2 - 1].end(), m_k);
      // If the element has been found, for l2:
      if (it != nList[l1 - 1].end()) {
        l2_neighbour = true;
        compare1 = (*basal1)[3];  // l4
        compare2 = (*basal1)[5];  // l6
      }  // End of check for l2 being a neighbour of m_k
    }    // end of checking to see if l2 is a neighbour or not
    // =================================
    if (l1_neighbour == false && l2_neighbour == false) {
      continue;
    }  // Continue the loop!
    // If l1 or l2 is a neighbour of m_k,
    // Find the next six elements: m_k, m_{k+1}, m_{k+2}, m_{k+3}, m_{k+4},
    // m_{k+5} evenElements will contain m_k, m_{k+2}, m_{k+4} oddElements will
    // contain m_{k+1}, m_{k+3}, m_{k+5}
    evenElements.clear();
    oddElements.clear();
    // Get the next six elements
    for (int m = k; m < k + 6; m++) {
      index = m;
      if (index >= ringSize) {
        index = m - ringSize;
      }
      // For m_k, m_{k+2}, m_{k+4}
      if ((m - k) % 2 == 0) {
        evenElements.push_back((*basal2)[index]);
      }  // end of updating m_k, m_{k+2}, m_{k+4}
      // For m_{k+1}, m_{k+3}, m_{k+5}
      else {
        oddElements.push_back((*basal2)[index]);
      }  // end of updating m_{k+1}, m_{k+3}, m_{k+5}
    }    // end of loop to get even and odd elements; the next six elements
    // =================================
    // Init
    isNeigh = false;
    notNeigh = false;
    // ---------------------------------
    // If m_k is a neighbour of l1:
    if (l1_neighbour) {
      isNeigh = ring::basalNeighbours(nList, &evenElements, compare1, compare2);
      if (isNeigh == false) {
        continue;
      }
      notNeigh = ring::notNeighboursOfRing(nList, &oddElements, basal1);
      if (notNeigh == false) {
        continue;
      }
      // All the conditions are true.
      return true;
    }  // end of checks and conditions for l1_neighbour
    // ---------------------------------
    // If m_k is a neighbour of l2:
    if (l2_neighbour) {
      isNeigh = ring::basalNeighbours(nList, &oddElements, compare1, compare2);
      if (isNeigh == false) {
        continue;
      }
      notNeigh = ring::notNeighboursOfRing(nList, &evenElements, basal1);
      if (notNeigh == false) {
        continue;
      }
      // All the conditions are true
      return true;
    }  // end of checks and conditions for l1_neighbour
    // ---------------------------------
    // =================================
  }  // End of loop (k) through all the elements in basal2
  // ---------------------------------------------
  return false;
}

/********************************************/ /**
                                                *  Tests whether the last two
                                                *elements of a triplet are
                                                *neighbours of two atom IDs
                                                *passed in
                                                ***********************************************/
bool ring::basalNeighbours(std::vector<std::vector<int>> nList,
                           std::vector<int> *triplet, int atomOne,
                           int atomTwo) {
  // Search for needles in a haystack :)
  int needle1 = (*triplet)[1];
  int needle2 = (*triplet)[2];
  bool neighbourFound = false;
  bool neighOne = false,
       neighTwo = false;  // Neighbour for atomOne, or neighbour for atomTwo
  // ----------------------------
  // For first element needle1, which must belong to either atomOne's or
  // atomTwo's neighbour list Search atomOne's neighbours
  auto it = std::find(nList[atomOne - 1].begin() + 1, nList[atomOne - 1].end(),
                      needle1);
  if (it != nList[atomOne - 1].end()) {
    neighbourFound = true;
    neighOne = true;
  }  // atomOne's neighbour
  // If it is not atomOne's neighbour, it might be atomTwo's neighbour
  if (!neighOne) {
    it = std::find(nList[atomTwo - 1].begin() + 1, nList[atomTwo - 1].end(),
                   needle1);
    if (it != nList[atomTwo - 1].end()) {
      neighbourFound = true;
      neighTwo = true;
    }  // end of check to see if neighbour was found
  }    // End of check to see if needle1 is atomTwo's neighbour
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
    it = std::find(nList[atomTwo - 1].begin() + 1, nList[atomTwo - 1].end(),
                   needle2);
    // It is a neighbour
    if (it != nList[atomTwo - 1].end()) {
      return true;
    }
    // It is not a neighbour
    else {
      return false;
    }
  }  // End of check for neighbour of atomTwo
  // ===============
  // if atomTwo was a neighbour of needle1, needle2 must be a neighbour of
  // atomOne
  else {
    it = std::find(nList[atomOne - 1].begin() + 1, nList[atomOne - 1].end(),
                   needle2);
    // It is a neighbour
    if (it != nList[atomOne - 1].end()) {
      return true;
    }
    // It is not a neighbour
    else {
      return false;
    }
  }
  // ===============
}

/********************************************/ /**
 *  Checks to make sure that the elements of the triplet are NOT
 neighbours of any elements inside a vector (ring) passed in (false)
 If any of them are neighbours, this function returns false
 ***********************************************/
bool ring::notNeighboursOfRing(std::vector<std::vector<int>> nList,
                               std::vector<int> *triplet,
                               std::vector<int> *ring) {
  int iatom;  // AtomID of the atom to be searched for inside the neighbour
              // lists
  int jatom;  // AtomID of in whose neighbour list iatom will be searched for
  std::vector<int>::iterator it;

  for (int i = 0; i < (*triplet).size(); i++) {
    iatom = (*triplet)[i];  // AtomID to be searched for
    // iatom must be searched for in the neighbour lists of all elements
    // of the ring vector
    for (int j = 0; j < (*ring).size(); j++) {
      jatom = (*ring)[j];
      // ------------------
      // Search for iatom in the neighbour list of jatom
      it = std::find(nList[jatom - 1].begin() + 1, nList[jatom - 1].end(),
                     iatom);
      // It is a neighbour!
      if (it != nList[jatom - 1].end()) {
        return false;
      }
      // ------------------
    }  // end of loop through every element of ring
  }    // end of loop through all triplet elements

  return true;
}

/********************************************/ /**
                                                *  Finding prismatic rings
                                                ***********************************************/
int ring::findPrismatic(std::vector<std::vector<int>> rings,
                        std::vector<int> *listHC,
                        std::vector<ring::strucType> *ringType, int iring,
                        int jring, std::vector<int> *prismaticRings) {
  int j1, j2, j3;    // Used for making rings to be searched
  int ringSize = 6;  // This is 6 for hexagons
  bool isEqual;      // Tests to see if two rings have the same elements or not
  std::vector<int> iTriplet;    // triplet formed from iring
  std::vector<int> jTriplet;    // triplet formed from jring
  std::vector<int> sortedRing;  // triplet formed from jring
  int j;

  // Make all possible alternate triplets from iring
  for (int k = 0; k < ringSize; k += 2) {
    // Clear vector
    iTriplet.clear();
    // -----------------------
    // Get a vector out of indices in iring and jring
    // Indices for making findRing out of iring and jring
    j1 = k;
    j2 = k + 1;
    j3 = k + 2;
    if (j3 >= ringSize) {
      j3 -= ringSize;
    }
    // Make a triplet out of iring (alternate!)
    // Replace this ugly code later!
    iTriplet.push_back(rings[iring][j1]);
    iTriplet.push_back(rings[iring][j2]);
    iTriplet.push_back(rings[iring][j3]);
    // -----------------------
    // Once the triplet for iring has been made,
    // make every possible triplet out of jring
    // and combine the two to create a new vector
    // Loop through every three to make a triplet
    for (int m = 0; m < ringSize; m++) {
      jTriplet.clear();  // Clear the triplet
      // ------
      // Get a triplet
      for (int i = m; i < m + 3; i++) {
        j = i;
        if (i >= ringSize) {
          j = i - ringSize;
        }
        jTriplet.push_back(rings[jring][j]);
      }  // end of getting a triplet from k
      // ------
      // Append iTriplet to jTriplet
      jTriplet.reserve(jTriplet.size() + iTriplet.size());
      jTriplet.insert(jTriplet.end(), iTriplet.begin(), iTriplet.end());
      // ------
      // Now search through all the rings to find a match for
      // jTriplet, which contains a triplet from iring and jring
      for (int kring = 0; kring < rings.size(); kring++) {
        //
        if (kring == iring || kring == jring) {
          continue;
        }
        //
        // Test to see if findRing and kring are equal or not
        //
        isEqual = ring::compareRings(jTriplet, rings[kring]);
        //
        // If they are equal, this is a prismatic ring!
        if (isEqual) {
          (*listHC).push_back(kring);          // Update listHC vector
          (*prismaticRings).push_back(kring);  // Update prismatic rings
          // Update the type inside ringType
          // If the ring is already a DDC ring, it is a mixed ring
          if ((*ringType)[kring] == ring::DDC) {
            (*ringType)[kring] = ring::bothPrismatic;
          }
          // If it is unclassified, it is just a prismatic ring
          if ((*ringType)[kring] == ring::unclassified) {
            (*ringType)[kring] = ring::HCprismatic;
          }
          break;
        }  // end of check for prismatic ring
      }    // end of search for match in rings for jTriplet
      // ------
    }  // end of loop through every three for jring
    // -----------------------
  }  // end of triplet making for iring

  return 0;
}

/********************************************/ /**
 *  Checks to see if two vectors (ring1, ring2) have the same
 elements (disordered). So the sequence is not important here.
 Returns true if the rings have the same elements
 ***********************************************/
bool ring::compareRings(std::vector<int> ring1, std::vector<int> ring2) {
  // Sort the rings first
  sort(ring1.begin(), ring1.end());  // Sort ring1 by ID
  sort(ring2.begin(), ring2.end());  // Sort ring2 by ID
  bool result;

  (ring1 == ring2) ? result = true : result = false;

  return result;
}

/********************************************/ /**
 *  Determines which hexagonal rings are both DDCs and HCs. This function
 returns a vector which contains the ring IDs of all the rings which are both.
 The ring IDs correspond to the index of the rings inside the vector of vector
 rings, starting from 0. Rings which are both are of enum type bothBasal OR
 bothPrismatic. Reassign rings which are mixed in listDDC and listHC as the
 dummy value -10
 ***********************************************/
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
    if ((*ringType)[iring] == ring::bothBasal ||
        (*ringType)[iring] == ring::bothPrismatic) {
      listMixed.push_back(iring);

      //-----------------
      // Search for iring in listDDC
      std::sort((*listDDC).begin(), (*listDDC).end());
      auto iter = std::find((*listDDC).begin(), (*listDDC).end(), iring);
      if (iter != (*listDDC).end()) {
        *iter = dummyValue;  // Assign dummy value to the mixed ring
      }                      // found in listDDC
      //-----------------
      //-----------------
      // Search for iring in listHC
      std::sort((*listHC).begin(), (*listHC).end());
      auto itr = std::find((*listHC).begin(), (*listHC).end(), iring);
      if (itr != (*listHC).end()) {
        *itr = dummyValue;  // Assign dummy value to the mixed ring
      }                     // found in listHC
      //-----------------

    }  // end of check for type
  }    // end of loop through all every ring

  return listMixed;
}

/********************************************/ /**
 *  Determines which complete cages are 'mixed'. For a cage to be classified as
 mixed, all rings must be of enum type bothBasal OR bothPrismatic.
 ***********************************************/
int ring::findMixedCages(std::vector<strucType> *ringType,
                         std::vector<cage::Cage> *cageList, int *numDDC,
                         int *numHC, int *numMC) {
  int iring;        // Ring index (starting from 0) comprising a particular cage
  int nrings;       // Number of rings in a particular cage
  int nMixedRings;  // Number of mixed rings
  cage::cageType currentType;  // Type of icage

  // Init
  *numDDC = 0;
  *numHC = 0;
  *numMC = 0;

  // Loop through every cage in cageList
  for (int icage = 0; icage < (*cageList).size(); icage++) {
    nrings = (*cageList)[icage].rings.size();
    nMixedRings = 0;  // init to zero
    // Loop through every ring inside icage
    for (int i = 0; i < nrings; i++) {
      iring = (*cageList)[icage].rings[i];  // Ring index of ith ring in icage
      // Check if iring is a mixed ring or not
      if ((*ringType)[iring] == ring::bothBasal ||
          (*ringType)[iring] == ring::bothPrismatic) {
        nMixedRings++;
      }  // end of check to see if iring is mixed
    }    // end of loop through every ring in icage
    // Check to see if all the rings in icage are mixed or not
    if (nMixedRings == nrings) {
      (*cageList)[icage].type = cage::Mixed;
      *numMC += 1;  // Add to the number of mixed cages
    }               // icage is a mixed cage
    else {
      // Add to the count of DDCs and HCs
      currentType = (*cageList)[icage].type;
      if (currentType == cage::DoubleDiaC) {
        *numDDC += 1;
      } else if (currentType == cage::HexC) {
        *numHC += 1;
      } else {
        std::cerr << "The cage is the wrong type!\n";
      }
    }
  }  // end of loop through every cage

  return 0;
}