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

#include <ring.hpp>

/**
 * @details Deletes the memory of a
 * vector of vectors (specifically for rings).
 * @param[in, out] rings The vector of vectors to be cleared.
 */
int ring::clearRingList(std::vector<std::vector<int>> &rings) {
  //
  std::vector<std::vector<int>> tempEmpty;

  rings.swap(tempEmpty);

  return 0;
}

/**
 * @details Assign an atomType (equal to the number of nodes in the ring) given
 * the rings vector.
 * @param[in] rings The vector of vectors containing the primitive rings, of a
 *  particular ring size.
 * @param[in, out] atomTypes A vector which contains a type for each atom,
 *  depending on it's type as classified by the prism identification scheme.
 * @param[in] nRings Number of rings.
 */
int ring::assignPolygonType(std::vector<std::vector<int>> rings,
                            std::vector<int> *atomTypes,
                            std::vector<int> nRings) {
  // Every value in listPrism corresponds to an index in rings.
  // Every ring contains atom indices, corresponding to the indices (not atom
  // IDs) in rings
  int iring;        // Index of current ring
  int iatom;        // Index of current atom
  int ringSize;     // Ring size of the current ring
  int prevRingSize; // Ring size previously assigned to a point

  // Dummy value corresponds to a value of 1.
  // If an atom is shared by more than one ring type, it is assigned the
  // value 2.

  // Loop through every ring in rings
  for (int iring = 0; iring < rings.size(); iring++) {
    ringSize = rings[iring].size();
    // Loop through every element in iring
    for (int j = 0; j < ringSize; j++) {
      iatom = rings[iring][j]; // Atom index
      // Update the atom type
      if ((*atomTypes)[iatom] == 1) {
        (*atomTypes)[iatom] = ringSize;
      } // The atom is unclassified
      else {
        // Only update the ring type if the number is higher
        prevRingSize = (*atomTypes)[iatom]; // Previously assigned ring size
        if (ringSize > prevRingSize) {
          (*atomTypes)[iatom] = ringSize;
        } // end of assigning the new ring size
      }   // only update if the number is higher
    }     // end of loop through every atom in iring
  }       // end of loop through every ring

  return 0;
} // end of function

/**
 *  @details Function that finds and returns a vector containing the elements
 *   shared by two input rings (each ring is a vector).
 *  @param[in] ring1 The first ring.
 *  @param[in] ring2 The second ring.
 *  @return A vector containing the common elements between the input rings.
 */
std::vector<int> ring::findsCommonElements(std::vector<int> ring1,
                                           std::vector<int> ring2) {
  //
  std::vector<int> common;
  int iatom; // Index to search for

  for (int i = 0; i < ring1.size(); i++) {
    iatom = ring1[i];
    // Search for iatom in ring2

    auto it = std::find(ring2.begin(), ring2.end(), iatom);

    if (it != ring2.end()) {
      common.push_back(iatom);
    } // iatom was found!
  }   // end of loop through every element of ring1

  return common;
}

/**
 * @details Function that finds the common elements
 *  in three input rings
 * @param[in] ring1 The first ring.
 * @param[in] ring2 The second ring.
 * @param[in] ring3 The third ring.
 * @return A value which is true if the three rings have at least one common
 * element, and false if the three rings have no elements in common.
 */
bool ring::commonElementsInThreeRings(std::vector<int> ring1,
                                      std::vector<int> ring2,
                                      std::vector<int> ring3) {
  std::vector<int>
      common1; // Vector containing the common elements of the first two rings
  std::vector<int>
      common2; // Vector containing the common elements of the three rings

  // Common elements among the first two rings
  common1 = ring::findsCommonElements(ring1, ring2);
  if (common1.size() == 0) {
    return false;
  } // no common elements in ring1 and ring2

  // Common elements among all three rings
  common2 = ring::findsCommonElements(common1, ring3);

  // If no common elements were found:
  if (common2.size() == 0) {
    return false;
  } // no common elements between ring1, ring2, and ring3

  return true; // Common elements found!
}

/**
 * @details Function that finds out if a given triplet is
 *  is present (in the same order or in reversed order) in a given ring.
 * @param[in] ring The input ring containing the indices of atoms.
 * @param[in] triplet Vector containing the triplet, for whose presence the
 *  input ring vector will be checked.
 * @return A bool value which is true if the triplet is present in the ring,
 *  and is false if the triplet is not in the ring.
 */
bool ring::findTripletInRing(std::vector<int> ring, std::vector<int> triplet) {
  //
  int ringSize = ring.size();   // should be 6
  std::vector<int> ringTriplet; // triplet from the ring to be searched
  int kIndex;                   // Used for making the triplet

  // Loop through every possible triplet in the ring to be searched
  for (int i = 0; i < ringSize; i++) {
    ringTriplet.clear(); // init
    // Get the first element of the ring triplet
    ringTriplet.push_back(ring[i]);
    //
    // Get the next two elements
    for (int k = 1; k < 3; k++) {
      kIndex = i + k;
      // Wrap-around
      if (kIndex >= ringSize) {
        kIndex -= ringSize;
      } // end of wrap-around
      ringTriplet.push_back(ring[kIndex]);
    } // next two elements of ringTriplet
    //
    // Obtained ringTriplet!
    // Check equality
    if (triplet == ringTriplet) {
      return true;
    } // triplet matches!
    //
    // Check the reversed triplet too
    std::reverse(ringTriplet.begin(), ringTriplet.end());

    // Check for equality
    if (triplet == ringTriplet) {
      return true;
    } // reversed triplet matches
  }   // first element of ringTriplet

  return false;
}

/**
 * @details Function that gets rings of a single ring size (i.e. a particular
 * number of nodes) from all primitive rings, and returns a vector of vectors
 * containing the rings of the specified size.
 * @param[in] rings The vector of vectors containing the primitive rings of all
 *  sizes.
 * @param[in] ringSize The desired ring size or number of nodes in each ring to
 *  be saved.
 * @return A vector of vectors containing primitive rings of one ring size or
 *  length.
 */
std::vector<std::vector<int>>
ring::getSingleRingSize(std::vector<std::vector<int>> rings, int ringSize) {
  //
  std::vector<std::vector<int>> ringSingleSize; // rings of one size

  // rings contains primitive rings of all sizes
  // Only save rings of a given size (ringSize) to the new
  // vector of vectors, ringSingleSize
  for (int iring = 0; iring < rings.size(); iring++) {
    // Check the size of the current ring
    // If it is the correct size, save it in ringSingleSize
    if (rings[iring].size() == ringSize) {
      ringSingleSize.push_back(rings[iring]);
    } // End of check of the size of iring
  }   // end of loop through all rings in rings

  return ringSingleSize;
}

/**
 * @details For two vectors, checks to see if there are common elements (true)
 *  or not (false).
 * @param[in] ring1 The vector of the first ring.
 * @param[in] ring2 The vector of the second ring.
 * @return A bool which is true if the input vectors have at least one common
 * element, and false if there are no common elements.
 */
bool ring::hasCommonElements(std::vector<int> ring1, std::vector<int> ring2) {
  std::vector<int> commonElements; // Vector containing common elements

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

/**
 * @details Checks to see if two vectors (ring1, ring2) have the same
 *  elements (which may or may not be disordered). The order or sequence of
 *  elements is not important, so the rings are sorted. Returns true if the
 *  rings have the same elements
 * @param[in] ring1 The first ring.
 * @param[in] ring2 The second ring.
 * @return A bool; true if the rings contain the same elements (not necessarily
 *  in the same sequence) and false if they do not have the same elements.
 */
bool ring::compareRings(std::vector<int> ring1, std::vector<int> ring2) {
  // Sort the rings first
  sort(ring1.begin(), ring1.end()); // Sort ring1 by ID
  sort(ring2.begin(), ring2.end()); // Sort ring2 by ID
  bool result;

  (ring1 == ring2) ? result = true : result = false;

  return result;
}
