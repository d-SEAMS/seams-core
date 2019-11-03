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
                                                *  Finds and returns all the
                                                *elements common between the two
                                                *vectors
                                                ***********************************************/
std::vector<int> ring::findsCommonElements(std::vector<int> ring1,
                                           std::vector<int> ring2) {
  //
  std::vector<int> common;
  int iatom;  // Index to search for

  for (int i = 0; i < ring1.size(); i++) {
    iatom = ring1[i];
    // Search for iatom in ring2

    auto it = std::find(ring2.begin(), ring2.end(), iatom);

    if (it != ring2.end()) {
      common.push_back(iatom);
    }  // iatom was found!
  }    // end of loop through every element of ring1

  return common;
}

/********************************************/ /**
                                                *  Finds the common elements in
                                                *three rings
                                                ***********************************************/
bool ring::commonElementsInThreeRings(std::vector<int> ring1,
                                      std::vector<int> ring2,
                                      std::vector<int> ring3) {
  //
  std::vector<int>
      common1;  // Vector containing the common elements of the first two rings
  std::vector<int>
      common2;  // Vector containing the common elements of the three rings

  // Common elements among the first two rings
  common1 = ring::findsCommonElements(ring1, ring2);
  if (common1.size() == 0) {
    return false;
  }  // no common elements in ring1 and ring2

  // Common elements among all three rings
  common2 = ring::findsCommonElements(common1, ring3);

  // If no common elements were found:
  if (common2.size() == 0) {
    return false;
  }  // no common elements between ring1, ring2, and ring3

  return true;  // Common elements found!
}

/********************************************/ /**
                                                *  Finds out if a triplet is
                                                *present (in the same order or
                                                *reversed) in a given ring
                                                ***********************************************/
bool ring::findTripletInRing(std::vector<int> ring, std::vector<int> triplet) {
  //
  int ringSize = ring.size();    // should be 6
  std::vector<int> ringTriplet;  // triplet from the ring to be searched
  int kIndex;                    // Used for making the triplet

  // Loop through every possible triplet in the ring to be searched
  for (int i = 0; i < ringSize; i++) {
    ringTriplet.clear();  // init
    // Get the first element of the ring triplet
    ringTriplet.push_back(ring[i]);
    //
    // Get the next two elements
    for (int k = 1; k < 3; k++) {
      kIndex = i + k;
      // Wrap-around
      if (kIndex >= ringSize) {
        kIndex -= ringSize;
      }  // end of wrap-around
      ringTriplet.push_back(ring[kIndex]);
    }  // next two elements of ringTriplet
    //
    // Obtained ringTriplet!
    // Check equality
    if (triplet == ringTriplet) {
      return true;
    }  // triplet matches!
    //
    // Check the reversed triplet too
    std::reverse(ringTriplet.begin(), ringTriplet.end());

    // Check for equality
    if (triplet == ringTriplet) {
      return true;
    }  // reversed triplet matches
  }    // first element of ringTriplet

  return false;
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