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