#ifndef __RINGS_H_
#define __RINGS_H_

#include <math.h>
#include <sys/stat.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <cage.hpp>
#include <mol_sys.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

namespace ring {

// General enum used throughout this program (Prism is for our prism
// classification) {per-ring classification}
enum strucType {
  unclassified,
  DDC,
  HCbasal,
  HCprismatic,
  bothBasal,
  bothPrismatic,
  Prism
};

// Returns a vector of vectors of rings of a single size
std::vector<std::vector<int>> getSingleRingSize(
    std::vector<std::vector<int>> rings, int ringSize);

// Check to see if two vectors have common elements or not
// True, if common elements are present and false if there are no common
// elements
bool hasCommonElements(std::vector<int> ring1, std::vector<int> ring2);

// Compares two disordered vectors and checks to see if they contain the same
// elements
bool compareRings(std::vector<int> ring1, std::vector<int> ring2);

// Erases memory for a vector of vectors for a list of rings
int clearRingList(std::vector<std::vector<int>> &rings);

}  // namespace ring

#endif  // __RINGS_H_
