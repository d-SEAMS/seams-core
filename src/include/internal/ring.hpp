#ifndef __RINGS_H_
#define __RINGS_H_

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <math.h>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
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

// Function for determining whether a ring has more than
// 3 consecutive water molecules or not.
bool checkRing(std::vector<int> ring, std::vector<std::vector<int>> nList);

// Compare a quadruplet with the neighbour list
bool compareQuad(std::vector<int> quad, std::vector<std::vector<int>> nList);

// Find out which hexagonal rings are DDC (Double Diamond Cages) rings.
// Returns a vector containing all the ring IDs which are DDC rings
std::vector<int> findDDC(std::vector<std::vector<int>> rings,
                         std::vector<strucType> *ringType,
                         std::vector<cage::Cage> *cageList);

// Find out which hexagonal rings are both DDCs (Double Diamond Cages) and HCs
// (Hexagonal Cages). Returns a vector containing all the ring IDs which are of
// this type
std::vector<int> findMixedRings(std::vector<std::vector<int>> rings,
                                std::vector<strucType> *ringType,
                                std::vector<int> *listDDC,
                                std::vector<int> *listHC);

// Find out which hexagonal rings are HC rings.
// Returns a vector containing all the ring IDs which are HC rings
std::vector<int> findHC(std::vector<std::vector<int>> rings,
                        std::vector<strucType> *ringType,
                        std::vector<std::vector<int>> nList,
                        std::vector<cage::Cage> *cageList);

// Find out which cages are mixed cages
int findMixedCages(std::vector<strucType> *ringType,
                   std::vector<cage::Cage> *cageList, int *numDDC, int *numHC,
                   int *numMC);

// Find out which rings are prisms.
// Returns a vector containing all the ring IDs which are prisms
std::vector<int>
findPrisms(std::vector<std::vector<int>> rings,
           std::vector<strucType> *ringType, int *nPrisms,
           std::vector<std::vector<int>> nList,
           molSys::PointCloud<molSys::Point<double>, double> *yCloud);

// Tests whether two rings are basal rings (true) or not (false) for a prism
// (strict criterion)
bool basalPrismConditions(std::vector<std::vector<int>> nList,
                          std::vector<int> *basal1, std::vector<int> *basal2);

// Reduced criterion: Two candidate basal rings of a prism block should have at
// least one bond between them
bool relaxedPrismConditions(std::vector<std::vector<int>> nList,
                            std::vector<int> *basal1, std::vector<int> *basal2);

// Checks whether two 4-membered rings are parallel in one dimension or not to
// prevent overcounting
bool discardExtraTetragonBlocks(
    std::vector<int> *basal1, std::vector<int> *basal2,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud);

// First condition for the DDC: There must be at least 3 other
// rings in which each element of the equatorial  ring is present
bool conditionOneDDC(std::vector<std::vector<int>> rings,
                     std::vector<int> *peripheralRings, int iring);

// Second condition for the DDC: There must be at least 1 other
// ring for every triplet in the equatorial  ring
bool conditionTwoDDC(std::vector<std::vector<int>> rings,
                     std::vector<int> *peripheralRings, int iring);

// Third condition for the DDC: Even (by vector index) numbered index triplets
// and odd triplets must have at least one element in common
bool conditionThreeDDC(std::vector<std::vector<int>> rings,
                       std::vector<int> *peripheralRings, int iring);

// Check to see if two vectors have common elements or not
// True, if common elements are present and false if there are no common
// elements
bool hasCommonElements(std::vector<int> ring1, std::vector<int> ring2);

// Tests whether two rings are basal rings (true) or not (false)
bool basalConditions(std::vector<std::vector<int>> nList,
                     std::vector<int> *basal1, std::vector<int> *basal2);

// Tests whether the last two elements of a triplet are neighbours of two atom
// IDs passed in
bool basalNeighbours(std::vector<std::vector<int>> nList,
                     std::vector<int> *triplet, int atomOne, int atomTwo);

// Tests to check that elements of a triplet are not neighbours of a ring
// (vector) passed
bool notNeighboursOfRing(std::vector<std::vector<int>> nList,
                         std::vector<int> *triplet, std::vector<int> *ring);

// Finds the prismatic rings from basal rings iring and jring
int findPrismatic(std::vector<std::vector<int>> rings, std::vector<int> *listHC,
                  std::vector<strucType> *ringType, int iring, int jring,
                  std::vector<int> *prismaticRings);

// Compares two disordered vectors and checks to see if they contain the same
// elements
bool compareRings(std::vector<int> ring1, std::vector<int> ring2);
} // namespace ring

#endif // __RINGS_H_
