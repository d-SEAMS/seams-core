#ifndef __TOPO_ONE_DIM_H_
#define __TOPO_ONE_DIM_H_

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

#include <mol_sys.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

namespace ring {

// Find out which rings are prisms.
// Returns a vector containing all the ring IDs which are prisms
std::vector<int> findPrisms(
    std::vector<std::vector<int>> rings, std::vector<strucType> *ringType,
    int *nPrisms, std::vector<std::vector<int>> nList,
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

}  // namespace ring

#endif  // __TOPOCONFINED_H_
