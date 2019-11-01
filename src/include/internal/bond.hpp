#ifndef __BONDING_H_
#define __BONDING_H_

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

// Internal
#include <cage.hpp>
#include <mol_sys.hpp>
#include <seams_input.hpp>

namespace bond {

// Create a vector of vectors (similar to the neighbour list conventions)
// containing the hydrogen bond connectivity information. Decides the existence
// of the hydrogen bond depending on the O--O and O--H vectors from the
// neighbour list already constructed
std::vector<std::vector<int>> populateHbonds(
    std::string filename,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, int targetFrame, int Htype);

// Calculates the distance of the hydrogen bond between O and H (of different
// atoms), given the respective pointClouds and the indices to each atom
double getHbondDistanceOH(
    molSys::PointCloud<molSys::Point<double>, double> *oCloud,
    molSys::PointCloud<molSys::Point<double>, double> *hCloud, int oAtomIndex,
    int hAtomIndex);

// Create a vector of vectors containing bond connectivity information. May
// contain duplicates! Gets the bond information from the vector of vectors
// containing the rings
std::vector<std::vector<int>> populateBonds(
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud);

// Creates a vector of vectors containing bond connectivity information from the
// rings vector of vectors and cage information
std::vector<std::vector<int>> createBondsFromCages(
    std::vector<std::vector<int>> rings, std::vector<cage::Cage> *cageList,
    cage::cageType type, int *nRings);

// Remove duplicate bonds
std::vector<std::vector<int>> trimBonds(std::vector<std::vector<int>> bonds);

// Remove bonds which are diagonals within rings
int rmDiagBonds(std::vector<std::vector<int>> rings,
                std::vector<std::vector<int>> bonds, std::vector<bool> *flag);

// Remove bonds which are longer than the cutoff
int rmLongBonds(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                std::vector<std::vector<int>> bonds, std::vector<bool> *flag,
                double cutoff);

// Search and set flag to false if a given 'bond' matches a bond in the vector
// of vectors
int searchBondMatch(std::vector<int> matchBond,
                    std::vector<std::vector<int>> bonds,
                    std::vector<bool> *flag);
}  // namespace bond

#endif  // __BONDING_H_
