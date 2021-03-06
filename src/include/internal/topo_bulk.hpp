//-----------------------------------------------------------------------------------
// d-SEAMS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef __TOPO_BULK_H_
#define __TOPO_BULK_H_

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
#include <order_parameter.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>
#include <shapeMatch.hpp>

/** @file topo_bulk.hpp
 *  @brief File containing functions used specific to bulk topological network
 * critera.
 */

/**
 *  @addtogroup ring
 *  @{
 */

namespace ring {

//! Find out which rings are DDCs or HCs, which are comprised of 6-membered
//! primitive rings. Start with a neighbour list (by index) and a vector of
//! vectors of rings (also by index). TODO: try 'square' ice and ice0
int topoBulkAnalysis(std::string path, std::vector<std::vector<int>> rings,
                     std::vector<std::vector<int>> nList,
                     molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                     int firstFrame, bool onlyTetrahedral = true);

//! Find out which hexagonal rings are DDC (Double Diamond Cages) rings.
//! Returns a vector containing all the ring IDs which are DDC rings
std::vector<int> findDDC(std::vector<std::vector<int>> rings,
                         std::vector<strucType> *ringType,
                         std::vector<int> listHC,
                         std::vector<cage::Cage> *cageList);

//! Find out which hexagonal rings are both DDCs (Double Diamond Cages) and HCs
//! (Hexagonal Cages). Returns a vector containing all the ring IDs which are
//! of this type
std::vector<int> findMixedRings(std::vector<std::vector<int>> rings,
                                std::vector<strucType> *ringType,
                                std::vector<int> *listDDC,
                                std::vector<int> *listHC);

//! Find out which hexagonal rings are HC rings.
//! Returns a vector containing all the ring IDs which are HC rings
std::vector<int> findHC(std::vector<std::vector<int>> rings,
                        std::vector<strucType> *ringType,
                        std::vector<std::vector<int>> nList,
                        std::vector<cage::Cage> *cageList);

//! First condition for the DDC: There must be at least 3 other
//! rings in which each element of the equatorial  ring is present
bool conditionOneDDC(std::vector<std::vector<int>> rings,
                     std::vector<int> *peripheralRings, int iring);

//! Second condition for the DDC: There must be at least 1 other
//! ring for every triplet in the equatorial  ring
bool conditionTwoDDC(std::vector<std::vector<int>> rings,
                     std::vector<int> *peripheralRings, int iring);

//! Third condition for the DDC: Even (by vector index) numbered index triplets
//! and odd triplets must have at least one element in common
bool conditionThreeDDC(std::vector<std::vector<int>> rings,
                       std::vector<int> *peripheralRings);

//! Tests whether two rings are basal rings (true) or not (false)
bool basalConditions(std::vector<std::vector<int>> nList,
                     std::vector<int> *basal1, std::vector<int> *basal2);

//! Tests whether the last two elements of a triplet are neighbours of two atom
//! IDs passed in
bool basalNeighbours(std::vector<std::vector<int>> nList,
                     std::vector<int> *triplet, int atomOne, int atomTwo);

//! Tests to check that elements of a triplet are not neighbours of a ring
//! (vector) passed
bool notNeighboursOfRing(std::vector<std::vector<int>> nList,
                         std::vector<int> *triplet, std::vector<int> *ring);

//! Finds the prismatic rings from basal rings iring and jring
int findPrismatic(std::vector<std::vector<int>> rings, std::vector<int> *listHC,
                  std::vector<strucType> *ringType, int iring, int jring,
                  std::vector<int> *prismaticRings);

//! Assigns a type of enum iceType, to every atom, using information from
//! ringType, which has the information of every ring
int getAtomTypesTopoBulk(std::vector<std::vector<int>> rings,
                         std::vector<ring::strucType> ringType,
                         std::vector<cage::iceType> *atomTypes);

//! Determines the number of HCs, DDCs, Mixed rings, prismatic and basal rings
int getStrucNumbers(std::vector<ring::strucType> ringType,
                    std::vector<cage::Cage> cageList, int *numHC, int *numDDC,
                    int *mixedRings, int *prismaticRings, int *basalRings);

} // namespace ring

/**
 *  @addtogroup prism3
 *  @{
 */

namespace prism3 {

//! Find out which rings are prisms.
int findBulkPrisms(std::vector<std::vector<int>> rings,
                   std::vector<ring::strucType> *ringType,
                   std::vector<std::vector<int>> nList,
                   molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                   std::vector<double> *rmsdPerAtom, double heightCutoff = 8);

//! Tests whether two rings are basal rings (true) or not (false) for a prism
//! (strict criterion)
bool basalPrismConditions(std::vector<std::vector<int>> nList,
                          std::vector<int> *basal1, std::vector<int> *basal2);

//! Reduced criterion: Two candidate basal rings of a prism block should have
//! at least one bond between them
bool relaxedPrismConditions(std::vector<std::vector<int>> nList,
                            std::vector<int> *basal1, std::vector<int> *basal2);

//! Check to see that candidate basal prisms are not really far from each other
bool basalRingsSeparation(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2, double heightCutoff = 8);
} // namespace prism3

#endif // __TOPO_BULK_H_
