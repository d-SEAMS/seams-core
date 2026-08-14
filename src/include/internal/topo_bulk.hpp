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

#ifndef SEAMS_TOPO_BULK_H_
#define SEAMS_TOPO_BULK_H_

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cmath>
#include <memory>
#include <sstream>
#include <string>
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

//! Find out rings in the bulk, looping through all ring sizes upto the
//! maxDepth The input ringsAllSizes array has rings of every size.
[[nodiscard]] int bulkPolygonRingAnalysis(
    std::string path, const std::vector<std::vector<int>> &rings,
    const std::vector<std::vector<int>> &nList,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, int maxDepth,
    int firstFrame);

// DDC HC Ring functions

//! Find out which rings are DDCs or HCs, which are comprised of 6-membered
//! primitive rings. Start with a neighbour list (by index) and a vector of
//! vectors of rings (also by index). TODO: try 'square' ice and ice0
[[nodiscard]] int topoBulkAnalysis(std::string path, const std::vector<std::vector<int>> &rings,
                     const std::vector<std::vector<int>> &nList,
                     molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     int firstFrame, bool onlyTetrahedral = true);

/** @struct ring::RingSearchIndex
 * @brief Inverted index from an atom to the rings that contain it.
 * @details The cage searches repeatedly ask which rings a given atom belongs
 *  to. Answering that by scanning every ring makes the searches quadratic in
 *  the ring count. One pass over the ring network answers it in constant time
 *  thereafter, and the ring indices in each row stay ascending, which is the
 *  order the searches previously visited them in.
 */
struct RingSearchIndex {
  std::vector<std::vector<int>> ringsContainingAtom;
};

//! Builds the inverted atom-to-rings index used by the cage searches
[[nodiscard]] RingSearchIndex
buildRingSearchIndex(const std::vector<std::vector<int>> &rings, int numAtoms);

//! Find out which hexagonal rings are DDC (Double Diamond Cages) rings.
//! Returns a vector containing all the ring IDs which are DDC rings
std::vector<int> findDDC(const std::vector<std::vector<int>> &rings,
                         std::vector<strucType> &ringType,
                         const std::vector<int> &listHC,
                         std::vector<cage::Cage> &cageList);

//! As findDDC, reusing an index the caller has already built
std::vector<int> findDDC(const std::vector<std::vector<int>> &rings,
                         std::vector<strucType> &ringType,
                         const std::vector<int> &listHC,
                         std::vector<cage::Cage> &cageList,
                         const RingSearchIndex &index);

//! Find out which hexagonal rings are both DDCs (Double Diamond Cages) and HCs
//! (Hexagonal Cages). Returns a vector containing all the ring IDs which are
//! of this type
std::vector<int> findMixedRings(const std::vector<std::vector<int>> &rings,
                                std::vector<strucType> &ringType,
                                std::vector<int> &listDDC,
                                std::vector<int> &listHC);

//! Find out which hexagonal rings are HC rings.
//! Returns a vector containing all the ring IDs which are HC rings
std::vector<int> findHC(const std::vector<std::vector<int>> &rings,
                        std::vector<strucType> &ringType,
                        const std::vector<std::vector<int>> &nList,
                        std::vector<cage::Cage> &cageList);

//! As findHC, reusing an index the caller has already built
std::vector<int> findHC(const std::vector<std::vector<int>> &rings,
                        std::vector<strucType> &ringType,
                        const std::vector<std::vector<int>> &nList,
                        std::vector<cage::Cage> &cageList,
                        const RingSearchIndex &index);

//! First condition for the DDC: There must be at least 3 other
//! rings in which each element of the equatorial  ring is present
bool conditionOneDDC(const std::vector<std::vector<int>> &rings,
                     std::vector<int> &peripheralRings, int iring);

//! As conditionOneDDC, reusing an index the caller has already built
bool conditionOneDDC(const std::vector<std::vector<int>> &rings,
                     std::vector<int> &peripheralRings, int iring,
                     const RingSearchIndex &index);

//! Second condition for the DDC: There must be at least 1 other
//! ring for every triplet in the equatorial  ring
bool conditionTwoDDC(const std::vector<std::vector<int>> &rings,
                     std::vector<int> &peripheralRings, int iring);

//! As conditionTwoDDC, answering each triplet from the inverted index
bool conditionTwoDDC(const std::vector<std::vector<int>> &rings,
                     std::vector<int> &peripheralRings, int iring,
                     const RingSearchIndex &index);

//! Third condition for the DDC: Even (by vector index) numbered index triplets
//! and odd triplets must have at least one element in common
bool conditionThreeDDC(const std::vector<std::vector<int>> &rings,
                       std::vector<int> &peripheralRings);

//! Tests whether two rings are basal rings (true) or not (false)
bool basalConditions(const std::vector<std::vector<int>> &nList,
                     std::vector<int> &basal1, std::vector<int> &basal2);

//! Tests whether the last two elements of a triplet are neighbours of two atom
//! IDs passed in
bool basalNeighbours(const std::vector<std::vector<int>> &nList,
                     std::vector<int> &triplet, int atomOne, int atomTwo);

//! Tests to check that elements of a triplet are not neighbours of a ring
//! (vector) passed
bool notNeighboursOfRing(const std::vector<std::vector<int>> &nList,
                         std::vector<int> &triplet, std::vector<int> &ring);

//! Finds the prismatic rings from basal rings iring and jring
[[nodiscard]] int findPrismatic(const std::vector<std::vector<int>> &rings, std::vector<int> &listHC,
                  std::vector<strucType> &ringType, int iring, int jring,
                  std::vector<int> &prismaticRings);

//! As findPrismatic, drawing candidate rings from the inverted index
[[nodiscard]] int findPrismatic(const std::vector<std::vector<int>> &rings,
                                std::vector<int> &listHC,
                                std::vector<strucType> &ringType, int iring,
                                int jring, std::vector<int> &prismaticRings,
                                const RingSearchIndex &index);

//! Assigns a type of enum class iceType, to every atom, using information from
//! ringType, which has the information of every ring
[[nodiscard]] int getAtomTypesTopoBulk(const std::vector<std::vector<int>> &rings,
                         const std::vector<ring::strucType> &ringType,
                         std::vector<cage::iceType> &atomTypes);

//! Determines the number of HCs, DDCs, Mixed rings, prismatic and basal rings
[[nodiscard]] int getStrucNumbers(const std::vector<ring::strucType> &ringType,
                    const std::vector<cage::Cage> &cageList, int &numHC, int &numDDC,
                    int &mixedRings, int &prismaticRings, int &basalRings);

} // namespace ring

/**
 *  @addtogroup prism3
 *  @{
 */

namespace prism3 {

//! Find out which rings are prisms.
[[nodiscard]] int findBulkPrisms(const std::vector<std::vector<int>> &rings,
                   std::vector<ring::strucType> &ringType,
                   const std::vector<std::vector<int>> &nList,
                   molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                   std::vector<double> &rmsdPerAtom, double heightCutoff = 8);

//! Tests whether two rings are basal rings (true) or not (false) for a prism
//! (strict criterion)
bool basalPrismConditions(const std::vector<std::vector<int>> &nList,
                          std::vector<int> &basal1, std::vector<int> &basal2);

//! Reduced criterion: Two candidate basal rings of a prism block should have
//! at least one bond between them
bool relaxedPrismConditions(const std::vector<std::vector<int>> &nList,
                            std::vector<int> &basal1, std::vector<int> &basal2);

//! Check to see that candidate basal prisms are not really far from each other
bool basalRingsSeparation(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &basal1, const std::vector<int> &basal2, double heightCutoff = 8);
} // namespace prism3

#endif // SEAMS_TOPO_BULK_H_
