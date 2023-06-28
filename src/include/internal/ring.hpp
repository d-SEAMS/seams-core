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

/** @file ring.hpp
 *  @brief File containing common functions used by bulk and confined
 * topological network critera.
 */

/**
 *  @addtogroup ring
 *  @{
 */

/** @brief Topological network criteria functions.
 *  @details This namespace contains functions for the topological network
 * criteria for bulk and confined systems both.
 *
 *  Although the namespace is shared by bulk and confined topological network
 *  criteria, the functions are split into files specific to each set of
 * criteria.
 *
 *  All the topological network criteria are based on the identification of
 *  primitive rings (using the Franzblau algorithm for primitive rings),
 * following which specific rules for combinations of primitive rings are used
 * to identify quasi-one-dimensional, quasi-two-dimensional and bulk ices.
 *
 *   ### Changelog ###
 *
 *  - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 13, 2019
 *  - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

namespace ring {

// General enum class used throughout this program (Prism is for our prism
// classification) {per-ring classification}
/** @enum class ring::strucType Qualifier for each ring, based on the classification
 * type determined by the bulk or confined ice topological network criteria.
 * @var ring::strucType unclassified
 * @brief The ring is unclassified, which may be either water or a deformed type
 * which cannot be classified by the criteria.
 * @var ring::strucType DDC
 * @brief The ring belongs to a double-diamond cage (DDC).
 *
 * @var ring::strucType HCbasal
 * @brief The ring belongs only to a hexagonal cage (HC). Specifically, the ring
 * is purely a basal ring of an HC.
 *
 * @var ring::strucType HCprismatic
 * @brief The ring belongs only to a hexagonal cage (HC); specifically the ring
 * is purely a prismatic ring of an HC. It is not shared by a DDC.
 *
 * @var ring::strucType bothBasal
 * @brief The ring belongs to both a DDC and HC. It is a 'mixed' ring. The ring
 * is also one of the basal rings of the HC of which it is part. A mixed ring
 * must be a peripheral ring of the DDC of which it is part by definition.
 *
 * @var ring::strucType bothPrismatic
 * @brief The ring belongs to both a DDC and HC and is, thus, a 'mixed' ring.
 * The ring is also one of the prismatic rings of the HC of which it is part. A
 * mixed ring must be a peripheral ring of the DDC of which it is part by
 * definition (can never be an equatorial ring of a DDC and also be part of an
 * HC).
 *
 * @var ring::strucType Prism
 * @brief The ring belongs to a prism block, classified according to the prism
 * identification scheme.
 *
 * @var ring::strucType bothBasal
 * @brief The ring belongs to both a DDC and HC. It is a 'mixed' ring. The ring
 * is also one of the basal rings of the HC of which it is part. A mixed ring
 * must be a peripheral ring of the DDC of which it is part by definition.
 *
 * @var ring::strucType bothPrismatic
 * @brief The ring belongs to both a DDC and HC and is, thus, a 'mixed' ring.
 * The ring is also one of the prismatic rings of the HC of which it is part. A
 * mixed ring must be a peripheral ring of the DDC of which it is part by
 * definition (can never be an equatorial ring of a DDC and also be part of an
 * HC).
 *
 * @var ring::strucType Prism
 * @brief The ring belongs to a prism block, classified according to the prism
 * identification scheme.
 */
enum class strucType {
  unclassified,
  DDC,
  HCbasal,
  HCprismatic,
  bothBasal,
  bothPrismatic,
  Prism,
  deformedPrism,
  mixedPrismRing
};

//! Returns a vector of vectors of rings of a single size
std::vector<std::vector<int>>
getSingleRingSize(std::vector<std::vector<int>> rings, int ringSize);

//! Check to see if two vectors have common elements or not
//! True, if common elements are present and false if there are no common
//! elements
bool hasCommonElements(std::vector<int> ring1, std::vector<int> ring2);

//! Compares two disordered vectors and checks to see if they contain the same
//! elements
bool compareRings(std::vector<int> ring1, std::vector<int> ring2);

//! Searches a particular ring for a triplet
bool findTripletInRing(std::vector<int> ring, std::vector<int> triplet);

//! Common elements in 3 rings
bool commonElementsInThreeRings(std::vector<int> ring1, std::vector<int> ring2,
                                std::vector<int> ring3);

//! Returns the common elements of two rings
std::vector<int> findsCommonElements(std::vector<int> ring1,
                                     std::vector<int> ring2);

//! Erases memory for a vector of vectors for a list of rings
int clearRingList(std::vector<std::vector<int>> &rings);

//! Assign an atomType (equal to the number of nodes in the ring)
//! given n-membered rings.
int assignPolygonType(std::vector<std::vector<int>> rings,
                      std::vector<int> *atomTypes, std::vector<int> nRings);

} // namespace ring

#endif // __RINGS_H_
