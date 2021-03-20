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

/*! \file ring.hpp
    \brief File containing common functions used by bulk and confined
   topological network critera.

    Details.
*/

/*!
 *  \addtogroup ring
 *  @{
 */

/*! \brief Topological network criteria functions.
 *         This namespace contains functions for the topological network
 criteria for bulk and confined systems both.
 *
 Although the namespace is shared by bulk and confined topological network
 criteria, the functions are split into files specific to each set of criteria.

 All the topological network criteria are based on the identification of
 primitive rings (using the Franzblau algorithm for primitive rings), following
 which specific rules for combinations of primitive rings are used to identify
 quasi-one-dimensional, quasi-two-dimensional and bulk ices.

  ### Changelog ###

  - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 13, 2019
 */

namespace ring {

// General enum used throughout this program (Prism is for our prism
// classification) {per-ring classification}
/*! \enum ring::strucType
 * Qualifier for each ring, based on the classification type determined by the
 * bulk or confined ice topological network criteria.
 *
 * \var ring::strucType unclassified
 * The ring is unclassified, which may be either water or a deformed type which
 * cannot be classified by the criteria.
 *
 * \var ring::strucType DDC
 * The ring belongs to a double-diamond cage (DDC).
 *
 * \var ring::strucType HCbasal
 * The ring belongs only to a hexagonal cage (HC). Specifically, the ring is
 * purely a basal ring of an HC.
 *
 * \var ring::strucType HCprismatic
 * The ring belongs only to a hexagonal cage (HC); specifically the ring is
 * purely a prismatic ring of an HC. It is not shared by a DDC.
 *
 * \var ring::strucType bothBasal
 * The ring belongs to both a DDC and HC. It is a 'mixed' ring. The ring is also
 * one of the basal rings of the HC of which it is part. A mixed ring must be a
 * peripheral ring of the DDC of which it is part by definition.
 *
 * \var ring::strucType bothPrismatic
 * The ring belongs to both a DDC and HC and is, thus, a 'mixed' ring. The ring
 * is also one of the prismatic rings of the HC of which it is part. A mixed
 * ring must be a peripheral ring of the DDC of which it is part by definition
 * (can never be an equatorial ring of a DDC and also be part of an HC).
 *
 * \var ring::strucType Prism
 * The ring belongs to a prism block, classified according to the prism
 * identification scheme.
 */
enum strucType {
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

// Returns a vector of vectors of rings of a single size
std::vector<std::vector<int>>
getSingleRingSize(std::vector<std::vector<int>> rings, int ringSize);

// Check to see if two vectors have common elements or not
// True, if common elements are present and false if there are no common
// elements
bool hasCommonElements(std::vector<int> ring1, std::vector<int> ring2);

// Compares two disordered vectors and checks to see if they contain the same
// elements
bool compareRings(std::vector<int> ring1, std::vector<int> ring2);

// Searches a particular ring for a triplet
bool findTripletInRing(std::vector<int> ring, std::vector<int> triplet);

// Common elements in 3 rings
bool commonElementsInThreeRings(std::vector<int> ring1, std::vector<int> ring2,
                                std::vector<int> ring3);

// Returns the common elements of two rings
std::vector<int> findsCommonElements(std::vector<int> ring1,
                                     std::vector<int> ring2);

// Erases memory for a vector of vectors for a list of rings
int clearRingList(std::vector<std::vector<int>> &rings);

} // namespace ring

#endif // __RINGS_H_
