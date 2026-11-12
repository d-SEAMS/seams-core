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
#include <cstdint>
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
 *
 *  findHC visits rings in the order of a key that does not use absolute
 *  coordinates: the sorted atom indices, then the unwrapped shift of the
 *  cycle. findDDC reports a maximum set of cages under the same key, so the
 *  cage count is independent of the order in which rings are claimed and
 *  does not change when the frame is translated.
 */
struct RingSearchIndex {
  std::vector<std::vector<int>> ringsContainingAtom;
  //! Same rows as ringsContainingAtom, packed. ringCol[ringOff[a] .. ringOff[a+1])
  //! is the ascending ring list of atom a. Empty when the index was not built.
  std::vector<int> ringOff;
  std::vector<int> ringCol;
  std::vector<int> atomRank;
  //! Integer minimum-image winding of each ring. Empty if no box was supplied.
  std::vector<std::array<int, 3>> ringShift;
};

//! Builds the inverted atom-to-rings index used by the cage searches
[[nodiscard]] RingSearchIndex
buildRingSearchIndex(const std::vector<std::vector<int>> &rings, int numAtoms);

//! As buildRingSearchIndex, ranking the atoms by position (ring::positionRank)
[[nodiscard]] RingSearchIndex buildRingSearchIndex(
    const std::vector<std::vector<int>> &rings,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud);

//! Rank of every atom of @a yCloud when the atoms are sorted by position
//! wrapped into the box, compared by x, then y, then z. Throws
//! std::invalid_argument when two atoms share a wrapped position
[[nodiscard]] std::vector<int>
positionRank(const molSys::PointCloud<molSys::Point<double>, double> &yCloud);

//! The rings of @a rings whose minimum-image bonds have integer lattice
//! images summing to zero. Any other ring winds through the periodic cell
//! and bounds no face. A cell narrower than a ring's span has such rings.
//! Without a box every ring is kept
[[nodiscard]] std::vector<std::vector<int>> contractibleRings(
    const std::vector<std::vector<int>> &rings,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud);

//! Faces for a cage search: contractibleRings, after refusing a cutoff that
//! can contain more than one periodic image of a pair. The narrowest face
//! separation must be longer than twice @a cutoff. A cloud without a box
//! keeps every ring. Throws std::invalid_argument when the cell is too thin
//! for @a cutoff
[[nodiscard]] std::vector<std::vector<int>> cageRings(
    const std::vector<std::vector<int>> &rings,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double cutoff);

//! Find out which hexagonal rings are DDC (Double Diamond Cages) rings.
//! Returns a vector containing all the ring IDs which are DDC rings
std::vector<int> findDDC(const std::vector<std::vector<int>> &rings,
                         std::vector<strucType> &ringType,
                         const std::vector<int> &listHC,
                         std::vector<cage::Cage> &cageList);

/** Proof that the reported number of equatorial rings is maximum.
 *  method 1: the mutual-conflict graph is the line graph of a bipartite
 *  graph, and a matching equals a König cover.
 *  method 2: each connected component was solved by branch-and-bound. The
 *  cover is a clique partition (an upper bound). When it equals the set,
 *  that partition is the certificate; otherwise the exhausted search is.
 */
struct PackingCertificate {
  int cages = 0;
  int upper = -1;
  int matching = 0;
  int cover = 0;
  int isolates = 0;
  int vertices = 0;
  int edges = 0;
  int components = 0;
  int dropped = 0;
  int method = 0;
  bool tight = false;
  double milliseconds = 0;
};

//! Peripheral rings of one equatorial candidate. n == 0 when the ring is
//! not a candidate. The first six sit inline. A defect that names more
//! conflicts stores the rest beside the row, so a crystal of candidates
//! stays a plain block of integers.
int equatorChoiceOverflow(const struct EquatorChoice *self, int at);
void equatorChoiceOverflowPut(struct EquatorChoice *self, int at, int value);

struct EquatorChoice {
  int ring[6]{};
  int n = 0;

  [[nodiscard]] int idAt(int k) const {
    return k < 6 ? ring[k] : equatorChoiceOverflow(this, k - 6);
  }

  void put(int k, int value) {
    if (k < 6) {
      ring[k] = value;
      return;
    }
    equatorChoiceOverflowPut(this, k - 6, value);
  }
};

void releaseEquatorOverflow();

//! Maximum set of mutually non-conflicting equatorial rings. @a periph[i].n
//! is 0 when ring i is not a candidate. The returned ids are a maximum
//! packing; @a certificate is filled when non-null.
[[nodiscard]] std::vector<int>
maximumEquators(const std::vector<EquatorChoice> &periph,
                const std::vector<int> &order, PackingCertificate *certificate);

//! As findDDC, reusing an index the caller has already built. When
//! @a certificate is given it receives the matching or search bound.
std::vector<int> findDDC(const std::vector<std::vector<int>> &rings,
                         std::vector<strucType> &ringType,
                         const std::vector<int> &listHC,
                         std::vector<cage::Cage> &cageList,
                         const RingSearchIndex &index,
                         PackingCertificate *certificate = nullptr);

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

//! Equatorial test of a DDC. Holds when every atom of ring iring lies in at
//! least three other rings and, for some starting atom and direction of
//! iring, six distinct rings P_k through its consecutive triplets satisfy
//! the even/odd common-atom and alternate-pair conditions. When
//! @a peripheralUnion is given it receives every ring that serves as a P_k
//! in some valid assignment, sorted. When @a chosen and @a rank are given,
//! @a chosen receives the six rings of the valid assignment whose ranks,
//! sorted, are lexicographically least, in ascending rank.
bool equatorialConditions(const std::vector<std::vector<int>> &rings,
                          int iring, const RingSearchIndex &index,
                          std::vector<int> *peripheralUnion,
                          std::vector<int> *chosen = nullptr,
                          const std::vector<int> *rank = nullptr,
                          EquatorChoice *raw = nullptr);

//! Membership of the neighbour rows after the leading self entry.
class NeighBond {
public:
  void build(const std::vector<std::vector<int>> &nList);

  // Row i, sorted and unique. Empty when i is outside the cloud.
  void row(int i, const int *&begin, const int *&end) const {
    if (static_cast<unsigned>(i) >= static_cast<unsigned>(n_) || off_.empty()) {
      begin = nullptr;
      end = nullptr;
      return;
    }
    begin = cols_.data() + off_[static_cast<std::size_t>(i)];
    end = cols_.data() + off_[static_cast<std::size_t>(i + 1)];
  }

  [[nodiscard]] bool operator()(int i, int j) const {
    const int *begin = nullptr;
    const int *end = nullptr;
    row(i, begin, end);
    if (begin == nullptr) {
      return false;
    }
    return std::binary_search(begin, end, j);
  }

private:
  int n_ = 0;
  // Row i is cols_[off_[i] .. off_[i + 1]], sorted and unique.
  std::vector<int> off_;
  std::vector<int> cols_;
};

//! Tests whether two rings are basal rings (true) or not (false): some
//! alternate triplet of basal1, starting at its first or second atom, is
//! bonded one-to-one to an alternate triplet of basal2, and no atom of the
//! other triplet of basal2 is bonded to basal1
bool basalConditions(const std::vector<std::vector<int>> &nList,
                     const std::vector<int> &basal1,
                     const std::vector<int> &basal2);

//! As basalConditions, answering bonds from a prebuilt NeighBond
bool basalConditions(const NeighBond &bonds, const std::vector<int> &basal1,
                     const std::vector<int> &basal2);

//! Tests whether the last two elements of a triplet are neighbours of two atom
//! IDs passed in
bool basalNeighbours(const std::vector<std::vector<int>> &nList,
                     std::vector<int> &triplet, int atomOne, int atomTwo);

//! Tests to check that elements of a triplet are not neighbours of a ring
//! (vector) passed
bool notNeighboursOfRing(const std::vector<std::vector<int>> &nList,
                         std::vector<int> &triplet,
                         const std::vector<int> &ring);

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

struct DdcSplit {
  double order_ms = 0;
  double equator_ms = 0;
  double pack_ms = 0;
  double conf_ms = 0;
  double split_ms = 0;
  double clique_ms = 0;
  double match_ms = 0;
  double search_ms = 0;
  double write_ms = 0;
};

struct PackParts {
  double conf_ms = 0;
  double split_ms = 0;
  double clique_ms = 0;
  double match_ms = 0;
  double search_ms = 0;
};

//! Conflict build, clique ids, matching and branch-and-bound from the last pack.
[[nodiscard]] PackParts lastPackParts();

//! Order, equatorial test and packing from the last findDDC call.
[[nodiscard]] DdcSplit lastDdcSplit();

} // namespace ring
/** @} */

/** @defgroup prism3 Bulk prisms
 *  Candidate basal rings and prism blocks in bulk ice.
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
/** @} */

#endif // SEAMS_TOPO_BULK_H_
