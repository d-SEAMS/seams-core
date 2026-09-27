//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------
#ifndef SEAMS_CAGE_ENUM_H_
#define SEAMS_CAGE_ENUM_H_

#include <cage.hpp>

#include <array>
#include <map>
#include <string>
#include <vector>

/** @file cage_enum.hpp
 *  @brief Signature-guided enumeration of polyhedral cages.
 *
 *  A cage is a connected set of primitive rings (faces) whose size
 *  census matches a Signature and whose edges form a closed
 *  polyhedron: every edge of every face is shared by exactly two
 *  faces of the set. Growth walks the ring adjacency graph (two
 *  rings are adjacent when they share an edge) and stays inside the
 *  signature budget. Distinct cages are the distinct sorted vertex
 *  sets; the nauty certificate, when linked, names the isomorphism
 *  class of each cage.
 */
namespace cage {

/** One cage matching a signature, closed or incomplete. */
struct FoundCage {
  Signature signature;
  std::vector<int> faces;     ///< ring indices into the input vector
  std::vector<int> vertices;  ///< sorted unique atom indices
  std::string certificate;    ///< nauty hex, empty when nauty is off
  bool closed = true;         ///< every face edge used twice
  int danglingEdges = 0;      ///< edges used once (cups / incomplete)
};

/** True when every edge of the listed faces is used by exactly two
 *  of those faces. An empty face list is not closed. */
bool isClosedPolyhedron(const std::vector<std::vector<int>> &rings,
                        const std::vector<int> &faces);

/** Face-sharing rings whose size census equals `signature` and whose
 *  edges close. The input may hold rings of every size; only sizes
 *  that appear in the signature take part. Named `hc` and `ddc`
 *  without a neighbour list use the geometric census
 *  (`4:6,6:2` and `6:7`). */
std::vector<FoundCage>
findBySignature(const std::vector<std::vector<int>> &rings,
                const Signature &signature);

/** True when `cycle` is `pattern` rotated or reversed. */
bool speciesCycleMatches(const std::vector<int> &cycle,
                         const std::vector<int> &pattern);

/** As findBySignature. A ring is a face only when its species sequence
 *  matches one pattern of the same length, up to rotation and reversal.
 *  An empty pattern list keeps every ring whose size is in the census.
 *  `species` is one class per atom index. A vertex past the end of
 *  `species` keeps its ring out. */
std::vector<FoundCage>
findBySignature(const std::vector<std::vector<int>> &rings,
                const Signature &signature, const std::vector<int> &species,
                const std::vector<std::vector<int>> &patterns);

/** As above. Named `hc` and `ddc` call findHC / findDDC on the
 *  six-membered rings so the vertex sets match those finders. */
std::vector<FoundCage>
findBySignature(const std::vector<std::vector<int>> &rings,
                const std::vector<std::vector<int>> &nList,
                const Signature &signature);

/** Connected face sets that stay inside the signature budget, have at
 *  least `minFaces` faces, and are not closed. Cups and incomplete
 *  cages during hydrate nucleation. Closed polyhedra are omitted.
 *  `minFaces <= 0` uses `max(1, signature.faceCount() / 2)`. A cup
 *  whose every face already belongs to a closed cage is omitted. */
std::vector<FoundCage>
findIncompleteBySignature(const std::vector<std::vector<int>> &rings,
                          const Signature &signature, int minFaces);

/** IRA/SOFI result for one cage. status 0 means the library ran on
 *  these vertices. status 1 means it is absent or the cloud is empty.
 *  nVertices is the cage, not the frame. */
struct CageShape {
  int nVertices = 0;
  int status = 1;
  double rmsd = -1.0;
  std::string pointGroup;
};

/** Coordinates of `vertices` only. An index past `all` is skipped. */
std::vector<std::array<double, 3>>
coordsOfVertices(const std::vector<std::array<double, 3>> &all,
                 const std::vector<int> &vertices);

/** Point group of this vertex set. Does not see any other atom. */
CageShape shapeOfVertices(const std::vector<std::array<double, 3>> &cageXyz);

/** Overlay `cageXyz` on `ref`. Row counts must agree. Does not see
 *  any atom outside the two sets. */
CageShape overlayVertices(const std::vector<std::array<double, 3>> &ref,
                          const std::vector<std::array<double, 3>> &cageXyz);

/** One network-former atom: coordination, same-species bonds, and the
 *  primitive rings that pass through it. `rings` counts are through
 *  this atom, so a ring of size n contributes to n rows. */
struct FormerRow {
  int index = -1;
  int species = 0;
  int coord = 0;
  int homopolar = 0;
  std::map<int, int> rings;
};

/** `formerSpecies < 0` keeps every atom. */
std::vector<FormerRow>
formerRows(const std::vector<std::vector<int>> &nList,
           const std::vector<std::vector<int>> &rings,
           const std::vector<int> &species, int formerSpecies);

/** True when the two tables have the same per-atom species, coordination,
 *  homopolar count, and ring-size multiset. */
bool sameNetwork(const std::vector<FormerRow> &early,
                 const std::vector<FormerRow> &late);

} // namespace cage

#endif // SEAMS_CAGE_ENUM_H_
