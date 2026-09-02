//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------
#ifndef SEAMS_CAGE_ENUM_H_
#define SEAMS_CAGE_ENUM_H_

#include <cage.hpp>

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

/** One closed polyhedron matching a signature. */
struct FoundCage {
  Signature signature;
  std::vector<int> faces;     ///< ring indices into the input vector
  std::vector<int> vertices;  ///< sorted unique atom indices
  std::string certificate;    ///< nauty hex, empty when nauty is off
};

/** True when every edge of the listed faces is used by exactly two
 *  of those faces. An empty face list is not closed. */
bool isClosedPolyhedron(const std::vector<std::vector<int>> &rings,
                        const std::vector<int> &faces);

/** Face-sharing rings whose size census equals `signature` and whose
 *  edges close. The input may hold rings of every size; only sizes
 *  that appear in the signature take part. */
std::vector<FoundCage>
findBySignature(const std::vector<std::vector<int>> &rings,
                const Signature &signature);

} // namespace cage

#endif // SEAMS_CAGE_ENUM_H_
