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

#ifndef SEAMS_CAGE_H_
#define SEAMS_CAGE_H_
#include <map>
#include <string>
#include <string_view>
#include <vector>

/** @file cage.hpp
 *  @brief File for cage types for topological network criteria.
 *
 */

/**
 *  @addtogroup cage
 *  @{
 */

/** @brief Functions for topological network criteria cage types.
 *         This namespace contains structs and enums for cage types
 *
 * Type definitions for atoms, rings and cages. HexC and DoubleDiaC are
 * the TUM ice cages. A Signature names any polyhedron by its ring-size
 * census (sodalite is {4:6, 6:8}); the enumerator in cage_enum.hpp
 * grows face-sharing ring sets that match that census.
 *
 * ### Changelog ###
 *
 * - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 13, 2019
 * - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

// Namespace for cages
namespace cage {

// Type of a cage (a group of rings)
/** @enum class cage::cageType
 * Qualifier for a cage.
 * According to the topological network criterion for DDCs and HCs
 *
 * @var cage::cageType HexC
 * The type for a hexagonal cage
 *
 * @var molSys::cageType DoubleDiaC
 * The type for a double-diamond cage
 */
enum class cageType { HexC, DoubleDiaC };

/** Ring-size census of a polyhedral cage: size -> number of faces.
 *
 *  Parse a comma list (`4:6,6:8`) or a named table entry
 *  (`sodalite`, `alpha`, `512`, `51262`, `hc`, `ddc`). Counts of a
 *  repeated size add. Every size and count must be positive.
 *  The names `hc` and `ddc` keep the TUM ice finders (five hexagons
 *  and seven hexagons); a raw list `4:6,6:2` is the hexagonal prism.
 */
struct Signature {
  enum class Kind { Census, HexC, DoubleDiaC };

  std::map<int, int> counts;
  Kind kind = Kind::Census;

  static Signature parse(std::string_view spec);
  std::string str() const;
  int faceCount() const;
  int maxRingSize() const;
  bool containsSize(int size) const;
  bool operator==(const Signature &other) const {
    return kind == other.kind && counts == other.counts;
  }
  bool operator!=(const Signature &other) const { return !(*this == other); }
};

// Type of ice for a particular atom. Dummy means that the atom is unclassified
// and is most probably water
/** @enum class cage::iceType
 * Qualifier for an atom, based on whether it is part of cages,
 * according to the topological network criterion for DDCs and HCs
 *
 * @var cage::iceType dummy
 * Type for an atom which does not belong to any kind of cage
 *
 * @var molSys::iceType hc
 * Type for an atom which belongs to an HC
 *
 * @var molSys::iceType ddc
 * Type for an atom which belongs to a DDC
 *
 * @var molSys::iceType mixed
 * Type for an atom which is part of a mixed ring, shared by both a DDC and an
 * HC
 */
enum class iceType { dummy, hc, ddc, mixed, pnc, mixed2 };

// Each DDC has one equatorial ring and 6 peripheral rings
// Each HC has two basal planes and 3 prismatic planes
/** @struct Cage
 * @brief This contains a cage, with the constituent rings
 *
 * Contains specifically the members:
 * - Cage classifier or qualifier, for each cage (can be a DDC or HC)
 * - Vector of rings in the cage
 */
struct Cage {
  cageType type;          //! type of the cage : can be DDC or HC
  std::vector<int> rings; //! coordinates
};

} // namespace cage

#endif // SEAMS_CAGE_H_
