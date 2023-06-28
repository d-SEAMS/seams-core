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

#ifndef __CAGE_H_
#define __CAGE_H_
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
 * Type definitions for atoms, rings and cages which are DDCs and HCs.
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

#endif // __CAGE_H_
