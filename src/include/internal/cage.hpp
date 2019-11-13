#ifndef __CAGE_H_
#define __CAGE_H_
#include <vector>

/*! \file cage.hpp
    \brief File for cage types for topological network criteria.

    Details.
*/

/*!
 *  \addtogroup cage
 *  @{
 */

/*! \brief Functions for topological network criteria cage types
 *         This namespace contains structs and enums for cage types
 *
Type definitions for atoms, rings and cages which are DDCs and HCs.

  ### Changelog ###

  - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 13, 2019
 */

// Namespace for cages
namespace cage {

// Type of a cage (a group of rings)
/*! \enum cage::cageType
 * Qualifier for a cage.
 * According to the topological network criterion for DDCs and HCs
 *
 * \var cage::cageType HexC
 * The type for a hexagonal cage
 *
 * \var molSys::cageType DoubleDiaC
 * The type for a double-diamond cage
 */
enum cageType { HexC, DoubleDiaC };

// Type of ice for a particular atom. Dummy means that the atom is unclassified
// and is most probably water
/*! \enum cage::iceType
 * Qualifier for an atom, based on whether it is part of cages,
 * according to the topological network criterion for DDCs and HCs
 *
 * \var cage::iceType dummy
 * Type for an atom which does not belong to any kind of cage
 *
 * \var molSys::iceType hc
 * Type for an atom which belongs to an HC
 *
 * \var molSys::iceType ddc
 * Type for an atom which belongs to a DDC
 *
 * \var molSys::iceType mixed
 * Type for an atom which is part of a mixed ring, shared by both a DDC and an HC
 */
enum iceType { dummy, hc, ddc, mixed };

// Each DDC has one equatorial ring and 6 peripheral rings
// Each HC has two basal planes and 3 prismatic planes
/*! \struct Cage
 * \brief This contains a cage, with the constituent rings
 *
 * Contains specifically the members:
 * - Cage classifier or qualifier, for each cage (can be a DDC or HC)
 * - Vector of rings in the cage
 */
struct Cage {
  cageType type;           // type of the cage : can be DDC or HC
  std::vector<int> rings;  // coordinates
};

}  // namespace cage

#endif  // __CAGE_H_
