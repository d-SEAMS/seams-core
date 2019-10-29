#ifndef __MOL_SYS_H_
#define __MOL_SYS_H_

#include <sys/stat.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include "boost/multi_array.hpp"

/// For debugging, instantiate the unordered map
template class std::unordered_map<int, int>;

/*! \file mol_sys.hpp
    \brief The main molecular system handler.

    Details.
*/

/*!
 *  \addtogroup molSys
 *  @{
 */

/*! \brief Bare-bones structs used throughout the architecture.
 *         This namespace defines Point and PointCloud structs, alongwith other
 basic functions and enums.
 *
 PointCloud is a struct that contains the information of a collection of Point
 structs, at every frame. Point contains the coordinates and types of each
 point, along with information about whether the particular particle is within a
 slice or not.
  ### Changelog ###

  - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Sept 19, 2019
 */

namespace molSys {

// Enum type for bond type

/*! \enum molSys::bond_type
 * Qualifier for the bond type between two nearest-neighbours, according to the
 * CHILL or CHILL+ classification scheme.
 *
 * \var molSys::bond_type staggered
 * The bond is a staggered bond, according to the \f$a(i,j)\f$ or \f$c(i,j)\f$
 * value.
 *
 * \var molSys::bond_type eclipsed
 * The bond is an eclipsed bond.
 *
 * \var molSys::bond_type out_of_range
 * The bond cannot be classified as either staggered or eclipsed.
 */
enum bond_type { staggered, eclipsed, out_of_range };

/*! \enum molSys::atom_state_type
 * Qualifier for the per-particle phase state, according to the CHILL, CHILL+,
 * or \f$q_6\f$ order parameter.
 *
 * \var molSys::atom_state_type cubic
 * Ic, or particle type signifying Cubic Ice.
 *
 * \var molSys::atom_state_type hexagonal
 * Ih, or particle type signifying Hexagonal Ice.
 *
 * \var molSys::atom_state_type water
 * Liquid/amorphous phase.
 *
 * \var molSys::atom_state_type interfacial
 * Interfacial ice: ice-like molecules which do not fulfill the strict criteria
 * of the Ic or Ih phases.
 *
 * \var molSys::atom_state_type clathrate
 * Clathrate ice phase.
 *
 * \var molSys::atom_state_type interClathrate
 * Interfacial clathrate ice phase.
 *
 * \var molSys::atom_state_type unclassified
 * Not classified into any other category.
 *
 * \var molSys::atom_state_type reCubic
 * Reclassified as cubic ice, according to the \f$q_6\f$ order parameter.
 *
 * \var molSys::atom_state_type reHex
 * Reclassified as hexagonal ice, according to the \f$q_6\f$ order parameter.
 */
enum atom_state_type {
  cubic,
  hexagonal,
  water,
  interfacial,
  clathrate,
  interClathrate,
  unclassified,
  reCubic,
  reHex
};

/*! \struct Result
 * \brief This contains the bond classifier of enum type #bond_type, and the
 * bond correlation factor.
 *
 * Contains specifically the members:
 * - Bond classifier or the type of the bond (staggered, eclipsed, out-of-range)
 * - Bond correlation factor
 */
struct Result {
  bond_type classifier;  // Classifier according to CHILL, CHILL+ etc
  double c_value;        // Bond correlation factor
};

/*! \struct Point
 * \brief This contains per-particle information.
 *
 * Specifically
 * - Type ID
 * - Molecular ID
 * - Atom ID
 * - Coordinates
 * - Neighbourlist
 * - Bond correlation type of type Result
 * - Type of #atom_state_type iceType
 * - In slice bool
 */
// Struct that contains per-particle information
template <typename T>
struct Point {
  int type, molID, atomID;   // type ID, molID, atomID
  T x, y, z;                 // coordinates
  std::vector<Result> c_ij;  // Results (contains bond correlation type)
  atom_state_type iceType =
      molSys::unclassified;  // Type of ice/water etc based on cij
  bool inSlice = true;       // Is the point inside the slice or not?
};

// Struct for a collection of points; contains information for a particular
// frame
/*! \struct PointCloud
 * \brief This contains a collection of points; contains information for a
 * particular frame.
 *
 * Specifically
 * - A vector of Point structs
 * - The current frame number
 * - The number of particles in the current frame
 * - A vector for the simulation box lengths in each dimension
 * - A vector containing the absolute lower box coordinates
 */
template <typename S, typename T>
struct PointCloud {
  std::vector<S> pts;     // Collection of points
  int currentFrame;       // Current frame number
  int nop;                // Number of atoms
  std::vector<T> box;     // Periodic box lengths
  std::vector<T> boxLow;  // xlo, ylo, zlo
  std::unordered_map<int, int> idIndexMap;
};

// Creates an unordered map, with the atomIDs as keys and molecular IDs as the
// values
std::unordered_map<int, int> createIDMolIDmap(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//// Function for clearing vectors in PointCloud after multiple usage
molSys::PointCloud<molSys::Point<double>, double> clearPointCloud(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud);
}  // namespace molSys

#endif  // __MOL_SYS_H_
