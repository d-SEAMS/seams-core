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

#ifndef __MOL_SYS_H_
#define __MOL_SYS_H_

#include "boost/multi_array.hpp"
#include <algorithm>
#include <array>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>
#include<unordered_map>

// For debugging, instantiate the unordered map [consider removal for
// production]
template class std::unordered_map<int, int>;
template class std::vector<std::vector<int>>;

/** @file mol_sys.hpp
 *  @brief The main molecular system handler.
 */

/**
 *  @addtogroup molSys
 *  @{
 */

/** @brief Bare-bones structs used throughout the architecture.
 *  @detials This namespace defines Point and PointCloud structs, alongwith
 * other basic functions and enums.
 *
 * PointCloud is a struct that contains the information of a collection of Point
 * structs, at every frame. Point contains the coordinates and types of each
 * point, along with information about whether the particular particle is within
 * a slice or not.
 *  ### Changelog ###
 *
 * - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Sept 19, 2019
 * - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

namespace molSys {

// Enum type for bond type

/** @enum class molSys::bond_type
 * @brief Qualifier for the bond type between two nearest-neighbours, according
 *  to the CHILL or CHILL+ classification scheme.
 *
 * @var molSys::bond_type staggered
 * @brief The bond is a staggered bond, according to the @f$a(i,j)@f$ or
 *  @f$c(i,j)@f$ value.
 *
 * @var molSys::bond_type eclipsed
 * @brief The bond is an eclipsed bond.
 *
 * @var molSys::bond_type out_of_range
 * @brief The bond cannot be classified as either staggered or eclipsed.
 */
enum class bond_type { staggered, eclipsed, out_of_range };

/** \enum class molSys::atom_state_type
 * @brief Qualifier for the per-particle phase state, according to the CHILL,
 *  CHILL+, or @f$q_6@f$ order parameter.
 *
 * @var molSys::atom_state_type cubic
 * @brief Ic, or particle type signifying Cubic Ice.
 *
 * @var molSys::atom_state_type hexagonal
 * @brief Ih, or particle type signifying
 *  Hexagonal Ice.
 *
 * @var molSys::atom_state_type water
 * @brief Liquid/amorphous phase.
 *
 * @var molSys::atom_state_type interfacial
 * @brief Interfacial ice: ice-like molecules
 *  which do not fulfill the strict criteria of the Ic or Ih phases.
 *
 * @var molSys::atom_state_type clathrate
 * @brief Clathrate ice phase.
 *
 * @var molSys::atom_state_type interClathrate
 * @brief Interfacial clathrate ice phase.
 *
 * @var molSys::atom_state_type unclassified
 * @brief Not classified into any other
 *  category.
 *
 * @var molSys::atom_state_type reCubic
 * @brief Reclassified as cubic ice, according to
 *  the @f$q_6@f$ order parameter.
 *
 * @var molSys::atom_state_type reHex
 * @brief Reclassified as hexagonal ice, according
 *  to the @f$q_6@f$ order parameter.
 */
enum class atom_state_type {
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

/** @struct Result
 * @brief This contains the bond classifier of enum class type #bond_type, and the
 *  bond correlation factor.
 * @details Contains specifically the members:
 * - Bond classifier or the type of the bond (staggered, eclipsed, out-of-range)
 * - Bond correlation factor
 */
struct Result {
  bond_type classifier; //! Classifier according to CHILL, CHILL+ etc
  double c_value;       //! Bond correlation factor
};

/** @struct Point
 * @brief This contains per-particle information.
 * @details Specifically
 *  - Type ID
 *  - Molecular ID
 *  - Atom ID
 *  - Coordinates
 *  - Neighbourlist
 *  - Bond correlation type of type Result
 *  - Type of #atom_state_type iceType
 *  - In slice bool
 */
template <typename T> struct Point {
  int type, molID, atomID;  //! type ID, molID, atomID
  T x, y, z;                //! coordinates
  std::vector<Result> c_ij; //! Results (contains bond correlation type)
  atom_state_type iceType =
      molSys::atom_state_type::unclassified; //! Type of ice/water etc based on cij
  bool inSlice = true;      //! Is the point inside the slice or not?
};

// Struct for a collection of points; contains information for a particular
// frame
/** @struct PointCloud
 *  @brief This contains a collection of points; contains information for a
 *  particular frame.
 * @details Specifically
 *  - A vector of Point structs
 *  - The current frame number
 *  - The number of particles in the current frame
 *  - A vector for the simulation box lengths in each dimension
 *  - A vector containing the absolute lower box coordinates
 */
template <typename S, typename T> struct PointCloud {
  std::vector<S> pts;    //! Collection of points
  int currentFrame;      //! Current frame number
  int nop;               //! Number of atoms
  std::vector<T> box;    //! Periodic box lengths
  std::vector<T> boxLow; //! xlo, ylo, zlo
  std::unordered_map<int, int> idIndexMap;
};

//! Creates an unordered map, with the atomIDs as keys and molecular IDs as the
//! values
std::unordered_map<int, int>
createIDMolIDmap(molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Creates an multimap with molecule IDs of the atoms as the keys and the
//! atom IDs as the values. More than one atom can have the same molecule ID
std::unordered_multimap<int, int>
createMolIDAtomIDMultiMap(molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Returns a vector of vectors, which contains the molIDs in the first column,
//! and the hydrogen atom indices (not atom IDs) in the row
std::vector<std::vector<int>>
hAtomMolList(molSys::PointCloud<molSys::Point<double>, double> *hCloud,
             molSys::PointCloud<molSys::Point<double>, double> *oCloud);

//! This function searches a vector of vectors molList, for a particular
//! molecular ID, and returns the index in molList
int searchMolList(std::vector<std::vector<int>> molList, int molIDtoFind);

//!//! Function for clearing vectors in PointCloud after multiple usage
molSys::PointCloud<molSys::Point<double>, double>
clearPointCloud(molSys::PointCloud<molSys::Point<double>, double> *yCloud);
} // namespace molSys

#endif // __MOL_SYS_H_
