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

#ifndef SEAMS_NEIGHBOURS_H_
#define SEAMS_NEIGHBOURS_H_

#include <array>
#include <set>
#include <utility>

#include <generic.hpp>
#include <mol_sys.hpp>

/** @file neighbours.hpp
 *  @brief Header file for neighbour list generation.
 */

/**
 *  @addtogroup nneigh
 *  @{
 */

/** @brief Functions for building neighbour lists.
 * This namespace contains functions that build neighbour lists (using
 * brute-force), saving either the atom IDs or atom indices (according to a
 * PointCloud) in a row-ordered vector of vectors.
 * Whether the atom IDs or atom indices (i.e. the indices of the elements in the
 * vector pts inside the PointCloud) are saved, the neighbour lists are
 * constructed such that the first element is the 'central atom', whose
 * neighbours are being saved on that particular line. The central atom is
 * followed by the atom IDs or indices of the nearest neighbours.
 *
 * In a 'full' neighbour list, if 1 is a neighbour of 2 then 1 is saved in the
 * neighbour list of 2 AND 2 is also saved in the neighbour list of 1.
 *
 * In a 'half' neighbour list, if 1 is a neighbour of 2, then 2 is saved in the
 * neighbour list of 1 but not vice versa.
 *
 *   ### Changelog ###
 *
 * - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 14, 2019
 * - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

namespace nneigh {
//! All these functions use atom IDs and not indices

//! Inefficient @f$O(n^2)@f$ implementation of neighbour lists when there are
//! two different types of atoms The neighbour list does not differentiate
//! between the types of atoms
std::vector<std::vector<int>> neighList(
    double rcutoff, const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    int typeI, int typeJ);

//! Inefficient @f$O(n^2)@f$ implementation of neighbour lists
//! You can only use this for neighbour lists with one type
std::vector<std::vector<int>> neighListO(
    double rcutoff, const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    int typeI);

//! Inefficient @f$O(n^2)@f$ implementation of neighbour lists
//! You can only use this for neighbour lists with one type
std::vector<std::vector<int>> halfNeighList(
    double rcutoff, const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    int typeI = 1);

//! The following function outputs a neighbour list using indices and NOT atom
//! IDs

//! Converts the neighbour list build with atom IDs into a neighbour list of
//! atom indices, according to the pointCloud
std::vector<std::vector<int>> neighbourListByIndex(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList);

//! Gets a neighbour list by index, according to a pointCloud given as the
//! input. Assume no slices or other skullduggery
std::vector<std::vector<int>> getNewNeighbourListByIndex(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff);

/** Bonded graph from the k nearest neighbours of each particle rather than a
 *  distance cutoff. The default intersection symmetrization (mutual = true)
 *  bonds i and j only when each lists the other among its k nearest; the
 *  union alternative bonds on either nomination. In a crystal the first
 *  shell is mutual and the two coincide; on disordered packings the mutual
 *  graph is sparser, which starves accidental ring structure -- measured on
 *  the dense null, mutual scores zero false crystal where union reaches
 *  2.5%. Candidates come from the cell list at
 *  candidateCutoff, which must comfortably exceed the k-th neighbour
 *  distance. On an undistorted tetrahedral lattice with k = 4 this graph
 *  equals the first-shell cutoff graph; under thermal distortion it keeps
 *  the neighbour identities a hard cutoff loses, which is what the cage
 *  predicates need. Rows are by atom ID with the leading self entry, like
 *  neighListO.
 */
std::vector<std::vector<int>> kNearestNeighbourList(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, int typeI, bool mutual = true);

/** The shell-separation certificate for the exact reduction of the k-nearest
 *  graph to the cutoff graph: returns {max_i d_k(i), min_i d_{k+1}(i)} over
 *  particles of the type. When max_i d_k(i) <= rcutoff <= min_i d_{k+1}(i),
 *  the two graphs coincide edge for edge and every downstream graph
 *  predicate is identical. Brute force, intended for validation.
 */
std::pair<double, double> shellSeparation(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    int typeI);

//! Erases memory for a vector of vectors for the neighbour list
[[nodiscard]] int clearNeighbourList(std::vector<std::vector<int>> &nList);

/** Persistent neighbour list with a LAMMPS skin.
 *  Vesin (or the brute-force fallback) builds candidates at cutoff+skin,
 *  which is the ghost halo: periodic images already sit in that shell.
 *  The cell list is rebuilt only when some atom has moved more than
 *  skin/2 from the last rebuild, the Verlet trigger. The bond graph
 *  for rings is the candidate pairs whose current distance is inside
 *  the analysis cutoff.
 */
class SkinNeighborList {
public:
  SkinNeighborList(double cutoff, double skin, int typeI);

  //! Refresh from a new frame. The returned list is ID-keyed with a
  //! leading self entry, the same shape as neighListO.
  const std::vector<std::vector<int>> &
  update(const molSys::PointCloud<molSys::Point<double>, double> &yCloud);

  [[nodiscard]] bool lastRebuilt() const { return rebuilt_; }
  //! Atoms whose cutoff bond set changed on the last update.
  [[nodiscard]] int lastChangedAtoms() const { return changedAtoms_; }
  [[nodiscard]] const std::vector<std::vector<int>> &bonds() const {
    return nList_;
  }

private:
  double cutoff_;
  double skin_;
  double cutoffSq_;
  double triggerSq_;
  int typeI_;
  bool rebuilt_{true};
  int changedAtoms_{0};
  std::vector<double> x0_;
  std::vector<double> y0_;
  std::vector<double> z0_;
  std::array<double, 3> box0_{};
  std::vector<std::pair<int, int>> candidates_;
  std::set<std::pair<int, int>> bonded_;
  std::vector<std::vector<int>> nList_;

  bool mustRebuild(
      const molSys::PointCloud<molSys::Point<double>, double> &yCloud) const;
  void rebuildCandidates(
      const molSys::PointCloud<molSys::Point<double>, double> &yCloud);
  void refreshBonds(
      const molSys::PointCloud<molSys::Point<double>, double> &yCloud);
};

}  // namespace nneigh

#endif  // SEAMS_NEIGHBOURS_H_
