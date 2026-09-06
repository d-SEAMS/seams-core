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

#include <algorithm>
#include <array>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <generic.hpp>
#include <mol_sys.hpp>

#ifdef SEAMS_HAS_LINKCELL
#include <linkcell.hpp>
#endif

/** @file neighbours.hpp
 *  @brief Header file for neighbour list generation.
 */

/**
 *  @addtogroup nneigh
 *  @{
 */

/** @brief Functions for building neighbour lists.
 * This namespace contains functions that build neighbour lists (vesin
 * cell list, with a brute-force fallback), saving either the atom IDs
 * or atom indices (according to a PointCloud) in a row-ordered vector
 * of vectors.
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

/// LAMMPS dump bound spans to restricted triclinic H (rows a, b, c).
///
/// `box[0..2]` are `xhi_bound - xlo_bound` etc. Optional `box[3..5]`
/// are tilt `xy, xz, yz`. `boxLow` is the bound lo. Inverse of
/// `xlo_bound = xlo + min(0, xy, xz, xy+xz)`.
inline void dumpBoundsToH(const std::vector<double> &box,
                          const std::vector<double> &boxLow, double H[3][3],
                          double origin[3]) {
  const double xspan = box.size() > 0 ? box[0] : 0.0;
  const double yspan = box.size() > 1 ? box[1] : 0.0;
  const double zspan = box.size() > 2 ? box[2] : 0.0;
  const double xlo_b = boxLow.size() > 0 ? boxLow[0] : 0.0;
  const double ylo_b = boxLow.size() > 1 ? boxLow[1] : 0.0;
  const double zlo_b = boxLow.size() > 2 ? boxLow[2] : 0.0;
  const double xy = box.size() >= 6 ? box[3] : 0.0;
  const double xz = box.size() >= 6 ? box[4] : 0.0;
  const double yz = box.size() >= 6 ? box[5] : 0.0;
  const double xmin = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
  const double xmax = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
  const double ymin = std::min(0.0, yz);
  const double ymax = std::max(0.0, yz);
  H[0][0] = xspan - xmax + xmin;
  H[0][1] = 0.0;
  H[0][2] = 0.0;
  H[1][0] = xy;
  H[1][1] = yspan - ymax + ymin;
  H[1][2] = 0.0;
  H[2][0] = xz;
  H[2][1] = yz;
  H[2][2] = zspan;
  origin[0] = xlo_b - xmin;
  origin[1] = ylo_b - ymin;
  origin[2] = zlo_b;
}

/// Recovered restricted-triclinic lengths lx, ly, lz (H diagonal).
inline void dumpCellLengths(const std::vector<double> &box,
                            const std::vector<double> &boxLow,
                            double lengths[3]) {
  double H[3][3];
  double origin[3];
  dumpBoundsToH(box, boxLow, H, origin);
  lengths[0] = H[0][0];
  lengths[1] = H[1][1];
  lengths[2] = H[2][2];
}

/// Longest recovered length: 0 = x, 1 = y, 2 = z.
inline int dumpAxialDim(const std::vector<double> &box,
                        const std::vector<double> &boxLow) {
  double lengths[3];
  dumpCellLengths(box, boxLow, lengths);
  int axial = 0;
  if (lengths[1] > lengths[axial]) {
    axial = 1;
  }
  if (lengths[2] > lengths[axial]) {
    axial = 2;
  }
  return axial;
}

inline int dumpAxialDim(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  return dumpAxialDim(yCloud.box, yCloud.boxLow);
}

/// Cartesian to fractional coordinates via dump H.
inline void dumpToFrac(const double H[3][3], const double origin[3], double x,
                       double y, double z, double s[3]) {
  const double lx = H[0][0];
  const double ly = H[1][1];
  const double lz = H[2][2];
  const double xy = H[1][0];
  const double xz = H[2][0];
  const double yz = H[2][1];
  s[2] = (z - origin[2]) / lz;
  s[1] = (y - origin[1] - yz * s[2]) / ly;
  s[0] = (x - origin[0] - xy * s[1] - xz * s[2]) / lx;
}

/// Fractional to cartesian coordinates via dump H.
inline void dumpFromFrac(const double H[3][3], const double origin[3],
                         const double s[3], double r[3]) {
  r[0] = origin[0] + H[0][0] * s[0] + H[1][0] * s[1] + H[2][0] * s[2];
  r[1] = origin[1] + H[1][1] * s[1] + H[2][1] * s[2];
  r[2] = origin[2] + H[2][2] * s[2];
}

/// Triclinic dump-cell volume |det(H)| from dumpBoundsToH.
/// Bound spans Lx*Ly*Lz are not the cell volume when tilt is present.
inline double dumpVolume(const std::vector<double> &box,
                         const std::vector<double> &boxLow) {
  double H[3][3];
  double origin[3];
  dumpBoundsToH(box, boxLow, H, origin);
  const double vol = H[0][0] * H[1][1] * H[2][2];
  return vol < 0.0 ? -vol : vol;
}

inline double dumpVolume(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  return dumpVolume(yCloud.box, yCloud.boxLow);
}

#ifdef SEAMS_HAS_LINKCELL
/// Map a LAMMPS dump box onto lc_cell.
inline lc_cell lammpsBoxToLcCell(const std::vector<double> &box,
                                 const std::vector<double> &boxLow) {
  double H[3][3], o[3];
  dumpBoundsToH(box, boxLow, H, o);
  lc_cell c = lc_cell_ortho(H[0][0], H[1][1], H[2][2]);
  c.ax = H[0][0];
  c.ay = H[0][1];
  c.az = H[0][2];
  c.bx = H[1][0];
  c.by = H[1][1];
  c.bz = H[1][2];
  c.cx = H[2][0];
  c.cy = H[2][1];
  c.cz = H[2][2];
  c.ox = o[0];
  c.oy = o[1];
  c.oz = o[2];
  return c;
}

inline linkcell::Cell lammpsBoxToLinkcell(const std::vector<double> &box,
                                          const std::vector<double> &boxLow) {
  return linkcell::Cell(lammpsBoxToLcCell(box, boxLow));
}

/// Cell and optional per-frame ortho lengths for analyzeResident.
/// boxLow set or nBox >= 6: one dump box via lammpsBoxToLinkcell,
/// frameLens is null. Otherwise Cell::ortho of box[0..2] and
/// frameLens is box.
inline void residentFrameCell(const double *box, const double *boxLow,
                              int nBox, linkcell::Cell &cell,
                              const double *&frameLens) {
  cell = linkcell::Cell::ortho(box[0], box[1], box[2]);
  frameLens = box;
  if (boxLow != nullptr || nBox >= 6) {
    std::vector<double> dump(static_cast<std::size_t>(std::max(nBox, 3)));
    for (int i = 0; i < nBox && i < static_cast<int>(dump.size()); ++i) {
      dump[static_cast<std::size_t>(i)] = box[i];
    }
    std::vector<double> lo;
    if (boxLow != nullptr) {
      lo = {boxLow[0], boxLow[1], boxLow[2]};
    }
    cell = lammpsBoxToLinkcell(dump, lo);
    frameLens = nullptr;
  }
}
#endif

//! All these functions use atom IDs and not indices

//! Full neighbour list for pairs of type I and type J (vesin cell list,
//! brute-force fallback). The neighbour list does not differentiate
//! between the types of atoms
std::vector<std::vector<int>> neighList(
    double rcutoff, const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    int typeI, int typeJ);

//! I-J neighbour list (I==J is like-type and reuses neighListO).
//! Unlike-type pairs use the same dump MIC as neighListO.
std::vector<std::vector<int>> neighListPair(
    double rcutoff, const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    int typeI, int typeJ);

//! Threaded cell-list neighbour rows. `rows[k]` holds the indices within
//! `rcutoff` of `subset[k]` under the minimum image convention, ascending
//! and without the atom itself. Atoms are binned in fractional coordinates
//! of the recovered triclinic cell H; every axis needs at least three cells
//! of perpendicular width `rcutoff`, so that the nearest image of every
//! neighbour sits in one of the 27 surrounding cells. Returns false, with
//! `rows` untouched, when the cell is too small for that, and the caller
//! takes the vesin or brute-force path. Each row is built by one thread.
bool cellListRowsThreaded(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &subset, double rcutoff,
    std::vector<std::vector<int>> &rows);

//! Full neighbour list for one type (threaded cell list from 2048 atoms
//! when OpenMP is compiled, vesin cell list, brute-force fallback)
//! You can only use this for neighbour lists with one type
std::vector<std::vector<int>> neighListO(
    double rcutoff, const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    int typeI);

//! Half neighbour list for one type (vesin cell list, brute-force fallback)
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
 *  2.5%. Nominations are a periodic linked-cell k-nearest search via
 *  linkcell (Allen and Tildesley): vesin is cutoff-only and KD-trees
 *  have no minimum-image convention. candidateCutoff is only a
 *  cell-size hint.
 *  On an undistorted tetrahedral lattice with k = 4 this graph equals
 *  the first-shell cutoff graph. Rows are by atom ID with the leading
 *  self entry, like neighListO.
 */
std::vector<std::vector<int>> kNearestNeighbourList(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, int typeI, bool mutual = true);

/** k-nearest graph on a type set (water oxygens). Substrate, ions and
 *  other species never appear as neighbours. An empty set is no atoms. */
std::vector<std::vector<int>> kNearestNeighbourList(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, const std::vector<int> &types, bool mutual = true);

/** Mutual and union k-nearest graphs from one candidate search. Seeded
 *  affiliation needs both; building them separately repeats the cell list.
 */
std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
kNearestNeighbourPair(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, int typeI);

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
kNearestNeighbourPair(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, const std::vector<int> &types);

/** The shell-separation certificate for the exact reduction of the k-nearest
 *  graph to the cutoff graph: returns {max_i d_k(i), min_i d_{k+1}(i)} over
 *  particles of the type. When max_i d_k(i) <= rcutoff <= min_i d_{k+1}(i),
 *  the two graphs coincide edge for edge and every downstream graph
 *  predicate is identical. Brute force, intended for validation.
 */
std::pair<double, double> shellSeparation(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    int typeI);

/** Nearest unlike image of each typeI particle among typeJ particles.
 *  Returns cloud-index pairs (i, j) and the MIC distance
 *  (sqrt of gen::periodicDistSq). On an ion cloud, typeI = 1
 *  (cation vertex) and typeJ = 2 (anion vertex). This is the
 *  nearest unlike neighbour, not a coordination number.
 */
std::vector<std::tuple<int, int, double>> nearestUnlike(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeI,
    int typeJ);

/** Subset of nearestUnlike where j's nearest typeI is i (mutual).
 *  A contact pair is this mutual edge, not first-shell membership.
 */
std::vector<std::pair<int, int>> mutualNearestUnlike(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeI,
    int typeJ);

//! Erases memory for a vector of vectors for the neighbour list
[[nodiscard]] int clearNeighbourList(std::vector<std::vector<int>> &nList);

/** Bond graph for TUM. Chosen at runtime so the three published
 *  assignments can be compared on the same frames.
 *  cutoff: pairs inside the distance cutoff (2020 graph).
 *  knn: mutual k-nearest (TUM v2 without hysteresis).
 *  knn-union: union k-nearest (the completion graph of the seeded rule).
 */
enum class BondGraph { Cutoff, KnnMutual, KnnUnion };

BondGraph bondGraphFromName(const std::string &name);
const char *bondGraphName(BondGraph graph);

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
  //! graph selects cutoff vs k-nearest (mutual or union). k is the
  //! neighbour count for the knn graphs (default 4).
  SkinNeighborList(double cutoff, double skin, int typeI,
                   BondGraph graph = BondGraph::KnnMutual, int k = 4);

  [[nodiscard]] BondGraph graph() const { return graph_; }

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
  int k_;
  BondGraph graph_{BondGraph::KnnMutual};
  bool mutual_{true};
  bool rebuilt_{true};
  int changedAtoms_{0};
  std::vector<double> x0_;
  std::vector<double> y0_;
  std::vector<double> z0_;
  std::vector<double> box0_;
  std::vector<std::pair<int, int>> candidates_;
  std::vector<std::pair<int, int>> bonded_;
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
