//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#ifndef SEAMS_VORONOI_QLM_H_
#define SEAMS_VORONOI_QLM_H_

#include <bop.hpp>
#include <mol_sys.hpp>

#include <vector>

/** @file voronoi_qlm.hpp
 *  @brief Voronoi facet-area weighted bond order parameters.
 */

namespace chill {

/** @struct VoronoiWeights
 * @brief The Voronoi facet neighbours of one particle and their area weights.
 * @details Neighbours are those sharing a facet with the particle's Voronoi
 *  cell; weights are facet areas normalised to sum to one. Cutoff-based
 *  parameters change discontinuously when a neighbour crosses the cutoff,
 *  which is the shortcoming Mickel, Kapfer, Schr\"oder-Turk and Mecke
 *  identified (J. Chem. Phys. 138, 044501, 2013): facet areas shrink to zero
 *  continuously as a neighbour leaves the first shell.
 */
struct VoronoiWeights {
  std::vector<int> neighbours;  //! Cloud indices of facet-sharing neighbours
  std::vector<double> weights;  //! Facet areas, normalised to sum to one
  bool certified = false;       //! The exactness certificate held (see below)
};

/** Facet neighbours and area weights for every particle, with an a
 *  posteriori exactness certificate. candidateCutoff seeds the bisector
 *  search; the cutoff is enlarged automatically until the certificate holds
 *  or the growth cap is reached.
 *
 *  Certificate. The cell of particle i clipped against every candidate
 *  within cutoff c is the true Voronoi cell whenever every vertex of the
 *  clipped cell lies within c/2 of i: a particle k beyond the cutoff has
 *  |d_k| > c, its bisector plane sits at distance |d_k|/2 > c/2 from i, and
 *  a half-space whose boundary is farther than c/2 cannot cut a region
 *  contained in the ball of radius c/2. An open or under-clipped cell keeps
 *  vertices on the seeding square at distance ~c and fails the certificate,
 *  which is exactly the failure the enlargement retries.
 */
[[nodiscard]] std::vector<VoronoiWeights> voronoiFacetWeights(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double candidateCutoff);

//! Steinhardt parameters of degree orderL (3, 4, 6 or 8) with the bond sum
//! weighted by Voronoi facet areas. qlBar averages the weighted q_lm over the
//! particle and its facet neighbours, after Lechner and Dellago.
[[nodiscard]] SteinhardtQl steinhardtQlVoronoi(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double candidateCutoff, int orderL);

} // namespace chill

#endif // SEAMS_VORONOI_QLM_H_
