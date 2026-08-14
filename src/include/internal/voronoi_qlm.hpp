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
};

//! Facet neighbours and area weights for every particle. candidateCutoff
//! bounds the bisector-plane search and must exceed the largest facet
//! neighbour distance; twice the first-shell distance is ample.
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
