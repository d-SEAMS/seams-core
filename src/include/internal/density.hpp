//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#ifndef SEAMS_DENSITY_H_
#define SEAMS_DENSITY_H_

#include <mol_sys.hpp>
#include <site.hpp>

#include <vector>

/** @file density.hpp
 *  @brief Type- or kind-resolved Cartesian number-density histogram.
 *
 *  This is not a 2-periodic MIC. Each slab volume is A_perp * dz,
 *  with A_perp from dump H (ly*lz for x, lx*lz for y, lx*ly for z)
 *  using recovered lx, ly, lz, not bound spans.
 */

namespace site {

struct DensityZ {
  std::vector<double> z;   //! Bin centres along the chosen axis
  std::vector<double> rho; //! Number density in each slab
  int type{0};             //! LAMMPS type; 0 = all atoms
};

//! Histogram Point::{x,y,z} for typeI (0 = every atom). axis is 0, 1, or 2.
DensityZ densityZ(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeI,
    int nbin, int axis);

//! Same histogram, restricted to atoms whose Table kind matches.
DensityZ densityZ(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const Table &table, Kind kind, int nbin, int axis);

} // namespace site

#endif // SEAMS_DENSITY_H_
