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

#ifndef SEAMS_RDF_H_
#define SEAMS_RDF_H_

#include <vector>

#include <mol_sys.hpp>

/** @file rdf.hpp
 *  @brief Partial 3D radial distribution functions under the dump MIC.
 */

/**
 *  @addtogroup rdf
 *  @{
 */

/** @brief Three-dimensional partial RDF g_IJ(r).
 *
 *  Pair distances use gen::periodicDist, so a sheared dump cell uses the
 *  same minimum image as the neighbour lists. The ideal-gas denominator
 *  uses the triclinic volume from nneigh::dumpVolume (det H), not the
 *  product of bound spans.
 */
namespace rdf {

struct PartialRdf {
  std::vector<double> r;    //! Bin centres
  std::vector<double> g;    //! g_IJ at each bin
  std::vector<int> count;   //! Unordered pair counts
  double rmax{0.0};
  double binwidth{0.0};
  double volume{0.0};
  int typeI{0};
  int typeJ{0};
  int nI{0};
  int nJ{0};
};

//! Histogram I-J pair distances and normalize to the dump-cell ideal gas.
PartialRdf partialRdf(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeI,
    int typeJ, double rmax, int nbins);

//! Running site-site CN: 4 pi rho_J integral_0^{r} s^2 g_IJ(s) ds per bin.
std::vector<double> runningCN(const PartialRdf &h, double rhoJ);

//! First minimum of g after the first maximum. Returns the bin index, or -1.
int firstMinimumBin(const PartialRdf &h);

//! Site-site CN integrated to rMax. rhoJ is the number density of type J.
double coordinationNumber(const PartialRdf &h, double rMax, double rhoJ);

} // namespace rdf

#endif // SEAMS_RDF_H_
