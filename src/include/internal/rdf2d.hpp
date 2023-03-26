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

#ifndef __RDF2D_H_
#define __RDF2D_H_

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <math.h>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include <mol_sys.hpp>
#include <order_parameter.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

/** @file rdf2d.hpp
 *  @brief File containing functions used to calculate the in-plane radial
 * distribution functions.
 */

/**
 *  @addtogroup rdf2
 *  @{
 */

/** @brief Defines RDF-specific functions.
 *
 * @details The radial distribution function @f$g(r)@f$, or pair distribution
 *  function. It can be used to illuminate features of short-range and
 * long-range order. The RDF is the probability of finding a particle at a
 * distance of
 *  @f$r@f$ from a tagged reference particle, relative to that of an ideal gas.
 *  For a system of @f$N@f$ particles, the pair correlation function for
 *  @f$N(N-1)@f$ pairs is:
 *
 *   @f[
 *   \rho_N^{(2)}(r,r') = ⟨\sum_{i=1}^{N} \sum_{j=1,j \neq i}^{N} \delta
 *  (r-r_i) \delta (r'-r_j) \langle @f]
 *
 *  The code essentially bins distances between pairs of particles, and
 *   normalizes the resulting histogram is normalized with respect to an ideal
 *   gas. The algorithm for the calculation of @f$g(r)@f$ may be divided into
 * the following steps:
 *
 *  1. @b Initialization: The @f$g(r)@f$ array is initialized to zero.
 *  2. @b Sampling: The histogram is added to for a particular bin, if the
 *  distance of a pair of atoms falls within the @f$r@f$ associated with the
 * bin.
 *  3. @b Normalization: Every bin of the @f$g(r)@f$ array is normalized by
 * the product of the number of ideal gas particles in that bin, and the number
 * of particles and number of frames.
 *
 *  To account for the geometry of the slab-like quasi-two-dimensional system,
 * the RDF must be additionally normalized by a form factor @f$f@f$ such that
 * @f$f = \frac{h}{2r}@f$ for @f$r>h@f$ and @f$f = \frac{h}{2r}@f$ for @f$r \le
 * h@f$; where @f$h@f$ is the height of the quasi-two-dimensional system.
 *
 *    ### Changelog ###
 *
 *   - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Dec 10, 2019
 *   - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

namespace rdf2 {

//! Main function for calculating the RDF for the same type of particle: calls
//! other functions for initializing, sampling and normalizing the RDF
int rdf2Danalysis_AA(std::string path, std::vector<double> *rdfValues,
                     molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                     double cutoff, double binwidth, int firstFrame,
                     int finalFrame);

//! Samples the RDF histogram at every step
std::vector<int>
sampleRDF_AA(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
             double cutoff, double binwidth, int nbin);

//! Normalize the histogram
int normalizeRDF(int nopA, std::vector<double> *rdfValues,
                 std::vector<int> histogram, double binwidth, int nbin,
                 std::vector<double> volumeLengths, int nIter);

//! Gets the lengths of the volume slice of the quasi-two-dimensional system
std::vector<double>
getSystemLengths(molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Gets the plane area from the volume lengths vector
double getPlaneArea(std::vector<double> volumeLengths);

} // namespace rdf2

#endif // __RDF2D_H_
