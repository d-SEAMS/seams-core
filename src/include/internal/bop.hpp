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

#ifndef __BOP_H_
#define __BOP_H_

#include <array>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cmath>
#include <complex>
#include <generic.hpp>
#include <math.h>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_output.hpp>

/** @file bop.hpp
 *  @brief File for the bond order parameter analysis.
 */

/**
 *  @addtogroup chill
 *  @{
 */

/** @brief CHILL and CHILL+ structure classification.
 *         This namespace contains functions that are used in the
CHILL/CHILL+ classification scheme, as well as a yodaCloud struct to hold
 all the information.
 *
 Both the <a
href="https://pubs.rsc.org/en/content/articlehtml/2010/cp/b919724a">CHILL</a>
and <a href="https://pubs.acs.org/doi/abs/10.1021/jp510289t">CHILL+</a> methods
are based on the <a href="https://aip.scitation.org/doi/10.1063/1.471721">local
bond order parameter method</a>, developed by ten Wolde et al. for the
identification of crystal nuclei in Lennard-Jones systems. The local order
 parameter method is based on <a
href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.47.1297">order
parameters</a> introduced by Steinhardt et al.

 The local environment of the water molecules is classified using an algorithm
based on spherical harmonics, which is independent of the specific crystal
structure and does not require the definition of a reference frame. The local
order around each water molecule @f$i@f$ is described by a local orientational
bond order parameter @f$q_{lm}(i)@f$, with @f$ 2l+1 @f$ complex components:

  @f[
  q_{lm}(i) = @frac{1}{N_b(i)} \Sigma_{j=1}^{N_b(i)} Y_{lm}(r_{ij})
  @f]

  Here, @f$N_b(i)=4@f$ is the number of nearest neighbours for water molecule
@f$i@f$, @f$l@f$ is a free integer parameter, @f$m@f$ is an integer such that
@f$m \in [-l,l]@f$. The functions @f$Y_{lm}(r_{ij})@f$ are the spherical
harmonics and @f$r_{ij}@f$ is the distance vector between molecule @f$i@f$ and
each one of its four nearest neighbours @f$j@f$.

  The alignment of the orientation of the local structures is measured by the
normalized dot product of the local orientational bond order parameter,
\f$q_l(i)\f$, between each pair of neighbor molecules \f$i\f$ and \f$j\f$, given
by:

  @f[
  a(i,j) = \frac{q_l(i).q_l(j)}{|q_l(i)||q_l(j)|} = \frac{ \Sigma_{m=-l}^{l}
q_l(i) q_l^*(j)}{ (\Sigma_{m=-l}^{l} q_l(j) q_l^*(i))^{1/2} (\Sigma_{m=-l}^{l}
q_l(j) q_l^*(i))^{1/2} } @f]

  where, @f$q_{lm}^*@f$ is the complex conjugate of @f$q_{lm}@f$.

  Depending on the values of @f$a(i,j)@f$, the bond between each pair of
molecules @f$i@f$ and @f$j@f$ is classified as being either staggered or
eclipsed.

  | Algorithm | Eclipsed Bonds                | Staggered Bonds       |
  | :----:    | :----:                        | :----:                |
  | CHILL     | @f$-0.05 >= a(i,j) >= -0.2@f$ | @f$a(i,j) <= -0.8@f$  |
  | CHILL+    | @f$0.25 >= a(i,j) >= -0.35@f$ | @f$a(i,j) <= -0.8@f$  |

  Since each molecule @f$i@f$ in deeply supercooled water has four nearest
neighbours, the type of the resultant four bonds is used to identify phases
according to the CHILL and CHILL+ algorithms.

  The CHILL algorithm can classify water molecules as belonging to the cubic,
hexagonal, interfacial or liquid (amorphous) phase:

  | Phase       | E   | S   | Neighbours | Description |
  |-------------|-----|-----|------------|---------------------------------------------------------|
  | Cubic       | 0   | 4   | 4          | All 4 neighbours are staggered. | |
Hexagonal   | 1   | 3   | 4          | 3 staggered and 1 eclipsed bond. | |
Interfacial | any | 2   | 4          | Belongs to the ice phase, but does not
satisfy strict criteria of cubic and hexagonal | |             | 0   | 3   | 4
|                                                         | | Liquid      | N/A
| N/A | any        | Classification by exclusion.                            |

  Here, E and S refer to eclipsed and staggered bonds. The neighbours column
contains the number of nearest neighbours.

  The CHILL+ algorithm, which is a modified version of CHILL, additionally
identifies clathrate and interfacial clathrate phases. The criteria are
  enumerated in the table below:

  | Phase                 | E   | S   | Neighbours | Description |
  |-----------------------|-----|-----|------------|---------------------------------------------------------------------------------------------------|
  | Cubic                 | 0   | 4   | 4          | no change was made to
staggered bond criterion from CHILL                                         | |
Hexagonal             | 1   | 3   | 4          | wider eclipsed range compared
to CHILL; identifies 99% hexagonal ice up to 270 K                  | |
Interfacial           | any | 2   | 4          | must have at least one first
neighbor water with more than two staggered bonds                    | | | 0   |
3   | 4          | must have at least one first neighbor water with more than
one staggered bond                     | | Clathrate             | 4   | 0   | 4
| bulk clathrate and part of the interface of clathrates can be identified with
four eclipsed bonds | | Interfacial clathrate | 3   | any | 4          | partial
clathrate cages and threads of clathrate-like order | | Liquid                |
N/A | N/A | any        | classifies as liquid if none of the above criteria are
fulfilled                                  |

  Although both the CHILL and CHILL+ classification schemes take into account
the local order of the environment of each particle, both schemes output a
per-particle classification.

  ### Changelog ###

  - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Sept 19, 2019
 */

namespace chill {

// 2*l+1 length complex vector
/** \struct YlmAtom
 * \brief This contains a complex vector of length \f$2l+1\f$.
 *
 * Contains specifically:
 * - A complex vector, projecting the orientational structure of a single pair
 * with each of the four nearest neighbours, on the basis of spherical harmonics
 */
struct YlmAtom {
  std::vector<std::complex<double>> ylm;
};

// Vector of 2*l+1 averaged over 4 nearest neighbours
/** @struct QlmAtom
 * @brief This is the local orientational bond order parameter @f$q_{lm}@f$, of
 * length @f$2l+1@f$.
 *
 * This complex vector is averaged over the four nearest neighbours, according
 * to the following equation:
 *
 * @f[
 * q_{lm}(i) = \frac{1}{N_b(i)} \Sigma_{j=1}^{N_b(i)} Y_{lm}(r_{ij})
 * @f]
 *
 * Here, @f$N_b(i)=4@f$ is the number of nearest neighbours for the molecule
 * @f$i@f$. This struct contains specifically:
 * - A complex vector of length @f$2l+1@f$, calculated according to the equation
 * above
 */
struct QlmAtom {
  std::vector<YlmAtom> ptq; // Averaged over neighbours
};

/**
 *  Function for getting the bond order correlations @f$c_{ij}@f$  (or
 @f$a_{ij}@f$ in some treatments) according to the CHILL algorithm.
 *  @param[in,out] yCloud The output molSys::PointCloud
 *  @param[in] nList The row-ordered neighbour list, by ID.
 The first element of each row is the particle ID, followed by the IDs of the
 neighbours
 *  @param[in] isSlice This decides whether there is a slice or not
 */
molSys::PointCloud<molSys::Point<double>, double>
getCorrel(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
          std::vector<std::vector<int>> nList, bool isSlice = false);

/**
 *  Function that classifies every particle's #molSys::atom_state_type ice
 type, according to the CHILL algorithm. Does not print out the information.
 *  @param[in,out] yCloud The output molSys::PointCloud
 *  @param[in] isSlice This decides whether there is a slice or not
 *  @param[in] nList Row-ordered neighbour list by atom ID
 */
molSys::PointCloud<molSys::Point<double>, double>
getIceTypeNoPrint(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                  std::vector<std::vector<int>> nList, bool isSlice = false);

// Classifies each atom according to the CHILL algorithm
/**
 *  Function that classifies every particle's #molSys::atom_state_type ice
 type, according to the CHILL algorithm.
 *  @param[in,out] yCloud The output molSys::PointCloud
 *  @param[in] nList Row-ordered neighbour list by atom ID
 *  @param[in] path Path to the output directory to which ice types are written
 out to
 *  @param[in] firstFrame First frame to be analyzed
 *  @param[in] isSlice This decides whether there is a slice or not
 *  @param[in] outputFileName Name of the output file, to which the ice types
 will be written out.
 */
molSys::PointCloud<molSys::Point<double>, double>
getIceType(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
           std::vector<std::vector<int>> nList, std::string path,
           int firstFrame, bool isSlice = false,
           std::string outputFileName = "chill.txt");

//! Gets c_ij and then classifies bond types according to the CHILL+ algorithm
molSys::PointCloud<molSys::Point<double>, double>
getCorrelPlus(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
              std::vector<std::vector<int>> nList, bool isSlice = false);

//! Classifies each atom according to the CHILL+ algorithm
molSys::PointCloud<molSys::Point<double>, double>
getIceTypePlus(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
               std::vector<std::vector<int>> nList, std::string path,
               int firstFrame, bool isSlice = false,
               std::string outputFileName = "chillPlus.txt");

//! q6 can distinguish between water and ice. Use this for the largest ice
//! cluster
std::vector<double>
getq6(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
      std::vector<std::vector<int>> nList, bool isSlice = false);

//! 'Test' condition for classifying hexagonal ice using averaged q6 and q3
//! Checks water
//! According to https://!pubs.rsc.org/en/content/articlehtml/2011/cp/c1cp22167a
//! Gets c_ij and then classifies bond types according to the CHILL+ algorithm
molSys::PointCloud<molSys::Point<double>, double>
reclassifyWater(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                std::vector<double> *q6);

//! Prints out the iceType for a particular frame onto the terminal
int printIceType(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 std::string path, int firstFrame, bool isSlice = false,
                 std::string outputFileName = "superChill.txt");

//! Checks if a given iatom is interfacial ice or not, according to the CHILL
//! algorithm
bool isInterfacial(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                   std::vector<std::vector<int>> nList, int iatom,
                   int num_staggrd, int num_eclipsd);

//! Finds the number of staggered bonds for a given atom of index jatom
int numStaggered(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 std::vector<std::vector<int>> nList, int jatom);

} // namespace chill

/** \brief Functions used for spherical harmonics
 *
 * The recommended function for calculating the spherical-harmonics based
 * complex vector is sph::spheriHarmo. This function uses the
 * [Boost](https://www.boost.org/) libraries.
 *
 * ### Changelog ###
 *
 * - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Sept 19, 2019
 * - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */
namespace sph {

// 7 is for Q3, orderL=3

std::vector<std::complex<double>>
spheriHarmo(int orderL, std::array<double, 2> radialCoord);

std::array<double, 2> radialCoord(std::array<double, 3> cartCoord);

//! Lookup table for Q3
std::vector<std::complex<double>>
lookupTableQ3Vec(std::array<double, 2> angles);

//! Lookup table for Q3 (m=0 to m=6)
std::complex<double> lookupTableQ3(int m, std::array<double, 2> angles);

//! Lookup table for Q6
std::vector<std::complex<double>>
lookupTableQ6Vec(std::array<double, 2> angles);

//! Lookup table for Q6 (m=0 to m=12)
std::complex<double> lookupTableQ6(int m, std::array<double, 2> angles);

} // namespace sph

#endif // __BOP_H_
