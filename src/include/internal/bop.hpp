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

#ifndef SEAMS_BOP_H_
#define SEAMS_BOP_H_

#include <array>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_MPI
#include <mpi.h>
#endif
#include <generic.hpp>
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
  a(i,j) = \frac{\sum_{m=-l}^{l} q_{lm}(i)\, q_{lm}^*(j)}
  {|q_l(i)|\,|q_l(j)|}
  @f]

  where, @f$q_{lm}^*@f$ is the complex conjugate of @f$q_{lm}@f$.

  Depending on the values of @f$a(i,j)@f$, the bond between each pair of
molecules @f$i@f$ and @f$j@f$ is classified as being either staggered or
eclipsed.

  @code
  Algorithm   Eclipsed bonds              Staggered bonds
  CHILL       -0.05 >= a(i,j) >= -0.2     a(i,j) <= -0.8
  CHILL+      0.25 >= a(i,j) >= -0.35     a(i,j) <= -0.8
  @endcode

  Since each molecule @f$i@f$ in deeply supercooled water has four nearest
neighbours, the type of the resultant four bonds is used to identify phases
according to the CHILL and CHILL+ algorithms.

  The CHILL algorithm can classify water molecules as belonging to the cubic,
hexagonal, interfacial or liquid (amorphous) phase:

  @code
  Phase        E    S    Neighbours  Description
  Cubic        0    4    4           All 4 neighbours are staggered.
  Hexagonal    1    3    4           3 staggered and 1 eclipsed bond.
  Interfacial  any  2    4           Belongs to the ice phase, but does not satisfy strict criteria of cubic and hexagonal.
               0    3    4
  Liquid       N/A  N/A  any         Classification by exclusion.
  @endcode

  Here, E and S refer to eclipsed and staggered bonds. The neighbours column
contains the number of nearest neighbours.

  The CHILL+ algorithm, which is a modified version of CHILL, additionally
identifies clathrate and interfacial clathrate phases. The criteria are
  enumerated in the table below:

  @code
  Phase                  E    S    Neighbours  Description
  Cubic                  0    4    4           no change was made to the staggered bond criterion from CHILL
  Hexagonal              1    3    4           wider eclipsed range compared to CHILL; identifies 99% hexagonal ice up to 270 K
  Interfacial            any  2    4           must have at least one first neighbor water with more than two staggered bonds
                         0    3    4           must have at least one first neighbor water with more than one staggered bond
  Clathrate              4    0    4           bulk clathrate and part of the interface of clathrates can be identified with four eclipsed bonds
  Interfacial clathrate  3    any  4           partial clathrate cages and threads of clathrate-like order
  Liquid                 N/A  N/A  any         classifies as liquid if none of the above criteria are fulfilled
  @endcode

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

/** @struct BondClassifier
 * @brief One rule set for classifying bond correlations @f$c_{ij}@f$.
 * @details A bond is staggered when @f$c_{ij} \le@f$ staggeredMax, eclipsed
 *  when it falls in [eclipsedMin, eclipsedMax], and out of range otherwise.
 *  The CHILL and CHILL+ water rules are two instances; other hydrogen-bonded
 *  or tetrahedral materials register their own thresholds and coordination.
 */
struct BondClassifier {
  double staggeredMax;    //! Staggered edge. Inclusive when closed is true
  double eclipsedMin;     //! Lower edge of the eclipsed band
  double eclipsedMax;     //! Upper edge of the eclipsed band
  int coordinationNumber; //! Nearest neighbours in the bond sum; <= 0 keeps
                          //! each atom's whole neighbour row
  bool closed = true;     //! Inclusive edges. OVITO Chill+ uses exclusive edges
};

//! Staggered, eclipsed, or out of range. Closed edges match Nguyen and
//! Molinero (c <= staggeredMax, eclipsedMin <= c <= eclipsedMax). Open
//! edges match OVITO 3.17.2 (c < staggeredMax, eclipsedMin < c < eclipsedMax).
inline molSys::bond_type classifyBond(double cij, const BondClassifier &rule) {
  const bool staggered =
      rule.closed ? cij <= rule.staggeredMax : cij < rule.staggeredMax;
  if (staggered) {
    return molSys::bond_type::staggered;
  }
  const bool eclipsed =
      rule.closed ? (cij >= rule.eclipsedMin && cij <= rule.eclipsedMax)
                  : (cij > rule.eclipsedMin && cij < rule.eclipsedMax);
  if (eclipsed) {
    return molSys::bond_type::eclipsed;
  }
  return molSys::bond_type::out_of_range;
}

//! The CHILL rule set (Moore et al., PCCP 12, 4124, 2010): strict staggered
//! bound, narrow eclipsed band, four nearest neighbours
[[nodiscard]] BondClassifier chillRule();

//! The CHILL+ rule set (Nguyen and Molinero, JPCB 119, 9369, 2015): wider
//! eclipsed band for interfacial and clathrate recognition
[[nodiscard]] BondClassifier chillPlusRule();

//! Look up a registered rule set by name. "CHILL" and "CHILL+" are always
//! present; registerBondClassifier adds more. Throws std::out_of_range for
//! unknown names.
[[nodiscard]] BondClassifier bondClassifier(const std::string &name);

//! Register (or replace) a named rule set at runtime, making the material's
//! thresholds available to every front end by name
void registerBondClassifier(const std::string &name,
                            const BondClassifier &rule);

//! Names of every registered rule set
[[nodiscard]] std::vector<std::string> bondClassifierNames();

//! Compute and classify the bond correlations c_ij under an arbitrary rule
//! set, filling yCloud.pts[i].c_ij. getCorrel and getCorrelPlus are this
//! engine under their canonical water rules.
void classifyBonds(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                   const std::vector<std::vector<int>> &nList,
                   const BondClassifier &rule, bool isSlice = false);

//! Unnormalised real l = 3 bond-order vector of the atoms @a which selects.
//! @a which null computes every atom. Ghost slots stay zero so a neighbour
//! rank can write the vector it computed for that atom.
struct BondOrder {
  int keepMax = 0;
  std::vector<double> q;    //! atom-major, 7 channels
  std::vector<double> norm; //! squared norm of q, with the m != 0 factor
  std::vector<int> nkeep;
  std::vector<int> chosen;             //! keepMax > 0, local indices, -1 empty
  std::vector<std::vector<int>> rows;  //! keepMax <= 0, local indices
};

void bondOrder(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               const std::vector<std::vector<int>> &nList, const BondClassifier &rule,
               BondOrder &order, const std::vector<char> *which);

//! c_ij for atoms @a which selects, from a bond-order vector that already
//! holds every neighbour the bond reads. Atoms outside @a which are left
//! untouched, including a ghost's bond list.
void correlFromOrder(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     const BondOrder &order, const BondClassifier &rule, bool isSlice,
                     const std::vector<char> *which);

#ifdef SEAMS_HAS_MPI
//! Published or OVITO CHILL+ on the owned atoms. Coordinates move one
//! cutoff; the bond-order vector of the boundary moves in one exchange,
//! overlapped with the interior correlations. The published neighbour test
//! then exchanges the staggered count of that same boundary. iceType is
//! written back onto the owned atoms of @a full.
void distributedChillPlus(molSys::PointCloud<molSys::Point<double>, double> &full,
                          double cutoff, BondClassifier rule, bool neighbourTest,
                          MPI_Comm comm);
#endif

/**
 *  Function for getting the bond order correlations @f$c_{ij}@f$  (or
 @f$a_{ij}@f$ in some treatments) according to the CHILL algorithm.
 *  @param[in,out] yCloud The output molSys::PointCloud
 *  @param[in] nList The row-ordered neighbour list, by ID.
 The first element of each row is the particle ID, followed by the IDs of the
 neighbours
 *  @param[in] isSlice This decides whether there is a slice or not
 *  @param[in] coordinationNumber How many nearest neighbours the bond sum
 runs over. The default of four is the coordination CHILL is defined and
 validated against (tetrahedral water); other values reuse the c_ij machinery
 for differently coordinated systems, and a non-positive value gives each
 atom its own count from its neighbour-list row. The ice-classification
 tables in getIceType remain four-bond rules
 */
void
getCorrel(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
          const std::vector<std::vector<int>> &nList, bool isSlice = false,
          int coordinationNumber = 4);

/**
 *  Function that classifies every particle's #molSys::atom_state_type ice
 type, according to the CHILL algorithm. Does not print out the information.
 *  @param[in,out] yCloud The output molSys::PointCloud
 *  @param[in] isSlice This decides whether there is a slice or not
 *  @param[in] nList Row-ordered neighbour list by atom ID
 */
void
getIceTypeNoPrint(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                  const std::vector<std::vector<int>> &nList, bool isSlice = false);

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
void
getIceType(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           const std::vector<std::vector<int>> &nList, std::string path,
           int firstFrame, bool isSlice = false,
           std::string outputFileName = "chill.txt");

//! Gets c_ij and then classifies bond types according to the CHILL+
//! algorithm. coordinationNumber as in getCorrel: four is the validated
//! CHILL+ water scheme, non-positive keeps each atom's whole neighbour row
void
getCorrelPlus(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
              const std::vector<std::vector<int>> &nList, bool isSlice = false,
              int coordinationNumber = 4, bool closedBounds = true);

//! Published CHILL+ on the cutoff shell. The pair rows are read inside
//! this call; the shared neighbour-list builder is not modified.
//! Returns false when that search is not available.
bool chillPlusPublishedCutoff(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
    int typeI);

//! One cutoff list and one bond pass for local q4, local q6, their
//! Lechner-Dellago averages, and the published CHILL+ labels. The four
//! q arrays match steinhardtQlCutoff with averaging on. Returns false
//! when the atoms are not one type or the cutoff list is not available,
//! and then writes nothing.
bool q4q6ChillPlus(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
    int typeI, std::vector<double> &q4, std::vector<double> &q6,
    std::vector<double> &q4bar, std::vector<double> &q6bar);

//! Cubic pre-filter from discovery/pysr_prefilter.py. The selected
//! expression is cmean^4. It writes cubic and returns true only when
//! every atom has a cutoff shell of four and a negative mean correlation
//! whose fourth power sits above the highest non-cubic score in that
//! search. That forces every bond below -0.84, which the published rule
//! labels staggered, so the atom is cubic. Otherwise it returns false
//! and leaves the cloud unchanged. The default published call does not
//! use it.
bool chillPlusSymbolicPrefilter(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
    int typeI);

//! Classifies each atom according to the CHILL+ algorithm
void
getIceTypePlus(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               const std::vector<std::vector<int>> &nList, std::string path,
               int firstFrame, bool isSlice = false,
               std::string outputFileName = "chillPlus.txt",
               bool neighbourTest = true);

//! CHILL+ ice types on the cloud. Does not write a file.
//! neighbourTest false is the OVITO interfacial rule (two staggered
//! bonds, or three staggered and no eclipsed bond, with no neighbour
//! lookup). Either rule still requires a cutoff shell of four.
void
getIceTypePlusNoPrint(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, bool isSlice = false,
    bool neighbourTest = true);

//! Set iceType to unclassified for atoms whose LAMMPS type is not in types.
//! Empty types leaves the cloud unchanged. Used after a water-type k-NN so
//! substrate and ions are not scored as water.
void clearIceTypesExcept(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &types);

//! q6 can distinguish between water and ice. Use this for the largest ice
//! cluster
std::vector<double>
getq6(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
      const std::vector<std::vector<int>> &nList, bool isSlice = false);

//! 'Test' condition for classifying hexagonal ice using averaged q6 and q3
//! Checks water
//! According to https://!pubs.rsc.org/en/content/articlehtml/2011/cp/c1cp22167a
//! Gets c_ij and then classifies bond types according to the CHILL+ algorithm
void
reclassifyWater(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                std::vector<double> &q6);

//! Prints out the iceType for a particular frame onto the terminal
[[nodiscard]] int printIceType(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                 std::string path, int firstFrame, bool isSlice = false,
                 std::string outputFileName = "superChill.txt");

//! Interfacial ice. Moore CHILL uses a neighbour with exactly two staggered
//! bonds; Nguyen CHILL+ uses a neighbour with more than one. Pass
//! chillPlus=true for the latter.
bool isInterfacial(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                   const std::vector<std::vector<int>> &nList, int iatom,
                   int num_staggrd, int num_eclipsd, bool chillPlus = false);

//! Finds the number of staggered bonds for a given atom of index jatom
[[nodiscard]] int numStaggered(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                 const std::vector<std::vector<int>> &nList, int jatom);

/** @struct SteinhardtQl
 * @brief Per-particle Steinhardt order parameters of a single degree @f$l@f$.
 * @details Holds both the local parameter and its neighbour-averaged
 *  counterpart. The averaged form of Lechner and Dellago folds in the second
 *  coordination shell, which separates the distributions of competing
 *  structures far more cleanly than the local form alone.
 */
struct SteinhardtQl {
  std::vector<double> ql;    //! Local @f$q_l(i)@f$
  std::vector<double> qlBar; //! Averaged @f$\bar{q}_l(i)@f$
  std::vector<double> wl;    //! Normalised @f$\hat{w}_l@f$, empty unless requested
};

//! Local and neighbour-averaged Steinhardt parameters of degree 3, 4, 6,
//! 8, or 12. Local ql is Steinhardt, Nelson and Ronchetti; qlBar is the
//! Lechner-Dellago average of q_lm over the particle and its neighbours.
//! l=12 is host-only (associated Legendre, or sphericart when linked)
//! for the Zeron hydrate pair with l=3; the device Ylm path has no
//! l=12. Other degrees return zeros. The compute
//! path flattens the neighbour list to CSR so the same kernel can run on
//! the host (OpenMP), across MPI ranks (atom split plus Allgatherv of
//! q_lm), or under OpenMP target offload when the build provides a device.
[[nodiscard]] SteinhardtQl
steinhardtQl(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
             const std::vector<std::vector<int>> &nList, int orderL);

//! Local @f$q_l@f$ on the cutoff list.
//!
//! Neighbours are orthoFullCsr, not a four-nearest rebuild, and the
//! kernel reads that CSR directly. Lechner-Dellago @f$\bar{q}_l@f$ is
//! filled only when @a average is true. Otherwise `qlBar` stays zero.
//! @a invariants fills @c wl with the normalised Steinhardt
//! @f$\hat{w}_l@f$ for l = 4 or 6.
[[nodiscard]] SteinhardtQl steinhardtQlCutoff(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double cutoff, int orderL, int typeI, bool average, bool invariants = false);

//! Mutual k-nearest Steinhardt parameters. The coordinate buffer built
//! for the nomination is the one the harmonics read, and the neighbour
//! rows stay in one packed array. l is 3, 4 or 6; other degrees return
//! the result of building the mutual graph and calling steinhardtQl.
[[nodiscard]] SteinhardtQl steinhardtQlMutual(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, int typeI, int orderL);

//! Whether steinhardtQl splits atoms over the ranks of an initialised
//! MPI_COMM_WORLD (true unless set). A driver that gives each rank whole
//! frames of its own turns it off.
void setSteinhardtAtomSplit(bool split);

//! Width of one atom's interleaved q_lm: real, imaginary, m = -l..l.
inline int steinhardtQlmWidth(int orderL) { return 2 * (2 * orderL + 1); }

//! q_lm for atoms @a which selects. Other atoms are left zero.
void writeSteinhardtQlm(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                        const std::vector<std::vector<int>> &nList, int orderL,
                        std::vector<double> &qlm, const std::vector<char> *which);

//! ql and qlBar for atoms @a which selects, reading q_lm of their neighbours.
void readSteinhardtQl(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                      const std::vector<std::vector<int>> &nList, int orderL,
                      const std::vector<double> &qlm, std::vector<double> &ql,
                      std::vector<double> &qlBar, const std::vector<char> *which);

#ifdef SEAMS_HAS_MPI
//! q_l and q_l bar of the owned atoms, frame-indexed. One cutoff of
//! coordinates and one exchange of q_lm.
SteinhardtQl distributedSteinhardt(
    const molSys::PointCloud<molSys::Point<double>, double> &full, double cutoff,
    int orderL, MPI_Comm comm);
#endif

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

//! Fill `out` with Y_lm, m = -l .. +l. Reuses `out` capacity.
void spheriHarmoInto(int orderL, std::array<double, 2> radialCoord,
                     std::vector<std::complex<double>> &out);

std::array<double, 2> radialCoord(std::array<double, 3> cartCoord);

//! Lookup table for Q3
std::vector<std::complex<double>>
lookupTableQ3Vec(std::array<double, 2> angles);

//! Lookup table for Q3 (m=0 to m=6)
std::complex<double> lookupTableQ3(int m, std::array<double, 2> angles);

//! Lookup table for Q4
std::vector<std::complex<double>>
lookupTableQ4Vec(std::array<double, 2> angles);

//! Lookup table for Q4 (m=0 to m=8)
std::complex<double> lookupTableQ4(int m, std::array<double, 2> angles);

//! Lookup table for Q6
std::vector<std::complex<double>>
lookupTableQ6Vec(std::array<double, 2> angles);

//! Lookup table for Q8
std::vector<std::complex<double>>
lookupTableQ8Vec(std::array<double, 2> angles);

//! Lookup table for Q8 (m=0 to m=16)
std::complex<double> lookupTableQ8(int m, std::array<double, 2> angles);

//! Lookup table for Q6 (m=0 to m=12)
std::complex<double> lookupTableQ6(int m, std::array<double, 2> angles);

} // namespace sph

#endif // SEAMS_BOP_H_
