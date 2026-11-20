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

#ifndef SEAMS_BONDING_H_
#define SEAMS_BONDING_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// Internal
#include <cage.hpp>
#include <mol_sys.hpp>
#include <seams_input.hpp>

/** @file bond.hpp
 *  @brief File for bond-related analyses (hydrogen bonds, bonded atoms for data
 * file write-outs etc.).
 */

/**
 *  @addtogroup bond
 *  @{
 */

/** @brief Functions for bond-related analyses
 *
 *  @details This namespace contains functions that are used for determining
 the hydrogen bonds, or lists of bonded atoms for write-outs to data files and
 visualization.
 *
 * These functions are distinct from the neighbour-list building functions.
 Bonds can be built from cages, ring lists etc.
 *
 * Hydrogen bonds are determined using a strict geometric criterion.
 *
 * A hydrogen bond between two water molecules exists when:
 *
 * 1. The distance between the acceptor oxygen atom and the donor hydrogen atom
 is less than 2.42 Angstrom
 * 2. The angle between the O-O vector and the O-H vector should be less than
 30 degrees
 *
 * Those numbers are the water defaults on populateHbonds and
 * populateHbondsWithInputClouds. populateHbondsFromDonors takes the same
 * defaults so a water donor-H set stays comparable; other chemistries pass
 * their own cutoffs. They are not an ionic-liquid default.
 *
 * ### Changelog ###
 *
 * - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 13, 2019
 * - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

namespace bond {

//! Create a vector of vectors (similar to the neighbour list conventions)
//! containing the hydrogen bond connectivity information. Decides the
//! existence of the hydrogen bond depending on the O--O and O--H vectors from
//! the neighbour list already constructed. The geometric criterion defaults
//! to the water values (acceptor--H distance below 2.42 Angstrom, O-O-H angle
//! below 30 degrees); other hydrogen-bonded systems pass their own thresholds
std::vector<std::vector<int>>
populateHbonds(std::string filename,
               molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               const std::vector<std::vector<int>> &nList, int targetFrame,
               int Htype, double distCutoff = 2.42, double angleCutoff = 30.0);

//! Create a vector of vectors (similar to the neighbour list conventions)
//! containing the hydrogen bond connectivity information. Decides the
//! existence of the hydrogen bond depending on the O--O and O--H vectors from
//! the neighbour list already constructed, taking a PointCloud for the H atoms
//! as input
// ! The H atom PointCloud should be for the entire system
std::vector<std::vector<int>> populateHbondsWithInputClouds(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    molSys::PointCloud<molSys::Point<double>, double> &hCloud,
    const std::vector<std::vector<int>> &nList, double distCutoff = 2.42,
    double angleCutoff = 30.0);

//! Geometric test for one donor-acceptor assignment. donorHs are hCloud
//! indices on the donor. O-O is gen::relDist(yCloud, acceptor, donor).
//! O-H uses gen::relDistFromPoint (dump MIC).
bool donatedHydrogenBond(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const molSys::PointCloud<molSys::Point<double>, double> &hCloud,
    int acceptorIndex, int donorIndex, const std::vector<int> &donorHs,
    double distCutoff = 2.42, double angleCutoff = 30.0);

//! Hydrogen-bond network from an explicit donor-H set. donorHs are hCloud
//! indices; each is paired with the heavy atom that shares its molID.
//! Water callers build donorHs from hAtomMolList (two H per molecule).
std::vector<std::vector<int>> populateHbondsFromDonors(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    molSys::PointCloud<molSys::Point<double>, double> &hCloud,
    const std::vector<std::vector<int>> &nList,
    const std::vector<int> &donorHs, double distCutoff = 2.42,
    double angleCutoff = 30.0);

//! One directed hydrogen bond. Indices are positions in the single cloud
//! passed to directedHydrogenBonds. Angle is degrees, acceptor-first:
//! the two vectors are the minimum image of acceptor minus donor and the
//! minimum image of acceptor minus the hydrogen. A linear
//! donor-hydrogen-acceptor geometry has angle 0.
struct DirectedHydrogenBond {
  int donor = -1;
  int hydrogen = -1;
  int acceptor = -1;
  double d_AH = 0.0;
  double d_DA = 0.0;
  double angle = 0.0;
};

//! Directed donor-hydrogen-acceptor edges on one cloud.
//!
//! Atoms whose type equals hydrogenType are hydrogens. Every other atom is
//! a heavy atom. Each hydrogen is given to the nearest heavy atom of a
//! donor type that shares its molID; the lower index wins an equal
//! distance. A hydrogen with no such heavy atom donates nothing. An edge is
//! kept when the acceptor-hydrogen minimum image is strictly shorter than
//! distCutoff and the acceptor-first angle is less than or equal to
//! angleCutoff. The neighbour list keeps a donor-acceptor pair only when
//! its minimum image is strictly shorter than pairCutoff. excludeSameMol
//! drops a pair that shares a molID. An empty donor or acceptor type list
//! returns no edges. Results are sorted by donor, hydrogen, acceptor.
std::vector<DirectedHydrogenBond> directedHydrogenBonds(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    double pairCutoff, double distCutoff, double angleCutoff,
    const std::vector<int> &donorTypes, const std::vector<int> &acceptorTypes,
    int hydrogenType, bool excludeSameMol = true);

//! Calculates the distance of the hydrogen bond between O and H (of different
//! atoms), given the respective pointClouds and the indices to each atom
double getHbondDistanceOH(
    const molSys::PointCloud<molSys::Point<double>, double> &oCloud,
    const molSys::PointCloud<molSys::Point<double>, double> &hCloud,
    int oAtomIndex, int hAtomIndex);

//! Create a vector of vectors containing bond connectivity information. May
//! contain duplicates! Gets the bond information from the vector of vectors
//! containing the neighbour list by index
std::vector<std::vector<int>>
populateBonds(const std::vector<std::vector<int>> &nList,
              const molSys::PointCloud<molSys::Point<double>, double> &yCloud);

//! Create a vector of vectors containing bond connectivity information
//! Gets the bond information from the vector of vectors
//! containing the neighbour list by index. Bonds between dummy atoms are not
//! filled.
std::vector<std::vector<int>>
populateBonds(const std::vector<std::vector<int>> &nList,
              const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
              const std::vector<cage::iceType> &atomTypes);

//! Creates a vector of vectors containing bond connectivity information from
//! the rings vector of vectors and cage information
std::vector<std::vector<int>>
createBondsFromCages(const std::vector<std::vector<int>> &rings,
                     std::vector<cage::Cage> &cageList, cage::cageType type,
                     int &nRings);

//! Remove duplicate bonds
std::vector<std::vector<int>> trimBonds(std::vector<std::vector<int>> bonds);

} // namespace bond

/** @} */
#endif // SEAMS_BONDING_H_
