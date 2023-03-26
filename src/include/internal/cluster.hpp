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

#ifndef __CLUSTER_H_
#define __CLUSTER_H_

#include <array>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <bop.hpp>
#include <cmath>
#include <complex>
#include <generic.hpp>
#include <iostream>
#include <math.h>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_output.hpp>

/** @file cluster.hpp
 *   @brief File for the bond order parameter analysis.
 *
 */

/**
 *  @addtogroup clump
 *  @{
 */

/*! @brief Clustering functions.
 *         This namespace contains functions that are used in the
 *  for clustering ice-like particles.
 *
 * Particles are clustered according to Stoddard's algorithm
 * [(Stoddard J. Comp. Phys., 27, 291,
 * 1977)](https://www.sciencedirect.com/science/article/pii/0021999178900116).
 *
 * The largest ice cluster determines the largest cluster of ice-like particles,
 * using a "linked list" created by the clustering algorithm. The determination
 * of ice-like molecules can be done using the @f$ q_6 @f$  parameter, CHILL or
 * CHILL+ parameters.
 *
 * The largest ice cluster thus determined can be re-centered for ease of
 * visualization.
 *
 *  ### Changelog ###
 *
 * - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 13, 2019
 * - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

namespace clump {

//! Finds the largest ice cluster
int largestIceCluster(
    std::string path, molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
    std::vector<std::vector<int>> nList, std::vector<bool> *isIce,
    std::vector<int> *clusterID, std::vector<int> *nClusters,
    std::unordered_map<int, int> *indexNumber, int firstFrame);

//! Get the linked list of a cluster, given by iceCloud, for a single cluster.
//! Required for cluster re-centering
int singleClusterLinkedList(
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
    std::vector<std::vector<int>> nList, std::vector<int> *linkedList);

//! Does the cluster analysis of ice particles in the system. Returns a
//! pointCloud of the largest ice cluster. The neighbour list returned is BY
//! INDEX of the largest ice cluster pointCloud.
int clusterAnalysis(std::string path,
                    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
                    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::vector<std::vector<int>> nList,
                    std::vector<std::vector<int>> &iceNeighbourList,
                    double cutoff, int firstFrame,
                    std::string bopAnalysis = "q6");

//! Recenters the coordinates of a pointCloud
int recenterClusterCloud(
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
    std::vector<std::vector<int>> nList);

} // namespace clump

#endif // __CLUSTER_H_
