#ifndef __CLUSTER_H_
#define __CLUSTER_H_

#include <math.h>
#include <array>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <bop.hpp>
#include <cmath>
#include <complex>
#include <generic.hpp>
#include <iostream>
#include <mol_sys.hpp>
#include <neighbours.hpp>

/*! \file cluster.hpp
    \brief File for the bond order parameter analysis.

    Details.
*/

/*!
 *  \addtogroup clump
 *  @{
 */

/*! \brief Clustering functions.
 *         This namespace contains functions that are used in the
for clustering ice-like particles.
 *
 Cluster analysis

  ### Changelog ###

  - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 4, 2019
 */

namespace clump {

// Finds the largest ice cluster
int largestIceCluster(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
    std::vector<std::vector<int>> nList, std::vector<bool> *isIce,
    std::vector<int> *clusterID, std::vector<int> *nClusters,
    std::unordered_map<int, int> *indexNumber);

// Does the cluster analysis of ice particles in the system. Returns a
// pointCloud of the largest ice cluster. The neighbour list returned is BY
// INDEX of the largest ice cluster pointCloud.
int clusterAnalysis(molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
                    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::vector<std::vector<int>> nList,
                    std::vector<std::vector<int>> &iceNeighbourList,
                    double cutoff, std::string bopAnalysis = "q6");

// Recenters the coordinates of a pointCloud
int recenterClusterCloud(
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud);

}  // namespace clump

#endif  // __CLUSTER_H_
