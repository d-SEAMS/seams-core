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

#ifndef __PNTCORRESPONDENCE_H_
#define __PNTCORRESPONDENCE_H_

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

// External
#include <Eigen/Core>

#include <cage.hpp>
#include <mol_sys.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

namespace pntToPnt {

//! Fills up an eigen matrix point set a reference ring, which is a regular
//! n-gonal polygon; where n is the number of nodes in the ring
Eigen::MatrixXd getPointSetRefRing(int n, int axialDim);

//! Creates an eigen matrix for the points of a prism block, constructed from
//! the points of a perfect polygon of radius 1, given the basal rings and
//! axial dimension
Eigen::MatrixXd
createPrismBlock(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 const Eigen::MatrixXd &refPoints, int ringSize,
                 std::vector<int> basal1, std::vector<int> basal2);

//! Calculate the average radial distance for the basal rings, calculated from
//! the centroid of each basal ring
double
getRadiusFromRings(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                   std::vector<int> basal1, std::vector<int> basal2);

//! Calculate the average height of the prism block, calculated using the basal
//! rings of the prism and the axial dimension
double getAvgHeightPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2);

//! Get the relative ordering of a pair of basal rings for a deformed
//! prism/perfect prism. Outputs a vector of vectors of indices, such that the
//! first vector is for the first basal ring, and the second vector is for the
//! second basal ring. The input neighbour list is with respect to indices, not
//! IDs
int relOrderPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2,
    std::vector<std::vector<int>> nList, std::vector<int> *outBasal1,
    std::vector<int> *outBasal2);

//! Get the relative ordering of a pair of basal rings for a deformed
//! prism/perfect prism. Outputs a vector of vectors of indices, such that the
//! first vector is for the first basal ring, and the second vector is for the
//! second basal ring.
int relOrderPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2,
    std::vector<int> *outBasal1, std::vector<int> *outBasal2);

//! Fill up an Eigen Matrix of a prism basal ring from an input vector of atom
//! indices
Eigen::MatrixXd
fillPointSetPrismRing(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                      std::vector<int> basalRing, int startingIndex);

//! Fill up an Eigen matrix of a prism block (two basal rings) from input
//! vectors for the basal rings
Eigen::MatrixXd fillPointSetPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2, int startingIndex);

//! Fills up an eigen matrix point set a reference ring, which is a regular
//! n-gonal polygon, constructed with radius 1 by default; where n is the
//! number of nodes in the ring
Eigen::MatrixXd getPointSetCage(ring::strucType type);

//! Matches the order of the basal rings of an HC or a potential HC
int relOrderHC(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
               std::vector<int> basal1, std::vector<int> basal2,
               std::vector<std::vector<int>> nList,
               std::vector<int> *matchedBasal1,
               std::vector<int> *matchedBasal2);

//! Matches the order of the basal rings of an DDC or a potential HC
std::vector<int> relOrderDDC(int index, std::vector<std::vector<int>> rings,
                             std::vector<cage::Cage> cageList);

//! Fills up an eigen matrix point set using the basal rings basal1 and basal2,
//! changing the order of the point set by filling up from the startingIndex
//! (starting from 0 to 5)
Eigen::MatrixXd
changeHexCageOrder(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                   std::vector<int> basal1, std::vector<int> basal2,
                   int startingIndex = 0);

//! Fills up an eigen matrix point set using information of the equatorial ring
//! and peripheral rings, embedded in a vector, already filled in relOrderDDC.
Eigen::MatrixXd
changeDiaCageOrder(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                   std::vector<int> ddcOrder, int startingIndex = 0);

} // namespace pntToPnt

#endif // __PNTCORRESPONDENCE_H_
