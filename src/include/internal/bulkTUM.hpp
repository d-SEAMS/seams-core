#ifndef __BULKTUM_H_
#define __BULKTUM_H_

#include <math.h>
#include <sys/stat.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <franzblau.hpp>
#include <neighbours.hpp>
#include <pntCorrespondence.hpp>
#include <topo_bulk.hpp>

/*! \file bulkTUM.hpp
    \brief File containing functions used specific to bulk topological
   unit matching (TUM) criterion
*/

/*!
 *  \addtogroup ring
 *  @{
 */

namespace tum3 {

// Topological unit matching for bulk water. If printClusters is true,
// individual clusters of connected cages are printed.
int topoUnitMatchingBulk(
    std::string path, std::vector<std::vector<int>> rings,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int firstFrame,
    bool printClusters, bool onlyTetrahedral);

// Build a reference Hexagonal cage, reading in from a template XYZ file
Eigen::MatrixXd buildRefHC(std::string fileName);

// Build a reference Double-Diamond cage, reading in from a template XYZ file
Eigen::MatrixXd buildRefDDC(std::string fileName);

// Shape-matching for a target HC
int shapeMatchHC(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 const Eigen::MatrixXd &refPoints, cage::Cage cageUnit,
                 std::vector<std::vector<int>> rings,
                 std::vector<std::vector<int>> nList, std::vector<double> *quat,
                 double *rmsd);

// Shape-matching for a target DDC
int shapeMatchDDC(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                  const Eigen::MatrixXd &refPoints,
                  std::vector<cage::Cage> cageList, int cageIndex,
                  std::vector<std::vector<int>> rings,
                  std::vector<double> *quat, double *rmsd);

// Calulate the RMSD for each ring, using RMSD values (rmsd) obtained from the
// shape-matching of each cage
int updateRMSDatom(std::vector<std::vector<int>> rings, cage::Cage cageUnit,
                   double rmsd, std::vector<double> *rmsdPerAtom,
                   std::vector<int> *noOfCommonAtoms,
                   std::vector<cage::iceType> atomTypes);

// Average the RMSD per atom
int averageRMSDatom(std::vector<double> *rmsdPerAtom,
                    std::vector<int> *noOfCommonAtoms);

// Topological network methods
// Finds the HCs and DDCs for the system
std::vector<cage::Cage> topoBulkCriteria(
    std::string path, std::vector<std::vector<int>> rings,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int firstFrame,
    int *numHC, int *numDDC, std::vector<ring::strucType> *ringType);

// Clustering
// Clusters cages using the Stillinger algorithm and prints out individual XYZ
// files of clusters.
int clusterCages(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 std::string path, std::vector<std::vector<int>> rings,
                 std::vector<cage::Cage> cageList, int numHC, int numDDC);

// Gets the atoms in the cages of a given cluster
std::vector<int> atomsFromCages(std::vector<std::vector<int>> rings,
                                std::vector<cage::Cage> cageList,
                                std::vector<int> clusterCages);

}  // namespace tum3

#endif  // __BULKTUM_H_
