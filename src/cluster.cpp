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

#include <cluster.hpp>
#include <iostream>

namespace bg = boost::geometry;

/**
 * @details Finds the number of particles in the largest ice cluster, for a
 * given frame, using Stoddard's clustering algorithm (Stoddard J. Comp. Phys.,
 * 27, 291,
 * 1977)](https://www.sciencedirect.com/science/article/pii/0021999178900116)
 * @param[in] path Output directory path, specified by the user
 * @param[in] yCloud The molSys::PointCloud for all particles, regardless of
 *  type classification
 * @param[in, out] iceCloud The molSys::PointCloud for the largest ice cluster
 *  of the ice-like molecules
 * @param[in] nList Row-ordered neighbour list by atom ID
 * @param[in] isIce Holds a bool value for each particle in yCloud. This is
 *  true for ice-like molecules and false otherwise
 * @param[in] list Linked list created by the clustering algorithm.
 * @param[in] nClusters Contains the number of particles in every ice-like
 *  cluster found
 * @param[in] indexNumber Unordered map for mapping the cluster ID indices to
 * the number in each cluster
 * @param[in] firstFrame First frame to be analyzed
 */
int clump::largestIceCluster(
    std::string path, molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
    std::vector<std::vector<int>> nList, std::vector<bool> *isIce,
    std::vector<int> *list, std::vector<int> *nClusters,
    std::unordered_map<int, int> *indexNumber, int firstFrame) {
  //
  int kAtomID;                    // Atom ID of the nearest neighbour
  int iClusterNumber;             // Number in the current cluster
  int nLargestCluster;            // Number in the largest cluster
  std::vector<int> linkedList;    // Linked list
  int j;                          // Index
  std::vector<int> startingIndex; // Contains the vector index in list
                                  // corresponding to a particular cluster
  int temp;                       // For swapping
  std::vector<bool>
      visited; // To make sure you don't go through the same atoms again.
  molSys::Point<double> iPoint; // Current point
  int currentIndex;             // Current index

  // -----------------------------------------------------------
  // INITIALIZATION
  linkedList.resize(yCloud->nop, -1); // init to dummy value
  // Initial values of the list. -1 is a dummy value if the molecule is
  // water or not in the slice
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    // Skip if the molecule is water or if it is not in the slice
    if ((*isIce)[iatom] == false || yCloud->pts[iatom].inSlice == false) {
      continue;
    } // skip for water or not in slice
    // Otherwise, assign the index as the ID
    linkedList[iatom] = iatom;
  } // init of cluster IDs
  // -----------------------------------------------------------
  // Get the linked list
  for (int i = 0; i < yCloud->nop - 1; i++) {
    //
    // Skip if the molecule is water or if it is not in the slice
    if (linkedList[i] == -1) {
      continue;
    } // skip for water or not in slice
    //
    // If iatom is already in a cluster, skip it
    if (linkedList[i] != i) {
      continue;
    } // Already part of a cluster
    //
    j = i; // Init of j
    // Execute the next part of the loop while j is not equal to i
    do {
      //
      // Go through the rest of the atoms (KLOOP)
      for (int k = i + 1; k < yCloud->nop; k++) {
        // Skip if not ice
        if ((*isIce)[k] == false) {
          continue;
        } // not ice
        // Skip if already part of a cluster
        if (linkedList[k] != k) {
          continue;
        } // Already part of a cluster
        //
        // Check to see if k is a nearest neighbour of j
        kAtomID = yCloud->pts[k].atomID; // Atom ID
        auto it = std::find(nList[j].begin() + 1, nList[j].end(), kAtomID);
        if (it != nList[j].end()) {
          // Swap!
          temp = linkedList[j];
          linkedList[j] = linkedList[k];
          linkedList[k] = temp;
        } // j and k are nearest neighbours
      }   // end of loop through k (KLOOP)
      //
      j = linkedList[j];
    } while (j != i); // end of control for j!=i
    //

  } // end of looping through every i

  // -----------------------------------------------------------
  // Get the starting index (which is the vector index in list) corresponding to
  // clusters with more than one molecule in them
  int nextElement; // value in list
  int index;       // starting index value
  // init
  visited.resize(yCloud->nop);

  for (int i = 0; i < yCloud->nop; i++) {
    //
    if (visited[i]) {
      continue;
    }                  // Already counted
    visited[i] = true; // Visited
    if (linkedList[i] == -1) {
      continue;
    } // not ice
    if (linkedList[i] == i) {
      continue;
    } // only one element
    //
    currentIndex = i;
    nextElement = linkedList[currentIndex];
    index = i;
    iClusterNumber = 1; // at least one value
    while (nextElement != index) {
      iClusterNumber++;
      currentIndex = nextElement;
      visited[currentIndex] = true;
      nextElement = linkedList[currentIndex];
    } // get number
    // Update startingIndex with index
    startingIndex.push_back(index);
    // Update the number of molecules in the cluster
    (*nClusters).push_back(iClusterNumber);
  } // end of loop through
  // -----------------------------------------------------------
  // Get the largest cluster and save the atoms to the iceCloud pointCloud
  nLargestCluster = *std::max_element((*nClusters).begin(), (*nClusters).end());
  int lClusIndex = distance(
      (*nClusters).begin(),
      max_element(
          (*nClusters).begin(),
          (*nClusters).end())); // index of the cluster with the largest number
  // -----------------------------------------------------------
  // Loop through the linked list from the starting element
  // startingIndex[lClusIndex].
  int startCluster =
      startingIndex[lClusIndex]; // index in linkedList to start from

  // L[i]->j->L[j]->k->L[k]->i
  index = startCluster;
  // Update iceCloud
  iPoint = yCloud->pts[index];
  iceCloud->pts.push_back(iPoint);
  //
  nextElement = linkedList[startCluster];
  while (nextElement != index) {
    currentIndex = nextElement;
    // update with currentIndex
    iPoint = yCloud->pts[currentIndex];
    iceCloud->pts.push_back(iPoint);
    // end update
    nextElement = linkedList[currentIndex];
  } // get the largest cluster
  // -----------------------------------------------------------
  // Update other variables in iceCloud

  iceCloud->currentFrame = yCloud->currentFrame;
  iceCloud->nop = iceCloud->pts.size();
  iceCloud->box = yCloud->box;
  iceCloud->boxLow = yCloud->boxLow;

  // Update idIndexMap
  for (int iatom = 0; iatom < iceCloud->nop; iceCloud++) {
    iceCloud->idIndexMap[iceCloud->pts[iatom].atomID] = iatom;
  } // end of loop through iceCloud

  // -----------------------------------------------------------
  // Write out the cluster statistics
  int totalClusters = (*nClusters).size(); // Total number of clusters
  int smallestCluster = nLargestCluster =
      *std::min_element((*nClusters).begin(),
                        (*nClusters).end()); // Size of the smallest cluster
  double avgClusterSize = 0.0;

  // Get the average cluster size
  for (int iCluster = 0; iCluster < totalClusters; iCluster++) {
    avgClusterSize += (*nClusters)[iCluster];
  } // Loop through the clusters
  // Normalize by the number
  if (totalClusters == 0) {
    avgClusterSize = 0;
  } else {
    avgClusterSize /= totalClusters;
  }

  iceCloud->currentFrame = yCloud->currentFrame;
  nLargestCluster = (*nClusters)[lClusIndex];

  // Write out to the file
  sout::writeClusterStats(path, yCloud->currentFrame, nLargestCluster,
                          totalClusters, smallestCluster, avgClusterSize,
                          firstFrame);

  // -----------------------------------------------------------
  return 0;
}

/**
 * @details Get the linked list of a cluster, given by iceCloud,
 *  for a single cluster. Required for cluster re-centering
 * @param[in] iceCloud The molSys::PointCloud for the largest ice cluster of
 *  the ice-like molecules
 * @param[in] nList Row-ordered neighbour list by atom index, not ID
 * @param[in, out] linkedList Linked list created by the clustering algorithm
 *  for the largest ice cluster
 */
int clump::singleClusterLinkedList(
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
    std::vector<std::vector<int>> nList, std::vector<int> *linkedList) {
  //
  int j;
  int temp; // For swapping indices
  //
  // -----------------------------------------------------------
  // INITIALIZATION
  (*linkedList).resize(iceCloud->nop);
  // Initial values of the list.
  for (int iatom = 0; iatom < iceCloud->nop; iatom++) {
    // Assign the index as the ID
    (*linkedList)[iatom] = iatom;
  } // init of cluster IDs
  // -----------------------------------------------------------
  // Get the linked list
  for (int i = 0; i < iceCloud->nop - 1; i++) {
    //
    // If iatom is already in a cluster, skip it
    if ((*linkedList)[i] != i) {
      continue;
    } // Already part of a cluster
    //
    j = i; // Init of j
    // Execute the next part of the loop while j is not equal to i
    do {
      //
      // Go through the rest of the atoms (KLOOP)
      for (int k = i + 1; k < iceCloud->nop; k++) {
        // Skip if already part of a cluster
        if ((*linkedList)[k] != k) {
          continue;
        } // Already part of a cluster
        //
        // Check to see if k is a nearest neighbour of j
        auto it = std::find(nList[j].begin() + 1, nList[j].end(), k);
        if (it != nList[j].end()) {
          // Swap!
          temp = (*linkedList)[j];
          (*linkedList)[j] = (*linkedList)[k];
          (*linkedList)[k] = temp;
        } // j and k are nearest neighbours
      }   // end of loop through k (KLOOP)
      //
      j = (*linkedList)[j];
    } while (j != i); // end of control for j!=i
    //

  } // end of looping through every i

  // Return 0
  return 0;
} // end of function

/**
 * @details Does the cluster analysis of ice particles in the system. Returns a
 *  molSys::PointCloud of the largest ice cluster (using the @f$ q_6 @f$
 * parameter by default). Uses the full neighbour list (by ID) according to the
 * full PointCloud yCloud. Returns a neighbour list by index, according to the
 * largest ice cluster.
 * @param[in] path Output directory path, specified by the user
 * @param[in, out] iceCloud The molSys::PointCloud for the largest ice cluster
 *  of the ice-like molecules
 * @param[in] yCloud The molSys::PointCloud for all the particles in the frame,
 *  regardless of ice type
 * @param[in] nList Row-ordered neighbour list by atom ID
 * @param[in] iceNeighbourList Row-ordered neighbour list by atom index, not
 *  ID, according to the iceCloud atoms
 * @param[in] cutoff Cutoff for the nearest neighbours
 * @param[in] firstFrame First frame to be analyzed
 * @param[in] bopAnalysis This determines which method to use for determining
 *  the ice-like nature of the particles. This can be "q6" or "chill", for using
 *  the @f$ q_6 @f$ parameter or CHILL algorithm, respectively
 */
int clump::clusterAnalysis(
    std::string path,
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList,
    std::vector<std::vector<int>> &iceNeighbourList, double cutoff,
    int firstFrame, std::string bopAnalysis) {
  //
  std::vector<bool> isIce;    // For every particle in yCloud, has a value
  int nTotalIce;              // Total number of ice-like molecules
  std::vector<int> clusterID; // Vector containing the cluster IDs of all the
                              // ice-like molecules.
  std::vector<int> nClusters; // Number of clusters for each index
  std::unordered_map<int, int>
      indexNumber; // Map for cluster index and index in the number vector

  // -------------------------------------------------------
  // Init
  // Clear the largest ice cluster pointCloud.
  *iceCloud = molSys::clearPointCloud(iceCloud);
  // Clear the neighbour list by index
  nneigh::clearNeighbourList(iceNeighbourList);
  // Init the vector of bools for every particle in yCloud
  isIce.resize(yCloud->nop);
  nTotalIce = 0; // Total number of ice-like molecules
  // -------------------------------------------------------
  // Use a bond-orientational parameter to find ice-like particles
  // Q6
  if (bopAnalysis == "q6") {
    //
    std::vector<double> q6Values;
    // Calculate the Q6 parameters for every point in yCloud
    q6Values = chill::getq6(yCloud, nList);
    // Assign values to isIce according to the values
    // of the q6 parameter. If q6 is greater than 0.5, it is ice-like.
    for (int iatom = 0; iatom < yCloud->nop; iatom++) {
      // If q6 is greater than 0.5, it is ice-like
      if (q6Values[iatom] > 0.5) {
        isIce[iatom] = true; // is ice-like; by default false
        nTotalIce++;         // Add to the number of ice-like molecules
      }                      // ice-like molecule found
    }                        // end of loop through every atom in yCloud
  } // end of getting a vector of bools for ice-like particles
  //
  // CHILL
  // Q6
  if (bopAnalysis == "chill") {
    //
    *yCloud = chill::getCorrel(yCloud, nList, false);
    // Get the ice types
    *yCloud = chill::getIceTypeNoPrint(yCloud, nList, false);
    // Assign values to isIce according to the CHILL algorithm
    for (int iatom = 0; iatom < yCloud->nop; iatom++) {
      // If it is an ice-like molecule, add it, otherwise skip
      if (yCloud->pts[iatom].iceType == molSys::atom_state_type::water) {
        continue;
      } // water
      if (yCloud->pts[iatom].iceType == molSys::atom_state_type::unclassified) {
        continue;
      }                    // unclassified
      isIce[iatom] = true; // is ice-like; by default false
      nTotalIce++;         // Add to the number of ice-like molecules

    } // end of loop through every atom in yCloud
  }   // end of getting a vector of bools for ice-like particles
  // -------------------------------------------------------
  // Get the largest ice cluster and other data
  clump::largestIceCluster(path, yCloud, iceCloud, nList, &isIce, &clusterID,
                           &nClusters, &indexNumber, firstFrame);

  // -------------------------------------------------------
  // Get the neighbour list by index according to the largest ice cluster
  // pointCloud
  iceNeighbourList = nneigh::getNewNeighbourListByIndex(iceCloud, cutoff);

  return 0;
} // end of function

/**
 * @details Recenters the largest ice cluster, by applying a transformation on
 * the largest ice cluster coordinates. Requires the neighbour list BY INDEX.
 * @param[in] iceCloud The molSys::PointCloud for the largest ice cluster of
 *  the ice-like molecules
 * @param[in] nList Row-ordered neighbour list by atom index, for the
 *  molSys::PointCloud iceCloud
 */
int clump::recenterClusterCloud(
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
    std::vector<std::vector<int>> nList) {
  //
  int dim = 3; // Dimensions
  std::vector<double> box = iceCloud->box;
  std::vector<double> boxLow = iceCloud->boxLow;
  std::vector<double> boxHigh;
  double xBoxCenter, yBoxCenter, zBoxCenter; // Centroid of the simulation box
  double x_centroid, y_centroid, z_centroid; // Centroid of the cluster
  // Variables for the linked list
  std::vector<int> linkedList; // Contains the linked list for the cluster
  std::vector<bool> visited; // Records whether an item has been visited or not

  // To avoid long confusing lines, fill boxHigh
  for (int k = 0; k < dim; k++) {
    boxHigh.push_back(boxLow[k] + box[k]);
  } // end of filling up boxHigh

  // --------------------------------------------------------------------------
  // Get the linked list of the cluster
  clump::singleClusterLinkedList(iceCloud, nList, &linkedList);
  // --------------------------------------------------------------------------
  // Loop through the entire looped list
  // init
  visited.resize(iceCloud->nop);

  // The starting value is the first atom
  int iatom = 0;    // Atom index of the 'starting value'
  int currentIndex; // Current atom
  int nextElement;  // Next linked atom
  int index;        // Keeps track of the first element in the linked list
  double x_ij, y_ij, z_ij; // Relative distance between the two atoms
  double xPBC, yPBC, zPBC; // Actual distance

  // Loop through the entire linked list
  for (int i = 0; i < iceCloud->nop; i++) {
    //
    if (visited[i]) {
      continue;
    }                  // Already counted
    visited[i] = true; // Visited
    // Should never execute
    if (linkedList[i] == i) {
      continue;
    } // only one element
    //
    currentIndex = i;
    nextElement = linkedList[currentIndex];
    index = i;
    // Keep looping
    while (nextElement != index) {
      currentIndex = nextElement;
      visited[currentIndex] = true;
      nextElement = linkedList[currentIndex];
      // -----------------------------------
      // Get the relative distance between the central atom (iatom)
      // and the next element
      // Coordinates
      // if (nextElement != index) {
      x_ij = iceCloud->pts[currentIndex].x - iceCloud->pts[nextElement].x;
      y_ij = iceCloud->pts[currentIndex].y - iceCloud->pts[nextElement].y;
      z_ij = iceCloud->pts[currentIndex].z - iceCloud->pts[nextElement].z;
      // Shift the nextElement if it's on the other side of the box
      // Shift x
      if (fabs(x_ij) > 0.5 * box[0]) {
        // Get the actual distance
        xPBC = box[0] - fabs(x_ij);
        if (x_ij < 0) {
          iceCloud->pts[nextElement].x = iceCloud->pts[currentIndex].x - xPBC;
        } // To the -x side of currentIndex
        else {
          iceCloud->pts[nextElement].x = iceCloud->pts[currentIndex].x + xPBC;
        } // Add to the + side
      }   // Shift nextElement
      //
      // Shift y
      if (fabs(y_ij) > 0.5 * box[1]) {
        // Get the actual distance
        yPBC = box[1] - fabs(y_ij);
        if (y_ij < 0) {
          iceCloud->pts[nextElement].y = iceCloud->pts[currentIndex].y - yPBC;
        } // To the -y side of currentIndex
        else {
          iceCloud->pts[nextElement].y = iceCloud->pts[currentIndex].y + yPBC;
        } // Add to the + side
      }   // Shift nextElement
      //
      // Shift z
      if (fabs(z_ij) > 0.5 * box[2]) {
        // Get the actual distance
        zPBC = box[2] - fabs(z_ij);
        if (z_ij < 0) {
          iceCloud->pts[nextElement].z = iceCloud->pts[currentIndex].z - zPBC;
        } // To the -z side of currentIndex
        else {
          iceCloud->pts[nextElement].z = iceCloud->pts[currentIndex].z + zPBC;
        } // Add to the + side
      }   // Shift nextElement
          // -----------------------------------
      // }  // don't shift the last atom!

    } // End of going through linked atoms
    //
  } // end of loop through atoms
  // --------------------------------------------------------------------------
  // Center of the simulation box
  xBoxCenter = 0.5 * (boxLow[0] + boxHigh[0]);
  yBoxCenter = 0.5 * (boxLow[1] + boxHigh[1]);
  zBoxCenter = 0.5 * (boxLow[2] + boxHigh[2]);
  // --------------------------------------------------------------------------
  // Get the centroid of the ice cluster.
  x_centroid = 0.0;
  y_centroid = 0.0;
  z_centroid = 0.0;

  for (int i = 0; i < iceCloud->nop; i++) {
    x_centroid += iceCloud->pts[i].x;
    y_centroid += iceCloud->pts[i].y;
    z_centroid += iceCloud->pts[i].z;
  } // end of loop through the particles

  if (iceCloud->nop == 0) {
    std::cerr << "There are no particles in the cluster.\n";
    return 1;
  } // error

  // Normalize by the number of particles.
  x_centroid /= iceCloud->nop;
  y_centroid /= iceCloud->nop;
  z_centroid /= iceCloud->nop;

  // --------------------------------------------------------------------------
  // Get the distance to shift by:
  double xShift = x_centroid - xBoxCenter;
  double yShift = y_centroid - yBoxCenter;
  double zShift = z_centroid - zBoxCenter;

  // Loop through all atoms and shift them
  for (int i = 0; i < iceCloud->nop; i++) {
    iceCloud->pts[i].x -= xShift;
    iceCloud->pts[i].y -= yShift;
    iceCloud->pts[i].z -= zShift;
  } // end of loop through the particles

  return 0;
} // end of function
