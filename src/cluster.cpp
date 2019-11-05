#include <cluster.hpp>
#include <iostream>

namespace bg = boost::geometry;

/********************************************/ /**
                                                *  Finds the number of particles
                                                *in the largest ice cluster, for
                                                *a given frame, using Stoddard's
                                                *1978 algorithm
                                                ***********************************************/
int clump::largestIceCluster(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
    std::vector<std::vector<int>> nList, std::vector<bool> *isIce,
    std::vector<int> *list, std::vector<int> *nClusters,
    std::unordered_map<int, int> *indexNumber) {
  //
  int kAtomID;                     // Atom ID of the nearest neighbour
  int iClusterNumber;              // Number in the current cluster
  int nLargestCluster;             // Number in the largest cluster
  std::vector<int> linkedList;     // Linked list
  int j;                           // Index
  std::vector<int> startingIndex;  // Contains the vector index in list
                                   // corresponding to a particular cluster
  int temp;                        // For swapping
  std::vector<bool>
      visited;  // To make sure you don't go through the same atoms again.
  molSys::Point<double> iPoint;  // Current point
  int currentIndex;              // Current index

  // -----------------------------------------------------------
  // INITIALIZATION
  linkedList.resize(yCloud->nop, -1);  // init to dummy value
  // Initial values of the list. -1 is a dummy value if the molecule is
  // water or not in the slice
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    // Skip if the molecule is water or if it is not in the slice
    if ((*isIce)[iatom] == false || yCloud->pts[iatom].inSlice == false) {
      continue;
    }  // skip for water or not in slice
    // Otherwise, assign the index as the ID
    linkedList[iatom] = iatom;
  }  // init of cluster IDs
  // -----------------------------------------------------------
  // Get the linked list
  for (int i = 0; i < yCloud->nop - 1; i++) {
    //
    // Skip if the molecule is water or if it is not in the slice
    if (linkedList[i] == -1) {
      continue;
    }  // skip for water or not in slice
    //
    // If iatom is already in a cluster, skip it
    if (linkedList[i] != i) {
      continue;
    }  // Already part of a cluster
    //
    j = i;  // Init of j
    // Execute the next part of the loop while j is not equal to i
    do {
      //
      // Go through the rest of the atoms (KLOOP)
      for (int k = i + 1; k < yCloud->nop; k++) {
        // Skip if not ice
        if ((*isIce)[k] == false) {
          continue;
        }  // not ice
        // Skip if already part of a cluster
        if (linkedList[k] != k) {
          continue;
        }  // Already part of a cluster
        //
        // Check to see if k is a nearest neighbour of j
        kAtomID = yCloud->pts[k].atomID;  // Atom ID
        auto it = std::find(nList[j].begin() + 1, nList[j].end(), kAtomID);
        if (it != nList[j].end()) {
          // Swap!
          temp = linkedList[j];
          linkedList[j] = linkedList[k];
          linkedList[k] = temp;
        }  // j and k are nearest neighbours
      }    // end of loop through k (KLOOP)
      //
      j = linkedList[j];
    } while (j != i);  // end of control for j!=i
    //

  }  // end of looping through every i

  // -----------------------------------------------------------
  // Get the starting index (which is the vector index in list) corresponding to
  // clusters with more than one molecule in them
  int nextElement;  // value in list
  int index;        // starting index value
  // init
  visited.resize(yCloud->nop);

  for (int i = 0; i < yCloud->nop; i++) {
    //
    if (visited[i]) {
      continue;
    }                   // Already counted
    visited[i] = true;  // Visited
    if (linkedList[i] == -1) {
      continue;
    }  // not ice
    if (linkedList[i] == i) {
      continue;
    }  // only one element
    //
    currentIndex = i;
    nextElement = linkedList[currentIndex];
    index = i;
    iClusterNumber = 1;  // at least one value
    while (nextElement != index) {
      iClusterNumber++;
      currentIndex = nextElement;
      visited[currentIndex] = true;
      nextElement = linkedList[currentIndex];
    }  // get number
    // Update startingIndex with index
    startingIndex.push_back(index);
    // Update the number of molecules in the cluster
    (*nClusters).push_back(iClusterNumber);
  }  // end of loop through
  // -----------------------------------------------------------
  // Get the largest cluster and save the atoms to the iceCloud pointCloud
  nLargestCluster = *std::max_element((*nClusters).begin(), (*nClusters).end());
  int lClusIndex = distance(
      (*nClusters).begin(),
      max_element(
          (*nClusters).begin(),
          (*nClusters).end()));  // index of the cluster with the largest number
  // -----------------------------------------------------------
  // Loop through the linked list from the starting element
  // startingIndex[lClusIndex].
  int startCluster =
      startingIndex[lClusIndex];  // index in linkedList to start from

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
  }  // get the largest cluster
  // -----------------------------------------------------------
  // Update other variables in iceCloud

  iceCloud->currentFrame = yCloud->currentFrame;
  iceCloud->nop = iceCloud->pts.size();
  iceCloud->box = yCloud->box;
  iceCloud->boxLow = yCloud->boxLow;

  // Update idIndexMap
  for (int iatom = 0; iatom < iceCloud->nop; iceCloud++) {
    iceCloud->idIndexMap[iceCloud->pts[iatom].atomID] = iatom;
  }  // end of loop through iceCloud

  return 0;
}

/********************************************/ /**
 *  Does the cluster analysis of ice particles in the system. Returns a
 pointCloud of the largest ice cluster (using the q6 parameter by default).
 Uses the full neighbour list (by ID) according to the full
 PointCloud yCloud. Returns a neighbour list by index, according to the largest
 ice cluster.
 ***********************************************/
int clump::clusterAnalysis(
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList,
    std::vector<std::vector<int>> &iceNeighbourList, double cutoff,
    std::string bopAnalysis) {
  //
  std::vector<bool> isIce;     // For every particle in yCloud, has a value
  int nTotalIce;               // Total number of ice-like molecules
  std::vector<int> clusterID;  // Vector containing the cluster IDs of all the
                               // ice-like molecules.
  std::vector<int> nClusters;  // Number of clusters for each index
  std::unordered_map<int, int>
      indexNumber;  // Map for cluster index and index in the number vector

  // -------------------------------------------------------
  // Init
  // Clear the largest ice cluster pointCloud.
  *iceCloud = molSys::clearPointCloud(iceCloud);
  // Clear the neighbour list by index
  nneigh::clearNeighbourList(iceNeighbourList);
  // Init the vector of bools for every particle in yCloud
  isIce.resize(yCloud->nop);
  nTotalIce = 0;  // Total number of ice-like molecules
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
        isIce[iatom] = true;  // is ice-like; by default false
        nTotalIce++;          // Add to the number of ice-like molecules
      }                       // ice-like molecule found
    }                         // end of loop through every atom in yCloud
  }  // end of getting a vector of bools for ice-like particles
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
      if (yCloud->pts[iatom].iceType == molSys::water) {
        continue;
      }  // water
      if (yCloud->pts[iatom].iceType == molSys::unclassified) {
        continue;
      }                     // unclassified
      isIce[iatom] = true;  // is ice-like; by default false
      nTotalIce++;          // Add to the number of ice-like molecules

    }  // end of loop through every atom in yCloud
  }    // end of getting a vector of bools for ice-like particles
  // -------------------------------------------------------
  // Get the largest ice cluster and other data
  clump::largestIceCluster(yCloud, iceCloud, nList, &isIce, &clusterID,
                           &nClusters, &indexNumber);

  // -------------------------------------------------------
  // Get the neighbour list by index according to the largest ice cluster
  // pointCloud
  iceNeighbourList = nneigh::getNewNeighbourListByIndex(iceCloud, cutoff);

  return 0;
}  // end of function

/********************************************/ /**
 *  Recenters the largest ice cluster, by applying a transformation on the
 largest ice cluster coordinates.
 ***********************************************/
int clump::recenterClusterCloud(
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud) {
  //
  int dim = 3;       // Dimensions
  bool nearBoxHigh;  // iatom is close to the higher box boundary
  bool nearBoxLow;   // iatom is close to the lower box boundary
  int iatom = 0;     // Atom index
  double dr_ij;      // Relative unwrapped distance
  std::vector<double> distFromBoxHi;  // Distance from the box boundaries
  std::vector<double> distFromBoxLo;
  std::vector<double> box = iceCloud->box;
  std::vector<double> boxLow = iceCloud->boxLow;
  std::vector<double> boxHigh;
  double xBoxCenter, yBoxCenter, zBoxCenter;  // Centroid of the simulation box
  double x_centroid, y_centroid, z_centroid;  // Centroid of the cluster

  // boxHigh filling
  boxHigh.push_back(box[0] + boxLow[0]);
  boxHigh.push_back(box[1] + boxLow[1]);
  boxHigh.push_back(box[2] + boxLow[2]);

  // Init
  // x coordinates
  distFromBoxLo.push_back(fabs(iceCloud->pts[iatom].x - boxLow[0]));
  distFromBoxHi.push_back(fabs(iceCloud->pts[iatom].x - boxHigh[0]));
  // y coordinates
  distFromBoxLo.push_back(fabs(iceCloud->pts[iatom].y - boxLow[1]));
  distFromBoxHi.push_back(fabs(iceCloud->pts[iatom].y - boxHigh[1]));
  // z coordinates
  distFromBoxLo.push_back(fabs(iceCloud->pts[iatom].z - boxLow[2]));
  distFromBoxHi.push_back(fabs(iceCloud->pts[iatom].z - boxHigh[2]));

  // --------------------------------------------------------------------------
  // Loop through every dimension (k=0 or x to k=2 or z)
  for (int k = 0; k < dim; k++) {
    // init
    nearBoxHigh = false;
    nearBoxLow = false;
    //
    // Figure out whether iatom is close to the higher or lower box boundary
    if (distFromBoxHi[k] <= distFromBoxLo[k]) {
      nearBoxHigh = true;
    }  // close to the rHi boundary
    else {
      nearBoxLow = true;
    }  // close to rLo boundary

    // Get atom pairs
    for (int jatom = 1; jatom < iceCloud->nop; jatom++) {
      // Get the actual relative distance without PBCs
      dr_ij = gen::relativeDist(iceCloud, iatom, jatom, k);
      // If dr_ij is greater than half the box length, then shift it,
      if (dr_ij > 0.5 * box[k]) {
        dr_ij -= iceCloud->box[k] * round(dr_ij / iceCloud->box[k]);
        dr_ij = fabs(dr_ij);
        //
        if (k == 0) {
          if (nearBoxHigh) {
            iceCloud->pts[jatom].x = iceCloud->pts[iatom].x + dr_ij;
          }  // xHi
          else {
            iceCloud->pts[jatom].x = iceCloud->pts[iatom].x - dr_ij;
          }  // xLo
        }    // x
        else if (k == 1) {
          if (nearBoxHigh) {
            iceCloud->pts[jatom].y = iceCloud->pts[iatom].y + dr_ij;
          }  // yHi
          else {
            iceCloud->pts[jatom].y = iceCloud->pts[iatom].y - dr_ij;
          }  // yLo
        }    // y
        else {
          if (nearBoxHigh) {
            iceCloud->pts[jatom].z = iceCloud->pts[iatom].z + dr_ij;
          }  // zHi
          else {
            iceCloud->pts[jatom].z = iceCloud->pts[iatom].z - dr_ij;
          }  // zLo
        }    // z
      }      // end of shift
    }        // end loop through atom pairs
  }          // end of loop through every dimension

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
  }  // end of loop through the particles

  if (iceCloud->nop == 0) {
    std::cerr << "There are no particles in the cluster.\n";
    return 1;
  }  // error

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
  }  // end of loop through the particles

  return 0;
}  // end of function
