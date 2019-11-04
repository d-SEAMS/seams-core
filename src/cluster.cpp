#include <cluster.hpp>
#include <iostream>

namespace bg = boost::geometry;

/********************************************/ /**
                                                *  Finds the number of particles
                                                *in the largest ice cluster, for
                                                *a given frame.
                                                *  @param[in] iceCloud The input
                                                *molSys::PointCloud for all the
                                                *ice-like particles
                                                *  @param[in] cutoff The cut-off
                                                *distance for determining
                                                *nearest-neighbours. For water,
                                                *the cut-off value is typically
                                                *taken to be \f$3.2\f$ or
                                                *\f$3.5\f$ Angstrom,
                                                *encompassing the
                                                *first-neighbour shell
                                                *molecules.
                                                *  @param[in] printCluster
                                                *Decides whether the cluster
                                                *should be printed out to a file
                                                *as a lammpstrj or not
                                                *  @param[in] isSlice Decides
                                                *whether there is a slice (true)
                                                *or not (false) \return an int
                                                *value holding the number of
                                                *particles in the largest ice
                                                *cluster
                                                ***********************************************/
int clump::largestIceCluster(
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud, double cutoff,
    bool printCluster, bool isSlice) {
  std::vector<int>
      clusterFlag;  // This will contain flags for each solid atom. If
                    // 'flagged', move onto the next unflagged element
  std::vector<int> linkedListCluster;  // This will contain the ID of the
                                       // cluster it belongs to
  std::vector<int>
      clusterID;  // This will contain the starting point of a cluster
  std::vector<int> nCluster;  // No. of particles in each cluster
  int nLargestCluster = 0;    // No. of particles in the largest cluster
  int nnumNeighbours;         // Number of nearest neighbours
  int j;
  int noc = 0;      // no. of particles in a cluster
  double r_jk = 0;  // Distance between j^th and k^th atom

  clusterFlag.resize(iceCloud->nop);  // Init cluster flag vector to 0. If added
                                      // to a cluster, then skip
  linkedListCluster.reserve(iceCloud->nop);

  // Initialize the linked list
  for (int iatom = 0; iatom < iceCloud->nop; iatom++) {
    linkedListCluster[iatom] = iatom;
  }

  // -------------------------------------------------
  // Construct the linked list
  for (int i = 0; i < iceCloud->nop - 1; i++) {
    if (linkedListCluster[i] != i) {
      continue;
    }  // i is already in a cluster
    j = i;
    do {
      for (int k = i + 1; k < iceCloud->nop; k++) {
        // If L[k] == k and r_jk < cutoff then
        if (linkedListCluster[k] == k) {
          r_jk = gen::periodicDist(iceCloud, j, k);
          if (r_jk <= cutoff) {
            iter_swap(linkedListCluster.begin() + j,
                      linkedListCluster.begin() + k);
          }  // end of check that j&k are neighbours
        }    // check for L[k] and k
      }
      j = linkedListCluster[j];
    } while (j != i);

  }  // Loop through N-1 atoms
  // -------------------------------------------------

  // Get the number of particles per cluster
  int jatom, iatom;

  for (int i = 0; i < iceCloud->nop; i++) {
    if (clusterFlag[i] == 1) {
      continue;
    }  // The particle has already been tagged as part of a previous cluster
    // Start a new cluster
    clusterID.push_back(i);
    iatom = i;
    noc = 0;
    for (int j = 0; j < iceCloud->nop; j++) {
      jatom = linkedListCluster[iatom];
      clusterFlag[jatom] = 1;
      iatom = jatom;
      noc++;
      if (jatom == i) {
        break;
      }
    }  // The maximum no. of particles in a cluster is N
    nCluster.push_back(noc);
  }  // Loop through all solid particles

  auto largest = std::max_element(nCluster.begin(), nCluster.end());
  int indexOfCluster = std::distance(std::begin(nCluster), largest);
  nLargestCluster = nCluster[indexOfCluster];
  // -----------------------------------------------
  // If the user wants, print out the largest ice cluster to a lammpstrj

  if (printCluster) {
    // Header of traj file
    std::ofstream outputFile;
    // Create a new file in the output directory
    outputFile.open("largestCluster.lammpstrj", std::ios_base::app);
    outputFile << "ITEM: TIMESTEP\n";
    outputFile << iceCloud->currentFrame << "\n";
    outputFile << "ITEM: NUMBER OF ATOMS\n";
    outputFile << nLargestCluster << "\n";
    outputFile << "ITEM: BOX BOUNDS pp pp pp\n";
    for (int k = 0; k < iceCloud->boxLow.size(); k++) {
      outputFile << iceCloud->boxLow[k] << " "
                 << iceCloud->boxLow[k] + iceCloud->box[k];
      // for triclinic boxes
      if (iceCloud->box.size() == 2 * iceCloud->boxLow.size()) {
        // The tilt factors are saved after the box lengths; so add 3
        outputFile
            << " "
            << iceCloud->box[k + iceCloud->boxLow.size()];  // this would be +2
                                                            // for a 2D box
      }
      outputFile << "\n";
    }  // end of printing box lengths
    outputFile << "ITEM: ATOMS id mol type x y z\n";

    iatom = clusterID[indexOfCluster];
    for (int i = 0; i < nLargestCluster; i++) {
      jatom = linkedListCluster[iatom];  // Get index of the atom to print
      iatom = jatom;
      outputFile << iceCloud->pts[jatom].atomID << " "
                 << iceCloud->pts[jatom].molID << " "
                 << iceCloud->pts[jatom].iceType << " "
                 << iceCloud->pts[jatom].x << " " << iceCloud->pts[jatom].y
                 << " " << iceCloud->pts[jatom].z << "\n";
    }  // End of traversing through the largest ice cluster
  }

  // -------------------------------------------------

  return nLargestCluster;
}

/********************************************/ /**
 *  Does the cluster analysis of ice particles in the system. Returns a
 pointCloud of the largest ice cluster (using the q6 parameter by default).
 Uses the full neighbour list (by ID) according to the full
 PointCloud yCloud.
 ***********************************************/
molSys::PointCloud<molSys::Point<double>, double> clump::clusterAnalysis(
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, std::string bopAnalysis) {
  //
  std::vector<bool> isIce;  // For every particle in yCloud, has a value
  int nTotalIce;            // Total number of ice-like molecules

  // -------------------------------------------------------
  // Init
  // Clear the largest ice cluster pointCloud.
  *iceCloud = molSys::clearPointCloud(iceCloud);
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
  // -------------------------------------------------------

  return *iceCloud;
}  // end of function
