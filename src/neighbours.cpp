#include <iostream>
#include <math.h>
#include <neighbours.hpp>

/********************************************/ /**
 *  Function for building neighbour lists for each
 particle. Inefficient O(n^2) implementation
 Full neighbour list. Use this when iatom and jatom
 are different
 ***********************************************/
molSys::PointCloud<molSys::Point<double>, double>
nneigh::neighList(double rcutoff,
                  molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                  int typeI, int typeJ) {
  nneigh::NeighbourList nList;
  nneigh::Jatom tempJatom;
  double r_ij;
  int indexYay;
  std::vector<int> tempListIatom;

  // Resize the neighbour list
  nList.iVector.resize(yCloud->nop);

  // Loop through every iatom and find nearest neighbours within rcutoff
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    if (yCloud->pts[iatom].type != typeI) {
      continue;
    }
    *yCloud =
        nneigh::clearNeighList(yCloud, iatom); // Clear neighbour list if full
    for (int jatom = 0; jatom < yCloud->nop; jatom++) {
      if (yCloud->pts[jatom].type != typeJ) {
        continue;
      }
      // If the distance is greater than rcutoff, continue
      r_ij = gen::periodicDist(yCloud, iatom, jatom);
      if (r_ij > rcutoff) {
        continue;
      }
      tempJatom.r = r_ij;
      tempJatom.index = jatom;

      // Increase number
      nList.iVector[iatom].nearest_neighbours += 1;
      nList.iVector[jatom].nearest_neighbours += 1;
      // Update the neighbour indices
      nList.iVector[iatom].n.push_back(tempJatom);
      tempJatom.index = iatom;
      nList.iVector[jatom].n.push_back(tempJatom);

    } // End of loop through jatom
  }   // End of loop for iatom

  // Now sort in ascending order
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    if (yCloud->pts[iatom].type != typeI) {
      continue;
    } // check atom type
    // clear nlist if full TODO
    std::sort(nList.iVector[iatom].n.begin(), nList.iVector[iatom].n.end(),
              nneigh::compareByLength);

    // Update neighbour lists in yCloud
    for (int t = 0; t < nList.iVector[iatom].nearest_neighbours; t++) {
      indexYay = nList.iVector[iatom].n[t].index;
      yCloud->pts[iatom].neighList.push_back(indexYay);
    }
  }

  return *yCloud;
}

/********************************************/ /**
 *  Function for building neighbour lists for each
 particle. Inefficient O(n^2) implementation. This will only work with one type of atom
 ***********************************************/
molSys::PointCloud<molSys::Point<double>, double>
nneigh::neighListO(double rcutoff,
                   molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                   int typeI) {
  nneigh::NeighbourList nList;
  nneigh::Jatom tempJatom;
  double r_ij;
  int indexYay;
  std::vector<int> tempListIatom;

  // Resize the neighbour list
  nList.iVector.resize(yCloud->nop);

  // Loop through every iatom and find nearest neighbours within rcutoff
  for (int iatom = 0; iatom < yCloud->nop - 1; iatom++) {
    if (yCloud->pts[iatom].type != typeI) {
      continue;
    }
    *yCloud =
        nneigh::clearNeighList(yCloud, iatom); // Clear neighbour list if full
    for (int jatom = iatom + 1; jatom < yCloud->nop; jatom++) {
      if (yCloud->pts[jatom].type != typeI) {
        continue;
      }
      // If the distance is greater than rcutoff, continue
      r_ij = gen::periodicDist(yCloud, iatom, jatom);
      if (r_ij > rcutoff) {
        continue;
      }
      tempJatom.r = r_ij;
      tempJatom.index = jatom;

      // Increase number
      nList.iVector[iatom].nearest_neighbours += 1;
      nList.iVector[jatom].nearest_neighbours += 1;
      // Update the neighbour indices
      nList.iVector[iatom].n.push_back(tempJatom);
      tempJatom.index = iatom;
      nList.iVector[jatom].n.push_back(tempJatom);

    } // End of loop through jatom
  }   // End of loop for iatom

  // Now sort in ascending order
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    if (yCloud->pts[iatom].type != typeI) {
      continue;
    } // check atom type
    // clear nlist if full TODO
    std::sort(nList.iVector[iatom].n.begin(), nList.iVector[iatom].n.end(),
              nneigh::compareByLength);

    // Update neighbour lists in yCloud
    for (int t = 0; t < nList.iVector[iatom].nearest_neighbours; t++) {
      indexYay = nList.iVector[iatom].n[t].index;
      yCloud->pts[iatom].neighList.push_back(indexYay);
    }
  }

  return *yCloud;
}

/********************************************/ /**
 *  Function for building neighbour lists for each
 particle. Inefficient O(n^2) implementation. This will only work with one type of atom
 Half neighbour list
 ***********************************************/
molSys::PointCloud<molSys::Point<double>, double>
nneigh::halfNeighList(double rcutoff,
                      molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                      int typeI) {
  nneigh::NeighbourList nList;
  nneigh::Jatom tempJatom;
  double r_ij;
  int indexYay;
  std::vector<int> tempListIatom;

  // Resize the neighbour list
  nList.iVector.resize(yCloud->nop);

  // Loop through every iatom and find nearest neighbours within rcutoff
  for (int iatom = 0; iatom < yCloud->nop - 1; iatom++) {
    if (yCloud->pts[iatom].type != typeI) {
      continue;
    }
    *yCloud =
        nneigh::clearNeighList(yCloud, iatom); // Clear neighbour list if full
    for (int jatom = iatom + 1; jatom < yCloud->nop; jatom++) {
      if (yCloud->pts[jatom].type != typeI) {
        continue;
      }
      // If the distance is greater than rcutoff, continue
      r_ij = gen::periodicDist(yCloud, iatom, jatom);
      if (r_ij > rcutoff) {
        continue;
      }
      tempJatom.r = r_ij;
      tempJatom.index = jatom;

      // Increase number
      nList.iVector[iatom].nearest_neighbours += 1;

      // Update the neighbour indices
      nList.iVector[iatom].n.push_back(tempJatom);

    } // End of loop through jatom
  }   // End of loop for iatom

  // // Now sort in ascending order
  // for(int iatom=0; iatom<yCloud->nop; iatom++){
  // 	if(yCloud->pts[iatom].type!=typeI){continue;} // check atom type
  // 	// clear nlist if full TODO
  // 	std::sort(nList.iVector[iatom].n.begin(), nList.iVector[iatom].n.end(), nneigh::compareByLength);

  // 	// Update neighbour lists in yCloud
  // 	for(int t=0; t<nList.iVector[iatom].nearest_neighbours; t++){
  // 		indexYay = nList.iVector[iatom].n[t].index;
  // 		yCloud->pts[iatom].neighList.push_back(indexYay);
  // 	}
  // }

  // No need to sort!
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    if (yCloud->pts[iatom].type != typeI) {
      continue;
    } // check atom type
    // clear nlist if full TODO
    // std::sort(nList.iVector[iatom].n.begin(), nList.iVector[iatom].n.end(), nneigh::compareByLength);

    // Update neighbour lists in yCloud
    for (int t = 0; t < nList.iVector[iatom].nearest_neighbours; t++) {
      indexYay = nList.iVector[iatom].n[t].index;
      yCloud->pts[iatom].neighList.push_back(indexYay);
    }
  }

  return *yCloud;
}

// Clear neighbour list for the i^th atom if it is full TODO: Gives segfault.
molSys::PointCloud<molSys::Point<double>, double> nneigh::clearNeighList(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatom) {

  if (yCloud->pts[iatom].neighList.size() != 0) {
    yCloud->pts[iatom].neighList.resize(0);
    yCloud->pts[iatom].neighList.shrink_to_fit();
  }

  return *yCloud;
}
