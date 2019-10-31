#include <math.h>
#include <iostream>
#include <neighbours.hpp>

/********************************************/ /**
 *  Function for building neighbour lists for each
 particle. Inefficient O(n^2) implementation
 Full neighbour list. Use this when iatom and jatom
 are different
 ***********************************************/
std::vector<std::vector<int>> nneigh::neighList(
    double rcutoff, molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    int typeI, int typeJ) {
  std::vector<std::vector<int>> nList;  // Vector of vector of ints
  int jatomIndex;                       // Atom ID corresponding to jatom
  int iatomIndex;                       // Atom ID corresponding to iatom
  double r_ij;                          // cutoff

  // Initialize with nop (irrespective of type)
  // Initialize and fill the first element with the current atom ID whose
  // neighbour list will be filled
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    // Find the atom ID (key) given the index or iatom (value)
    auto itr = std::find_if(
        yCloud->idIndexMap.begin(), yCloud->idIndexMap.end(),
        [&iatom](const std::pair<int, int> &p) { return p.second == iatom; });
    // If found:
    if (itr == yCloud->idIndexMap.end()) {
      std::cerr << "Something is wrong with your idIndexMap!\n";
      continue;
    } else {
      iatomIndex = itr->first;
    }  // End of finding the atom ID to fill as the first element in the
       // neighbour list
    nList.push_back(std::vector<int>());  // Empty vector for the index iatom
    // Fill the first element with the atom ID of iatom itself
    nList[iatom].push_back(iatomIndex);
  }  // end of init

  // pairs of atoms of type I and J
  // Loop through every iatom and find nearest neighbours within rcutoff
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    if (yCloud->pts[iatom].type != typeI) {
      continue;
    }
    // Loop through the other atoms
    for (int jatom = 0; jatom < yCloud->nop; jatom++) {
      if (yCloud->pts[jatom].type != typeJ) {
        continue;
      }
      // If the distance is greater than rcutoff, continue
      r_ij = gen::periodicDist(yCloud, iatom, jatom);
      if (r_ij > rcutoff) {
        continue;
      }

      // Get the atom IDs for iatom and jatom
      auto gotI = std::find_if(
          yCloud->idIndexMap.begin(), yCloud->idIndexMap.end(),
          [&iatom](const std::pair<int, int> &p) { return p.second == iatom; });
      if (gotI == yCloud->idIndexMap.end()) {
        std::cerr << "Something is wrong with your idIndexMap!\n";
        return nList;
      } else {
        iatomIndex = gotI->first;
      }  // End of finding the atom ID for iatom
      // Find the atom ID of jatom
      auto gotJ = std::find_if(
          yCloud->idIndexMap.begin(), yCloud->idIndexMap.end(),
          [&jatom](const std::pair<int, int> &p) { return p.second == jatom; });
      if (gotJ == yCloud->idIndexMap.end()) {
        std::cerr << "Something is wrong with your idIndexMap!\n";
        return nList;
      } else {
        jatomIndex = gotJ->first;
      }  // End of finding the atom ID for jatom
      // Update the neighbour indices with atom IDs for iatom and jatom both
      // (full list)
      nList[iatom].push_back(jatomIndex);
      nList[jatom].push_back(iatomIndex);

    }  // End of loop through jatom
  }    // End of loop for iatom

  return nList;
}

/********************************************/ /**
 *  Function for building neighbour lists for each
 particle. Inefficient O(n^2) implementation. This will only work with one type
 of atom
 ***********************************************/
std::vector<std::vector<int>> nneigh::neighListO(
    double rcutoff, molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    int typeI) {
  std::vector<std::vector<int>>
      nList;       // Vector of vectors of the neighbour list
  double r_ij;     // Distance between iatom and jatom
  int iatomIndex;  // Atomic ID of the atom with index iatom
  int jatomIndex;  // Atomic ID of the atom with index jatom
  int indexYay;
  std::vector<int> tempListIatom;

  // Initialize and fill the first element with the current atom ID whose
  // neighbour list will be filled
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    // Find the atom ID (key) given the index or iatom (value)
    auto itr = std::find_if(
        yCloud->idIndexMap.begin(), yCloud->idIndexMap.end(),
        [&iatom](const std::pair<int, int> &p) { return p.second == iatom; });
    // If found:
    if (itr == yCloud->idIndexMap.end()) {
      std::cerr << "Something is wrong with your idIndexMap!\n";
      continue;
    } else {
      iatomIndex = itr->first;
    }  // End of finding the atom ID to fill as the first element in the
       // neighbour list
    nList.push_back(std::vector<int>());  // Empty vector for the index iatom
    // Fill the first element with the atom ID of iatom itself
    nList[iatom].push_back(iatomIndex);
  }  // end of init

  // Loop through every iatom and find nearest neighbours within rcutoff
  for (int iatom = 0; iatom < yCloud->nop - 1; iatom++) {
    if (yCloud->pts[iatom].type != typeI) {
      continue;
    }
    // Loop through the other atoms
    for (int jatom = iatom + 1; jatom < yCloud->nop; jatom++) {
      if (yCloud->pts[jatom].type != typeI) {
        continue;
      }
      // If the distance is greater than rcutoff, continue
      r_ij = gen::periodicDist(yCloud, iatom, jatom);
      if (r_ij > rcutoff) {
        continue;
      }

      // Get the atom IDs for iatom and jatom
      auto gotI = std::find_if(
          yCloud->idIndexMap.begin(), yCloud->idIndexMap.end(),
          [&iatom](const std::pair<int, int> &p) { return p.second == iatom; });
      if (gotI == yCloud->idIndexMap.end()) {
        std::cerr << "Something is wrong with your idIndexMap!\n";
        return nList;
      } else {
        iatomIndex = gotI->first;
      }  // End of finding the atom ID for iatom
      // Find the atom ID of jatom
      auto gotJ = std::find_if(
          yCloud->idIndexMap.begin(), yCloud->idIndexMap.end(),
          [&jatom](const std::pair<int, int> &p) { return p.second == jatom; });
      if (gotJ == yCloud->idIndexMap.end()) {
        std::cerr << "Something is wrong with your idIndexMap!\n";
        return nList;
      } else {
        jatomIndex = gotJ->first;
      }  // End of finding the atom ID for jatom
      // Update the neighbour indices with atom IDs for iatom and jatom both
      // (full list)
      nList[iatom].push_back(jatomIndex);
      nList[jatom].push_back(iatomIndex);

    }  // End of loop through jatom
  }    // End of loop for iatom

  return nList;
}

/********************************************/ /**
 *  Function for building neighbour lists for each
 particle. Inefficient O(n^2) implementation. This will only work with one type
 of atom Half neighbour list
 ***********************************************/
std::vector<std::vector<int>> nneigh::halfNeighList(
    double rcutoff, molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    int typeI) {
  std::vector<std::vector<int>>
      nList;       // Vector of vectors of the neighbour list
  double r_ij;     // Distance between iatom and jatom
  int iatomIndex;  // Atomic ID of the atom with index iatom
  int jatomIndex;  // Atomic ID of the atom with index jatom
  int indexYay;
  std::vector<int> tempListIatom;

  // Initialize and fill the first element with the current atom ID whose
  // neighbour list will be filled
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    // Find the atom ID (key) given the index or iatom (value)
    auto itr = std::find_if(
        yCloud->idIndexMap.begin(), yCloud->idIndexMap.end(),
        [&iatom](const std::pair<int, int> &p) { return p.second == iatom; });
    // If found:
    if (itr == yCloud->idIndexMap.end()) {
      std::cerr << "Something is wrong with your idIndexMap!\n";
      continue;
    } else {
      iatomIndex = itr->first;
    }  // End of finding the atom ID to fill as the first element in the
       // neighbour list
    nList.push_back(std::vector<int>());  // Empty vector for the index iatom
    // Fill the first element with the atom ID of iatom itself
    nList[iatom].push_back(iatomIndex);
  }  // end of init

  // Loop through every iatom and find nearest neighbours within rcutoff
  for (int iatom = 0; iatom < yCloud->nop - 1; iatom++) {
    if (yCloud->pts[iatom].type != typeI) {
      continue;
    }
    // Loop through the other atoms
    for (int jatom = iatom + 1; jatom < yCloud->nop; jatom++) {
      if (yCloud->pts[jatom].type != typeI) {
        continue;
      }
      // If the distance is greater than rcutoff, continue
      r_ij = gen::periodicDist(yCloud, iatom, jatom);
      if (r_ij > rcutoff) {
        continue;
      }

      // Get the atom IDs for iatom and jatom
      auto gotI = std::find_if(
          yCloud->idIndexMap.begin(), yCloud->idIndexMap.end(),
          [&iatom](const std::pair<int, int> &p) { return p.second == iatom; });
      if (gotI == yCloud->idIndexMap.end()) {
        std::cerr << "Something is wrong with your idIndexMap!\n";
        return nList;
      } else {
        iatomIndex = gotI->first;
      }  // End of finding the atom ID for iatom
      // Find the atom ID of jatom
      auto gotJ = std::find_if(
          yCloud->idIndexMap.begin(), yCloud->idIndexMap.end(),
          [&jatom](const std::pair<int, int> &p) { return p.second == jatom; });
      if (gotJ == yCloud->idIndexMap.end()) {
        std::cerr << "Something is wrong with your idIndexMap!\n";
        return nList;
      } else {
        jatomIndex = gotJ->first;
      }  // End of finding the atom ID for jatom
      // Update the neighbour indices with atom IDs for iatom and jatom both
      // (full list)
      nList[iatom].push_back(jatomIndex);

    }  // End of loop through jatom
  }    // End of loop for iatom

  return nList;
}

/********************************************/ /**
 *  Function for getting the neighbour list by index instead of by atom ID.
 The ordering is with respect to the pointCloud with the coordinates.The first
 element is the atom for which the other atom indices are neighbours For
 example, if the neighbours of 1 are 2, 3, 4 the sub-vector would have 1 2 3 4
 ***********************************************/
std::vector<std::vector<int>> nneigh::neighbourListByIndex(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList) {
  //
  std::vector<std::vector<int>>
      indexNlist;              // Desired neighbour list of indices
  int iatomID, jatomID;        // Atom IDs
  int iatomIndex, jatomIndex;  // Indices of iatom and jatom
  int nnumNeighbours;          // Number of nearest neighbours

  // Loop through every atom whose neighbours are contained in the neighbour
  // list
  for (int iatom = 0; iatom < nList.size(); iatom++) {
    iatomID = nList[iatom][0];  // Atom ID
    // Get the index of iatom
    auto gotI = yCloud->idIndexMap.find(iatomID);
    if (gotI != yCloud->idIndexMap.end()) {
      iatomIndex = gotI->second;
    }  // found iatomIndex
    //
    nnumNeighbours = nList[iatomIndex].size() - 1;
    // Update the new neighbour list
    indexNlist.push_back(
        std::vector<int>());  // Empty vector for the index iatom
    // Fill the first element with the atom ID of iatom itself
    indexNlist[iatom].push_back(iatomIndex);
    //
    // Loop through the neighbours of iatom
    for (int jatom = 1; jatom <= nnumNeighbours; jatom++) {
      jatomID = nList[iatomIndex][jatom];  // Atom ID of neighbour
      //
      // Get the index of the j^th atom
      auto gotJ = yCloud->idIndexMap.find(jatomID);
      if (gotJ != yCloud->idIndexMap.end()) {
        jatomIndex = gotJ->second;
      }  // found jatomIndex
      // Add to the neighbour list
      indexNlist[iatom].push_back(jatomIndex);
    }  // end of loop through neighbours
  }    // end of loop through every atom

  // Return the new neighbour list
  return indexNlist;
}