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
