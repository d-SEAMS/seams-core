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

#include <pntCorrespondence.hpp>

/**
 * @details Fills up an eigen matrix point set a reference ring, which is a
 * regular n-gonal polygon, constructed with radius 1 by default; where n is the
 * number of nodes in the ring
 */
Eigen::MatrixXd pntToPnt::getPointSetRefRing(int n, int axialDim) {
  //
  Eigen::MatrixXd pointSet(n, 3); // Output point set for a regular polygon
  std::vector<int> dims;          // Apart from the axial dimension

  // Set the axial dimension to 0, and fill up the vector of the other
  // dimensions
  for (int k = 0; k < 3; k++) {
    if (k != axialDim) {
      dims.push_back(k);
    } // other dimensions
  }   // end of setting the axial dimension

  // To get the vertices of a regular polygon, use the following formula:
  // x[i] = r * cos(2*pi*i/n)
  // y[i] = r * sin(2*pi*i/n)
  // where n is the number of points and the index i is 0<=i<=n

  // Loop through every particle
  for (int i = 0; i < n; i++) {
    pointSet(i, dims[0]) = cos((2.0 * gen::pi * i) / n); // x
    pointSet(i, dims[1]) = sin((2.0 * gen::pi * i) / n); // y
    // Set the axial dimension to zero
    pointSet(i, axialDim) = 0.0;
  } // end of loop through all the points

  return pointSet;
} // end of function

/**
 * @details Creates an eigen matrix for the points of a prism block, constructed
 * from the points of a perfect polygon of radius 1, given the basal rings and
 * axial dimension
 */
Eigen::MatrixXd pntToPnt::createPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    const Eigen::MatrixXd &refPoints, int ringSize, std::vector<int> basal1,
    std::vector<int> basal2) {
  //
  int nop = ringSize * 2;           // Number of particles in the prism block
  Eigen::MatrixXd pointSet(nop, 3); // Output point set for a regular prism
  int axialDim;                     // Has a value of 0, 1 or 2 for x, y or z
  double avgRadius; // Average radius within which the basal ring points lie
  double avgHeight; // Get the average height of the prism
  int iBasalOne,
      iBasalTwo; // Indices for the basal1 and basal2 points in the point set
  // --------------------------------------
  // Get the axial dimension
  // The axial dimension will have the largest box length
  // Index -> axial dimension
  // 0 -> x dim
  // 1 -> y dim
  // 2 -> z dim
  axialDim = std::max_element(yCloud->box.begin(), yCloud->box.end()) -
             yCloud->box.begin();
  // --------------------------------------
  // Get the average radius
  avgRadius = pntToPnt::getRadiusFromRings(yCloud, basal1, basal2);
  // --------------------------------------
  // Get the average height of the prism
  avgHeight = pntToPnt::getAvgHeightPrismBlock(yCloud, basal1, basal2);
  // --------------------------------------
  // Fill in the matrix of the point set
  //
  // Fill basal1 first, and then basal2
  for (int i = 0; i < ringSize; i++) {
    iBasalOne = i;            // index in point set for basal1 point
    iBasalTwo = i + ringSize; // index in point set for basal2 point
    // Fill up the dimensions
    for (int k = 0; k < 3; k++) {
      // For the axial dimension
      if (k == axialDim) {
        pointSet(iBasalOne, k) = 0.0;       // basal1
        pointSet(iBasalTwo, k) = avgHeight; // basal2
      } // end of filling up the axial dimension
      else {
        pointSet(iBasalOne, k) = avgRadius * refPoints(i, k); // basal1
        pointSet(iBasalTwo, k) = avgRadius * refPoints(i, k); // basal2
      } // fill up the non-axial dimension coordinate dimension
    }   // Fill up all the dimensions
  }     // end of filling in the point set
  // --------------------------------------
  return pointSet;
} // end of function

/**
 * @details Calculate the average radial distance for the basal rings,
 * calculated from the centroid of each basal ring
 */
double pntToPnt::getRadiusFromRings(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2) {
  //
  double avgRadius = 0.0;
  std::vector<double> centroid1, centroid2;
  std::vector<double> dist;     // Distances
  int ringSize = basal1.size(); // Number of nodes in the basal rings
  double r1, r2;                // Distance from the centroid
  int iatom; // Atom index in yCloud for the current point in basal1
  int jatom; // Atom index in yCloud for the current point in basal2
  // -----------------------
  // Calculate the centroid for basal1 and basal2
  centroid1.resize(3);
  centroid2.resize(3); // init
  // Loop through basal1 and basal2
  for (int i = 0; i < ringSize; i++) {
    // For basal1
    centroid1[0] += yCloud->pts[basal1[i]].x;
    centroid1[1] += yCloud->pts[basal1[i]].y;
    centroid1[2] += yCloud->pts[basal1[i]].z;
    //
    // For basal2
    centroid2[0] += yCloud->pts[basal2[i]].x;
    centroid2[1] += yCloud->pts[basal2[i]].y;
    centroid2[2] += yCloud->pts[basal2[i]].z;
    //
  } // end of loop through basal1 and basal2
  // Normalize by the number of nodes
  for (int k = 0; k < 3; k++) {
    centroid1[k] /= ringSize;
    centroid2[k] /= ringSize;
  } // end of dividing by the number
  // Calculated the centroid for basal1 and basal2!
  // -----------------------
  // Calculate the distances of each point from the respective centroid
  for (int i = 0; i < ringSize; i++) {
    // Get the current point in basal1 and basal2
    iatom = basal1[i];
    jatom = basal2[i];
    // Get the distance from the respective centroid
    // basal1
    r1 = gen::unWrappedDistFromPoint(yCloud, iatom, centroid1);
    // basal2
    r2 = gen::unWrappedDistFromPoint(yCloud, jatom, centroid2);
    // Update the distances
    dist.push_back(r1);
    dist.push_back(r2);
  } // Loop through every point in basal1 and basal2
  // -----------------------
  // Calculate the average, excluding the outliers
  avgRadius = gen::getAverageWithoutOutliers(dist);
  // -----------------------
  return avgRadius;
} // end of function

/**
 * @details Calculate the average height of the prism block, calculated using
 * the basal rings of the prism and the axial dimension
 */
double pntToPnt::getAvgHeightPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2) {
  //
  double avgHeight = 0.0;
  double r_ij;                  // Distance between a point in basal1 and basal2
  std::vector<double> dist;     // Distances
  int ringSize = basal1.size(); // Number of nodes in the basal rings
  int iatom; // Atom index in yCloud for the current point in basal1
  int jatom; // Atom index in yCloud for the current point in basal2
  // -----------------------
  // Calculate the distances of each point from the respective centroid
  for (int i = 0; i < ringSize; i++) {
    // Get the current point in basal1 and basal2
    iatom = basal1[i];
    jatom = basal2[i];
    // Get the distance between a point in basal1 and the corresponding point in
    // basal2
    r_ij = gen::periodicDist(yCloud, iatom, jatom);
    // Update the distances
    dist.push_back(r_ij);
  } // Loop through every point in basal1 and basal2
  // -----------------------
  // Calculate the average, excluding the outliers
  avgHeight = gen::getAverageWithoutOutliers(dist);
  // -----------------------
  return avgHeight;
} // end of function

/**
 * @details Get the relative ordering of a pair of basal rings for a deformed
 * prism/perfect prism. Outputs a vector of vectors of indices, such that the
 * first vector is for the first basal ring, and the second vector is for the
 * second basal ring. The input neighbour list is with respect to indices, not
 * IDs
 */
int pntToPnt::relOrderPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2,
    std::vector<std::vector<int>> nList, std::vector<int> *outBasal1,
    std::vector<int> *outBasal2) {
  //
  int ringSize = basal1.size(); // Number of nodes in basal1 and basal2
  int nBonds;       // Number of bonds between two parallel basal rings
  bool isNeighbour; // Bool for checking if two atoms are neighbours or not
  int l_k, m_k;     // Elements in basal1 and basal2
  bool isClock, isAntiClock; // Clockwise and anti-clockwise ordering of basal1
                             // and basal2
  int iatom, jatom; // Index in the ring to start from, for basal1 and basal2
  int currentIatom, currentJatom;
  // ---------------
  // Find the nearest neighbours of basal1 elements in basal2
  nBonds = 0;
  isNeighbour = false;
  // Loop through every element of basal1
  for (int l = 0; l < ringSize; l++) {
    l_k = basal1[l]; // This is the atom particle C++ index

    // Search for the nearest neighbour of l_k in basal2
    // Loop through basal2 elements
    for (int m = 0; m < ringSize; m++) {
      m_k = basal2[m]; // Atom index to find in the neighbour list of iatom

      // Find m_k inside l_k neighbour list
      auto it = std::find(nList[l_k].begin() + 1, nList[l_k].end(), m_k);

      // If the element has been found, for l1
      if (it != nList[l_k].end()) {
        isNeighbour = true;
        iatom = l; // index of basal1
        jatom = m; // index of basal2
        break;
      } // found element

    } // end of loop through all atomIDs in basal2

    if (isNeighbour) {
      break;
    } // nearest neighbour found
  }   // end of loop through all the atomIDs in basal1

  if (!isNeighbour) {
    std::cerr << "Something is wrong with your deformed prism.\n";
    // Error handling
    return 1;
  }
  // ---------------------------------------------------
  // Find out if the order of basal2 is 'clockwise' or 'anticlockwise'
  isClock = false; // init
  isAntiClock = false;

  // atom index in the ring
  int tempJfor, tempJback;

  tempJfor = jatom + 1;
  tempJback = jatom - 1;

  if (jatom == ringSize - 1) {
    tempJfor = 0;
    tempJback = ringSize - 2;
  }
  if (jatom == 0) {
    tempJfor = 1;
    tempJback = ringSize - 1;
  }

  int forwardJ = basal2[tempJfor];
  int backwardJ = basal2[tempJback];
  int currentI = basal1[iatom];

  // Check clockwise
  double distClock = gen::periodicDist(yCloud, currentI, forwardJ);
  double distAntiClock = gen::periodicDist(yCloud, currentI, backwardJ);

  // Clockwise
  if (distClock < distAntiClock) {
    isClock = true;
  } // end of clockwise check
  // Anti-clockwise
  if (distAntiClock < distClock) {
    isAntiClock = true;
  } // end of anti-clockwise check
  // Some error
  if (isClock == false && isAntiClock == false) {
    // std::cerr << "The points are equidistant.\n";
    // Error handling
    return 1;
  } // end of error handling
  // ---------------------------------------------------
  // Get the order of basal1 and basal2
  for (int i = 0; i < ringSize; i++) {
    currentIatom = iatom + i;
    if (currentIatom >= ringSize) {
      currentIatom -= ringSize;
    } // end of basal1 element wrap-around

    // In clockwise order
    if (isClock) {
      currentJatom = jatom + i;
      if (currentJatom >= ringSize) {
        currentJatom -= ringSize;
      } // wrap around
    }   // end of clockwise update
    else {
      currentJatom = jatom - i;
      if (currentJatom < 0) {
        currentJatom += ringSize;
      } // wrap around
    }   // end of anti-clockwise update

    // Add to outBasal1 and outBasal2 now
    (*outBasal1).push_back(basal1[currentIatom]);
    (*outBasal2).push_back(basal2[currentJatom]);
  } //
  //

  return 0;
} // end of function

/**
 * @details Get the relative ordering of a pair of basal rings for a deformed
 * prism/perfect prism. Outputs a vector of vectors of indices, such that the
 * first vector is for the first basal ring, and the second vector is for the
 * second basal ring. The input neighbour list is with respect to indices, not
 * IDs
 */
int pntToPnt::relOrderPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2,
    std::vector<int> *outBasal1, std::vector<int> *outBasal2) {
  //
  int ringSize = basal1.size(); // Number of nodes in basal1 and basal2
  int nBonds;       // Number of bonds between two parallel basal rings
  bool isNeighbour; // Bool for checking if two atoms are neighbours or not
  int l_k, m_k;     // Elements in basal1 and basal2
  bool isClock, isAntiClock; // Clockwise and anti-clockwise ordering of basal1
                             // and basal2
  int iatom, jatom; // Index in the ring to start from, for basal1 and basal2
  int currentIatom, currentJatom;
  // ---------------
  // Find the nearest neighbours of basal1 elements in basal2
  nBonds = 0;
  double infHugeNumber = 100000;
  double leastDist = infHugeNumber;
  int index = -1; // starting index
  // For the first element of basal1:

  l_k = basal1[0]; // This is the atom particle C++ index

  // Search for the nearest neighbour of l_k in basal2
  // Loop through basal2 elements
  for (int m = 0; m < ringSize; m++) {
    m_k = basal2[m]; // Atom index to find in the neighbour list of iatom

    // Calculate the distance
    double dist = gen::periodicDist(yCloud, l_k, m_k);

    // Update the least distance
    if (leastDist > dist) {
      leastDist = dist; // This is the new least distance
      index = m;
    } // end of update of the least distance

  } // found element

  // If the element has been found, for l1
  if (leastDist < infHugeNumber) {
    isNeighbour = true;
    iatom = 0;     // index of basal1
    jatom = index; // index of basal2
  }                // end of check
  else {
    std::cerr << "Something is wrong with your deformed prism.\n";
    // Error handling
    return 1;
  }
  // ---------------------------------------------------
  // Find out if the order of basal2 is 'clockwise' or 'anticlockwise'
  isClock = false; // init
  isAntiClock = false;

  // atom index in the ring
  int tempJfor, tempJback;

  tempJfor = jatom + 1;
  tempJback = jatom - 1;

  if (jatom == ringSize - 1) {
    tempJfor = 0;
    tempJback = ringSize - 2;
  }
  if (jatom == 0) {
    tempJfor = 1;
    tempJback = ringSize - 1;
  }

  int forwardJ = basal2[tempJfor];
  int backwardJ = basal2[tempJback];
  int currentI = basal1[iatom];

  // Check clockwise
  double distClock = gen::periodicDist(yCloud, currentI, forwardJ);
  double distAntiClock = gen::periodicDist(yCloud, currentI, backwardJ);

  // Clockwise
  if (distClock < distAntiClock) {
    isClock = true;
  } // end of clockwise check
  // Anti-clockwise
  if (distAntiClock < distClock) {
    isAntiClock = true;
  } // end of anti-clockwise check
  // Some error
  if (isClock == false && isAntiClock == false) {
    // std::cerr << "The points are equidistant.\n";
    // Error handling
    return 1;
  } // end of error handling
  // ---------------------------------------------------
  // Get the order of basal1 and basal2
  for (int i = 0; i < ringSize; i++) {
    currentIatom = iatom + i;
    if (currentIatom >= ringSize) {
      currentIatom -= ringSize;
    } // end of basal1 element wrap-around

    // In clockwise order
    if (isClock) {
      currentJatom = jatom + i;
      if (currentJatom >= ringSize) {
        currentJatom -= ringSize;
      } // wrap around
    }   // end of clockwise update
    else {
      currentJatom = jatom - i;
      if (currentJatom < 0) {
        currentJatom += ringSize;
      } // wrap around
    }   // end of anti-clockwise update

    // Add to outBasal1 and outBasal2 now
    (*outBasal1).push_back(basal1[currentIatom]);
    (*outBasal2).push_back(basal2[currentJatom]);
  } //
  //

  return 0;
} // end of function

/**
 * @details Fill up an Eigen Matrix for a prism basal ring from an input vector
 * of atom indices
 */
Eigen::MatrixXd pntToPnt::fillPointSetPrismRing(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basalRing, int startingIndex) {
  //
  //
  int dim = 3;                        // Number of dimensions (3)
  int ringSize = basalRing.size();    // Number of nodes in the ring
  int nop = basalRing.size();         // Number of particles in the basal ring
  Eigen::MatrixXd pointSet(nop, dim); // Output set of 3D points
  // Indices for the first and second basal rings being filled
  int currentPosition; // Current index in basalRing
  int index;           // Index in yCloud

  // Check
  if (startingIndex >= ringSize || startingIndex < 0) {
    startingIndex = 0;
  } // end of check

  // Beginning from the starting index, get the points in the basal ring
  for (int i = 0; i < nop; i++) {
    //
    // Getting currentIndex
    currentPosition = startingIndex + i;
    // Wrapping around for the ring
    if (currentPosition >= ringSize) {
      currentPosition -= ringSize;
    } // end of wrap-around
    //
    // -------------------
    // Basal ring points
    index = basalRing[currentPosition];    // Index of the current point
    pointSet(i, 0) = yCloud->pts[index].x; // x coord
    pointSet(i, 1) = yCloud->pts[index].y; // y coord
    pointSet(i, 2) = yCloud->pts[index].z; // z coord
  } // end of point filling from the relative ordering vector of vectors

  // Return the set of points
  return pointSet;
} // end of function

/**
 * @details Fill up an Eigen Matrix for a prism block from input vectors of atom
 * indices of the basal rings
 */
Eigen::MatrixXd pntToPnt::fillPointSetPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2, int startingIndex) {
  //
  //
  int dim = 3;                        // Number of dimensions (3)
  int ringSize = basal1.size();       // Number of nodes in the ring
  int nop = ringSize * 2;             // Number of particles in prism block
  Eigen::MatrixXd pointSet(nop, dim); // Output set of 3D points
  // Indices for the first and second basal rings being filled
  int currentPosition;        // Current position in the basal ring vectors
  int iatom, jatom;           // Current index in the point set
  int iatomIndex, jatomIndex; // Index in yCloud corresponding to the basal1
                              // and basal2 points

  // Fill up basal1 points, followed by basal2 points

  // Check
  if (startingIndex >= ringSize || startingIndex < 0) {
    startingIndex = 0;
  } // end of check

  // Beginning from the starting index, get the points in the basal ring
  for (int i = 0; i < ringSize; i++) {
    //
    // Getting currentIndex
    currentPosition = startingIndex + i;
    // Wrapping around for the ring
    if (currentPosition >= ringSize) {
      currentPosition -= ringSize;
    } // end of wrap-around
    //
    // -------------------
    // Basal1 ring points
    iatomIndex =
        basal1[currentPosition]; // Index of the current point in basal1
    iatom = i;                   // index in the point set being filled
    pointSet(iatom, 0) = yCloud->pts[iatomIndex].x; // x coord
    pointSet(iatom, 1) = yCloud->pts[iatomIndex].y; // y coord
    pointSet(iatom, 2) = yCloud->pts[iatomIndex].z; // z coord
    // -------------------
    // Basal2 ring points
    jatomIndex =
        basal2[currentPosition]; // Index of the current point in basal2
    jatom = i + ringSize;        // index in the point set being filled
    pointSet(jatom, 0) = yCloud->pts[jatomIndex].x; // x coord
    pointSet(jatom, 1) = yCloud->pts[jatomIndex].y; // y coord
    pointSet(jatom, 2) = yCloud->pts[jatomIndex].z; // z coord
  } // end of point filling from the relative ordering vector of vectors

  // Return the set of points
  return pointSet;
} // end of function

// ---------------------------------------------------
/**
 * @details REFERENCE POINT SETS FOR BULK ICES
 * Fills up an eigen matrix point set for a reference cage, saved in the
 * templates folder relative to the top-level directory NOT NEEDED MAYBE
 */
Eigen::MatrixXd pntToPnt::getPointSetCage(ring::strucType type) {
  //
  Eigen::MatrixXd pointSet(12, 3); // Reference point set (Eigen matrix)
  molSys::PointCloud<molSys::Point<double>, double>
      setCloud; // PointCloud for holding the reference point values

  if (type == ring::strucType::HCbasal) {
    // Read in the XYZ file
    std::string fileName = "templates/hc.xyz";
    //
    sinp::readXYZ(fileName);
    int n = setCloud.nop; // Number of points in the reference point set
    Eigen::MatrixXd pointSet(n, 3); // Output point set for a regular polygon

    // Loop through every particle
    for (int i = 0; i < n; i++) {
      pointSet(i, 0) = setCloud.pts[i].x; // x
      pointSet(i, 1) = setCloud.pts[i].y; // y
      pointSet(i, 2) = setCloud.pts[i].z; // z
    } // end of looping through every particle

    // Return the filled point set
    return pointSet;
  } // end of reference point set for HCs
  // else {
  //   // ddc
  //   Eigen::MatrixXd pointSet(14, 3);  // Output point set for a regular
  //   polygon return pointSet;
  // }  // DDCs

  // Eigen::MatrixXd pointSet(14, 3);  // Output point set for a regular polygon
  // // Never happens
  return pointSet;

} // end of function

/**
 * @details Matches the order of the basal rings of an HC
 * or a potential HC. The goal is to find which elements are bonded to which and
 * the relative order
 */
int pntToPnt::relOrderHC(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2,
    std::vector<std::vector<int>> nList, std::vector<int> *matchedBasal1,
    std::vector<int> *matchedBasal2) {
  //
  int l1 = basal1[0];           // First element of basal1
  int l2 = basal1[1];           // Second element of basal1
  int ringSize = basal1.size(); // Number of nodes in the basal rings
  bool neighOne, neighTwo;      // Basal2 element is the neighbour of l1 or l2
  bool neighbourFound;          // neighbour found
  int m_k;               // Element of basal2 which is a neighbour of l1 or l2
  int m_kIndex;          // Index of basal2 which is a neighbour of l1 or l2
  int iatom;             // Current element of basal2 being searched for
  int nextBasal1Element; // Next bonded element of basal1
  int nextBasal2Element; // Next bonded element of basal2
  bool isReversedOrder;  // basal2 is reversed wrt basal1
  int index;

  (*matchedBasal1).resize(ringSize);
  (*matchedBasal2).resize(ringSize);

  // Search to see if l1 or l2 is a neighbour
  // -------------------
  // Searching for the neighbours of l1 or l2
  for (int i = 0; i < ringSize; i++) {
    iatom = basal2[i]; // Element of basal2
    // Search for the current element in the neighbour list of l1
    auto it = std::find(nList[l1].begin() + 1, nList[l1].end(), iatom);
    // If iatom is the neighbour of l1
    if (it != nList[l1].end()) {
      neighbourFound = true;
      neighOne = true;
      m_k = iatom;  // Element found
      m_kIndex = i; // Index in basal2
      nextBasal1Element = basal1[2];
      break;
    } // iatom is the neighbour of l1
    // l2 neighbour check
    else {
      auto it1 = std::find(nList[l2].begin() + 1, nList[l2].end(), iatom);
      // If iatom is the neighbour of l2
      if (it1 != nList[l2].end()) {
        neighbourFound = true;
        neighOne = false;
        neighTwo = true;
        nextBasal1Element = basal1[3];
        m_k = iatom;  // Element found
        m_kIndex = i; // Index in basal2
        break;
      } // iatom is the neighbour of l2
    }   // Check for the neighbour of l2
  }     // end of search through basal2 elements
  //
  // -------------------
  // If a neighbour was not found, then there is some mistake
  if (!neighbourFound) {
    // std::cerr << "This is not an HC\n";
    return 1;
  }

  // ------------------------------
  //
  // This element should be the nearest neighbour of the corresponding element
  // in basal2
  //
  // Testing the original order
  index = m_kIndex + 2;
  // wrap-around
  if (index >= ringSize) {
    index -= ringSize;
  } // wrap-around
  nextBasal2Element = basal2[index];
  // Search for the next basal element
  auto it = std::find(nList[nextBasal1Element].begin() + 1,
                      nList[nextBasal1Element].end(), nextBasal2Element);
  // If this element is found, then the original order is correct
  if (it != nList[nextBasal1Element].end()) {
    // Fill up the temporary vector with basal2 elements
    for (int i = 0; i < ringSize; i++) {
      index = m_kIndex + i; // index in basal2
      // wrap-around
      if (index >= ringSize) {
        index -= ringSize;
      }                                    // end of wrap-around
      (*matchedBasal2)[i] = basal2[index]; // fill up values
    }                                      // end of filling up tempBasal2
  }                                        // the original order is correct
  else {
    //
    index = m_kIndex - 2;
    // wrap-around
    if (index < 0) {
      index += ringSize;
    } // wrap-around
    nextBasal2Element = basal2[index];
    // Search for the next basal element
    auto it = std::find(nList[nextBasal1Element].begin() + 1,
                        nList[nextBasal1Element].end(), nextBasal2Element);
    // If this element is found, then the original order is correct
    if (it != nList[nextBasal1Element].end()) {
      // Fill up the temporary vector with basal2 elements
      for (int i = 0; i < ringSize; i++) {
        index = m_kIndex - i; // index in basal2
        // wrap-around
        if (index < 0) {
          index += ringSize;
        }                                    // end of wrap-around
        (*matchedBasal2)[i] = basal2[index]; // fill up values
      }                                      // end of filling up tempBasal2

    }
    //
    else {
      // std::cerr << "not an HC\n";
      return 1;
    }
    //
  } // the reversed order is correct!
  // ------------------------------
  // Fill up basal1
  (*matchedBasal1) = basal1;
  return 0;
} // end of the function

/**
 * @details Fills up an eigen matrix point set using the basal rings basal1 and
 * basal2,
 * changing the order of the point set by filling up from the startingIndex
 * (starting from 0 to 5)
 */
Eigen::MatrixXd pntToPnt::changeHexCageOrder(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2, int startingIndex) {
  Eigen::MatrixXd pointSet(12, 3);
  int iatomIndex, jatomIndex; // Current atom index in yCloud, according to
                              // basal1 and basal2 respectively
  int iPnt;                   // Current index in the Eigen matrix pointSet
  int cageSize = 12;          // Number of points in the cage
  std::vector<int> newBasal1, newBasal2;
  std::array<double, 3> dr; // Components of the distance
  int iatomOne;             // Index of the first atom

  // Checks and balances
  //
  if (startingIndex > 5 || startingIndex < 0) {
    startingIndex = 0;
  } // no invalid starting index

  // Change the order
  if (startingIndex > 0) {
    for (int k = 0; k < 6; k++) {
      iPnt = k + startingIndex;
      if (iPnt >= 6) {
        iPnt -= 6;
      } // wrap-around
      newBasal1.push_back(basal1[iPnt]);
      newBasal2.push_back(basal2[iPnt]);
    } // change the order
  }   // end of filling for startingIndex>0
  else {
    newBasal1 = basal1;
    newBasal2 = basal2;
  } // end of filling up the reordered basal rings

  //
  // FIRST POINT
  // basal1
  iatomOne = newBasal1[0];
  pointSet(0, 0) = yCloud->pts[iatomOne].x;
  pointSet(0, 1) = yCloud->pts[iatomOne].y;
  pointSet(0, 2) = yCloud->pts[iatomOne].z;
  // basal2
  jatomIndex = newBasal2[0];
  // Get the distance from basal1
  dr = gen::relDist(yCloud, iatomOne, jatomIndex);

  // basal2
  pointSet(6, 0) = yCloud->pts[iatomOne].x + dr[0];
  pointSet(6, 1) = yCloud->pts[iatomOne].y + dr[1];
  pointSet(6, 2) = yCloud->pts[iatomOne].z + dr[2];
  //
  // Loop through the rest of the points
  for (int i = 1; i < 6; i++) {
    // basal1
    iatomIndex = newBasal1[i]; // Atom index to be filled for basal1
    jatomIndex = newBasal2[i]; // Atom index to be filled for basal2
    //
    // Get the distance from the first atom
    dr = gen::relDist(yCloud, iatomOne, iatomIndex);
    //
    pointSet(i, 0) = yCloud->pts[iatomOne].x + dr[0];
    pointSet(i, 1) = yCloud->pts[iatomOne].y + dr[1];
    pointSet(i, 2) = yCloud->pts[iatomOne].z + dr[2];
    //
    // Get the distance from the first atom
    dr = gen::relDist(yCloud, iatomOne, jatomIndex);
    // basal2
    pointSet(i + 6, 0) = yCloud->pts[iatomOne].x + dr[0];
    pointSet(i + 6, 1) = yCloud->pts[iatomOne].y + dr[1];
    pointSet(i + 6, 2) = yCloud->pts[iatomOne].z + dr[2];
    //
  } // end of loop

  return pointSet;
}

/**
 * @details Reorders the particles of a DDC
 * into a vector, which contains the atom indices of the DDC particles. Here,
 * index is the index in the cageList vector of structs, referring to a specific
 * DDC.
 * The order created is as follows:
 * 1. The first 6 atoms are the atoms of the equatorial ring, such that adjacent
 * particles are bonded to each other (l1, l2, l3, l4, l5, l6).
 * 2. The particles (l1, l3, l5) and (l2, l4, l6) are nearest neighbours of
 * atoms in two different sets of peripheral rings. Each set of three peripheral
 * rings have a second-shell neighbour of the equatorial ring particles (called
 * apex1 and apex2 respectively.)
 * 3. The order is such that the next three atoms in the vector are the first
 * nearest neighbours in the peripheral rings of one triplet (peripheral1),
 * followed by apex1, and three nearest neighbours of the other triplet,
 * followed by apex2. Basically: ((6 equatorial ring atoms), (3 nearest
 * neighbours of [l1,l3,l5] in the peripheral rings), (1 second shell neighbour
 * of [l1,l3,l5]), (first nearest neighbours of [l2,l4,l6]), (second nearest
 * neighbour of [l2,l4,l6]) )
 *
 * Thus, when you want to change the order of the DDC, this is how the order
 * should be wrapped:
 * 1. The first 6 atoms should be wrapped around (i.e. the sixth atom should
 * become the first and so on)
 * 2. The next 3 atoms should be wrapped but only when moved by multiples of 2.
 * 3. The 9th element should not be changed.
 * 4. The next three elements should be wrapped around (multiples of 2), since
 * alternate elements of the equatorial ring are bonded.
 */
std::vector<int> pntToPnt::relOrderDDC(int index,
                                       std::vector<std::vector<int>> rings,
                                       std::vector<cage::Cage> cageList) {
  //
  std::vector<int> ddcOrder; // Order of the particles in the DDC.
  int nop = 14;              // Number of elements in the DDC
  int ringSize = 6;          // Number of nodes in each ring
  int iring, jring;          // Ring indices of the DDC rings
  int iatom;                 // Atom index in the equatorial ring
  std::vector<int> peripheral1, peripheral2; // Vectors holding the neighbours
                                             // of (l1,l3,l5) and (l2,l4,l6)
  int apex1, apex2;
  bool neighbourFound;            // Neighbour found
  int nextI, prevI, nextJ, prevJ; // elements around iatom and jatom
  int jatomIndex, atomIndex;

  // Add the equatorial ring particles
  //
  iring = cageList[index]
              .rings[0]; // Equatorial ring index in rings vector of vectors
  //
  // Loop through all the atoms of the equatorial ring
  for (int i = 0; i < ringSize; i++) {
    ddcOrder.push_back(rings[iring][i]);
  } // end of adding equatorial ring particles
  //
  // ------------------------------
  // Find the neighbouring atoms of (l1, l3 and l5) and (l2,l4,l6) by looping
  // through the peripheral rings
  //
  for (int i = 0; i < ringSize; i++) {
    // Init
    iatom = ddcOrder[i]; // element of the equatorial ring
    // Get the next and previous elements of iatom
    // next element
    if (i == ringSize - 1) {
      nextI = ddcOrder[0];
      prevI = ddcOrder[i - 1];
    } // if last element
    else if (i == 0) {
      nextI = ddcOrder[i + 1];
      prevI = ddcOrder[5]; // wrapped
    } else {
      nextI = ddcOrder[i + 1];
      prevI = ddcOrder[i - 1];
    }
    // ------------------
    // Search all the peripheral rings for iatom, get the nearest neighbours and
    // apex1 and apex2, which are bonded to the nearest neighbours of the
    // equatorial ring (thus they are second nearest neighbours)
    for (int j = 1; j <= ringSize; j++) {
      //
      jring = cageList[index].rings[j]; // Peripheral ring
      // Search for iatom
      //
      auto it = std::find(rings[jring].begin(), rings[jring].end(), iatom);
      // If iatom was found in jring peripheral ring
      if (it != rings[jring].end()) {
        // ----------------
        // Atom index for iatom found
        jatomIndex = std::distance(rings[jring].begin(), it);
        // ----------------
        // Get the next and previous element
        //
        // Next atom
        atomIndex = jatomIndex + 1;
        // wrap-around
        if (atomIndex >= ringSize) {
          atomIndex -= ringSize;
        } // end of wrap-around
        nextJ = rings[jring][atomIndex];
        // Previous atom
        atomIndex = jatomIndex - 1;
        if (atomIndex < 0) {
          atomIndex += ringSize;
        } // end of wrap-around
        prevJ = rings[jring][atomIndex];
        // ----------------
        // Check to see if nextJ or prevJ are different from nextI and prevI
        //
        // Checking prevJ
        if (prevJ != nextI && prevJ != prevI) {
          // Add to the vector
          // if even peripheral1
          if (i % 2 == 0) {
            peripheral1.push_back(prevJ);
            // Get apex1 for i=0
            if (i == 0) {
              // Go back two elements from jatomIndex
              atomIndex = jatomIndex - 2;
              // wrap-around
              if (atomIndex < 0) {
                atomIndex += ringSize;
              } // end of wrap-around
              apex1 = rings[jring][atomIndex];
            } // apex1 for i=0
          }   // peripheral1
          // if odd peripheral2
          else {
            peripheral2.push_back(prevJ);
            // Get apex2 for i=0
            if (i == 1) {
              // Go back two elements from jatomIndex
              atomIndex = jatomIndex - 2;
              // wrap-around
              if (atomIndex < 0) {
                atomIndex += ringSize;
              } // end of wrap-around
              apex2 = rings[jring][atomIndex];
            } // apex2 for i=1
          }   // peripheral 2
          // Get apex1 or apex2 for i=0 or i=1
          break;
        } // check if prevJ is the atom to be added, bonded to iatom
        // ----------------------
        //
        // Checking nextJ
        if (nextJ != nextI && nextJ != prevI) {
          // Add to the vector
          // if even peripheral1
          if (i % 2 == 0) {
            peripheral1.push_back(nextJ);
            // Get apex1 for i=0
            if (i == 0) {
              // Go foward two elements from jatomIndex
              atomIndex = jatomIndex + 2;
              // wrap-around
              if (atomIndex >= ringSize) {
                atomIndex -= ringSize;
              } // end of wrap-around
              apex1 = rings[jring][atomIndex];
            } // apex1 for i=0
          }   // peripheral1
          // if odd peripheral2
          else {
            peripheral2.push_back(prevJ);
            // Get apex2 for i=0
            if (i == 1) {
              // Go foward two elements from jatomIndex
              atomIndex = jatomIndex + 2;
              // wrap-around
              if (atomIndex >= ringSize) {
                atomIndex -= ringSize;
              } // end of wrap-around
              apex2 = rings[jring][atomIndex];
            } // apex2 for i=1
          }   // peripheral 2
          // Get apex1 or apex2 for i=0 or i=1
          break;
        } // check if prevJ is the atom to be added, bonded to iatom
        // ----------------
      } // end of check for iatom in jring
    }   // end of search for the peripheral rings
    // ------------------
  } // end of looping through the elements of the equatorial ring
  // ------------------------------
  // Update ddcOrder with peripheral1, followed by apex1, peripheral2 and apex2
  //
  // peripheral1
  for (int i = 0; i < peripheral1.size(); i++) {
    ddcOrder.push_back(peripheral1[i]);
  } // end of adding peripheral1
  //
  // apex1
  ddcOrder.push_back(apex1);
  //
  // peripheral2
  for (int i = 0; i < peripheral2.size(); i++) {
    ddcOrder.push_back(peripheral2[i]);
  } // end of adding peripheral2
  //
  // apex2
  ddcOrder.push_back(apex2);
  //
  return ddcOrder;
} // end of the function

/**
 * @details Fills up an eigen matrix point set using information of the
 * equatorial ring and peripheral rings, embedded in a vector, already filled in
 * relOrderDDC. The order is such that the next three atoms in the vector are
 * the first nearest neighbours in the peripheral rings of one triplet
 * (peripheral1), followed by apex1, and three nearest neighbours of the other
 * triplet, followed by apex2. Basically: ((6 equatorial ring atoms), (3 nearest
 * neighbours of [l1,l3,l5] in the peripheral rings), (1 second shell neighbour
 * of [l1,l3,l5]), (first nearest neighbours of [l2,l4,l6]), (second nearest
 * neighbour of [l2,l4,l6]) )
 *
 * Thus, when you want to change the order of the DDC, this is how the order
 * should be wrapped:
 * 1. The first 6 atoms should be wrapped around (i.e. the sixth atom should
 * become the first and so on)
 * 2. The next 3 atoms should be wrapped but only when moved by multiples of 2.
 * 3. The 9th element should not be changed.
 * 4. The next three elements should be wrapped around (multiples of 2), since
 * alternate elements of the equatorial ring are bonded.
 */
Eigen::MatrixXd pntToPnt::changeDiaCageOrder(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> ddcOrder, int startingIndex) {
  int nop = 14;                     // Number of elements in the DDC
  int ringSize = 6;                 // Six nodes in the rings
  Eigen::MatrixXd pointSet(nop, 3); // Output point set
  int peripheralStartingIndex; // Index using which the elements of peripheral1
                               // and peripheral2 will be wrapped
  std::vector<int> wrappedDDC; // Changed order of the DDC: should be the same
                               // size as the original ddcOrder vector
  int currentIndex;            // Current index
  // Variables for filling the point set
  int iatomIndex, jatomIndex;
  std::array<double, 3> dr; // Components of the distance

  if (startingIndex == 0) {
    wrappedDDC = ddcOrder;
  } // if the order does not have to be changed
  else {
    wrappedDDC.resize(nop); // Should be the same size as ddcOrder
    // ------------------
    // EQUATORIAL RING
    // Change the order of the equatorial ring
    for (int i = 0; i < ringSize; i++) {
      currentIndex = startingIndex + i;
      // Wrap-around
      if (currentIndex >= ringSize) {
        currentIndex -= ringSize;
      } // end of wrap-around
      // Update the equatorial ring
      wrappedDDC[i] = ddcOrder[currentIndex];
    } // end of wrapping the equatorial ring
    // ------------------
    // PERIPHERAL RINGS
    //
    // The index of the peripheral ring starting indices
    if (startingIndex <= 1) {
      peripheralStartingIndex = 0;
    } else if (startingIndex > 1 && startingIndex <= 3) {
      peripheralStartingIndex = 1;
    } else {
      peripheralStartingIndex = 2;
    }
    //
    // Update the portions of the wrappedDDC vector
    for (int i = 6; i < 9; i++) {
      currentIndex = i + peripheralStartingIndex;
      // wrap-around
      if (currentIndex <= 9) {
        currentIndex -= 3;
      } // end of wrap-around
      wrappedDDC[i] = ddcOrder[currentIndex];
    } // peripheral1
    //
    // Update the peripheral2 portions of the wrappedDDC vector
    for (int i = 10; i < 13; i++) {
      currentIndex = i + peripheralStartingIndex;
      // wrap-around
      if (currentIndex <= 13) {
        currentIndex -= 3;
      } // end of wrap-around
      wrappedDDC[i] = ddcOrder[currentIndex];
    } // peripheral2
    // ------------------
    // SECOND-NEAREST NEIGHBOURS
    wrappedDDC[9] = ddcOrder[9];   // apex1
    wrappedDDC[13] = ddcOrder[13]; // apex2
    // ------------------
  } // the order of the DDC has to be changed

  // FILL UP THE EIGEN MATRIX
  iatomIndex = wrappedDDC[0]; // first point
  pointSet(0, 0) = yCloud->pts[iatomIndex].x;
  pointSet(0, 1) = yCloud->pts[iatomIndex].y;
  pointSet(0, 2) = yCloud->pts[iatomIndex].z;
  // Loop through the rest of the equatorial ring points
  for (int i = 1; i < 14; i++) {
    //
    jatomIndex = wrappedDDC[i]; // Atom index to be filled
    // Get the distance from iatomIndex
    dr = gen::relDist(yCloud, iatomIndex, jatomIndex);

    // basal2
    pointSet(i, 0) = yCloud->pts[iatomIndex].x + dr[0];
    pointSet(i, 1) = yCloud->pts[iatomIndex].y + dr[1];
    pointSet(i, 2) = yCloud->pts[iatomIndex].z + dr[2];
  } // end of loop

  return pointSet;
}
