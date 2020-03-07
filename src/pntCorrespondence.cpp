#include <pntCorrespondence.hpp>

// Fills up an eigen matrix point set for an HC, according to an input
// pointCloud, the relative order given by the basal rings, beginning from the
// startingIndex
Eigen::MatrixXd pntToPnt::fillPointSetHC(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> relOrder, int startingIndex) {
  //
  int dim = 3;                          // Number of dimensions (3)
  int ringSize = relOrder[0].size();    // Number of elements in each basal ring
  int nop = ringSize * relOrder.size(); // Number of particles in the HC
  Eigen::MatrixXd pointSet(nop, dim);   // Output set of 3D points
  // Indices for the first and second basal rings being filled
  int iBasalOne, iBasalTwo; // Indices being filled
  int currentPosition;      // Current index in the first or second basal rings
                            // inside relOrder
  int index;                // Index in yCloud

  // Check
  if (startingIndex >= ringSize || startingIndex < 0) {
    startingIndex = 0;
  } // end of check

  // Beginning from the starting index, get the points such that the first basal
  // ring is filled up first, followed by the connected second basal ring
  for (int i = 0; i < ringSize; i++) {
    //
    iBasalOne = i;            // basal one ring is filled first
    iBasalTwo = i + ringSize; // basal two ring is filled after basal one is
                              // filled in sequential order
    // Getting currentIndex
    currentPosition = startingIndex + i;
    // Wrapping around for the ring
    if (currentPosition >= ringSize) {
      currentPosition -= ringSize;
    } // end of wrap-around
    //
    // -------------------
    // basal one points
    index = relOrder[0][currentPosition];          // Index of the current point
    pointSet(iBasalOne, 0) = yCloud->pts[index].x; // x coord
    pointSet(iBasalOne, 1) = yCloud->pts[index].y; // y coord
    pointSet(iBasalOne, 2) = yCloud->pts[index].z; // z coord
    // -------------------
    // basal two points
    index = relOrder[1][currentPosition];          // Index of the current point
    pointSet(iBasalTwo, 0) = yCloud->pts[index].x; // x coord
    pointSet(iBasalTwo, 1) = yCloud->pts[index].y; // y coord
    pointSet(iBasalTwo, 2) = yCloud->pts[index].z; // z coord
  } // end of point filling from the relative ordering vector of vectors

  // Return the set of points
  return pointSet;
} // end function

// Fills up an eigen matrix point set a reference ring, which is a regular
// n-gonal polygon, constructed with radius 1 by default; where n is the number
// of nodes in the ring
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

// Creates an eigen matrix for the points of a prism block, constructed from the
// points of a perfect polygon of radius 1, given the basal rings and axial
// dimension
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

// Calculate the average radial distance for the basal rings, calculated from
// the centroid of each basal ring
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

// Calculate the average height of the prism block, calculated using the basal
// rings of the prism and the axial dimension
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

// Get the relative ordering of a pair of basal rings for a deformed
// prism/perfect prism. Outputs a vector of vectors of indices, such that the
// first vector is for the first basal ring, and the second vector is for the
// second basal ring. The input neighbour list is with respect to indices, not
// IDs
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

// Get the relative ordering of a pair of basal rings for a deformed
// prism/perfect prism. Outputs a vector of vectors of indices, such that the
// first vector is for the first basal ring, and the second vector is for the
// second basal ring. The input neighbour list is with respect to indices, not
// IDs
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

// Fill up an Eigen Matrix for a prism basal ring from an input vector of atom
// indices
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

// Fill up an Eigen Matrix for a prism block from input vectors of atom
// indices of the basal rings
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