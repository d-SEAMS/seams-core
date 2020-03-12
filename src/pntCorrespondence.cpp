#include <pntCorrespondence.hpp>

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

// ---------------------------------------------------
// REFERENCE POINT SETS FOR BULK ICES
// Fills up an eigen matrix point set for a reference cage, saved in the
// templates folder relative to the top-level directory NOT NEEDED MAYBE
Eigen::MatrixXd pntToPnt::getPointSetCage(ring::strucType type) {
  //
  Eigen::MatrixXd pointSet(12, 3); // Reference point set (Eigen matrix)
  molSys::PointCloud<molSys::Point<double>, double>
      setCloud; // PointCloud for holding the reference point values

  if (type == ring::HCbasal) {
    // Read in the XYZ file
    std::string fileName = "templates/hc.xyz";
    //
    sinp::readXYZ(fileName, &setCloud);
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

/********************************************/ /**
 *  Matches the order of the basal rings of an HC
 or a potential HC. The goal is to find which elements are bonded to which and
 the relative order
 ***********************************************/
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
  // Now that the nearest neighbours have been found, find the reversed order
  std::vector<int> tempBasal2; // temporary vector for basal2 matched elements
  tempBasal2.resize(ringSize); // Init to 6 elements
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
      }                              // end of wrap-around
      tempBasal2[i] = basal2[index]; // fill up values
    }                                // end of filling up tempBasal2
  }                                  // the original order is correct
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
        }                              // end of wrap-around
        tempBasal2[i] = basal2[index]; // fill up values
      }                                // end of filling up tempBasal2

    }
    //
    else {
      // std::cerr << "not an HC\n";
      return 1;
    }
    //
  } // the reversed order is correct!
  // ------------------------------
  // If l1 is a neighbour, both basal rings can start from the next element
  if (neighOne) {
    for (int i = 0; i < ringSize; i++) {
      index = 1 + i;
      if (index >= ringSize) {
        index -= ringSize;
      }
      (*matchedBasal1)[i] = basal1[i];
      (*matchedBasal2)[i] = tempBasal2[i];
    } // end of filling of basal1 and basal2
    return 0;
  } // offset by 1 for neighOne
  // No offset
  else if (neighTwo) {
    //
    (*matchedBasal1) = basal1;
    (*matchedBasal2) = tempBasal2;
    return 0;
  } // end of filling
  //

  // std::cerr << "Function should not reach this point.\n";
  return 1;
} // end of the function

/********************************************/ /**
 *  Fills up an eigen matrix point set using the basal rings basal1 and
 basal2,
 changing the order of the point set by filling up from the startingIndex
 (starting from 0 to 5)
 ***********************************************/
Eigen::MatrixXd pntToPnt::changeHexCageOrder(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> basal1, std::vector<int> basal2, int startingIndex) {
  //
  Eigen::MatrixXd pointSet(12, 3);
  int iatomIndex, jatomIndex; // Current atom index in yCloud, according to
                              // basal1 and basal2 respectively
  int iPnt;                   // Current index in the Eigen matrix pointSet
  int cageSize = 12;          // Number of points in the cage
  std::vector<int> newBasal1, newBasal2;
  //
  // Checks and balances
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

  // Loop through the points
  for (int i = 0; i < 6; i++) {
    // basal1
    iatomIndex = newBasal1[i]; // Atom index to be filled
    pointSet(i, 0) = yCloud->pts[iatomIndex].x;
    pointSet(i, 1) = yCloud->pts[iatomIndex].y;
    pointSet(i, 2) = yCloud->pts[iatomIndex].z;
    //
    // basal2
    jatomIndex = newBasal2[i];
    pointSet(i + 6, 0) = yCloud->pts[jatomIndex].x;
    pointSet(i + 6, 1) = yCloud->pts[jatomIndex].y;
    pointSet(i + 6, 2) = yCloud->pts[jatomIndex].z;
    //
  } // end of loop

  return pointSet;
}