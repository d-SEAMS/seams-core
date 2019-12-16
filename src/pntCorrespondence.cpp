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
Eigen::MatrixXd pntToPnt::getPointSetRefRing(int n) {
  //
  Eigen::MatrixXd pointSet(n, 3); // Output point set for a regular polygon

  // To get the vertices of a regular polygon, use the following formula:
  // x[i] = r * cos(2*pi*i/n)
  // y[i] = r * sin(2*pi*i/n)
  // where n is the number of points and the index i is 0<=i<=n

  // Loop through every particle
  for (int i = 0; i < n; i++) {
    pointSet(i, 0) = cos((2.0 * gen::pi * i) / n); // x
    pointSet(i, 1) = sin((2.0 * gen::pi * i) / n); // y
    pointSet(i, 2) = 0.0;                          // z
  } // end of loop through all the points

  return pointSet;
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

// Fill up an Eigen Matrix from an input vector of atom indices
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