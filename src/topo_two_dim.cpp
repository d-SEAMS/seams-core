#include <topo_two_dim.hpp>

// -----------------------------------------------------------------------------------------------------
// MONOLAYER ALGORITHMS
// -----------------------------------------------------------------------------------------------------

/**
 * @details Function that loops through the primitive rings (which is a vector
 * of vectors) of all sizes, upto maxDepth (the largest ring size). The function
 * returns a vector which contains the ring indices of all the rings which
 * constitute prism blocks. The ring IDs correspond to the index of the rings
 * inside the vector of vector rings, starting from 0. This is registered as a
 * Lua function, and is exposed to the user directly. The function calls the
 * following functions internally:
 * - ring::clearRingList (Clears the vector of vectors for rings of a single
 * type, to prevent excessive memory being blocked).
 * - ring::getSingleRingSize (Fill a vector of vectors for rings of a particular
 * ring size).
 * - ring::findPrisms (Now that rings of a particular size have been obtained,
 * prism blocks are found and saved).
 * - topoparam::normHeightPercent (Gets the height% for the prism blocks).
 * - ring::assignPrismType (Assigns a prism type to each atom type).
 * - sout::writePrismNum (Write out the prism information for the current
 * frame).
 * - sout::writeLAMMPSdataAllPrisms (Writes out a LAMMPS data file for the
 * current frame, which can be visualized in OVITO).
 * @param[in] path The string to the output directory, in which files will be
 *  written out.
 * @param[in] rings Row-ordered vector of vectors for rings of a single type.
 * @param[in] nList Row-ordered neighbour list by index.
 * @param[in] yCloud The input PointCloud.
 * @param[in] maxDepth The maximum possible size of the primitive rings.
 * @param[in] sheetArea Area calculated using the two significant dimensions of
 *  the quasi-two-dimensional sheet.
 * @param[in] firstFrame The first frame to be analyzed
 */
int ring::polygonRingAnalysis(
    std::string path, std::vector<std::vector<int>> rings,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int maxDepth,
    double sheetArea, int firstFrame) {
  //
  std::vector<std::vector<int>>
      ringsOneType;           // Vector of vectors of rings of a single size
  int nRings;                 // Number of rings of the current type
  std::vector<int> nRingList; // Vector of the values of the number of prisms
                              // for a particular frame
  std::vector<double> coverageAreaXY; // Coverage area percent on the XY plane
                                      // for a particular n and frame
  std::vector<double> coverageAreaXZ; // Coverage area percent on the XZ plane
                                      // for a particular n and frame
  std::vector<double> coverageAreaYZ; // Coverage area percent on the YZ plane
                                      // for a particular n and frame
  std::vector<double> singleAreas; // Coverage area on the XY, XZ and YZ planes
                                   // for a single ring
  std::vector<int>
      atomTypes; // contains int values for each prism type considered
  // -------------------------------------------------------------------------------
  // Init
  nRingList.resize(
      maxDepth -
      2); // Has a value for every value of ringSize from 3, upto maxDepth
  coverageAreaXY.resize(maxDepth - 2);
  coverageAreaXZ.resize(maxDepth - 2);
  coverageAreaYZ.resize(maxDepth - 2);
  // The atomTypes vector is the same size as the pointCloud atoms
  atomTypes.resize(yCloud->nop, 1); // The dummy or unclassified value is 1
  // -------------------------------------------------------------------------------
  // Run this loop for rings of sizes upto maxDepth
  // The smallest possible ring is of size 3
  for (int ringSize = 3; ringSize <= maxDepth; ringSize++) {
    // Clear ringsOneType
    ring::clearRingList(ringsOneType);
    // Get rings of the current ring size
    ringsOneType = ring::getSingleRingSize(rings, ringSize);
    //
    // Continue if there are zero rings of ringSize
    if (ringsOneType.size() == 0) {
      nRingList[ringSize - 3] = 0; // Update the number of prisms
      // Update the coverage area%
      coverageAreaXY[ringSize - 3] = 0.0;
      coverageAreaXZ[ringSize - 3] = 0.0;
      coverageAreaYZ[ringSize - 3] = 0.0;
      continue;
    } // skip if there are no rings
    //
    // -------------
    // Number of rings with n nodes
    nRings = ringsOneType.size();
    // -------------
    // Now that you have rings of a certain size:
    nRingList[ringSize - 3] = nRings; // Update the number of n-membered rings
    // Update the coverage area% for the n-membered ring phase.
    singleAreas = topoparam::calcCoverageArea(yCloud, ringsOneType, sheetArea);
    // XY plane
    coverageAreaXY[ringSize - 3] = singleAreas[0];
    // XZ plane
    coverageAreaXZ[ringSize - 3] = singleAreas[1];
    // YZ plane
    coverageAreaYZ[ringSize - 3] = singleAreas[2];
    // -------------
  } // end of loop through every possible ringSize

  // Get the atom types for all the ring types
  ring::assignPolygonType(rings, &atomTypes, nRingList);

  // Write out the ring information
  sout::writeRingNum(path, yCloud->currentFrame, nRingList, coverageAreaXY,
                     coverageAreaXZ, coverageAreaYZ, maxDepth, firstFrame);
  // Write out the lammps data file for the particular frame
  sout::writeLAMMPSdataAllRings(yCloud, nList, atomTypes, maxDepth, path);

  return 0;
}