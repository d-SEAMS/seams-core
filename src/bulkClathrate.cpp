//-----------------------------------------------------------------------------------
// d-SEAMS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <bulkClathrate.hpp>

// -----------------------------------------------------------------------------------------------------
// ALGORITHMS FOR CLATHRATES 
// -----------------------------------------------------------------------------------------------------

/**
 * @details Shape-matching algorithms for S2 clathrate. This is registered as a Lua
 * function and is accessible to the user. Internally, this function calls the
 * following functions:
 *   - ring::getSingleRingSize (Saves rings of a single ring size into a new
 * vector of vectors, which is subsequently used for finding DDCs, HCs etc).
 *   - ring::findDDC (Finds the DDCs).
 *   - ring::findHC (Finds the HCs).
 *   - ring::findMixedRings (Finds the mixed rings, which are shared by DDCs and
 * HCs both).
 *   - ring::getStrucNumbers (Gets the number of structures (DDCs, HCs, mixed
 * rings, basal rings, prismatic rings, to be used for write-outs).
 *   - sout::writeTopoBulkData (Writes out the numbers and data obtained for the
 *    current frame).
 *   - ring::getAtomTypesTopoBulk (Gets the atom type for every atom, to be used
 * for printing out the ice types found).
 *   - sout::writeLAMMPSdataTopoBulk (Writes out the atoms, with the classified
 *    types, into a LAMMPS data file, which can be visualized in OVITO).
 *  @param[in] path The file path of the output directory to which output files
 *   will be written.
 *  @param[in] rings Vector of vectors containing the primitive rings. This
 *   contains rings of all sizes.
 *  @param[in] nList Row-ordered neighbour list, by index.
 *  @param[in] yCloud The input PointCloud, with respect to which the indices in
 *   the rings and nList vector of vectors have been saved.
 *  @param[in] templateFileO Path to the LAMMPS trajectory file containing the template
 *   reference structure of the S2 clathrate. 
  *  @param[in] oxygenAtomType Type ID for the oxygen atoms. 
 */
void clath::shapeMatchS2ClathrateSystem(std::string path, std::vector<std::vector<int>> nList, molSys::PointCloud<molSys::Point<double>, double> yCloud, 
  std::string filename, int targetFrame,
  int atomTypeI, bool isSlice, std::array<double, 3> coordLow, std::array<double, 3> coordHigh,
  std::string templateFileO, int oxygenAtomType){
  // Reference structure stuff 
  int dim = 3;   // Number of dimensions
  int nOxy = 28; // Number of O atoms in the clathrate hexadecahedron
  Eigen::MatrixXdRowMajor refPntsO(nOxy, dim); // Reference point set of just O atoms (Eigen matrix)
  molSys::PointCloud<molSys::Point<double>, double>
        refCloud;  // pointCloud for the reference template S2 cage 
  std::vector<std::vector<int>> ringsRef;  // Rings for the reference structure 
  // Candidate structure stuff
  molSys::PointCloud<molSys::Point<double>, double>
        targetCloud;  // pointCloud for the target/candidate S2 cage
  //
  std::vector<std::vector<int>>
      ringsHex;    // Vector of vectors of rings of just 6-membered rings
  std::vector<std::vector<double>> centroidCoord; // Coordinates of the centroid of the molecules
  std::vector<int> cageOAtomIndices; // List of indices of water O atoms which are closest to a THF molecule 
  double maxCutoff = 30.0; // in Angstroms; TODO change 
  // -----------------------------------
  // Build the reference point set
  // This is row-ordered 
  std::tie(refCloud, ringsRef, refPntsO) = clath::buildRefS2Cage(templateFileO, oxygenAtomType);
  // -----------------------------------
  // Get the COM of the encapsulated molecules in the region 
  // where you want to do shape-matching 
  centroidCoord = misc::getCentroidMolecules(filename, targetFrame, atomTypeI, isSlice, coordLow, coordHigh);

  // Loop through the THF centroids
  for (auto& centroidPnt : centroidCoord){
    // Find the 28 closest water molecules 
    cageOAtomIndices = misc::kClosestPoints(yCloud, oxygenAtomType, 
      centroidPnt, nOxy, maxCutoff);
  } // loop through the THF centroid points 
  // Find a test template structure (28 closest water molecules)
  // Find the primitive rings for the candidate cage 

  // Get just the 6-membered rings for a given test template structure 
  // ringsHex = ring::getSingleRingSize(rings, 6);

  // Shape-matching of test cage to the reference cage 

  return;
} // end of function 

/**
 * @details Build a reference SII cage, consisting of 28 water molecules, reading it in from a template
 * file saved in the templates directory 
 * In the reference structure with O H H atoms, H atoms follow the O atoms in the same molecule
 */
std::pair<Eigen::MatrixXdRowMajor, Eigen::MatrixXdRowMajor> 
clath::buildRefS2CageLammpsTrj(std::string filename, std::string filenameO, int oxygenAtomType) {
  //
  // Get the reference point set
  // Row major 
  Eigen::MatrixXdRowMajor refPntsO(28, 3); // Reference point set of just O atoms (Eigen matrix)
  Eigen::MatrixXdRowMajor refPntsWat(84, 3); // Reference point set of O H H water atoms (Eigen matrix)
  molSys::PointCloud<molSys::Point<double>, double>
      waterCloud; // PointCloud for holding the reference point values for all the water molecules
      molSys::PointCloud<molSys::Point<double>, double>
      oCloud; // PointCloud for holding the reference point values for just the O atoms
  std::vector<int> oxyIDorder, watIndexOrder; // Vectors containing atom IDs of the O atoms
  // and the atom indices of the water molecules (H H following the O atom in each molecule; with the same molecule ID)
  int iatomID, jatomID; // LAMMPS atom ID
  int iatomMolID; // molecule ID
  int iatomIndex, jatomIndex; // LAMMPS atom index in waterCloud   
  std::unordered_multimap<int, int> molIDAtomIDmap; // Unordered multimap with molecule IDs of the atoms as the keys and the
                  // atom IDs as the values. More than one atom can have the same molecule ID
  std::vector<int> atomIndexListMolID; // Vector with atom indices with the same molecule ID  
  //
  // read in all the water molecules into the pointCloud setCloud
  //
  // All water molecules 
  waterCloud = sinp::readLammpsTrj(filename, 1, &waterCloud);
  // Just the oxygen atoms 
  oCloud = sinp::readLammpsTrjO(filenameO, 1, &oCloud, oxygenAtomType);
  // --------------
  // Get the unordered map of the atom IDs (keys) and the molecular IDs
  // (values)
  molIDAtomIDmap = molSys::createMolIDAtomIDMultiMap(&waterCloud);
  // --------------
  // Get the order of O atoms in the reference PointCloud 
  for (int iatom = 0; iatom < oCloud.nop; iatom++)
  {
    iatomID = oCloud.pts[iatom].atomID; // current atom ID of the O atom
    // Add to the vector 
    oxyIDorder.push_back(iatomID);
    // Update the Eigen matrix for the reference points of O 
    refPntsO(iatom, 0) = oCloud.pts[iatom].x;
    refPntsO(iatom, 1) = oCloud.pts[iatom].y;
    refPntsO(iatom, 2) = oCloud.pts[iatom].z;
  } // order of oxygen atoms; update of oxyIDorder

  // --------------
  // Go through the O atoms, and add H atoms which belong 
  // to the same molecule 
  for (auto it = oxyIDorder.begin(); it < oxyIDorder.end(); ++it)
  {
    // -------
    // Find the index in waterCloud
    iatomID = *it; // Oxygen atom ID 
    // Get the index from idIndexMap in waterCloud
    auto gotI = waterCloud.idIndexMap.find(iatomID);
    iatomIndex = gotI->second;
    // Update the watIndexOrder vector (according to waterCloud)
    watIndexOrder.push_back(iatomIndex);
    // -------
    // Get the molecule ID and find the H atoms 
    iatomMolID = waterCloud.pts[iatomIndex].molID; 
    // ---
    // Find all atoms with iatomMolID
    atomIndexListMolID = gen::atomIndexWithMolID(waterCloud, iatomMolID, molIDAtomIDmap);
    
    // Go through the vector
    for (int i = 0; i < atomIndexListMolID.size(); ++i)
    {
      jatomIndex = atomIndexListMolID[i]; // atom index 
      if (iatomIndex!=jatomIndex){
        watIndexOrder.push_back(jatomIndex);
      } // skip for O atoms
    } // end of loop through atomIndexListMolID
  }
  // --------------
  // Build the Eigen vector refPntsWat using watIndexOrder
  // which is according to the PointCloud waterCloud
  for (int i = 0; i < watIndexOrder.size(); ++i)
  {
    iatomIndex = watIndexOrder[i]; // Atom index in waterCloud 
    // Update the Eigen matrix for the reference points of water 
    refPntsWat(i, 0) = waterCloud.pts[iatomIndex].x;
    refPntsWat(i, 1) = waterCloud.pts[iatomIndex].y;
    refPntsWat(i, 2) = waterCloud.pts[iatomIndex].z;
  } // end of building the Eigen vector for refPntsWat
  // --------------
  return std::make_pair (refPntsO, refPntsWat);
}

/**
 * @details Build a reference SII cage, consisting of 28 water molecules, reading it in from a template
 * file saved in the templates directory. Only the O atoms are read in.
 * Give the file path for the O atom trajectory 
 * Returns a tuple of a PointCloud of the reference structure, vector of vectors of six-membered rings
 * and the 4 other molecules not part of the 6-membered rings, and an Eigen row matrix of (28,3) dimensions
 * with the coordinates. The 6-membered rings are mutually exclusive and share no common elements.
 * The rings contain atom index values and not atom IDs. 
 * The order in the reference matrix is: 
 * Ring 1, Ring 2, Ring 3, Ring 4, 4 elements not part of the 6-membered rings  
 */
std::tuple<molSys::PointCloud<molSys::Point<double>, double>,std::vector<std::vector<int>> , Eigen::MatrixXdRowMajor> 
clath::buildRefS2Cage(std::string filename, int oxygenAtomType) {
  //
  int nop = 28; // No. of O atoms 
  int dim = 3;
  // Get the reference point set
  // Row major   
  Eigen::MatrixXdRowMajor refPntsO(nop, dim); // Reference point set of just O atoms (Eigen matrix)
  molSys::PointCloud<molSys::Point<double>, double>
      refCloud; // PointCloud for holding the reference point values for the O atoms 
  std::vector<std::vector<int> > rings{ {2, 26, 22, 21, 5, 1},
   {3, 4, 9, 17, 16, 25}, {8, 14, 13, 19, 6, 7}, {12, 20, 23, 27, 15, 11}, {24, 18, 0, 10} }; // Rings (predefined for the reference structure)
  // --------------
  // Get the reference PointCloud 
  refCloud = sinp::readLammpsTrjO(filename, 1, &refCloud, oxygenAtomType);
  // --------------
  // Fill in the Eigen matrix, given the rings vector of vector 
  refPntsO = pntToPnt::fillTargetEigenPointSetFromRings(
    refCloud, rings, nop, dim);
  // --------------
  return std::make_tuple (refCloud, rings, refPntsO);
}

/**
 * @details Add a given vector to the current rings vector of vectors for a 
 * target clathrate structure. Match the reference and target structures. 
 */
void clath::matchClathrateLastRing(std::vector<std::vector<int>> targetRings, std::vector<int> lastTargetVec,
  std::vector<std::vector<int>> refRings, molSys::PointCloud<molSys::Point<double>, double> targetCloud,
  molSys::PointCloud<molSys::Point<double>, double> refCloud,
  std::vector<double> *quat, double *rmsd,
  std::vector<double> *rmsdList, double *scale) {
  //
  int numRings; // Number of rings/vectors to match
  int nop=0; // No. of water molecules/ no. of O atoms; points to match
  int dim = 3; // No. of dimensions 
  // Contains only the number of rings corresponding to the target rings
  std::vector<std::vector<int>> reducedRefRings; 
  // --------------
  // Add the last vector to the target rings vector of vectors
  targetRings.push_back(lastTargetVec);
  // Get the number of elements and rings to match
  numRings = targetRings.size(); // No. of rings to match
  for (auto& iring : targetRings){
    // Number of elements in this ring
    nop += iring.size();
  } // loop through the vector of vectors 
  // 
  // Get the reduced reference rings 
  for (int i = 0; i < numRings; ++i)
  {
    reducedRefRings.push_back(refRings[i]);
  } // building reduced reference rings 
  // --------------
  // Get the Eigen row-major matrices of the target and reference point set 
  Eigen::MatrixXdRowMajor refPnts(nop, dim); // Reference point set (Eigen matrix) 
  Eigen::MatrixXdRowMajor targetPnts(nop, dim); // Target point set (Eigen matrix)
  // Fill in the Eigen matrices, given the rings vectors of vector
  // Reference point set  
  refPnts = pntToPnt::fillTargetEigenPointSetFromRings(refCloud, reducedRefRings, nop, dim);
  // Target point set 
  targetPnts = pntToPnt::fillTargetEigenPointSetFromRings(targetCloud, targetRings, nop, dim);
  // --------------
  // SHAPE-MATCHING 
  absor::hornAbsOrientationRowMajor(refPnts, targetPnts, quat,
            rmsd, rmsdList, scale);
  // --------------
  return;
}

/**
 * @details Function for returning a vector of vectors of the centroid x y z values for molecules, 
 * given a constituent atom type. 
 * @param[in] yCloud The given PointCloud
 * @param[in] molID The molecule ID for which atom index values will be returned
 * @param[in] molIDAtomIDmap Unordered multimap mapping the molecule IDs to the atom IDs for the PointCloud
 * @return A vector containing atom indices with the molID molecule ID, corresponding to the given PointCloud
 */
std::vector<std::vector<double>> 
misc::getCentroidMolecules(std::string filename, int targetFrame,
  int atomTypeI, bool isSlice, std::array<double, 3> coordLow, std::array<double, 3> coordHigh) {
  //
  std::vector<std::vector<double>> centroidCoord; // Coordinates of the centroid of the molecules
  molSys::PointCloud<molSys::Point<double>, double> yCloud; // PointCloud of all atoms 
  std::vector<double> currentCoord(3,0.0); // Coordinates of the current COM
  int numAtomsMol; // Number of atoms in the molecule 
  std::vector<bool> calcCentroidFlag; // Flag to make sure COM is not calculated twice 
  int iatomMolID; // Molecule ID of current molecule 
  int jatomID, jatomIndex; // Atoms belonging to the molecule iatomMolID
  std::unordered_multimap<int, int>
      molIDAtomIDmap; // atom IDs as keys and mol IDs as values

  // Get the pointCloud for all atoms in the system 
  yCloud = sinp::readLammpsTrj(filename, targetFrame, &yCloud, isSlice,coordLow,coordHigh);
  // Flags for calculating the COM for every atom
  // Set to true if the centroid has been calculated 
  calcCentroidFlag.resize(yCloud.nop, false);
  
  // Get the unordered map of the atom IDs (keys) and the molecular IDs
  // (values)
  molIDAtomIDmap = molSys::createMolIDAtomIDMultiMap(&yCloud);

  // Loop through all the atoms of atom type atomTypeI 
  // Calculate the centroid for those molecules 
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // Skip if the atom is not of type or if COM flag is true 
    if (yCloud.pts[iatom].type != atomTypeI || calcCentroidFlag[iatom])
    {
      continue;
    } // skip for other molecules or if COM is true 
    // --------
    // Find the COM for the molecule containing iatom
    calcCentroidFlag[iatom] = true; // set flag to true for COM calculation 
    // Find mol ID and atoms with that particular molecular ID 
    iatomMolID = yCloud.pts[iatom].molID; // molecule ID
    // Set currentCoord to iatom coordinates 
    currentCoord[0] = yCloud.pts[iatom].x; // x 
    currentCoord[1] = yCloud.pts[iatom].y; // y 
    currentCoord[2] = yCloud.pts[iatom].z; // z 
    numAtomsMol = 1;
    // --------
    // Find all atoms with iatomMolID and calculate the COM
    auto range = molIDAtomIDmap.equal_range(iatomMolID); 
    // Loop through all atoms with iatomMolID
    for (auto it = range.first; it != range.second; it++)
    {
      // it->second gives the value (in this case, the atom ID)
      jatomID = it->second; // Atom ID with molecule ID equal to iatomMolID
      auto gotJ = yCloud.idIndexMap.find(jatomID);
      jatomIndex = gotJ->second;
      // If jatomIndex has already been added, skip it
      if (calcCentroidFlag[jatomIndex])
      {
        continue;
      }
      // Add jatom coordinates to centroid
      currentCoord[0] += yCloud.pts[jatomIndex].x; // x 
      currentCoord[1] += yCloud.pts[jatomIndex].y; // y 
      currentCoord[2] += yCloud.pts[jatomIndex].z; // z 
      // Add to the number of atoms in the molecule 
      numAtomsMol++; 
      // Set the jatom calcCentroidFlag bool to true
      calcCentroidFlag[jatomIndex] = true; // jatomIndex COM has been calculated 
    } // end of loop with all atoms with iatomMolID
    // Divide the COM coordinates by numAtomsMol
    for (int k = 0; k < 3; ++k)
    {
      currentCoord[k] /= numAtomsMol; 
    } // divide by numAtomsMol
    // TODO: Do some error handling to avoid dividing by 0
    // Add to the vector of vectors
    centroidCoord.push_back(currentCoord); 
    // -----------
  } // end of loop through all atoms in yCloud 

  // Return the coordinates of the centers of masses of the molecules
  // with the same molecule ID, for the given atom type 
  return centroidCoord;
}

/**
 * @details Function for finding the k closest points (of type atomType) in a pointCloud, 
 * from a given target point (x y z coordinates).
 * @param[in] yCloud The given PointCloud
 * @param[in] molID The molecule ID for which atom index values will be returned
 * @param[in] maxCutoff Maximum cutoff distance in which to calculate the k closest points 
 * @return A vector of k atom indices in yCloud corresponding to the k closest points 
 */
std::vector<int> misc::kClosestPoints(molSys::PointCloud<molSys::Point<double>, double> yCloud, 
  int atomType, std::vector<double> targetPointCoord, int k, double maxCutoff) {
  //
  std::vector<int> pntIndices; // Vector of the closest points indices
  double dist; // Distance of iatom from the target point 
  std::vector< std::pair<double, int> > dIndVec; // Vector of distances and indices   

  // ----------------
  // Loop through all atoms in yCloud 
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // Skip for atoms which are not of type atomType 
    if (yCloud.pts[iatom].type != atomType)
    {
      continue;
    } // exclude atoms which are not of atomType
    // ----
    // Get the unwrapped distance from the target point 
    dist = gen::unWrappedDistFromPoint(&yCloud, iatom, targetPointCoord);
    //
    // Skip this if the distance is larger than maxCutoff
    if (dist > maxCutoff)
    {
      continue;
    } // skip for dist>maxCutoff
    //
    // ------
    // Otherwise, add to the vector of distances and indices 
    dIndVec.push_back(
      std::make_pair(dist, iatom)
      );
  } // end of loop through iatom
  // ----------------

  // Sort through the dIndVec using the distance 
  std::sort(dIndVec.begin(), dIndVec.end());

  // Error handling: if k is greater than dIndVec then something is 
  // wrong (TODO)

  // Loop through dIndVec and add the first k to pntIndices
  for (int i = 0; i < k; i++) {
    pntIndices.push_back(dIndVec[i].second);
  } // end of loop through first k points  

  // Return the indices of the k closest particles to the target point
  return pntIndices;
}