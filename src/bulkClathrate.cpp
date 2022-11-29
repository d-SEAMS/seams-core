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
 *  @param[in] rcutoff Distance cutoff for building neighbour lists for the target cages. 
 *   (3.5 Angstrom is generally a good value)
 */
void clath::shapeMatchS2ClathrateSystem(std::string path, molSys::PointCloud<molSys::Point<double>, double> yCloud, 
  std::string filename, int targetFrame,
  int atomTypeI, bool isSlice, std::array<double, 3> coordLow, std::array<double, 3> coordHigh,
  std::string templateFileO, int oxygenAtomType, double rcutoff){
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
  Eigen::MatrixXdRowMajor targetPnts(nOxy, dim); // Target point set of just O atoms (Eigen matrix)
  //
  std::vector<int> atomIndices; // List of atom indices in yCloud 
  std::vector<std::vector<int>>
      rings;    // Vector of vectors of rings of 4 6-membered rings and 4 other atoms 
  std::vector<std::vector<double>> centroidCoord; // Coordinates of the centroid of the molecules
  double maxCutoff = 30.0; // in Angstroms; TODO change 
  // Shape-matching 
  std::vector<double> quaternionRot;         // quaternion rotation
  double rmsd;                               // least RMSD
  std::vector<double> rmsdList;              // List of RMSD per atom
  double scale;                              // Scale factor
  // Vector for the RMSD per atom and the RMSD per ring:
  std::vector<double> rmsdPerAtom;
  std::vector<int> noOfCommonElements; // An atom can belong to more than one ring or shape
  std::vector<cage::iceType>
      atomTypes; // This vector will have an ice type for every atom in yCloud 
  // -----------------------------------
  // INIT 
  // Init the atom type vector
  atomTypes.resize(yCloud.nop); // Has a value for each atom
  // Init the rmsd per atom
  rmsdPerAtom.resize(yCloud.nop, -1);    // Has a value for each atom
  noOfCommonElements.resize(yCloud.nop); // Has a value for each atom
  // -----------------------------------
  // Build the reference point set
  // This is row-ordered 
  std::tie(refCloud, ringsRef, refPntsO) = clath::buildRefS2Cage(templateFileO, oxygenAtomType);
  // -----------------------------------
  // Get the centroid of the encapsulated molecules in the region 
  // where you want to do shape-matching 
  centroidCoord = misc::getCentroidMolecules(filename, targetFrame, atomTypeI, isSlice, coordLow, coordHigh);

  // Loop through the THF centroids
  for (auto& centroidPnt : centroidCoord){
    // Find the 28 closest water molecules
    // output it into a vector of atom indices
    // and save it in a pointCloud
    std::tie(atomIndices,targetCloud) = misc::kClosestPoints(yCloud, oxygenAtomType, 
      centroidPnt, nOxy, maxCutoff);
    // The shape-matching with the S2 clathrate reference structure
    // is achieved by generating a vector of vectors corresponding to 
    // 4 6-membered rings with no elements in common, and 4 other 
    // O atoms corresponding to nexus points connecting the 6-membered 
    // rings via 5-membered rings.
    // These ordered arrangements are obtained by shape-matching individual
    // rings, added up together tp get the closest permutation 
    // TODO: periodic distance check??
    rings = clath::getShapeMatchedRingsInCage(targetCloud, refCloud, ringsRef, 
      rcutoff, oxygenAtomType);
    // Now get a row-ordered Eigen matrix of points (n x 3)
    // from the rings vector of vectors
    targetPnts = pntToPnt::fillTargetEigenPointSetFromRings(targetCloud, rings, nOxy, dim); 

    // Shape-matching with the reference point set 
    absor::hornAbsOrientationRowMajor(refPntsO, targetPnts, &quaternionRot, &rmsd,
                              &rmsdList, &scale);
    // Compare rmsd with some parametrized value (TODO: find that reference)
    // If the cage qualifies as a clathrate cage, change the atomTypes vector
    if (rmsd<10.0)
     {
       clath::assignAtomTypes(atomIndices, &atomTypes, cage::clathS2);
       // Reorder the atomIndices vector according to the order in the 
       // rings vector of vectors 
       // Assign RMSD per atom values and update noOfCommonElements
     } // the cage qualifies as a clathrate cage 

    // Put the bonds generated by the rings into an unordered multimap
  } // loop through the THF centroid points 

  // Write out the THF centroid points as single points 
  // Write out the oxygen atoms with the classification 
  //
  // Write out the oxygen atoms and THF centroid points 
  // and the RMSD per atom to a LAMMPS trajectory file
  int j=0; 

  return;
} // end of function

/**
 * @details Builds a vector of vectors with atom indices, corresponding to the arrangements
 * of 28 water molecules. 
 * This algorithm takes advantage of the connectivity information encoded in the 
 * constituent 6-membered and 5-membered rings in the candidate cage, thereby significantly reducing 
 * the number of possible point-to-point correspondence arrangements.
 * A typical $5^{12} 6^4$ cage, has 4 non-overlapping 6-membered rings and 4 water molecules which
 *  span between the rings, only participating in 5-membered rings.
 * 1. In the candidate $5^{12} 6^4$ cage, first find all the 6-membered primitive rings. If 4 6-membered rings are present, then proceed to the next steps. 
 * 2. The first ring should be matched with the first ring of the reference cage, which corresponds to finding the least RMSD of 12 arrangements. 
 * 3. Once the optimal point-to-point correspondence of the first ring has been obtained, add the second ring to the list of indices, and change the order of the last ring until the optimal arrangement has been found (12 arrangements). 
 * 4. Repeat this process for the third and fourth rings (24 arrangements). 
 * 5. The 4 remaining water molecules which are not part of the 6-membered rings can be tested with all possible permutations (24 arrangements). 

 * Therefore, this entire procedure only requires 72 arrangements to be tested, instead of $28!$ or $3E+29$ permutations. 
 * TODO: This can possibly be modified if $3$ or $2$ $6$-membered rings are found instead of $4$, but this would increase the number of arrangements to test. 
 */
std::vector<std::vector<int>>
clath::getShapeMatchedRingsInCage(molSys::PointCloud<molSys::Point<double>, double> targetCloud, molSys::PointCloud<molSys::Point<double>, double> refCloud, 
  std::vector<std::vector<int>> ringsRef, double rcutoff, int oxygenAtomType) {
  //
  std::vector<std::vector<int>> rings; // Output vector of vectors
  // which will contain the indices of targetCloud matched in the best arrangement
  // Target cage stuff
  std::vector<std::vector<int>> nList;  // Neighbour list for targetCloud
  std::vector<std::vector<int>> ringsHex, ringsAll;    // Vector of vectors of 
  // 28 atom indices in targetCloud
  std::vector<int> targetAtomIndices; // atom indices in targetCloud (0-27)
  // Vectors and variables for ring shape-matching 
  std::vector<int> tempRing, revRing; // 6-membered rings 
  std::vector<int> currentLastRing; // Current last ring for targetCloud
  std::vector<int> lastTargetRing; // The last ring for targetCloud with 4 elements 
  // RMSD stuff 
  std::vector<double> quaternionRot;         // quaternion rotation
  double rmsd1, rmsd2;                       // least RMSD
  std::vector<double> rmsdList1, rmsdList2;  // List of RMSD per atom
  double scale;                              // Scale factor
  // Variables for looping through possible permutations
  std::vector<double> currentQuat;      // quaternion rotation
  double currentRmsd;                   // least RMSD
  std::vector<double> currentRmsdList;  // List of RMSD per atom
  double currentScale;
  // --------------
  // Get 6-membered rings and neighbour list for the target cage 
  //
  // For the target point set 
  // Calculate a neighbour list
  nList = nneigh::neighListO(rcutoff, &targetCloud, oxygenAtomType);
  // Neighbour list by index
  nList = nneigh::neighbourListByIndex(&targetCloud, nList);
  // Find the vector of vector of rings
  ringsAll = primitive::ringNetwork(nList, 6); 
  // Get just the 6-membered rings 
  ringsHex = ring::getSingleRingSize(ringsAll, 6);
  // --------------
  // Shape-matching rings 
  //
  // Loop through the targetCloud rings 
  for (auto& iring : ringsHex)
  {
    tempRing = iring;
    revRing = iring;
    // Reverse the ring
    std::reverse(revRing.begin(), revRing.end());
    // No reversal 
    for (int i = 0; i < ringsHex[0].size(); ++i)
    {
      if (i==0)
      {
        // -------------------
        clath::matchClathrateLastRing(rings, tempRing, ringsRef, 
            targetCloud, refCloud,&currentQuat,
            &currentRmsd, &currentRmsdList, &currentScale);
        // Update for the first time 
        quaternionRot = currentQuat;
        rmsd1 = currentRmsd;
        rmsdList1 = currentRmsdList;
        scale = currentScale;
        currentLastRing = tempRing;
        // -------------------
      } // first step 
      else
      {
        // Change the order of the ring 
        rotate(tempRing.begin(), tempRing.begin()+1, tempRing.end());
        // Shape-matching 
        clath::matchClathrateLastRing(rings, tempRing, ringsRef, 
            targetCloud, refCloud,&currentQuat,
            &currentRmsd, &currentRmsdList, &currentScale);
        // Update if currentRmsd is less than rmsd1
        if (currentRmsd < rmsd1)
        {
          quaternionRot = currentQuat;
          rmsd1 = currentRmsd;
          rmsdList1 = currentRmsdList;
          scale = currentScale;
          currentLastRing = tempRing;
        } // end of update 
      } // all other steps apart from the first 
    } // go through 12 arrangements of the ring 

    // Reversed ring  
    for (int i = 0; i < ringsHex[0].size(); ++i)
    {
      if (i==0)
      {
        // -------------------
        clath::matchClathrateLastRing(rings, revRing, ringsRef, 
            targetCloud, refCloud,&currentQuat,
            &currentRmsd, &currentRmsdList, &currentScale);
        if (currentRmsd < rmsd1)
        {
          quaternionRot = currentQuat;
          rmsd1 = currentRmsd;
          rmsdList1 = currentRmsdList;
          scale = currentScale;
          currentLastRing = revRing;
        } // end of update 
        // -------------------
      } // first step 
      else{
        // Change the order of the ring 
        rotate(revRing.begin(), revRing.begin()+1, revRing.end());
        // Shape-matching 
        clath::matchClathrateLastRing(rings, revRing, ringsRef, 
            targetCloud, refCloud,&currentQuat,
            &currentRmsd, &currentRmsdList, &currentScale);
        // Update if currentRmsd is less than rmsd1
        if (currentRmsd < rmsd1)
        {
          quaternionRot = currentQuat;
          rmsd1 = currentRmsd;
          rmsdList1 = currentRmsdList;
          scale = currentScale;
          currentLastRing = revRing;
        } // end of update 
      } // all other steps apart from the first 
    } // go through 12 arrangements of the ring 

    rings.push_back(currentLastRing);

  } // end of loop through ringsHex
  // --------------
  // Get the last 4 elements, corresponding to the 5th ring in ringsRef
  // ringsRef[4]
  // Get all the 28 indices of the atoms in targetCloud
  for (int iatom = 0; iatom < targetCloud.nop; iatom++)
  {
      targetAtomIndices.push_back(iatom);
  } // overkill!
  // Flattened rings vector, which will contain all the indices
  // in the rings vector for the targetCloud
  auto flattenedRingsVec = std::accumulate(rings.begin(), rings.end(), decltype(rings)::value_type{},
          [](auto &x, auto &y) {
      x.insert(x.end(), y.begin(), y.end());
      return x;
  });
  // Sort the flattened rings vector and targetAtomIndices
  std::sort(targetAtomIndices.begin(), targetAtomIndices.end());
  std::sort(flattenedRingsVec.begin(), flattenedRingsVec.end());
  // --------------
  // Elements not in common 
  std::set_symmetric_difference(
      targetAtomIndices.begin(), targetAtomIndices.end(),
      flattenedRingsVec.begin(), flattenedRingsVec.end(),
      std::back_inserter(lastTargetRing));

  // Get all possible permutations and do shape-matching
  // of the last ring 
  std::sort(lastTargetRing.begin(), lastTargetRing.end()); // should be in ascending order
  int count=0; // for looping through all permutations
  //
  // Go through all the permutations 
  do {
    // reordered lastTargetRing
    if (count==0)
    {
      // -------------------
      clath::matchClathrateLastRing(rings, lastTargetRing, ringsRef, 
          targetCloud, refCloud,&currentQuat,
          &currentRmsd, &currentRmsdList, &currentScale);
      // Update for the first time 
      quaternionRot = currentQuat;
      rmsd1 = currentRmsd;
      rmsdList1 = currentRmsdList;
      scale = currentScale;
      currentLastRing = lastTargetRing;
      // -------------------
    } // first permutation
    else
    {
      // Shape-matching 
      clath::matchClathrateLastRing(rings, lastTargetRing, ringsRef, 
          targetCloud, refCloud,&currentQuat,
          &currentRmsd, &currentRmsdList, &currentScale);
      // Update if currentRmsd is less than rmsd1
      if (currentRmsd < rmsd1)
      {
        quaternionRot = currentQuat;
        rmsd1 = currentRmsd;
        rmsdList1 = currentRmsdList;
        scale = currentScale;
        currentLastRing = lastTargetRing;
      } // end of update
    } // subsequent permutations after the first  

    count++; // update the number of times the loop has run 
  } while (std::next_permutation(lastTargetRing.begin(), lastTargetRing.end()));

  // Presumably currentLastRing is the best match 
  rings.push_back(currentLastRing);
  // --------------
  // Return the rings for targetCloud, corresponding to the
  // best match to the reference S2 cage structure. 
  return rings;
} 

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
 * @details Assigns a type (cage::iceType) to each atom, according to the
 * previous classification, producing a vector of indices. 
 * @param[in] atomIndices Vector of atom indices.
 * @param[in] clathType Enum ice type to assign
 * @param[in] atomTypes Structural type for each atom.
 */
void clath::assignAtomTypes(std::vector<int> atomIndices,
                               std::vector<cage::iceType> *atomTypes,
                               cage::iceType clathType) {
  //
  int iatomIndex; // current atom index 

  for (auto it = atomIndices.begin(); it < atomIndices.end(); ++it)
  {
    iatomIndex = *it; // current atom index
    // Set the atom type to the clathrate type  
    (*atomTypes)[iatomIndex] = clathType;
  } // end of loop through atomIndices 

  return;
}

/**
 * @details Update the calculated RMSD per ring using the RMSD values of each
 * cage, and also update the values in the noOfCommonRings vector, which will be
 * used for averaging the RMSD per atom depending on the number of cages that
 * share that particular ring.
 */
void misc::updateRMSDatom(std::vector<int> atomIndices, 
  std::vector<std::vector<int>> rings,
  std::vector<double> rmsdList, std::vector<double> *rmsdPerAtom,
  std::vector<int> *noOfCommonAtoms) {
  //
  std::vector<int> reorderAtomIndices; // atom indices reordered according to rings

  // Reorder the atom indices according to the arrangement in the rings 
  // vector of vectors, which is how the shape-matching is done 
  reorderAtomIndices = misc::reorderAtomIndices(atomIndices, rings);

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
 * @return Vector of atom indices in yCloud and a pointCloud of k atoms in yCloud corresponding to the k closest points 
 */
std::pair<std::vector<int>, molSys::PointCloud<molSys::Point<double>, double>> misc::kClosestPoints(molSys::PointCloud<molSys::Point<double>, double> yCloud, 
  int atomType, std::vector<double> targetPointCoord, int k, double maxCutoff) {
  //
  molSys::PointCloud<molSys::Point<double>, double>
        outCloud;  // pointCloud for the closest points
  double dist; // Distance of iatom from the target point 
  std::vector< std::pair<double, int> > dIndVec; // Vector of distances and indices  
  std::vector<int> atomIndices; // atom indices according to yCloud 
  int iatom; // Atom index in yCloud  

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
    iatom = dIndVec[i].second;
    // Add to the atom index vector (according to yCloud)
    atomIndices.push_back(iatom);
    // Add to the output pointCloud 
    outCloud.pts.push_back(yCloud.pts[iatom]);
    outCloud.idIndexMap[yCloud.pts[iatom].atomID] = outCloud.pts.size() - 1;
  } // end of loop through first k points  

  // Set the number of particles and the current frame number
  outCloud.nop = outCloud.pts.size(); // This should be the same as k 
  // TODO: error handling here.  
  outCloud.currentFrame = yCloud.currentFrame; // Current frame number
  outCloud.box = yCloud.box; // Set the box lengths so that the periodic distance works  

  // Return the indices of the k closest particles to the target point
  return std::make_pair(atomIndices, outCloud);
}

/**
 * @details Function for reordering atom indices (corresponding to atom indices
 * in a larger pointCloud, according to a vector of vectors of the rings (according 
 * to a smaller pointCloud) 
 */
std::vector<int> misc::reorderAtomIndices(std::vector<int> atomIndices, 
    std::vector<std::vector<int>> rings) {
  // 
  std::vector<int> outAtomIndices; // atom indices according to the larger pointCloud 
  int targetCloudIndex; // Index according to the smaller targetCloud
  int yCloudIndex; // Index according to the larger yCloud  

  // Loop through the rings vector of vectors
  // Here the atom indices correspond to the smaller pointCloud 
  for (auto& iring : rings){
    // Loop through the current iring 
    for (auto it = iring.begin(); it < iring.end(); ++it)
    {
      targetCloudIndex = *it; // index in the rings
      yCloudIndex = atomIndices[targetCloudIndex];
      // Update the output vector 
      outAtomIndices.push_back(yCloudIndex);  
    } // end of loop through iring 
  } // end of loop through rings vector of vectors 

  // Return the atom indices ordered according 
  // to the arrangement in the clathrate cage, 
  // corresponding to the atom indices in the larger pointCloud 
  return outAtomIndices;
}