#include <bulkTUM.hpp>

// -----------------------------------------------------------------------------------------------------
// TOPOLOGICAL UNIT MATCHING ALGORITHMS
// -----------------------------------------------------------------------------------------------------

/********************************************/ /**
 *  Finds out if rings constitute double-diamond cages or hexagonal cages.
 Requires a neighbour list (by index) and a vector of vectors of primitive rings
 which should also be by index. This is registered as a Lua function and is
 accessible to the user.
 Internally, this function calls the following functions:
 - ring::getSingleRingSize (Saves rings of a single ring size into a new vector
 of vectors, which is subsequently used for finding DDCs, HCs etc).
 - ring::findDDC (Finds the DDCs).
 - ring::findHC (Finds the HCs).
 - ring::findMixedRings (Finds the mixed rings, which are shared by DDCs and HCs
 both).
 - ring::getStrucNumbers (Gets the number of structures (DDCs, HCs, mixed rings,
 basal rings, prismatic rings, to be used for write-outs).
 - sout::writeTopoBulkData (Writes out the numbers and data obtained for the
 current frame).
 - ring::getAtomTypesTopoBulk (Gets the atom type for every atom, to be used for
 printing out the ice types found).
 - sout::writeLAMMPSdataTopoBulk (Writes out the atoms, with the classified
 types, into a LAMMPS data file, which can be visualized in OVITO).
 *  @param[in] path The file path of the output directory to which output files
 will be written.
 *  @param[in] rings Vector of vectors containing the primitive rings. This
 contains rings of all sizes.
 *  @param[in] nList Row-ordered neighbour list, by index.
 *  @param[in] yCloud The input PointCloud, with respect to which the indices in
 the rings and nList vector of vectors have been saved.
 *  @param[in] firstFrame First frame to be analyzed
 *  @param[in] onlyTetrahedral Flag for only finding DDCs and HCs (true) or also
 finding PNCs (false)
 ***********************************************/
int tum3::topoUnitMatchingBulk(
    std::string path, std::vector<std::vector<int>> rings,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int firstFrame,
    bool printClusters, bool onlyTetrahedral) {
  //
  // The input rings vector has rings of all sizes
  //
  // ringType has a value for rings of a particular size
  std::vector<ring::strucType>
      ringType;  // This vector will have a value for each ring inside
                 // ringList, of type enum strucType in gen.hpp
  // Make a list of all the DDCs and HCs
  std::vector<cage::Cage> cageList;
  std::vector<std::vector<int>>
      ringsOneType;     // Vector of vectors of rings of a single size
  int initRingSize;     // Todo or not: calculate the PNCs or not
  int maxRingSize = 6;  // DDCs and HCs are for 6-membered rings
  std::vector<cage::iceType>
      atomTypes;  // This vector will have a value for every atom
  // Number of types
  int numHC, numDDC, mixedRings, prismaticRings, basalRings;
  // Shape-matching variables ----
  double rmsd;               // RMSD value for a particular cage type
  std::vector<double> quat;  // Quaternion obtained from shape-matching
  std::vector<std::vector<double>> quatList;  // List of quaternions
  // Reference points
  Eigen::MatrixXd refPntsDDC(
      14, 3);  // Reference point set (Eigen matrix) for a DDC
  Eigen::MatrixXd refPntsHC(12,
                            3);  // Reference point set (Eigen matrix) for a DDC

  // Vector for the RMSD per atom and the RMSD per ring:
  std::vector<double> rmsdPerAtom, rmsdPerRing;
  std::vector<int>
      noOfCommonRings;  // An atom can belong to more than one ring or shape
  //

  if (onlyTetrahedral) {
    initRingSize = 6;
  } else {
    initRingSize = 5;
  }

  // Init the atom type vector
  atomTypes.resize(yCloud->nop);  // Has a value for each atom
  // Init the rmsd per atom
  rmsdPerAtom.resize(yCloud->nop, -1);  // Has a value for each atom

  // ----------------------------------------------
  // Init
  // Get rings of size 5 or 6.
  for (int ringSize = initRingSize; ringSize <= maxRingSize; ringSize++) {
    // Clear ringsOneType
    ring::clearRingList(ringsOneType);
    // Get rings of the current ring size
    ringsOneType = ring::getSingleRingSize(rings, ringSize);
    // Skip for zero rings
    if (ringsOneType.size() == 0) {
      continue;
    }  // skip for no rings of ringSize
    //
    // Init the ringType vector
    ringType.resize(
        ringsOneType.size());  // Has a value for each ring. init to zero.
    // ----------------------------------------------
    if (ringSize == 6) {
      // Get the cages

      // Find the DDCs and HCs
      cageList = tum3::topoBulkCriteria(path, ringsOneType, nList, yCloud,
                                        firstFrame, &numHC, &numDDC, &ringType);

      // Gets the atom type for every atom, to be used for printing out the ice
      // types found
      ring::getAtomTypesTopoBulk(ringsOneType, ringType, &atomTypes);
    }
    // Finding the 5-membered prismatic blocks
    else {
      // Get the prism block classifications
      prism3::findBulkPrisms(ringsOneType, &ringType, nList, yCloud,
                             &rmsdPerAtom);
      // Gets the atom type for every atom, to be used for printing out the ice
      // types found
      ring::getAtomTypesTopoBulk(ringsOneType, ringType, &atomTypes);
    }
  }

  // --------------------------------------------------
  // SHAPE-MATCHING
  //
  // Get the reference point sets
  //
  std::string filePathHC = "templates/hc.xyz";
  std::string filePathDDC = "templates/ddc.xyz";
  refPntsHC = tum3::buildRefHC(filePathHC);     // HC
  refPntsDDC = tum3::buildRefDDC(filePathDDC);  // DDC
  //
  // Init
  noOfCommonRings.resize(ringsOneType.size());  // Initialized to zero
  rmsdPerRing.resize(
      ringsOneType.size(),
      -1);  // Initialized to -1 which is the value for dummy atoms
  //
  // Loop through the entire cageList vector of cages to match the HCs and DDCs
  //
  // HCs are first, followed by DDCs.
  //
  // Go through all the HCs
  for (int icage = 0; icage < numHC; icage++) {
    // Match against a perfect HC
    tum3::shapeMatchHC(yCloud, refPntsHC, cageList[icage], ringsOneType, nList,
                       &quat, &rmsd);
    // Update the vector of quaternions
    quatList.push_back(quat);
    // Update the RMSD per ring
    tum3::updateRMSDring(cageList[icage], rmsd, &rmsdPerRing, &noOfCommonRings);
  }  // end of looping through all HCs
  // --------------------------------------------------
  // Go through all the DDCs
  for (int icage = numHC; icage < cageList.size(); icage++) {
    // Match against a perfect DDC
    tum3::shapeMatchDDC(yCloud, refPntsDDC, cageList, icage, ringsOneType,
                        &quat, &rmsd);
    // Update the vector of quaternions
    quatList.push_back(quat);
    // Update the RMSD per ring
    tum3::updateRMSDring(cageList[icage], rmsd, &rmsdPerRing, &noOfCommonRings);
  }  // end of looping through all HCs

  // --------------------------------------------------
  // Getting the RMSD per atom
  // Update the RMSD per atom
  tum3::updateRMSDatom(rmsdPerRing, noOfCommonRings, &rmsdPerAtom, atomTypes,
                       ringsOneType);
  // --------------------------------------------------

  // Print out the lammps data file with the bonds and types
  sout::writeLAMMPSdataTopoBulk(yCloud, nList, atomTypes, path);
  // To output the bonds between dummy atoms, uncomment the following line
  // sout::writeLAMMPSdataTopoBulk(yCloud, nList, atomTypes, path, true);

  // --------------------------------------------------
  // Writes out the dumpfile
  //
  std::vector<int> dumpAtomTypes;  // Must be typecast to int
  dumpAtomTypes.resize(atomTypes.size());
  // Change from enum to int
  for (int i = 0; i < atomTypes.size(); i++) {
    if (atomTypes[i] == cage::hc) {
      dumpAtomTypes[i] = 1;
    }  // HC
    else if (atomTypes[i] == cage::ddc) {
      dumpAtomTypes[i] = 2;
    }  // DDC
    else if (atomTypes[i] == cage::mixed) {
      dumpAtomTypes[i] = 3;
    }  // mixed rings
    else if (atomTypes[i] == cage::pnc) {
      dumpAtomTypes[i] = 4;
    }  // pnc
    else if (atomTypes[i] == cage::mixed2) {
      dumpAtomTypes[i] = 5;
    }  // shared by pnc and ddc/hc
    else {
      dumpAtomTypes[i] = 0;
    }  // dummy atoms
  }    // end of changing the type to int
  sout::writeLAMMPSdumpCages(yCloud, rmsdPerAtom, dumpAtomTypes, path,
                             firstFrame);
  // --------------------------------------------------
  return 0;
}

// Shape-matching for an HC
/********************************************/ /**
 *  Match the input cage with a perfect HC.
 The quaternion for the rotation and the RMSD is outputted.
 *  @param[in] refPoints The point set of the reference system (or right
system). This is a \f$ (n \times 3) \f$ Eigen matrix. Here, \f$ n \f$ is the
number of particles.
 *  @param[in] targetPoints \f$ (n \times 3) \f$ Eigen matrix of the
candidate/test system (or left system).
 *  @param[in, out] quat The quaternion for the optimum rotation of the left
system into the right system.
 *  @param[in, out] scale The scale factor obtained from Horn's algorithm.
 ***********************************************/
int tum3::shapeMatchHC(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    const Eigen::MatrixXd &refPoints, cage::Cage cageUnit,
    std::vector<std::vector<int>> rings, std::vector<std::vector<int>> nList,
    std::vector<double> *quat, double *rmsd) {
  //
  int iring,
      jring;  // Indices in the ring vector of vectors for the basal rings
  std::vector<int> basal1,
      basal2;        // Re-ordered basal rings matched to each other
  int ringSize = 6;  // Each ring has 6 nodes
  Eigen::MatrixXd targetPointSet(12, 3);  // Target point set (Eigen matrix)
  //
  std::vector<double> rmsdList;  // List of RMSD per atom
  // Variables for looping through possible permutations
  //
  std::vector<double> currentQuat;  // quaternion rotation
  double currentRmsd;               // least RMSD
  double currentScale;              // scale
  double scale;                     // Final scale
  int index;  // Int for describing which permutation matches the best

  // Init
  iring = cageUnit.rings[0];     // Index of basal1
  jring = cageUnit.rings[1];     // Index of basal2
  rmsdList.resize(yCloud->nop);  // Not actually updated here
  //
  // ----------------
  // Re-order the basal rings so that they are matched
  pntToPnt::relOrderHC(yCloud, rings[iring], rings[jring], nList, &basal1,
                       &basal2);
  // ----------------
  // Loop through all possible permutations
  //
  for (int i = 0; i < ringSize; i++) {
    // Change the order of the target points somehow!
    //
    targetPointSet = pntToPnt::changeHexCageOrder(yCloud, basal1, basal2, i);
    // Shape-matching
    absor::hornAbsOrientation(refPoints, targetPointSet, &currentQuat,
                              &currentRmsd, &rmsdList, &currentScale);
    if (i == 0) {
      *quat = currentQuat;
      *rmsd = currentRmsd;
      scale = currentScale;
      index = 0;
    } else {
      if (currentRmsd < *rmsd) {
        *quat = currentQuat;
        *rmsd = currentRmsd;
        scale = currentScale;
        index = i;
      }  // update
    }    // Update if this is a better match
  }      // Loop through possible permutations
  // ----------------
  return 0;
}

// Shape-matching for a DDC
/********************************************/ /**
 *  Match the input cage with a perfect DDC.
 The quaternion for the rotation and the RMSD is outputted.
 *  @param[in] refPoints The point set of the reference system (or right
system). This is a \f$ (n \times 3) \f$ Eigen matrix. Here, \f$ n \f$ is the
number of particles.
 *  @param[in] targetPoints \f$ (n \times 3) \f$ Eigen matrix of the
candidate/test system (or left system).
 *  @param[in, out] quat The quaternion for the optimum rotation of the left
system into the right system.
 *  @param[in, out] scale The scale factor obtained from Horn's algorithm.
 ***********************************************/
int tum3::shapeMatchDDC(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    const Eigen::MatrixXd &refPoints, std::vector<cage::Cage> cageList,
    int cageIndex, std::vector<std::vector<int>> rings,
    std::vector<double> *quat, double *rmsd) {
  //
  std::vector<int> ddcOrder;              // Connectivity of the DDC
  int ringSize = 6;                       // Each ring has 6 nodes
  Eigen::MatrixXd targetPointSet(14, 3);  // Target point set (Eigen matrix)
  //
  std::vector<double> rmsdList;  // List of RMSD per atom
  // Variables for looping through possible permutations
  //
  std::vector<double> currentQuat;  // quaternion rotation
  double currentRmsd;               // least RMSD
  double currentScale;              // scale
  double scale;                     // Final scale
  int index;  // Int for describing which permutation matches the best

  // ----------------
  // Save the order of the DDC in a vector
  ddcOrder = pntToPnt::relOrderDDC(cageIndex, rings, cageList);
  // ----------------
  // Loop through all possible permutations
  //
  for (int i = 0; i < ringSize; i++) {
    // Change the order of the target points somehow!
    //
    targetPointSet = pntToPnt::changeDiaCageOrder(yCloud, ddcOrder, i);
    // Shape-matching
    absor::hornAbsOrientation(refPoints, targetPointSet, &currentQuat,
                              &currentRmsd, &rmsdList, &currentScale);
    if (i == 0) {
      *quat = currentQuat;
      *rmsd = currentRmsd;
      scale = currentScale;
      index = 0;
    } else {
      if (currentRmsd < *rmsd) {
        *quat = currentQuat;
        *rmsd = currentRmsd;
        scale = currentScale;
        index = i;
      }  // update
    }    // Update if this is a better match
  }      // Loop through possible permutations
  // ----------------
  return 0;
}

/********************************************/ /**
 *  Build a reference Hexagonal cage,
 reading it in from a templates directory
 ***********************************************/
Eigen::MatrixXd tum3::buildRefHC(std::string fileName) {
  //
  Eigen::MatrixXd refPnts(12, 3);  // Reference point set (Eigen matrix)
  // Get the reference HC point set
  molSys::PointCloud<molSys::Point<double>, double>
      setCloud;  // PointCloud for holding the reference point values
  // Variables for rings
  std::vector<std::vector<int>> nList;  // Neighbour list
  std::vector<std::vector<int>> rings;  // Rings
  std::vector<ring::strucType>
      ringType;  // This vector will have a value for each ring inside
  std::vector<int> listHC;  // Contains atom indices of atoms making up HCs
  // Make a list of all the DDCs and HCs
  std::vector<cage::Cage> cageList;
  int iring, jring;
  //
  // read in the XYZ file into the pointCloud setCloud
  //
  sinp::readXYZ(fileName, &setCloud);
  // box lengths
  for (int i = 0; i < 3; i++) {
    setCloud.box[i] = 50;
  }  // end of setting box lengths
  //

  nList = nneigh::neighListO(3.5, &setCloud, 1);
  // Neighbour list by index
  nList = nneigh::neighbourListByIndex(&setCloud, nList);
  // Find the vector of vector of rings
  rings = primitive::ringNetwork(nList, 6);
  // init the ringType vector
  ringType.resize(rings.size());
  // Find the HCs
  listHC = ring::findHC(rings, &ringType, nList, &cageList);
  // Get the basal rings from cageList
  iring = cageList[0].rings[0];
  jring = cageList[0].rings[1];
  //
  std::vector<int> matchedBasal1,
      matchedBasal2;  // Re-ordered basal rings 1 and 2
  // Reordered basal rings
  // Getting the target Eigen vectors
  // Get the re-ordered matched basal rings, ordered with respect to each
  // other

  pntToPnt::relOrderHC(&setCloud, rings[iring], rings[jring], nList,
                       &matchedBasal1, &matchedBasal2);
  // Get the reference point set
  refPnts =
      pntToPnt::changeHexCageOrder(&setCloud, matchedBasal1, matchedBasal2, 0);
  //
  return refPnts;
}

/********************************************/ /**
 *  Build a reference Double-Diamond cage,
 reading it in from a templates directory
 ***********************************************/
Eigen::MatrixXd tum3::buildRefDDC(std::string fileName) {
  //
  Eigen::MatrixXd refPnts(14, 3);  // Reference point set (Eigen matrix)
  // Get the reference HC point set
  molSys::PointCloud<molSys::Point<double>, double>
      setCloud;  // PointCloud for holding the reference point values
  // Variables for rings
  std::vector<std::vector<int>> nList;  // Neighbour list
  std::vector<std::vector<int>> rings;  // Rings
  std::vector<ring::strucType>
      ringType;  // This vector will have a value for each ring inside
  std::vector<int> listDDC,
      listHC;  // Contains atom indices of atoms making up DDCs and HCs
  std::vector<int> ddcOrder;  // Atom indices of particles in the DDC
  // Make a list of all the DDCs and HCs
  std::vector<cage::Cage> cageList;
  int iring, jring;
  //
  // read in the XYZ file into the pointCloud setCloud
  //
  sinp::readXYZ(fileName, &setCloud);

  // box lengths
  for (int i = 0; i < 3; i++) {
    setCloud.box[i] = 50;
  }  // end of setting box lengths
  //

  nList = nneigh::neighListO(3.5, &setCloud, 1);
  // Neighbour list by index
  nList = nneigh::neighbourListByIndex(&setCloud, nList);
  // Find the vector of vector of rings
  rings = primitive::ringNetwork(nList, 6);
  // init the ringType vector
  ringType.resize(rings.size());
  // Find the DDCs
  listDDC = ring::findDDC(rings, &ringType, listHC, &cageList);
  // Save the order of the DDC in a vector
  ddcOrder = pntToPnt::relOrderDDC(0, rings, cageList);
  // Get the reference point set
  refPnts = pntToPnt::changeDiaCageOrder(&setCloud, ddcOrder, 0);
  //
  return refPnts;
}

/********************************************/ /**
 *  Update the calculated RMSD per ring using the RMSD values of each cage,
 and also update the values in the noOfCommonRings vector, which will be used
 for averaging the RMSD per atom depending on the number of cages that share
 that particular ring.
 ***********************************************/
int tum3::updateRMSDring(cage::Cage cageUnit, double rmsd,
                         std::vector<double> *rmsdPerRing,
                         std::vector<int> *noOfCommonRings) {
  //
  int nRings = cageUnit.rings.size();  // Number of rings in the current cage
  int iring;  // Index according to the rings vector of vector, for the current
              // ring inside the cage

  // Loop through the rings in each cage
  for (int i = 0; i < nRings; i++) {
    iring = cageUnit.rings[i];  // Current ring

    if ((*rmsdPerRing)[iring] == -1) {
      (*rmsdPerRing)[iring] = rmsd;
      (*noOfCommonRings)[iring] = 1;
    }  // initialization
    else {
      (*rmsdPerRing)[iring] += rmsd;
      (*noOfCommonRings)[iring] += 1;
    }  // add to the number of shared rings

  }  // end of looping through the rings in the cage

  return 0;
}

/********************************************/ /**
 *  Update the calculated RMSD per atom using the RMSD values of each ring,
 averaged by, the values in the noOfCommonRings vector.
 ***********************************************/
int tum3::updateRMSDatom(std::vector<double> rmsdPerRing,
                         std::vector<int> noOfCommonRings,
                         std::vector<double> *rmsdPerAtom,
                         std::vector<cage::iceType> atomTypes,
                         std::vector<std::vector<int>> rings) {
  //
  int iatom;                       // Atom index in iring
  int ringSize = rings[0].size();  // This should be 6
  double rmsdVal;                  // Current RMSD value
  //
  // Loop through the calculated RMSD
  for (int iring = 0; iring < rings.size(); iring++) {
    //
    if (noOfCommonRings[iring] == 0) {
      rmsdVal = -1;
    }  // for a dummy atom
    else {
      rmsdVal = rmsdPerRing[iring] / noOfCommonRings[iring];
    }  // Get the averaged rmsd Value
    //
    // Update the atoms
    // Loop through every atom
    for (int i = 0; i < ringSize; i++) {
      iatom = rings[iring][i];  // Atom index
      // Do not update for PNC atoms
      if (atomTypes[iatom] == cage::pnc || atomTypes[iatom] == cage::mixed2) {
        continue;
      }  // do not update for pnc
      else {
        (*rmsdPerAtom)[iatom] = rmsdVal;
      }
    }  // end of looping through ring atoms
  }    // Loop through every ring

  return 0;
}

// -----------------------------------------------------------------
//
// TOPOLOGICAL NETWORK ALGORITHMS
//
// -----------------------------------------------------------------
/********************************************/ /**
 *  Finds out if rings constitute double-diamond cages or hexagonal cages.
 Requires a neighbour list (by index) and a vector of vectors of primitive rings
 which should also be by index. This is only for rings of size 6!
 Internally, this function calls the following functions:
 - ring::getSingleRingSize (Saves rings of a single ring size into a new vector
 of vectors, which is subsequently used for finding DDCs, HCs etc).
 - ring::findDDC (Finds the DDCs).
 - ring::findHC (Finds the HCs).
 - ring::findMixedRings (Finds the mixed rings, which are shared by DDCs and HCs
 both).
 - ring::getStrucNumbers (Gets the number of structures (DDCs, HCs, mixed rings,
 basal rings, prismatic rings, to be used for write-outs).
 - sout::writeTopoBulkData (Writes out the numbers and data obtained for the
 current frame).
 - ring::getAtomTypesTopoBulk (Gets the atom type for every atom, to be used for
 printing out the ice types found).
 - sout::writeLAMMPSdataTopoBulk (Writes out the atoms, with the classified
 types, into a LAMMPS data file, which can be visualized in OVITO).
 *  @param[in] path The file path of the output directory to which output files
 will be written.
 *  @param[in] rings Vector of vectors containing the primitive rings. This
 should contain rings of only size 6.
 *  @param[in] nList Row-ordered neighbour list, by index.
 *  @param[in] yCloud The input PointCloud, with respect to which the indices in
 the rings and nList vector of vectors have been saved.
 *  @param[in] firstFrame First frame to be analyzed
 *  @param[in] onlyTetrahedral Flag for only finding DDCs and HCs (true) or also
 finding PNCs (false)
 ***********************************************/
std::vector<cage::Cage> tum3::topoBulkCriteria(
    std::string path, std::vector<std::vector<int>> rings,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int firstFrame,
    int *numHC, int *numDDC, std::vector<ring::strucType> *ringType) {
  //
  // Ring IDs of each type will be saved in these vectors
  std::vector<int> listDDC;  // Vector for ring indices of DDC
  std::vector<int> listHC;   // Vector for ring indices of HC
  std::vector<int>
      listMixed;  // Vector for ring indices of rings that are both DDC and HC
                  // ringList, of type enum strucType in gen.hpp
  // Make a list of all the DDCs and HCs
  std::vector<cage::Cage> cageList;
  // Number of types
  int mixedRings, prismaticRings, basalRings;

  // ----------------------------------------------
  // Init
  //
  *numHC = 0;   // Number of hexagonal cages
  *numDDC = 0;  // Init the number of DDCs
  // Quit the function for zero rings
  if (rings.size() == 0) {
    return cageList;
  }  // skip for no rings of ringSize
  //
  // Init the ringType vector
  (*ringType).resize(rings.size());  // Has a value for each ring. init to zero.
  // ----------------------------------------------
  // Get the cages

  // Find HC rings, saving the ring IDs (starting from 0) to listHC
  listHC = ring::findHC(rings, ringType, nList, &cageList);

  // Find DDC rings, saving the IDs to listDDC
  listDDC = ring::findDDC(rings, ringType, listHC, &cageList);

  // Find rings which are both DDCs and HCs (mixed)
  // A dummy value of -10 in the listDDC and listHC vectors for mixed rings
  listMixed = ring::findMixedRings(rings, ringType, &listDDC, &listHC);

  // Get the number of structures (DDCs, HCs, mixed rings, basal rings,
  // prismatic rings)
  ring::getStrucNumbers(*ringType, cageList, numHC, numDDC, &mixedRings,
                        &prismaticRings, &basalRings);

  // Write out to a file
  sout::writeTopoBulkData(path, yCloud->currentFrame, *numHC, *numDDC,
                          mixedRings, basalRings, prismaticRings, firstFrame);

  return cageList;
}