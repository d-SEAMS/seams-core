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
 *  @param[in] firstFrame First frame to be analyzed
 *  @param[in] onlyTetrahedral Flag for only finding DDCs and HCs (true) or also
 *   finding PNCs (false)
 */
void clath::shapeMatchS2ClathrateSystem(std::string templateFileName, std::string templateFileO, int oxygenAtomType){
  //
  Eigen::MatrixXdRowMajor refPntsO(28, 3); // Reference point set of just O atoms (Eigen matrix)
  Eigen::MatrixXdRowMajor refPntsWat(84, 3); // Reference point set of O H H water atoms (Eigen matrix)
  // -----------------------------------
  // Build the reference point sets 
  std::tie(refPntsO, refPntsWat) = clath::buildRefS2CageLammpsTrj(templateFileName, templateFileO, oxygenAtomType);
  // -----------------------------------

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
