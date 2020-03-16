#include <bulkTUM.hpp>

// -----------------------------------------------------------------------------------------------------
// DDC / HC ALGORITHMS
// -----------------------------------------------------------------------------------------------------

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