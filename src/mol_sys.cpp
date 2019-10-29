#include <iostream>
#include <memory>
#include <mol_sys.hpp>

/********************************************/ /**
 *  Function for clearing PointCloud if it is already
 filled. This should be called before every frame is read in.
 *  @param[out] yCloud The cleared PointCloud
 ***********************************************/
molSys::PointCloud<molSys::Point<double>, double> molSys::clearPointCloud(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  //
  std::vector<molSys::Point<double>> tempPts;
  std::vector<double> tempBox;
  //
  std::vector<double> tempBox1;

  tempPts.swap(yCloud->pts);
  tempBox.swap(yCloud->box);
  tempBox1.swap(yCloud->boxLow);
  yCloud->idIndexMap.clear();

  return *yCloud;
}

/********************************************/ /**
                                                *  Function for creating an
                                                *unordered map with the atomIDs
                                                *in the pointCloud as the keys
                                                *and the molecular IDs as the
                                                *values
                                                ***********************************************/
std::unordered_map<int, int> molSys::createIDMolIDmap(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  std::unordered_map<int, int>
      idMolIDmap;  // atom IDs as keys and mol IDs as values
  int iatomMolID;  // molID of the current iatom
  int iatomID;     // atom ID of the current iatom

  // Loop through the atoms in yCloud
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    iatomID = yCloud->pts[iatom].atomID;    // atom ID
    iatomMolID = yCloud->pts[iatom].molID;  // molecular ID
    // Update the unordered map
    idMolIDmap[iatomID] = iatomMolID;
  }  // end of loop through every iatom in pointCloud

  return idMolIDmap;
}
