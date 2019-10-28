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
