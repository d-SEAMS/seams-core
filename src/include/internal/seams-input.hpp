#ifndef __SEAMS_INPUT_H_
#define __SEAMS_INPUT_H_

#include <iostream>
#include <memory>
#include <mol_sys.hpp>

//// Boost
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"

namespace sinp {

//// Get file list inside the input folder
std::vector<std::string> getInpFileList(std::string inputFolder);

//// Function for reading in a specified frame (frame number and not timestep value)
molSys::PointCloud<molSys::Point<double>, double>
readLammpsTrj(std::string filename, int targetFrame,
              molSys::PointCloud<molSys::Point<double>, double> *yCloud,
              bool isSlice = false,
              std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
              std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

//// Function for reading in a specified frame (frame number and not timestep value)
//// This only reads in oxygen atoms
molSys::PointCloud<molSys::Point<double>, double> readLammpsTrjO(
    std::string filename, int targetFrame,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int typeO,
    bool isSlice = false,
    std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
    std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

//// Function for reading in atom coordinates from an XYZ file
int readXYZ(std::string filename,
            molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//// Function for clearing vectors in PointCloud after multiple usage
molSys::PointCloud<molSys::Point<double>, double>
clearPointCloud(molSys::PointCloud<molSys::Point<double>, double> *yCloud);

inline bool atomInSlice(double x, double y, double z,
                        std::array<double, 3> coordLow,
                        std::array<double, 3> coordHigh) {
  int flag = 0; // If this is 3 then the particle is inside the volume slice

  if (x >= coordLow[0] && x <= coordHigh[0]) {
    flag++;
  }
  if (y >= coordLow[1] && y <= coordHigh[1]) {
    flag++;
  }
  if (z >= coordLow[2] && z <= coordHigh[2]) {
    flag++;
  }

  if (flag == 3) {
    return true;
  } else {
    return false;
  }
}

} // namespace sinp

#endif //// __SEAMS_INPUT_H_
