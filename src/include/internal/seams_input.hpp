//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#ifndef __SEAMS_INPUT_H_
#define __SEAMS_INPUT_H_

#include <iostream>
#include <memory>
#include <mol_sys.hpp>
#include <ring.hpp>
#include <string>

//// Boost
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"

/** @file seams_input.hpp
 *  @brief File for functions that read in files).
 */

/**
 *  @addtogroup sinp
 *  @{
 */

/** @brief Functions for the d-SEAMS readers.
 *  @details This namespace contains functions that are used for reading in the
 *   formats supported by d-SEAMS. LAMMPS trajectory files and XYZ files are
 *   currently supported, though it is recommended to use LAMMPS trajectory
 * files, since the simulation box size etc. are not normally present in XYZ
 * files, and many analyses may fail without the correct box dimensions.
 *
 * ### Changelog ###
 *
 * - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Dec 26, 2019
 * - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

namespace sinp {

//! Get file list inside the input folder
std::vector<std::string> getInpFileList(std::string inputFolder);

//! Function for reading in a specified frame (frame number and not timestep
//! value)
molSys::PointCloud<molSys::Point<double>, double>
readLammpsTrj(std::string filename, int targetFrame,
              bool isSlice = false,
              std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
              std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

//! Function for reading in a specified frame (frame number and not timestep
//! value) / This only reads in oxygen atoms
molSys::PointCloud<molSys::Point<double>, double> readLammpsTrjO(
    std::string filename, int targetFrame, int typeO,
    bool isSlice = false,
    std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
    std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

//! Function that reads in only atoms pf the desired type and ignores all atoms
//! which are not in the slice as well
molSys::PointCloud<molSys::Point<double>, double> readLammpsTrjreduced(
    std::string filename, int targetFrame, int typeI,
    bool isSlice = false,
    std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
    std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

//! Function for reading in atom coordinates from an XYZ file
molSys::PointCloud<molSys::Point<double>, double> readXYZ(std::string filename);

//! Reads bonds into a vector of vectors from a file with a specific format
std::vector<std::vector<int>> readBonds(std::string filename);

inline bool atomInSlice(double x, double y, double z,
                        std::array<double, 3> coordLow,
                        std::array<double, 3> coordHigh) {
  int flag = 0; //! If this is 3 then the particle is inside the volume slice

  if (((x >= coordLow[0]) && (x <= coordHigh[0])) ||
      coordLow[0] == coordHigh[0]) {
    flag++;
  }
  if (((y >= coordLow[1]) && (y <= coordHigh[1])) ||
      coordLow[1] == coordHigh[1]) {
    flag++;
  }
  if (((z >= coordLow[2]) && (z <= coordHigh[2])) ||
      coordLow[2] == coordHigh[2]) {
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
