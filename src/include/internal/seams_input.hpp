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

#ifndef SEAMS_SEAMS_INPUT_H_
#define SEAMS_SEAMS_INPUT_H_

#include <functional>
#include <iostream>
#include <memory>
#include <mol_sys.hpp>
#include <ring.hpp>
#include <string>

#include <filesystem>

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

//! Number of ITEM: TIMESTEP frames in a LAMMPS dump, or 0 if the file
//! cannot be read. The first full count per path walks the file once
//! and later reads seek, matching chemfiles Trajectory::nsteps and
//! readcon's frame-offset table. Sequential load_frame walks reuse a
//! live cursor the way LAMMPS ReaderNative does on rerun.
int nLammpsFrames(const std::string &filename);

//! Drop a cached dump session (tests that rewrite a path in place).
void dropLammpsDumpIndex(const std::string &filename);

//! Call fn(frame, cloud) for each frame in [first, last] (1-based,
//! inclusive). last <= 0 means every ITEM: TIMESTEP. typeFilter <= 0
//! keeps every atom. nThreads <= 0 uses the OpenMP default; 1 is
//! serial. Each worker opens its own handle and seeks the shared
//! offset table. Incremental RingUpdater / AffiliationUpdater state
//! cannot be shared across workers: use the batch classifiers.
void forEachLammpsFrame(
    const std::string &filename, int first, int last, int typeFilter,
    const std::function<void(
        int, molSys::PointCloud<molSys::Point<double>, double> &)> &fn,
    int nThreads = 0);

//! Function for reading in a specified frame (frame number and not timestep
//! value)
molSys::PointCloud<molSys::Point<double>, double>
readLammpsTrj(std::string filename, int targetFrame,
              molSys::PointCloud<molSys::Point<double>, double> &yCloud,
              bool isSlice = false,
              std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
              std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

//! One LAMMPS dump frame, one atom type. The type argument is any
//! LAMMPS type (the O in the name is historical). If isSlice is true,
//! each kept point gets inSlice set; atoms outside the box stay in
//! the cloud. nop is the type-filtered count.
molSys::PointCloud<molSys::Point<double>, double> readLammpsTrjO(
    std::string filename, int targetFrame,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeO,
    bool isSlice = false,
    std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
    std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

//! One LAMMPS dump frame, one atom type, dropping atoms outside the
//! slice when isSlice is true. nop is the kept count. An axis with
//! lo == hi is unconstrained.
molSys::PointCloud<molSys::Point<double>, double> readLammpsTrjreduced(
    std::string filename, int targetFrame,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeI,
    bool isSlice = false,
    std::array<double, 3> coordLow = std::array<double, 3>{0, 0, 0},
    std::array<double, 3> coordHigh = std::array<double, 3>{0, 0, 0});

//! Function for reading in atom coordinates from an XYZ file
molSys::PointCloud<molSys::Point<double>, double> readXYZ(std::string filename);

//! Reads bonds into a vector of vectors from a file with a specific format
std::vector<std::vector<int>> readBonds(std::string filename);

#ifdef SEAMS_HAS_CHEMFILES
//! Read any trajectory format supported by chemfiles (PDB, GRO, DCD, etc.)
molSys::PointCloud<molSys::Point<double>, double>
readChemfiles(std::string filename, int targetFrame,
              molSys::PointCloud<molSys::Point<double>, double> &yCloud,
              int typeFilter = -1);
#endif

#ifdef SEAMS_HAS_READCON
//! Read a .con format file (eOn saddle point search trajectories)
molSys::PointCloud<molSys::Point<double>, double>
readCon(std::string filename, int targetFrame,
        molSys::PointCloud<molSys::Point<double>, double> &yCloud);
#endif

//! True when each component lies in [lo, hi], or that axis has lo == hi
//! (unconstrained). ([0,0,0], [50,0,0]) is x in [0, 50], y and z open.
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

#endif //// SEAMS_SEAMS_INPUT_H_
