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

#ifndef __SEAMS_OUTPUT_H_
#define __SEAMS_OUTPUT_H_

#include <bond.hpp>
#include <cage.hpp>
#include <errno.h> // errno, ENOENT, EEXIST
#include <generic.hpp>
#include <iostream>
#include <memory>
#include <mol_sys.hpp>
#include <sys/stat.h> // stat
#if defined(_WIN32)
#include <direct.h> // _mkdir
#endif

// Boost
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
namespace fs = boost::filesystem;
// #include <filesystem>
// namespace fs = std::filesystem;

namespace sout {

/**
 *  Inline function for checking if the directory exists or not
 *  @param[in] path The path of the directory
 *  @return True or false
 */
inline bool isDirExist(const std::string &path) {
#if defined(_WIN32)
  struct _stat info;
  if (_stat(path.c_str(), &info) != 0) {
    return false;
  }
  return (info.st_mode & _S_IFDIR) != 0;
#else
  struct stat info;
  if (stat(path.c_str(), &info) != 0) {
    return false;
  }
  return (info.st_mode & S_IFDIR) != 0;
#endif
}

/**
 *  Inline function for creating the desried directory.
 *  @param[in] path The path of the directory
 */
inline int makePath(const std::string &path) {
#if defined(_WIN32)
  int ret = _mkdir(path.c_str());
#else
  mode_t mode = 0755;
  int ret = mkdir(path.c_str(), mode);
#endif
  if (ret == 0)
    return 0;

  switch (errno) {
  case ENOENT:
    // parent didn't exist, try to create it
    {
      int pos = path.find_last_of('/');
      if (pos == std::string::npos)
#if defined(_WIN32)
        pos = path.find_last_of('\\');
      if (pos == std::string::npos)
#endif
        return 1;
      if (!makePath(path.substr(0, pos)))
        return 1;
    }
// now, try to create again
#if defined(_WIN32)
    return 0 == _mkdir(path.c_str());
#else
    return 0 == mkdir(path.c_str(), mode);
#endif

  case EEXIST:
    // done!
    if (isDirExist(path)) {
      return 0;
    } else {
      return 1;
    }

  default:
    return 1;
  }
}

//! Function for printing out ring info, when there is no volume slice
int writeRings(std::vector<std::vector<int>> rings,
               std::string filename = "rings.dat");

//! Function for printing out the number of prism blocks, with or without
//! slices. Be careful when using slices!
int writePrismNum(std::string path, std::vector<int> nPrisms,
                  std::vector<int> nDefPrisms,
                  std::vector<double> heightPercent, int maxDepth,
                  int currentFrame, int firstFrame);

//! Function for printing out the coverage area and the number of rings of each
//! type
int writeRingNum(std::string path, int currentFrame, std::vector<int> nRings,
                 std::vector<double> coverageAreaXY,
                 std::vector<double> coverageAreaXZ,
                 std::vector<double> coverageAreaYZ, int maxDepth,
                 int firstFrame);

//! Function for printing out the number of rings of each
//! type in a bulk system 
int writeRingNumBulk(std::string path, int currentFrame, std::vector<int> nRings, int maxDepth,
                 int firstFrame);

//! Function for printing out the RDF, given the filename
int printRDF(std::string fileName, std::vector<double> *rdfValues,
             double binwidth, int nbin);

//! Function for printing out the number of DDCs, HCs, mixed rings, basal and
//! prismatic rings
int writeTopoBulkData(std::string path, int currentFrame, int numHC, int numDDC,
                      int mixedRings, int basalRings, int prismaticRings,
                      int firstFrame);

//! Function for writing out each prism
int writePrisms(std::vector<int> *basal1, std::vector<int> *basal2,
                int prismNum,
                molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Function for writing out cluster statistics
int writeClusterStats(std::string path, int currentFrame, int largestCluster,
                      int numOfClusters, int smallestCluster,
                      double avgClusterSize, int firstFrame);

//! Function for printing out the molecule IDs present in the slice (compatible with 
//! the LAMMPS group command 
int writeMoleculeIDsInSlice(std::string path, molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Function for printing out the molecule IDs present in the slice (compatible with 
//! the OVITO Expression Select command 
int writeMoleculeIDsExpressionSelectOVITO(std::string path, molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Write a data file for rings
int writeLAMMPSdata(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::vector<std::vector<int>> rings,
                    std::vector<std::vector<int>> bonds,
                    std::string filename = "system-rings.data");

//! Write out a LAMMPS dump file containing the RMSD per atom
int writeLAMMPSdumpINT(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<double> rmsdPerAtom, std::vector<int> atomTypes, int maxDepth,
    std::string path);

//! Write out a LAMMPS dump file containing the inSlice value for each atom
//! for a user-defined slice 
int writeLAMMPSdumpSlice(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, std::string path);

//! Write out a LAMMPS dump file containing the RMSD per atom for bulk ice
int writeLAMMPSdumpCages(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<double> rmsdPerAtom, std::vector<int> atomTypes,
    std::string path, int firstFrame);

//! Write a data file for prisms of every type
int writeLAMMPSdataAllPrisms(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, std::vector<int> atomTypes,
    int maxDepth, std::string path, bool doShapeMatching = false);

//! Write a data file for rings of every type for a monolayer
int writeLAMMPSdataAllRings(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, std::vector<int> atomTypes,
    int maxDepth, std::string path, bool isMonolayer = true);

//! Write a data file for a particular frame, writing out topological bulk ice
//! structures (DDCs/HCs)
int writeLAMMPSdataTopoBulk(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, std::vector<cage::iceType> atomTypes,
    std::string path, bool bondsBetweenDummy = false);

//! Write a data file for prisms of a single type
int writeLAMMPSdataPrisms(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> rings, bool useBondFile, std::string bondFile,
    std::vector<int> listPrism, std::vector<std::vector<int>> nList,
    std::string filename = "system-prisms.data");

//! Write out a lammps data file for DDCs or HCs, assuming that there is no
//! slice
int writeLAMMPSdataCages(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> rings, std::vector<cage::Cage> *cageList,
    cage::cageType type, int numCages,
    std::string filename = "system-cages.data");

//! Write out all cages of all types into a folder called cages inside the
//! output directory
int writeAllCages(std::string path, std::vector<cage::Cage> *cageList,
                  std::vector<std::vector<int>> rings,
                  std::vector<std::vector<int>> nList,
                  molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                  int currentFrame);

//! Write out a particular cage to a file
int writeEachCage(std::vector<int> currentCage, int cageNum,
                  cage::cageType type, std::vector<std::vector<int>> rings,
                  molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Write out the basal rings of a particular Hexagonal cage
int writeBasalRingsHex(std::vector<int> currentCage, int cageNum,
                       std::vector<std::vector<int>> nList,
                       std::vector<std::vector<int>> rings);

//! Write out the basal rings for a particular prism
int writeBasalRingsPrism(
    std::vector<int> *basal1, std::vector<int> *basal2, int prismNum,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, bool isDeformed);

//! Generic function for writing out to a dump file
int writeDump(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
              std::string path, std::string outFile);

//! Function for printing out Q6, Cij and averaged Q3 values as single columns
//! to text files The file names are cij, q6, q3
int writeHisto(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
               std::vector<std::vector<int>> nList, std::vector<double> avgQ6);

//! Function for printing the largest ice cluster
int writeCluster(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 std::string fileName = "cluster.txt", bool isSlice = false,
                 int largestIceCluster = 0);

//! Function for writing out the XYZ files for each cluster
int writeXYZcluster(std::string path,
                    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::vector<int> atoms, int clusterID, cage::cageType type);
} // namespace sout
#endif // __SEAMS_OUTPUT_H_
