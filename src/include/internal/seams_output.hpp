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
#include <generic.hpp>
#include <iostream>
#include <memory>
#include <mol_sys.hpp>

#include <filesystem>
namespace fs = std::filesystem;

namespace sout {

/**
 *  Inline function for creating the desired directory (and parents).
 *  Uses std::filesystem for portability.
 *  @param[in] path The path of the directory
 *  @return 0 on success, 1 on failure
 */
[[nodiscard]] inline int makePath(const std::string &path) {
  std::error_code ec;
  if (fs::is_directory(path)) {
    return 0;
  }
  if (fs::create_directories(path, ec)) {
    return 0;
  }
  return ec ? 1 : 0;
}

//! Function for printing out ring info, when there is no volume slice
[[nodiscard]] int writeRings(const std::vector<std::vector<int>> &rings,
               std::string filename = "rings.dat");

//! Function for printing out the number of prism blocks, with or without
//! slices. Be careful when using slices!
[[nodiscard]] int writePrismNum(std::string path, const std::vector<int> &nPrisms,
                  const std::vector<int> &nDefPrisms,
                  const std::vector<double> &heightPercent, int maxDepth,
                  int currentFrame, int firstFrame);

//! Function for printing out the coverage area and the number of rings of each
//! type
[[nodiscard]] int writeRingNum(std::string path, int currentFrame, const std::vector<int> &nRings,
                 const std::vector<double> &coverageAreaXY,
                 const std::vector<double> &coverageAreaXZ,
                 const std::vector<double> &coverageAreaYZ, int maxDepth,
                 int firstFrame);

//! Function for printing out the number of rings of each
//! type in a bulk system
[[nodiscard]] int writeRingNumBulk(std::string path, int currentFrame, const std::vector<int> &nRings, int maxDepth,
                 int firstFrame);

//! Function for printing out the RDF, given the filename
[[nodiscard]] int printRDF(std::string fileName, std::vector<double> *rdfValues,
             double binwidth, int nbin);

//! Function for printing out the number of DDCs, HCs, mixed rings, basal and
//! prismatic rings
[[nodiscard]] int writeTopoBulkData(std::string path, int currentFrame, int numHC, int numDDC,
                      int mixedRings, int basalRings, int prismaticRings,
                      int firstFrame);

//! Function for writing out each prism
[[nodiscard]] int writePrisms(std::vector<int> *basal1, std::vector<int> *basal2,
                int prismNum,
                molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Function for writing out cluster statistics
[[nodiscard]] int writeClusterStats(std::string path, int currentFrame, int largestCluster,
                      int numOfClusters, int smallestCluster,
                      double avgClusterSize, int firstFrame);

//! Function for printing out the molecule IDs present in the slice (compatible with
//! the LAMMPS group command
[[nodiscard]] int writeMoleculeIDsInSlice(std::string path, molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Function for printing out the molecule IDs present in the slice (compatible with
//! the OVITO Expression Select command
[[nodiscard]] int writeMoleculeIDsExpressionSelectOVITO(std::string path, molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Write a data file for rings
[[nodiscard]] int writeLAMMPSdata(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::vector<std::vector<int>> rings,
                    std::vector<std::vector<int>> bonds,
                    std::string filename = "system-rings.data");

//! Write out a LAMMPS dump file containing the RMSD per atom
[[nodiscard]] int writeLAMMPSdumpINT(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    const std::vector<double> &rmsdPerAtom, const std::vector<int> &atomTypes, int maxDepth,
    std::string path);

//! Write out a LAMMPS dump file containing the inSlice value for each atom
//! for a user-defined slice
[[nodiscard]] int writeLAMMPSdumpSlice(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, std::string path);

//! Write out a LAMMPS dump file containing the RMSD per atom for bulk ice
[[nodiscard]] int writeLAMMPSdumpCages(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    const std::vector<double> &rmsdPerAtom, const std::vector<int> &atomTypes,
    std::string path, int firstFrame);

//! Write a data file for prisms of every type
[[nodiscard]] int writeLAMMPSdataAllPrisms(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    const std::vector<std::vector<int>> &nList, const std::vector<int> &atomTypes,
    int maxDepth, std::string path, bool doShapeMatching = false);

//! Write a data file for rings of every type for a monolayer
[[nodiscard]] int writeLAMMPSdataAllRings(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    const std::vector<std::vector<int>> &nList, const std::vector<int> &atomTypes,
    int maxDepth, std::string path, bool isMonolayer = true);

//! Write a data file for a particular frame, writing out topological bulk ice
//! structures (DDCs/HCs)
[[nodiscard]] int writeLAMMPSdataTopoBulk(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    const std::vector<std::vector<int>> &nList, const std::vector<cage::iceType> &atomTypes,
    std::string path, bool bondsBetweenDummy = false);

//! Write a data file for prisms of a single type
[[nodiscard]] int writeLAMMPSdataPrisms(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    const std::vector<std::vector<int>> &rings, bool useBondFile, std::string bondFile,
    const std::vector<int> &listPrism, const std::vector<std::vector<int>> &nList,
    std::string filename = "system-prisms.data");

//! Write out a lammps data file for DDCs or HCs, assuming that there is no
//! slice
[[nodiscard]] int writeLAMMPSdataCages(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    const std::vector<std::vector<int>> &rings, std::vector<cage::Cage> *cageList,
    cage::cageType type, int numCages,
    std::string filename = "system-cages.data");

//! Write out all cages of all types into a folder called cages inside the
//! output directory
[[nodiscard]] int writeAllCages(std::string path, std::vector<cage::Cage> *cageList,
                  const std::vector<std::vector<int>> &rings,
                  const std::vector<std::vector<int>> &nList,
                  molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                  int currentFrame);

//! Write out a particular cage to a file
[[nodiscard]] int writeEachCage(const std::vector<int> &currentCage, int cageNum,
                  cage::cageType type, const std::vector<std::vector<int>> &rings,
                  molSys::PointCloud<molSys::Point<double>, double> *yCloud);

//! Write out the basal rings of a particular Hexagonal cage
[[nodiscard]] int writeBasalRingsHex(const std::vector<int> &currentCage, int cageNum,
                       const std::vector<std::vector<int>> &nList,
                       const std::vector<std::vector<int>> &rings);

//! Write out the basal rings for a particular prism
[[nodiscard]] int writeBasalRingsPrism(
    std::vector<int> *basal1, std::vector<int> *basal2, int prismNum,
    const std::vector<std::vector<int>> &nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, bool isDeformed);

//! Generic function for writing out to a dump file
[[nodiscard]] int writeDump(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
              std::string path, std::string outFile);

//! Function for printing out Q6, Cij and averaged Q3 values as single columns
//! to text files The file names are cij, q6, q3
[[nodiscard]] int writeHisto(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
               const std::vector<std::vector<int>> &nList, const std::vector<double> &avgQ6);

//! Function for printing the largest ice cluster
[[nodiscard]] int writeCluster(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 std::string fileName = "cluster.txt", bool isSlice = false,
                 int largestIceCluster = 0);

//! Function for writing out the XYZ files for each cluster
[[nodiscard]] int writeXYZcluster(std::string path,
                    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    const std::vector<int> &atoms, int clusterID, cage::cageType type);
} // namespace sout
#endif // __SEAMS_OUTPUT_H_
