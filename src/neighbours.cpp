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

#include <cmath>
#include <iostream>
#include <numeric>
#include <span>
#include <utility>

#include <neighbours.hpp>
#include <simd_distance.hpp>

#ifdef SEAMS_HAS_VESIN
#include <vesin.h>
#endif

namespace {

/**
 * @details Reports whether a point cloud carries the three box lengths that
 *  the minimum image convention needs. A cloud whose trajectory frame failed
 *  to load leaves molSys::PointCloud::box empty, and every distance routine
 *  indexes it unconditionally.
 * @param[in] yCloud The input molSys::PointCloud.
 * @return True when all three box lengths are present.
 */
bool hasPeriodicBox(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  if (yCloud.box.size() >= 3) {
    return true;
  }
  std::cerr << "The point cloud has no simulation box; cannot build a "
               "neighbour list.\n";
  return false;
}

/**
 * @details Builds the reverse of molSys::PointCloud::idIndexMap: a dense table
 *  mapping a particle index to its atom ID. The forward map answers
 *  ID -> index in constant time; the reverse direction is needed once per
 *  particle when emitting neighbour lists keyed by atom ID, and a single
 *  linear pass over the map amortises it.
 * @param[in] yCloud The input molSys::PointCloud.
 * @return Table of length molSys::PointCloud::nop, holding the atom ID at each
 *  index, or -1 where the map has no entry for that index.
 */
std::vector<int> indexToIDTable(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  std::vector<int> indexToID(yCloud.nop, -1);
  for (const auto &[atomID, index] : yCloud.idIndexMap) {
    if (index >= 0 && index < yCloud.nop) {
      indexToID[index] = atomID;
    }
  }
  return indexToID;
}

#ifdef SEAMS_HAS_VESIN
/**
 * @details Runs the vesin cell list over the particles named by @a subset and
 *  reports the neighbour pairs as cloud indices. Cell-list construction is
 *  linear in the particle count, against the quadratic cost of comparing every
 *  pair.
 * @param[in] yCloud The input molSys::PointCloud.
 * @param[in] subset Cloud indices of the particles to search over.
 * @param[in] rcutoff Distance cutoff, within which two atoms are neighbours.
 * @param[out] pairs Neighbour pairs, as cloud indices, in both directions.
 * @return True when vesin produced a neighbour list, false when the caller
 *  should fall back to the brute-force path.
 */
bool cellListPairs(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                   const std::vector<int> &subset, double rcutoff,
                   std::vector<std::pair<int, int>> &pairs) {
  const size_t nSubset = subset.size();
  std::vector<std::array<double, 3>> positions(nSubset);
  for (size_t i = 0; i < nSubset; i++) {
    const int idx = subset[i];
    positions[i] = {yCloud.pts[idx].x, yCloud.pts[idx].y, yCloud.pts[idx].z};
  }

  // Box matrix (row-major, diagonal for orthorhombic)
  double box[3][3] = {{yCloud.box[0], 0.0, 0.0},
                      {0.0, yCloud.box[1], 0.0},
                      {0.0, 0.0, yCloud.box[2]}};
  bool periodic[3] = {true, true, true};

  VesinOptions options;
  options.cutoff = rcutoff;
  options.full = true; // full neighbor list (both i->j and j->i)
  options.sorted = false;
  options.algorithm = VesinAutoAlgorithm;
  options.return_shifts = false;
  options.return_distances = false;
  options.return_vectors = false;

  VesinNeighborList neighbors;
  const char *error_message = nullptr;
  VesinDevice device = {VesinCPU, 0};

  const int status = vesin_neighbors(
      reinterpret_cast<const double (*)[3]>(positions.data()), nSubset, box,
      periodic, device, options, &neighbors, &error_message);

  if (status != 0) {
    std::cerr << "Vesin failed: " << (error_message ? error_message : "unknown")
              << "; falling back to brute force.\n";
    return false;
  }

  pairs.clear();
  pairs.reserve(neighbors.length);
  for (size_t k = 0; k < neighbors.length; k++) {
    pairs.emplace_back(subset[neighbors.pairs[k][0]],
                       subset[neighbors.pairs[k][1]]);
  }

  vesin_free(&neighbors);
  return true;
}
#endif

/**
 * @details Seeds a neighbour list so that row @a i begins with the atom ID of
 *  particle @a i, which is the row-header convention used throughout d-SEAMS.
 * @param[in] indexToID Reverse index-to-ID table.
 * @param[in] nop Number of particles.
 * @return Neighbour list of @a nop rows, each holding its own atom ID.
 */
std::vector<std::vector<int>> seedWithSelfIDs(const std::vector<int> &indexToID,
                                              int nop) {
  std::vector<std::vector<int>> nList(nop);
  for (int iatom = 0; iatom < nop; iatom++) {
    if (indexToID[iatom] == -1) {
      std::cerr << "Something is wrong with your idIndexMap!\n";
      continue;
    }
    nList[iatom].push_back(indexToID[iatom]);
  }
  return nList;
}

} // namespace

/**
 * @details Function for building neighbour lists for each
 *  particle. Inefficient brute-force \f$ O(n^2) \f$ implementation.
 *  This generates the full neighbour list, by ID.
 * @param[in] rcutoff Distance cutoff, within which two atoms are neighbours.
 * @param[in] yCloud The input molSys::PointCloud
 * @param[in] typeI Type ID of particles of type I.
 * @param[in] typeJ Type ID of particles of type J.
 * @return Row-ordered full neighbour list, by atom ID.
 */
std::vector<std::vector<int>>
nneigh::neighList(double rcutoff,
                  const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                  int typeI, int typeJ) {
  if (!hasPeriodicBox(yCloud)) {
    return {};
  }

  const std::vector<int> indexToID = indexToIDTable(yCloud);
  std::vector<std::vector<int>> nList = seedWithSelfIDs(indexToID, yCloud.nop);

  // Compare squared distances so that the per-pair square root is avoided
  const double rcutoffSq = rcutoff * rcutoff;

  // pairs of atoms of type I and J
  // Loop through every iatom and find nearest neighbours within rcutoff
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    if (yCloud.pts[iatom].type != typeI) {
      continue;
    }
    const int iatomIndex = indexToID[iatom];
    // Loop through the other atoms
    for (int jatom = 0; jatom < yCloud.nop; jatom++) {
      if (yCloud.pts[jatom].type != typeJ) {
        continue;
      }
      // If the distance is greater than rcutoff, continue
      if (gen::periodicDistSq(yCloud, iatom, jatom) > rcutoffSq) {
        continue;
      }

      // Update the neighbour indices with atom IDs for iatom and jatom both
      // (full list)
      nList[iatom].push_back(indexToID[jatom]);
      nList[jatom].push_back(iatomIndex);

    } // End of loop through jatom
  }   // End of loop for iatom

  return nList;
}

/**
 * @details Function for building neighbour lists for each
 *  particle of only one type. Inefficient brute-force \f$ O(n^2) \f$
 * implementation. This generates the full neighbour list, by ID. This function
 * will only work for building a neighbour list between one type of particles.
 * @param[in] rcutoff Distance cutoff, within which two atoms are neighbours.
 * @param[in] yCloud The input molSys::PointCloud
 * @param[in] typeI Type ID of the \f$ i^{th} \f$ particle type.
 * @return Row-ordered full neighbour list, by atom ID.
 */
std::vector<std::vector<int>>
nneigh::neighListO(double rcutoff,
                   const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                   int typeI) {
  if (!hasPeriodicBox(yCloud)) {
    return {};
  }

  std::vector<std::vector<int>>
      nList; // Vector of vectors of the neighbour list

  // Build index-to-atomID lookup once; every path below needs it
  const std::vector<int> indexToID = indexToIDTable(yCloud);

  // Collect indices of typeI atoms
  std::vector<int> typeIIndices;
  for (int i = 0; i < yCloud.nop; i++) {
    if (yCloud.pts[i].type == typeI) {
      typeIIndices.push_back(i);
    }
  }

#ifdef SEAMS_HAS_VESIN
  // O(n) cell-list neighbor search via vesin
  {
    std::vector<std::pair<int, int>> pairs;
    if (cellListPairs(yCloud, typeIIndices, rcutoff, pairs)) {
      // Initialize nList for ALL atoms (not just typeI)
      nList = seedWithSelfIDs(indexToID, yCloud.nop);

      // Fill neighbor pairs from vesin output, which is already bidirectional
      for (const auto &[iatom, jatom] : pairs) {
        nList[iatom].push_back(indexToID[jatom]);
      }

      return nList;
    }
    // If vesin failed, fall through to brute-force
  }
#endif

  // Initialize and fill the first element with the current atom ID whose
  // neighbour list will be filled
  nList = seedWithSelfIDs(indexToID, yCloud.nop);

  const double rcutoffSq = rcutoff * rcutoff;
  const double bx = yCloud.box[0];
  const double by = yCloud.box[1];
  const double bz = yCloud.box[2];

  // Scratch buffers for the batched distance kernel, sized once for the
  // largest batch and reused across iatom
  const size_t nTypeIAtoms = typeIIndices.size();
  std::vector<double> dx(nTypeIAtoms), dy(nTypeIAtoms), dz(nTypeIAtoms);
  std::vector<double> distSq(nTypeIAtoms);

  // Loop through every iatom and find nearest neighbours within rcutoff
  for (size_t ii = 0; ii < nTypeIAtoms; ii++) {
    const int iatom = typeIIndices[ii];
    const int iatomIndex = indexToID[iatom];

    // Collect coordinate differences for all j > i of typeI
    const size_t remaining = nTypeIAtoms - ii - 1;
    if (remaining == 0) {
      continue;
    }

    const double xi = yCloud.pts[iatom].x;
    const double yi = yCloud.pts[iatom].y;
    const double zi = yCloud.pts[iatom].z;
    for (size_t jj = 0; jj < remaining; jj++) {
      const int jatom = typeIIndices[ii + 1 + jj];
      dx[jj] = xi - yCloud.pts[jatom].x;
      dy[jj] = yi - yCloud.pts[jatom].y;
      dz[jj] = zi - yCloud.pts[jatom].z;
    }

    // Batch compute squared periodic distances (SIMD when available)
    seams::BatchPeriodicDistSq(std::span(dx).first(remaining),
                               std::span(dy).first(remaining),
                               std::span(dz).first(remaining), bx, by, bz,
                               std::span(distSq).first(remaining));

    // Filter by cutoff and update neighbour lists
    for (size_t jj = 0; jj < remaining; jj++) {
      if (distSq[jj] > rcutoffSq) {
        continue;
      }
      const int jatom = typeIIndices[ii + 1 + jj];
      nList[iatom].push_back(indexToID[jatom]);
      nList[jatom].push_back(iatomIndex);
    }
  }   // End of loop for iatom

  return nList;
}

/**
 * @details Function for building neighbour lists for each
 *  particle of only one type. Inefficient brute-force \f$ O(n^2) \f$
 *  implementation. This generates the half neighbour list, by ID. This function
 *  will only work for building a neighbour list between one type of particles.
 * @param[in] rcutoff Distance cutoff, within which two atoms are neighbours.
 * @param[in] yCloud The input molSys::PointCloud
 * @param[in] typeI Type ID of the \f$ i^{th} \f$ particle type.
 * @return Row-ordered half neighbour list, by atom ID.
 */
std::vector<std::vector<int>>
nneigh::halfNeighList(double rcutoff,
                      const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                      int typeI) {
  if (!hasPeriodicBox(yCloud)) {
    return {};
  }

  const std::vector<int> indexToID = indexToIDTable(yCloud);
  std::vector<std::vector<int>> nList = seedWithSelfIDs(indexToID, yCloud.nop);

  // Compare squared distances so that the per-pair square root is avoided
  const double rcutoffSq = rcutoff * rcutoff;

  // Loop through every iatom and find nearest neighbours within rcutoff
  for (int iatom = 0; iatom < yCloud.nop - 1; iatom++) {
    if (yCloud.pts[iatom].type != typeI) {
      continue;
    }
    // Loop through the other atoms
    for (int jatom = iatom + 1; jatom < yCloud.nop; jatom++) {
      if (yCloud.pts[jatom].type != typeI) {
        continue;
      }
      // If the distance is greater than rcutoff, continue
      if (gen::periodicDistSq(yCloud, iatom, jatom) > rcutoffSq) {
        continue;
      }

      // Update the neighbour indices with the atom ID of jatom (half list)
      nList[iatom].push_back(indexToID[jatom]);

    } // End of loop through jatom
  }   // End of loop for iatom

  return nList;
}

/**
 * @details Function for creating a neighbour list by index (from scratch)
 * instead of by atom ID. The ordering is with respect to the pointCloud with
 * the coordinates.The first element is the atom for which the other atom
 * indices are neighbours For example, if the neighbours of 1 are 2, 3, 4 the
 * sub-vector would have 1 2 3 4
 * @param[in] yCloud The input molSys::PointCloud
 * @param[in] cutoff Distance cutoff, within which two atoms are neighbours.
 * @return Row-ordered full neighbour list, by index, NOT atom ID.
 */
std::vector<std::vector<int>> nneigh::getNewNeighbourListByIndex(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff) {
  //
  if (!hasPeriodicBox(yCloud)) {
    return {};
  }

  std::vector<std::vector<int>> nList(yCloud.nop);

  // Initialize and fill the first element with the index of the atom whose
  // neighbour list will be filled
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    nList[iatom].push_back(iatom);
  } // end of init
  // -------------------------------------------------------
#ifdef SEAMS_HAS_VESIN
  // O(n) cell-list neighbour search via vesin, over every particle
  {
    std::vector<int> allIndices(yCloud.nop);
    std::iota(allIndices.begin(), allIndices.end(), 0);

    std::vector<std::pair<int, int>> pairs;
    if (cellListPairs(yCloud, allIndices, cutoff, pairs)) {
      // vesin returns both directions, so each row is filled from its own pairs
      for (const auto &[iatom, jatom] : pairs) {
        nList[iatom].push_back(jatom);
      }
      return nList;
    }
    // If vesin failed, fall through to brute-force
  }
#endif
  // -------------------------------------------------------
  // Compare squared distances so that the per-pair square root is avoided
  const double cutoffSq = cutoff * cutoff;
  const double bx = yCloud.box[0];
  const double by = yCloud.box[1];
  const double bz = yCloud.box[2];

  // Scratch buffers for the batched distance kernel, sized once for the
  // largest batch and reused across iatom
  std::vector<double> dx(yCloud.nop), dy(yCloud.nop), dz(yCloud.nop);
  std::vector<double> distSq(yCloud.nop);

  // Loop through every iatom and find nearest neighbours within cutoff
  for (int iatom = 0; iatom < yCloud.nop - 1; iatom++) {
    const size_t remaining = static_cast<size_t>(yCloud.nop - iatom - 1);

    const double xi = yCloud.pts[iatom].x;
    const double yi = yCloud.pts[iatom].y;
    const double zi = yCloud.pts[iatom].z;
    for (size_t jj = 0; jj < remaining; jj++) {
      const int jatom = iatom + 1 + static_cast<int>(jj);
      dx[jj] = xi - yCloud.pts[jatom].x;
      dy[jj] = yi - yCloud.pts[jatom].y;
      dz[jj] = zi - yCloud.pts[jatom].z;
    }

    // Batch compute squared periodic distances (SIMD when available)
    seams::BatchPeriodicDistSq(std::span(dx).first(remaining),
                               std::span(dy).first(remaining),
                               std::span(dz).first(remaining), bx, by, bz,
                               std::span(distSq).first(remaining));

    // Update the neighbour indices for iatom and jatom both (full list)
    for (size_t jj = 0; jj < remaining; jj++) {
      if (distSq[jj] > cutoffSq) {
        continue;
      }
      const int jatom = iatom + 1 + static_cast<int>(jj);
      nList[iatom].push_back(jatom);
      nList[jatom].push_back(iatom);
    } // End of loop through jatom
  }   // End of loop for iatom

  return nList;
} // end of function

/**
 * @details Function for getting the neighbour list by index instead of by atom
 *  ID from a previously constructed input neighbour list by ID. The ordering is
 *  with respect to the pointCloud with the coordinates.The first element is the
 *  atom for which the other atom indices are neighbours For example, if the
 *  neighbours of 1 are 2, 3, 4 the sub-vector would have 1 2 3 4
 * @param[in] yCloud The input molSys::PointCloud
 * @param[in] nList Full neighbour list, by atom ID.
 * @return Row-ordered full neighbour list, by index, NOT atom ID.
 */
std::vector<std::vector<int>> nneigh::neighbourListByIndex(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList) {
  //
  std::vector<std::vector<int>> indexNlist; // Desired neighbour list of indices
  int iatomID, jatomID;                     // Atom IDs
  int iatomIndex, jatomIndex;               // Indices of iatom and jatom
  int nnumNeighbours;                       // Number of nearest neighbours

  // Loop through every atom whose neighbours are contained in the neighbour
  // list
  for (int iatom = 0; iatom < nList.size(); iatom++) {
    iatomID = nList[iatom][0]; // Atom ID
    // Get the index of iatom
    auto gotI = yCloud.idIndexMap.find(iatomID);
    if (gotI != yCloud.idIndexMap.end()) {
      iatomIndex = gotI->second;
    } // found iatomIndex
    //
    nnumNeighbours = nList[iatomIndex].size() - 1;
    // Update the new neighbour list
    indexNlist.push_back(
        std::vector<int>()); // Empty vector for the index iatom
    // Fill the first element with the atom ID of iatom itself
    indexNlist[iatom].push_back(iatomIndex);
    //
    // Loop through the neighbours of iatom
    for (int jatom = 1; jatom <= nnumNeighbours; jatom++) {
      jatomID = nList[iatomIndex][jatom]; // Atom ID of neighbour
      //
      // Get the index of the j^th atom
      auto gotJ = yCloud.idIndexMap.find(jatomID);
      if (gotJ != yCloud.idIndexMap.end()) {
        jatomIndex = gotJ->second;
      } // found jatomIndex
      // Add to the neighbour list
      indexNlist[iatom].push_back(jatomIndex);
    } // end of loop through neighbours
  }   // end of loop through every atom

  // Return the new neighbour list
  return indexNlist;
}

/**
 * @details Deletes the memory of a
 *  vector of vectors. Call this before creating the neighbour list for a new
 *  frame.
 *  @param[in, out] nList Vector of vectors, of the neighbour list to be erased.
 */
int nneigh::clearNeighbourList(std::vector<std::vector<int>> &nList) {
  nList.clear();
  nList.shrink_to_fit();
  return 0;
}
