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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <set>
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
    vesin_free(&neighbors);
    return false;
  }

  pairs.clear();
  pairs.reserve(neighbors.length);
  for (size_t k = 0; k < neighbors.length; k++) {
    const int iatom = subset[neighbors.pairs[k][0]];
    const int jatom = subset[neighbors.pairs[k][1]];
    // A cell list enumerates periodic images, so a particle can appear as its
    // own neighbour through an image, and one neighbour can arrive through
    // several images at once. Both happen as soon as the box stops being
    // larger than twice the cutoff. The minimum image convention that the
    // brute-force path applies admits each ordered pair once and never the
    // self pair, so reduce to that here rather than letting box size change
    // the meaning of a neighbour list.
    if (iatom == jatom) {
      continue;
    }
    pairs.emplace_back(iatom, jatom);
  }

  vesin_free(&neighbors);

  std::sort(pairs.begin(), pairs.end());
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());

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

//! True when index has a real atom ID in the reverse table.
bool hasAtomID(const std::vector<int> &indexToID, int index) {
  return index >= 0 && static_cast<size_t>(index) < indexToID.size() &&
         indexToID[index] != -1;
}

//! Appends jatom's ID to iatom's row only when both particles have IDs.
//! An unmapped index is left as an empty row, which callers treat as
//! "no self header, skip this particle".
void appendNeighbourID(std::vector<std::vector<int>> &nList,
                       const std::vector<int> &indexToID, int iatom,
                       int jatom) {
  if (!hasAtomID(indexToID, iatom) || !hasAtomID(indexToID, jatom)) {
    return;
  }
  nList[iatom].push_back(indexToID[jatom]);
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

  // Pairs of type I and type J. When the types coincide the full i x j
  // product would write each unordered pair twice and accept iatom == jatom
  // (distance 0). Walk j > i in that case, matching halfNeighList.
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    if (yCloud.pts[iatom].type != typeI) {
      continue;
    }
    const int jStart = (typeI == typeJ) ? iatom + 1 : 0;
    for (int jatom = jStart; jatom < yCloud.nop; jatom++) {
      if (yCloud.pts[jatom].type != typeJ) {
        continue;
      }
      if (gen::periodicDistSq(yCloud, iatom, jatom) > rcutoffSq) {
        continue;
      }
      appendNeighbourID(nList, indexToID, iatom, jatom);
      appendNeighbourID(nList, indexToID, jatom, iatom);
    }
  }

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
        appendNeighbourID(nList, indexToID, iatom, jatom);
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
      appendNeighbourID(nList, indexToID, iatom, jatom);
      appendNeighbourID(nList, indexToID, jatom, iatom);
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
      appendNeighbourID(nList, indexToID, iatom, jatom);

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
  std::vector<std::vector<int>> indexNlist(nList.size());
  int iatomID, jatomID;
  int iatomIndex;

  for (size_t iatom = 0; iatom < nList.size(); iatom++) {
    if (nList[iatom].empty()) {
      continue;
    }
    iatomID = nList[iatom][0];
    auto gotI = yCloud.idIndexMap.find(iatomID);
    if (gotI == yCloud.idIndexMap.end()) {
      continue;
    }
    iatomIndex = gotI->second;
    indexNlist[iatom].push_back(iatomIndex);
    for (size_t j = 1; j < nList[iatom].size(); j++) {
      jatomID = nList[iatom][j];
      auto gotJ = yCloud.idIndexMap.find(jatomID);
      if (gotJ == yCloud.idIndexMap.end()) {
        continue;
      }
      indexNlist[iatom].push_back(gotJ->second);
    }
  }

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

/**
 * @details Bonded graph from the k nearest neighbours of each particle,
 *  union-symmetrized: i and j are bonded when either lists the other among
 *  its k nearest. Candidates come from the same cell-list search as
 *  neighListO at @a candidateCutoff, which must comfortably exceed the k-th
 *  neighbour distance. On an undistorted tetrahedral lattice with k = 4 this
 *  equals the first-shell cutoff graph; under thermal distortion it keeps
 *  the neighbour identities a hard cutoff loses, which is what the ring and
 *  cage predicates consume.
 * @param[in] yCloud The input molSys::PointCloud.
 * @param[in] k Number of nearest neighbours each particle nominates.
 * @param[in] candidateCutoff Cell-list search radius for the candidates.
 * @param[in] typeI Type ID of the particles in the graph.
 * @return Row-ordered full neighbour list by atom ID, one row per atom with
 *  the leading self entry, like neighListO.
 */
std::vector<std::vector<int>>
nneigh::kNearestNeighbourList(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, int typeI, bool mutual) {
  std::vector<std::vector<int>> candidateRows =
      nneigh::neighListO(candidateCutoff, yCloud, typeI);
  if (candidateRows.empty() || k <= 0) {
    return candidateRows;
  }

  // Directed nominations: each atom's k nearest among the candidates
  std::vector<std::vector<int>> nominated(yCloud.nop);
  std::vector<std::pair<double, int>> byDist; // (distance^2, neighbour index)
  for (int i = 0; i < yCloud.nop; i++) {
    const auto &row = candidateRows[i];
    if (row.size() <= 1) {
      continue;
    }
    byDist.clear();
    for (size_t m = 1; m < row.size(); m++) {
      const auto it = yCloud.idIndexMap.find(row[m]);
      if (it == yCloud.idIndexMap.end()) {
        continue;
      }
      byDist.emplace_back(gen::periodicDistSq(yCloud, i, it->second),
                          it->second);
    }
    const size_t keep = std::min(static_cast<size_t>(k), byDist.size());
    std::partial_sort(byDist.begin(), byDist.begin() + keep, byDist.end());
    nominated[i].reserve(keep);
    for (size_t m = 0; m < keep; m++) {
      nominated[i].push_back(byDist[m].second);
    }
  }

  // Exactness fallback: an atom whose k-th neighbour lies beyond the
  // candidate cutoff got a truncated nomination; recompute those few by a
  // brute-force scan so the k nearest are exact regardless of the cutoff
  for (int i = 0; i < yCloud.nop; i++) {
    if (yCloud.pts[i].type != typeI ||
        static_cast<int>(nominated[i].size()) >= k) {
      continue;
    }
    byDist.clear();
    for (int j = 0; j < yCloud.nop; j++) {
      if (j == i || yCloud.pts[j].type != typeI) {
        continue;
      }
      byDist.emplace_back(gen::periodicDistSq(yCloud, i, j), j);
    }
    const size_t keep = std::min(static_cast<size_t>(k), byDist.size());
    std::partial_sort(byDist.begin(), byDist.begin() + keep, byDist.end());
    nominated[i].clear();
    for (size_t m = 0; m < keep; m++) {
      nominated[i].push_back(byDist[m].second);
    }
  }

  // Union symmetrization into ID rows with the leading self entry
  std::vector<std::vector<int>> out(yCloud.nop);
  for (int i = 0; i < yCloud.nop; i++) {
    out[i].push_back(candidateRows[i].empty() ? yCloud.pts[i].atomID
                                              : candidateRows[i][0]);
  }
  std::set<std::pair<int, int>> bonds;
  if (mutual) {
    // Intersection symmetrization: a bond requires both nominations. In a
    // crystal the first shell is mutual, so this equals the union graph
    // there; on disordered packings one-sided nominations vanish, which is
    // what starves the accidental ring complexes
    std::set<std::pair<int, int>> directed;
    for (int i = 0; i < yCloud.nop; i++) {
      for (const int j : nominated[i]) {
        directed.emplace(i, j);
      }
    }
    for (const auto &[i, j] : directed) {
      if (directed.count({j, i})) {
        bonds.emplace(std::min(i, j), std::max(i, j));
      }
    }
  } else {
    for (int i = 0; i < yCloud.nop; i++) {
      for (const int j : nominated[i]) {
        bonds.emplace(std::min(i, j), std::max(i, j));
      }
    }
  }
  for (const auto &[i, j] : bonds) {
    out[i].push_back(yCloud.pts[j].atomID);
    out[j].push_back(yCloud.pts[i].atomID);
  }
  return out;
}

/**
 * @details The shell-separation certificate for the exact reduction of the
 *  k-nearest bonded graph to the cutoff graph. Returns the largest k-th
 *  neighbour distance and the smallest (k+1)-th neighbour distance over all
 *  particles of the type. Whenever
 *  max_i d_k(i) <= rcutoff <= min_i d_{k+1}(i), every particle's cutoff
 *  neighbourhood is exactly its k nearest, so the union-symmetrized
 *  k-nearest graph and the cutoff graph coincide edge for edge and every
 *  graph predicate downstream -- rings, cages, affiliation -- is identical.
 * @return {max over i of d_k(i), min over i of d_{k+1}(i)}; zeros when the
 *  system has too few particles.
 */
std::pair<double, double> nneigh::shellSeparation(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    int typeI) {
  double maxKth = 0.0;
  double minNext = 0.0;
  bool haveNext = false;
  std::vector<double> dists;
  for (int i = 0; i < yCloud.nop; i++) {
    if (yCloud.pts[i].type != typeI) {
      continue;
    }
    dists.clear();
    for (int j = 0; j < yCloud.nop; j++) {
      if (j == i || yCloud.pts[j].type != typeI) {
        continue;
      }
      dists.push_back(gen::periodicDistSq(yCloud, i, j));
    }
    if (static_cast<int>(dists.size()) < k + 1) {
      continue;
    }
    std::partial_sort(dists.begin(), dists.begin() + k + 1, dists.end());
    maxKth = std::max(maxKth, std::sqrt(dists[k - 1]));
    const double next = std::sqrt(dists[k]);
    minNext = haveNext ? std::min(minNext, next) : next;
    haveNext = true;
  }
  return {maxKth, haveNext ? minNext : 0.0};
}

nneigh::SkinNeighborList::SkinNeighborList(double cutoff, double skin,
                                           int typeI)
    : cutoff_(cutoff), skin_(skin), cutoffSq_(cutoff * cutoff),
      triggerSq_((0.5 * skin) * (0.5 * skin)), typeI_(typeI) {}

bool nneigh::SkinNeighborList::mustRebuild(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) const {
  if (static_cast<int>(x0_.size()) != yCloud.nop || yCloud.box.size() < 3) {
    return true;
  }
  for (int k = 0; k < 3; k++) {
    if (std::fabs(yCloud.box[k] - box0_[static_cast<std::size_t>(k)]) >
        1e-8) {
      return true;
    }
  }
  std::vector<double> old(3);
  for (int i = 0; i < yCloud.nop; i++) {
    if (yCloud.pts[i].type != typeI_) {
      continue;
    }
    old[0] = x0_[static_cast<std::size_t>(i)];
    old[1] = y0_[static_cast<std::size_t>(i)];
    old[2] = z0_[static_cast<std::size_t>(i)];
    const double r = gen::unWrappedDistFromPoint(yCloud, i, old);
    if (r * r > triggerSq_) {
      return true;
    }
  }
  return false;
}

void nneigh::SkinNeighborList::rebuildCandidates(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  candidates_.clear();
  const double wide = cutoff_ + skin_;
#ifdef SEAMS_HAS_VESIN
  std::vector<int> subset;
  subset.reserve(static_cast<std::size_t>(yCloud.nop));
  for (int i = 0; i < yCloud.nop; i++) {
    if (yCloud.pts[i].type == typeI_) {
      subset.push_back(i);
    }
  }
  std::vector<std::pair<int, int>> pairs;
  if (cellListPairs(yCloud, subset, wide, pairs)) {
    for (const auto &[i, j] : pairs) {
      if (i < j) {
        candidates_.emplace_back(i, j);
      }
    }
  } else
#endif
  {
    const auto wideList = nneigh::neighListO(wide, yCloud, typeI_);
    for (int i = 0; i < yCloud.nop; i++) {
      if (static_cast<int>(wideList.size()) <= i || wideList[i].size() <= 1) {
        continue;
      }
      for (std::size_t n = 1; n < wideList[i].size(); n++) {
        const auto it = yCloud.idIndexMap.find(wideList[i][n]);
        if (it == yCloud.idIndexMap.end()) {
          continue;
        }
        const int j = it->second;
        if (i < j) {
          candidates_.emplace_back(i, j);
        }
      }
    }
  }
  x0_.resize(static_cast<std::size_t>(yCloud.nop));
  y0_.resize(static_cast<std::size_t>(yCloud.nop));
  z0_.resize(static_cast<std::size_t>(yCloud.nop));
  for (int i = 0; i < yCloud.nop; i++) {
    x0_[static_cast<std::size_t>(i)] = yCloud.pts[i].x;
    y0_[static_cast<std::size_t>(i)] = yCloud.pts[i].y;
    z0_[static_cast<std::size_t>(i)] = yCloud.pts[i].z;
  }
  box0_ = {yCloud.box[0], yCloud.box[1], yCloud.box[2]};
}

void nneigh::SkinNeighborList::refreshBonds(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  std::set<std::pair<int, int>> next;
  for (const auto &[i, j] : candidates_) {
    if (i < 0 || j < 0 || i >= yCloud.nop || j >= yCloud.nop) {
      continue;
    }
    const double r2 = gen::periodicDistSq(yCloud, i, j);
    if (r2 <= cutoffSq_) {
      next.emplace(i, j);
    }
  }
  std::set<int> touched;
  for (const auto &p : bonded_) {
    if (next.count(p) == 0) {
      touched.insert(p.first);
      touched.insert(p.second);
    }
  }
  for (const auto &p : next) {
    if (bonded_.count(p) == 0) {
      touched.insert(p.first);
      touched.insert(p.second);
    }
  }
  changedAtoms_ = static_cast<int>(touched.size());
  bonded_.swap(next);

  const auto indexToID = indexToIDTable(yCloud);
  nList_ = seedWithSelfIDs(indexToID, yCloud.nop);
  for (const auto &[i, j] : bonded_) {
    appendNeighbourID(nList_, indexToID, i, j);
    appendNeighbourID(nList_, indexToID, j, i);
  }
}

const std::vector<std::vector<int>> &nneigh::SkinNeighborList::update(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  rebuilt_ = mustRebuild(yCloud);
  if (rebuilt_) {
    rebuildCandidates(yCloud);
  }
  refreshBonds(yCloud);
  return nList_;
}
