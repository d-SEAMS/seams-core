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
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>

#include <neighbours.hpp>
#include <simd_distance.hpp>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

#ifdef SEAMS_HAS_VESIN
#include <vesin.h>
#endif
#ifdef SEAMS_HAS_LINKCELL
#include <linkcell.hpp>
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

// Ortho SIMD wrap on box[0..2] is wrong for a tilt dump: those
// entries are bound spans, not lx, ly, lz. gen::periodicDistSq
// already inverts the dump 3x3 when box.size() >= 6.
void fillPairDistSq(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int iatom,
    const int *jatom, std::size_t n, double *dx, double *dy, double *dz,
    double *distSq) {
  if (yCloud.box.size() >= 6) {
    gen::batchPeriodicDistSq(yCloud, iatom, jatom, n, distSq);
    return;
  }
  const double xi = yCloud.pts[static_cast<std::size_t>(iatom)].x;
  const double yi = yCloud.pts[static_cast<std::size_t>(iatom)].y;
  const double zi = yCloud.pts[static_cast<std::size_t>(iatom)].z;
  for (std::size_t k = 0; k < n; k++) {
    const auto &pj = yCloud.pts[static_cast<std::size_t>(jatom[k])];
    dx[k] = xi - pj.x;
    dy[k] = yi - pj.y;
    dz[k] = zi - pj.z;
  }
  seams::BatchPeriodicDistSq(dx, dy, dz, yCloud.box[0], yCloud.box[1],
                             yCloud.box[2], distSq, n);
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
  if (nSubset == 0) {
    pairs.clear();
    return true;
  }
  std::vector<std::array<double, 3>> positions(nSubset);
  for (size_t i = 0; i < nSubset; i++) {
    const int idx = subset[i];
    positions[i] = {yCloud.pts[idx].x, yCloud.pts[idx].y, yCloud.pts[idx].z};
  }

  double box[3][3];
  double origin[3];
  nneigh::dumpBoundsToH(yCloud.box, yCloud.boxLow, box, origin);
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
 *  particle. Vesin cell list when built with SEAMS_HAS_VESIN, else
 *  brute-force \f$ O(n^2) \f$. This generates the full neighbour list, by ID.
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
  return neighListPair(rcutoff, yCloud, typeI, typeJ);
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

namespace {
constexpr int kThreadedCellRowsMinAtoms = 2048;
} // namespace

bool nneigh::cellListRowsThreaded(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &subset, double rcutoff,
    std::vector<std::vector<int>> &rows) {
  const std::size_t n = subset.size();
  if (n == 0) {
    rows.clear();
    return true;
  }
  if (rcutoff <= 0.0 || !hasPeriodicBox(yCloud)) {
    return false;
  }
  double H[3][3];
  double origin[3];
  nneigh::dumpBoundsToH(yCloud.box, yCloud.boxLow, H, origin);
  const double lx = H[0][0];
  const double ly = H[1][1];
  const double lz = H[2][2];
  const double xy = H[1][0];
  const double xz = H[2][0];
  const double yz = H[2][1];
  if (lx <= 0.0 || ly <= 0.0 || lz <= 0.0) {
    return false;
  }
  // perpendicular widths of the cell along a, b and c: V / |b x c| etc.
  const double volume = lx * ly * lz;
  const auto norm3 = [](double x, double y, double z) {
    return std::sqrt(x * x + y * y + z * z);
  };
  // a = (lx, 0, 0), b = (xy, ly, 0), c = (xz, yz, lz)
  const double bxc = norm3(ly * lz, -xy * lz, xy * yz - ly * xz);
  const double cxa = norm3(0.0, xz * 0.0 + lz * lx, -yz * lx);
  const double axb = norm3(0.0, 0.0, lx * ly);
  const double width[3] = {volume / bxc, volume / cxa, volume / axb};
  int ncell[3];
  for (int k = 0; k < 3; k++) {
    ncell[k] = static_cast<int>(std::floor(width[k] / rcutoff));
    if (ncell[k] < 3) {
      return false;
    }
  }
  // fractional coordinates of the subset, wrapped into [0, 1)
  std::vector<int> cellOf(n);
  const auto wrap = [](double f) {
    f -= std::floor(f);
    return f >= 1.0 ? 0.0 : f;
  };
  for (std::size_t k = 0; k < n; k++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(subset[k])];
    const double rx = p.x - origin[0];
    const double ry = p.y - origin[1];
    const double rz = p.z - origin[2];
    const double fc = rz / lz;
    const double fb = (ry - fc * yz) / ly;
    const double fa = (rx - fb * xy - fc * xz) / lx;
    const double f[3] = {wrap(fa), wrap(fb), wrap(fc)};
    int c[3];
    for (int d = 0; d < 3; d++) {
      c[d] = std::min(ncell[d] - 1, static_cast<int>(f[d] * ncell[d]));
    }
    cellOf[k] = (c[0] * ncell[1] + c[1]) * ncell[2] + c[2];
  }
  const int nCells = ncell[0] * ncell[1] * ncell[2];
  std::vector<int> cellStart(static_cast<std::size_t>(nCells) + 1, 0);
  for (std::size_t k = 0; k < n; k++) {
    ++cellStart[static_cast<std::size_t>(cellOf[k]) + 1];
  }
  for (int cidx = 0; cidx < nCells; cidx++) {
    cellStart[static_cast<std::size_t>(cidx) + 1] += cellStart[static_cast<std::size_t>(cidx)];
  }
  std::vector<int> cellAtoms(n);
  {
    std::vector<int> fill(cellStart.begin(), cellStart.end() - 1);
    for (std::size_t k = 0; k < n; k++) {
      cellAtoms[static_cast<std::size_t>(fill[static_cast<std::size_t>(cellOf[k])]++)] =
          static_cast<int>(k);
    }
  }
  rows.assign(n, {});
  const double rc2 = rcutoff * rcutoff;
  // each row belongs to one thread; the 27 cells around an atom hold the
  // nearest image of every neighbour because every axis has at least
  // three cells of width rcutoff
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(dynamic, 256) if (n >= 4096)
#endif
  for (std::int64_t kk = 0; kk < static_cast<std::int64_t>(n); kk++) {
    const std::size_t k = static_cast<std::size_t>(kk);
    const int i = subset[k];
    const int cid = cellOf[k];
    const int cz = cid % ncell[2];
    const int cy = (cid / ncell[2]) % ncell[1];
    const int cx = cid / (ncell[1] * ncell[2]);
    auto &row = rows[k];
    for (int dx = -1; dx <= 1; dx++) {
      const int nx = (cx + dx + ncell[0]) % ncell[0];
      for (int dy = -1; dy <= 1; dy++) {
        const int ny = (cy + dy + ncell[1]) % ncell[1];
        for (int dz = -1; dz <= 1; dz++) {
          const int nz = (cz + dz + ncell[2]) % ncell[2];
          const int ncid = (nx * ncell[1] + ny) * ncell[2] + nz;
          for (int m = cellStart[static_cast<std::size_t>(ncid)];
               m < cellStart[static_cast<std::size_t>(ncid) + 1]; m++) {
            const int j = subset[static_cast<std::size_t>(cellAtoms[static_cast<std::size_t>(m)])];
            if (j == i) {
              continue;
            }
            if (gen::periodicDistSq(yCloud, i, j) <= rc2) {
              row.push_back(j);
            }
          }
        }
      }
    }
    std::sort(row.begin(), row.end());
  }
  return true;
}

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

#ifdef SEAMS_HAS_OPENMP
  // threaded cell list: one row per thread, the same minimum-image set
  if (typeIIndices.size() >= static_cast<std::size_t>(kThreadedCellRowsMinAtoms)) {
    std::vector<std::vector<int>> rows;
    if (cellListRowsThreaded(yCloud, typeIIndices, rcutoff, rows)) {
      nList = seedWithSelfIDs(indexToID, yCloud.nop);
      for (std::size_t k = 0; k < rows.size(); k++) {
        auto &dest = nList[static_cast<std::size_t>(typeIIndices[k])];
        for (const int j : rows[k]) {
          dest.push_back(indexToID[static_cast<std::size_t>(j)]);
        }
      }
      return nList;
    }
  }
#endif

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

    fillPairDistSq(yCloud, iatom, typeIIndices.data() + ii + 1, remaining,
                   dx.data(), dy.data(), dz.data(), distSq.data());

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
 * @details I-J neighbour list. Like-type (I==J) reuses neighListO,
 *  including its dump MIC. Unlike-type pairs use the same vesin /
 *  fillPairDistSq path, then keep only mixed I-J contacts.
 * @param[in] rcutoff Distance cutoff, within which two atoms are neighbours.
 * @param[in] yCloud The input molSys::PointCloud
 * @param[in] typeI Type ID of particles of type I.
 * @param[in] typeJ Type ID of particles of type J.
 * @return Row-ordered full neighbour list, by atom ID.
 */
std::vector<std::vector<int>>
nneigh::neighListPair(
    double rcutoff,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeI,
    int typeJ) {
  if (typeI == typeJ) {
    return neighListO(rcutoff, yCloud, typeI);
  }
  if (!hasPeriodicBox(yCloud)) {
    return {};
  }

  const std::vector<int> indexToID = indexToIDTable(yCloud);
  std::vector<std::vector<int>> nList = seedWithSelfIDs(indexToID, yCloud.nop);

#ifdef SEAMS_HAS_VESIN
  {
    std::vector<int> subset;
    subset.reserve(static_cast<size_t>(yCloud.nop));
    for (int i = 0; i < yCloud.nop; i++) {
      const int t = yCloud.pts[i].type;
      if (t == typeI || t == typeJ) {
        subset.push_back(i);
      }
    }
    std::vector<std::pair<int, int>> pairs;
    if (cellListPairs(yCloud, subset, rcutoff, pairs)) {
      for (const auto &[iatom, jatom] : pairs) {
        const int ti = yCloud.pts[iatom].type;
        const int tj = yCloud.pts[jatom].type;
        const bool mixed = (ti == typeI && tj == typeJ) ||
                           (ti == typeJ && tj == typeI);
        if (!mixed) {
          continue;
        }
        appendNeighbourID(nList, indexToID, iatom, jatom);
      }
      return nList;
    }
  }
#endif

  const double rcutoffSq = rcutoff * rcutoff;
  std::vector<int> iIdx;
  std::vector<int> jIdx;
  iIdx.reserve(static_cast<size_t>(yCloud.nop));
  jIdx.reserve(static_cast<size_t>(yCloud.nop));
  for (int i = 0; i < yCloud.nop; i++) {
    const int t = yCloud.pts[i].type;
    if (t == typeI) {
      iIdx.push_back(i);
    }
    if (t == typeJ) {
      jIdx.push_back(i);
    }
  }

  std::vector<double> dx(jIdx.size()), dy(jIdx.size()), dz(jIdx.size());
  std::vector<double> distSq(jIdx.size());
  for (const int iatom : iIdx) {
    if (jIdx.empty()) {
      break;
    }
    fillPairDistSq(yCloud, iatom, jIdx.data(), jIdx.size(), dx.data(),
                   dy.data(), dz.data(), distSq.data());
    for (size_t jj = 0; jj < jIdx.size(); jj++) {
      if (distSq[jj] > rcutoffSq) {
        continue;
      }
      const int jatom = jIdx[jj];
      if (iatom == jatom) {
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
 *  particle of only one type. Vesin cell list when built with
 *  SEAMS_HAS_VESIN, else brute-force \f$ O(n^2) \f$. This generates the
 *  half neighbour list, by ID. This function will only work for building a
 *  neighbour list between one type of particles.
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

#ifdef SEAMS_HAS_VESIN
  {
    std::vector<int> typeIIndices;
    for (int i = 0; i < yCloud.nop; i++) {
      if (yCloud.pts[i].type == typeI) {
        typeIIndices.push_back(i);
      }
    }
    std::vector<std::pair<int, int>> pairs;
    if (cellListPairs(yCloud, typeIIndices, rcutoff, pairs)) {
      // Half list: keep the lower cloud index as the owner, matching the
      // brute walk that only writes j into i when i < j.
      for (const auto &[iatom, jatom] : pairs) {
        if (iatom < jatom) {
          appendNeighbourID(nList, indexToID, iatom, jatom);
        }
      }
      return nList;
    }
  }
#endif

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
#ifdef SEAMS_HAS_OPENMP
  if (yCloud.nop >= kThreadedCellRowsMinAtoms) {
    std::vector<int> allIndices(static_cast<std::size_t>(yCloud.nop));
    std::iota(allIndices.begin(), allIndices.end(), 0);
    std::vector<std::vector<int>> rows;
    if (cellListRowsThreaded(yCloud, allIndices, cutoff, rows)) {
      for (std::size_t k = 0; k < rows.size(); k++) {
        nList[k].insert(nList[k].end(), rows[k].begin(), rows[k].end());
      }
      return nList;
    }
  }
#endif
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

  // Scratch buffers for the batched distance kernel, sized once for the
  // largest batch and reused across iatom
  std::vector<double> dx(yCloud.nop), dy(yCloud.nop), dz(yCloud.nop);
  std::vector<double> distSq(yCloud.nop);
  std::vector<int> allIdx(static_cast<std::size_t>(yCloud.nop));
  std::iota(allIdx.begin(), allIdx.end(), 0);

  // Loop through every iatom and find nearest neighbours within cutoff
  for (int iatom = 0; iatom < yCloud.nop - 1; iatom++) {
    const size_t remaining = static_cast<size_t>(yCloud.nop - iatom - 1);

    fillPairDistSq(yCloud, iatom, allIdx.data() + iatom + 1, remaining,
                   dx.data(), dy.data(), dz.data(), distSq.data());

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

  // rows are independent and the map is only read, so each row is one task
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nList.size() >= 4096)
#endif
  for (std::int64_t ii = 0; ii < static_cast<std::int64_t>(nList.size()); ii++) {
    const std::size_t iatom = static_cast<std::size_t>(ii);
    if (nList[iatom].empty()) {
      continue;
    }
    const int iatomID = nList[iatom][0];
    auto gotI = yCloud.idIndexMap.find(iatomID);
    if (gotI == yCloud.idIndexMap.end()) {
      continue;
    }
    auto &out = indexNlist[iatom];
    out.reserve(nList[iatom].size());
    out.push_back(gotI->second);
    for (size_t j = 1; j < nList[iatom].size(); j++) {
      auto gotJ = yCloud.idIndexMap.find(nList[iatom][j]);
      if (gotJ == yCloud.idIndexMap.end()) {
        continue;
      }
      out.push_back(gotJ->second);
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

// k-nearest nominations: linkcell mesh, not vesin occupancy.
// candidateCutoff is a cell-size hint for knearest_into.
namespace {

#ifdef SEAMS_HAS_LINKCELL
// Writes n*k cloud indices, nearest first, unused slots -1.
bool nominateByLinkcell(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    int typeI, double cellHint, std::vector<int> &nom) {
  const int n = yCloud.nop;
  if (n <= 0 || k <= 0 || yCloud.box.size() < 3) {
    return true;
  }
  const double bx = yCloud.box[0];
  const double by = yCloud.box[1];
  const double bz = yCloud.box[2];
  if (!(bx > 0.0 && by > 0.0 && bz > 0.0)) {
    return true;
  }
  std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
  std::vector<int> mask(static_cast<std::size_t>(n), 0);
  int nSrc = 0;
  for (int i = 0; i < n; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(i) * 3 + 0] = p.x;
    xyz[static_cast<std::size_t>(i) * 3 + 1] = p.y;
    xyz[static_cast<std::size_t>(i) * 3 + 2] = p.z;
    if (p.type == typeI) {
      mask[static_cast<std::size_t>(i)] = 1;
      ++nSrc;
    }
  }
  if (nSrc == 0) {
    return true;
  }
  const linkcell::Cell cell =
      nneigh::lammpsBoxToLinkcell(yCloud.box, yCloud.boxLow);
  try {
    linkcell::knearest_into(xyz.data(), static_cast<std::size_t>(n), cell,
                            static_cast<std::size_t>(k), nom.data(),
                            nom.size(), mask.data(), cellHint);
  } catch (const linkcell::Error &e) {
    std::cerr << "linkcell failed: " << e.what()
              << "; falling back to the in-tree cell list.\n";
    return false;
  }
  return true;
}
#endif

// Packed n*k nominations, nearest first, unused slots -1.
std::vector<int> nominatePacked(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    int typeI, double cellHint) {
  const int n = yCloud.nop;
  std::vector<int> nom(static_cast<std::size_t>(n) * static_cast<std::size_t>(k),
                       -1);
  if (n <= 0 || k <= 0 || yCloud.box.size() < 3) {
    return nom;
  }
#ifdef SEAMS_HAS_LINKCELL
  if (nominateByLinkcell(yCloud, k, typeI, cellHint, nom)) {
    return nom;
  }
#endif
  std::vector<int> owners;
  owners.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    if (yCloud.pts[static_cast<std::size_t>(i)].type == typeI) {
      owners.push_back(i);
    }
  }
  if (yCloud.box.size() >= 6) {
    for (const int i : owners) {
      std::priority_queue<std::pair<double, int>> heap;
      for (const int j : owners) {
        if (j == i) {
          continue;
        }
        const double d2 = gen::periodicDistSq(yCloud, i, j);
        if (static_cast<int>(heap.size()) < k) {
          heap.push({d2, j});
        } else if (d2 < heap.top().first) {
          heap.pop();
          heap.push({d2, j});
        }
      }
      const int take = std::min(k, static_cast<int>(heap.size()));
      for (int t = take - 1; t >= 0; t--) {
        nom[static_cast<std::size_t>(i) * static_cast<std::size_t>(k) +
            static_cast<std::size_t>(t)] = heap.top().second;
        heap.pop();
      }
    }
    return nom;
  }
  const double bx = yCloud.box[0];
  const double by = yCloud.box[1];
  const double bz = yCloud.box[2];
  if (!(bx > 0.0 && by > 0.0 && bz > 0.0)) {
    return nom;
  }
  const double xlo = yCloud.boxLow.size() >= 3 ? yCloud.boxLow[0] : 0.0;
  const double ylo = yCloud.boxLow.size() >= 3 ? yCloud.boxLow[1] : 0.0;
  const double zlo = yCloud.boxLow.size() >= 3 ? yCloud.boxLow[2] : 0.0;
  double edge = cellHint > 0.0 ? cellHint : 3.0;
  edge = std::min({edge, bx, by, bz});
  int nx = std::max(1, static_cast<int>(std::floor(bx / edge)));
  int ny = std::max(1, static_cast<int>(std::floor(by / edge)));
  int nz = std::max(1, static_cast<int>(std::floor(bz / edge)));
  const int ncell = nx * ny * nz;
  const double invx = static_cast<double>(nx) / bx;
  const double invy = static_cast<double>(ny) / by;
  const double invz = static_cast<double>(nz) / bz;
  const double cellMin = std::min({bx / nx, by / ny, bz / nz});

  std::vector<int> head(static_cast<std::size_t>(ncell), -1);
  std::vector<int> next(static_cast<std::size_t>(yCloud.nop), -1);

  auto wrap = [](double x, double L) {
    double t = x / L;
    t -= std::floor(t);
    return t * L;
  };
  auto cellOf = [&](int iatom, int &ix, int &iy, int &iz) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(iatom)];
    const double x = wrap(p.x - xlo, bx);
    const double y = wrap(p.y - ylo, by);
    const double z = wrap(p.z - zlo, bz);
    ix = std::min(nx - 1, std::max(0, static_cast<int>(x * invx)));
    iy = std::min(ny - 1, std::max(0, static_cast<int>(y * invy)));
    iz = std::min(nz - 1, std::max(0, static_cast<int>(z * invz)));
  };
  auto cellIndex = [&](int ix, int iy, int iz) {
    ix %= nx;
    if (ix < 0) {
      ix += nx;
    }
    iy %= ny;
    if (iy < 0) {
      iy += ny;
    }
    iz %= nz;
    if (iz < 0) {
      iz += nz;
    }
    return (iz * ny + iy) * nx + ix;
  };

  for (const int i : owners) {
    int ix = 0;
    int iy = 0;
    int iz = 0;
    cellOf(i, ix, iy, iz);
    const int c = cellIndex(ix, iy, iz);
    next[static_cast<std::size_t>(i)] = head[static_cast<std::size_t>(c)];
    head[static_cast<std::size_t>(c)] = i;
  }

  const int maxReach = std::max({nx, ny, nz}) / 2 + 1;
  std::vector<unsigned> seen(static_cast<std::size_t>(ncell), 0);
  unsigned stamp = 1;
  for (const int i : owners) {
    std::priority_queue<std::pair<double, int>> heap;
    int ix = 0;
    int iy = 0;
    int iz = 0;
    cellOf(i, ix, iy, iz);
    ++stamp;
    if (stamp == 0) {
      std::fill(seen.begin(), seen.end(), 0);
      stamp = 1;
    }
    int reach = 1;
    while (reach <= maxReach) {
      for (int dx = -reach; dx <= reach; dx++) {
        for (int dy = -reach; dy <= reach; dy++) {
          for (int dz = -reach; dz <= reach; dz++) {
            const bool shell = (reach == 1) || (std::abs(dx) == reach ||
                                                std::abs(dy) == reach ||
                                                std::abs(dz) == reach);
            if (!shell && reach > 1) {
              continue;
            }
            const int c = cellIndex(ix + dx, iy + dy, iz + dz);
            if (seen[static_cast<std::size_t>(c)] == stamp) {
              continue;
            }
            seen[static_cast<std::size_t>(c)] = stamp;
            int j = head[static_cast<std::size_t>(c)];
            while (j >= 0) {
              if (j != i) {
                const double d2 = gen::periodicDistSq(yCloud, i, j);
                if (static_cast<int>(heap.size()) < k) {
                  heap.push({d2, j});
                } else if (d2 < heap.top().first) {
                  heap.pop();
                  heap.push({d2, j});
                }
              }
              j = next[static_cast<std::size_t>(j)];
            }
          }
        }
      }
      if (static_cast<int>(heap.size()) >= k) {
        const double bound = static_cast<double>(reach) * cellMin;
        if (heap.top().first <= bound * bound) {
          break;
        }
      }
      ++reach;
    }
    const int take = std::min(k, static_cast<int>(heap.size()));
    for (int t = take - 1; t >= 0; t--) {
      nom[static_cast<std::size_t>(i) * static_cast<std::size_t>(k) +
          static_cast<std::size_t>(t)] = heap.top().second;
      heap.pop();
    }
  }
  return nom;
}

bool listsPacked(const std::vector<int> &nom, int n, int k, int i, int j) {
  if (i < 0 || i >= n) {
    return false;
  }
  const std::size_t base =
      static_cast<std::size_t>(i) * static_cast<std::size_t>(k);
  for (int t = 0; t < k; t++) {
    if (nom[base + static_cast<std::size_t>(t)] == j) {
      return true;
    }
  }
  return false;
}

void addUndirectedIDs(
    std::vector<std::vector<int>> &out,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int i,
    int j) {
  out[static_cast<std::size_t>(i)].push_back(
      yCloud.pts[static_cast<std::size_t>(j)].atomID);
  out[static_cast<std::size_t>(j)].push_back(
      yCloud.pts[static_cast<std::size_t>(i)].atomID);
}

void fillGraphFromNom(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &nom, int k, bool mutual,
    std::vector<std::vector<int>> &out) {
  const int n = yCloud.nop;
  std::vector<int> deg(static_cast<std::size_t>(n), 0);
  for (int i = 0; i < n; i++) {
    for (int t = 0; t < k; t++) {
      const int j =
          nom[static_cast<std::size_t>(i) * static_cast<std::size_t>(k) +
              static_cast<std::size_t>(t)];
      if (j < 0 || j >= n || j == i) {
        continue;
      }
      const bool listed = listsPacked(nom, n, k, j, i);
      if (mutual) {
        if (i < j && listed) {
          deg[static_cast<std::size_t>(i)]++;
          deg[static_cast<std::size_t>(j)]++;
        }
      } else if (i < j || !listed) {
        deg[static_cast<std::size_t>(i)]++;
        deg[static_cast<std::size_t>(j)]++;
      }
    }
  }
  out.assign(static_cast<std::size_t>(n), {});
  for (int i = 0; i < n; i++) {
    out[static_cast<std::size_t>(i)].reserve(
        static_cast<std::size_t>(1 + deg[static_cast<std::size_t>(i)]));
    out[static_cast<std::size_t>(i)].push_back(
        yCloud.pts[static_cast<std::size_t>(i)].atomID);
  }
  for (int i = 0; i < n; i++) {
    for (int t = 0; t < k; t++) {
      const int j =
          nom[static_cast<std::size_t>(i) * static_cast<std::size_t>(k) +
              static_cast<std::size_t>(t)];
      if (j < 0 || j >= n || j == i) {
        continue;
      }
      const bool listed = listsPacked(nom, n, k, j, i);
      if (mutual) {
        if (i < j && listed) {
          addUndirectedIDs(out, yCloud, i, j);
        }
      } else if (i < j) {
        addUndirectedIDs(out, yCloud, i, j);
      } else if (!listed) {
        addUndirectedIDs(out, yCloud, i, j);
      }
    }
  }
}

void fillBothGraphsFromNom(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &nom, int k, std::vector<std::vector<int>> &mutual,
    std::vector<std::vector<int>> &uni) {
  fillGraphFromNom(yCloud, nom, k, true, mutual);
  fillGraphFromNom(yCloud, nom, k, false, uni);
}

} // namespace

std::vector<std::vector<int>>
nneigh::kNearestNeighbourList(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, int typeI, bool mutual) {
  if (yCloud.nop <= 0 || k <= 0) {
    return {};
  }
  const double hint = candidateCutoff > 0.0 ? 0.75 * candidateCutoff : 3.0;
  const auto nom = nominatePacked(yCloud, k, typeI, hint);
  std::vector<std::vector<int>> out;
  fillGraphFromNom(yCloud, nom, k, mutual, out);
  return out;
}

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
nneigh::kNearestNeighbourPair(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, int typeI) {
  if (yCloud.nop <= 0 || k <= 0) {
    return {{}, {}};
  }
  const double hint = candidateCutoff > 0.0 ? 0.75 * candidateCutoff : 3.0;
  const auto nom = nominatePacked(yCloud, k, typeI, hint);
  std::vector<std::vector<int>> mutual;
  std::vector<std::vector<int>> uni;
  fillBothGraphsFromNom(yCloud, nom, k, mutual, uni);
  return {std::move(mutual), std::move(uni)};
}

/**
 * @details Nearest unlike neighbour of every typeI index among typeJ
 *  indices. Each typeI walks every typeJ with gen::periodicDistSq
 *  (the dump MIC). Ties keep the lowest typeJ index. Particles with
 *  no unlike partner are omitted. Results are cloud indices and the
 *  MIC distance, not a coordination number.
 */
std::vector<std::tuple<int, int, double>> nneigh::nearestUnlike(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeI,
    int typeJ) {
  std::vector<std::tuple<int, int, double>> out;
  if (yCloud.box.size() < 3 || yCloud.nop <= 0) {
    return out;
  }

  std::vector<int> jIdx;
  jIdx.reserve(static_cast<std::size_t>(yCloud.nop));
  for (int j = 0; j < yCloud.nop; j++) {
    if (yCloud.pts[static_cast<std::size_t>(j)].type == typeJ) {
      jIdx.push_back(j);
    }
  }
  if (jIdx.empty()) {
    return out;
  }

  for (int i = 0; i < yCloud.nop; i++) {
    if (yCloud.pts[static_cast<std::size_t>(i)].type != typeI) {
      continue;
    }
    int bestJ = -1;
    double bestD2 = std::numeric_limits<double>::infinity();
    for (const int j : jIdx) {
      if (j == i) {
        continue;
      }
      const double d2 = gen::periodicDistSq(yCloud, i, j);
      if (d2 < bestD2) {
        bestD2 = d2;
        bestJ = j;
      }
    }
    if (bestJ >= 0) {
      out.emplace_back(i, bestJ, std::sqrt(bestD2));
    }
  }
  return out;
}

/**
 * @details Contact pairs: those nearestUnlike edges whose typeJ
 *  partner also lists the typeI index as its nearest unlike
 *  neighbour. Not first-shell membership.
 */
std::vector<std::pair<int, int>> nneigh::mutualNearestUnlike(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeI,
    int typeJ) {
  const auto forward = nearestUnlike(yCloud, typeI, typeJ);
  const auto reverse = nearestUnlike(yCloud, typeJ, typeI);
  std::vector<int> partnerOfJ(static_cast<std::size_t>(yCloud.nop), -1);
  for (const auto &row : reverse) {
    const int j = std::get<0>(row);
    if (j >= 0 && j < yCloud.nop) {
      partnerOfJ[static_cast<std::size_t>(j)] = std::get<1>(row);
    }
  }
  std::vector<std::pair<int, int>> out;
  for (const auto &row : forward) {
    const int i = std::get<0>(row);
    const int j = std::get<1>(row);
    if (j >= 0 && j < yCloud.nop &&
        partnerOfJ[static_cast<std::size_t>(j)] == i) {
      out.emplace_back(i, j);
    }
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

nneigh::BondGraph nneigh::bondGraphFromName(const std::string &name) {
  if (name == "cutoff") {
    return BondGraph::Cutoff;
  }
  if (name == "knn" || name == "knn-mutual" || name == "mutual") {
    return BondGraph::KnnMutual;
  }
  if (name == "knn-union" || name == "union") {
    return BondGraph::KnnUnion;
  }
  throw std::invalid_argument(
      "unknown bond graph '" + name +
      "' (use cutoff, knn, knn-union)");
}

const char *nneigh::bondGraphName(BondGraph graph) {
  switch (graph) {
  case BondGraph::Cutoff:
    return "cutoff";
  case BondGraph::KnnMutual:
    return "knn";
  case BondGraph::KnnUnion:
    return "knn-union";
  }
  return "cutoff";
}

nneigh::SkinNeighborList::SkinNeighborList(double cutoff, double skin,
                                           int typeI, BondGraph graph, int k)
    : cutoff_(cutoff), skin_(skin), cutoffSq_(cutoff * cutoff),
      triggerSq_((0.5 * skin) * (0.5 * skin)), typeI_(typeI),
      k_(graph == BondGraph::Cutoff ? 0 : k), graph_(graph),
      mutual_(graph != BondGraph::KnnUnion) {}

bool nneigh::SkinNeighborList::mustRebuild(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) const {
  if (static_cast<int>(x0_.size()) != yCloud.nop || yCloud.box.size() < 3) {
    return true;
  }
  if (yCloud.box.size() != box0_.size()) {
    return true;
  }
  for (std::size_t k = 0; k < yCloud.box.size(); k++) {
    if (std::fabs(yCloud.box[k] - box0_[k]) > 1e-8) {
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
  box0_ = yCloud.box;
}

void nneigh::SkinNeighborList::refreshBonds(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  std::vector<std::pair<int, int>> next;
  next.reserve(candidates_.size());
  const std::size_t m = candidates_.size();
  std::vector<double> dx(m), dy(m), dz(m), r2(m);
  for (std::size_t t = 0; t < m; t++) {
    const int i = candidates_[t].first;
    const int j = candidates_[t].second;
    if (i < 0 || j < 0 || i >= yCloud.nop || j >= yCloud.nop) {
      dx[t] = dy[t] = dz[t] = 0.0;
      continue;
    }
    dx[t] = yCloud.pts[static_cast<std::size_t>(j)].x -
            yCloud.pts[static_cast<std::size_t>(i)].x;
    dy[t] = yCloud.pts[static_cast<std::size_t>(j)].y -
            yCloud.pts[static_cast<std::size_t>(i)].y;
    dz[t] = yCloud.pts[static_cast<std::size_t>(j)].z -
            yCloud.pts[static_cast<std::size_t>(i)].z;
  }
  if (yCloud.box.size() >= 6 && m > 0) {
    for (std::size_t t = 0; t < m; t++) {
      const int i = candidates_[t].first;
      const int j = candidates_[t].second;
      if (i < 0 || j < 0 || i >= yCloud.nop || j >= yCloud.nop) {
        r2[t] = 0.0;
        continue;
      }
      r2[t] = gen::periodicDistSq(yCloud, i, j);
    }
  } else if (yCloud.box.size() >= 3 && m > 0) {
    seams::BatchPeriodicDistSq(dx.data(), dy.data(), dz.data(), yCloud.box[0],
                               yCloud.box[1], yCloud.box[2], r2.data(), m);
  }
  if (k_ <= 0) {
    for (std::size_t t = 0; t < m; t++) {
      const int i = candidates_[t].first;
      const int j = candidates_[t].second;
      if (i < 0 || j < 0 || i >= yCloud.nop || j >= yCloud.nop) {
        continue;
      }
      if (r2[t] <= cutoffSq_) {
        next.emplace_back(i, j);
      }
    }
  } else {
    std::vector<std::vector<std::pair<double, int>>> byDist(
        static_cast<std::size_t>(yCloud.nop));
    for (std::size_t t = 0; t < m; t++) {
      const int i = candidates_[t].first;
      const int j = candidates_[t].second;
      if (i < 0 || j < 0 || i >= yCloud.nop || j >= yCloud.nop) {
        continue;
      }
      byDist[static_cast<std::size_t>(i)].emplace_back(r2[t], j);
      byDist[static_cast<std::size_t>(j)].emplace_back(r2[t], i);
    }
    std::vector<std::vector<int>> nominated(static_cast<std::size_t>(yCloud.nop));
    for (int i = 0; i < yCloud.nop; i++) {
      auto &row = byDist[static_cast<std::size_t>(i)];
      if (row.empty()) {
        continue;
      }
      const int take = std::min(k_, static_cast<int>(row.size()));
      std::partial_sort(row.begin(), row.begin() + take, row.end());
      nominated[static_cast<std::size_t>(i)].reserve(static_cast<std::size_t>(take));
      for (int n = 0; n < take; n++) {
        nominated[static_cast<std::size_t>(i)].push_back(row[static_cast<std::size_t>(n)].second);
      }
    }
    for (int i = 0; i < yCloud.nop; i++) {
      for (const int j : nominated[static_cast<std::size_t>(i)]) {
        const int a = std::min(i, j);
        const int b = std::max(i, j);
        if (mutual_) {
          const auto &other = nominated[static_cast<std::size_t>(j)];
          if (std::find(other.begin(), other.end(), i) == other.end()) {
            continue;
          }
        }
        next.emplace_back(a, b);
      }
    }
    std::sort(next.begin(), next.end());
    next.erase(std::unique(next.begin(), next.end()), next.end());
  }
  std::sort(next.begin(), next.end());
  next.erase(std::unique(next.begin(), next.end()), next.end());
  std::vector<int> touched;
  {
    std::size_t a = 0;
    std::size_t b = 0;
    while (a < bonded_.size() || b < next.size()) {
      if (b == next.size() || (a < bonded_.size() && bonded_[a] < next[b])) {
        touched.push_back(bonded_[a].first);
        touched.push_back(bonded_[a].second);
        ++a;
      } else if (a == bonded_.size() || next[b] < bonded_[a]) {
        touched.push_back(next[b].first);
        touched.push_back(next[b].second);
        ++b;
      } else {
        ++a;
        ++b;
      }
    }
    std::sort(touched.begin(), touched.end());
    touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
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
