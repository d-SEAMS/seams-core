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

#include <index_range.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include <neighbours.hpp>
#include <region_errors.hpp>
#include <simd_distance.hpp>
#include <thread_team.hpp>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

#ifdef SEAMS_HAS_MINIMAGE
#include <minimage.h>
#endif
#include <linkcell.hpp>

namespace {

// Set while a caller wants the one-thread cutoff list and the OpenMP
// team left untouched for a later parallel region.
thread_local bool gSerialCutoff = false;

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
    double H[3][3];
    double origin[3];
    nneigh::dumpBoundsToH(yCloud.box, yCloud.boxLow, H, origin);
    if (H[0][0] > 0.0 && H[1][1] > 0.0 && H[2][2] > 0.0) {
      return true;
    }
    std::cerr << "The simulation box has a non-positive edge; cannot build a "
                 "neighbour list.\n";
    return false;
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
  std::vector<int> indexToID(static_cast<std::size_t>(yCloud.nop), -1);
  for (const auto &[atomID, index] : yCloud.idIndexMap) {
    if (index >= 0 && index < yCloud.nop) {
      indexToID[static_cast<std::size_t>(index)] = atomID;
    }
  }
  return indexToID;
}

// Below this many atoms a row-by-row pass is cheaper than waking a team.
constexpr int kParallelRowsMinAtoms = 4096;

// Cells of about this many atoms keep a tiny cutoff in a large box under
// linkcell's mesh limit instead of one cell per cutoff width.
constexpr double kAtomsPerLinkcellCell = 20.0;

//! Volume of the Cartesian extent linkcell bins on a cell with an open
//! axis: the kept atoms, folded into the box along the periodic axes
//! only, so an open axis spans wherever its atoms are. Each side is at
//! least @a floor.
double binnedVolume(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    const std::vector<double> &xyz, const std::vector<int> &mask,
                    const std::array<bool, 3> &periodic, double floor) {
  const bool anyPeriodic = periodic[0] || periodic[1] || periodic[2];
  double H[3][3] = {};
  double origin[3] = {};
  if (anyPeriodic) {
    nneigh::dumpBoundsToH(yCloud.box, yCloud.boxLow, H, origin);
  }
  double lo[3];
  double hi[3];
  for (std::size_t axis = 0; axis < 3; ++axis) {
    lo[axis] = std::numeric_limits<double>::infinity();
    hi[axis] = -lo[axis];
  }
  for (std::size_t i = 0; 3 * i < xyz.size(); ++i) {
    if (!mask.empty() && mask[i] == 0) {
      continue;
    }
    double r[3] = {xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]};
    if (anyPeriodic) {
      double s[3];
      nneigh::dumpToFrac(H, origin, r[0], r[1], r[2], s);
      for (std::size_t axis = 0; axis < 3; ++axis) {
        if (periodic[axis]) {
          s[axis] -= std::floor(s[axis]);
        }
      }
      nneigh::dumpFromFrac(H, origin, s, r);
    }
    for (std::size_t axis = 0; axis < 3; ++axis) {
      lo[axis] = std::min(lo[axis], r[axis]);
      hi[axis] = std::max(hi[axis], r[axis]);
    }
  }
  double volume = 1.0;
  for (std::size_t axis = 0; axis < 3; ++axis) {
    volume *= std::max(hi[axis] - lo[axis], floor);
  }
  return volume;
}

constexpr std::array<bool, 3> kPeriodic{true, true, true};

//! Width of the cutoff search and of the row fill for this call: one
//! thread inside a parallel region or under preferSerialCutoff, otherwise
//! the OpenMP team.
int cutoffWidth() {
#ifdef SEAMS_HAS_OPENMP
  if (gSerialCutoff || omp_in_parallel()) {
    return 1;
  }
  return std::max(1, omp_get_max_threads());
#else
  return 1;
#endif
}

//! fill(i) for every row, split over @a width threads once the cloud is
//! large enough to pay for the team. Each row is written by one thread.
template <typename Fill> void forEachRow(int nop, int width, Fill &&fill) {
  seams::RegionErrors errors;
  // nvc++ -mp=gpu SIGSEGVs host parallel-fors in libnvomp
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
#pragma omp parallel for schedule(static) num_threads(width) if (width > 1 && nop >= kParallelRowsMinAtoms)
#else
  (void)width;
#endif
  for (int i = 0; i < nop; i++) {
    errors.run([&] { fill(i); });
  }
  errors.rethrow();
}

//! One flag per atom, nonzero where keep(type) holds, for the linkcell
//! mask. Empty when every atom is kept, so the search takes no mask.
template <typename Keep>
std::vector<int> typeMask(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                          Keep &&keep) {
  std::vector<int> mask(static_cast<std::size_t>(std::max(yCloud.nop, 0)), 0);
  bool every = true;
  for (int i = 0; i < yCloud.nop; i++) {
    const bool kept = keep(yCloud.pts[static_cast<std::size_t>(i)].type);
    mask[static_cast<std::size_t>(i)] = kept ? 1 : 0;
    every = every && kept;
  }
  if (every) {
    mask.clear();
  }
  return mask;
}

//! Row @a i of @a csr by atom ID, self first, into @a row. An atom without
//! an ID leaves the row empty, and a partner without one is skipped.
void idRow(const linkcell::PairCsr &csr, const std::vector<int> &indexToID, int i,
           std::vector<int> &row) {
  row.clear();
  if (indexToID[static_cast<std::size_t>(i)] == -1) {
    return;
  }
  const int a = csr.offsets[static_cast<std::size_t>(i)];
  const int b = csr.offsets[static_cast<std::size_t>(i) + 1];
  row.reserve(static_cast<std::size_t>(1 + b - a));
  row.push_back(indexToID[static_cast<std::size_t>(i)]);
  for (int k = a; k < b; k++) {
    const int id = indexToID[static_cast<std::size_t>(csr.cols[static_cast<std::size_t>(k)])];
    if (id != -1) {
      row.push_back(id);
    }
  }
}

//! Row @a i of @a csr by cloud index, self first, into @a row.
void indexRow(const linkcell::PairCsr &csr, int i, std::vector<int> &row) {
  const int a = csr.offsets[static_cast<std::size_t>(i)];
  const int b = csr.offsets[static_cast<std::size_t>(i) + 1];
  row.resize(static_cast<std::size_t>(1 + b - a));
  row[0] = i;
  std::copy(csr.cols.begin() + a, csr.cols.begin() + b, row.begin() + 1);
}

//! Rows by atom ID for every atom of the cloud, each filled on the thread
//! that allocates it.
std::vector<std::vector<int>> idRows(const linkcell::PairCsr &csr,
                                     const std::vector<int> &indexToID, int nop,
                                     int width) {
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(nop));
  forEachRow(nop, width,
             [&](int i) { idRow(csr, indexToID, i, rows[static_cast<std::size_t>(i)]); });
  return rows;
}

/**
 * @details Seeds a neighbour list so that row @a i begins with the atom ID of
 *  particle @a i, which is the row-header convention used throughout d-SEAMS.
 * @param[in] indexToID Reverse index-to-ID table.
 * @param[in] nop Number of particles.
 * @return Neighbour list of @a nop rows, each holding its own atom ID.
 */
std::vector<std::vector<int>> seedWithSelfIDs(const std::vector<int> &indexToID,
                                              int nop) {
  std::vector<std::vector<int>> nList(static_cast<std::size_t>(nop));
  for (int iatom = 0; iatom < nop; iatom++) {
    if (indexToID[static_cast<std::size_t>(iatom)] == -1) {
      std::cerr << "Something is wrong with your idIndexMap!\n";
      continue;
    }
    nList[static_cast<std::size_t>(iatom)].push_back(indexToID[static_cast<std::size_t>(iatom)]);
  }
  return nList;
}

//! True when index has a real atom ID in the reverse table.
bool hasAtomID(const std::vector<int> &indexToID, int index) {
  return index >= 0 && static_cast<size_t>(index) < indexToID.size() &&
         indexToID[static_cast<std::size_t>(index)] != -1;
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
  nList[static_cast<std::size_t>(iatom)].push_back(indexToID[static_cast<std::size_t>(jatom)]);
}


} // namespace

void nneigh::linkcellCutoffCsr(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double rcutoff,
    bool half, const std::array<bool, 3> &periodic, const std::vector<int> &mask,
    int threads, linkcell::PairCsr &out) {
  if (!std::isfinite(rcutoff) || !(rcutoff > 0.0)) {
    throw std::invalid_argument("the neighbour cutoff must be finite and positive");
  }
  const int nop = std::max(yCloud.nop, 0);
  if (!mask.empty() && mask.size() != static_cast<std::size_t>(nop)) {
    throw std::invalid_argument("the neighbour mask needs one flag per atom");
  }
  const bool anyPeriodic = periodic[0] || periodic[1] || periodic[2];
  if (anyPeriodic && yCloud.box.size() < 3) {
    throw std::invalid_argument("a periodic axis needs the cloud's box");
  }
  const int width = std::max(1, threads);
  thread_local std::vector<double> scratch;
  std::vector<double> &xyz = scratch;
  xyz.resize(3 * static_cast<std::size_t>(nop));
  forEachRow(nop, width, [&](int i) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    xyz[3 * static_cast<std::size_t>(i)] = p.x;
    xyz[3 * static_cast<std::size_t>(i) + 1] = p.y;
    xyz[3 * static_cast<std::size_t>(i) + 2] = p.z;
  });
  // Every axis open needs no cell; linkcell ignores the one passed.
  const linkcell::Cell cell = anyPeriodic
                                  ? nneigh::lammpsBoxToLinkcell(yCloud.box, yCloud.boxLow)
                                  : linkcell::Cell::ortho(1.0, 1.0, 1.0);
  const std::size_t kept =
      mask.empty() ? static_cast<std::size_t>(nop)
                   : static_cast<std::size_t>(
                         std::count_if(mask.begin(), mask.end(), [](int v) { return v != 0; }));
  double hint = rcutoff;
  if (kept > 0) {
    const bool allPeriodic = periodic[0] && periodic[1] && periodic[2];
    const double volume = allPeriodic ? nneigh::dumpVolume(yCloud.box, yCloud.boxLow)
                                      : binnedVolume(yCloud, xyz, mask, periodic, rcutoff);
    if (std::isfinite(volume) && volume > 0.0) {
      hint = std::max(rcutoff,
                      std::cbrt(kAtomsPerLinkcellCell * volume / static_cast<double>(kept)));
    }
  }
  const double none = 0.0;
  lc_set_threads(width);
  seams::withTeamAsleep(width > 1 ? seams::teamWidth() : 1, [&] {
    linkcell::pairs_within_csr(xyz.empty() ? &none : xyz.data(),
                               static_cast<std::size_t>(nop), cell, rcutoff, out, half,
                               periodic, mask.empty() ? nullptr : mask.data(), hint);
  });
}

/**
 * @details Function for building neighbour lists for each
 *  particle, from linkcell's cutoff rows. This generates the full
 *  neighbour list, by ID.
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

namespace {
constexpr int kThreadedCellRowsParallelAtoms = 4096;

// One thread bins one slab of an orthorhombic box plus the cutoff halo,
// periodic images folded next to that slab. Partner indices are cloud
// indices. Returns false unless several threads are available.
bool orthoThreadSlabs(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                      const std::vector<int> &subset, double cutoff,
                      std::vector<std::vector<int>> &rows) {
  // A private grid per thread, covering one slab plus the cutoff halo,
  // matched the fused rows and still took longer than the shared cell
  // list on a filled box. The shared list is the threaded path.
  (void)yCloud;
  (void)subset;
  (void)cutoff;
  (void)rows;
  return false;
}


} // namespace

bool nneigh::cellListRowsThreaded(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &subset, double rcutoff,
    std::vector<std::vector<int>> &rows) {
  if (orthoThreadSlabs(yCloud, subset, rcutoff, rows)) {
    return true;
  }
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
  seams::coarsenUntilFits(ncell[0], ncell[1], ncell[2], 3);
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
  int nCells = 0;
  if (!seams::gridFits(ncell[0], ncell[1], ncell[2], nCells)) {
    return false;
  }
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
  for (auto &row : rows) {
    row.reserve(16);
  }
  const double rc2 = rcutoff * rcutoff;
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  // each row belongs to one thread; the 27 cells around an atom hold the
  // nearest image of every neighbour because every axis has at least
  // three cells of width rcutoff
  seams::RegionErrors errors;
#ifdef SEAMS_HAS_OPENMP
  // nvc++ -mp=gpu (SEAMS_HAS_OFFLOAD) SIGSEGVs this host parallel-for
  // in libnvomp (__kmpc_fork_call).
#if !defined(SEAMS_HAS_OFFLOAD)
#pragma omp parallel if (n >= static_cast<std::size_t>(kThreadedCellRowsParallelAtoms))
#endif
#endif
  {
    std::vector<int> cand;
#ifdef SEAMS_HAS_OPENMP
#if !defined(SEAMS_HAS_OFFLOAD)
#pragma omp for schedule(static)
#endif
#endif
    for (std::int64_t kk = 0; kk < static_cast<std::int64_t>(n); kk++) {
      errors.run([&] {
        const std::size_t k = static_cast<std::size_t>(kk);
        const int i = subset[k];
        const int cid = cellOf[k];
        const int cz = cid % ncell[2];
        const int cy = (cid / ncell[2]) % ncell[1];
        const int cx = cid / (ncell[1] * ncell[2]);
        cand.clear();
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
                if (j != i) {
                  cand.push_back(j);
                }
              }
            }
          }
        }
        auto &row = rows[k];
        for (std::size_t t = 0; t < cand.size(); t++) {
          const double r2 = gen::periodicDistSq(frame, yCloud, i, cand[t]);
          if (r2 <= rc2) {
            row.push_back(cand[t]);
          }
        }
        std::sort(row.begin(), row.end());
      });
    }
  }
  errors.rethrow();
  return true;
}

namespace {

struct alignas(64) PairBag {
  std::vector<std::pair<int, int>> pairs;
};

// Pairs inside one cutoff bin of the unwrapped window. Coordinates are
// stored in cell order so each thread walks a contiguous slab. The
// periodic axes are compile-time flags: a runtime branch in the inner
// loop kept the hybrid teams off the pure-MPI pace.
template <bool PX, bool PY, bool PZ>
void searchSlots(const double *cx, const double *cy, const double *cz, const int *slotAtom,
                 const int *cellStart, const int *cellOf, int slot0, int slot1, int nc0, int nc1,
                 int nc2, double cut2, double ex, double ey, double ez,
                 std::vector<std::pair<int, int>> &mine) {
  mine.clear();
  for (int slot = slot0; slot < slot1; ++slot) {
      const int i = slotAtom[slot];
      const double ixv = cx[slot];
      const double iyv = cy[slot];
      const double izv = cz[slot];
      const int home = cellOf[i];
      const int hx = home % nc0;
      const int hy = (home / nc0) % nc1;
      const int hz = home / (nc0 * nc1);
      // A non-periodic window has no image of a pair in another cell, so
      // the forward half of the 27-cell stencil lists each pair once.
      const int dz0 = (!PX && !PY && !PZ) ? 0 : -1;
      for (int dz = dz0; dz <= 1; ++dz) {
        int jz = hz + dz;
        if constexpr (PZ) {
          if (jz < 0) {
            jz += nc2;
          } else if (jz >= nc2) {
            jz -= nc2;
          }
        } else if (jz < 0 || jz >= nc2) {
          continue;
        }
        const int dy0 = (!PX && !PY && !PZ && dz == 0) ? 0 : -1;
        for (int dy = dy0; dy <= 1; ++dy) {
          int jy = hy + dy;
          if constexpr (PY) {
            if (jy < 0) {
              jy += nc1;
            } else if (jy >= nc1) {
              jy -= nc1;
            }
          } else if (jy < 0 || jy >= nc1) {
            continue;
          }
          const int dx0 = (!PX && !PY && !PZ && dz == 0 && dy == 0) ? 0 : -1;
          for (int dx = dx0; dx <= 1; ++dx) {
            int jx = hx + dx;
            if constexpr (PX) {
              if (jx < 0) {
                jx += nc0;
              } else if (jx >= nc0) {
                jx -= nc0;
              }
            } else if (jx < 0 || jx >= nc0) {
              continue;
            }
            const int c = (jz * nc1 + jy) * nc0 + jx;
            const int begin = cellStart[c];
            const int end = cellStart[c + 1];
            for (int t = begin; t < end; ++t) {
              const int j = slotAtom[t];
              if (j == i) {
                continue;
              }
              // Full stencil: one direction. Half stencil: the same cell
              // still needs that cut, and a forward cell is seen once.
              if ((PX || PY || PZ || (dx == 0 && dy == 0 && dz == 0)) && j <= i) {
                continue;
              }
              double dxv = ixv - cx[t];
              double dyv = iyv - cy[t];
              double dzv = izv - cz[t];
              if constexpr (PX) {
                if (dxv > 0.5 * ex) {
                  dxv -= ex;
                } else if (dxv < -0.5 * ex) {
                  dxv += ex;
                }
              }
              if constexpr (PY) {
                if (dyv > 0.5 * ey) {
                  dyv -= ey;
                } else if (dyv < -0.5 * ey) {
                  dyv += ey;
                }
              }
              if constexpr (PZ) {
                if (dzv > 0.5 * ez) {
                  dzv -= ez;
                } else if (dzv < -0.5 * ez) {
                  dzv += ez;
                }
              }
              if (dxv * dxv + dyv * dyv + dzv * dzv < cut2) {
                mine.emplace_back(i, j);
              }
            }
          }
        }
      }
    }
}

template <bool PX, bool PY, bool PZ>
void searchWindow(const double *cx, const double *cy, const double *cz, const int *slotAtom,
                  const int *cellStart, const int *cellOf, int n, int nc0, int nc1, int nc2,
                  double cut2, double ex, double ey, double ez, int nThreads,
                  std::vector<PairBag> &partial) {
  // Contiguous slabs, one per thread already pinned by OMP_PLACES.
  seams::RegionErrors errors;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (nThreads > 1)
#endif
  {
    int tid = 0;
    int nt = 1;
#ifdef SEAMS_HAS_OPENMP
    tid = omp_get_thread_num();
    nt = omp_get_num_threads();
#endif
    const int slot0 = static_cast<int>((static_cast<long long>(tid) * n) / nt);
    const int slot1 = static_cast<int>((static_cast<long long>(tid + 1) * n) / nt);
    errors.run([&] {
      searchSlots<PX, PY, PZ>(cx, cy, cz, slotAtom, cellStart, cellOf, slot0, slot1, nc0, nc1,
                              nc2, cut2, ex, ey, ez, partial[static_cast<std::size_t>(tid)].pairs);
    });
  }
  errors.rethrow();
}

} // namespace

bool nneigh::blockWindowRows(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double rcutoff,
    const std::vector<char> *source, std::vector<std::vector<int>> &rows, bool byId,
    bool everyRow) {
  if (yCloud.box.size() != 3 || yCloud.nop <= 0 || rcutoff <= 0.0) {
    return false;
  }
  const double lenBox[3] = {yCloud.box[0], yCloud.box[1], yCloud.box[2]};
  if (!(lenBox[0] > 0.0 && lenBox[1] > 0.0 && lenBox[2] > 0.0)) {
    return false;
  }
  const double origin[3] = {yCloud.boxLow.size() > 0 ? yCloud.boxLow[0] : 0.0,
                            yCloud.boxLow.size() > 1 ? yCloud.boxLow[1] : 0.0,
                            yCloud.boxLow.size() > 2 ? yCloud.boxLow[2] : 0.0};
  const int n = yCloud.nop;
  const auto windowT0 = std::chrono::steady_clock::now();
#ifdef SEAMS_HAS_OPENMP
  const int nThreads = std::max(1, omp_get_max_threads());
#else
  const int nThreads = 1;
#endif
  // One pass over the anchors. Each thread keeps its own bounds; a
  // reduction on the shared min raced under this libgomp.
  struct Bound {
    double lo[3];
    double hi[3];
    int n = 0;
  };
  std::vector<Bound> bounds(static_cast<std::size_t>(nThreads));
  for (auto &b : bounds) {
    b.lo[0] = b.lo[1] = b.lo[2] = 1.0e300;
    b.hi[0] = b.hi[1] = b.hi[2] = -1.0e300;
    b.n = 0;
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (nThreads > 1)
#endif
  {
    int tid = 0;
#ifdef SEAMS_HAS_OPENMP
    tid = omp_get_thread_num();
#endif
    Bound &b = bounds[static_cast<std::size_t>(tid)];
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int i = 0; i < n; i++) {
      if (source != nullptr &&
          (static_cast<int>(source->size()) <= i || (*source)[static_cast<std::size_t>(i)] == 0)) {
        continue;
      }
      const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
      b.lo[0] = std::min(b.lo[0], p.x);
      b.lo[1] = std::min(b.lo[1], p.y);
      b.lo[2] = std::min(b.lo[2], p.z);
      b.hi[0] = std::max(b.hi[0], p.x);
      b.hi[1] = std::max(b.hi[1], p.y);
      b.hi[2] = std::max(b.hi[2], p.z);
      ++b.n;
    }
  }
  double lo0 = 1.0e300, lo1 = 1.0e300, lo2 = 1.0e300;
  double hi0 = -1.0e300, hi1 = -1.0e300, hi2 = -1.0e300;
  int nAnchor = 0;
  for (const auto &b : bounds) {
    lo0 = std::min(lo0, b.lo[0]);
    lo1 = std::min(lo1, b.lo[1]);
    lo2 = std::min(lo2, b.lo[2]);
    hi0 = std::max(hi0, b.hi[0]);
    hi1 = std::max(hi1, b.hi[1]);
    hi2 = std::max(hi2, b.hi[2]);
    nAnchor += b.n;
  }
  if (nAnchor == 0) {
    return false;
  }
  // A domain that crosses the periodic boundary has atoms at both ends of
  // the box. The largest empty gap is the outside of that domain. Folding
  // the anchors onto the other side keeps the window shorter than the box.
  auto foldSpan = [&](double &lo, double &hi, int axis) {
    const double edge = lenBox[axis];
    if (!(edge > 0.0) || !(hi - lo > 0.5 * edge) || nAnchor < 2) {
      return;
    }
    constexpr int kBins = 64;
    int hist[kBins] = {};
    for (int i = 0; i < n; i++) {
      if (source != nullptr &&
          (static_cast<int>(source->size()) <= i || (*source)[static_cast<std::size_t>(i)] == 0)) {
        continue;
      }
      const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
      const double xyz[3] = {p.x, p.y, p.z};
      double f = xyz[axis] / edge;
      f -= std::floor(f);
      int b = static_cast<int>(f * kBins);
      if (b < 0) {
        b = 0;
      }
      if (b >= kBins) {
        b = kBins - 1;
      }
      hist[b]++;
    }
    int best = 0;
    int startBin = 0;
    for (int i = 0; i < kBins; i++) {
      if (hist[i] != 0) {
        continue;
      }
      int run = 1;
      while (run < kBins && hist[(i + run) % kBins] == 0) {
        run++;
      }
      if (run > best) {
        best = run;
        startBin = i;
      }
    }
    if (best <= 0 || best >= kBins) {
      return;
    }
    const double width = edge / kBins;
    const int cluster = (startBin + best) % kBins;
    lo = cluster * width - width;
    hi = startBin * width + width;
    if (hi <= lo) {
      hi += edge;
    }
  };
  if (source != nullptr) {
    foldSpan(lo0, hi0, 0);
    foldSpan(lo1, hi1, 1);
    foldSpan(lo2, hi2, 2);
  }
  const double alo = lo0, ahi0 = hi0, alo1 = lo1, ahi1 = hi1, alo2 = lo2, ahi2 = hi2;
  const double box0 = lenBox[0], box1 = lenBox[1], box2 = lenBox[2];
  // Owned atoms of one block do not wrap. Ghosts that do are shifted by a
  // box length so they sit against that block. The grid is the resulting
  // window, and it does not wrap.
  std::vector<double> local(static_cast<std::size_t>(n) * 3);
  double *const lptr = local.data();
  std::vector<Bound> span(static_cast<std::size_t>(nThreads));
  for (auto &b : span) {
    b.lo[0] = b.lo[1] = b.lo[2] = 1.0e300;
    b.hi[0] = b.hi[1] = b.hi[2] = -1.0e300;
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (nThreads > 1)
#endif
  {
    int tid = 0;
#ifdef SEAMS_HAS_OPENMP
    tid = omp_get_thread_num();
#endif
    Bound &b = span[static_cast<std::size_t>(tid)];
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int i = 0; i < n; i++) {
      const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
      const double c[3] = {p.x, p.y, p.z};
      const double boxes[3] = {box0, box1, box2};
      const double los[3] = {alo, alo1, alo2};
      const double his[3] = {ahi0, ahi1, ahi2};
      for (int d = 0; d < 3; d++) {
        double best = c[d];
        double bestDist = best < los[d] ? los[d] - best : (best > his[d] ? best - his[d] : 0.0);
        for (int wrap = -1; wrap <= 1; wrap += 2) {
          const double s = c[d] + wrap * boxes[d];
          const double dist = s < los[d] ? los[d] - s : (s > his[d] ? s - his[d] : 0.0);
          if (dist < bestDist) {
            bestDist = dist;
            best = s;
          }
        }
        lptr[static_cast<std::size_t>(i) * 3 + static_cast<std::size_t>(d)] = best;
        b.lo[d] = std::min(b.lo[d], best);
        b.hi[d] = std::max(b.hi[d], best);
      }
    }
  }
  double plow0 = 1.0e300, plow1 = 1.0e300, plow2 = 1.0e300;
  double plen0 = -1.0e300, plen1 = -1.0e300, plen2 = -1.0e300;
  for (const auto &b : span) {
    plow0 = std::min(plow0, b.lo[0]);
    plow1 = std::min(plow1, b.lo[1]);
    plow2 = std::min(plow2, b.lo[2]);
    plen0 = std::max(plen0, b.hi[0]);
    plen1 = std::max(plen1, b.hi[1]);
    plen2 = std::max(plen2, b.hi[2]);
  }
  const double plow[3] = {plow0, plow1, plow2};
  const double plen[3] = {plen0, plen1, plen2};
  bool opened = false;
  double len[3];
  for (int d = 0; d < 3; d++) {
    len[d] = std::max(plen[d] - plow[d], 0.0);
    if (lenBox[d] - len[d] > rcutoff) {
      opened = true;
    }
  }
  const double shift0 = plow[0], shift1 = plow[1], shift2 = plow[2];
  const double len0 = len[0], len1 = len[1], len2 = len[2];
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nThreads > 1)
#endif
  for (int i = 0; i < n; i++) {
    const double shifts[3] = {shift0, shift1, shift2};
    const double lens[3] = {len0, len1, len2};
    for (int d = 0; d < 3; d++) {
      double &s = lptr[static_cast<std::size_t>(i) * 3 + static_cast<std::size_t>(d)];
      s -= shifts[d];
      if (!(s < lens[d]) && lens[d] > 0.0) {
        s = std::nextafter(lens[d], 0.0);
      }
    }
  }
  if (!opened) {
    return false;
  }
  (void)origin;
  bool periodic[3];
  double edge[3];
  double originW[3];
  for (int d = 0; d < 3; d++) {
    // A gap wider than the cutoff holds no pair, so that axis does not wrap.
    // One extra cutoff keeps the periodic image of the window outside the shell.
    periodic[d] = !(lenBox[d] - len[d] > rcutoff);
    if (periodic[d]) {
      edge[d] = lenBox[d];
      originW[d] = origin[d];
    } else {
      edge[d] = std::max(len[d], 0.0) + rcutoff + 1.0e-6;
      originW[d] = 0.0;
    }
  }
  const seams::PinSession placed;
  seams::adoptThreads();
  // Not thread_local: the fill below is shared by the team. A thread-local
  // buffer would be empty on every worker and the write would fault.
  std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
  std::vector<char> anchor(static_cast<std::size_t>(n), 1);
  if (source != nullptr && !everyRow) {
    for (int i = 0; i < n; i++) {
      anchor[static_cast<std::size_t>(i)] =
          static_cast<int>(source->size()) > i && (*source)[static_cast<std::size_t>(i)] != 0;
    }
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nThreads > 1)
#endif
  for (int i = 0; i < n; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    const double raw[3] = {p.x, p.y, p.z};
    for (int d = 0; d < 3; d++) {
      xyz[static_cast<std::size_t>(i) * 3 + static_cast<std::size_t>(d)] =
          periodic[d] ? raw[d]
                      : local[static_cast<std::size_t>(i) * 3 + static_cast<std::size_t>(d)];
    }
  }
  // A bin of one cutoff, wrapped on the axes that still see the whole box.
  // The previous periodic walk built that full box in linkcell and did not
  // get faster as the rank count grew, so two or four ranks stayed near the
  // one-rank time. Threads of one rank share this loop.
  const double cut2 = rcutoff * rcutoff;
  int nc[3];
  for (int d = 0; d < 3; d++) {
    nc[d] = std::max(1, static_cast<int>(std::floor(edge[d] / rcutoff)));
  }
  seams::coarsenUntilFits(nc[0], nc[1], nc[2], 1);
  int nCells = 0;
  if (!seams::gridFits(nc[0], nc[1], nc[2], nCells)) {
    return false;
  }
  std::vector<int> cellOf(static_cast<std::size_t>(n), 0);
  std::vector<int> cellCount(static_cast<std::size_t>(nCells), 0);
  auto cellIndex = [&](double s, int axis) {
    int ix = static_cast<int>(std::floor(s / rcutoff));
    if (ix < 0) {
      ix = 0;
    }
    if (ix >= nc[axis]) {
      ix = nc[axis] - 1;
    }
    return ix;
  };
  for (int i = 0; i < n; i++) {
    const double *p = xyz.data() + static_cast<std::size_t>(i) * 3;
    const int c = (cellIndex(p[2], 2) * nc[1] + cellIndex(p[1], 1)) * nc[0] + cellIndex(p[0], 0);
    cellOf[static_cast<std::size_t>(i)] = c;
    cellCount[static_cast<std::size_t>(c)]++;
  }
  std::vector<int> cellStart(static_cast<std::size_t>(nCells) + 1, 0);
  for (int c = 0; c < nCells; c++) {
    cellStart[static_cast<std::size_t>(c) + 1] =
        cellStart[static_cast<std::size_t>(c)] + cellCount[static_cast<std::size_t>(c)];
  }
  std::vector<int> cellAtom(static_cast<std::size_t>(n));
  {
    std::vector<int> fill = cellStart;
    for (int i = 0; i < n; i++) {
      const int c = cellOf[static_cast<std::size_t>(i)];
      cellAtom[static_cast<std::size_t>(fill[static_cast<std::size_t>(c)]++)] = i;
    }
  }
  // Cell order, one coordinate array per axis. The pair walk then reads
  // neighbours of a bin sequentially instead of gathering the frame order.
  std::vector<double> cx(static_cast<std::size_t>(n));
  std::vector<double> cy(static_cast<std::size_t>(n));
  std::vector<double> cz(static_cast<std::size_t>(n));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nThreads > 1)
#endif
  for (int slot = 0; slot < n; slot++) {
    const int i = cellAtom[static_cast<std::size_t>(slot)];
    const double *p = xyz.data() + static_cast<std::size_t>(i) * 3;
    cx[static_cast<std::size_t>(slot)] = p[0];
    cy[static_cast<std::size_t>(slot)] = p[1];
    cz[static_cast<std::size_t>(slot)] = p[2];
  }
  const auto pairT0 = std::chrono::steady_clock::now();
  if (const char *split = std::getenv("SEAMS_SPLIT")) {
    if (split[0] == '1') {
      const double prepMs =
          std::chrono::duration<double, std::milli>(pairT0 - windowT0).count();
      std::fprintf(stderr, "prep ms %.3f\n", prepMs);
    }
  }
  // One cache line per bag. Adjacent vector headers were sharing a
  // line, and every emplace on a hybrid team invalidated the others.
  std::vector<PairBag> partial(static_cast<std::size_t>(nThreads));
  for (auto &bag : partial) {
    bag.pairs.reserve(static_cast<std::size_t>(n) * 4 / static_cast<std::size_t>(nThreads) + 16);
  }
  const int mask = (periodic[0] ? 1 : 0) | (periodic[1] ? 2 : 0) | (periodic[2] ? 4 : 0);
  const auto runSearch = [&](auto &&go) {
    go(cx.data(), cy.data(), cz.data(), cellAtom.data(), cellStart.data(), cellOf.data(), n, nc[0],
       nc[1], nc[2], cut2, edge[0], edge[1], edge[2], nThreads, partial);
  };
  switch (mask) {
  case 0:
    runSearch(searchWindow<false, false, false>);
    break;
  case 1:
    runSearch(searchWindow<true, false, false>);
    break;
  case 2:
    runSearch(searchWindow<false, true, false>);
    break;
  case 3:
    runSearch(searchWindow<true, true, false>);
    break;
  case 4:
    runSearch(searchWindow<false, false, true>);
    break;
  case 5:
    runSearch(searchWindow<true, false, true>);
    break;
  case 6:
    runSearch(searchWindow<false, true, true>);
    break;
  default:
    runSearch(searchWindow<true, true, true>);
    break;
  }
  std::size_t nPairs = 0;
  for (const auto &part : partial) {
    nPairs += part.pairs.size();
  }
  if (const char *split = std::getenv("SEAMS_SPLIT")) {
    if (split[0] == '1') {
      const double pairMs =
          std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - pairT0)
              .count();
      std::fprintf(stderr, "pairs threads %d n %d cells %d hits %zu ms %.3f\n", nThreads, n, nCells,
                   nPairs, pairMs);
    }
  }
  const std::vector<int> indexToID = byId ? indexToIDTable(yCloud) : std::vector<int>{};
  std::vector<int> threadDegree(static_cast<std::size_t>(nThreads) * static_cast<std::size_t>(n), 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (nThreads > 1)
#endif
  {
    int tid = 0;
#ifdef SEAMS_HAS_OPENMP
    tid = omp_get_thread_num();
#endif
    int *deg = threadDegree.data() + static_cast<std::size_t>(tid) * static_cast<std::size_t>(n);
    for (const auto &hit : partial[static_cast<std::size_t>(tid)].pairs) {
      if (anchor[static_cast<std::size_t>(hit.first)]) {
        ++deg[hit.first];
      }
      if (anchor[static_cast<std::size_t>(hit.second)]) {
        ++deg[hit.second];
      }
    }
  }
  std::vector<int> degree(static_cast<std::size_t>(n), 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nThreads > 1)
#endif
  for (int i = 0; i < n; i++) {
    int sum = 0;
    for (int t = 0; t < nThreads; t++) {
      sum += threadDegree[static_cast<std::size_t>(t) * static_cast<std::size_t>(n) +
                          static_cast<std::size_t>(i)];
    }
    degree[static_cast<std::size_t>(i)] = sum;
  }
  std::vector<int> cursor(static_cast<std::size_t>(n) + 1, 0);
  for (int i = 0; i < n; i++) {
    cursor[static_cast<std::size_t>(i) + 1] =
        cursor[static_cast<std::size_t>(i)] + degree[static_cast<std::size_t>(i)];
  }
  std::vector<int> partner(static_cast<std::size_t>(cursor[static_cast<std::size_t>(n)]));
  std::vector<int> threadAt(static_cast<std::size_t>(nThreads) * static_cast<std::size_t>(n), 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nThreads > 1)
#endif
  for (int i = 0; i < n; i++) {
    int at = cursor[static_cast<std::size_t>(i)];
    for (int t = 0; t < nThreads; t++) {
      threadAt[static_cast<std::size_t>(t) * static_cast<std::size_t>(n) +
               static_cast<std::size_t>(i)] = at;
      at += threadDegree[static_cast<std::size_t>(t) * static_cast<std::size_t>(n) +
                         static_cast<std::size_t>(i)];
    }
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (nThreads > 1)
#endif
  {
    int tid = 0;
#ifdef SEAMS_HAS_OPENMP
    tid = omp_get_thread_num();
#endif
    int *at = threadAt.data() + static_cast<std::size_t>(tid) * static_cast<std::size_t>(n);
    for (const auto &hit : partial[static_cast<std::size_t>(tid)].pairs) {
      if (anchor[static_cast<std::size_t>(hit.first)]) {
        partner[static_cast<std::size_t>(at[hit.first]++)] = hit.second;
      }
      if (anchor[static_cast<std::size_t>(hit.second)]) {
        partner[static_cast<std::size_t>(at[hit.second]++)] = hit.first;
      }
    }
  }
  if (static_cast<int>(rows.size()) != n) {
    rows.assign(static_cast<std::size_t>(n), {});
  } else {
    for (auto &row : rows) {
      row.clear();
    }
  }
  auto idOf = [&](int i) {
    if (!byId) {
      return i;
    }
    const int mapped = indexToID[static_cast<std::size_t>(i)];
    if (mapped >= 0) {
      return mapped;
    }
    return yCloud.pts[static_cast<std::size_t>(i)].atomID;
  };
  seams::RegionErrors rowErrors;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nThreads > 1)
#endif
  for (int i = 0; i < n; i++) {
    rowErrors.run([&] {
      const int idi = idOf(i);
      if (idi < 0) {
        return;
      }
      auto &row = rows[static_cast<std::size_t>(i)];
      row.push_back(idi);
      if (!anchor[static_cast<std::size_t>(i)]) {
        return;
      }
      const int begin = cursor[static_cast<std::size_t>(i)];
      const int end = cursor[static_cast<std::size_t>(i) + 1];
      std::sort(partner.begin() + begin, partner.begin() + end);
      int prev = -1;
      for (int t = begin; t < end; t++) {
        const int j = partner[static_cast<std::size_t>(t)];
        if (j == prev) {
          continue;
        }
        prev = j;
        const int idj = idOf(j);
        if (idj >= 0) {
          row.push_back(idj);
        }
      }
    });
  }
  rowErrors.rethrow();
  if (const char *split = std::getenv("SEAMS_SPLIT")) {
    if (split[0] == '1') {
      const double rowMs =
          std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - pairT0)
              .count();
      std::fprintf(stderr, "rows ms %.3f\n", rowMs);
    }
  }
  return true;
}

const std::vector<std::vector<int>> &nneigh::reusedShellList(
    double rcutoff, const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<char> *source) {
  // One scratch set per calling thread. Two overlapping calls used to
  // refill the same rows. The team writes them through these references.
  thread_local std::vector<std::vector<int>> rowScratch;
  thread_local linkcell::PairCsr csrScratch;
  std::vector<std::vector<int>> &rows = rowScratch;
  linkcell::PairCsr &csr = csrScratch;
  const int n = yCloud.nop;
  if (n <= 0 || rcutoff <= 0.0) {
    rows.clear();
    return rows;
  }
  if (blockWindowRows(yCloud, rcutoff, source, rows)) {
    return rows;
  }
  const int width = cutoffWidth();
  linkcellCutoffCsr(yCloud, rcutoff, false, kPeriodic, {}, width, csr);
  const std::vector<int> indexToID = indexToIDTable(yCloud);
  rows.resize(static_cast<std::size_t>(n));
  forEachRow(n, width, [&](int i) {
    auto &row = rows[static_cast<std::size_t>(i)];
    if (source != nullptr && (static_cast<int>(source->size()) <= i ||
                              (*source)[static_cast<std::size_t>(i)] == 0)) {
      row.clear();
      if (indexToID[static_cast<std::size_t>(i)] != -1) {
        row.push_back(indexToID[static_cast<std::size_t>(i)]);
      }
      return;
    }
    idRow(csr, indexToID, i, row);
    if (!row.empty()) {
      std::sort(row.begin() + 1, row.end());
    }
  });
  return rows;
}

void nneigh::preferSerialCutoff(bool on) { gSerialCutoff = on; }

bool nneigh::refillNeighListO(double rcutoff,
                             const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                             int typeI, std::vector<std::vector<int>> &rows) {
  if (!hasPeriodicBox(yCloud)) {
    return false;
  }
  const int width = cutoffWidth();
  thread_local linkcell::PairCsr scratch;
  linkcell::PairCsr &csr = scratch;
  linkcellCutoffCsr(yCloud, rcutoff, false, kPeriodic,
                    typeMask(yCloud, [typeI](int t) { return t == typeI; }), width, csr);
  const std::vector<int> indexToID = indexToIDTable(yCloud);
  rows.resize(static_cast<std::size_t>(yCloud.nop));
  forEachRow(yCloud.nop, width,
             [&](int i) { idRow(csr, indexToID, i, rows[static_cast<std::size_t>(i)]); });
  return true;
}

bool nneigh::refillIndexRows(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                            double cutoff, std::vector<std::vector<int>> &rows) {
  if (!hasPeriodicBox(yCloud)) {
    return false;
  }
  const int width = cutoffWidth();
  thread_local linkcell::PairCsr scratch;
  linkcell::PairCsr &csr = scratch;
  linkcellCutoffCsr(yCloud, cutoff, false, kPeriodic, {}, width, csr);
  rows.resize(static_cast<std::size_t>(yCloud.nop));
  forEachRow(yCloud.nop, width,
             [&](int i) { indexRow(csr, i, rows[static_cast<std::size_t>(i)]); });
  return true;
}

std::vector<std::vector<int>>
nneigh::neighListLinkcell(
    double rcutoff,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeI) {
  if (!hasPeriodicBox(yCloud) || !(rcutoff > 0.0) || yCloud.nop <= 0) {
    return {};
  }
  const seams::PinSession placed;
  seams::prepareCompute();
  int pool = cutoffWidth();
  int wide = seams::poolWidth(pool);
  // Eight threads loses to four below a few thousand molecules: the
  // list and the cage walk do not have a full team of work, and the
  // extra workers show up as a slow sample. Larger crystals keep the
  // width the caller asked for.
  if (yCloud.nop < 4096) {
    pool = std::min(pool, 4);
    wide = std::min(wide, 4);
  }
#ifdef SEAMS_HAS_OPENMP
  if (!omp_in_parallel()) {
    omp_set_num_threads(pool);
  }
#endif
  const auto mask = typeMask(yCloud, [typeI](int t) { return t == typeI; });
  if (!mask.empty() && std::find(mask.begin(), mask.end(), 1) == mask.end()) {
    return {};
  }
  thread_local linkcell::PairCsr scratch;
  linkcell::PairCsr &csr = scratch;
  linkcellCutoffCsr(yCloud, rcutoff, false, kPeriodic, mask, wide, csr);
  return idRows(csr, indexToIDTable(yCloud), yCloud.nop, pool);
}

std::vector<std::vector<int>>
nneigh::neighListO(double rcutoff,
                   const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                   int typeI) {
  // The row fill is an OpenMP walk. Pausing that team here makes the
  // next region pay a wake-up, and an unsized Rayon pool would spin beside
  // it. Record the OpenMP width for linkcell and leave the team running.
  const seams::PinSession placed;
  seams::prepareCompute();
  if (!hasPeriodicBox(yCloud)) {
    return {};
  }
  const int width = cutoffWidth();
  thread_local linkcell::PairCsr scratch;
  linkcell::PairCsr &csr = scratch;
  linkcellCutoffCsr(yCloud, rcutoff, false, kPeriodic,
                    typeMask(yCloud, [typeI](int t) { return t == typeI; }),
                    seams::poolWidth(width), csr);
  return idRows(csr, indexToIDTable(yCloud), yCloud.nop, width);
}

/**
 * @details I-J neighbour list. Like-type (I==J) reuses neighListO,
 *  including its dump MIC. Unlike-type pairs come from the same linkcell
 *  rows over both types, keeping only the mixed I-J contacts.
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
  const seams::PinSession placed;
  seams::prepareCompute();
  if (!hasPeriodicBox(yCloud)) {
    return {};
  }
  const int width = cutoffWidth();
  thread_local linkcell::PairCsr scratch;
  linkcell::PairCsr &csr = scratch;
  linkcellCutoffCsr(yCloud, rcutoff, false, kPeriodic,
                    typeMask(yCloud, [typeI, typeJ](int t) { return t == typeI || t == typeJ; }),
                    seams::poolWidth(width), csr);
  const std::vector<int> indexToID = indexToIDTable(yCloud);
  std::vector<std::vector<int>> nList(static_cast<std::size_t>(yCloud.nop));
  forEachRow(yCloud.nop, width, [&](int i) {
    auto &row = nList[static_cast<std::size_t>(i)];
    idRow(csr, indexToID, i, row);
    if (row.empty()) {
      return;
    }
    const int ti = yCloud.pts[static_cast<std::size_t>(i)].type;
    const int other = ti == typeI ? typeJ : typeI;
    std::size_t w = 1;
    for (int k = csr.offsets[static_cast<std::size_t>(i)];
         k < csr.offsets[static_cast<std::size_t>(i) + 1]; k++) {
      const int j = csr.cols[static_cast<std::size_t>(k)];
      if (indexToID[static_cast<std::size_t>(j)] == -1) {
        continue;
      }
      if (yCloud.pts[static_cast<std::size_t>(j)].type == other) {
        row[w++] = indexToID[static_cast<std::size_t>(j)];
      }
    }
    row.resize(w);
  });
  return nList;
}

/**
 * @details Function for building neighbour lists for each
 *  particle of only one type, from linkcell's half CSR rows. This
 *  generates the half neighbour list, by ID: the lower cloud index owns
 *  the pair. This function will only work for building a neighbour list
 *  between one type of particles.
 * @param[in] rcutoff Distance cutoff, within which two atoms are neighbours.
 * @param[in] yCloud The input molSys::PointCloud
 * @param[in] typeI Type ID of the \f$ i^{th} \f$ particle type.
 * @return Row-ordered half neighbour list, by atom ID.
 */
std::vector<std::vector<int>>
nneigh::halfNeighList(double rcutoff,
                      const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                      int typeI) {
  const seams::PinSession placed;
  seams::prepareCompute();
  if (!hasPeriodicBox(yCloud)) {
    return {};
  }
  const int width = cutoffWidth();
  thread_local linkcell::PairCsr scratch;
  linkcell::PairCsr &csr = scratch;
  linkcellCutoffCsr(yCloud, rcutoff, true, kPeriodic,
                    typeMask(yCloud, [typeI](int t) { return t == typeI; }),
                    seams::poolWidth(width), csr);
  return idRows(csr, indexToIDTable(yCloud), yCloud.nop, width);
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
  const seams::PinSession placed;
  seams::prepareCompute();
  if (!hasPeriodicBox(yCloud)) {
    return {};
  }
  const int width = cutoffWidth();
  thread_local linkcell::PairCsr scratch;
  linkcell::PairCsr &csr = scratch;
  linkcellCutoffCsr(yCloud, cutoff, false, kPeriodic, {}, seams::poolWidth(width), csr);
  std::vector<std::vector<int>> nList(static_cast<std::size_t>(yCloud.nop));
  forEachRow(yCloud.nop, width,
             [&](int i) { indexRow(csr, i, nList[static_cast<std::size_t>(i)]); });
  return nList;
}

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
  seams::RegionErrors errors;
  // nvc++ -mp=gpu SIGSEGVs this host parallel-for in libnvomp
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
#pragma omp parallel for schedule(static) if (nList.size() >= 4096)
#endif
  for (std::int64_t ii = 0; ii < static_cast<std::int64_t>(nList.size()); ii++) {
    errors.run([&] {
      const std::size_t iatom = static_cast<std::size_t>(ii);
      if (nList[iatom].empty()) {
        return;
      }
      const int iatomID = nList[iatom][0];
      auto gotI = yCloud.idIndexMap.find(iatomID);
      if (gotI == yCloud.idIndexMap.end()) {
        return;
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
    });
  }
  errors.rethrow();

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

// k-nearest nominations from the linkcell mesh.
// candidateCutoff is a cell-size hint for knearest_into.
namespace {

bool typeAllowed(int type, int typeI, const std::vector<int> *types) {
  if (types == nullptr || types->empty()) {
    return type == typeI;
  }
  return std::find(types->begin(), types->end(), type) != types->end();
}

// Packed positions and, when some atom is not a source, the mask the
// linkcell search takes. False when there is nothing to search.
bool stageLinkcell(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                   int typeI, const std::vector<int> *types, std::vector<double> &xyz,
                   std::vector<int> &mask, const int *&maskPtr) {
  const int n = yCloud.nop;
  maskPtr = nullptr;
  if (n <= 0 || yCloud.box.size() < 3) {
    return false;
  }
  const double bx = yCloud.box[0];
  const double by = yCloud.box[1];
  const double bz = yCloud.box[2];
  if (!(bx > 0.0 && by > 0.0 && bz > 0.0)) {
    return false;
  }
  if (static_cast<int>(xyz.size()) != n * 3) {
    xyz.resize(static_cast<std::size_t>(n) * 3);
  }
  int nSrc = 0;
  bool uniform = true;
  for (int i = 0; i < n; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(i) * 3 + 0] = p.x;
    xyz[static_cast<std::size_t>(i) * 3 + 1] = p.y;
    xyz[static_cast<std::size_t>(i) * 3 + 2] = p.z;
    if (typeAllowed(p.type, typeI, types)) {
      ++nSrc;
    } else {
      uniform = false;
    }
  }
  if (nSrc == 0) {
    return false;
  }
  // An all-true mask takes the filtered walk. The unmasked walk is the
  // same set when every atom is a source.
  if (!uniform) {
    mask.assign(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; i++) {
      if (typeAllowed(yCloud.pts[static_cast<std::size_t>(i)].type, typeI, types)) {
        mask[static_cast<std::size_t>(i)] = 1;
      }
    }
    maskPtr = mask.data();
  }
  return true;
}

// Writes n*k cloud indices, nearest first, unused slots -1. Throws
// linkcell::Error when the search fails.
void nominateByLinkcell(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    int typeI, double cellHint, std::vector<int> &nom,
    const std::vector<int> *types = nullptr,
    std::vector<double> *keepXyz = nullptr) {
  const int n = yCloud.nop;
  if (n <= 0 || k <= 0) {
    return;
  }
  std::vector<double> localXyz;
  std::vector<double> &xyz = keepXyz != nullptr ? *keepXyz : localXyz;
  std::vector<int> mask;
  const int *maskPtr = nullptr;
  if (!stageLinkcell(yCloud, typeI, types, xyz, mask, maskPtr)) {
    return;
  }
  const linkcell::Cell cell =
      nneigh::lammpsBoxToLinkcell(yCloud.box, yCloud.boxLow);
  seams::withTeamAsleep(seams::teamWidth(), [&] {
    linkcell::knearest_into(xyz.data(), static_cast<std::size_t>(n), cell,
                            static_cast<std::size_t>(k), nom.data(), nom.size(),
                            maskPtr, cellHint);
  });
}

// Packed n*k nominations, nearest first, unused slots -1.
std::vector<int> nominatePacked(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    int typeI, double cellHint, const std::vector<int> *types = nullptr) {
  const int n = yCloud.nop;
  std::vector<int> nom(static_cast<std::size_t>(std::max(n, 0)) *
                           static_cast<std::size_t>(std::max(k, 0)),
                       -1);
  nominateByLinkcell(yCloud, k, typeI, cellHint, nom, types);
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
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
  if (n >= 4096 && k > 0 && omp_get_max_threads() > 1 && !omp_in_parallel()) {
    // The serial fill appends edge (i, j) of nomination slot i * k + t to
    // rows i and j in slot order. A row therefore holds its own accepted
    // slots and the accepted slots naming it, merged by slot. Each thread
    // builds whole rows from those two sorted lists.
    const std::size_t slots = static_cast<std::size_t>(n) * static_cast<std::size_t>(k);
    std::vector<char> accepted(slots, 0);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
      for (int t = 0; t < k; t++) {
        const std::size_t s = static_cast<std::size_t>(i) * static_cast<std::size_t>(k) +
                              static_cast<std::size_t>(t);
        const int j = nom[s];
        if (j < 0 || j >= n || j == i) {
          continue;
        }
        const bool listed = listsPacked(nom, n, k, j, i);
        accepted[s] = mutual ? (i < j && listed) : (i < j || !listed);
      }
    }
    std::vector<int> inOff(static_cast<std::size_t>(n) + 1, 0);
    for (std::size_t s = 0; s < slots; s++) {
      if (accepted[s] != 0) {
        ++inOff[static_cast<std::size_t>(nom[s]) + 1];
      }
    }
    for (int i = 0; i < n; i++) {
      inOff[static_cast<std::size_t>(i) + 1] += inOff[static_cast<std::size_t>(i)];
    }
    std::vector<int> inSlot(static_cast<std::size_t>(inOff[static_cast<std::size_t>(n)]));
    {
      std::vector<int> cursor(inOff.begin(), inOff.end() - 1);
      for (std::size_t s = 0; s < slots; s++) {
        if (accepted[s] != 0) {
          inSlot[static_cast<std::size_t>(cursor[static_cast<std::size_t>(nom[s])]++)] =
              static_cast<int>(s);
        }
      }
    }
    out.assign(static_cast<std::size_t>(n), {});
    seams::RegionErrors errors;
#pragma omp parallel for schedule(static)
    for (int r = 0; r < n; r++) {
      errors.run([&] {
        const std::size_t own = static_cast<std::size_t>(r) * static_cast<std::size_t>(k);
        int nOwn = 0;
        for (int t = 0; t < k; t++) {
          nOwn += accepted[own + static_cast<std::size_t>(t)] != 0 ? 1 : 0;
        }
        const int b = inOff[static_cast<std::size_t>(r)];
        const int e = inOff[static_cast<std::size_t>(r) + 1];
        auto &row = out[static_cast<std::size_t>(r)];
        row.reserve(static_cast<std::size_t>(1 + nOwn + (e - b)));
        row.push_back(yCloud.pts[static_cast<std::size_t>(r)].atomID);
        int t = 0;
        int c = b;
        while (t < k || c < e) {
          const std::size_t mine = own + static_cast<std::size_t>(t);
          if (t < k && accepted[mine] == 0) {
            ++t;
            continue;
          }
          if (t < k &&
              (c >= e || mine < static_cast<std::size_t>(inSlot[static_cast<std::size_t>(c)]))) {
            row.push_back(yCloud.pts[static_cast<std::size_t>(nom[mine])].atomID);
            ++t;
          } else {
            const int from = inSlot[static_cast<std::size_t>(c)] / k;
            row.push_back(yCloud.pts[static_cast<std::size_t>(from)].atomID);
            ++c;
          }
        }
      });
    }
    errors.rethrow();
    return;
  }
#endif
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
  const int n = yCloud.nop;
  // The two graphs read the same nominations. Building them together asks
  // listsPacked once per candidate. Each thread keeps the candidates of its
  // contiguous atom block, and the blocks are replayed in atom order, so
  // the neighbour rows match the serial fills.
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
  const bool scale = n >= 4096 && k > 0 && omp_get_max_threads() > 1 && !omp_in_parallel();
#else
  const bool scale = false;
#endif
  if (!scale) {
    fillGraphFromNom(yCloud, nom, k, true, mutual);
    fillGraphFromNom(yCloud, nom, k, false, uni);
    return;
  }
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
  struct Edge {
    int i;
    int j;
  };
  const int team = omp_get_max_threads();
  std::vector<std::vector<Edge>> mutualParts(static_cast<std::size_t>(team));
  std::vector<std::vector<Edge>> uniParts(static_cast<std::size_t>(team));
  seams::RegionErrors errors;
#pragma omp parallel
  {
    std::vector<Edge> mutualLocal;
    std::vector<Edge> uniLocal;
#pragma omp for schedule(static)
    for (int i = 0; i < n; i++) {
      errors.run([&] {
        for (int t = 0; t < k; t++) {
          const int j =
              nom[static_cast<std::size_t>(i) * static_cast<std::size_t>(k) +
                  static_cast<std::size_t>(t)];
          if (j < 0 || j >= n || j == i) {
            continue;
          }
          const bool listed = listsPacked(nom, n, k, j, i);
          if (i < j && listed) {
            mutualLocal.push_back(Edge{i, j});
          }
          if (i < j || !listed) {
            uniLocal.push_back(Edge{i, j});
          }
        }
      });
    }
    mutualParts[static_cast<std::size_t>(omp_get_thread_num())] = std::move(mutualLocal);
    uniParts[static_cast<std::size_t>(omp_get_thread_num())] = std::move(uniLocal);
  }
  errors.rethrow();
  auto replay = [&](const std::vector<std::vector<Edge>> &parts,
                    std::vector<std::vector<int>> &out) {
    out.assign(static_cast<std::size_t>(n), {});
    for (int i = 0; i < n; i++) {
      auto &row = out[static_cast<std::size_t>(i)];
      row.reserve(8);
      row.push_back(yCloud.pts[static_cast<std::size_t>(i)].atomID);
    }
    for (const auto &part : parts) {
      for (const Edge &edge : part) {
        addUndirectedIDs(out, yCloud, edge.i, edge.j);
      }
    }
  };
  replay(mutualParts, mutual);
  replay(uniParts, uni);
#endif
}

} // namespace

bool nneigh::directedKNearestCsr(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, std::vector<int> &offsets, std::vector<int> &cols) {
  const int n = yCloud.nop;
  offsets.assign(static_cast<std::size_t>(std::max(n, 0)) + 1, 0);
  cols.clear();
  if (n <= 0 || k <= 0 || yCloud.box.size() < 3) {
    return false;
  }
  const double hint = candidateCutoff > 0.0 ? 0.75 * candidateCutoff : 3.0;
  const auto nom = nominatePacked(yCloud, k, 1, hint);
  std::vector<int> counts(static_cast<std::size_t>(n), 0);
  for (int i = 0; i < n; i++) {
    for (int t = 0; t < k; t++) {
      const int j = nom[static_cast<std::size_t>(i) * static_cast<std::size_t>(k) +
                        static_cast<std::size_t>(t)];
      if (j >= 0 && j < n && j != i) {
        counts[static_cast<std::size_t>(i)]++;
      }
    }
  }
  if (!seams::prefixCounts(counts.data(), n, offsets)) {
    return false;
  }
  cols.assign(static_cast<std::size_t>(offsets[static_cast<std::size_t>(n)]), 0);
  std::vector<int> slot(offsets.begin(), offsets.end() - 1);
  for (int i = 0; i < n; i++) {
    for (int t = 0; t < k; t++) {
      const int j = nom[static_cast<std::size_t>(i) * static_cast<std::size_t>(k) +
                        static_cast<std::size_t>(t)];
      if (j >= 0 && j < n && j != i) {
        cols[static_cast<std::size_t>(slot[static_cast<std::size_t>(i)]++)] = j;
      }
    }
  }
  return true;
}

std::vector<std::vector<int>>
nneigh::kNearestNeighbourList(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, int typeI, bool mutual) {
  const seams::PinSession placed;
  seams::prepareCompute();
  [[maybe_unused]] const seams::TeamWidth team(yCloud.nop);
  lc_set_threads(seams::poolFor(yCloud.nop));
  if (yCloud.nop <= 0 || k <= 0) {
    return {};
  }
  const double hint = candidateCutoff > 0.0 ? 0.75 * candidateCutoff : 3.0;
  const bool timing = std::getenv("SEAMS_KNN_TIME") != nullptr;
  const auto t0 = std::chrono::steady_clock::now();
  const auto nom = nominatePacked(yCloud, k, typeI, hint);
  const auto tNom = std::chrono::steady_clock::now();
  std::vector<std::vector<int>> out;
  fillGraphFromNom(yCloud, nom, k, mutual, out);
  if (timing) {
    const auto ms = [](auto a, auto b) {
      return std::chrono::duration<double, std::milli>(b - a).count();
    };
    std::cerr << "knn-split n=" << yCloud.nop << " nom_ms=" << ms(t0, tNom)
              << " graph_ms=" << ms(tNom, std::chrono::steady_clock::now()) << "\n";
  }
  return out;
}

void nneigh::kNearestCsr(const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
                         double candidateCutoff, int typeI, bool mutual,
                         std::vector<int> &offsets, std::vector<int> &cols, bool park) {
  seams::prepareCompute();
  [[maybe_unused]] const seams::TeamWidth team(yCloud.nop);
  // The search and the graph build share one linkcell region. An OpenMP
  // team left in active wait by the previous region sits on those cores,
  // so park it for this call when the caller asks. A caller
  // that classifies in OpenMP immediately afterwards leaves the team
  // running: parking it makes that region pay the wake inside the same
  // call. A search already inside a parallel region stays on the caller.
#ifdef SEAMS_HAS_OPENMP
  if (omp_in_parallel()) {
    lc_set_threads(1);
  } else {
    if (park) {
      omp_pause_resource_all(omp_pause_soft);
    }
    lc_set_threads(seams::poolFor(yCloud.nop));
  }
#else
  (void)park;
  lc_set_threads(seams::poolFor(yCloud.nop));
#endif
  const int n = yCloud.nop;
  offsets.assign(static_cast<std::size_t>(std::max(n, 0)) + 1, 0);
  cols.clear();
  if (n <= 0 || k <= 0) {
    return;
  }
  const double hint = candidateCutoff > 0.0 ? 0.75 * candidateCutoff : 3.0;
  // Linkcell builds the graph from its nominations in the search's own
  // region: the rows are those fillGraphFromNom gives (mutual: lower
  // partners ascending, then higher ones nearest first; otherwise each
  // edge in the order of the nomination slot that adds it). No OpenMP
  // team is started for the call.
  std::vector<double> xyz;
  std::vector<int> mask;
  const int *maskPtr = nullptr;
  if (!stageLinkcell(yCloud, typeI, nullptr, xyz, mask, maskPtr)) {
    return;
  }
  linkcell::knearest_csr(xyz.data(), static_cast<std::size_t>(n),
                         nneigh::lammpsBoxToLinkcell(yCloud.box, yCloud.boxLow),
                         static_cast<std::size_t>(k), mutual, offsets, cols, maskPtr, hint);
}

void nneigh::kNearestPacked(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, int typeI, std::vector<double> &xyz,
    std::vector<int> &nom) {
  // The q path calls this and then an OpenMP harmonic. Pausing the team
  // here makes the next parallel region pay a wake-up, which is the
  // several-millisecond swing on a few thousand atoms.
  const seams::PinSession placed;
  seams::prepareCompute();
  lc_set_threads(seams::poolFor(yCloud.nop));
  const int n = yCloud.nop;
  nom.assign(static_cast<std::size_t>(std::max(n, 0)) * static_cast<std::size_t>(std::max(k, 0)),
             -1);
  if (n <= 0 || k <= 0) {
    return;
  }
  const double hint = candidateCutoff > 0.0 ? 0.75 * candidateCutoff : 3.0;
  nominateByLinkcell(yCloud, k, typeI, hint, nom, nullptr, &xyz);
  if (static_cast<int>(xyz.size()) == n * 3) {
    return;
  }
  xyz.resize(static_cast<std::size_t>(n) * 3);
  for (int i = 0; i < n; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(3 * i)] = p.x;
    xyz[static_cast<std::size_t>(3 * i + 1)] = p.y;
    xyz[static_cast<std::size_t>(3 * i + 2)] = p.z;
  }
}

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
nneigh::kNearestNeighbourPair(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, int typeI) {
  if (yCloud.nop <= 0 || k <= 0) {
    return {{}, {}};
  }
  const seams::PinSession placed;
  seams::prepareCompute();
  lc_set_threads(seams::poolFor(yCloud.nop));
  const double hint = candidateCutoff > 0.0 ? 0.75 * candidateCutoff : 3.0;
  const auto nom = nominatePacked(yCloud, k, typeI, hint);
  std::vector<std::vector<int>> mutual;
  std::vector<std::vector<int>> uni;
  fillBothGraphsFromNom(yCloud, nom, k, mutual, uni);
  return {std::move(mutual), std::move(uni)};
}

std::vector<std::vector<int>>
nneigh::kNearestNeighbourList(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, const std::vector<int> &types, bool mutual) {
  const seams::PinSession placed;
  seams::prepareCompute();
  [[maybe_unused]] const seams::TeamWidth team(yCloud.nop);
  lc_set_threads(seams::poolFor(yCloud.nop));
  if (yCloud.nop <= 0 || k <= 0 || types.empty()) {
    return {};
  }
  const double hint = candidateCutoff > 0.0 ? 0.75 * candidateCutoff : 3.0;
  const auto nom = nominatePacked(yCloud, k, types.front(), hint, &types);
  std::vector<std::vector<int>> out;
  fillGraphFromNom(yCloud, nom, k, mutual, out);
  return out;
}

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
nneigh::kNearestNeighbourPair(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, const std::vector<int> &types) {
  if (yCloud.nop <= 0 || k <= 0 || types.empty()) {
    return {{}, {}};
  }
  const seams::PinSession placed;
  seams::prepareCompute();
  lc_set_threads(seams::poolFor(yCloud.nop));
  const double hint = candidateCutoff > 0.0 ? 0.75 * candidateCutoff : 3.0;
  const auto nom = nominatePacked(yCloud, k, types.front(), hint, &types);
  std::vector<std::vector<int>> mutual;
  std::vector<std::vector<int>> uni;
  fillBothGraphsFromNom(yCloud, nom, k, mutual, uni);
  return {std::move(mutual), std::move(uni)};
}

/**
 * @details Nearest unlike neighbour of every typeI index among typeJ
 *  indices. Each typeI walks every typeJ with the fractional wrap, and
 *  calls the Euclidean image only when that nearest wrap reaches half
 *  the shortest edge. Ties keep the lowest typeJ index. Particles with
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

  const gen::FracBox frame = gen::makeFracBox(yCloud);
  std::vector<int> iIdx;
  iIdx.reserve(static_cast<std::size_t>(yCloud.nop));
  for (int i = 0; i < yCloud.nop; i++) {
    if (yCloud.pts[static_cast<std::size_t>(i)].type == typeI) {
      iIdx.push_back(i);
    }
  }
  out.resize(iIdx.size());
  const int nI = static_cast<int>(iIdx.size());
  const std::size_t nJ = jIdx.size();
#ifdef SEAMS_HAS_MINIMAGE
  // minimage's Euclidean images start the exact search past the Smith ball,
  // which replaces any that are not lattice images. Its engine wrap is not
  // used: minimage 0.1.4 gets it wrong on needle cells even inside the ball.
  const gen::CellFrame cells = gen::makeCellFrame(yCloud);
#endif
  seams::RegionErrors errors;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (nI >= 256)
#endif
  {
    std::vector<double> row;
    errors.run([&] { row.resize(nJ); });
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int t = 0; t < nI; t++) {
      errors.run([&] {
        const int i = iIdx[static_cast<std::size_t>(t)];
        const auto &pi = yCloud.pts[static_cast<std::size_t>(i)];
        for (std::size_t k = 0; k < nJ; k++) {
          const auto &pj = yCloud.pts[static_cast<std::size_t>(jIdx[k])];
          row[k] = frame.ok ? gen::fracDistSq(frame, pi.x, pi.y, pi.z, pj.x, pj.y, pj.z)
                            : gen::periodicDistSq(frame, yCloud, i, jIdx[k]);
        }
        int bestJ = -1;
        double bestD2 = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0; k < nJ; k++) {
          if (jIdx[k] != i && row[k] < bestD2) {
            bestD2 = row[k];
            bestJ = jIdx[k];
          }
        }
        if (bestJ >= 0 && frame.ok && !gen::smithInside(frame, bestD2)) {
          bestJ = -1;
          bestD2 = std::numeric_limits<double>::infinity();
          for (std::size_t k = 0; k < nJ; k++) {
            const int j = jIdx[k];
            if (j == i) {
              continue;
            }
            double d2 = row[k];
            if (!gen::smithInside(frame, d2)) {
              d2 = -1.0;
#ifdef SEAMS_HAS_MINIMAGE
              if (cells.ok) {
                d2 = gen::euclideanDistSq(cells, i, j);
              }
#endif
              if (d2 < 0.0) {
                d2 = gen::periodicDistSq(frame, yCloud, i, j);
              }
            }
            if (d2 < bestD2) {
              bestD2 = d2;
              bestJ = j;
            }
          }
        }
        if (bestJ >= 0) {
          out[static_cast<std::size_t>(t)] = {i, bestJ, std::sqrt(bestD2)};
        } else {
          out[static_cast<std::size_t>(t)] = {-1, -1, 0.0};
        }
      });
    }
  }
  errors.rethrow();
  out.erase(std::remove_if(out.begin(), out.end(),
                           [](const std::tuple<int, int, double> &row) {
                             return std::get<0>(row) < 0;
                           }),
            out.end());
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
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  std::vector<double> dists;
  for (int i = 0; i < yCloud.nop; i++) {
    if (yCloud.pts[static_cast<std::size_t>(i)].type != typeI) {
      continue;
    }
    const auto &pi = yCloud.pts[static_cast<std::size_t>(i)];
    dists.clear();
    for (int j = 0; j < yCloud.nop; j++) {
      if (j == i || yCloud.pts[static_cast<std::size_t>(j)].type != typeI) {
        continue;
      }
      const auto &pj = yCloud.pts[static_cast<std::size_t>(j)];
      dists.push_back(frame.ok ? gen::fracDistSq(frame, pi.x, pi.y, pi.z, pj.x,
                                                 pj.y, pj.z)
                               : gen::periodicDistSq(frame, yCloud, i, j));
    }
    // A site with fewer than k others has no k-th neighbour, so no cutoff
    // reproduces its k-nearest row. The certificate fails rather than
    // skipping that site.
    const int nOther = static_cast<int>(dists.size());
    if (nOther < k) {
      return {std::numeric_limits<double>::infinity(), 0.0};
    }
    // With exactly k others the k-th distance still bounds the cutoff from
    // below; there is no (k+1)-th distance to bound it from above.
    const int want = nOther < k + 1 ? k : k + 1;
    std::partial_sort(dists.begin(), dists.begin() + want, dists.end());
    if (frame.ok && !gen::smithInside(frame, dists[static_cast<std::size_t>(want - 1)])) {
      dists.clear();
      for (int j = 0; j < yCloud.nop; j++) {
        if (j == i || yCloud.pts[static_cast<std::size_t>(j)].type != typeI) {
          continue;
        }
        dists.push_back(gen::periodicDistSq(frame, yCloud, i, j));
      }
      std::partial_sort(dists.begin(), dists.begin() + want, dists.end());
    }
    if (want == k) {
      maxKth = std::max(maxKth, std::sqrt(dists[static_cast<std::size_t>(k - 1)]));
      continue;
    }
    maxKth = std::max(maxKth, std::sqrt(dists[static_cast<std::size_t>(k - 1)]));
    const double next = std::sqrt(dists[static_cast<std::size_t>(k)]);
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
  // Same trigger as an atom: a box edge or tilt that has changed by at
  // most skin/2 cannot bring a pair from outside the candidate shell
  // inside the cutoff.
  const double trigger = 0.5 * skin_;
  for (std::size_t k = 0; k < yCloud.box.size(); k++) {
    if (std::fabs(yCloud.box[k] - box0_[k]) > trigger) {
      return true;
    }
  }
  std::vector<double> old(3);
  for (int i = 0; i < yCloud.nop; i++) {
    if (yCloud.pts[static_cast<std::size_t>(i)].type != typeI_) {
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
  if (!hasPeriodicBox(yCloud)) {
    throw std::invalid_argument("SkinNeighborList needs the cloud's box");
  }
  thread_local linkcell::PairCsr scratch;
  linkcell::PairCsr &csr = scratch;
  linkcellCutoffCsr(yCloud, wide, true, kPeriodic,
                    typeMask(yCloud, [this](int t) { return t == typeI_; }),
                    seams::poolWidth(cutoffWidth()), csr);
  candidates_.reserve(csr.cols.size());
  for (int i = 0; i < yCloud.nop; i++) {
    for (int k = csr.offsets[static_cast<std::size_t>(i)];
         k < csr.offsets[static_cast<std::size_t>(i) + 1]; k++) {
      candidates_.emplace_back(i, csr.cols[static_cast<std::size_t>(k)]);
    }
  }
  x0_.resize(static_cast<std::size_t>(yCloud.nop));
  y0_.resize(static_cast<std::size_t>(yCloud.nop));
  z0_.resize(static_cast<std::size_t>(yCloud.nop));
  for (int i = 0; i < yCloud.nop; i++) {
    x0_[static_cast<std::size_t>(i)] = yCloud.pts[static_cast<std::size_t>(i)].x;
    y0_[static_cast<std::size_t>(i)] = yCloud.pts[static_cast<std::size_t>(i)].y;
    z0_[static_cast<std::size_t>(i)] = yCloud.pts[static_cast<std::size_t>(i)].z;
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
    const gen::FracBox frame = gen::makeFracBox(yCloud);
    for (std::size_t t = 0; t < m; t++) {
      const int i = candidates_[t].first;
      const int j = candidates_[t].second;
      if (i < 0 || j < 0 || i >= yCloud.nop || j >= yCloud.nop) {
        r2[t] = 0.0;
        continue;
      }
      r2[t] = gen::periodicDistSq(frame, yCloud, i, j);
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
      if (r2[t] < cutoffSq_) {
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
