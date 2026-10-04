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

#include <neighbours.hpp>
#include <rdf2d.hpp>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

#ifdef SEAMS_HAS_VESIN
#include <vesin.h>
#endif

#ifdef SEAMS_HAS_HWY
#include "hwy/highway.h"
#endif

namespace {

#ifdef SEAMS_HAS_HWY
namespace hn = hwy::HWY_NAMESPACE;

HWY_ATTR void accumulatePackedCell(const double *px, const double *py,
                                   const double *pz, int begin, int end,
                                   double sxi, double syi, double szi,
                                   double lx, double ly, double lz, double xy,
                                   double xz, double yz, double cut2,
                                   double cutoff, double binwidth, int nbin,
                                   int *local) {
  const hn::ScalableTag<double> d;
  const size_t lanes = hn::Lanes(d);
  const auto vsxi = hn::Set(d, sxi);
  const auto vsyi = hn::Set(d, syi);
  const auto vszi = hn::Set(d, szi);
  const auto vlx = hn::Set(d, lx);
  const auto vly = hn::Set(d, ly);
  const auto vlz = hn::Set(d, lz);
  const auto vxy = hn::Set(d, xy);
  const auto vxz = hn::Set(d, xz);
  const auto vyz = hn::Set(d, yz);
  alignas(64) double r2buf[16];
  int t = begin;
  if (lanes <= 16) {
  for (; t + static_cast<int>(lanes) <= end; t += static_cast<int>(lanes)) {
    const auto dsx0 =
        hn::Sub(vsxi, hn::LoadU(d, px + static_cast<std::size_t>(t)));
    const auto dsy0 =
        hn::Sub(vsyi, hn::LoadU(d, py + static_cast<std::size_t>(t)));
    const auto dsz0 =
        hn::Sub(vszi, hn::LoadU(d, pz + static_cast<std::size_t>(t)));
    const auto dsx = hn::NegMulAdd(hn::Set(d, 1.0), hn::Round(dsx0), dsx0);
    const auto dsy = hn::NegMulAdd(hn::Set(d, 1.0), hn::Round(dsy0), dsy0);
    const auto dsz = hn::NegMulAdd(hn::Set(d, 1.0), hn::Round(dsz0), dsz0);
    const auto rx =
        hn::MulAdd(vxz, dsz, hn::MulAdd(vxy, dsy, hn::Mul(vlx, dsx)));
    const auto ry = hn::MulAdd(vyz, dsz, hn::Mul(vly, dsy));
    const auto rz = hn::Mul(vlz, dsz);
    const auto r2 =
        hn::MulAdd(rx, rx, hn::MulAdd(ry, ry, hn::Mul(rz, rz)));
    hn::StoreU(r2, d, r2buf);
    for (size_t k = 0; k < lanes; k++) {
      if (r2buf[k] > cut2) {
        continue;
      }
      const double r = std::sqrt(r2buf[k]);
      if (r > cutoff) {
        continue;
      }
      int b = static_cast<int>(r / binwidth);
      if (b < 0) {
        continue;
      }
      if (b >= nbin) {
        b = nbin - 1;
      }
      local[b] += 2;
    }
  }
  }
  for (; t < end; t++) {
    double dsx = sxi - px[static_cast<std::size_t>(t)];
    double dsy = syi - py[static_cast<std::size_t>(t)];
    double dsz = szi - pz[static_cast<std::size_t>(t)];
    dsx -= std::round(dsx);
    dsy -= std::round(dsy);
    dsz -= std::round(dsz);
    const double rx = lx * dsx + xy * dsy + xz * dsz;
    const double ry = ly * dsy + yz * dsz;
    const double rz = lz * dsz;
    const double r2 = rx * rx + ry * ry + rz * rz;
    if (r2 > cut2) {
      continue;
    }
    const double r = std::sqrt(r2);
    if (r > cutoff) {
      continue;
    }
    int b = static_cast<int>(r / binwidth);
    if (b < 0) {
      continue;
    }
    if (b >= nbin) {
      b = nbin - 1;
    }
    local[b] += 2;
  }
}
#endif

}  // namespace

namespace {

// Cutoff-sized periodic grid. Below half the shortest edge each pair has one
// image, so a neighbour-cell hit is binned directly.
bool histogramPackedGrid(
    const gen::FracBox &frame,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double cutoff, double binwidth, int nbin, std::vector<int> &histogram) {
  const int n = yCloud.nop;
  if (!frame.ok || n < 2 || !(cutoff > 0.0) || nbin <= 0) {
    return false;
  }
  const double span = cutoff + std::max(cutoff, 1.0) * 1e-12;
  if (!(frame.lx > span && frame.ly > span && frame.lz > span)) {
    return false;
  }
  const int nx = std::max(1, static_cast<int>(std::floor(frame.lx / span)));
  const int ny = std::max(1, static_cast<int>(std::floor(frame.ly / span)));
  const int nz = std::max(1, static_cast<int>(std::floor(frame.lz / span)));
  const long long ncell64 =
      static_cast<long long>(nx) * ny * static_cast<long long>(nz);
  if (ncell64 < 27 || ncell64 > static_cast<long long>(n) * 8 ||
      ncell64 > 4000000LL) {
    return false;
  }
  const int ncell = static_cast<int>(ncell64);

  std::vector<double> sx(static_cast<std::size_t>(n));
  std::vector<double> sy(static_cast<std::size_t>(n));
  std::vector<double> sz(static_cast<std::size_t>(n));
  std::vector<int> cellOf(static_cast<std::size_t>(n));
  std::vector<int> count(static_cast<std::size_t>(ncell), 0);
  const double invLx = 1.0 / frame.lx;
  const double invLy = 1.0 / frame.ly;
  const double invLz = 1.0 / frame.lz;
  for (int i = 0; i < n; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    const double fz = (p.z - frame.oz) * invLz;
    const double fy = (p.y - frame.oy - frame.yz * fz) * invLy;
    const double fx = (p.x - frame.ox - frame.xy * fy - frame.xz * fz) * invLx;
    double wx = fx - std::floor(fx);
    double wy = fy - std::floor(fy);
    double wz = fz - std::floor(fz);
    if (wx >= 1.0) {
      wx = 0.0;
    }
    if (wy >= 1.0) {
      wy = 0.0;
    }
    if (wz >= 1.0) {
      wz = 0.0;
    }
    int ix = static_cast<int>(wx * nx);
    int iy = static_cast<int>(wy * ny);
    int iz = static_cast<int>(wz * nz);
    if (ix >= nx) {
      ix = nx - 1;
    }
    if (iy >= ny) {
      iy = ny - 1;
    }
    if (iz >= nz) {
      iz = nz - 1;
    }
    if (ix < 0) {
      ix = 0;
    }
    if (iy < 0) {
      iy = 0;
    }
    if (iz < 0) {
      iz = 0;
    }
    sx[static_cast<std::size_t>(i)] = wx;
    sy[static_cast<std::size_t>(i)] = wy;
    sz[static_cast<std::size_t>(i)] = wz;
    const int cell = (ix * ny + iy) * nz + iz;
    cellOf[static_cast<std::size_t>(i)] = cell;
    count[static_cast<std::size_t>(cell)] += 1;
  }

  std::vector<int> offsets(static_cast<std::size_t>(ncell) + 1, 0);
  for (int c = 0; c < ncell; c++) {
    offsets[static_cast<std::size_t>(c) + 1] =
        offsets[static_cast<std::size_t>(c)] + count[static_cast<std::size_t>(c)];
  }
  std::vector<int> cursor = offsets;
  std::vector<int> ids(static_cast<std::size_t>(n));
  std::vector<double> px(static_cast<std::size_t>(n));
  std::vector<double> py(static_cast<std::size_t>(n));
  std::vector<double> pz(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    const int slot = cursor[static_cast<std::size_t>(cellOf[static_cast<std::size_t>(i)])]++;
    ids[static_cast<std::size_t>(slot)] = i;
    px[static_cast<std::size_t>(slot)] = sx[static_cast<std::size_t>(i)];
    py[static_cast<std::size_t>(slot)] = sy[static_cast<std::size_t>(i)];
    pz[static_cast<std::size_t>(slot)] = sz[static_cast<std::size_t>(i)];
  }

  std::vector<int> neighOf(static_cast<std::size_t>(ncell) + 1, 0);
  std::vector<int> neigh;
  neigh.reserve(static_cast<std::size_t>(ncell) * 27);
  std::vector<int> stamp(static_cast<std::size_t>(ncell), 0);
  int stampId = 0;
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        const int home = (ix * ny + iy) * nz + iz;
        neighOf[static_cast<std::size_t>(home)] = static_cast<int>(neigh.size());
        ++stampId;
        for (int dx = -1; dx <= 1; dx++) {
          int jx = ix + dx;
          jx %= nx;
          if (jx < 0) {
            jx += nx;
          }
          for (int dy = -1; dy <= 1; dy++) {
            int jy = iy + dy;
            jy %= ny;
            if (jy < 0) {
              jy += ny;
            }
            for (int dz = -1; dz <= 1; dz++) {
              int jz = iz + dz;
              jz %= nz;
              if (jz < 0) {
                jz += nz;
              }
              const int nb = (jx * ny + jy) * nz + jz;
              if (stamp[static_cast<std::size_t>(nb)] == stampId) {
                continue;
              }
              stamp[static_cast<std::size_t>(nb)] = stampId;
              neigh.push_back(nb);
            }
          }
        }
      }
    }
  }
  neighOf[static_cast<std::size_t>(ncell)] = static_cast<int>(neigh.size());

  const double cut2 = cutoff * cutoff;
  const double lx = frame.lx;
  const double ly = frame.ly;
  const double lz = frame.lz;
  const double xy = frame.xy;
  const double xz = frame.xz;
  const double yz = frame.yz;
  int nthreads = 1;
#ifdef SEAMS_HAS_OPENMP
  nthreads = omp_get_max_threads();
#endif
  std::vector<int> locals(static_cast<std::size_t>(nthreads) *
                              static_cast<std::size_t>(nbin),
                          0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (n >= 512)
#endif
  {
    int tid = 0;
#ifdef SEAMS_HAS_OPENMP
    tid = omp_get_thread_num();
#endif
    int *local = locals.data() + static_cast<std::size_t>(tid) *
                                     static_cast<std::size_t>(nbin);
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int i = 0; i < n; i++) {
      const double sxi = sx[static_cast<std::size_t>(i)];
      const double syi = sy[static_cast<std::size_t>(i)];
      const double szi = sz[static_cast<std::size_t>(i)];
      const int home = cellOf[static_cast<std::size_t>(i)];
      const int nb0 = neighOf[static_cast<std::size_t>(home)];
      const int nb1 = neighOf[static_cast<std::size_t>(home) + 1];
      for (int u = nb0; u < nb1; u++) {
        const int nb = neigh[static_cast<std::size_t>(u)];
        const int begin = offsets[static_cast<std::size_t>(nb)];
        const int end = offsets[static_cast<std::size_t>(nb) + 1];
        const int *idbase = ids.data() + begin;
        const int len = end - begin;
        const int *upper = std::upper_bound(idbase, idbase + len, i);
        const int start = begin + static_cast<int>(upper - idbase);
#ifdef SEAMS_HAS_HWY
        accumulatePackedCell(px.data(), py.data(), pz.data(), start, end, sxi,
                             syi, szi, lx, ly, lz, xy, xz, yz, cut2, cutoff,
                             binwidth, nbin, local);
#else
        for (int t = start; t < end; t++) {
          double dsx = sxi - px[static_cast<std::size_t>(t)];
          double dsy = syi - py[static_cast<std::size_t>(t)];
          double dsz = szi - pz[static_cast<std::size_t>(t)];
          dsx -= std::round(dsx);
          dsy -= std::round(dsy);
          dsz -= std::round(dsz);
          const double rx = lx * dsx + xy * dsy + xz * dsz;
          const double ry = ly * dsy + yz * dsz;
          const double rz = lz * dsz;
          const double r2 = rx * rx + ry * ry + rz * rz;
          if (r2 > cut2) {
            continue;
          }
          const double r = std::sqrt(r2);
          if (r > cutoff) {
            continue;
          }
          int b = static_cast<int>(r / binwidth);
          if (b < 0) {
            continue;
          }
          if (b >= nbin) {
            b = nbin - 1;
          }
          local[b] += 2;
        }
#endif
      }
    }
  }
  for (int tid = 0; tid < nthreads; tid++) {
    const int *local = locals.data() + static_cast<std::size_t>(tid) *
                                           static_cast<std::size_t>(nbin);
    for (int b = 0; b < nbin; b++) {
      histogram[static_cast<std::size_t>(b)] += local[b];
    }
  }
  return true;
}

}  // namespace

// -----------------------------------------------------------------------------------------------------
// IN-PLANE RDF
// -----------------------------------------------------------------------------------------------------

/**
 * @details Calculates the in-plane RDF for quasi-two-dimensional water, when
 * both the atoms are of the same type. The input PointCloud only has particles
 * of type A in it.
 * This is registered as a Lua function and is
 * accessible to the user.
 * Internally, this function calls the following functions:
 * - rdf2::getSystemLengths (Gets the dimensions of the quasi-two-dimensional
 * system).
 * - ring::rdf2::sampleRDF_AA (Samples the current frame, binning the
 * coordinates).
 * - ring::rdf2::normalizeRDF (Normalizes the RDF).
 * - ring::sout::printRDF (Writes out the RDF to the desired output directory,
 * in the form of an ASCII file)
 *  @param[in] path The file path of the output directory to which output files
 * will be written.
 *  @param[in] rdfValues Vector containing the RDF values.
 *  @param[in] yCloud The input PointCloud.
 *  @param[in] cutoff Cutoff for the RDF. This should not be greater than half
 * the box length.
 *  @param[in] binwidth Width of the bin for histogramming.
 *  @param[in] firstFrame The first frame for RDF binning.
 *  @param[in] finalFrame The final frame for RDF binning.
 */
int rdf2::rdf2Danalysis_AA(
    std::string path, std::vector<double> &rdfValues,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
    double binwidth, int firstFrame, int finalFrame) {
  //
  int nopA = yCloud.nop; // Number of particles of type A in the current frame.
  int currentFrame = yCloud.currentFrame; // The current frame
  int nbin;                                // Number of bins
  std::vector<double> volumeLengths;       // Lengths of the volume of the
                                           // quasi-two-dimensional system
  std::vector<int> histogram; // Histogram for the RDF, at every step
  int nIter = finalFrame - firstFrame + 1; // Number of iterations
  // ----------------------------------------------
  // INITIALIZATION
  if (currentFrame == firstFrame) {
    // -----------------
    // Checks and balances?
    // -----------------
    nbin = static_cast<int>(cutoff / binwidth);
    if (nbin < 1) {
      nbin = 1;
    }
    rdfValues.resize(nbin);
  }
  volumeLengths = rdf2::getSystemLengths(yCloud);
  nbin = static_cast<int>(cutoff / binwidth);
  if (nbin < 1) {
    nbin = 1;
  }
  if (static_cast<int>(rdfValues.size()) != nbin) {
    rdfValues.assign(nbin, 0.0);
  }
  // Sample for the current frame!!
  histogram = rdf2::sampleRDF_AA(yCloud, cutoff, binwidth, nbin);
  // ----------------------------------------------
  // NORMALIZATION
  rdf2::normalizeRDF(nopA, rdfValues, histogram, binwidth, nbin, volumeLengths,
                     nIter);
  // ----------------------------------------------
  // PRINT OUT
  if (currentFrame == finalFrame) {
    // Create folder if required
    sout::makePath(path);
    std::string outputDirName = path + "topoMonolayer";
    sout::makePath(outputDirName);
    //
    std::string fileName = path + "topoMonolayer/rdf.dat";
    //
    // Comment line
    std::ofstream outputFile; // For the output file
    outputFile.open(fileName, std::ios_base::app | std::ios_base::out);
    outputFile << "# r  g(r)\n";
    outputFile.close();
    //
    //
    sout::printRDF(fileName, rdfValues, binwidth, nbin);
  } // end of print out
  // ----------------------------------------------

  return 0;
} // end of function

/**
 * @details Samples the RDF for a particular frame
 *  The input PointCloud only has particles
 *  of type A in it.
 *  - gen::periodicDist (Periodic distance between a pair of atoms).
 * @param[in] yCloud The input PointCloud.
 * @param[in] cutoff Cutoff for the RDF calculation.
 * @param[in] binwidth Width of the bin.
 * @param[in] nbin Number of bins.
 * @return RDF histogram for the current frame. Each pair contributes
 *  once, at its minimum-image distance.
 */
std::vector<int>
rdf2::sampleRDF_AA(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                   double cutoff, double binwidth, int nbin) {
  //
  std::vector<int> histogram;

  // Init the histogram to 0
  histogram.resize(nbin);

  const gen::FracBox frame = gen::makeFracBox(yCloud);
  double edge[3] = {frame.lx, frame.ly, frame.lz};
  std::sort(edge, edge + 3);
  const double pad = std::max(cutoff, 1.0) * 1e-8;
  // Below half the shortest edge a pair has one image, and that image is
  // the distance vesin already returns. A thinner cell can still use the
  // cell list: the histogram keeps one minimum-image distance. Past the
  // middle edge the candidate list is most of the pairs, so the direct
  // loop is faster.
  const bool oneImage =
      frame.ok && edge[0] > 0.0 && cutoff + pad < frame.halfMin;
  const bool sparse = oneImage && edge[1] > 0.0 && cutoff < 0.35 * edge[1];
  if (oneImage && yCloud.nop > 1 &&
      histogramPackedGrid(frame, yCloud, cutoff, binwidth, nbin, histogram)) {
    return histogram;
  }
#ifdef SEAMS_HAS_VESIN
  if (sparse && yCloud.nop > 1 && yCloud.box.size() >= 3) {
    std::vector<std::array<double, 3>> positions(
        static_cast<size_t>(yCloud.nop));
    for (int i = 0; i < yCloud.nop; i++) {
      positions[static_cast<size_t>(i)] = {
          yCloud.pts[i].x, yCloud.pts[i].y, yCloud.pts[i].z};
    }
    double box[3][3];
    double origin[3];
    nneigh::dumpBoundsToH(yCloud.box, yCloud.boxLow, box, origin);
    bool periodic[3] = {true, true, true};
    VesinOptions options{};
    options.cutoff = cutoff + pad;
    options.full = true;
    options.sorted = false;
    options.algorithm = VesinAutoAlgorithm;
    options.return_shifts = false;
    options.return_distances = true;
    options.return_vectors = false;
    VesinNeighborList neighbors;
    const char *error_message = nullptr;
    VesinDevice device = {VesinCPU, 0};
    const int status = vesin_neighbors(
        reinterpret_cast<const double (*)[3]>(positions.data()),
        static_cast<size_t>(yCloud.nop), box, periodic, device, options,
        &neighbors, &error_message);
    if (status == 0) {
      if (oneImage && neighbors.distances != nullptr) {
        for (size_t k = 0; k < neighbors.length; k++) {
          const int iatom = static_cast<int>(neighbors.pairs[k][0]);
          const int jatom = static_cast<int>(neighbors.pairs[k][1]);
          if (iatom >= jatom) {
            continue;
          }
          const double r = neighbors.distances[k];
          if (r > cutoff) {
            continue;
          }
          int b = static_cast<int>(r / binwidth);
          if (b < 0) {
            continue;
          }
          if (b >= nbin) {
            b = nbin - 1;
          }
          histogram[static_cast<size_t>(b)] += 2;
        }
      } else {
        std::vector<std::pair<int, int>> uniq;
        uniq.reserve(neighbors.length / 2 + 1);
        for (size_t k = 0; k < neighbors.length; k++) {
          int iatom = static_cast<int>(neighbors.pairs[k][0]);
          int jatom = static_cast<int>(neighbors.pairs[k][1]);
          if (iatom == jatom) {
            continue;
          }
          if (iatom > jatom) {
            std::swap(iatom, jatom);
          }
          uniq.emplace_back(iatom, jatom);
        }
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        for (const auto &pair : uniq) {
          const double r = std::sqrt(
              gen::periodicDistSq(frame, yCloud, pair.first, pair.second));
          if (r > cutoff) {
            continue;
          }
          int b = static_cast<int>(r / binwidth);
          if (b < 0) {
            continue;
          }
          if (b >= nbin) {
            b = nbin - 1;
          }
          histogram[static_cast<size_t>(b)] += 2;
        }
      }
      vesin_free(&neighbors);
      return histogram;
    }
    std::cerr << "Vesin failed: "
              << (error_message ? error_message : "unknown")
              << "; falling back to brute force.\n";
    vesin_free(&neighbors);
  }
#endif

  const double cut2 = cutoff * cutoff;
  const int nop = yCloud.nop;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (nop >= 512)
  {
    std::vector<int> local(static_cast<std::size_t>(nbin), 0);
#pragma omp for schedule(static)
    for (int iatom = 0; iatom < nop - 1; iatom++) {
      const auto &pi = yCloud.pts[static_cast<std::size_t>(iatom)];
      for (int jatom = iatom + 1; jatom < nop; jatom++) {
        const auto &pj = yCloud.pts[static_cast<std::size_t>(jatom)];
        double rij = 0.0;
        if (oneImage && frame.ok) {
          const double r2 =
              gen::fracDistSq(frame, pi.x, pi.y, pi.z, pj.x, pj.y, pj.z);
          if (r2 > cut2) {
            continue;
          }
          rij = std::sqrt(r2);
        } else {
          rij = std::sqrt(gen::periodicDistSq(frame, yCloud, iatom, jatom));
          if (rij > cutoff) {
            continue;
          }
        }
        int b = static_cast<int>(rij / binwidth);
        if (b < 0) {
          continue;
        }
        if (b >= nbin) {
          b = nbin - 1;
        }
        local[static_cast<std::size_t>(b)] += 2;
      }
    }
#pragma omp critical
    {
      for (int b = 0; b < nbin; b++) {
        histogram[static_cast<std::size_t>(b)] +=
            local[static_cast<std::size_t>(b)];
      }
    }
  }
#else
  for (int iatom = 0; iatom < nop - 1; iatom++) {
    const auto &pi = yCloud.pts[static_cast<std::size_t>(iatom)];
    for (int jatom = iatom + 1; jatom < nop; jatom++) {
      const auto &pj = yCloud.pts[static_cast<std::size_t>(jatom)];
      double rij = 0.0;
      if (oneImage && frame.ok) {
        const double r2 =
            gen::fracDistSq(frame, pi.x, pi.y, pi.z, pj.x, pj.y, pj.z);
        if (r2 > cut2) {
          continue;
        }
        rij = std::sqrt(r2);
      } else {
        rij = std::sqrt(gen::periodicDistSq(frame, yCloud, iatom, jatom));
        if (rij > cutoff) {
          continue;
        }
      }
      int b = static_cast<int>(rij / binwidth);
      if (b < 0) {
        continue;
      }
      if (b >= nbin) {
        b = nbin - 1;
      }
      histogram[static_cast<std::size_t>(b)] += 2;
    }
  }
#endif

  // Return the histogram
  return histogram;

} // end of function

/**
 * @details Normalizes the histogram and adds it to the RDF.
 * The normalization requires the plane area and height.
 * @param[in] nopA The number of particles of type A.
 * @param[in] rdfValues Radial distribution function values for all the frames,
 *  in the form of a vector.
 * @param[in] histogram The histogram for the current frame.
 * @param[in] binwidth The width of each bin for the RDF histogram.
 * @param[in] nbin The number of bins for the RDF.
 * @param[in] volumeLengths The confining dimensions of the
 *  quasi-two-dimensional system, which may be the slice dimensions or the
 *  dimensions of the box.
 * @param[in] nIter The number of iterations for which the coordinates will be
 *  binned. This is basically equivalent to the number of frames over which the
 * RDF will be calculated.
 */
int rdf2::normalizeRDF(int nopA, std::vector<double> &rdfValues,
                       std::vector<int> histogram, double binwidth, int nbin,
                       std::vector<double> volumeLengths, int nIter) {
  // volumeLengths is the particle AABB (box.size()==3) or dump H
  // {lx, ly, lz} (box.size()>=6). Product of the dump lengths is
  // nneigh::dumpVolume; Lx*Ly is the restricted-triclinic in-plane
  // area and Lz is the confined height. A vanishing in-plane product
  // (AABB of a line or point) keeps the smallest extent as height.
  const double planeArea = volumeLengths[0] * volumeLengths[1];
  const double volume =
      volumeLengths[0] * volumeLengths[1] * volumeLengths[2];
  const double height =
      planeArea > 0.0
          ? volume / planeArea
          : *std::min_element(volumeLengths.begin(), volumeLengths.end());
  double r;         // Distance for the current bin
  double factor;    // Factor for accounting for the two-dimensional slab
  double binVolume; // Volume of the current bin
  double volumeDensity = nopA / volume;

  // Loop through every bin
  for (int ibin = 0; ibin < nbin; ibin++) {
    //
    r = binwidth * (ibin + 0.5);
    // ----------------------
    // Factor calculation
    factor = 1;
    if (r > height) {
      // FACTOR = HEIGHT/(2*R)
      factor = height / (2 * r);
    } // r>height
    else {
      // FACTOR = 1 - R/(2*HEIGHT)
      factor = 1 - r / (2 * height);
    } // r<=height
    // ----------------------
    // binVolume = 4.0*PI*(DELTA_R**3)*((I_BIN+1)**3-I_BIN**3)/3.0
    binVolume = 4 * gen::pi * pow(binwidth, 3.0) *
                (pow(ibin + 1, 3.0) - pow(ibin, 3.0)) / 3.0;
    // Update the RDF
    rdfValues[ibin] +=
        histogram[ibin] / (nIter * binVolume * nopA * volumeDensity * factor);
  } // end of loop through every bin

  return 0;

} // end of function

/**
 * @details Calculates the lengths of the quasi-two-dimensional
 *  system. A tilt dump (box.size()>=6) returns dump H {lx, ly, lz}.
 *  A length-3 box returns the particle AABB; the smallest length is
 *  then the slab height.
 * @param[in] yCloud The molSys::PointCloud struct for the system.
 * @return {lx, ly, lz} from dumpBoundsToH, or the AABB extents.
 */
std::vector<double> rdf2::getSystemLengths(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  // Dump tilt: lattice edges from dumpBoundsToH, not the particle
  // AABB and not the bound-span product.
  if (yCloud.box.size() >= 6) {
    double H[3][3];
    double origin[3];
    nneigh::dumpBoundsToH(yCloud.box, yCloud.boxLow, H, origin);
    return {H[0][0], H[1][1], H[2][2]};
  }

  std::vector<double> lengths; // Volume lengths
  std::vector<double>
      rMax; // Max of the coordinates {0 is for x, 1 is for y and 2 is for z}
  std::vector<double>
      rMin; // Min of the coordinates {0 is for x, 1 is for y and 2 is for z}
  std::vector<double> r_iatom; // Current point coordinates
  int dim = 3;

  // Init
  r_iatom.push_back(yCloud.pts[0].x);
  r_iatom.push_back(yCloud.pts[0].y);
  r_iatom.push_back(yCloud.pts[0].z);
  rMin = r_iatom;
  rMax = r_iatom;

  // Loop through the entire PointCloud
  for (int iatom = 1; iatom < yCloud.nop; iatom++) {
    r_iatom[0] = yCloud.pts[iatom].x; // x coordinate of iatom
    r_iatom[1] = yCloud.pts[iatom].y; // y coordinate of iatom
    r_iatom[2] = yCloud.pts[iatom].z; // z coordinate of iatom
    // Loop through every dimension
    for (int k = 0; k < dim; k++) {
      // Update rMin
      if (r_iatom[k] < rMin[k]) {
        rMin[k] = r_iatom[k];
      } // end of rMin update
      // Update rMax
      if (r_iatom[k] > rMax[k]) {
        rMax[k] = r_iatom[k];
      } // end of rMax update
    }   // end of looping through every dimension
  }     // end of loop through all the atoms

  // Get the lengths
  for (int k = 0; k < dim; k++) {
    lengths.push_back(rMax[k] - rMin[k]);
  } // end of updating lengths

  return lengths;
} // end of function

/**
 * @details Calculates the plane area from the volume lengths.
 *  This is the product of the two largest dimensions of the
 * quasi-two-dimensional system.
 * @param[in] volumeLengths A vector of the lengths of the volume slice or
 *  simulation box
 * @return The plane area of the two significant dimensions
 */
double rdf2::getPlaneArea(std::vector<double> volumeLengths) {
  //
  // Sort the vector in descending order
  std::sort(volumeLengths.begin(), volumeLengths.end());
  std::reverse(volumeLengths.begin(), volumeLengths.end());

  return volumeLengths[0] * volumeLengths[1];

} // end of function
