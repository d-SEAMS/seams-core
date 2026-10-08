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
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

// A pair inside the cutoff adds 2 at its distance.
inline void binPair(double r, double cutoff, double binwidth, int nbin,
                    int *histogram) {
  if (r > cutoff) {
    return;
  }
  const int b = std::min(static_cast<int>(r / binwidth), nbin - 1);
  if (b >= 0) {
    histogram[b] += 2;
  }
}

// Rapaport's cell pairs (The Art of Molecular Dynamics Simulation, 2004).
// Cells span at least the cutoff across each face separation, so inside the
// one-image ball a pair sits in one cell or in neighbouring ones. The half
// stencil visits each neighbouring pair of cells once, and one lattice
// shift serves the whole pair, so the inner loop does no wrap.
constexpr int kHalfStencil[13][3] = {
    {0, 0, 1},  {0, 1, -1}, {0, 1, 0},  {0, 1, 1},  {1, -1, -1},
    {1, -1, 0}, {1, -1, 1}, {1, 0, -1}, {1, 0, 0},  {1, 0, 1},
    {1, 1, -1}, {1, 1, 0},  {1, 1, 1}};

bool histogramCellPairs(
    const gen::FracBox &frame,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double cutoff, double binwidth, int nbin, std::vector<int> &histogram) {
  const int n = yCloud.nop;
  const double span = cutoff + std::max(cutoff, 1.0) * 1e-12;
  if (!frame.ok || n < 2 || !(cutoff > 0.0) || nbin <= 0 ||
      !(frame.wx > span && frame.wy > span && frame.wz > span)) {
    return false;
  }
  int dims[3] = {static_cast<int>(std::min(frame.wx / span, 1.0e6)),
                 static_cast<int>(std::min(frame.wy / span, 1.0e6)),
                 static_cast<int>(std::min(frame.wz / span, 1.0e6))};
  // A dilute frame takes coarser cells rather than walking empty ones.
  while (1LL * dims[0] * dims[1] * dims[2] > 2LL * n + 27) {
    int *widest = std::max_element(dims, dims + 3);
    *widest = std::max(1, *widest / 2);
  }
  const int ny = dims[1];
  const int nz = dims[2];
  const int ncell = dims[0] * ny * nz;
  const double lx = frame.lx, ly = frame.ly, lz = frame.lz;
  const double xy = frame.xy, xz = frame.xz, yz = frame.yz;

  std::vector<int> cellOf(static_cast<std::size_t>(n));
  std::vector<double> wrapped(3 * static_cast<std::size_t>(n));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 4096)
#endif
  for (int i = 0; i < n; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    const double fz = (p.z - frame.oz) / lz;
    const double fy = (p.y - frame.oy - yz * fz) / ly;
    const double fx = (p.x - frame.ox - xy * fy - xz * fz) / lx;
    double s[3] = {fx - std::floor(fx), fy - std::floor(fy),
                   fz - std::floor(fz)};
    int idx[3];
    for (int k = 0; k < 3; k++) {
      if (s[k] >= 1.0) {
        s[k] = 0.0;
      }
      idx[k] = std::min(static_cast<int>(s[k] * dims[k]), dims[k] - 1);
    }
    cellOf[static_cast<std::size_t>(i)] = (idx[0] * ny + idx[1]) * nz + idx[2];
    wrapped[3 * static_cast<std::size_t>(i)] = lx * s[0] + xy * s[1] + xz * s[2];
    wrapped[3 * static_cast<std::size_t>(i) + 1] = ly * s[1] + yz * s[2];
    wrapped[3 * static_cast<std::size_t>(i) + 2] = lz * s[2];
  }
  std::vector<int> start(static_cast<std::size_t>(ncell) + 1, 0);
  for (const int c : cellOf) {
    start[static_cast<std::size_t>(c) + 1]++;
  }
  for (int c = 0; c < ncell; c++) {
    start[static_cast<std::size_t>(c) + 1] += start[static_cast<std::size_t>(c)];
  }
  std::vector<int> fill(start.begin(), start.end() - 1);
  std::vector<double> px(static_cast<std::size_t>(n));
  std::vector<double> py(static_cast<std::size_t>(n));
  std::vector<double> pz(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    const auto s = static_cast<std::size_t>(
        fill[static_cast<std::size_t>(cellOf[static_cast<std::size_t>(i)])]++);
    px[s] = wrapped[3 * static_cast<std::size_t>(i)];
    py[s] = wrapped[3 * static_cast<std::size_t>(i) + 1];
    pz[s] = wrapped[3 * static_cast<std::size_t>(i) + 2];
  }

  const double *X = px.data();
  const double *Y = py.data();
  const double *Z = pz.data();
  const int *first = start.data();
  const double cut2 = cutoff * cutoff;
  const std::size_t stride = static_cast<std::size_t>(nbin + 15) / 16 * 16;
  int nthreads = 1;
#ifdef SEAMS_HAS_OPENMP
  if (n >= 512) {
    nthreads = omp_get_max_threads();
  }
#endif
  std::vector<int> locals(static_cast<std::size_t>(nthreads) * stride, 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
  {
    int tid = 0;
#ifdef SEAMS_HAS_OPENMP
    tid = omp_get_thread_num();
#endif
    int *local = locals.data() + static_cast<std::size_t>(tid) * stride;
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int c = 0; c < ncell; c++) {
      const int a0 = first[c];
      const int a1 = first[c + 1];
      if (a0 == a1) {
        continue;
      }
      for (int i = a0; i < a1; i++) {
        for (int j = i + 1; j < a1; j++) {
          const double dx = X[j] - X[i];
          const double dy = Y[j] - Y[i];
          const double dz = Z[j] - Z[i];
          const double r2 = dx * dx + dy * dy + dz * dz;
          if (r2 <= cut2) {
            binPair(std::sqrt(r2), cutoff, binwidth, nbin, local);
          }
        }
      }
      const int home[3] = {c / (ny * nz), (c / nz) % ny, c % nz};
      for (const auto &d : kHalfStencil) {
        int nb[3];
        int w[3];
        for (int k = 0; k < 3; k++) {
          nb[k] = home[k] + d[k];
          w[k] = nb[k] < 0 ? -1 : (nb[k] >= dims[k] ? 1 : 0);
          nb[k] -= w[k] * dims[k];
        }
        const int cell = (nb[0] * ny + nb[1]) * nz + nb[2];
        const int b0 = first[cell];
        const int b1 = first[cell + 1];
        const double sx = lx * w[0] + xy * w[1] + xz * w[2];
        const double sy = ly * w[1] + yz * w[2];
        const double sz = lz * w[2];
        for (int i = a0; i < a1; i++) {
          const double xi = X[i] - sx;
          const double yi = Y[i] - sy;
          const double zi = Z[i] - sz;
          for (int j = b0; j < b1; j++) {
            const double dx = X[j] - xi;
            const double dy = Y[j] - yi;
            const double dz = Z[j] - zi;
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 <= cut2) {
              binPair(std::sqrt(r2), cutoff, binwidth, nbin, local);
            }
          }
        }
      }
    }
  }
  for (int t = 0; t < nthreads; t++) {
    for (int b = 0; b < nbin; b++) {
      histogram[static_cast<std::size_t>(b)] +=
          locals[static_cast<std::size_t>(t) * stride + static_cast<std::size_t>(b)];
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

  // Inside the one-image ball every pair has a single image within the
  // cutoff, and the cell pairs find it. Past the ball the direct loop takes
  // the minimum image of each pair.
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  const double pad = std::max(cutoff, 1.0) * 1e-8;
  if (frame.ok && cutoff + pad < frame.halfMin &&
      histogramCellPairs(frame, yCloud, cutoff, binwidth, nbin, histogram)) {
    return histogram;
  }

  const int nop = yCloud.nop;
#ifdef SEAMS_HAS_MINIMAGE
  const gen::CellFrame cells = gen::makeCellFrame(yCloud);
#endif
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (nop >= 512)
#endif
  {
    std::vector<int> local(histogram.size(), 0);
    std::vector<double> row(static_cast<std::size_t>(std::max(nop, 1)));
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
    for (int iatom = 0; iatom < nop - 1; iatom++) {
      const auto rest = static_cast<std::size_t>(nop - iatom - 1);
      bool batched = false;
#ifdef SEAMS_HAS_MINIMAGE
      const std::size_t at = 3 * static_cast<std::size_t>(iatom);
      batched = cells.ok && frame.ok &&
                mi_dist2_euclidean_many(&cells.cell, &cells.xyz[at],
                                        &cells.xyz[at + 3], rest,
                                        row.data()) == 0;
#endif
      for (std::size_t k = 0; k < rest; k++) {
        const int jatom = iatom + 1 + static_cast<int>(k);
        const double r2 =
            batched ? row[k] : gen::periodicDistSq(frame, yCloud, iatom, jatom);
        binPair(std::sqrt(r2), cutoff, binwidth, nbin, local.data());
      }
    }
#ifdef SEAMS_HAS_OPENMP
#pragma omp critical
#endif
    for (std::size_t b = 0; b < local.size(); b++) {
      histogram[b] += local[b];
    }
  }

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
