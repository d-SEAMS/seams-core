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

#include <rdf2d.hpp>

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
    std::string path, std::vector<double> *rdfValues,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, double cutoff,
    double binwidth, int firstFrame, int finalFrame) {
  //
  int nopA = yCloud->nop; // Number of particles of type A in the current frame.
  int currentFrame = yCloud->currentFrame; // The current frame
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
    nbin = 1 + int(cutoff / binwidth); // Number of bins
    (*rdfValues).resize(nbin);         // RDF initialized to 0
  }                                    // end of initialization
  // ----------------------------------------------
  // SAMPLING
  // Get the volume lengths of the quasi-two-dimensional system
  volumeLengths = rdf2::getSystemLengths(yCloud);
  // Update cutoff if required? later
  nbin = int(cutoff / binwidth); // Number of bins
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
 * @param[in] cutoff Cutoff for the RDF calculation, which should be less than
 *  or equal to half the box length.
 * @param[in] binwidth Width of the bin.
 * @param[in] nbin Number of bins.
 * @return RDF histogram for the current frame.
 */
std::vector<int>
rdf2::sampleRDF_AA(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                   double cutoff, double binwidth, int nbin) {
  //
  std::vector<int> histogram; // Histogram for the RDF
  double r_ij;                // Unwrapped distance between iatom and jatom
  int ibin; // Index of the bin in which the distance r_ij falls

  // Init the histogram to 0
  histogram.resize(nbin);

  // Loop through pairs of atoms
  for (int iatom = 0; iatom < yCloud->nop - 1; iatom++) {
    //
    // Loop through the next atom in the atom pair
    for (int jatom = iatom + 1; jatom < yCloud->nop; jatom++) {
      // Calculate the distance between iatom and jatom
      r_ij = gen::periodicDist(yCloud, iatom, jatom);
      // Check if the distance is within the cutoff. Update the histogram if
      // r_ij is within the cutoff
      if (r_ij <= cutoff) {
        ibin = int(r_ij / binwidth); // Bin in which r_ij falls
        histogram[ibin] += 2;        // Update for iatom and jatom both
      }                              // end of histogram update
    }                                // end of loop through jatom
  }                                  // end of loop through iatom

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
int rdf2::normalizeRDF(int nopA, std::vector<double> *rdfValues,
                       std::vector<int> histogram, double binwidth, int nbin,
                       std::vector<double> volumeLengths, int nIter) {
  //
  auto height = *std::min_element(
      volumeLengths.begin(),
      volumeLengths.end()); // The height or the smallest dimension
  // double planeArea = rdf2::getPlaneArea(
  //     volumeLengths);  // The area of the quasi-two-dimensional water
  double r;         // Distance for the current bin
  double factor;    // Factor for accounting for the two-dimensional slab
  double binVolume; // Volume of the current bin
  double volumeDensity =
      nopA / (volumeLengths[0] * volumeLengths[1] * volumeLengths[2]);

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
    (*rdfValues)[ibin] +=
        histogram[ibin] / (nIter * binVolume * nopA * volumeDensity * factor);
  } // end of loop through every bin

  return 0;

} // end of function

/**
 * @details Calculates the lengths of the quasi-two-dimensional
 *  system. The smallest length is the 'height'.
 * @param[in] yCloud The molSys::PointCloud struct for the system.
 * @return The length (i.e. the 'height') of the smallest dimension of
 *  quasi-two-dimensional system.
 */
std::vector<double> rdf2::getSystemLengths(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  //
  std::vector<double> lengths; // Volume lengths
  std::vector<double>
      rMax; // Max of the coordinates {0 is for x, 1 is for y and 2 is for z}
  std::vector<double>
      rMin; // Min of the coordinates {0 is for x, 1 is for y and 2 is for z}
  std::vector<double> r_iatom; // Current point coordinates
  int dim = 3;

  // Init
  r_iatom.push_back(yCloud->pts[0].x);
  r_iatom.push_back(yCloud->pts[0].y);
  r_iatom.push_back(yCloud->pts[0].z);
  rMin = r_iatom;
  rMax = r_iatom;

  // Loop through the entire PointCloud
  for (int iatom = 1; iatom < yCloud->nop; iatom++) {
    r_iatom[0] = yCloud->pts[iatom].x; // x coordinate of iatom
    r_iatom[1] = yCloud->pts[iatom].y; // y coordinate of iatom
    r_iatom[2] = yCloud->pts[iatom].z; // z coordinate of iatom
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
