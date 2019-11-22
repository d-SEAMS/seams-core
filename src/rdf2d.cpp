#include <rdf2d.hpp>

// -----------------------------------------------------------------------------------------------------
// IN-PLANE RDF
// -----------------------------------------------------------------------------------------------------

/********************************************/ /**
 *  Calculates the in-plane RDF for quasi-two-dimensional water, when
 both the atoms are of the same type. The input PointCloud only has particles
 of type A in it.
 This is registered as a Lua function and is
 accessible to the user.
 Internally, this function calls the following functions:
 - ring::getSingleRingSize (Saves rings of a single ring size into a new vector
 of vectors, which is subsequently used for finding DDCs, HCs etc).
 - ring::findDDC (Finds the DDCs).
 - ring::findHC (Finds the HCs).
 - ring::findMixedRings (Finds the mixed rings, which are shared by DDCs and HCs
 both).
 - ring::getStrucNumbers (Gets the number of structures (DDCs, HCs, mixed rings,
 basal rings, prismatic rings, to be used for write-outs).
 - sout::writeTopoBulkData (Writes out the numbers and data obtained for the
 current frame).
 - ring::getAtomTypesTopoBulk (Gets the atom type for every atom, to be used for
 printing out the ice types found).
 - sout::writeLAMMPSdataTopoBulk (Writes out the atoms, with the classified
 types, into a LAMMPS data file, which can be visualized in OVITO).
 *  @param[in] path The file path of the output directory to which output files
 will be written.
 *  @param[in] rings Vector of vectors containing the primitive rings. This
 contains rings of all sizes.
 *  @param[in] nList Row-ordered neighbour list, by index.
 *  @param[in] yCloud The input PointCloud, with respect to which the indices in
 the rings and nList vector of vectors have been saved.
 *  @param[in] printEachCage Flag for printing the information of each cage in
 the frame (true) or not printing the coordinates/connectivity of each cage
 (false).
 ***********************************************/
int rdf2::rdf2Danalysis_AA(
    std::string path, std::vector<double> *rdfValues,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, double cutoff,
    double binwidth, int firstFrame, int finalFrame) {
  //
  int nopA =
      yCloud->nop;  // Number of particles of type A in the current frame.
  int currentFrame = yCloud->currentFrame;  // The current frame
  int nbin;                                 // Number of bins
  std::vector<double> volumeLengths;        // Lengths of the volume of the
                                            // quasi-two-dimensional system
  std::vector<int> histogram;  // Histogram for the RDF, at every step
  int nIter = finalFrame - firstFrame + 1;  // Number of iterations
  // ----------------------------------------------
  // INITIALIZATION
  if (currentFrame == firstFrame) {
    // -----------------
    // Checks and balances?
    // -----------------
    nbin = 1 + int(cutoff / binwidth);  // Number of bins
    (*rdfValues).resize(nbin);          // RDF initialized to 0
  }                                     // end of initialization
  // ----------------------------------------------
  // SAMPLING
  // Get the volume lengths of the quasi-two-dimensional system
  volumeLengths = rdf2::getSystemLengths(yCloud);
  // Update cutoff if required? later
  nbin = int(cutoff / binwidth);  // Number of bins
  // Sample for the current frame!!
  histogram = rdf2::sampleRDF_AA(yCloud, cutoff, binwidth, nbin);
  // ----------------------------------------------
  // NORMALIZATION
  rdf2::normalizeRDF(nopA, rdfValues, histogram, binwidth, nbin, volumeLengths,
                     nIter);
  // ----------------------------------------------
  // PRINT OUT
  if (currentFrame == finalFrame) {
    std::string fileName = path + "topoMonolayer/rdf.dat";
    sout::printRDF(fileName, rdfValues, binwidth, nbin);
  }  // end of print out
  // ----------------------------------------------

  return 0;
}  // end of function

/********************************************/ /**
 *  Samples the RDF for a particular frame
 The input PointCloud only has particles
 of type A in it.
 - sout::writeLAMMPSdataTopoBulk (Writes out the atoms, with the classified
 types, into a LAMMPS data file, which can be visualized in OVITO).
 *  @param[in] path The file path of the output directory to which output files
 will be written.
 *  @param[in] rings Vector of vectors containing the primitive rings. This
 contains rings of all sizes.
 *  @param[in] nList Row-ordered neighbour list, by index.
 *  @param[in] yCloud The input PointCloud, with respect to which the indices in
 the rings and nList vector of vectors have been saved.
 *  @param[in] printEachCage Flag for printing the information of each cage in
 the frame (true) or not printing the coordinates/connectivity of each cage
 (false).
 ***********************************************/
std::vector<int> rdf2::sampleRDF_AA(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, double cutoff,
    double binwidth, int nbin) {
  //
  std::vector<int> histogram;  // Histogram for the RDF
  double r_ij;                 // Unwrapped distance between iatom and jatom
  int ibin;  // Index of the bin in which the distance r_ij falls

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
        ibin = int(r_ij / binwidth);  // Bin in which r_ij falls
        histogram[ibin] += 2;         // Update for iatom and jatom both
      }                               // end of histogram update
    }                                 // end of loop through jatom
  }                                   // end of loop through iatom

  // Return the histogram
  return histogram;

}  // end of function

/********************************************/ /**
 *  Normalizes the histogram and adds it to the RDF.
 The normalization requires the plane area and height.
 *  @param[in] path The file path of the output directory to which output files
 will be written.
 *  @param[in] rings Vector of vectors containing the primitive rings. This
 contains rings of all sizes.
 *  @param[in] nList Row-ordered neighbour list, by index.
 *  @param[in] yCloud The input PointCloud, with respect to which the indices in
 the rings and nList vector of vectors have been saved.
 *  @param[in] printEachCage Flag for printing the information of each cage in
 the frame (true) or not printing the coordinates/connectivity of each cage
 (false).
 ***********************************************/
int rdf2::normalizeRDF(int nopA, std::vector<double> *rdfValues,
                       std::vector<int> histogram, double binwidth, int nbin,
                       std::vector<double> volumeLengths, int nIter) {
  //
  auto height = *std::min_element(
      volumeLengths.begin(),
      volumeLengths.end());  // The height or the smallest dimension
  // double planeArea = rdf2::getPlaneArea(
  //     volumeLengths);  // The area of the quasi-two-dimensional water
  double r;          // Distance for the current bin
  double factor;     // Factor for accounting for the two-dimensional slab
  double binVolume;  // Volume of the current bin
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
    }  // r>height
    else {
      // FACTOR = 1 - R/(2*HEIGHT)
      factor = 1 - r / (2 * height);
    }  // r<=height
    // ----------------------
    // binVolume = 4.0*PI*(DELTA_R**3)*((I_BIN+1)**3-I_BIN**3)/3.0
    binVolume = 4 * gen::pi * pow(binwidth, 3.0) *
                (pow(ibin + 1, 3.0) - pow(ibin, 3.0)) / 3.0;
    // Update the RDF
    (*rdfValues)[ibin] +=
        histogram[ibin] / (nIter * binVolume * nopA * volumeDensity * factor);
  }  // end of loop through every bin

  return 0;

}  // end of function

/********************************************/ /**
 *  Calculates the lengths of the quasi-two-dimensional
 system. The smallest length is the 'height'.
 *  @param[in] path The file path of the output directory to which output files
 will be written.
 *  @param[in] rings Vector of vectors containing the primitive rings. This
 contains rings of all sizes.
 *  @param[in] nList Row-ordered neighbour list, by index.
 *  @param[in] yCloud The input PointCloud, with respect to which the indices in
 the rings and nList vector of vectors have been saved.
 *  @param[in] printEachCage Flag for printing the information of each cage in
 the frame (true) or not printing the coordinates/connectivity of each cage
 (false).
 ***********************************************/
std::vector<double> rdf2::getSystemLengths(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  //
  std::vector<double> lengths;  // Volume lengths
  std::vector<double>
      rMax;  // Max of the coordinates {0 is for x, 1 is for y and 2 is for z}
  std::vector<double>
      rMin;  // Min of the coordinates {0 is for x, 1 is for y and 2 is for z}
  std::vector<double> r_iatom;  // Current point coordinates
  int dim = 3;

  // Init
  r_iatom.push_back(yCloud->pts[0].x);
  r_iatom.push_back(yCloud->pts[0].y);
  r_iatom.push_back(yCloud->pts[0].z);
  rMin = r_iatom;
  rMax = r_iatom;

  // Loop through the entire PointCloud
  for (int iatom = 1; iatom < yCloud->nop; iatom++) {
    r_iatom[0] = yCloud->pts[iatom].x;  // x coordinate of iatom
    r_iatom[1] = yCloud->pts[iatom].y;  // y coordinate of iatom
    r_iatom[2] = yCloud->pts[iatom].z;  // z coordinate of iatom
    // Loop through every dimension
    for (int k = 0; k < dim; k++) {
      // Update rMin
      if (r_iatom[k] < rMin[k]) {
        rMin[k] = r_iatom[k];
      }  // end of rMin update
      // Update rMax
      if (r_iatom[k] > rMax[k]) {
        rMax[k] = r_iatom[k];
      }  // end of rMax update
    }    // end of looping through every dimension
  }      // end of loop through all the atoms

  // Get the lengths
  for (int k = 0; k < dim; k++) {
    lengths.push_back(rMax[k] - rMin[k]);
  }  // end of updating lengths

  return lengths;
}  // end of function

/********************************************/ /**
 *  Calculates the plane area from the volume lengths.
 This is the product of the two largest dimensions of the quasi-two-dimensional
 water.
 *  @param[in] path The file path of the output directory to which output files
 will be written.
 *  @param[in] rings Vector of vectors containing the primitive rings. This
 contains rings of all sizes.
 *  @param[in] nList Row-ordered neighbour list, by index.
 *  @param[in] yCloud The input PointCloud, with respect to which the indices in
 the rings and nList vector of vectors have been saved.
 *  @param[in] printEachCage Flag for printing the information of each cage in
 the frame (true) or not printing the coordinates/connectivity of each cage
 (false).
 ***********************************************/
double rdf2::getPlaneArea(std::vector<double> volumeLengths) {
  //
  // Sort the vector in descending order
  std::sort(volumeLengths.begin(), volumeLengths.end());
  std::reverse(volumeLengths.begin(), volumeLengths.end());

  return volumeLengths[0] * volumeLengths[1];

}  // end of function