#include <rdf.hpp>

// Function for accumulating the RDF for a frame, when the central atom is of the same type as
// the distribution atom
correl::PairCorrel correl::accumulateRDFii(
    MolSys::PointCloud<MolSys::Point<double>, double> *yCloud,
    correl::PairCorrel *rdf, double binwidth, double cutoff, int typeI,
    bool isSlice, std::array<double, 3> coordLow,
    std::array<double, 3> coordHigh) {
  int n_iatoms = 0;                  // Number of atoms of type I
  int nbin = int(cutoff / binwidth); // No. of bins
  int nnumNeighbours = 0; // Number of nearest neighbours within the cutoff
  int ibin = 0;           // Index of current bin number
  int jatom = 0;          // Index of neighbour
  double r_ij = 0;        // Unwrapped distance r_ij

  rdf->histo.resize(nbin);  // Resize the histo vector
  rdf->nframes += 1;        // Update the frame number
  rdf->binwidth = binwidth; // Update binwidth

  for (int iatom = 0; iatom < yCloud->nop - 1; iatom++) {
    // Check for the right type
    if (yCloud->pts[iatom].type != typeI) {
      continue;
    }
    n_iatoms++; // Increment the number of iatoms
    nnumNeighbours =
        yCloud->pts[iatom].neighList.size(); // Get number of neighbours
    if (nnumNeighbours == 0) {
      continue;
    }
    // Loop over half neighbour list
    for (int j = 0; j < nnumNeighbours; j++) {
      jatom = yCloud->pts[iatom].neighList[j];        // Index of neighbour
      r_ij = gen::periodicDist(yCloud, iatom, jatom); // Get unwrapped distance
      ibin = int(r_ij / binwidth); // Find which bin the particle falls in
      // Now add to histogram (+2 because this for both iatom and jatom)
      rdf->histo[ibin] += 2;
    } // Loop over nearest neighbours

  } // End of loop through all atoms

  rdf->n_iatoms = n_iatoms; // n_iatoms
  rdf->n_jatoms = n_iatoms; // n_iatoms

  return *rdf;
}

/********************************************/ /**
 *  Function for normalizing RDF
 ***********************************************/
correl::PairCorrel
correl::normalizeRDF(MolSys::PointCloud<MolSys::Point<double>, double> *yCloud,
                     correl::PairCorrel *rdf, bool isSlice,
                     std::array<double, 3> coordLow,
                     std::array<double, 3> coordHigh) {
  double bin_volume = 0; // Bin volume
  double nideal = 0;     // No. of ideal gas particles in each bin_volume
  double volume = yCloud->box[0] * yCloud->box[1] * yCloud->box[2];
  double rho = rdf->n_jatoms / volume; // Number density of distribution atoms J

  // Loop over all bins
  for (int ibin = 0; ibin < rdf->histo.size(); ibin++) {
    // Volume between bin k+1 and k
    bin_volume = (pow(ibin + 2, 3) - pow(ibin + 1, 3)) * pow(rdf->binwidth, 3);
    // Assuming the total nop does not change with time
    // Number of ideal gas particles in bin_volume
    nideal = (4.0 / 3.0) * correl::PI * bin_volume * rho;
    // Normalization
    rdf->histo[ibin] /= (rdf->nframes * rdf->n_iatoms * nideal);
  }

  return *rdf;
}

/********************************************/ /**
 *  Function for printing out info in PairCorrel struct
 ***********************************************/
int correl::prettyPrintRDF(correl::PairCorrel *rdf, std::string outFile) {
  std::ofstream outputFile;
  // Create a new file in the output directory
  outputFile.open(outFile);

  if (outputFile.is_open()) {
    // First line
    outputFile << "# r \t rdfValue \t over " << rdf->nframes << " frames\n";
    // Write out the arrays x and y to a file in the output folder
    for (int ibin = 0; ibin < rdf->histo.size(); ibin++) {
      outputFile << (rdf->binwidth * (ibin + 1.5)) << "\t" << rdf->histo[ibin]
                 << "\n";
    }
  }
  // Close the file
  outputFile.close();
  return 1;
}

/********************************************/ /**
 *  Function for clearing RDF PairCorrel struct
 ***********************************************/
correl::PairCorrel correl::clearPairCorrel(correl::PairCorrel *rdf) {

  std::vector<double> tempPts;

  tempPts.swap(rdf->histo);
  rdf->nframes = 0;

  return *rdf;
}
