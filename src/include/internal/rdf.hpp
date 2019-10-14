#ifndef __RDF_H_
#define __RDF_H_

#include <fstream>
#include <iostream>
#include <memory>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <sstream>

namespace correl {

// Value of pi
const double PI = 3.14159;

// Struct for the RDF
struct PairCorrel {
  std::vector<double> histo; // vector for the histogram
  double binwidth;           // Binwidth
  int nframes = 0;           // Number of snapshots processed
  int n_iatoms = 1;          // Number of atoms of type I
  int n_jatoms = 1;          // Number of atoms of type J
};

// Function for accumulating the RDF for a frame, when the central atom is of the same type as
// the distribution atom
PairCorrel
accumulateRDFii(MolSys::PointCloud<MolSys::Point<double>, double> *yCloud,
                PairCorrel *rdf, double binwidth, double cutoff, int typeI,
                bool isSlice = false,
                std::array<double, 3> = std::array<double, 3>{0, 0, 0},
                std::array<double, 3> = std::array<double, 3>{0, 0, 0});

// Function for normalizing the RDF for a frame, when the central atom is of the same type as
// the distribution atom
PairCorrel
normalizeRDF(MolSys::PointCloud<MolSys::Point<double>, double> *yCloud,
             PairCorrel *rdf, bool isSlice = false,
             std::array<double, 3> = std::array<double, 3>{0, 0, 0},
             std::array<double, 3> = std::array<double, 3>{0, 0, 0});

// Function for clearing large vectors in the PairCorrel struct
PairCorrel clearPairCorrel(PairCorrel *rdf);

// Function for pretty printing? a struct of PairCorrel
int prettyPrintRDF(PairCorrel *rdf, std::string outFile);

} // namespace correl

#endif // __RDF_H_
