#ifndef __RDF2D_H_
#define __RDF2D_H_

#include <math.h>
#include <sys/stat.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <mol_sys.hpp>
#include <order_parameter.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

/*! \file rdf2d.hpp
    \brief File containing functions used to calculate the
   in-plane radial distribution functions.
*/

/*!
 *  \addtogroup rdf2
 *  @{
 */

namespace rdf2 {

// Main function for calculating the RDF for the same type of particle: calls
// other functions for initializing, sampling and normalizing the RDF
int rdf2Danalysis_AA(std::string path, std::vector<double> *rdfValues,
                     molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                     double cutoff, double binwidth, int firstFrame,
                     int finalFrame);

// Samples the RDF histogram at every step
std::vector<int> sampleRDF_AA(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, double cutoff,
    double binwidth, int nbin);

// Normalize the histogram
int normalizeRDF(int nopA, std::vector<double> *rdfValues,
                 std::vector<int> histogram, double binwidth, int nbin,
                 std::vector<double> volumeLengths, int nIter);

// Gets the lengths of the volume slice of the quasi-two-dimensional system
std::vector<double> getSystemLengths(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud);

// Gets the plane area from the volume lengths vector
double getPlaneArea(std::vector<double> volumeLengths);

}  // namespace rdf2

#endif  // __RDF2D_H_
