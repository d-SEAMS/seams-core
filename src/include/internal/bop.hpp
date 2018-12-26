#ifndef __BOP_H_
#define __BOP_H_

#include <array>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cmath>
#include <complex>
#include <generic.hpp>
#include <math.h>
#include <mol_sys.hpp>
#include <neighbours.hpp>

namespace chill {

// 2*l+1 length complex vector
struct YlmAtom {
  std::vector<std::complex<double>> ylm;
};

// Vector of 2*l+1 averaged over 4 nearest neighbours
struct QlmAtom {

  std::vector<YlmAtom> ptq; // Averaged over neighbours
};

// Uses Boost for spherical harmonics, and gets c_ij according to the CHILL algorithm
molSys::PointCloud<molSys::Point<double>, double>
getCorrel(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
          bool isSlice = false);

// Classifies each atom according to the CHILL algorithm
molSys::PointCloud<molSys::Point<double>, double>
getIceType(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
           bool isSlice = false, std::string outputFileName = "chill.txt");

// Gets c_ij and then classifies bond types according to the CHILL+ algorithm
molSys::PointCloud<molSys::Point<double>, double>
getCorrelPlus(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
              bool isSlice = false);

// Classifies each atom according to the CHILL+ algorithm
molSys::PointCloud<molSys::Point<double>, double>
getIceTypePlus(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
               bool isSlice = false,
               std::string outputFileName = "chillPlus.txt");

// q6 can distinguish between water and ice. Use this for the largest ice cluster
std::vector<double>
getq6(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
      bool isSlice = false);

// 'Test' condition for classifying hexagonal ice using averaged q6 and q3
// Checks water
// According to https://pubs.rsc.org/en/content/articlehtml/2011/cp/c1cp22167a
// Gets c_ij and then classifies bond types according to the CHILL+ algorithm
molSys::PointCloud<molSys::Point<double>, double>
reclassifyWater(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                std::vector<double> *q6);

// Gets a PointCloud struct of the ice particles in a given frame
molSys::PointCloud<molSys::Point<double>, double>
getIceCloud(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
            molSys::PointCloud<molSys::Point<double>, double> *iceCloud);

// Finds the largest ice cluster
int largestIceCluster(
    molSys::PointCloud<molSys::Point<double>, double> *iceCloud, double cutoff,
    bool printCluster = false, bool isSlice = false);

// Prints out the iceType for a particular frame onto the terminal
int printIceType(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 bool isSlice = false,
                 std::string outputFileName = "superChill.txt");

// Checks if a given iatom is interfacial ice or not, according to the CHILL algorithm
bool isInterfacial(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                   int iatom, int num_staggrd, int num_eclipsd);

// Finds the number of staggered bonds for a given atom of index jatom
int numStaggered(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 int jatom);

} // namespace chill

namespace sph {

// 7 is for Q3, orderL=3

std::vector<std::complex<double>> spheriHarmo(int orderL,
                                              std::array<double, 2> coordSph);

std::array<double, 2> radialCoord(std::array<double, 3>);

// Lookup table for Q3
std::vector<std::complex<double>>
lookupTableQ3Vec(std::array<double, 2> angles);

// Lookup table for Q3 (m=0 to m=6)
std::complex<double> lookupTableQ3(int m, std::array<double, 2> angles);

// Lookup table for Q6
std::vector<std::complex<double>>
lookupTableQ6Vec(std::array<double, 2> angles);

// Lookup table for Q6 (m=0 to m=12)
std::complex<double> lookupTableQ6(int m, std::array<double, 2> angles);

} // namespace sph

inline std::vector<std::complex<double>>
avgVector(std::vector<std::complex<double>> v, int l, int neigh) {
  if (neigh == 0) {
    return v;
  }
  for (int m = 0; m < 2 * l + 1; m++) {
    v[m] = (1.0 / (double)neigh) * v[m];
  }

  return v;
}

#endif // __BOP_H_
