#ifndef SEAMS_PHASE_H_
#define SEAMS_PHASE_H_

#include <mol_sys.hpp>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace phase {

enum class GlassKind { other = 0, ice, lda, mda, hda };

//! Open hexagonal channel along z: six-rings stacked, not a closed 512.
int openChannelCount(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList);

//! 64-bit key of donated H directions per oxygen (molID-shared hydrogens).
std::uint64_t protonKey(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    int oxygenType, int hydrogenType);

//! Mean squared displacement of hydrogenType between two frames.
double hydrogenMSD(
    const molSys::PointCloud<molSys::Point<double>, double> &frame0,
    const molSys::PointCloud<molSys::Point<double>, double> &frame1,
    int hydrogenType);

//! Ice XXI library: 152 sites in a BCT box a=b=20.197 A, c=7.891 A.
struct IceXXIHit {
  bool match = false;
  int nSites = 0;
  double a = 0.0;
  double c = 0.0;
  double density = 0.0;
};

IceXXIHit iceXXILibrary(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud);

//! Local number density (1/A^3) in a sphere of radius rcut.
std::vector<double> localDensity(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double rcut);

//! HDA/MDA/LDA/ice from local density. Ice Ih sits below iceMax.
//! The dense null sits at or above mdaMin.
GlassKind glassFromDensity(double rho, double iceMax, double mdaMin);

} // namespace phase

#endif
