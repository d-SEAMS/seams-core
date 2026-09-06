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

//! Ice XXI library: Lee et al., Nat. Mater. 25, 302 (2026), I-4 2d,
//! Z=152, a=b=20.197 A, c=7.891 A, 1.413 g/cm^3. A hit also requires a
//! tetrahedral four-nearest graph (mean coordination in [3.5, 4.5]) and at
//! least 50 primitive six-rings. A simple-cubic packing of 152 sites in that
//! cell is not a hit.
struct IceXXIHit {
  bool match = false;
  int nSites = 0;
  int nSix = 0;
  double a = 0.0;
  double c = 0.0;
  double density = 0.0;
  double meanCoord = 0.0;
};

IceXXIHit iceXXILibrary(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud);

//! Local number density (1/A^3) in a sphere of radius rcut.
std::vector<double> localDensity(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double rcut);

//! Bulk mass density in g/cm^3 from nop and the box (18.015 g/mol).
double frameDensity(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud);

//! HDA/MDA/LDA/ice from local number density (1/A^3). Ice first-shell
//! density is four neighbours in the 3.5 A sphere (~0.022 A^-3). HDA is
//! the dense shell (>= 0.070 A^-3). ldaMax splits LDA from MDA.
GlassKind glassFromDensity(double rho, double iceMax, double ldaMax,
                           double mdaMin);

//! Three-argument form: LDA/MDA split at the midpoint of iceMax and mdaMin.
GlassKind glassFromDensity(double rho, double iceMax, double mdaMin);

//! Literature local-density windows: ice < 0.028, LDA < 0.040, MDA < 0.070,
//! else HDA.
GlassKind glassFromDensity(double rho);

} // namespace phase

#endif
