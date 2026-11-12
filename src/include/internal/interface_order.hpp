//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#ifndef SEAMS_INTERFACE_ORDER_H_
#define SEAMS_INTERFACE_ORDER_H_

#include <mol_sys.hpp>

#include <cstdint>
#include <string>
#include <vector>

/** @file interface_order.hpp
 *  @brief The two published CHILL+ rules on an interface, and the ring
 *  labels of confined ice.
 *
 *  Thresholds are the 2015 values. Nothing here is fitted.
 *
 *  Nguyen and Molinero, J. Phys. Chem. B 119 (2015) 9369, equation 1
 *  and Table 1: q3 is the mean of Y_3m over the four nearest neighbours,
 *  a staggered bond is c(i,j) <= -0.8, an eclipsed bond is
 *  -0.35 <= c(i,j) <= 0.25, a crystal row has four neighbours inside
 *  the cutoff, and an interfacial row keeps the neighbour test in
 *  footnote a of that table.
 *
 *  OVITO 3.17.2 Chill+ (manual page chill_plus) and the mdapy kernel
 *  that was cross-checked against OVITO's ChillPlusModifier: q3 runs
 *  over every neighbour inside the cutoff, a staggered bond is
 *  c(i,j) < -0.8, an eclipsed bond is -0.35 < c(i,j) < 0.25, and an
 *  interfacial row is two staggered bonds, or three staggered and no
 *  eclipsed bond, with no neighbour test.
 */

namespace iface {

enum class Cause { agree = 0, neighbour = 1, q3Shell = 2, openEdge = 3, mixed = 4 };

enum class Kind { slab = 0, nucleus = 1 };

struct Front {
  double position = 0; //! 0.5 crossing, or the nucleus radius
  double width = 0;    //! 10-90 width along the same path
  bool defined = false;
};

struct Profile {
  double crystal = 0;      //! (cubic + hexagonal) / N
  double interfacial = 0;  //! interfacial / N
  double solid = 0;        //! crystal + interfacial, over N
  double thickness = 0;    //! ice span between the two 0.5 crossings, or 2 R
  std::vector<Front> fronts;
  bool defined = false;
};

struct Bin {
  double mid = 0;
  int n = 0;
  int disagree = 0;
  double crystalPub = 0;
  double solidPub = 0;
  double crystalOvi = 0;
  double solidOvi = 0;
};

struct Study {
  int n = 0;
  int nDisagree = 0;
  int causeCount[5] = {};
  int countPub[9] = {};
  int countOvi[9] = {};
  int transition[9][9] = {};
  int shellDisagree[8] = {}; //! 0..6 and 7 meaning 7 or more
  int shellAgree[8] = {};
  int ringDisagree[3] = {}; //! atoms in at least one ring of size 4, 5, 6
  int ringAgree[3] = {};
  int prismAtoms = 0;
  int prismBlocks = 0;
  int bilayerAtoms = 0;
  int polygonAtoms = 0;
  int thinSheets = 0;
  int diamondOther = 0;
  int diamondCubic = 0;
  int diamondHex = 0;
  int cnaOther = 0;
  int cnaCrystal = 0; //! FCC, HCP, icosahedral or BCC under adaptive CNA
  int ptmPass = 0;
  int ptmPassDisagree = 0;
  double q6Disagree = 0;
  double q6Agree = 0;
  double depthDisagree = 0; //! mean distance to the nearest front
  double depthAgree = 0;
  int depthDisagreeN = 0;
  int depthAgreeN = 0;
  Kind kind = Kind::nucleus;
  int axis = 2;
  Profile solidPub;    //! cubic + hexagonal + interfacial
  Profile crystalPub;  //! cubic + hexagonal
  Profile solidOvi;
  Profile crystalOvi;
  std::vector<Bin> bins;
};

//! Both rules, the ring census, and the interface or nucleus profile.
//! Cutoff and bond thresholds stay at the published values.
//! census false skips rings, sheets, diamond, CNA, q6 and the template.
//! The CHILL+ labels and the thickness profile are unchanged.
[[nodiscard]] Study studyCloud(
    molSys::PointCloud<molSys::Point<double>, double> cloud, double cutoff = 3.5,
    double binWidth = 1.0, bool census = true);

[[nodiscard]] Study studyDump(const std::string &path, int frame, double cutoff = 3.5,
                              double binWidth = 1.0, bool census = true);

[[nodiscard]] std::string formatStudy(const Study &study, const std::string &tag);

//! Ordinary least squares slope. Empty or a single point returns 0.
[[nodiscard]] double slope(const std::vector<double> &x, const std::vector<double> &y);

struct GrowthFrame {
  int frame = 0;
  double timePs = 0;
  Study study;
};

struct Growth {
  std::vector<GrowthFrame> frames;
  double rateSolidPub = 0; //! d(thickness)/dt in A/ps
  double rateSolidOvi = 0;
  double rateCrystalPub = 0;
  double rateCrystalOvi = 0;
};

//! Every frame of a dump. time of frame i is (i - 1) * dumpEveryPs.
[[nodiscard]] Growth studyTrajectory(const std::string &path, double dumpEveryPs,
                                     double cutoff = 3.5);

[[nodiscard]] std::string formatGrowth(const Growth &growth);

//! Lechner-Dellago q6 on the 3.5 A shell, the counting rule in Espinosa,
//! Navarro, Sanz, Valeriani and Vega, J. Chem. Phys. 145, 211922 (2016),
//! page 211922-3. The threshold is not fitted. The report is the fraction
//! of molecules at or above each listed q6bar and local q6, plus the
//! frame-1 histogram of the ice core and the liquid core.
[[nodiscard]] std::string formatQ6Growth(const std::string &path, double dumpEveryPs,
                                         double cutoff = 3.5);

//! Ideal oxygen frameworks at a nearest-neighbour distance of `bond`.
//! The box is orthorhombic and large enough that a 3.5 A shell does not
//! see a lateral image.
[[nodiscard]] molSys::PointCloud<molSys::Point<double>, double>
nanotube(int ring, int stacks, double bond = 2.75);

[[nodiscard]] molSys::PointCloud<molSys::Point<double>, double>
honeycombSheet(int nx, int ny, int layers, double bond = 2.75);

//! 10-90 width of a step smoothed by a Gaussian of standard deviation xi.
//! The profile is the normal cdf, so the width is 2 * 1.2815515655446004 * xi.
[[nodiscard]] double gaussianStepWidth(double xi);

struct BootstrapInterval {
  double mean = 0;
  double lo = 0; //! 2.5 percentile, Hyndman-Fan type 7
  double hi = 0; //! 97.5 percentile
  bool excludesZero = false;
};

//! Percentile bootstrap of the mean. nDraw samples with replacement.
[[nodiscard]] BootstrapInterval bootstrapMean(const std::vector<double> &samples, int nDraw,
                                              std::uint64_t seed);

//! Mean of a[i] - b[i], resampled by shared index.
[[nodiscard]] BootstrapInterval pairedDifference(const std::vector<double> &a,
                                                 const std::vector<double> &b, int nDraw,
                                                 std::uint64_t seed);

struct WillardChandler {
  double xi = 2.4;
  double meanRadius = 0;
  double surfaceRms = 0;     //! rms of the instantaneous surface about its mean
  double densityWidth = 0;   //! 10-90 of the coarse density normal to that surface
  double solidWidthPub = 0;  //! 10-90 of the 2015 solid fraction on the same coordinate
  double solidWidthOvi = 0;
  int nRays = 0;
  bool defined = false;
};

//! Instantaneous interface of Willard and Chandler, J. Phys. Chem. B 114,
//! 1954 (2010). xi is their Gaussian width. The surface is the set of
//! rays from the crystal centre of mass on which the coarse density
//! crosses the midpoint of the core and the outer liquid.
[[nodiscard]] WillardChandler willardChandler(
    molSys::PointCloud<molSys::Point<double>, double> cloud, double xi = 2.4);

[[nodiscard]] std::string formatWillardChandler(const WillardChandler &surface);

} // namespace iface

#endif
