//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//-----------------------------------------------------------------------------------

#ifndef SEAMS_SCIENCE_FRONT_H_
#define SEAMS_SCIENCE_FRONT_H_

#include <mol_sys.hpp>

#include <cstdint>
#include <string>
#include <vector>

/** @file science_front.hpp
 *  @brief Layer profiles and the single-trajectory critical size.
 *
 *  CHILL+ is the published Nguyen and Molinero rule (J. Phys. Chem. B
 *  119, 9369, 2015): q3 over the four nearest neighbours, a crystal row
 *  of four neighbours inside the cutoff, and the interfacial neighbour
 *  test in Table 1. The OVITO rule is applied separately and is not
 *  used to choose these labels. The shell label is the second-shell
 *  template: 0 hexagonal, 1 cubic, 2 stacking-disordered, and -1
 *  otherwise. q6-bar is the Lechner-Dellago average on the same 3.5 A
 *  cutoff. None of these thresholds is fitted.
 *
 *  The critical size on one spontaneous trajectory is the mean size of
 *  the largest cluster in the contiguous band around its minimum. The
 *  band stops where the size exceeds that minimum by the larger of 5
 *  molecules and 15 percent. The interval is a percentile bootstrap of
 *  that mean. A minimum in the first two frames, followed by a larger
 *  cluster, is an upper bound: the dump does not contain a smaller
 *  nucleus. It is not a committor.
 */

namespace science {

//! Published CHILL+ cutoff, in angstrom. Not a fitted parameter.
inline constexpr double kCutoffA = 3.5;

//! Lechner-Dellago q6-bar ice cut used for the cluster mask.
//! The same gap as the nucleation census, not refitted here.
inline constexpr double kQ6BarIceCut = 0.36;

struct Front {
  double position = 0; //! 0.5 crossing along the profile, angstrom
  double width = 0;    //! 10-90 width of that crossing
  bool defined = false;
};

struct LayerRow {
  double mid = 0;
  int n = 0;
  int cubic = 0;
  int hexagonal = 0;
  int interfacial = 0;
  int liquid = 0;
  int shellCubic = 0;
  int shellHex = 0;
  int shellStack = 0;
  int shellOther = 0;
  double q6 = 0;
  double q6bar = 0;
};

struct Layers {
  int axis = 2;
  double length = 0;
  int n = 0;
  std::vector<LayerRow> rows;
  std::vector<Front> fronts;
  double thickness = 0;
  bool defined = false;
};

//! 10-90 fronts of a binned solid fraction. `span` is the axis length.
//! A flat crystal, with no fall through 0.5, is undefined.
[[nodiscard]] Layers fitSolidFraction(const std::vector<double> &phi,
                                      const std::vector<double> &mid, double span);

//! CHILL+, the shell label and q6 binned along `axis` (0, 1 or 2).
//! A negative axis uses the longest box edge.
[[nodiscard]] Layers
layersOf(const molSys::PointCloud<molSys::Point<double>, double> &cloud, int axis = -1,
         double binWidth = 1.0);

struct SlopeInterval {
  double slope = 0;
  double lo = 0;
  double hi = 0;
  int n = 0;
};

//! Ordinary-least-squares slope, with a percentile bootstrap.
//! A draw whose x values have no spread is dropped.
[[nodiscard]] SlopeInterval bootstrapSlope(const std::vector<double> &x,
                                           const std::vector<double> &y, int nDraw,
                                           std::uint64_t seed);

struct CriticalSize {
  double nStar = 0;
  double sigma = 0; //! sample standard deviation of the window
  double lo = 0;
  double hi = 0;
  int onset = -1;
  int window = 0;
  bool upperBound = false;
  bool defined = false;
};

//! Critical size of one cluster-size series. See the file comment.
[[nodiscard]] CriticalSize criticalSize(const std::vector<double> &series, int nDraw,
                                        std::uint64_t seed);

struct ClusterPart {
  int size = 0;
  int cubic = 0;
  int hexagonal = 0;
  int interfacial = 0;
  int shellCubic = 0;
  int shellHex = 0;
  int shellStack = 0;
  int shellOther = 0;
};

struct NucleusFrame {
  int frame = 0;
  int nop = 0;
  int chillCubic = 0;
  int chillHex = 0;
  int chillInter = 0;
  int chillLiquid = 0;
  int shellCubic = 0;
  int shellHex = 0;
  int shellStack = 0;
  int shellOther = 0;
  int q6Ice = 0;
  double q6barMean = 0;
  //! Shell-stacking molecules, split by the published CHILL+ label.
  int stackAsCubic = 0;
  int stackAsHex = 0;
  int stackAsInter = 0;
  int stackAsLiquid = 0;
  //! Molecules that both the shell and CHILL+ call cubic or hexagonal,
  //! and the subset whose polytype agrees.
  int polytypeBoth = 0;
  int polytypeAgree = 0;
  ClusterPart core;
  ClusterPart extended;
  ClusterPart shell;
  ClusterPart q6;
};

//! One frame. The published CHILL+ label is left on the cloud.
[[nodiscard]] NucleusFrame
nucleusOf(molSys::PointCloud<molSys::Point<double>, double> &cloud, int frame,
          double q6Cut = kQ6BarIceCut);

struct AtomRow {
  int id = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  int chill = 2;
  int chillOvito = 2;
  int shell = -1;
  double q6 = 0;
  double q6bar = 0;
  int inExtended = 0;
};

//! Per-molecule labels. `inExtended` marks the largest CHILL+ cluster
//! of cubic, hexagonal and interfacial molecules.
[[nodiscard]] std::vector<AtomRow>
labelAtoms(molSys::PointCloud<molSys::Point<double>, double> cloud,
           double q6Cut = kQ6BarIceCut);

} // namespace science

#endif
