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

#ifndef SEAMS_NUCLEATION_PATHWAY_HPP_
#define SEAMS_NUCLEATION_PATHWAY_HPP_

#include <mol_sys.hpp>

#include <vector>

/** @file nucleation_pathway.hpp
 *  @brief One-frame ice census and the trajectory reductions used to
 *  read a nucleation pathway.
 *
 *  CHILL+ is applied twice, with the published Nguyen and Molinero
 *  edges and with the ovito-chill+ rule. Neither call changes those
 *  edges. Rings are the Franzblau primitive rings of the cutoff graph.
 *  q6-bar is the Lechner-Dellago average from chill::steinhardtQlCutoff
 *  on that same cutoff. The fused q4/q6/CHILL+ entry point on later
 *  neighbour work returns the local parameters and does not fill
 *  q-bar, so this census keeps the averaged call.
 *  Symmetry labels are the diamond common-neighbour test: 1 is cubic
 *  diamond (ice Ic) and 0 is hexagonal diamond (ice Ih).
 *
 *  The largest cluster of a mask is the largest connected component of
 *  the cutoff graph restricted to that mask. A tie takes the component
 *  whose root index is smaller.
 */

namespace nucleation {

//! Published CHILL+ cutoff. Not a fitted parameter.
inline constexpr double kCutoffA = 3.5;

//! Primitive rings through this size. Depth 7 is the adjacency-row search.
inline constexpr int kRingDepth = 7;

//! Lechner-Dellago \f$\bar{q}_6\f$ above this value is an ice-like
//! environment. Set from the gap between a perfect mW Ic crystal and a
//! liquid reference, not from the nucleation trajectory.
inline constexpr double kQ6BarIceCut = 0.36;

struct FrameRecord {
  int frame = 0;
  long long timestep = 0;
  int nop = 0;
  //! Primitive rings of each size. Index is the ring size. 0 and 1 stay 0.
  long long ringCount[8] = {};
  int chillCubic = 0;
  int chillHex = 0;
  int chillInter = 0;
  int chillClath = 0;
  int chillInterClath = 0;
  int chillWater = 0;
  int ovitoCubic = 0;
  int ovitoHex = 0;
  int ovitoInter = 0;
  int ovitoClath = 0;
  int ovitoInterClath = 0;
  int ovitoWater = 0;
  int symCubic = 0;
  int symHex = 0;
  double q6barMean = 0;
  double q6barP50 = 0;
  int q6barIce = 0;
  //! Largest cluster. c and h are cubic and hexagonal members.
  //! For the q6-bar cluster, c and h are diamond-symmetry labels,
  //! because q6-bar does not separate Ic from Ih.
  int nChillCore = 0;
  int cChillCore = 0;
  int hChillCore = 0;
  int nChillExt = 0;
  int cChillExt = 0;
  int hChillExt = 0;
  int nOvitoCore = 0;
  int cOvitoCore = 0;
  int hOvitoCore = 0;
  int nOvitoExt = 0;
  int cOvitoExt = 0;
  int hOvitoExt = 0;
  int nQ6 = 0;
  int cQ6 = 0;
  int hQ6 = 0;
  int nSym = 0;
  int cSym = 0;
  int hSym = 0;
};

//! Label one frame. Leaves iceType as the published CHILL+ label.
//! typeI selects the molecules passed to q6-bar. q6Cut is the ice cut
//! on \f$\bar{q}_6\f$.
[[nodiscard]] FrameRecord
analyzeFrame(molSys::PointCloud<molSys::Point<double>, double> &cloud,
             double cutoff, int ringDepth, int typeI, double q6Cut, int frame);

//! Cubic / (cubic + hexagonal). Quiet NaN when the denominator is 0.
[[nodiscard]] double stackFraction(int cubic, int hexagonal);

//! Index of the last sample before the strict climb that reaches the
//! maximum. A flat series returns 0 and climbed is false.
struct Ascent {
  int index = 0;
  int maxIndex = 0;
  double baseline = 0;
  double maximum = 0;
  bool climbed = false;
};

[[nodiscard]] Ascent finalAscent(const std::vector<double> &series);

//! First sample at which the series has covered half of its rise from
//! the first sample to the maximum. A series that never rises has
//! ok false. This is the operational marker used for induction on one
//! noisy growth window. It is not a committor.
struct HalfRise {
  bool ok = false;
  int index = -1;
  double start = 0;
  double level = 0;
  double maximum = 0;
};

[[nodiscard]] HalfRise halfRise(const std::vector<double> &series);

//! Positive lag means the ring series moves before the ice series.
//! halfFrames is the ice half-rise index minus the ring half-rise index.
struct Lead {
  bool ok = false;
  int ascentFrames = 0;
  int halfFrames = 0;
  int corrLag = 0;
  double corr = 0;
};

//! Ring signal is n6 / (n5 + n6). Ice signal is the ice-like count
//! divided by the molecule count. Both are intensive.
[[nodiscard]] Lead ringLead(const std::vector<double> &ringFraction,
                            const std::vector<double> &iceFraction);

struct Shift {
  int dSize = 0;
  int dFrame = 0;
};

//! How far definition b moves the onset size and the onset frame
//! relative to definition a. dSize = b.size - a.size at each onset.
[[nodiscard]] Shift definitionShift(int sizeA, int frameA, int sizeB, int frameB);

} // namespace nucleation

#endif
