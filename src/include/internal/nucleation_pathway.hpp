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
 *  on that same cutoff, with the average filled. chill::q4q6ChillPlus
 *  at 29b3c2b2 returns that average as well, with q4, q6 and the
 *  published CHILL+ labels. The call is built on a passFast rewrite
 *  that conflicts with the shared orthorhombic stencil, so this census
 *  keeps steinhardtQlCutoff.
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

//! Grow basin of one label. The barrier is the cluster size at the foot
//! of that label's own final strict climb to its maximum. A series that
//! never climbs has occupied false: it never enters the basin. The
//! barrier is not a shared multiple of the molecule count.
struct GrowBasin {
  bool occupied = false;
  double barrier = 0;
  int index = -1;
};

[[nodiscard]] GrowBasin growBasin(const std::vector<double> &series);

//! Outcome of one replica. committed means the final size is still at or
//! above that replica's own grow-basin barrier. truncated means the series
//! ends on the cap 2N while that cap is also its maximum: the run was cut
//! there and is not a committor sample. The barrier is never 2N itself.
enum class CommittorFate { uncommitted, committed, truncated };

struct CommittorCall {
  CommittorFate fate = CommittorFate::uncommitted;
  GrowBasin basin;
  bool cutAtTwiceN = false;
};

[[nodiscard]] CommittorCall committorOf(const std::vector<double> &series, int molecules);

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
