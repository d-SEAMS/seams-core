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

#ifndef SEAMS_BRINE_FRONT_H_
#define SEAMS_BRINE_FRONT_H_

#include <mol_sys.hpp>
#include <site.hpp>

#include <vector>

/** @file brine_front.hpp
 *  @brief In-situ measurements at a planar ice/brine front.
 *
 *  The ice label is an input. CHILL+ (or any other classifier) stays
 *  outside these functions, so a synthetic slab tests the arithmetic
 *  without a trajectory. Lengths are the cloud's length unit.
 */

namespace site {

inline constexpr double kWaterMolarMassKg = 0.01801528;

struct ConcentrationBin {
  double zLo = 0;
  double zHi = 0;
  double zMid = 0;
  int nWater = 0;
  int nIceWater = 0;
  int nIon = 0;
  double iceFraction = 0;
  double molality = 0; ///< formula units per kg of water in the bin
};

struct ConcentrationProfile {
  std::vector<ConcentrationBin> bins;
  double binWidth = 0;
  double zOrigin = 0;
  double length = 0;
};

/// Two planar fronts of one ice slab in a periodic box. `zEnter` is where
/// the ice fraction crosses `cut` on entering the ice along +z, and `zExit`
/// is the crossing on leaving it. `thickness` is the ice run along z.
struct FrontPair {
  double zEnter = 0;
  double zExit = 0;
  double thickness = 0;
  bool wrapped = false;
  bool found = false;
};

/// k = c_ice / c_liquid, both molalities of formula units. `defined` is
/// false when the liquid bulk contains no salt and the ice bulk does,
/// because the ratio is then unbounded.
struct RejectionCoefficient {
  double cIce = 0;
  double cLiquid = 0;
  double k = 0;
  int nIceWater = 0;
  int nLiquidWater = 0;
  int nIceIon = 0;
  int nLiquidIon = 0;
  bool defined = false;
};

struct LinearFit {
  double slope = 0;
  double intercept = 0;
  double slopeStdErr = 0; ///< standard error of the slope (not `stderr`, a macro on musl)
  int n = 0;
};

struct BrinePocket {
  std::vector<int> waters;
  int nIons = 0;
};

struct PocketCensus {
  std::vector<BrinePocket> pockets; ///< not including the liquid reservoir
  std::vector<int> sizes;           ///< pocket water counts, ascending
  int reservoirWaters = 0;
  int nPockets = 0;
};

/// Waters within `halfWidth` of either front, their primitive-ring census,
/// and the coordination defects (a water whose oxygen neighbour count is
/// not 4). `ringCensus[s]` is the number of primitive rings of size s that
/// lie entirely in the layer. `coordination[c]` counts waters with c
/// neighbours, and the last bin also holds every larger count.
struct FrontLayerReport {
  std::vector<int> waters;
  std::vector<int> ringCensus;
  std::vector<int> coordination;
  int nDefect = 0;
  int nFour = 0;
  int nHexagonal = 0;
  int nCubic = 0;
  int nInterfacial = 0;
  int nClathrate = 0;
  int nLiquidLabel = 0;
};

struct TrapEvent {
  int ion = -1; ///< index into the ion list, not the cloud
  IonState from = IonState::liquid;
  IonState to = IonState::liquid;
};

struct FrontOptions {
  int nBins = 20;
  int axis = 2;
  double frontCut = 0.5;
  double bulkIce = 0.85;
  double bulkLiquid = 0.15;
  int ionsPerFormula = 2;
  double shellCutoff = 3.5;
  double layerHalfWidth = 3.5;
  int maxRing = 6;
  int coordHist = 8;
};

struct FrontSnapshot {
  ConcentrationProfile profile;
  FrontPair fronts;
  RejectionCoefficient rejection;
  IonEnvironment ions;
  PocketCensus pockets;
  FrontLayerReport layer;
};

ConcentrationProfile concentrationProfile(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const std::vector<bool> &iceFlag, const std::vector<int> &waterIndices,
    const std::vector<int> &ionIndices, int nBins, int axis = 2,
    int ionsPerFormula = 2);

FrontPair locateFronts(const ConcentrationProfile &profile, double cut = 0.5);

RejectionCoefficient rejectionCoefficient(const ConcentrationProfile &profile,
                                          double bulkIce = 0.85,
                                          double bulkLiquid = 0.15);

LinearFit linearFit(const std::vector<double> &time,
                    const std::vector<double> &value);

PocketCensus brinePockets(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const std::vector<bool> &iceFlag, const std::vector<int> &waterIndices,
    const std::vector<int> &ionIndices, double cutoff);

FrontLayerReport frontLayer(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const std::vector<int> &waterIndices, const FrontPair &fronts, double halfWidth,
    double cutoff, int maxRing, int axis = 2);

/// A trapping event is a change of IonState into ice. A release is a
/// change out of ice. Other changes are recorded as well; the caller
/// filters on `to`.
std::vector<TrapEvent> trappingTransitions(const IonEnvironment &before,
                                           const IonEnvironment &after);

FrontSnapshot analyseFront(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const std::vector<bool> &iceFlag, const std::vector<int> &waterIndices,
    const std::vector<int> &ionIndices, const FrontOptions &options = {});

} // namespace site

#endif
