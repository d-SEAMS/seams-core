//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#ifndef SEAMS_ORDER_PARAMETER_H_
#define SEAMS_ORDER_PARAMETER_H_

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <cage.hpp>
#include <mol_sys.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

namespace topoparam {

//! Calculates the height%, an average measure of filled volume. The average
//! height of a prism can be taken to be 2.75-2.85 Angstrom. (Koga et. al.,
//! 2001)
double
normHeightPercent(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                  int nPrisms, double avgPrismHeight);

//! Calculates the coverage area%, an area-based measure of relative proportion
//! of monolayer ices.
std::vector<double>
calcCoverageArea(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                 const std::vector<std::vector<int>> &rings, double sheetArea);

//! Calculates the projected area on the XY, YZ and XZ planes
std::vector<double>
projAreaSingleRing(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                   const std::vector<int> &ring);

//! Per-atom Rodger F4 = mean of cos(3 phi) over H-O...O-H dihedrals of
//! this oxygen with its oxygen neighbours. Hydrogens share molID with
//! their oxygen. An oxygen with no hydrogens (mW) is quiet_NaN.
//! nList is by atom ID with the leading self entry, as neighListO.
std::vector<double>
rodgerF4(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
         const std::vector<std::vector<int>> &nList, int oxygenType,
         int hydrogenType);

//! Mean of the finite per-atom F4 values; quiet_NaN when none are finite.
double meanFinite(const std::vector<double> &values);

//! Jump-rotor tau_90 from two frames. Each water is the H-H vector of
//! the two hydrogens that share the oxygen molID. Returns dt when any
//! molecule rotates by at least 90 degrees, quiet_NaN otherwise.
double jumpRotorTau90(
    const molSys::PointCloud<molSys::Point<double>, double> &frame0,
    const molSys::PointCloud<molSys::Point<double>, double> &frame1, double dt,
    int oxygenType, int hydrogenType);

//! CHILL+ cubic/hexagonal molecules binned into basal layers along
//! `axis` (0=x, 1=y, 2=z). phiC is N_Ic / (N_Ic + N_Ih). sequence is
//! one character per layer: C cubic-majority, H hexagonal-majority,
//! M a tie with ice present, . empty.
struct LayerStack {
  double phiC = 0.0;
  std::string sequence;
  std::vector<int> cubicPerLayer;
  std::vector<int> hexPerLayer;
};

LayerStack
layerCubicity(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
              int axis = 2, double layerWidth = 3.7);

//! Literature I_sd reference: CHILL+ molecules binned by coordinate.
//! TUM/rings stacking: only HC-basal and DDC-equatorial six-rings vote.
//! Prismatic/peripheral rings and non-plane cages (clathrate 5^12) do
//! not write a letter. Ions stay off the ring graph the caller passes.
//! phiC is N_equatorial / (N_basal + N_equatorial) over plane rings.
LayerStack
tumLayerStack(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
              const std::vector<std::vector<int>> &rings,
              const std::vector<bool> &basal, const std::vector<bool> &equatorial,
              int axis = 2, double layerWidth = 3.7);

} // namespace topoparam

#endif // SEAMS_ORDER_PARAMETER_H_
