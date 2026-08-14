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

#ifndef SEAMS_SHAPEMATCH_H_
#define SEAMS_SHAPEMATCH_H_

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

#include <absOrientation.hpp>
#include <mol_sys.hpp>
#include <pntCorrespondence.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

namespace match {

//! Shape-matching for a pair of polygon basal rings. Returns true if the pair
//! of basal rings form a prism block.
bool matchPrism(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                const std::vector<std::vector<int>> &nList,
                const Eigen::MatrixXd &refPoints, std::vector<int> &basal1,
                std::vector<int> &basal2, std::vector<double> &rmsdPerAtom,
                bool isPerfect = true);

//! Shape-matching for a pair of polygon basal rings. Returns true if the pair
//! of basal rings form a prism block.
bool matchUntetheredPrism(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, const Eigen::MatrixXd &refPoints,
    std::vector<int> &basal1, std::vector<int> &basal2,
    std::vector<double> &rmsdPerAtom);

//! Shape-matching for a pair of polygon basal rings, comparing with a complete
//! prism block. Returns true if the pair of basal rings form a prism block.
bool matchPrismBlock(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     const std::vector<std::vector<int>> &nList,
                     const Eigen::MatrixXd &refPoints, std::vector<int> &basal1,
                     std::vector<int> &basal2, int &beginIndex);

//! Update the per-particle RMSD for a prism block basal ring.
[[nodiscard]] int updatePerAtomRMSDRing(const std::vector<int> &basalRing, int startingIndex,
                          const std::vector<double> &rmsdFromMatch,
                          std::vector<double> &rmsdPerAtom);

//! Update the RMSD of each particle in a prism block basal ring with the RMSD
//! of the ring.
[[nodiscard]] int updateRMSDRing(const std::vector<int> &basalRing, int startingIndex,
                   double rmsdVal, std::vector<double> &rmsdPerAtom);

} // namespace match

#endif // SEAMS_SHAPEMATCH_H_
