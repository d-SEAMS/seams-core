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

#ifndef SEAMS_CHILL_PLUS_SCHEME_H_
#define SEAMS_CHILL_PLUS_SCHEME_H_

#include <mol_sys.hpp>

#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

/** @file chill_plus_scheme.hpp
 *  @brief The two CHILL+ rules.
 *
 *  `chill+` is the published rule. Table 1 footnote a of Nguyen and
 *  Molinero (J. Phys. Chem. B 2015, 119, 9369) counts neighbours inside
 *  the cutoff shell, and a crystal row requires that count to be 4.
 *  Equation 1 averages q_3 over the four nearest of those neighbours.
 *  The interfacial rows keep their neighbour test.
 *
 *  `ovito-chill+` is the OVITO 3.17.2 Chill+ modifier: the same shell
 *  count, q_3 averaged over the whole cutoff shell, open c_ij edges,
 *  and no neighbour test.
 *
 *  Plain CHILL does not use this file.
 */

namespace chill {

enum class ChillPlusRule { Published, Ovito };

//! Neighbours in one cutoff row, ignoring the leading atom id and any
//! id that is not in the cloud. This is the Table 1 shell count.
[[nodiscard]] int cutoffShellCount(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, int iatom);

//! `chill+`, `ovito-chill+`, a YAML document with `rule:` or
//! `chill-plus:`, or a Lua table fragment `rule = "..."`.
//! An empty string is the published rule. Anything else throws.
[[nodiscard]] ChillPlusRule chillPlusRule(std::string_view spec);

//! Read `spec` as a rule name. If `spec` names a readable file, the
//! file is scanned for a `rule` or `chill-plus` key. Comments are
//! skipped. The file is not executed as Lua or YAML.
[[nodiscard]] ChillPlusRule chillPlusRuleFromFile(const std::string &spec);

[[nodiscard]] const char *chillPlusRuleName(ChillPlusRule rule);

//! Published or OVITO CHILL+ on one cloud. The neighbour list is the
//! cutoff shell the rule counts.
void applyChillPlus(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, ChillPlusRule rule,
    bool isSlice = false);

} // namespace chill

#endif
