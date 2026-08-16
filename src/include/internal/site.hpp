//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#ifndef SEAMS_SITE_H_
#define SEAMS_SITE_H_

#include <mol_sys.hpp>

#include <string_view>
#include <unordered_map>
#include <vector>

/** @file site.hpp
 *  @brief Per-atom site chemistry, separate from CHILL ice labels.
 *
 *  LAMMPS type IDs are not chemistry. A Table maps type (and optional
 *  atom-ID override) to Kind. Family is an input, not an inference.
 */

namespace site {

enum class Kind {
  unspecified,
  cationHead,
  anion,
  tail,
  donorH,
  acceptor,
  polar,
  apolar,
  waterO,
  waterH,
  solvent
};

enum class Family {
  waterIce,
  ionicLiquid,
  moltenSalt,
  des,
  electrolyte,
  confinedIL,
  confinedWater,
  networkFormer
};

struct Table {
  Family family = Family::waterIce;
  std::unordered_map<int, Kind> typeToKind;
  std::unordered_map<int, Kind> atomOverride;

  Kind of(const molSys::Point<double> &p) const;
  Kind ofType(int typeId) const;
};

std::vector<int>
indicesOf(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
          const Table &table, Kind kind);

int lammpsTypeOfKind(const Table &table, Kind kind);

molSys::PointCloud<molSys::Point<double>, double>
ionCloud(const molSys::PointCloud<molSys::Point<double>, double> &src,
         const Table &table);

Table parseSiteSpec(std::string_view spec);

} // namespace site

#endif // SEAMS_SITE_H_
