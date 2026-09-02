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
 *  Table::of does not write polar or apolar back onto the table;
 *  indicesOf applies those two derived unions.
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
  waterIce,      // default; CHILL/TUM allowed
  ionicLiquid,   // CHILL/TUM refused
  moltenSalt,    // CHILL/TUM refused
  des,           // CHILL/TUM refused
  electrolyte,   // CHILL/TUM only if the caller also names a waterIce subset
  confinedIL,    // CHILL/TUM refused
  confinedWater, // 2D RDF / monolayer rings; bulk CHILL refused
  networkFormer  // silica / BeF2; Franzblau yes, CHILL no
};

struct Table {
  Family family = Family::waterIce;
  std::unordered_map<int, Kind> typeToKind;   // LAMMPS type ID
  std::unordered_map<int, Kind> atomOverride; // atom ID wins over type
  Kind of(const molSys::Point<double> &p) const;
  Kind ofType(int typeId) const;
};

std::vector<int>
indicesOf(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
          const Table &table, Kind kind);

int lammpsTypeOfKind(const Table &table, Kind kind); // error if not unique

// One vertex per molID that carries ionKind (cationHead or anion).
// Coordinates: unweighted geometric COM of atoms of that molID whose
// kind is ionKind, unwrapped with gen::relDist to the first atom of
// the group (createMolIDAtomIDMultiMap, mol_sys.hpp:187-190). Copies
// box / boxLow from src. Point::type is 1 for cationHead molecules
// and 2 for anion molecules so neighList(merged, 1, 2) is legal.
// A molecule that is already one tagged site (UA, or one designated
// type per ion) is a no-op: the COM is that site.
molSys::PointCloud<molSys::Point<double>, double>
ionCloud(const molSys::PointCloud<molSys::Point<double>, double> &src,
         const Table &table);

Table parseSiteSpec(std::string_view spec);


/// Ions read against a per-atom ice assignment. Ions are not part of the
/// hydrogen-bond network; the assignment is computed on the water and each
/// ion is classed by its first water shell: every shell molecule labelled
/// is `ice`, none is `liquid`, otherwise `front`. An ion with an empty shell
/// is liquid.
enum class IonState { liquid = 0, front = 1, ice = 2 };
struct IonEnvironment {
  std::vector<int> ion;            ///< cloud indices of the ions, in input order
  std::vector<int> shell;          ///< water molecules within the cutoff
  std::vector<double> iceFraction; ///< labelled share of that shell
  std::vector<IonState> state;
  int nIce = 0;
  int nFront = 0;
  int nLiquid = 0;
};
/// `iceFlag` is indexed like `yCloud.pts`; `waterType` selects the water
/// oxygens (0 accepts every atom that is not an ion); `cutoff` is the first
/// shell radius in the cloud's length unit.
IonEnvironment
ionEnvironment(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               const std::vector<bool> &iceFlag, const std::vector<int> &ionIndices,
               int waterType, double cutoff);
} // namespace site

#endif // SEAMS_SITE_H_
