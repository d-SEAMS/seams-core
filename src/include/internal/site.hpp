//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#ifndef SEAMS_SITE_H_
#define SEAMS_SITE_H_

#include <mol_sys.hpp>

#include <array>
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
  std::vector<std::vector<int>> members;  ///< per ion, the shell molecules (cloud indices)
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

/// Rings of the water network that pass through a first shell: the census
/// by size of every ring with at least one vertex in `shell`. An ion is not
/// a vertex of the network, so the rings it would have closed are gone and
/// the rings its shell still carries measure how far the network survives
/// around it (the hydration-shell ring census of a brine or an ionic
/// solution). `census[s]` counts rings of size s up to `maxRingSize`.
std::vector<int> shellRingCensus(const std::vector<std::vector<int>> &rings,
                                 const std::vector<int> &shell, int maxRingSize);

/// Guests read against enumerated cages. A cage is a set of vertex atoms;
/// its centre is the periodic centroid of the vertices. A guest (methane,
/// THF, an ion, any atom outside the network) sits in the cage whose
/// centre is nearest to it when that centre is within `radius`; a guest
/// with no centre within `radius` is free. This is the occupancy of a
/// clathrate hydrate (fraction of 5^12 and 5^12 6^2 cages filled) and the
/// test of whether an ion sits inside the network rather than at a vertex.
struct GuestOccupancy {
  std::vector<int> guestsPerCage;        ///< per cage, guests inside
  std::vector<int> cageOfGuest;          ///< per guest in input order, cage index or -1
  std::vector<double> centreDistance;    ///< per guest, distance to its cage centre (-1 when free)
  int occupied = 0;                      ///< cages with at least one guest
  int multiply = 0;                      ///< cages with more than one guest
  int free = 0;                          ///< guests in no cage
};
/// `cages` are vertex index lists into `yCloud.pts`; `radius` in the
/// cloud's length unit (half the cage diameter, about 4 A for a 5^12 cage).
GuestOccupancy
guestOccupancy(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               const std::vector<std::vector<int>> &cages,
               const std::vector<int> &guestIndices, double radius);

/// Periodic centroid of a set of atoms: every atom is unwrapped to its
/// minimum image about the first, and the mean is taken there.
std::array<double, 3>
periodicCentroid(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                 const std::vector<int> &atoms);

/// Occupancy by ray-parity: a guest sits in a cage when a ray from the
/// guest, in the frame that unwraps the cage about its periodic
/// centroid, crosses an odd number of fan-triangulated faces. Among
/// cages that contain the guest the nearest centroid wins. `rings` is
/// the ring list the cage faces index; `cageFaces[c]` is the face
/// indices of cage c.
GuestOccupancy
guestOccupancyInside(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     const std::vector<std::vector<int>> &rings,
                     const std::vector<std::vector<int>> &cageFaces,
                     const std::vector<int> &guestIndices);

/// Connected components of ice atoms on the water graph, and the ion
/// count per component: each ion is assigned to the cluster of the
/// nearest ice atom within `cutoff`.
struct IceClusterIons {
  std::vector<int> clusterOf;      ///< per atom, cluster id or -1
  std::vector<int> ionsInCluster;  ///< per cluster
  std::vector<int> clusterOfIon;   ///< per ion in input order, or -1
  int nClusters = 0;
};

IceClusterIons
iceClusterIonCensus(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    const std::vector<bool> &ice,
                    const std::vector<std::vector<int>> &nListByIndex,
                    const std::vector<int> &ionIndices, double cutoff);
} // namespace site

#endif // SEAMS_SITE_H_
