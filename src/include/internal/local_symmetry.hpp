#ifndef SEAMS_LOCAL_SYMMETRY_H_
#define SEAMS_LOCAL_SYMMETRY_H_

#include <mol_sys.hpp>

#include <Eigen/Core>
#include <string>
#include <vector>

// Local symmetry of ice, hydrate cages and rings.
//
// Matching uses IRA (Gunde, Salles, Hemeryck, Martin-Samos, JCIM 2021,
// 10.1021/acs.jcim.1c00567) when libira is linked, and a Kabsch-Hungarian
// search otherwise. Point groups use SOFI (Gunde et al., JCP 2024,
// 10.1063/5.0215689). Space groups use spglib (Togo, Shinohara, Tanaka,
// STAM Methods 2024, 10.1080/27660400.2024.2384822) when libsymspg is linked.
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT

namespace symm {

enum class Kind {
  iceIh = 0,
  iceIc = 1,
  iceIsd = 2,
  cage512 = 3,
  cage51262 = 4,
  cage51264 = 5,
  ring5 = 6,
  ring6 = 7
};

inline constexpr int kTemplateCount = 8;

const char *name(Kind kind);

// Ideal oxygen coordinates, rows are atoms, angstrom. The first row of an
// ice shell is the central oxygen. Cage and ring templates are centred.
Eigen::MatrixXd points(Kind kind);

struct Overlay {
  int status = 1;
  Kind kind = Kind::iceIh;
  std::string name;
  double rmsd = 0;
  double hausdorff = 0;
  // R and t take the environment onto the template:
  // template[i] ≈ R * environment[permutation[i]] + t.
  Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
  Eigen::Vector3d translation = Eigen::Vector3d::Zero();
  std::vector<int> permutation;
  std::string backend;
};

// One row per template with the same atom count, lowest RMSD first.
// backend is "ira", "kabsch", or nullptr for whichever is linked.
std::vector<Overlay> score(const Eigen::MatrixXd &environment,
                           const char *backend = nullptr);
Overlay best(const Eigen::MatrixXd &environment, const char *backend = nullptr);

// In-tree matcher (3x3 SVD Kabsch, bond-graph permutations). backend "kernel".
Overlay matchKernel(const Eigen::MatrixXd &environment);

// Same matcher on many environments. Independent shells run in parallel.
void matchAll(const std::vector<Eigen::MatrixXd> &environments, std::vector<Overlay> &out);

// A bond is counted eclipsed when its azimuthal cosine exceeds this value.
// Both bonds are projected onto the plane perpendicular to the central bond.
// Ideal eclipsed is cos 0 = 1 and ideal staggered is cos 60 degrees = 0.5.
// The ideal-cell fast path does not use this cosine. It uses the raw dot
// product of the second-shell bond with another first-shell bond, which is
// 1/3 staggered and 7/9 eclipsed, and it compares the square of that number.
// At 220-270 K the staggered azimuthal cosines of mW ice reach about 0.85
// while the eclipsed bonds stay at 1, so 0.90 is the gap the second-shell
// template leaves between them. 0.75 counts those warmed staggered bonds
// as eclipsed.
inline constexpr double kEclipsedCosine = 0.90;

// Scale-free ceiling on the second-shell residual. Both clouds are centred,
// then the RMSD is taken after the optimal rotation and the optimal scale
// and divided by the template's root-mean-square radius (Horn 1987; Larsen,
// Schmidt, Schiotz 2016). Zero is a perfect shape. 0.12 is the quenched mW
// operating point. A TIP4P/2005 shell at the melting point still vibrates:
// raising the ceiling from 0.12 to 0.15 recovers quench-confirmed crystal,
// and the next step to 0.18 mostly adds atoms that neither PTM nor CHILL+
// call ice after the same minimisation.
inline constexpr double kShellRmsd = 0.15;

// Per-atom ice label. kind[i] is iceIh, iceIc, iceIsd, or -1. shellRmsd is
// the scale-free residual ceiling. The label is the template with the lowest
// residual after the bond-graph permutation of the tetrahedron.
void classifyIce(const Eigen::MatrixXd &cartesian, const Eigen::Matrix3d &vectors,
                 std::vector<int> &kind, double shellRmsd = kShellRmsd);

// Column SoA coordinates. The canon shell is built once. kind is written in full.
void classifyIceArrays(const double *x, const double *y, const double *z, int n, const double lat[3][3],
                       int *kind, double shellRmsd = kShellRmsd);

// Same labels as classifyIceArrays, packed from a point cloud. On ideal
// cubic ice the codes match diamondFromCsr (0 hexagonal, 1 cubic, -1 other).
void classifyIceCloud(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                      std::vector<int> &kind);

// Azimuthal cosine of each tetrahedral bond, largest first. The value is
// the best alignment of an outer neighbour with one of the other three
// bonds. Ideal eclipsed is 1 and ideal staggered is 0.5. Returns how many
// of the four bonds were measured, or -1 when the environment is not a shell.
int bondCosines(const Eigen::MatrixXd &environment, double out[4]);

// Bond-graph RMSD of a 17-atom shell against one ice template. Negative
// when the shell cannot be aligned. This does not consult the eclipsed count,
// so a hexagonal shell can be scored against the cubic template.
double iceTemplateRmsd(const Eigen::MatrixXd &environment, Kind kind);

// Alvarez continuous shape measure of environment against its best template.
// 0 is a perfect match. Scale, rotation and translation are free.
double shapeMeasure(const Eigen::MatrixXd &environment);

struct Group {
  int status = 1;
  std::string symbol;
  int nOperations = 0;
};

// status 2 when SOFI is not linked.
Group pointGroup(const Eigen::MatrixXd &xyz, double symThr = 0.05);

struct Space {
  int status = 1;
  int number = 0;
  std::string symbol;
};

// `vectors` has the three cell edges as rows. `cartesian` is n x 3.
// status 2 when spglib is not linked.
Space spaceGroup(const Eigen::Matrix3d &vectors, const Eigen::MatrixXd &cartesian,
                 double symprec = 1e-3);

// Space group of one conventional cell when every label is the same ideal
// ice and the box is a cubic integer supercell of that cell. Returns false
// when the frame is not that crystal; the caller then uses spaceGroup on
// the whole cell. The international symbol matches the full-cell call.
bool conventionalCrystalSpace(const Eigen::Matrix3d &vectors,
                              const Eigen::MatrixXd &cartesian,
                              const std::vector<int> &kind, double symprec,
                              Space &out);

enum class Crystal { iceIh, iceIc, fault };

struct Lattice {
  Eigen::Matrix3d vectors = Eigen::Matrix3d::Identity();
  Eigen::MatrixXd cartesian;
};

// Ideal ice Ih (P6_3/mmc), ice Ic (Fd-3m), and an Ih cell with one
// bilayer shifted onto the cubic stacking site.
Lattice crystal(Crystal which);

// Central oxygen plus its 16 nearest periodic neighbours, when the 16th
// neighbour lies inside 5 A. Empty when the atom is not in a bulk ice shell.
Eigen::MatrixXd iceShell(const Eigen::MatrixXd &cartesian,
                         const Eigen::Matrix3d &vectors, int atom);

bool iraLinked();
bool spglibLinked();
// "ira" when libira is linked, otherwise "kabsch".
const char *matcher();

} // namespace symm

#endif
