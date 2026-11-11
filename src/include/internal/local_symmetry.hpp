#ifndef SEAMS_LOCAL_SYMMETRY_H_
#define SEAMS_LOCAL_SYMMETRY_H_

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

// In-tree matcher (quaternion Kabsch, fixed permutations). backend "kernel".
Overlay matchKernel(const Eigen::MatrixXd &environment);

// Per-atom ice label for a periodic cell. kind[i] is iceIh, iceIc, iceIsd,
// or -1 when the oxygen is not a bulk ice shell. Labels match matchKernel
// on the shells the test cells produce. Parallel over atoms.
void classifyIce(const Eigen::MatrixXd &cartesian, const Eigen::Matrix3d &vectors,
                 std::vector<int> &kind);

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
