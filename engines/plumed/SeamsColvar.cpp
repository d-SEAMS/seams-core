// Copyright (c) d-SEAMS developers
// PLUMED collective variable. q6 carries finite-difference derivatives so a
// RESTRAINT can bias it. ice is the CHILL+ cubic fraction. ring is the
// fraction of atoms on a primitive six-ring.

#include "core/ActionAtomistic.h"
#include "core/ActionRegister.h"
#include "core/Colvar.h"
#include "seams.h"
#include "../triclinic_box.hpp"

#include <cmath>
#include <vector>

namespace PLMD {

class SeamsCV : public Colvar {
 public:
  static void registerKeywords(Keywords &keys);
  explicit SeamsCV(const ActionOptions &ao);
  void calculate() override;

 private:
  double cutoff_;
  bool score(const std::vector<double> &xyz, const double box9[9], double &q6, double &ice,
             double &ring) const;
};

void SeamsCV::registerKeywords(Keywords &keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms", "ATOMS", "atoms of the cell");
  keys.add("compulsory", "CUTOFF", "3.5", "neighbour cutoff in the current length unit");
  keys.addOutputComponent("q6", "default", "mean local q6, the continuous bias component");
  keys.addOutputComponent("ice", "default", "fraction of atoms labelled cubic by published CHILL+");
  keys.addOutputComponent("ring", "default", "fraction of atoms on a primitive six-ring");
}

SeamsCV::SeamsCV(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao), cutoff_(3.5) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);
  parse("CUTOFF", cutoff_);
  checkRead();
  requestAtoms(atoms);
  addComponentWithDerivatives("q6");
  componentIsNotPeriodic("q6");
  addComponentWithDerivatives("ice");
  componentIsNotPeriodic("ice");
  addComponentWithDerivatives("ring");
  componentIsNotPeriodic("ring");
}

bool SeamsCV::score(const std::vector<double> &xyz, const double box9[9], double &q6, double &ice,
                    double &ring) const {
  const int n = static_cast<int>(xyz.size() / 3);
  seams_system *sys = seams_system_create();
  if (sys == nullptr || seams_system_set_frame(sys, xyz.data(), n, box9, nullptr, nullptr) != SEAMS_OK) {
    seams_system_destroy(sys);
    return false;
  }
  std::vector<int> labels(static_cast<std::size_t>(n));
  int counts[9] = {};
  std::vector<double> q6v(static_cast<std::size_t>(n));
  if (seams_chill_plus(sys, cutoff_, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels.data(),
                       counts) != SEAMS_OK ||
      seams_q4q6(sys, cutoff_, 0, nullptr, q6v.data()) != SEAMS_OK) {
    seams_system_destroy(sys);
    return false;
  }
  double sum = 0.0;
  for (double value : q6v) sum += value;
  q6 = sum / static_cast<double>(n);
  ice = static_cast<double>(counts[SEAMS_LABEL_CUBIC]) / static_cast<double>(n);
  seams_rings *rings = nullptr;
  if (seams_rings_find(sys, cutoff_, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK) {
    seams_system_destroy(sys);
    return false;
  }
  if (seams_rings_fraction(rings, 6, &ring) != SEAMS_OK) {
    seams_rings_destroy(rings);
    seams_system_destroy(sys);
    return false;
  }
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  return true;
}

void SeamsCV::calculate() {
  // With UNITS LENGTH=A, getPosition and getBox are already Angstrom.
  // The driver still reads an xyz file as nm, so that deck is written in nm.
  const unsigned n = getNumberOfAtoms();
  std::vector<double> xyz(3 * n);
  for (unsigned i = 0; i < n; ++i) {
    const Vector p = getPosition(i);
    xyz[3 * i] = p[0];
    xyz[3 * i + 1] = p[1];
    xyz[3 * i + 2] = p[2];
  }
  const Tensor box = getBox();
  double raw[3][3];
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      raw[a][b] = box(a, b);
    }
  }
  double box9[9];
  plumed_box_to_lammps(raw, box9);
  double q6 = 0, ice = 0, ring = 0;
  if (!score(xyz, box9, q6, ice, ring)) {
    error("d-SEAMS rejected this frame");
  }
  // Central difference on q6 only. ice and ring stay finite but are not biased.
  const double eps = 1.0e-3;
  std::vector<Vector> dq(n);
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned d = 0; d < 3; ++d) {
      xyz[3 * i + d] += eps;
      double q_plus = 0, ice_plus = 0, ring_plus = 0;
      if (!score(xyz, box9, q_plus, ice_plus, ring_plus)) {
        error("d-SEAMS rejected a displaced frame");
      }
      xyz[3 * i + d] -= 2.0 * eps;
      double q_minus = 0, ice_minus = 0, ring_minus = 0;
      if (!score(xyz, box9, q_minus, ice_minus, ring_minus)) {
        error("d-SEAMS rejected a displaced frame");
      }
      xyz[3 * i + d] += eps;
      dq[i][d] = (q_plus - q_minus) / (2.0 * eps);
    }
  }
  Value *vq = getPntrToComponent("q6");
  Value *vi = getPntrToComponent("ice");
  Value *vr = getPntrToComponent("ring");
  vq->set(q6);
  vi->set(ice);
  vr->set(ring);
  const Vector zero;
  for (unsigned i = 0; i < n; ++i) {
    setAtomsDerivatives(vq, static_cast<int>(i), dq[i]);
    setAtomsDerivatives(vi, static_cast<int>(i), zero);
    setAtomsDerivatives(vr, static_cast<int>(i), zero);
  }
  setBoxDerivativesNoPbc(vq);
  setBoxDerivativesNoPbc(vi);
  setBoxDerivativesNoPbc(vr);
}

PLUMED_REGISTER_ACTION(SeamsCV, "SEAMS")

}  // namespace PLMD
