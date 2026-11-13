// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// PLUMED collective variable. q6 carries finite-difference derivatives so a
// RESTRAINT can bias it. ice is the CHILL+ cubic fraction and ring is the
// fraction of atoms on a primitive six-ring. Both are discrete, so they are
// registered without derivatives and cannot be biased.

#include "plumed/core/ActionAtomistic.h"
#include "plumed/core/ActionRegister.h"
#include "plumed/core/Colvar.h"
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
  mutable int hexagonal_ = 0;
  mutable int six_ = 0;
  bool score(const std::vector<double> &xyz, const double box9[9], double &q6, double &ice,
             double &ring) const;
};

void SeamsCV::registerKeywords(Keywords &keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms", "ATOMS", "atoms of the cell");
  keys.add("compulsory", "CUTOFF", "3.5", "neighbour cutoff in angstrom");
  keys.addOutputComponent("q6", "default", "mean local q6, the continuous bias component");
  keys.addOutputComponent("ice", "default", "CHILL+ cubic fraction; not biasable, no derivatives");
  keys.addOutputComponent("ring", "default", "primitive six-ring fraction; not biasable, no derivatives");
  keys.addOutputComponent("hexagonal", "default", "CHILL+ hexagonal count; not biasable, no derivatives");
  keys.addOutputComponent("six", "default", "primitive six-ring count; not biasable, no derivatives");
}

SeamsCV::SeamsCV(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao), cutoff_(3.5) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);
  parse("CUTOFF", cutoff_);
  checkRead();
  requestAtoms(atoms);
  addComponentWithDerivatives("q6");
  componentIsNotPeriodic("q6");
  addComponent("ice");
  componentIsNotPeriodic("ice");
  addComponent("ring");
  componentIsNotPeriodic("ring");
  addComponent("hexagonal");
  componentIsNotPeriodic("hexagonal");
  addComponent("six");
  componentIsNotPeriodic("six");
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
  hexagonal_ = counts[SEAMS_LABEL_HEXAGONAL];
  seams_rings *rings = nullptr;
  if (seams_rings_find(sys, cutoff_, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK) {
    seams_system_destroy(sys);
    return false;
  }
  int hist[7] = {};
  if (seams_rings_fraction(rings, 6, &ring) != SEAMS_OK ||
      seams_rings_histogram(rings, hist, 7) != SEAMS_OK) {
    seams_rings_destroy(rings);
    seams_system_destroy(sys);
    return false;
  }
  six_ = hist[6];
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  return true;
}

void SeamsCV::calculate() {
  // getPosition and getBox are PLUMED's internal nanometres. d-SEAMS
  // expects angstroms, and CUTOFF is parsed as a bare number in angstroms.
  const double to_angstrom = 10.0;
  const unsigned n = getNumberOfAtoms();
  std::vector<double> xyz(3 * static_cast<std::size_t>(n));
  for (unsigned i = 0; i < n; ++i) {
    const Vector p = getPosition(i);
    xyz[3 * i] = p[0] * to_angstrom;
    xyz[3 * i + 1] = p[1] * to_angstrom;
    xyz[3 * i + 2] = p[2] * to_angstrom;
  }
  const Tensor box = getBox();
  double raw[3][3];
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      raw[a][b] = box(a, b) * to_angstrom;
    }
  }
  double box9[9];
  plumed_box_to_lammps(raw, box9);
  double q6 = 0, ice = 0, ring = 0;
  if (!score(xyz, box9, q6, ice, ring)) {
    error("d-SEAMS rejected this frame");
  }
  const int hexagonal = hexagonal_;
  const int six = six_;
  Value *vq = getPntrToComponent("q6");
  Value *vi = getPntrToComponent("ice");
  Value *vr = getPntrToComponent("ring");
  Value *vh = getPntrToComponent("hexagonal");
  Value *vs = getPntrToComponent("six");
  vq->set(q6);
  vi->set(ice);
  vr->set(ring);
  vh->set(static_cast<double>(hexagonal));
  vs->set(static_cast<double>(six));
  // A PRINT of the value does not allocate derivative storage. Writing
  // into that empty array aborts. A RESTRAINT does allocate it.
  if (vq->getNumberOfDerivatives() == 0) {
    return;
  }
  // Central difference on q6 only. ice and ring have no derivatives.
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
      // eps is an angstrom step. PLUMED differentiates with respect to its
      // internal nanometres, and one nanometre is ten angstroms.
      dq[i][d] = (q_plus - q_minus) / (2.0 * eps) * to_angstrom;
    }
  }
  for (unsigned i = 0; i < n; ++i) {
    setAtomsDerivatives(vq, static_cast<int>(i), dq[i]);
  }
  setBoxDerivativesNoPbc(vq);
}

PLUMED_REGISTER_ACTION(SeamsCV, "SEAMS")

}  // namespace PLMD
