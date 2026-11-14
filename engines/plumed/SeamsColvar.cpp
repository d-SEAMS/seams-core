// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// PLUMED collective variable. q6 carries finite-difference derivatives so a
// RESTRAINT can bias it. ice is the CHILL+ cubic fraction and ring is the
// fraction of atoms on a primitive six-ring. Both are discrete, so they are
// registered without derivatives and cannot be biased.
//
// MODE=all (the default) computes every component. MODE=q6 computes only the
// mean local q6 and MODE=q6bar only the mean Lechner-Dellago q6-bar, so a
// timing against PLUMED's own Q6 or LOCAL_AVERAGE pays for that order
// parameter alone. DUMP=file appends one LAMMPS dump frame per call with
// the coordinates handed to d-SEAMS, in angstrom, and the per-atom values.

#include "plumed/config/version.h"
#include "plumed/core/ActionAtomistic.h"
#include "plumed/core/ActionRegister.h"
#include "plumed/core/Colvar.h"
#include "plumed/core/PlumedMain.h"
#include "plumed/tools/Units.h"
#if PLUMED_VERSION_MAJOR == 2 && PLUMED_VERSION_MINOR < 10
#include "plumed/core/Atoms.h"
#endif
#include "seams.h"
#include "../triclinic_box.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace PLMD {

class ActionWithValue;
void addPlainComponent(ActionWithValue *self, const std::string &name);
void addDerivativeComponent(ActionWithValue *self, const std::string &name);

class SeamsCV : public Colvar {
 public:
  static void registerKeywords(Keywords &keys);
  explicit SeamsCV(const ActionOptions &ao);
  ~SeamsCV() override;
  void calculate() override;

 private:
  enum class Mode { all, q6, q6bar };
  struct Score {
    double q6 = 0.0;
    double ice = 0.0;
    double ring = 0.0;
    int hexagonal = 0;
    int six = 0;
  };
  double cutoff_;
  Mode mode_ = Mode::all;
  std::string dump_;
  long calls_ = 0;
  double seconds_ = 0.0;
  std::vector<int> labels_;
  std::vector<double> peratom_;
  bool score(const std::vector<double> &xyz, const double box9[9], Score &out, bool keep);
  double lengthToAngstrom();
  void writeDump(const std::vector<double> &xyz, const double box9[9]);
};

void SeamsCV::registerKeywords(Keywords &keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms", "ATOMS", "atoms of the cell");
  keys.add("compulsory", "CUTOFF", "3.5", "neighbour cutoff in angstrom");
  keys.add("compulsory", "MODE", "all", "all, q6 (mean local q6 only) or q6bar (mean q6-bar only)");
  keys.add("optional", "DUMP", "append a LAMMPS dump frame of the per-atom values on each call");
  keys.addOutputComponent("q6", "default",
                          "mean local q6 (MODE=all or q6) or mean q6-bar (MODE=q6bar)");
  keys.addOutputComponent("ice", "default", "CHILL+ cubic fraction; not biasable, no derivatives");
  keys.addOutputComponent("ring", "default", "primitive six-ring fraction; not biasable, no derivatives");
  keys.addOutputComponent("hexagonal", "default", "CHILL+ hexagonal count; not biasable, no derivatives");
  keys.addOutputComponent("six", "default", "primitive six-ring count; not biasable, no derivatives");
}

SeamsCV::SeamsCV(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao), cutoff_(3.5) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);
  parse("CUTOFF", cutoff_);
  std::string mode;
  parse("MODE", mode);
  if (mode == "all") {
    mode_ = Mode::all;
  } else if (mode == "q6") {
    mode_ = Mode::q6;
  } else if (mode == "q6bar") {
    mode_ = Mode::q6bar;
  } else {
    error("MODE is all, q6 or q6bar");
  }
  parse("DUMP", dump_);
  checkRead();
  requestAtoms(atoms);
  // The one-argument library symbols. See plumed_components.cpp.
  // q6-bar mixes every neighbour's harmonics, so only q6 is biasable.
  if (mode_ == Mode::q6bar) {
    addPlainComponent(this, "q6");
  } else {
    addDerivativeComponent(this, "q6");
  }
  componentIsNotPeriodic("q6");
  if (mode_ == Mode::all) {
    for (const char *name : {"ice", "ring", "hexagonal", "six"}) {
      addPlainComponent(this, name);
      componentIsNotPeriodic(name);
    }
  }
  log.printf("  d-SEAMS %s, cutoff %f angstrom, mode %s\n", seams_version_string(), cutoff_,
             mode.c_str());
}

SeamsCV::~SeamsCV() {
  if (calls_ > 0) {
    log.printf("  SEAMS calls %ld seconds %.9e per_call %.9e\n", calls_, seconds_,
               seconds_ / static_cast<double>(calls_));
  }
}

double SeamsCV::lengthToAngstrom() {
  // getPosition and getBox are in PLUMED's length unit, in nanometres.
#if PLUMED_VERSION_MAJOR > 2 || (PLUMED_VERSION_MAJOR == 2 && PLUMED_VERSION_MINOR >= 10)
  const double nm = getUnits().getLength();
#else
  const double nm = plumed.getAtoms().getUnits().getLength();
#endif
  return nm * 10.0;
}

bool SeamsCV::score(const std::vector<double> &xyz, const double box9[9], Score &out, bool keep) {
  const int n = static_cast<int>(xyz.size() / 3);
  seams_system *sys = seams_system_create();
  if (sys == nullptr || seams_system_set_frame(sys, xyz.data(), n, box9, nullptr, nullptr) != SEAMS_OK) {
    seams_system_destroy(sys);
    return false;
  }
  peratom_.resize(static_cast<std::size_t>(n));
  bool ok = true;
  if (mode_ == Mode::q6bar) {
    double mean = 0.0;
    ok = seams_q6bar(sys, cutoff_, 0, 0.36, peratom_.data(), peratom_.size(), &mean, nullptr,
                     nullptr) == SEAMS_OK;
    double sum = 0.0;
    for (double value : peratom_) sum += value;
    out.q6 = sum / static_cast<double>(n);
  } else {
    ok = seams_q4q6(sys, cutoff_, 0, nullptr, peratom_.data()) == SEAMS_OK;
    double sum = 0.0;
    for (double value : peratom_) sum += value;
    out.q6 = sum / static_cast<double>(n);
  }
  if (ok && mode_ == Mode::all && keep) {
    labels_.resize(static_cast<std::size_t>(n));
    int counts[9] = {};
    seams_rings *rings = nullptr;
    int hist[7] = {};
    ok = seams_chill_plus(sys, cutoff_, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                          labels_.data(), counts) == SEAMS_OK &&
         seams_rings_find(sys, cutoff_, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK &&
         seams_rings_fraction(rings, 6, &out.ring) == SEAMS_OK &&
         seams_rings_histogram(rings, hist, 7) == SEAMS_OK;
    seams_rings_destroy(rings);
    out.ice = static_cast<double>(counts[SEAMS_LABEL_CUBIC]) / static_cast<double>(n);
    out.hexagonal = counts[SEAMS_LABEL_HEXAGONAL];
    out.six = hist[6];
  }
  seams_system_destroy(sys);
  return ok;
}

void SeamsCV::writeDump(const std::vector<double> &xyz, const double box9[9]) {
  if (dump_.empty() || comm.Get_rank() != 0) {
    return;
  }
  FILE *fp = std::fopen(dump_.c_str(), calls_ == 1 ? "w" : "a");
  if (fp == nullptr) {
    error("cannot write DUMP file " + dump_);
  }
  const unsigned n = getNumberOfAtoms();
  const bool tilted = box9[6] != 0.0 || box9[7] != 0.0 || box9[8] != 0.0;
  std::fprintf(fp, "ITEM: TIMESTEP\n%lld\nITEM: NUMBER OF ATOMS\n%u\n",
               static_cast<long long>(getStep()), n);
  const double xy = box9[6], xz = box9[7], yz = box9[8];
  const double xlo = box9[0], ylo = box9[1], zlo = box9[2];
  const double xhi = xlo + box9[3], yhi = ylo + box9[4], zhi = zlo + box9[5];
  if (tilted) {
    const double xmin = std::min({0.0, xy, xz, xy + xz}), xmax = std::max({0.0, xy, xz, xy + xz});
    const double ymin = std::min(0.0, yz), ymax = std::max(0.0, yz);
    std::fprintf(fp, "ITEM: BOX BOUNDS xy xz yz pp pp pp\n%.17g %.17g %.17g\n%.17g %.17g %.17g\n"
                     "%.17g %.17g %.17g\n",
                 xlo + xmin, xhi + xmax, xy, ylo + ymin, yhi + ymax, xz, zlo, zhi, yz);
  } else {
    std::fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n%.17g %.17g\n%.17g %.17g\n%.17g %.17g\n", xlo,
                 xhi, ylo, yhi, zlo, zhi);
  }
  const char *value = mode_ == Mode::q6bar ? "q6bar" : "q6";
  std::fprintf(fp, "ITEM: ATOMS id type x y z %s%s\n", value, mode_ == Mode::all ? " label" : "");
  for (unsigned i = 0; i < n; ++i) {
    std::fprintf(fp, "%d 1 %.17g %.17g %.17g %.17g", getAbsoluteIndex(i).serial(), xyz[3 * i],
                 xyz[3 * i + 1], xyz[3 * i + 2], peratom_[i]);
    if (mode_ == Mode::all) {
      std::fprintf(fp, " %d", labels_[i]);
    }
    std::fputc('\n', fp);
  }
  std::fclose(fp);
}

void SeamsCV::calculate() {
  const double to_angstrom = lengthToAngstrom();
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
  Score got;
  const auto t0 = std::chrono::steady_clock::now();
  if (!score(xyz, box9, got, true)) {
    error(std::string("d-SEAMS rejected this frame: ") + seams_last_error());
  }
  seconds_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  calls_ += 1;
  writeDump(xyz, box9);
  Value *vq = getPntrToComponent("q6");
  vq->set(got.q6);
  if (mode_ == Mode::all) {
    getPntrToComponent("ice")->set(got.ice);
    getPntrToComponent("ring")->set(got.ring);
    getPntrToComponent("hexagonal")->set(static_cast<double>(got.hexagonal));
    getPntrToComponent("six")->set(static_cast<double>(got.six));
  }
  // A PRINT of the value does not allocate derivative storage. Writing
  // into that empty array aborts. A RESTRAINT does allocate it.
  if (mode_ == Mode::q6bar || vq->getNumberOfDerivatives() == 0) {
    return;
  }
  // Central difference on q6 only. ice and ring have no derivatives.
  const double eps = 1.0e-3;
  std::vector<Vector> dq(n);
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned d = 0; d < 3; ++d) {
      Score plus, minus;
      xyz[3 * i + d] += eps;
      if (!score(xyz, box9, plus, false)) {
        error("d-SEAMS rejected a displaced frame");
      }
      xyz[3 * i + d] -= 2.0 * eps;
      if (!score(xyz, box9, minus, false)) {
        error("d-SEAMS rejected a displaced frame");
      }
      xyz[3 * i + d] += eps;
      // eps is an angstrom step. PLUMED differentiates with respect to its
      // own length unit, to_angstrom angstroms long.
      dq[i][d] = (plus.q6 - minus.q6) / (2.0 * eps) * to_angstrom;
    }
  }
  for (unsigned i = 0; i < n; ++i) {
    setAtomsDerivatives(vq, static_cast<int>(i), dq[i]);
  }
  setBoxDerivativesNoPbc(vq);
}

PLUMED_REGISTER_ACTION(SeamsCV, "SEAMS")

}  // namespace PLMD
