// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
#ifdef COMPUTE_CLASS
ComputeStyle(dseams, ComputeDseams)
#else
#ifndef LMP_COMPUTE_DSEAMS_H
#define LMP_COMPUTE_DSEAMS_H

#include "compute.h"

namespace LAMMPS_NS {

// Per-atom CHILL+ label (the seams_label integer). Atoms outside the group
// read -1. `type I` or `mol water` picks the water oxygens of the group and
// the other group atoms are unclassified (6); a group of several types
// needs one of them. Every rank sends its owned atoms to rank 0, which
// classifies the whole group and returns the labels, so the result does
// not depend on the domain decomposition.
class ComputeDseams : public Compute {
 public:
  ComputeDseams(class LAMMPS *, int, char **);
  ~ComputeDseams() override;
  void init() override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nmax_;
  double *label_;
  double cutoff_;
  int graph_;
  int rule_;
  int oxygen_type_;
  bool mol_water_;
  long calls_;
  double seconds_;
};

}  // namespace LAMMPS_NS
#endif
#endif
