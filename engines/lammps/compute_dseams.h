// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
#ifdef COMPUTE_CLASS
ComputeStyle(dseams, ComputeDseams)
#else
#ifndef LMP_COMPUTE_DSEAMS_H
#define LMP_COMPUTE_DSEAMS_H

#include "compute.h"

namespace LAMMPS_NS {

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
};

}  // namespace LAMMPS_NS
#endif
#endif
