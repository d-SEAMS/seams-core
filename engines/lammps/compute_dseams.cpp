// Copyright (c) d-SEAMS developers
#include "compute_dseams.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "update.h"

#include "seams.h"

#include <cstring>
#include <vector>

using namespace LAMMPS_NS;

ComputeDseams::ComputeDseams(LAMMPS *lmp, int narg, char **arg)
    : Compute(lmp, narg, arg), nmax_(0), label_(nullptr), cutoff_(3.5),
      graph_(SEAMS_GRAPH_CUTOFF), rule_(SEAMS_CHILL_PUBLISHED) {
  if (narg < 3) error->all(FLERR, "Illegal compute dseams command");
  for (int iarg = 3; iarg < narg; iarg += 2) {
    if (iarg + 1 >= narg) error->all(FLERR, "Illegal compute dseams command");
    if (strcmp(arg[iarg], "cutoff") == 0) {
      cutoff_ = strtod(arg[iarg + 1], nullptr);
    } else if (strcmp(arg[iarg], "graph") == 0) {
      graph_ = strcmp(arg[iarg + 1], "mutual") == 0 ? SEAMS_GRAPH_MUTUAL : SEAMS_GRAPH_CUTOFF;
    } else if (strcmp(arg[iarg], "rule") == 0) {
      rule_ = strcmp(arg[iarg + 1], "ovito") == 0 ? SEAMS_CHILL_OVITO : SEAMS_CHILL_PUBLISHED;
    } else {
      error->all(FLERR, "Illegal compute dseams command");
    }
  }
  peratom_flag = 1;
  size_peratom_cols = 0;
}

ComputeDseams::~ComputeDseams() {
  memory->destroy(label_);
}

void ComputeDseams::init() {}

void ComputeDseams::compute_peratom() {
  invoked_peratom = update->ntimestep;
  const int nlocal = atom->nlocal;
  if (nlocal > nmax_) {
    memory->destroy(label_);
    nmax_ = nlocal;
    memory->create(label_, nmax_, "dseams:label");
    vector_atom = label_;
  }
  std::vector<double> xyz(static_cast<std::size_t>(nlocal) * 3);
  std::vector<int> labels(static_cast<std::size_t>(nlocal));
  for (int i = 0; i < nlocal; ++i) {
    xyz[static_cast<std::size_t>(3 * i)] = atom->x[i][0];
    xyz[static_cast<std::size_t>(3 * i + 1)] = atom->x[i][1];
    xyz[static_cast<std::size_t>(3 * i + 2)] = atom->x[i][2];
  }
  const double box9[9] = {domain->boxlo[0],
                          domain->boxlo[1],
                          domain->boxlo[2],
                          domain->xprd,
                          domain->yprd,
                          domain->zprd,
                          domain->triclinic ? domain->xy : 0.0,
                          domain->triclinic ? domain->xz : 0.0,
                          domain->triclinic ? domain->yz : 0.0};
  seams_system *sys = seams_system_create();
  if (sys == nullptr ||
      seams_system_set_frame(sys, xyz.data(), nlocal, box9, nullptr, nullptr) != SEAMS_OK ||
      seams_chill_plus(sys, cutoff_, 0, graph_, rule_, labels.data(), nullptr) != SEAMS_OK) {
    seams_system_destroy(sys);
    error->one(FLERR, "compute dseams failed");
  }
  seams_system_destroy(sys);
  for (int i = 0; i < nlocal; ++i) label_[i] = static_cast<double>(labels[static_cast<std::size_t>(i)]);
}

double ComputeDseams::memory_usage() {
  return static_cast<double>(nmax_) * sizeof(double);
}
