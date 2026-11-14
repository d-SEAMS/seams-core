// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
#include "compute_dseams.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "update.h"
#include "utils.h"

#include "seams.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

using namespace LAMMPS_NS;

namespace {

// The water oxygen type of the group under `mol water`: the one type that
// occurs exactly once in every multi-atom molecule. With no multi-atom
// molecule the group must be one type. The rules and messages follow
// water::selectOxygens in the seams CLI.
std::string waterType(const std::vector<int> &types, const std::vector<tagint> &mols,
                      int &type) {
  std::map<tagint, std::vector<int>> groups;
  for (std::size_t k = 0; k < types.size(); ++k) groups[mols[k]].push_back(types[k]);
  std::set<int> votes;
  int multi = 0;
  int ambiguous = 0;
  for (const auto &item : groups) {
    if (item.second.size() < 2) continue;
    ++multi;
    std::map<int, int> count;
    for (int t : item.second) count[t] += 1;
    std::vector<int> singles;
    for (const auto &entry : count) {
      if (entry.second == 1) singles.push_back(entry.first);
    }
    if (singles.size() == 1) {
      votes.insert(singles[0]);
    } else if (singles.size() > 1) {
      ++ambiguous;
    }
  }
  if (multi > 0) {
    if (votes.size() == 1 && ambiguous == 0) {
      type = *votes.begin();
      return "";
    }
    if (votes.size() > 1 || ambiguous > 0) return "molecule has several singleton sites; pass type";
    return "molecule has no water oxygen; pass type";
  }
  const std::set<int> distinct(types.begin(), types.end());
  if (distinct.size() == 1) {
    type = *distinct.begin();
    return "";
  }
  return "mixed species need a water type or molecule";
}

}  // namespace

ComputeDseams::ComputeDseams(LAMMPS *lmp, int narg, char **arg)
    : Compute(lmp, narg, arg), nmax_(0), label_(nullptr), cutoff_(3.5),
      graph_(SEAMS_GRAPH_CUTOFF), rule_(SEAMS_CHILL_PUBLISHED), oxygen_type_(0),
      mol_water_(false), calls_(0), seconds_(0.0) {
  if (narg < 3) error->all(FLERR, "Illegal compute dseams command");
  for (int iarg = 3; iarg < narg; iarg += 2) {
    if (iarg + 1 >= narg) error->all(FLERR, "Illegal compute dseams command");
    if (strcmp(arg[iarg], "cutoff") == 0) {
      cutoff_ = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (!(cutoff_ > 0.0)) error->all(FLERR, "compute dseams: cutoff must be positive");
    } else if (strcmp(arg[iarg], "graph") == 0) {
      if (strcmp(arg[iarg + 1], "mutual") == 0) {
        graph_ = SEAMS_GRAPH_MUTUAL;
      } else if (strcmp(arg[iarg + 1], "cutoff") == 0) {
        graph_ = SEAMS_GRAPH_CUTOFF;
      } else {
        error->all(FLERR, "compute dseams: graph is cutoff or mutual");
      }
    } else if (strcmp(arg[iarg], "rule") == 0) {
      if (strcmp(arg[iarg + 1], "ovito") == 0) {
        rule_ = SEAMS_CHILL_OVITO;
      } else if (strcmp(arg[iarg + 1], "published") == 0) {
        rule_ = SEAMS_CHILL_PUBLISHED;
      } else {
        error->all(FLERR, "compute dseams: rule is published or ovito");
      }
    } else if (strcmp(arg[iarg], "type") == 0) {
      oxygen_type_ = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (oxygen_type_ <= 0 || oxygen_type_ > atom->ntypes)
        error->all(FLERR, "compute dseams: type must be an atom type of this system");
    } else if (strcmp(arg[iarg], "mol") == 0) {
      if (strcmp(arg[iarg + 1], "water") != 0) error->all(FLERR, "compute dseams: mol is water");
      if (!atom->molecule_flag)
        error->all(FLERR, "compute dseams: mol water needs an atom style with molecule IDs");
      mol_water_ = true;
    } else {
      error->all(FLERR, "Illegal compute dseams command");
    }
  }
  peratom_flag = 1;
  size_peratom_cols = 0;
}

ComputeDseams::~ComputeDseams() {
  memory->destroy(label_);
  if (comm->me == 0 && calls_ > 0) {
    char line[160];
    std::snprintf(line, sizeof(line),
                  "compute dseams calls %ld seconds %.9e per_call %.9e\n", calls_, seconds_,
                  seconds_ / static_cast<double>(calls_));
    utils::logmesg(lmp, line);
  }
}

void ComputeDseams::init() {}

void ComputeDseams::compute_peratom() {
  invoked_peratom = update->ntimestep;
  const double t0 = MPI_Wtime();
  if (atom->nmax > nmax_) {
    memory->destroy(label_);
    nmax_ = atom->nmax;
    memory->create(label_, nmax_, "dseams:label");
    vector_atom = label_;
  }
  const int nlocal = atom->nlocal;
  std::vector<tagint> tags;
  std::vector<tagint> mols;
  std::vector<int> types;
  std::vector<double> xyz;
  for (int i = 0; i < nlocal; ++i) {
    label_[i] = -1.0;
    if (!(atom->mask[i] & groupbit)) continue;
    tags.push_back(atom->tag[i]);
    types.push_back(atom->type[i]);
    mols.push_back(atom->molecule_flag ? atom->molecule[i] : atom->tag[i]);
    xyz.insert(xyz.end(), {atom->x[i][0], atom->x[i][1], atom->x[i][2]});
  }
  const int nprocs = comm->nprocs;
  const int me = comm->me;
  int mine = static_cast<int>(tags.size());
  std::vector<int> counts(nprocs), displs(nprocs), counts3(nprocs), displs3(nprocs);
  MPI_Gather(&mine, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, world);
  int total = 0;
  if (me == 0) {
    for (int p = 0; p < nprocs; ++p) {
      displs[p] = total;
      counts3[p] = 3 * counts[p];
      displs3[p] = 3 * total;
      total += counts[p];
    }
  }
  std::vector<tagint> all_tags(me == 0 ? total : 0);
  std::vector<tagint> all_mols(me == 0 ? total : 0);
  std::vector<int> all_types(me == 0 ? total : 0);
  std::vector<double> all_xyz(me == 0 ? 3 * static_cast<std::size_t>(total) : 0);
  MPI_Gatherv(tags.data(), mine, MPI_LMP_TAGINT, all_tags.data(), counts.data(), displs.data(),
              MPI_LMP_TAGINT, 0, world);
  MPI_Gatherv(mols.data(), mine, MPI_LMP_TAGINT, all_mols.data(), counts.data(), displs.data(),
              MPI_LMP_TAGINT, 0, world);
  MPI_Gatherv(types.data(), mine, MPI_INT, all_types.data(), counts.data(), displs.data(),
              MPI_INT, 0, world);
  MPI_Gatherv(xyz.data(), 3 * mine, MPI_DOUBLE, all_xyz.data(), counts3.data(), displs3.data(),
              MPI_DOUBLE, 0, world);
  MPI_Bcast(&total, 1, MPI_INT, 0, world);

  // Rank 0 classifies the water oxygens of the group in atom-ID order and
  // broadcasts the labels. Group atoms that are not water oxygens are
  // unclassified.
  std::vector<tagint> sorted_tags(total);
  std::vector<int> labels(total, SEAMS_LABEL_UNCLASSIFIED);
  int status = SEAMS_OK;
  std::string why;
  if (me == 0 && total > 0) {
    std::vector<int> order(total);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return all_tags[a] < all_tags[b]; });
    std::vector<int> sorted_types(total);
    std::vector<tagint> sorted_mols(total);
    for (int k = 0; k < total; ++k) {
      sorted_tags[k] = all_tags[order[k]];
      sorted_types[k] = all_types[order[k]];
      sorted_mols[k] = all_mols[order[k]];
    }
    int keep_type = oxygen_type_;
    if (keep_type > 0) {
      if (std::find(sorted_types.begin(), sorted_types.end(), keep_type) == sorted_types.end())
        why = "atom type " + std::to_string(keep_type) + " is not in the group";
    } else if (mol_water_) {
      why = waterType(sorted_types, sorted_mols, keep_type);
    } else if (std::any_of(sorted_types.begin(), sorted_types.end(),
                           [&](int t) { return t != sorted_types[0]; })) {
      why = "mixed species need a water type or molecule";
    }
    std::vector<int> picked;
    for (int k = 0; k < total && why.empty(); ++k) {
      if (keep_type == 0 || sorted_types[k] == keep_type) picked.push_back(k);
    }
    const int n = static_cast<int>(picked.size());
    std::vector<double> sorted_xyz(3 * static_cast<std::size_t>(n));
    std::vector<int> ids(n);
    for (int k = 0; k < n; ++k) {
      const int src = order[picked[k]];
      ids[k] = static_cast<int>(all_tags[src]);
      for (int d = 0; d < 3; ++d) sorted_xyz[3 * k + d] = all_xyz[3 * src + d];
    }
    std::vector<int> picked_labels(n);
    const double box9[9] = {domain->boxlo[0],
                            domain->boxlo[1],
                            domain->boxlo[2],
                            domain->xprd,
                            domain->yprd,
                            domain->zprd,
                            domain->triclinic ? domain->xy : 0.0,
                            domain->triclinic ? domain->xz : 0.0,
                            domain->triclinic ? domain->yz : 0.0};
    seams_system *sys = why.empty() ? seams_system_create() : nullptr;
    if (!why.empty()) {
      status = SEAMS_ERR_ARG;
    } else if (sys == nullptr) {
      status = SEAMS_ERR_NOMEM;
    } else {
      status = seams_system_set_frame(sys, sorted_xyz.data(), n, box9, nullptr, ids.data());
      if (status == SEAMS_OK) {
        status = seams_chill_plus(sys, cutoff_, 0, graph_, rule_, picked_labels.data(), nullptr);
      }
      seams_system_destroy(sys);
    }
    for (int k = 0; k < n && status == SEAMS_OK; ++k) labels[picked[k]] = picked_labels[k];
  }
  int selection_failed = why.empty() ? 0 : 1;
  MPI_Bcast(&status, 1, MPI_INT, 0, world);
  MPI_Bcast(&selection_failed, 1, MPI_INT, 0, world);
  // Only rank 0 has the message; error->all prints from rank 0.
  if (selection_failed) error->all(FLERR, "compute dseams: " + why);
  if (status != SEAMS_OK) {
    error->all(FLERR, std::string("compute dseams: ") + seams_status_string(
                          static_cast<seams_status>(status)) + ": " + seams_last_error());
  }
  if (total > 0) {
    MPI_Bcast(sorted_tags.data(), total, MPI_LMP_TAGINT, 0, world);
    MPI_Bcast(labels.data(), total, MPI_INT, 0, world);
  }
  for (int i = 0; i < nlocal; ++i) {
    if (!(atom->mask[i] & groupbit)) continue;
    const auto at = std::lower_bound(sorted_tags.begin(), sorted_tags.end(), atom->tag[i]);
    label_[i] = static_cast<double>(labels[static_cast<std::size_t>(at - sorted_tags.begin())]);
  }
  seconds_ += MPI_Wtime() - t0;
  calls_ += 1;
}

double ComputeDseams::memory_usage() {
  return static_cast<double>(nmax_) * sizeof(double);
}
