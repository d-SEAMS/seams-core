#include "compute_dseams.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "update.h"

#include <seams_c_api.h>

#include <cstdlib>
#include <cstring>
#include <vector>

using namespace LAMMPS_NS;

ComputeDseams::ComputeDseams(LAMMPS *lmp, int narg, char **arg)
    : Compute(lmp, narg, arg), nmax(0), chill(nullptr), cutoff(3.5),
      graph("cutoff"), oxygen_type(0), mol_mode(0), mol_id(0) {
  if (narg < 3) {
    error->all(FLERR, "Illegal compute dseams command");
  }
  for (int iarg = 3; iarg < narg; iarg += 2) {
    if (iarg + 1 >= narg) {
      error->all(FLERR, "Illegal compute dseams command");
    }
    if (std::strcmp(arg[iarg], "graph") == 0) {
      graph = arg[iarg + 1];
      if (graph != "cutoff" && graph != "mutual") {
        error->all(FLERR, "compute dseams: graph is cutoff or mutual");
      }
    } else if (std::strcmp(arg[iarg], "cutoff") == 0) {
      cutoff = std::strtod(arg[iarg + 1], nullptr);
      if (!(cutoff > 0.0)) {
        error->all(FLERR, "compute dseams: cutoff must be positive");
      }
    } else if (std::strcmp(arg[iarg], "type") == 0) {
      oxygen_type = std::atoi(arg[iarg + 1]);
      if (oxygen_type <= 0) {
        error->all(FLERR, "compute dseams: type must be a positive atom type");
      }
    } else if (std::strcmp(arg[iarg], "mol") == 0) {
      if (std::strcmp(arg[iarg + 1], "water") == 0) {
        mol_mode = 1;
      } else {
        mol_id = std::atoi(arg[iarg + 1]);
        if (mol_id == 0 && arg[iarg + 1][0] != '0') {
          error->all(FLERR, "compute dseams: mol is water or a molecule id");
        }
        mol_mode = 2;
      }
    } else {
      error->all(FLERR, "Illegal compute dseams command");
    }
  }
  peratom_flag = 1;
  size_peratom_cols = 0;
}

void ComputeDseams::init() {}

void ComputeDseams::compute_peratom() {
  invoked_peratom = update->ntimestep;
  const int nlocal = atom->nlocal;
  if (nlocal > nmax) {
    memory->destroy(chill);
    nmax = nlocal;
    memory->create(chill, nmax, "dseams:chill");
    vector_atom = chill;
  }
  std::vector<double> xyz(static_cast<std::size_t>(nlocal) * 3);
  std::vector<int> types(static_cast<std::size_t>(nlocal), 1);
  std::vector<int> mols(static_cast<std::size_t>(nlocal), 0);
  std::vector<int> labels(static_cast<std::size_t>(nlocal), 6);
  const int *molArg = nullptr;
  for (int i = 0; i < nlocal; i++) {
    xyz[static_cast<std::size_t>(i) * 3 + 0] = atom->x[i][0];
    xyz[static_cast<std::size_t>(i) * 3 + 1] = atom->x[i][1];
    xyz[static_cast<std::size_t>(i) * 3 + 2] = atom->x[i][2];
    types[static_cast<std::size_t>(i)] = atom->type[i];
    mols[static_cast<std::size_t>(i)] =
        atom->molecule != nullptr ? static_cast<int>(atom->molecule[i]) : i + 1;
  }
  if (atom->molecule != nullptr || mol_mode != 0) {
    molArg = mols.data();
  }
  // Edge origin and edge lengths, plus the tilt. Passing only xprd
  // as an orthorhombic length drops xy, xz and yz.
  const double box9[9] = {
      domain->boxlo[0],
      domain->boxlo[1],
      domain->boxlo[2],
      domain->xprd,
      domain->yprd,
      domain->zprd,
      domain->triclinic ? domain->xy : 0.0,
      domain->triclinic ? domain->xz : 0.0,
      domain->triclinic ? domain->yz : 0.0};
  const int *oneMol = mol_mode == 2 ? &mol_id : nullptr;
  if (seams_chill_plus_species(xyz.data(), nlocal, box9, types.data(), molArg,
                              oxygen_type, mol_mode, oneMol,
                              mol_mode == 2 ? 1 : 0, cutoff, graph.c_str(),
                              labels.data(), nullptr, nullptr, nullptr) != 0) {
    const char *why = seams_c_error();
    error->one(FLERR, why != nullptr && why[0] != '\0'
                          ? why
                          : "compute dseams: seams_chill_plus_species failed");
  }
  for (int i = 0; i < nlocal; i++) {
    chill[i] = static_cast<double>(labels[static_cast<std::size_t>(i)]);
  }
}

double ComputeDseams::memory_usage() {
  return static_cast<double>(nmax) * sizeof(double);
}
