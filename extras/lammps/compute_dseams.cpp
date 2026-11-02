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
      graph("cutoff") {
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
  std::vector<int> labels(static_cast<std::size_t>(nlocal), 6);
  for (int i = 0; i < nlocal; i++) {
    xyz[static_cast<std::size_t>(i) * 3 + 0] = atom->x[i][0];
    xyz[static_cast<std::size_t>(i) * 3 + 1] = atom->x[i][1];
    xyz[static_cast<std::size_t>(i) * 3 + 2] = atom->x[i][2];
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
  if (seams_chill_plus_graph_lammps(xyz.data(), nlocal, box9, cutoff,
                                   graph.c_str(), labels.data()) != 0) {
    error->one(FLERR, "compute dseams: seams_chill_plus_graph failed");
  }
  for (int i = 0; i < nlocal; i++) {
    chill[i] = static_cast<double>(labels[static_cast<std::size_t>(i)]);
  }
}

double ComputeDseams::memory_usage() {
  return static_cast<double>(nmax) * sizeof(double);
}
