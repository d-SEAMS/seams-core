#include "compute_dseams.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "update.h"

#include <seams_c_api.h>

#include <cstring>
#include <vector>

using namespace LAMMPS_NS;

ComputeDseams::ComputeDseams(LAMMPS *lmp, int narg, char **arg)
    : Compute(lmp, narg, arg), nmax(0), chill(nullptr) {
  if (narg < 3) {
    error->all(FLERR, "Illegal compute dseams command");
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
  const double box[6] = {domain->boxlo[0], domain->boxlo[1], domain->boxlo[2],
                         domain->xprd, domain->yprd, domain->zprd};
  if (seams_chill_plus(xyz.data(), nlocal, box, labels.data()) != 0) {
    error->one(FLERR, "compute dseams: seams_chill_plus failed");
  }
  for (int i = 0; i < nlocal; i++) {
    chill[i] = static_cast<double>(labels[static_cast<std::size_t>(i)]);
  }
}

double ComputeDseams::memory_usage() {
  return static_cast<double>(nmax) * sizeof(double);
}
