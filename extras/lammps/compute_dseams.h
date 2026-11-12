#ifdef COMPUTE_CLASS
ComputeStyle(dseams, ComputeDseams)
#else

#ifndef LMP_COMPUTE_DSEAMS_H
#define LMP_COMPUTE_DSEAMS_H

#include "compute.h"

#include <string>

namespace LAMMPS_NS {

class ComputeDseams : public Compute {
public:
  ComputeDseams(class LAMMPS *, int, char **);
  void init() override;
  void compute_peratom() override;
  double memory_usage() override;

private:
  int nmax;
  double *chill;
  double cutoff;
  std::string graph;
  int oxygen_type;
  int mol_mode;
  int mol_id;
};

} // namespace LAMMPS_NS

#endif
#endif
