// Copyright (c) d-SEAMS developers
#include "../engines/triclinic_box.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

static int near(double got, double want) {
  if (!(std::abs(got - want) < 1e-12)) {
    std::cerr << "got " << got << " want " << want << "\n";
    return 1;
  }
  return 0;
}

int main() {
  int fails = 0;
  double gromacs[3][3] = {};
  gromacs[0][0] = 10.0;
  gromacs[1][1] = 8.0;
  gromacs[2][2] = 9.0;
  gromacs[1][0] = 1.5;
  gromacs[2][0] = 0.4;
  gromacs[2][1] = -0.2;
  double box[9];
  gromacs_box_to_lammps(gromacs, box, 10.0);
  fails += near(box[3], 100.0);
  fails += near(box[4], 80.0);
  fails += near(box[5], 90.0);
  fails += near(box[6], 15.0);
  fails += near(box[7], 4.0);
  fails += near(box[8], -2.0);

  double plumed[3][3] = {};
  plumed[0][0] = 10.0;
  plumed[1][1] = 8.0;
  plumed[2][2] = 9.0;
  plumed[1][0] = 1.5;
  plumed[2][0] = 0.4;
  plumed[2][1] = -0.2;
  plumed_box_to_lammps(plumed, box);
  fails += near(box[3], 10.0);
  fails += near(box[4], 8.0);
  fails += near(box[5], 9.0);
  fails += near(box[6], 1.5);
  fails += near(box[7], 0.4);
  fails += near(box[8], -0.2);
  return fails == 0 ? 0 : 1;
}
