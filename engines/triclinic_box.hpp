// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// Lower-triangular cells as the nine LAMMPS numbers ox,oy,oz,lx,ly,lz,xy,xz,yz.
//
// GROMACS and PLUMED both store the lattice vectors as rows. A
// lower-triangular cell (a along x, b in the xy plane, c with all three
// components) keeps the upper triangle at zero, which GROMACS requires:
//   xy = box[1][0], xz = box[2][0], yz = box[2][1].

#pragma once

inline void gromacs_box_to_lammps(const double box[3][3], double out[9], double scale) {
  out[0] = 0.0;
  out[1] = 0.0;
  out[2] = 0.0;
  out[3] = box[0][0] * scale;
  out[4] = box[1][1] * scale;
  out[5] = box[2][2] * scale;
  out[6] = box[1][0] * scale;
  out[7] = box[2][0] * scale;
  out[8] = box[2][1] * scale;
}

inline void plumed_box_to_lammps(const double box[3][3], double out[9]) {
  out[0] = 0.0;
  out[1] = 0.0;
  out[2] = 0.0;
  out[3] = box[0][0];
  out[4] = box[1][1];
  out[5] = box[2][2];
  out[6] = box[1][0];
  out[7] = box[2][0];
  out[8] = box[2][1];
}
