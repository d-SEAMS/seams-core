// Copyright (c) d-SEAMS developers
// Lower-triangular cells as the nine LAMMPS numbers ox,oy,oz,lx,ly,lz,xy,xz,yz.
//
// GROMACS stores box vectors as columns of box[dimension][vector], so a
// lower-triangular cell (a along x, b in the xy plane, c with all three
// components) has its tilts in the upper triangle:
//   xy = box[0][1], xz = box[0][2], yz = box[1][2].
// PLUMED stores the same cell with the lattice vectors as rows, so the
// tilts sit in the lower triangle:
//   xy = box(1,0), xz = box(2,0), yz = box(2,1).

#pragma once

inline void gromacs_box_to_lammps(const double box[3][3], double out[9], double scale) {
  out[0] = 0.0;
  out[1] = 0.0;
  out[2] = 0.0;
  out[3] = box[0][0] * scale;
  out[4] = box[1][1] * scale;
  out[5] = box[2][2] * scale;
  out[6] = box[0][1] * scale;
  out[7] = box[0][2] * scale;
  out[8] = box[1][2] * scale;
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
