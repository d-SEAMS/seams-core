#ifndef SEAMS_SHELL_HOT_H_
#define SEAMS_SHELL_HOT_H_

// Hot path of the ice classifier. No Eigen types: this header is safe to
// compile with AVX2, which faults when it is mixed with Eigen's allocators.
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT

namespace symm {

// One ideal shell, three templates. first[k][a] is the a-th tetrahedron
// vertex of template k, relative to the centre. arm[k][a] is the three
// second-shell atoms on that bond, ordered by azimuth.
struct CanonShell {
  double first[3][4][3];
  double second[3][12][3];
  int parent[3][12];
  int arm[3][4][3];
  int nArm[3][4];
  double pts[3][17][3];
};

// Seventeen oxygens: the centre, the tetrahedron, and the twelve outer atoms.
struct Shell17 {
  double p[17][3];
};

// xyz is n by 3, row major. lat rows are the cell edges. kind[i] is written
// as ice Ih (0), ice Ic (1), ice Isd (2), or -1.
// x, y, z are length n. Column-major coordinates, no interleaving copy.
void classifyPrepared(const double *x, const double *y, const double *z, int n,
                      const double lat[3][3], const CanonShell &canon, double gate, int *kind);

} // namespace symm

#endif
