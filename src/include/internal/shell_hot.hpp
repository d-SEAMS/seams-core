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

// First-shell prefilter. It is not the crystal decision: that is the
// second-shell fingerprint and the scale-free fit. On TIP4P/2005 at its
// melting point a quench-majority crystal has tetrahedral-order 5th
// percentile 0.87 and the liquid median is 0.80, so 0.85 lets the thermal
// crystal reach the second shell and still sits above the liquid median.
// 0.94 is the floor of a quenched mW crystal and cuts that thermal body.
inline constexpr double kThermalQ = 0.85;
// Ratio of the longest to the shortest of the four nearest bonds. The
// same crystal's 95th percentile is 1.17.
inline constexpr double kThermalBondRatio = 1.20;

// xyz is n by 3, row major. lat rows are the cell edges. kind[i] is written
// as ice Ih (0), ice Ic (1), ice Isd (2), or -1.
// x, y, z are length n. Column-major coordinates, no interleaving copy.
void classifyPrepared(const double *x, const double *y, const double *z, int n,
                      const double lat[3][3], const CanonShell &canon, double gate, int *kind,
                      double thermalQ = kThermalQ, double bondRatio = kThermalBondRatio);

} // namespace symm

#endif
