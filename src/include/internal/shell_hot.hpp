#ifndef SEAMS_SHELL_HOT_H_
#define SEAMS_SHELL_HOT_H_

// Hot path of the ice classifier. No Eigen types: this header is safe to
// compile with AVX2, which faults when it is mixed with Eigen's allocators.
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT

namespace symm {

// One ideal shell, three templates. first[k][a] is the a-th tetrahedron
// vertex of template k, relative to the center. arm[k][a] is the three
// second-shell atoms on that bond, ordered by azimuth.
struct CanonShell {
  double first[3][4][3];
  double second[3][12][3];
  int parent[3][12];
  int arm[3][4][3];
  int nArm[3][4];
  double pts[3][17][3];
};

// Seventeen oxygens: the center, the tetrahedron, and the twelve outer atoms.
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
// openFace completes a short second shell on a bonded neighbour that has
// fewer than four bonds. A full coordination shell is never padded. Off
// by default, which keeps the labels the shell match gave before.
void classifyPrepared(const double *x, const double *y, const double *z, int n,
                      const double lat[3][3], const CanonShell &canon, double gate, int *kind,
                      double thermalQ = kThermalQ, double bondRatio = kThermalBondRatio,
                      bool openFace = false);

// Common-neighbour signatures of one atom's twelve second neighbours, given
// relative to the atom: 1 for twelve (4,2,1), 0 for six (4,2,1) and six
// (4,2,2), 3 or 4 for a mostly cubic or mostly hexagonal shell with eight
// good signatures, -2 for at least `ambiguous` good ones, -1 otherwise.
int diamondCna(const double *sx, const double *sy, const double *sz, int ambiguous,
               int *n421 = nullptr, int *n422 = nullptr, int *nGood = nullptr);

// True when the AVX2 build of shell_classify.cpp is linked in, the CPU
// reports AVX2 and FMA, and seams::cpu::allowSimdKernels has not turned it
// off. classifyPrepared and diamondCna call symm::avx2 only then, and
// symm::generic otherwise.
bool avx2Kernels();

// The two builds of shell_classify.cpp. Call them through the functions
// above; symm::avx2 exists only where avx2Kernels() can be true.
namespace generic {
void classifyPrepared(const double *x, const double *y, const double *z, int n,
                      const double lat[3][3], const CanonShell &canon, double gate, int *kind,
                      double thermalQ, double bondRatio, bool openFace);
int diamondCna(const double *sx, const double *sy, const double *sz, int ambiguous, int *n421,
               int *n422, int *nGood);
} // namespace generic

namespace avx2 {
void classifyPrepared(const double *x, const double *y, const double *z, int n,
                      const double lat[3][3], const CanonShell &canon, double gate, int *kind,
                      double thermalQ, double bondRatio, bool openFace);
int diamondCna(const double *sx, const double *sy, const double *sz, int ambiguous, int *n421,
               int *n422, int *nGood);
} // namespace avx2

} // namespace symm

#endif
