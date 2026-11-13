#ifndef SEAMS_YLM_BATCH_H_
#define SEAMS_YLM_BATCH_H_

#include <ylm_cartesian.hpp>

namespace seams {
namespace ylm {

// Add up to eight l = 12 bonds into one interleaved q_lm row (real, imag).
// n is 1..8. Shorter tails do not contribute. Same values as the scalar
// solid harmonic.
void ylm12AddBonds(const double *dx, const double *dy, const double *dz, int n,
                   double *qlmInterleaved);

// Bond-major output: component m of bond b is at b * (2 * orderL + 1) + m.
// Host only. l is 3, 4 or 6.
void ylmBatch(int orderL, const double *dx, const double *dy, const double *dz,
              int n, double *re, double *im);

// Eight bonds, AVX-512. The first four sum into accA and the next four
// into accB. Each accumulator is interleaved real, imag, length
// 2 * (2 * l + 1). l is 4 or 6. The recurrence is the Cartesian
// polynomial, not a trigonometric Y_lm.
void ylmAdd8(int orderL, const double *dx, const double *dy, const double *dz,
             double *accA, double *accB);

// Eight atoms, four bonds each. dx, dy and dz are 4 * 8 values, bond
// slot in the outer index, and each group of eight is 64-byte aligned.
// The solid harmonics stay in SoA while the four bonds are added, then
// one blocked transpose writes atom a at rows[a * stride + 2 * m].
// l is 4 or 6. Returns false if any bond has zero length.
bool ylmSumRows4(int orderL, const double *dx, const double *dy,
                 const double *dz, double *rows, int stride);

// Eight atoms, a variable shell. dx/dy/dz are bond-major, nbonds groups
// of eight. nUsed[a] is the number of leading live bonds of atom a;
// a later slot is zero and is not summed. Writes ql[8] for l = 4 or 6.
// The sum is the same Cartesian solid harmonic as ylmSumRows4, added in
// bond order, then divided by nUsed.
bool ylmLocalQl(int orderL, const double *dx, const double *dy, const double *dz,
                int nbonds, const int *nUsed, double prefactor, double *ql);

} // namespace ylm
} // namespace seams

#endif
