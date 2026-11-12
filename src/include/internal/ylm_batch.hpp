#ifndef SEAMS_YLM_BATCH_H_
#define SEAMS_YLM_BATCH_H_

#include <ylm_cartesian.hpp>

namespace seams {
namespace ylm {

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

} // namespace ylm
} // namespace seams

#endif
