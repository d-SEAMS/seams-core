#ifndef SEAMS_YLM_BATCH_H_
#define SEAMS_YLM_BATCH_H_

#include <ylm_cartesian.hpp>

namespace seams {
namespace ylm {

// Bond-major output: component m of bond b is at b * (2 * orderL + 1) + m.
// Host only. l is 3, 4 or 6.
void ylmBatch(int orderL, const double *dx, const double *dy, const double *dz,
              int n, double *re, double *im);

} // namespace ylm
} // namespace seams

#endif
