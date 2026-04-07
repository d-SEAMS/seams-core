#ifndef __SIMD_DISTANCE_H_
#define __SIMD_DISTANCE_H_

#include <cmath>

#ifdef SEAMS_HAS_HWY

#include "hwy/highway.h"

namespace seams {

namespace hn = hwy::HWY_NAMESPACE;

// Compute squared periodic distances for a batch of atom pairs.
// dx, dy, dz: coordinate differences (length n)
// bx, by, bz: periodic box dimensions
// out: output squared distances (length n)
inline void BatchPeriodicDistSq(const double* HWY_RESTRICT dx,
                                const double* HWY_RESTRICT dy,
                                const double* HWY_RESTRICT dz,
                                double bx, double by, double bz,
                                double* HWY_RESTRICT out, size_t n) {
  const hn::ScalableTag<double> d;
  const size_t N = hn::Lanes(d);

  const auto vbx = hn::Set(d, bx);
  const auto vby = hn::Set(d, by);
  const auto vbz = hn::Set(d, bz);

  size_t i = 0;
  for (; i + N <= n; i += N) {
    auto vdx = hn::Load(d, dx + i);
    auto vdy = hn::Load(d, dy + i);
    auto vdz = hn::Load(d, dz + i);

    // Absolute values
    vdx = hn::Abs(vdx);
    vdy = hn::Abs(vdy);
    vdz = hn::Abs(vdz);

    // Periodic wrap: dr -= box * round(dr / box)
    vdx = hn::Sub(vdx, hn::Mul(vbx, hn::Round(hn::Div(vdx, vbx))));
    vdy = hn::Sub(vdy, hn::Mul(vby, hn::Round(hn::Div(vdy, vby))));
    vdz = hn::Sub(vdz, hn::Mul(vbz, hn::Round(hn::Div(vdz, vbz))));

    // r2 = dx*dx + dy*dy + dz*dz
    auto r2 = hn::MulAdd(vdx, vdx, hn::MulAdd(vdy, vdy, hn::Mul(vdz, vdz)));
    hn::Store(r2, d, out + i);
  }

  // Scalar remainder
  for (; i < n; i++) {
    double ddx = std::fabs(dx[i]);
    double ddy = std::fabs(dy[i]);
    double ddz = std::fabs(dz[i]);
    ddx -= bx * std::round(ddx / bx);
    ddy -= by * std::round(ddy / by);
    ddz -= bz * std::round(ddz / bz);
    out[i] = ddx * ddx + ddy * ddy + ddz * ddz;
  }
}

}  // namespace seams

#else  // !SEAMS_HAS_HWY

// Scalar fallback
namespace seams {
inline void BatchPeriodicDistSq(const double* dx, const double* dy,
                                const double* dz, double bx, double by,
                                double bz, double* out, size_t n) {
  for (size_t i = 0; i < n; i++) {
    double ddx = std::fabs(dx[i]);
    double ddy = std::fabs(dy[i]);
    double ddz = std::fabs(dz[i]);
    ddx -= bx * std::round(ddx / bx);
    ddy -= by * std::round(ddy / by);
    ddz -= bz * std::round(ddz / bz);
    out[i] = ddx * ddx + ddy * ddy + ddz * ddz;
  }
}
}  // namespace seams

#endif  // SEAMS_HAS_HWY

#endif  // __SIMD_DISTANCE_H_
