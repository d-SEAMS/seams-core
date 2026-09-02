#ifndef SEAMS_SIMD_DISTANCE_H_
#define SEAMS_SIMD_DISTANCE_H_

#include <cmath>
#include <algorithm>
#include <cstddef>
#if defined(__has_include)
#if __has_include(<span>)
#include <span>
#define SEAMS_HAS_STD_SPAN 1
#endif
#endif

#ifdef SEAMS_HAS_HWY

#include "hwy/highway.h"

namespace seams {

namespace hn = hwy::HWY_NAMESPACE;

// Compute squared periodic distances for a batch of atom pairs.
// Orthorhombic lengths only: one independent wrap per axis.
// A LAMMPS dump with tilt (PointCloud.box.size() >= 6) stores bound
// spans in box[0..2], not lx, ly, lz. Call gen::periodicDistSq or
// gen::batchPeriodicDistSq for those clouds.
// dx, dy, dz: coordinate differences (length n)
// bx, by, bz: periodic box dimensions
// out: output squared distances (length n)
inline HWY_ATTR void BatchPeriodicDistSq(const double* HWY_RESTRICT dx,
                                const double* HWY_RESTRICT dy,
                                const double* HWY_RESTRICT dz,
                                double bx, double by, double bz,
                                double* HWY_RESTRICT out, size_t n) {
  const hn::ScalableTag<double> d;
  const size_t N = hn::Lanes(d);

  const auto vbx = hn::Set(d, bx);
  const auto vby = hn::Set(d, by);
  const auto vbz = hn::Set(d, bz);

  // The box is fixed across the batch, so the three divisions of the minimum
  // image convention collapse into one reciprocal each, hoisted out of the
  // loop. Vector division has several times the latency of multiplication.
  const double rbx = 1.0 / bx;
  const double rby = 1.0 / by;
  const double rbz = 1.0 / bz;
  const auto vrbx = hn::Set(d, rbx);
  const auto vrby = hn::Set(d, rby);
  const auto vrbz = hn::Set(d, rbz);

  size_t i = 0;
  for (; i + N <= n; i += N) {
    // Callers pass std::vector<double> scratch (16-byte typical). Load/Store
    // require native vector alignment and fault under AVX2/AVX-512.
    auto vdx = hn::LoadU(d, dx + i);
    auto vdy = hn::LoadU(d, dy + i);
    auto vdz = hn::LoadU(d, dz + i);

    // Absolute values
    vdx = hn::Abs(vdx);
    vdy = hn::Abs(vdy);
    vdz = hn::Abs(vdz);

    // Periodic wrap: dr -= box * round(dr * (1 / box))
    vdx = hn::NegMulAdd(vbx, hn::Round(hn::Mul(vdx, vrbx)), vdx);
    vdy = hn::NegMulAdd(vby, hn::Round(hn::Mul(vdy, vrby)), vdy);
    vdz = hn::NegMulAdd(vbz, hn::Round(hn::Mul(vdz, vrbz)), vdz);

    // r2 = dx*dx + dy*dy + dz*dz
    auto r2 = hn::MulAdd(vdx, vdx, hn::MulAdd(vdy, vdy, hn::Mul(vdz, vdz)));
    hn::StoreU(r2, d, out + i);
  }

  // Scalar remainder
  for (; i < n; i++) {
    double ddx = std::fabs(dx[i]);
    double ddy = std::fabs(dy[i]);
    double ddz = std::fabs(dz[i]);
    ddx -= bx * std::round(ddx * rbx);
    ddy -= by * std::round(ddy * rby);
    ddz -= bz * std::round(ddz * rbz);
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
  // Matches the vectorised path: one reciprocal per axis for the whole batch
  const double rbx = 1.0 / bx;
  const double rby = 1.0 / by;
  const double rbz = 1.0 / bz;
  for (size_t i = 0; i < n; i++) {
    double ddx = std::fabs(dx[i]);
    double ddy = std::fabs(dy[i]);
    double ddz = std::fabs(dz[i]);
    ddx -= bx * std::round(ddx * rbx);
    ddy -= by * std::round(ddy * rby);
    ddz -= bz * std::round(ddz * rbz);
    out[i] = ddx * ddx + ddy * ddy + ddz * ddz;
  }
}
}  // namespace seams

#endif  // SEAMS_HAS_HWY

namespace seams {

/**
 * @brief Squared minimum-image distances for a batch of coordinate
 *  differences, with the extents carried alongside the data.
 * @details The ranges must agree in length; the shortest one bounds the work,
 *  so a mismatched call computes fewer results rather than running off the end
 *  of an array. This overload forwards to whichever kernel the build selected.
 * @param[in] dx,dy,dz Coordinate differences.
 * @param[in] bx,by,bz Periodic box dimensions.
 * @param[out] out Squared distances.
 */
#ifdef SEAMS_HAS_STD_SPAN
inline void BatchPeriodicDistSq(std::span<const double> dx,
                                std::span<const double> dy,
                                std::span<const double> dz, double bx,
                                double by, double bz, std::span<double> out) {
  const size_t n = std::min({dx.size(), dy.size(), dz.size(), out.size()});
  BatchPeriodicDistSq(dx.data(), dy.data(), dz.data(), bx, by, bz, out.data(),
                      n);
}
#endif

}  // namespace seams

#endif  // SEAMS_SIMD_DISTANCE_H_
