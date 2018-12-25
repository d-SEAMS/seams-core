#ifndef _SPHERICAL_HARMONICS_H_
#define _SPHERICAL_HARMONICS_H_

#include <array>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cmath>
#include <complex>

namespace trans {

// 7 is for Q3, orderL=3

std::vector<std::complex<double>> spheriHarmo(int orderL,
                                              std::array<double, 2> coordSph);

std::array<double, 2> radialCoord(std::array<double, 3>);

} // namespace trans

#endif // _SPHERICAL_HARMONICS_H_
