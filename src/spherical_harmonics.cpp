#include "spherical_harmonics.h"

namespace bg = boost::geometry;

blaze::StaticVector<std::complex<double>, 7UL>
trans::spheriHarmo(int orderL, blaze::StaticVector<double, 2UL> radialCoord) {
  // Iterate over the index of order
  std::vector<int> v = {-3, -2, -1, 0, 1, 2, 3};
  blaze::StaticVector<std::complex<double>, 7UL> result;
  for (auto n : v) {
    auto theta = radialCoord[1];
    auto phi = radialCoord[0];
    // This is for l=3
    std::complex<double> b = boost::math::spherical_harmonic(3, n, theta, phi);
    // std::complex<double> b2 = boost::math::spherical_harmonic(3, n, 0.524, 0);
    // std::complex<double> b1 = boost::math::spherical_harmonic(3, n, 1.047, 0);
    result[n] = b;
  }
  return result;
}

std::array<double, 2> trans::radialCoord(std::array<double, 3> cartCoord) {
  // The output
  std::array<double, 2> result;
  // Point Definitions
  bg::model::point<long double, 3, bg::cs::cartesian> cartesianPoint;
  bg::model::point<long double, 3, bg::cs::spherical<bg::radian>> azuPoint;
  // Set Value (TODO: Recurse this)
  bg::set<0>(cartesianPoint, cartCoord[0]);
  bg::set<1>(cartesianPoint, cartCoord[1]);
  bg::set<2>(cartesianPoint, cartCoord[2]);
  // Transform
  bg::transform(cartesianPoint, azuPoint);
  result[0] = bg::get<0>(azuPoint);
  result[1] = bg::get<1>(azuPoint);
  return result;
}
