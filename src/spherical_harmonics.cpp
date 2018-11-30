#include "spherical_harmonics.h"

// Temporary
#include <cmath>
#include <iostream>
namespace bg = boost::geometry;

blaze::StaticVector<std::complex<double>, 7UL>
spheriHarmo(int orderL, blaze::StaticVector<double, 2UL> radialCoord) {
  // Iterate over the index of order
  std::vector<int> v = {-3, -2, -1, 0, 1, 2, 3};
  blaze::StaticVector<std::complex<double>, 7UL> result;
  for (auto n : v) {
    std::complex<double> b2 = boost::math::spherical_harmonic(3, n, 0.524, 0);
    std::complex<double> b1 = boost::math::spherical_harmonic(3, n, 1.047, 0);
    std::cout << "Boost " << b1 << " at " << n << std::endl;
    std::cout << "Boost2 " << b2 << " at " << n << std::endl;
    result[n] = b1;
  }
  for (auto n : v) {
    std::cout << result[n] << std::endl;
  }
  return result;
}

// std::array<double, 2> radialCoord(std::array<double, 3> cartCoord) {
//   std::cout << cartCoord[0];
//   std::cout << cartCoord[1];
//   std::cout << cartCoord[2];
//   // The output
//   std::array<double, 2> result;
//   // Point Definitions
//   bg::model::point<long double, 3, bg::cs::cartesian> cartesianPoint;
//   bg::model::point<long double, 2, bg::cs::spherical<bg::degree>> azuPoint;
//   // Set Value (TODO: Recurse this)
//   bg::set<0>(cartesianPoint, cartCoord[0]);
//   bg::set<1>(cartesianPoint, cartCoord[1]);
//   bg::set<2>(cartesianPoint, cartCoord[2]);

//   // Transform
//   bg::transform(cartesianPoint, azuPoint);
//   result[0] = bg::get<0>(azuPoint);
//   result[1] = bg::get<1>(azuPoint);
//   // Test
//   std::cout << "x=" << bg::get<0>(cartesianPoint) << " should be "
//             << cartCoord[0] << " is now " << result[0] << std::endl
//             << " y=" << bg::get<1>(cartesianPoint) << " should be "
//             << cartCoord[1] << " is now " << result[1] << std::endl
//             << " z=" << bg::get<2>(cartesianPoint) << " should be "
//             << cartCoord[2] << std::endl;
//   return result;
// }
