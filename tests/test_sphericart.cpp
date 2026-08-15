#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sphericart_ylm.hpp>
#include <steinhardt_device.hpp>

#include <cmath>

TEST_CASE("sphericart_ylm::available matches the compile-time flag",
          "[sphericart]") {
#ifdef SEAMS_HAS_SPHERICART
  REQUIRE(seams::sphericart_ylm::available());
#else
  REQUIRE_FALSE(seams::sphericart_ylm::available());
  double dummy = 0.0;
  REQUIRE(seams::sphericart_ylm::ylmCartesian(6, &dummy, 0, &dummy) == 0);
#endif
}

#ifdef SEAMS_HAS_SPHERICART

TEST_CASE("sphericart real Ylm have the same l=6 modulus as ylmAll",
          "[sphericart]") {
  const double xyz[3] = {0.3, 0.4, 0.8660254037844386};
  double cart[26];
  REQUIRE(seams::sphericart_ylm::ylmCartesian(6, xyz, 1, cart) == 0);
  const double r = std::sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
  const double phi = std::atan2(xyz[0], xyz[1]);
  const double theta = std::acos(xyz[2] / r);
  double closed[26];
  seams::steinhardt::ylmAll(6, theta, phi, closed);
  double mCart = 0.0;
  double mClosed = 0.0;
  for (int m = 0; m < 13; m++) {
    mCart += cart[2 * m] * cart[2 * m] + cart[2 * m + 1] * cart[2 * m + 1];
    mClosed +=
        closed[2 * m] * closed[2 * m] + closed[2 * m + 1] * closed[2 * m + 1];
  }
  REQUIRE_THAT(mCart, Catch::Matchers::WithinRel(mClosed, 1e-6));
}

#endif
