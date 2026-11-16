// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//
// A consumer of libyodaLib linked without --allow-shlib-undefined. When
// the build took a sphericart library, every symbol libyodaLib needs from
// it must resolve here, and the Ylm call must agree with the closed form.

#include <cmath>
#include <cstdio>
#include <sphericart_ylm.hpp>
#include <vector>

int main() {
  if (!seams::sphericart_ylm::available()) {
    std::puts("sphericart not linked");
    return 77;
  }
  // Y_6^0 on the z axis is sqrt(13/(4 pi)).
  const double xyz[3] = {0.0, 0.0, 1.0};
  std::vector<double> out(2 * 13, 0.0);
  if (seams::sphericart_ylm::ylmCartesian(6, xyz, 1, out.data()) != 0) {
    std::puts("ylmCartesian failed");
    return 1;
  }
  const double want = std::sqrt(13.0 / (4.0 * M_PI));
  const double got = out[2 * 6];
  std::printf("Y60 %.12f expected %.12f\n", got, want);
  return std::fabs(got - want) < 1e-10 ? 0 : 1;
}
