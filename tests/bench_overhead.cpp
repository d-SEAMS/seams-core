/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Measures where the spherical-harmonic time goes: the Cartesian to
** angle conversion against the harmonics themselves.
**
** Build target: bench_overhead.  Run: bcpp/tests/bench_overhead
*/

#include <bop.hpp>

#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

template <typename Fn> double bestMillis(Fn &&fn, int reps) {
  double best = std::numeric_limits<double>::max();
  for (int rep = 0; rep < reps; rep++) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    best = std::min(
        best, std::chrono::duration<double, std::milli>(end - start).count());
  }
  return best;
}

} // namespace

int main() {
  // ------------------------------------------------------------------
  // Where the spherical-harmonic time actually goes.
  //
  // The closed forms need sin and cos of the polar angle and the azimuthal
  // phase.  Those come from a Cartesian displacement, so the current path
  // pays acos and atan2 to build angles it then immediately turns back into
  // sines and cosines.  Libraries that evaluate the harmonics directly in
  // Cartesian coordinates skip that round trip entirely, so the question is
  // how much of the cost it represents.
  {
    constexpr size_t nDir = 2000000;
    std::vector<std::array<double, 3>> dirs(nDir);
    unsigned long long state = 99194853094755497ULL;
    auto uniform = [&state]() {
      state ^= state << 13; state ^= state >> 7; state ^= state << 17;
      return static_cast<double>(state >> 11) / 9007199254740992.0 - 0.5;
    };
    for (size_t i = 0; i < nDir; i++) {
      dirs[i] = {uniform(), uniform(), uniform()};
      if (dirs[i][0] == 0.0 && dirs[i][1] == 0.0 && dirs[i][2] == 0.0) {
        dirs[i] = {0.0, 0.0, 1.0};
      }
    }

    // Angles precomputed once, so the harmonics are timed on their own
    std::vector<std::array<double, 2>> angles(nDir);
    for (size_t i = 0; i < nDir; i++) {
      const auto &d = dirs[i];
      const double r = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
      angles[i] = {std::atan2(d[0], d[1]), std::acos(d[2] / r)};
    }

    double sink = 0.0;
    const double tConvert = bestMillis([&]() {
      for (size_t i = 0; i < nDir; i++) {
        const auto &d = dirs[i];
        const double r = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        sink += std::atan2(d[0], d[1]) + std::acos(d[2] / r);
      }
    }, 3);

    const double tHarmonic = bestMillis([&]() {
      for (size_t i = 0; i < nDir; i++) {
        sink += sph::lookupTableQ6Vec(angles[i])[6].real();
      }
    }, 3);

    const double tBoth = bestMillis([&]() {
      for (size_t i = 0; i < nDir; i++) {
        const auto &d = dirs[i];
        const double r = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        const std::array<double, 2> a = {std::atan2(d[0], d[1]),
                                         std::acos(d[2] / r)};
        sink += sph::lookupTableQ6Vec(a)[6].real();
      }
    }, 3);

    std::cout << "\nSpherical harmonics, " << nDir << " directions (l=6)\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::left << std::setw(34) << "Cartesian -> angles only"
              << std::fixed << std::setprecision(2) << tConvert << " ms\n";
    std::cout << std::left << std::setw(34) << "Y_6m from angles only"
              << tHarmonic << " ms\n";
    std::cout << std::left << std::setw(34) << "combined (current path)"
              << tBoth << " ms\n";
    std::cout << std::left << std::setw(34) << "conversion share"
              << std::setprecision(1) << (100.0 * tConvert / tBoth) << " %\n";
    if (sink == 12345.6789) std::cout << "";
  }

  return 0;
}
