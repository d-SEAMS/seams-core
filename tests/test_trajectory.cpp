//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "seams.hpp"
#include <trajectory.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <random>
#include <vector>

namespace {

std::array<double, 3> wrapTriclinic(const std::array<double, 3> &r,
                                    const double box9[9]) {
  const double ox = box9[0], oy = box9[1], oz = box9[2];
  const double lx = box9[3], ly = box9[4], lz = box9[5];
  const double xy = box9[6], xz = box9[7], yz = box9[8];
  double sz = (r[2] - oz) / lz;
  double sy = (r[1] - oy - yz * sz) / ly;
  double sx = (r[0] - ox - xy * sy - xz * sz) / lx;
  sx -= std::floor(sx);
  sy -= std::floor(sy);
  sz -= std::floor(sz);
  return {ox + lx * sx + xy * sy + xz * sz, oy + ly * sy + yz * sz,
          oz + lz * sz};
}

} // namespace

TEST_CASE("Triclinic multi-crossing unwrapping and straddling molecule centres",
          "[trajectory][unwrap]") {
  const double box9[9] = {-2.0, 1.5, -0.5, 11.0, 12.0, 13.0, 3.2, -1.8, 2.4};
  constexpr std::size_t nFrames = 35;
  constexpr std::size_t nMols = 3;
  constexpr std::size_t nAtoms = 3 * nMols;

  const std::vector<int> mols = {10, 10, 10, 20, 20, 20, 30, 30, 30};
  const std::vector<double> masses = {15.9994, 1.00794, 1.00794,
                                      15.9994, 1.00794, 1.00794,
                                      15.9994, 1.00794, 1.00794};
  const double mTot = 15.9994 + 2.0 * 1.00794;

  // Rigid water geometry relative to oxygen.
  const std::array<std::array<double, 3>, 3> rel = {{
      {0.0, 0.0, 0.0},
      {0.75695, 0.0, 0.58588},
      {-0.75695, 0.0, 0.58588},
  }};
  const std::array<double, 3> comShift = {
      (1.00794 * (rel[1][0] + rel[2][0])) / mTot,
      (1.00794 * (rel[1][1] + rel[2][1])) / mTot,
      (1.00794 * (rel[1][2] + rel[2][2])) / mTot,
  };

  std::vector<double> trueCom(nFrames * 3 * nMols, 0.0);
  std::vector<double> wrappedAtoms(nFrames * 3 * nAtoms, 0.0);

  const std::array<std::array<double, 3>, nMols> r0 = {{
      {-1.8, 1.7, -0.3}, // right near box corner so H atoms straddle PBC
      {5.0, 8.0, 11.9},
      {2.0, 3.0, 4.0},
  }};
  const std::array<std::array<double, 3>, nMols> vel = {{
      {1.7, -1.5, 1.9},
      {-2.1, 1.8, -1.6},
      {0.9, 2.2, 1.4},
  }};

  for (std::size_t f = 0; f < nFrames; ++f) {
    const double t = static_cast<double>(f);
    for (std::size_t m = 0; m < nMols; ++m) {
      const std::array<double, 3> oPos = {
          r0[m][0] + vel[m][0] * t,
          r0[m][1] + vel[m][1] * t,
          r0[m][2] + vel[m][2] * t,
      };
      // Frame 0 wrapped anchor determines global lattice image offset.
      const auto oWrap0 = wrapTriclinic(r0[m], box9);
      const std::array<double, 3> latticeShift = {
          oWrap0[0] - r0[m][0],
          oWrap0[1] - r0[m][1],
          oWrap0[2] - r0[m][2],
      };
      trueCom[(f * nMols + m) * 3 + 0] = oPos[0] + comShift[0] + latticeShift[0];
      trueCom[(f * nMols + m) * 3 + 1] = oPos[1] + comShift[1] + latticeShift[1];
      trueCom[(f * nMols + m) * 3 + 2] = oPos[2] + comShift[2] + latticeShift[2];

      for (std::size_t a = 0; a < 3; ++a) {
        const std::array<double, 3> rawAtom = {
            oPos[0] + rel[a][0],
            oPos[1] + rel[a][1],
            oPos[2] + rel[a][2],
        };
        const auto wAtom = wrapTriclinic(rawAtom, box9);
        const std::size_t idx = (f * nAtoms + 3 * m + a) * 3;
        wrappedAtoms[idx + 0] = wAtom[0];
        wrappedAtoms[idx + 1] = wAtom[1];
        wrappedAtoms[idx + 2] = wAtom[2];
      }
    }
  }

  const auto res = seams::unwrap_centres(
      wrappedAtoms, nFrames, nAtoms, std::span<const double>(box9, 9), mols,
      masses);
  REQUIRE(res.mol_ids == std::vector<int>{10, 20, 30});
  REQUIRE(res.centres.size() == trueCom.size());
  for (std::size_t i = 0; i < trueCom.size(); ++i) {
    CHECK(res.centres[i] == Catch::Approx(trueCom[i]).margin(1e-11));
  }
}

namespace {

std::array<double, 3> roundedFractionalImage(const std::array<double, 3> &d,
                                             const double box9[9]) {
  const double lx = box9[3], ly = box9[4], lz = box9[5];
  const double xy = box9[6], xz = box9[7], yz = box9[8];
  const double sz = d[2] / lz;
  const double sy = (d[1] - yz * sz) / ly;
  const double sx = (d[0] - xy * sy - xz * sz) / lx;
  const double tx = sx - std::nearbyint(sx);
  const double ty = sy - std::nearbyint(sy);
  const double tz = sz - std::nearbyint(sz);
  return {lx * tx + xy * ty + xz * tz, ly * ty + yz * tz, lz * tz};
}

std::array<double, 3> bruteForceImage(const std::array<double, 3> &d,
                                      const double box9[9], int range) {
  const double lx = box9[3], ly = box9[4], lz = box9[5];
  const double xy = box9[6], xz = box9[7], yz = box9[8];
  std::array<double, 3> best = d;
  double best2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  for (int i = -range; i <= range; ++i) {
    for (int j = -range; j <= range; ++j) {
      for (int k = -range; k <= range; ++k) {
        const std::array<double, 3> c = {d[0] - lx * i - xy * j - xz * k,
                                         d[1] - ly * j - yz * k,
                                         d[2] - lz * k};
        const double c2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
        if (c2 < best2) {
          best2 = c2;
          best = c;
        }
      }
    }
  }
  return best;
}

double dist3(const std::array<double, 3> &a, const std::array<double, 3> &b) {
  return std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

double norm3(const std::array<double, 3> &a) {
  return std::hypot(a[0], a[1], a[2]);
}

} // namespace

TEST_CASE("Unwrapping takes the Euclidean image on a skewed LAMMPS cell",
          "[trajectory][unwrap][triclinic]") {
  // LAMMPS-legal tilts: |xy| <= lx/2, |xz| <= lx/2, |yz| <= ly/2.
  const double box0[9] = {0.0,         0.0,        0.0,
                          18.03850526, 10.63623724, 11.40838967,
                          -8.77049545, 6.89505882,  2.90161048};
  REQUIRE(std::abs(box0[6]) <= 0.5 * box0[3]);
  REQUIRE(std::abs(box0[7]) <= 0.5 * box0[3]);
  REQUIRE(std::abs(box0[8]) <= 0.5 * box0[4]);

  const std::array<double, 3> rawDelta = {-16.98407757, -4.19565342,
                                          8.11028162};
  const std::array<double, 3> euclid = {3.42737868, 3.53897334, -3.29810805};
  const auto frac = roundedFractionalImage(rawDelta, box0);
  // Rounding each fractional coordinate lands one a-vector away.
  REQUIRE(dist3(frac, euclid) == Catch::Approx(box0[3]).margin(1e-7));
  REQUIRE(norm3(frac) > norm3(euclid) + 9.0);

  const auto cell0 = traj::BoxCell::fromBox9(box0);
  REQUIRE(cell0.ok);
  REQUIRE(cell0.triclinic);
  const auto mic = cell0.micDelta(rawDelta[0], rawDelta[1], rawDelta[2]);
  const auto search =
      cell0.micDeltaSearch(rawDelta[0], rawDelta[1], rawDelta[2]);
  CHECK(dist3({mic[0], mic[1], mic[2]}, euclid) < 1e-7);
  CHECK(dist3({search[0], search[1], search[2]}, euclid) < 1e-7);

  // Frame 2 is the same cell scaled by 3%; the pair scales with it, so the
  // rounded fractional step is again one (scaled) a-vector off.
  constexpr double scale = 1.03;
  double box1[9];
  std::copy_n(box0, 9, box1);
  for (int i = 3; i < 9; ++i) {
    box1[i] *= scale;
  }
  REQUIRE(std::abs(box1[6]) <= 0.5 * box1[3]);
  REQUIRE(std::abs(box1[7]) <= 0.5 * box1[3]);
  REQUIRE(std::abs(box1[8]) <= 0.5 * box1[4]);

  const std::array<double, 3> u0 = {4.1, 2.7, 5.3};
  const std::array<double, 3> u1 = {u0[0] + euclid[0], u0[1] + euclid[1],
                                    u0[2] + euclid[2]};
  const std::array<double, 3> u2 = {u1[0] + scale * euclid[0],
                                    u1[1] + scale * euclid[1],
                                    u1[2] + scale * euclid[2]};
  const double *frameBoxes[3] = {box0, box0, box1};
  const std::array<std::array<double, 3>, 3> truth = {u0, u1, u2};

  std::vector<double> boxes(27);
  std::vector<double> wrapped(9);
  for (int f = 0; f < 3; ++f) {
    std::copy_n(frameBoxes[f], 9, boxes.begin() + 9 * f);
    const auto w = wrapTriclinic(truth[f], frameBoxes[f]);
    std::copy(w.begin(), w.end(), wrapped.begin() + 3 * f);
  }
  const std::array<double, 3> raw2 = {wrapped[6] - u1[0], wrapped[7] - u1[1],
                                      wrapped[8] - u1[2]};
  const std::array<double, 3> step2 = {scale * euclid[0], scale * euclid[1],
                                       scale * euclid[2]};
  REQUIRE(dist3(roundedFractionalImage(raw2, box1), step2) ==
          Catch::Approx(box1[3]).margin(1e-7));

  std::vector<double> out(9, 0.0);
  REQUIRE(traj::unwrapPositions(wrapped.data(), 3, 1, boxes.data(), 3,
                                out.data()));
  // The first frame is already in the cell, so no lattice offset is needed.
  for (int f = 0; f < 3; ++f) {
    for (int k = 0; k < 3; ++k) {
      CHECK(out[3 * f + k] == Catch::Approx(truth[f][k]).margin(1e-9));
    }
  }

  const auto res = seams::unwrap_positions(wrapped, 3, 1, boxes);
  REQUIRE(res.size() == 9);
  for (std::size_t i = 0; i < 9; ++i) {
    CHECK(res[i] == Catch::Approx(out[i]).margin(1e-12));
  }
}

TEST_CASE("Minimum image matches brute force on random LAMMPS-legal cells",
          "[trajectory][unwrap][triclinic]") {
  std::mt19937_64 rng(20261120);
  std::uniform_real_distribution<double> len(1.0, 30.0);
  std::uniform_real_distribution<double> tilt(-0.5, 0.5);
  std::uniform_real_distribution<double> frac(-1.5, 1.5);
  int roundedMisses = 0;
  for (int n = 0; n < 2000; ++n) {
    const double lx = len(rng), ly = len(rng), lz = len(rng);
    const double box9[9] = {0.0,       0.0,       0.0,       lx,       ly,
                            lz,        tilt(rng) * lx, tilt(rng) * lx,
                            tilt(rng) * ly};
    const double a = frac(rng), b = frac(rng), c = frac(rng);
    const std::array<double, 3> d = {lx * a + box9[6] * b + box9[7] * c,
                                     ly * b + box9[8] * c, lz * c};
    const auto cell = traj::BoxCell::fromBox9(box9);
    REQUIRE(cell.ok);
    const double want = norm3(bruteForceImage(d, box9, 16));
    const auto mic = cell.micDelta(d[0], d[1], d[2]);
    const auto search = cell.micDeltaSearch(d[0], d[1], d[2]);
    CHECK(norm3({mic[0], mic[1], mic[2]}) ==
          Catch::Approx(want).margin(1e-9));
    CHECK(norm3({search[0], search[1], search[2]}) ==
          Catch::Approx(want).margin(1e-9));
    if (norm3(roundedFractionalImage(d, box9)) > want + 1e-9) {
      ++roundedMisses;
    }
  }
  // The sample must contain cells where rounding is not the minimum image.
  CHECK(roundedMisses > 100);
}

TEST_CASE("Ballistic and anisotropic MSD with arbitrary axis and negative controls",
          "[trajectory][msd]") {
  const std::array<double, 3> axis = {1.0, 2.0, 2.0}; // norm = 3
  const double dt = 0.25;
  constexpr std::size_t nFrames = 25;

  SECTION("Ballistic motion wrapped into periodic cell matches exact v^2 t^2") {
    const double box9[9] = {0, 0, 0, 9.0, 10.0, 11.0, 1.5, -0.5, 1.0};
    // v = 1.2 * e_hat + 0.8 * e_perp, where e_hat = (1,2,2)/3 and e_perp = (2,-1,0)/sqrt(5)
    const double vPar = 1.2;
    const double vPerp = 0.8;
    const double inv5 = 1.0 / std::sqrt(5.0);
    const std::array<double, 3> vel = {
        vPar * (1.0 / 3.0) + vPerp * (2.0 * inv5),
        vPar * (2.0 / 3.0) - vPerp * (1.0 * inv5),
        vPar * (2.0 / 3.0),
    };

    seams::Msd msdAll(1, 0, 1, dt, axis);
    seams::Msd msdBounded(1, 10, 1, dt, axis);
    for (std::size_t f = 0; f < nFrames; ++f) {
      const double t = static_cast<double>(f) * dt;
      const std::array<double, 3> pos = {
          1.0 + vel[0] * t, 2.0 + vel[1] * t, 3.0 + vel[2] * t};
      const auto w = wrapTriclinic(pos, box9);
      msdAll.push(w, std::span<const double, 9>(box9, 9));
      msdBounded.push(w, std::span<const double, 9>(box9, 9));
    }

    const auto cAll = msdAll.curve();
    const auto cBnd = msdBounded.curve();
    REQUIRE(cAll.lag_time.size() == nFrames);
    REQUIRE(cBnd.lag_time.size() == 11);

    for (std::size_t k = 0; k < cAll.lag_time.size(); ++k) {
      const double t = static_cast<double>(k) * dt;
      CHECK(cAll.axial[k] == Catch::Approx(vPar * vPar * t * t).margin(1e-11));
      CHECK(cAll.lateral[k] == Catch::Approx(vPerp * vPerp * t * t).margin(1e-11));
      CHECK(cAll.total[k] ==
            Catch::Approx((vPar * vPar + vPerp * vPerp) * t * t).margin(1e-11));
      if (k < cBnd.lag_time.size()) {
        CHECK(cBnd.total[k] == Catch::Approx(cAll.total[k]).margin(1e-12));
      }
    }

    const auto fit = msdAll.fit(0.25, 4.0);
    CHECK(fit.alpha_fitted);
    CHECK(fit.alpha == Catch::Approx(2.0).margin(1e-10));
  }

  SECTION("A skipped log-log fit reports NaN alpha, not zero") {
    const std::array<double, 3> xAxis = {1.0, 0.0, 0.0};
    // One origin: lag 1 moves 1e-11 (MSD 1e-22 < 1e-20), lag 2 moves 2.
    seams::Msd msd(1, 0, 10, 1.0, xAxis);
    for (const double x : {0.0, 1e-11, 2.0}) {
      msd.push(std::array<double, 3>{x, 0.0, 0.0});
    }
    const auto fit = msd.fit(0.0, 2.0);
    CHECK_FALSE(fit.alpha_fitted);
    CHECK(std::isnan(fit.alpha));
    CHECK(fit.d_total == Catch::Approx(1.0 / 3.0).margin(1e-12));
    CHECK(fit.d_axial == Catch::Approx(1.0).margin(1e-12));

    traj::DiffusionFitResult res;
    REQUIRE(traj::fitDiffusionCurve({0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0},
                                    {0.0, 0.0}, 0.0, 1.0, res));
    CHECK_FALSE(res.alphaFitted);
    CHECK(std::isnan(res.alpha));
    CHECK(res.dTotal == Catch::Approx(1.0 / 6.0).margin(1e-14));
  }

  SECTION("Negative controls: stationary and purely lateral motion") {
    seams::Msd msdStatic(2, 8, 1, 1.0, axis);
    const std::vector<double> fixed = {1.0, 2.0, 3.0, -1.0, 0.5, 4.0};
    for (int f = 0; f < 6; ++f) {
      msdStatic.push(fixed);
    }
    const auto cs = msdStatic.curve();
    for (double val : cs.total) {
      CHECK(val == Catch::Approx(0.0).margin(1e-15));
    }

    // Motion strictly perpendicular to axis (1, 2, 2): direction (2, -1, 0).
    seams::Msd msdLat(1, 8, 1, 1.0, axis);
    for (int f = 0; f < 6; ++f) {
      const std::array<double, 3> p = {2.0 * f, -1.0 * f, 0.0};
      msdLat.push(p);
    }
    const auto cl = msdLat.curve();
    for (std::size_t k = 1; k < cl.lag_time.size(); ++k) {
      CHECK(cl.axial[k] == Catch::Approx(0.0).margin(1e-24));
      CHECK(cl.lateral[k] == Catch::Approx(cl.total[k]).margin(1e-12));
    }
  }
}

TEST_CASE("Water vectors, generic dipoles, and axis order parameters",
          "[trajectory][dipole]") {
  const std::array<double, 9> box9 = {0, 0, 0, 10.0, 10.0, 10.0, 0, 0, 0};
  // Place O near z = 9.9 and H1, H2 wrapped across z = 0 boundary (+z bisector).
  const double dOhZ = 0.58588;
  const double dOhX = 0.75695;
  const std::vector<double> xyz = {
      5.0,          5.0, 9.9,                 // O
      5.0 + dOhX,   5.0, 9.9 + dOhZ - 10.0,   // H1 wrapped
      5.0 - dOhX,   5.0, 9.9 + dOhZ - 10.0,   // H2 wrapped
  };
  const std::vector<int> types = {2, 1, 1};
  const std::vector<int> ids = {1, 2, 3};
  const std::vector<int> mols = {7, 7, 7};

  seams::System sys;
  sys.set_frame(xyz, box9, types, ids);
  sys.set_mols(mols);

  const auto wv = sys.water_vectors(2, 1);
  REQUIRE(wv.mol_ids == std::vector<int>{7});
  CHECK(wv.oh_mid[0] == Catch::Approx(0.0).margin(1e-12));
  CHECK(wv.oh_mid[1] == Catch::Approx(0.0).margin(1e-12));
  CHECK(wv.oh_mid[2] == Catch::Approx(dOhZ).margin(1e-12));
  CHECK(wv.bisector[2] == Catch::Approx(1.0).margin(1e-12));
  // H1 -> H2 points along -x.
  CHECK(wv.hh[0] == Catch::Approx(-2.0 * dOhX).margin(1e-12));
  CHECK(wv.hh_unit[0] == Catch::Approx(-1.0).margin(1e-12));

  // Charge-weighted SPC/E dipole matches bisector direction and exact magnitude.
  const std::vector<double> qSpce = {-0.8476, 0.4238, 0.4238};
  const auto dipQ = sys.molecule_dipoles(mols, qSpce);
  CHECK(dipQ.dipoles[2] == Catch::Approx(2.0 * 0.4238 * dOhZ).margin(1e-12));
  CHECK(dipQ.unit_vectors[2] == Catch::Approx(1.0).margin(1e-12));

  // Charge-free generic orientation matches O -> H-midpoint!
  const auto dipGeom = sys.molecule_dipoles();
  CHECK(dipGeom.dipoles[2] == Catch::Approx(dOhZ).margin(1e-12));
  CHECK(dipGeom.unit_vectors[2] == Catch::Approx(1.0).margin(1e-12));

  // Axis order parameters: z-axis gives P1=1, P2=1 for bisector and P1=0, P2=-0.5 for HH.
  const std::array<double, 3> zAxis = {0.0, 0.0, 1.0};
  const auto ordBis = seams::orientation_order(wv.bisector, zAxis);
  CHECK(ordBis.mean_p1 == Catch::Approx(1.0).margin(1e-12));
  CHECK(ordBis.mean_p2 == Catch::Approx(1.0).margin(1e-12));

  const auto ordHh = seams::orientation_order(wv.hh_unit, zAxis);
  CHECK(ordHh.mean_p1 == Catch::Approx(0.0).margin(1e-12));
  CHECK(ordHh.mean_p2 == Catch::Approx(-0.5).margin(1e-12));

  // Negative control: symmetric octahedral 6-vector set has <P1> = 0 and <P2> = 0.
  const std::vector<double> octa = {
      1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1};
  const auto ordIso = seams::orientation_order(octa, zAxis);
  CHECK(ordIso.mean_p1 == Catch::Approx(0.0).margin(1e-14));
  CHECK(ordIso.mean_p2 == Catch::Approx(0.0).margin(1e-14));
}

TEST_CASE("Continuous vs intermittent H-bond lifetime and rigid-rotor reorientation",
          "[trajectory][tcf]") {
  SECTION("Exact blinking H-bond negative control: S_HB dies while C_HB recovers") {
    seams::HbondTcf tcf(0, 1, 1.0, true);
    const std::vector<int> d1 = {1}, a1 = {2};
    const std::vector<int> empty;
    // Frame 0: {(1,2)}, Frame 1: {}, Frame 2: {(1,2)}
    tcf.push(d1, a1);
    tcf.push(empty, empty);
    tcf.push(d1, a1);
    const auto c = tcf.curve();
    REQUIRE(c.lag_time.size() == 3);
    CHECK(c.c_intermittent[0] == Catch::Approx(1.0));
    CHECK(c.s_continuous[0] == Catch::Approx(1.0));
    // Lag 1: origins 0 ({(1,2)} -> {}) and 1 ({} -> {(1,2)}), denom = 1 + 0 = 1, num = 0.
    CHECK(c.c_intermittent[1] == Catch::Approx(0.0));
    CHECK(c.s_continuous[1] == Catch::Approx(0.0));
    // Lag 2: origin 0 ({(1,2)} -> {} -> {(1,2)}): intermittent = 1, continuous = 0!
    CHECK(c.c_intermittent[2] == Catch::Approx(1.0));
    CHECK(c.s_continuous[2] == Catch::Approx(0.0));
  }

  SECTION("Rigid rotor exact C1(t) = cos(w t) and C2(t) = 0.5*(3 cos^2(w t) - 1)") {
    const double dt = 0.05;
    const double omega = 1.3;
    constexpr std::size_t nFrames = 30;
    seams::Reorientation rot(2, 0, 1, dt);
    for (std::size_t f = 0; f < nFrames; ++f) {
      const double theta = omega * static_cast<double>(f) * dt;
      const std::vector<double> vecs = {
          std::cos(theta), std::sin(theta), 0.0,
          0.0, std::cos(theta), std::sin(theta),
      };
      rot.push(vecs);
    }
    const auto rc = rot.curve();
    REQUIRE(rc.lag_time.size() == nFrames);
    for (std::size_t k = 0; k < nFrames; ++k) {
      const double wt = omega * static_cast<double>(k) * dt;
      const double wantC1 = std::cos(wt);
      const double wantC2 = 0.5 * (3.0 * wantC1 * wantC1 - 1.0);
      CHECK(rc.c1[k] == Catch::Approx(wantC1).margin(1e-12));
      CHECK(rc.c2[k] == Catch::Approx(wantC2).margin(1e-12));
    }
  }
}

TEST_CASE("Accumulators and unwrappers refuse counts whose storage size wraps",
          "[traj][bounds]") {
  constexpr std::size_t most = std::numeric_limits<std::size_t>::max();
  // 3 * (most / 3 + 1) is 2 modulo 2^64, and most + 1 is 0.
  const std::size_t wraps3 = most / 3 + 1;
  const double axis[3] = {0.0, 0.0, 1.0};
  const double box9[9] = {0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0};

  CHECK_FALSE(traj::MsdAccumulator(wraps3, 4, 1, 1.0, axis).valid());
  CHECK_FALSE(traj::MsdAccumulator(traj::kMaxTriples + 1, 4, 1, 1.0, axis).valid());
  CHECK_FALSE(traj::MsdAccumulator(1, most, 1, 1.0, axis).valid());
  CHECK_FALSE(traj::MsdAccumulator(1, most - 1, 1, 1.0, axis).valid());
  CHECK_FALSE(traj::MsdAccumulator(1, traj::kMaxLag + 1, 1, 1.0, axis).valid());
  CHECK_FALSE(traj::ReorientationAccumulator(wraps3, 4, 1, 1.0).valid());
  CHECK_FALSE(traj::ReorientationAccumulator(1, most, 1, 1.0).valid());
  CHECK_FALSE(traj::HbondLifetimeAccumulator(most, 1, 1.0, true).valid());
  CHECK_FALSE(traj::HbondLifetimeAccumulator(traj::kMaxLag + 1, 1, 1.0, true).valid());

  // A refused accumulator stays refused: push neither reads nor writes.
  traj::MsdAccumulator msd(1, most, 1, 1.0, axis);
  const double one[3] = {1.0, 2.0, 3.0};
  CHECK_FALSE(msd.push(one, 1, nullptr));
  CHECK(msd.numLags() == 0);
  traj::HbondLifetimeAccumulator hb(most, 1, 1.0, true);
  const int d[1] = {0};
  const int a[1] = {1};
  CHECK_FALSE(hb.push(d, a, 1));
  traj::ReorientationAccumulator reo(1, most, 1, 1.0);
  CHECK_FALSE(reo.push(one, 1));

  traj::TrajectoryUnwrapper unwrapper(wraps3);
  CHECK(unwrapper.numParticles() == 0);
  const traj::BoxCell cell = traj::BoxCell::fromBox9(box9);
  double out[3] = {0.0, 0.0, 0.0};
  CHECK_FALSE(unwrapper.pushFrame(one, wraps3, cell, out));

  std::vector<double> xyz(6, 1.0);
  std::vector<double> unwrapped(6, 0.0);
  CHECK_FALSE(traj::unwrapPositions(xyz.data(), most / 3 + 2, 1, box9, 1, unwrapped.data()));
  CHECK_FALSE(traj::unwrapPositions(xyz.data(), 1, wraps3, box9, 1, unwrapped.data()));
  const int mols[2] = {1, 1};
  std::vector<int> ids;
  std::vector<double> centres;
  CHECK_FALSE(traj::unwrapMoleculeCentres(xyz.data(), most / 3 + 2, 1, box9, 1, mols, nullptr,
                                          ids, centres));

  // max_lag 0 keeps every lag and needs no lag bins up front.
  CHECK(traj::MsdAccumulator(1, 0, 1, 1.0, axis).valid());
  CHECK(traj::HbondLifetimeAccumulator(0, 1, 1.0, true).valid());
}
