//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#ifndef SEAMS_TRAJECTORY_HPP_
#define SEAMS_TRAJECTORY_HPP_

#include <mol_sys.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <vector>

/**
 * @brief Multi-frame PBC unwrapping, molecular orientation/dipole vectors,
 *        and streaming time-correlation accumulators.
 */
namespace traj {

/**
 * @brief Largest particle or vector count n whose 3 * n doubles one
 * std::vector can hold. Larger counts are refused, so 3 * n never wraps.
 */
inline constexpr std::size_t kMaxTriples =
    static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) /
    (3 * sizeof(double));

/**
 * @brief Largest box count n whose 9 * n doubles fit in a single array.
 * Larger counts are refused, so 9 * n and pointer offsets never overflow.
 */
inline constexpr std::size_t kMaxBoxes =
    static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) /
    (9 * sizeof(double));

/**
 * @brief Largest max_lag whose max_lag + 1 lag bins one std::vector can
 * hold. Larger lags are refused, so max_lag + 1 never wraps to zero.
 */
inline constexpr std::size_t kMaxLag =
    static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) /
        sizeof(double) -
    1;

/**
 * @brief out = a * b; false, leaving out alone, when that overflows size_t.
 */
[[nodiscard]] inline bool checkedMul(std::size_t a, std::size_t b,
                                     std::size_t &out) noexcept {
  if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) {
    return false;
  }
  out = a * b;
  return true;
}

/**
 * @brief Restricted-triclinic periodic simulation cell.
 *
 * Edge matrix H has columns a = (lx, 0, 0), b = (xy, ly, 0), c = (xz, yz, lz)
 * and origin (ox, oy, oz).
 */
struct BoxCell {
  double ox = 0.0;
  double oy = 0.0;
  double oz = 0.0;
  double lx = 0.0;
  double ly = 0.0;
  double lz = 0.0;
  double xy = 0.0;
  double xz = 0.0;
  double yz = 0.0;
  double halfMin = 0.0;
  bool triclinic = false;
  bool ok = false;

  static BoxCell fromBox9(const double box9[9]) noexcept;
  static BoxCell fromCloud(
      const molSys::PointCloud<molSys::Point<double>, double> &cloud) noexcept;
  void toBox9(double out[9]) const noexcept;

  /**
   * @brief Euclidean minimum-image Cartesian displacement for raw delta
   * (dx, dy, dz).
   *
   * Uses minimage's reduced-basis search when built with it, otherwise
   * micDeltaSearch. Independent rounding of fractional coordinates is not
   * the minimum image on a tilted cell.
   */
  [[nodiscard]] std::array<double, 3> micDelta(double dx, double dy,
                                               double dz) const noexcept;

  /**
   * @brief Euclidean minimum image by exhaustive search over the lattice
   * shifts that can beat the fractional image.
   */
  [[nodiscard]] std::array<double, 3> micDeltaSearch(double dx, double dy,
                                                     double dz) const noexcept;
};

/**
 * @brief Stateful streaming unwrapper across periodic frames.
 */
class TrajectoryUnwrapper {
public:
  explicit TrajectoryUnwrapper(std::size_t nParticles = 0);

  /// A count above kMaxTriples leaves the unwrapper empty, and pushFrame
  /// refuses such a count.
  void reset(std::size_t nParticles = 0);

  [[nodiscard]] std::size_t numParticles() const noexcept { return nParticles_; }
  [[nodiscard]] std::size_t framesSeen() const noexcept { return framesSeen_; }

  bool pushFrame(const double *xyz, std::size_t nParticles, const BoxCell &cell,
                 double *unwrappedOut);

private:
  std::size_t nParticles_ = 0;
  std::size_t framesSeen_ = 0;
  std::vector<double> currentUnwrapped_;
};

/**
 * @brief Unwrap a multi-frame coordinate array [nFrames * 3 * nAtoms].
 *
 * @param boxes Either 1 cell (9 doubles) or nFrames cells (9 * nFrames doubles).
 */
bool unwrapPositions(const double *xyzFrames, std::size_t nFrames,
                     std::size_t nAtoms, const double *boxes,
                     std::size_t nBoxes, double *unwrappedOut);

/**
 * @brief Compute intact PBC-assembled per-molecule centres for one frame.
 *
 * Molecules appear in first-seen order of molIds. If masses is nullptr,
 * equal atom weights are used.
 */
bool moleculeCentres(const double *xyz, std::size_t nAtoms, const BoxCell &cell,
                     const int *molIds, const double *masses,
                     std::vector<int> &uniqueMolIds,
                     std::vector<double> &centresOut);

/**
 * @brief Cloud overload using stored atom coordinates and molIDs.
 */
bool moleculeCentres(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const double *masses, std::vector<int> &uniqueMolIds,
    std::vector<double> &centresOut);

/**
 * @brief Compute and unwrap molecule centres across nFrames frames.
 */
bool unwrapMoleculeCentres(const double *xyzFrames, std::size_t nFrames,
                           std::size_t nAtoms, const double *boxes,
                           std::size_t nBoxes, const int *molIds,
                           const double *masses,
                           std::vector<int> &uniqueMolIds,
                           std::vector<double> &unwrappedCentresOut);

/**
 * @brief Per-water orientation vectors (O->H midpoint and H->H).
 */
struct WaterVectorsResult {
  std::vector<int> molIds;
  std::vector<int> oxygenIndices;
  std::vector<double> ohMidVectors;
  std::vector<double> bisectorUnit;
  std::vector<double> hhVectors;
  std::vector<double> hhUnit;
};

/**
 * @brief Extract per-water O->midpoint(H1,H2) and H1->H2 vectors under MIC.
 *
 * Pairs each oxygen (oxygenType) with its two hydrogens (hydrogenType) by
 * shared molID when multi-atom molecule IDs exist, falling back to the two
 * nearest hydrogens within covalentCutoff Angstrom when molIDs do not group O
 * and H atoms.
 */
bool waterOrientationVectors(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    int oxygenType, int hydrogenType, WaterVectorsResult &out,
    double covalentCutoff = 1.35);

/**
 * @brief Generic per-molecule dipole / orientation vectors.
 */
struct DipoleVectorsResult {
  std::vector<int> molIds;
  std::vector<double> dipoles;
  std::vector<double> unitVectors;
  std::vector<double> magnitudes;
};

/**
 * @brief Compute per-molecule dipole or geometric orientation vectors.
 *
 * When charges is non-null, computes sum_a q_a (r_a^intact - R_geom) in e*A
 * (origin-independent for neutral molecules and translation-invariant for ions).
 * When charges is null, computes the mean MIC bond vector from the first atom
 * of each molecule to its remaining atoms (reducing to O->H-midpoint for O-H-H).
 */
bool moleculeDipoleVectors(const double *xyz, std::size_t nAtoms,
                           const BoxCell &cell, const int *molIds,
                           const double *charges, DipoleVectorsResult &out);

bool moleculeDipoleVectors(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const int *molIdsOverride, const double *charges, DipoleVectorsResult &out);

/**
 * @brief Axis-projected orientational order parameters P1 and P2.
 */
struct OrientationOrderResult {
  std::vector<double> cosTheta;
  std::vector<double> p2;
  double meanP1 = 0.0;
  double meanP2 = 0.0;
  std::size_t nValid = 0;
};

bool axisOrderParameters(const double *vectors, std::size_t nVectors,
                         const double axis[3], OrientationOrderResult &out);

/**
 * @brief Linear Einstein fit of total, 1D axial, and 2D lateral MSD.
 */
struct DiffusionFitResult {
  double dTotal = 0.0;
  double dAxial = 0.0;
  double dLateral = 0.0;
  /// NaN unless alphaFitted: the window needs two lags with t > 0 and
  /// MSD_total > 1e-20 for the log-log fit.
  double alpha = std::numeric_limits<double>::quiet_NaN();
  bool alphaFitted = false;
  double r2Total = 0.0;
};

/**
 * @brief Integrated and single-exponential relaxation time summary.
 */
struct DecayFitResult {
  double tauIntegral = 0.0;
  double tauExp = 0.0;
  double amplitude = 1.0;
  double r2 = 0.0;
};

bool fitDiffusionCurve(const std::vector<double> &lagTime,
                       const std::vector<double> &msdTotal,
                       const std::vector<double> &msdAxial,
                       const std::vector<double> &msdLateral, double tStart,
                       double tEnd, DiffusionFitResult &out);

bool fitCorrelationDecay(const std::vector<double> &lagTime,
                         const std::vector<double> &corr, double tStart,
                         double tEnd, DecayFitResult &out);

/**
 * @brief Streaming multi-origin MSD accumulator with arbitrary-axis split.
 *
 * When maxLag > 0, memory is bounded to O(nParticles * maxLag / originStride).
 * When maxLag == 0, all lags up to framesSeen - 1 are retained.
 */
class MsdAccumulator {
public:
  MsdAccumulator(std::size_t nParticles, std::size_t maxLag,
                 std::size_t originStride, double dt, const double axis[3]);

  [[nodiscard]] bool valid() const noexcept { return valid_; }
  [[nodiscard]] std::size_t numParticles() const noexcept { return nParticles_; }
  [[nodiscard]] std::size_t framesSeen() const noexcept { return framesSeen_; }
  [[nodiscard]] std::size_t numLags() const noexcept;

  /**
   * @brief Push one frame of coordinates.
   *
   * If box9OrNull is non-null, coordinates are unwrapped incrementally across
   * consecutive frames before accumulating displacements. If null, xyz is used
   * directly as unwrapped coordinates.
   */
  bool push(const double *xyz, std::size_t nParticles,
            const double *box9OrNull = nullptr);

  void curve(std::vector<double> &lagTime, std::vector<double> &total,
             std::vector<double> &axial, std::vector<double> &lateral,
             std::vector<std::size_t> &nOrigins) const;

  bool fit(double tStart, double tEnd, DiffusionFitResult &out) const;

private:
  struct OriginSlot {
    std::size_t frameIndex = 0;
    std::vector<double> coords;
  };

  std::size_t nParticles_ = 0;
  std::size_t maxLag_ = 0;
  std::size_t originStride_ = 1;
  double dt_ = 1.0;
  std::array<double, 3> axisUnit_{0.0, 0.0, 1.0};
  bool valid_ = false;

  std::size_t framesSeen_ = 0;
  bool usedBox_ = false;
  TrajectoryUnwrapper unwrapper_;
  std::vector<double> scratchUnwrapped_;
  std::deque<OriginSlot> activeOrigins_;

  std::vector<double> sumTotal_;
  std::vector<double> sumAxial_;
  std::vector<double> sumLateral_;
  std::vector<std::size_t> originCounts_;
};

/**
 * @brief Streaming multi-origin H-bond lifetime accumulator (C_HB and S_HB).
 */
class HbondLifetimeAccumulator {
public:
  HbondLifetimeAccumulator(std::size_t maxLag, std::size_t originStride,
                           double dt, bool directed = true);

  [[nodiscard]] bool valid() const noexcept { return valid_; }
  [[nodiscard]] std::size_t framesSeen() const noexcept { return framesSeen_; }
  [[nodiscard]] std::size_t numLags() const noexcept;

  bool push(const int *donors, const int *acceptors, std::size_t nEdges);

  void curve(std::vector<double> &lagTime, std::vector<double> &cIntermittent,
             std::vector<double> &sContinuous,
             std::vector<std::size_t> &nOrigins) const;

  bool fit(double tStart, double tEnd, DecayFitResult &intermittentFit,
           DecayFitResult &continuousFit) const;

private:
  struct OriginSlot {
    std::size_t frameIndex = 0;
    std::vector<std::uint64_t> initialEdges;
    std::vector<std::uint64_t> survivingEdges;
  };

  std::size_t maxLag_ = 0;
  std::size_t originStride_ = 1;
  double dt_ = 1.0;
  bool directed_ = true;
  bool valid_ = false;

  std::size_t framesSeen_ = 0;
  std::deque<OriginSlot> activeOrigins_;

  std::vector<double> intermittentNum_;
  std::vector<double> continuousNum_;
  std::vector<double> denom_;
  std::vector<std::size_t> originCounts_;
};

/**
 * @brief Streaming multi-origin orientational time-correlation accumulator
 *        for C1(t) = <P1(u(0).u(t))> and C2(t) = <P2(u(0).u(t))>.
 */
class ReorientationAccumulator {
public:
  ReorientationAccumulator(std::size_t nVectors, std::size_t maxLag,
                           std::size_t originStride, double dt);

  [[nodiscard]] bool valid() const noexcept { return valid_; }
  [[nodiscard]] std::size_t numVectors() const noexcept { return nVectors_; }
  [[nodiscard]] std::size_t framesSeen() const noexcept { return framesSeen_; }
  [[nodiscard]] std::size_t numLags() const noexcept;

  bool push(const double *vectors, std::size_t nVectors);

  void curve(std::vector<double> &lagTime, std::vector<double> &c1,
             std::vector<double> &c2,
             std::vector<std::size_t> &nOrigins) const;

  bool fit(double tStart, double tEnd, DecayFitResult &c1Fit,
           DecayFitResult &c2Fit) const;

private:
  struct OriginSlot {
    std::size_t frameIndex = 0;
    std::vector<double> unitVectors;
    std::vector<std::uint8_t> mask;
  };

  std::size_t nVectors_ = 0;
  std::size_t maxLag_ = 0;
  std::size_t originStride_ = 1;
  double dt_ = 1.0;
  bool valid_ = false;

  std::size_t framesSeen_ = 0;
  std::deque<OriginSlot> activeOrigins_;

  std::vector<double> sumC1_;
  std::vector<double> sumC2_;
  std::vector<std::size_t> pairCounts_;
  std::vector<std::size_t> originCounts_;
};

} // namespace traj

#endif // SEAMS_TRAJECTORY_HPP_
