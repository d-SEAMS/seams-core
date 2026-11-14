//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <chrono>

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <domain.hpp>
#include <neighbours.hpp>
#include <thread_team.hpp>
#ifdef SEAMS_HAS_LINKCELL
#include <linkcell.h>
#endif
#include <ozaki.hpp>
#include <ylm_batch.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#define SEAMS_BOP_X86 1
#endif
#include <iostream>
#include <mutex>
#include <numbers>
#include <unordered_map>
#include <utility>
#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

/**
 * @brief Powers of @f$\sin\theta@f$ and @f$\cos\theta@f$ together with the
 *  unit-modulus phases @f$e^{im\phi}@f$, evaluated once per direction.
 * @details The closed forms for @f{3m}@f$, @f{4m}@f$ and @f{6m}@f$ need
 *  @f$\sin^k\theta@f$ and @f$\cos^k\theta@f$ up to @f$k=6@f$ and the phases for
 *  @f$|m| \leq 6@f$. Building them by recurrence costs two calls to the
 *  trigonometric library per direction, shared across every @f$m@f$.
 */
struct AngularTerms {
  std::array<double, 9> sinPow;             //! sinPow[k] = sin(theta)^k
  std::array<double, 9> cosPow;             //! cosPow[k] = cos(theta)^k
  std::array<std::complex<double>, 9> phase; //! phase[k] = exp(i k phi)

  // maxPow is the highest power a caller reads. l = 3 reads 3, l = 6
  // reads 6, l = 8 reads 8. Stopping there leaves the lower powers
  // identical to filling the whole table.
  AngularTerms(double theta, double phi, int maxPow = 8) {
    const double sinTheta = std::sin(theta);
    const double cosTheta = std::cos(theta);
    const std::complex<double> unitPhase(std::cos(phi), std::sin(phi));

    sinPow[0] = 1.0;
    cosPow[0] = 1.0;
    phase[0] = {1.0, 0.0};
    const int last = std::min(maxPow, static_cast<int>(sinPow.size()) - 1);
    for (int k = 1; k <= last; k++) {
      sinPow[static_cast<std::size_t>(k)] =
          sinPow[static_cast<std::size_t>(k - 1)] * sinTheta;
      cosPow[static_cast<std::size_t>(k)] =
          cosPow[static_cast<std::size_t>(k - 1)] * cosTheta;
      phase[static_cast<std::size_t>(k)] =
          phase[static_cast<std::size_t>(k - 1)] * unitPhase;
    }
  }
};

/**
 * @brief Normalisation constants @f$N_{lm}@f$ of the Condon-Shortley spherical
 *  harmonics, for @f$m \geq 0@f$.
 * @details Held as function-local statics so that the square roots are
 *  evaluated once for the lifetime of the process rather than per call.
 */
const std::array<double, 4> &normQ3() {
  static const std::array<double, 4> values = {
      0.25 * std::sqrt(7.0 / std::numbers::pi),           // m = 0
      0.125 * std::sqrt(21.0 / std::numbers::pi),         // |m| = 1
      0.25 * std::sqrt(105.0 / (2.0 * std::numbers::pi)), // |m| = 2
      0.125 * std::sqrt(35.0 / std::numbers::pi)};        // |m| = 3
  return values;
}

const std::array<double, 5> &normQ4() {
  static const std::array<double, 5> values = {
      0.1875 * std::sqrt(1.0 / std::numbers::pi),          // m = 0
      0.375 * std::sqrt(5.0 / std::numbers::pi),           // |m| = 1
      0.375 * std::sqrt(5.0 / (2.0 * std::numbers::pi)),   // |m| = 2
      0.375 * std::sqrt(35.0 / std::numbers::pi),          // |m| = 3
      0.1875 * std::sqrt(35.0 / (2.0 * std::numbers::pi))}; // |m| = 4
  return values;
}

const std::array<double, 9> &normQ8() {
  static const std::array<double, 9> values = {
      std::sqrt(17.0 / (4.0 * std::numbers::pi)),              // m = 0
      std::sqrt(17.0 / (288.0 * std::numbers::pi)),            // |m| = 1
      std::sqrt(17.0 / (20160.0 * std::numbers::pi)),          // |m| = 2
      std::sqrt(17.0 / (1330560.0 * std::numbers::pi)),        // |m| = 3
      std::sqrt(17.0 / (79833600.0 * std::numbers::pi)),       // |m| = 4
      std::sqrt(17.0 / (4151347200.0 * std::numbers::pi)),     // |m| = 5
      std::sqrt(17.0 / (174356582400.0 * std::numbers::pi)),   // |m| = 6
      std::sqrt(17.0 / (5230697472000.0 * std::numbers::pi)),  // |m| = 7
      std::sqrt(17.0 / (83691159552000.0 * std::numbers::pi))}; // |m| = 8
  return values;
}

const std::array<double, 7> &normQ6() {
  static const std::array<double, 7> values = {
      0.03125 * std::sqrt(13.0 / std::numbers::pi),         // m = 0
      0.0625 * std::sqrt(273.0 / (2.0 * std::numbers::pi)), // |m| = 1
      0.015625 * std::sqrt(1365.0 / std::numbers::pi),      // |m| = 2
      0.03125 * std::sqrt(1365.0 / std::numbers::pi),       // |m| = 3
      0.09375 * std::sqrt(91.0 / (2.0 * std::numbers::pi)), // |m| = 4
      0.09375 * std::sqrt(1001.0 / std::numbers::pi),       // |m| = 5
      0.015625 * std::sqrt(3003.0 / std::numbers::pi)};     // |m| = 6
  return values;
}

/**
 * @brief The associated-Legendre amplitude of @f$Y_{lm}@f$ stripped of its
 *  normalisation and azimuthal phase.
 * @param[in] orderL Degree @f@f$; 3, 4 or 6.
 * @param[in] absM Absolute order @f$|m|@f$.
 * @param[in] terms Precomputed powers of the polar angle.
 * @return The real amplitude @f$P_{l|m|}(\cos\theta)@f$ in the convention used
 *  by the hard-coded tables.
 */
double legendreAmplitude(int orderL, int absM, const AngularTerms &terms) {
  const auto &s = terms.sinPow;
  const auto &c = terms.cosPow;

  if (orderL == 3) {
    if (absM <= 0) {
      return 5.0 * c[3] - 3.0 * c[1];
    }
    if (absM == 1) {
      return s[1] * (5.0 * c[2] - 1.0);
    }
    if (absM == 2) {
      return s[2] * c[1];
    }
    // Callers pass absM from 0 through orderL, so 3 is the remaining value.
    return s[3];
  }

  if (orderL == 4) {
    if (absM <= 0) {
      return 35.0 * c[4] - 30.0 * c[2] + 3.0;
    }
    if (absM == 1) {
      return s[1] * (7.0 * c[3] - 3.0 * c[1]);
    }
    if (absM == 2) {
      return s[2] * (7.0 * c[2] - 1.0);
    }
    if (absM == 3) {
      return s[3] * c[1];
    }
    return s[4];
  }

  if (orderL == 6) {
    if (absM <= 0) {
      return 231.0 * c[6] - 315.0 * c[4] + 105.0 * c[2] - 5.0;
    }
    if (absM == 1) {
      return s[1] * (33.0 * c[5] - 30.0 * c[3] + 5.0 * c[1]);
    }
    if (absM == 2) {
      return s[2] * (33.0 * c[4] - 18.0 * c[2] + 1.0);
    }
    if (absM == 3) {
      return s[3] * (11.0 * c[3] - 3.0 * c[1]);
    }
    if (absM == 4) {
      return s[4] * (11.0 * c[2] - 1.0);
    }
    if (absM == 5) {
      return s[5] * c[1];
    }
    return s[6];
  }

  if (absM <= 0) {
    return (35.0 + c[2] * (-1260.0 + c[2] * (6930.0 + c[2] * (-12012.0 + c[2] * 6435.0)))) / 128.0;
  }
  if (absM == 1) {
    return s[1] * c[1] * (-315.0 + c[2] * (3465.0 + c[2] * (-9009.0 + c[2] * 6435.0))) / 16.0;
  }
  if (absM == 2) {
    return s[2] * (-315.0 + c[2] * (10395.0 + c[2] * (-45045.0 + c[2] * 45045.0))) / 16.0;
  }
  if (absM == 3) {
    return s[3] * c[1] * (10395.0 + c[2] * (-90090.0 + c[2] * 135135.0)) / 8.0;
  }
  if (absM == 4) {
    return s[4] * (10395.0 + c[2] * (-270270.0 + c[2] * 675675.0)) / 8.0;
  }
  if (absM == 5) {
    return s[5] * c[1] * (-135135.0 + c[2] * 675675.0) / 2.0;
  }
  if (absM == 6) {
    return s[6] * (-135135.0 + c[2] * 2027025.0) / 2.0;
  }
  if (absM == 7) {
    return s[7] * c[1] * 2027025.0;
  }
  return s[8] * 2027025.0;
}

/**
 * @brief Evaluates the conjugate pair @f$(Y_{l,-m}, Y_{l,+m})@f$ for
 *  @f$m \geq 0@f$.
 * @details The two members satisfy
 *  @f$Y_{l,m} = (-1)^m Y_{l,-m}^{*}@f$ exactly, since both are formed from the
 *  same amplitude and the same phase.
 * @param[in] orderL Degree @f@f$; 3, 4 or 6.
 * @param[in] absM Absolute order @f$|m|@f$.
 * @param[in] terms Precomputed powers and phases.
 * @return A pair holding @f$Y_{l,-m}@f$ then @f$Y_{l,+m}@f$.
 */
std::pair<std::complex<double>, std::complex<double>>
harmonicPair(int orderL, int absM, const AngularTerms &terms) {
  const double norm = (orderL == 3)   ? normQ3()[absM]
                      : (orderL == 4) ? normQ4()[absM]
                      : (orderL == 8) ? normQ8()[absM]
                                      : normQ6()[absM];
  const double amplitude = norm * legendreAmplitude(orderL, absM, terms);

  const std::complex<double> negative = amplitude * std::conj(terms.phase[absM]);
  const std::complex<double> positive =
      ((absM % 2 == 0) ? amplitude : -amplitude) * terms.phase[absM];

  return {negative, positive};
}

// Real l = 3 harmonics, one row per bond. Channels are
// (Y_30, A_1 cos phi, A_1 sin phi, A_2 cos 2phi, A_2 sin 2phi,
//  A_3 cos 3phi, A_3 sin 3phi), with A_m the normalised associated
//  Legendre value. The correlation multiplies every |m| > 0 channel
//  by two, which is the complex inner product of the conjugate pair.
//  Bonds are evaluated in structure-of-arrays groups so the polynomial
//  is independent across the group and contracts with fma.
constexpr int kRealL3 = 7;
// An odd harmonic of a centrosymmetric shell cancels. The leftover is
// roundoff, and its correlation is not a bond angle. Exact cancellation
// already returns 0; this floor does the same for a shell that only
// cancels to machine noise, so a translation cannot change the label.
constexpr double kCancelledShell = 1e-12;

double correlationOf(double dot, double normI, double normJ) {
  if (!(normI > kCancelledShell && normJ > kCancelledShell)) {
    return 0.0;
  }
  const double denom = std::sqrt(normI * normJ);
  return denom > 0.0 ? dot / denom : 0.0;
}

void addRealYlm3Scalar(double x, double y, double z, const std::array<double, 4> &norm, double *acc) {
  // Unnormalised solid harmonic. A zero bond is the north pole, which is
  // what the trigonometric form produced when r == 0.
  const double x2 = x * x;
  const double y2 = y * y;
  const double z2 = z * z;
  const double r2 = x2 + y2 + z2;
  if (!(r2 > 0.0)) {
    acc[0] += norm[0] * 2.0;
    return;
  }
  const double inv = 1.0 / (r2 * std::sqrt(r2));
  const double t1 = 5.0 * z2 - r2;
  acc[0] += norm[0] * (z * (5.0 * z2 - 3.0 * r2)) * inv;
  acc[1] += norm[1] * (y * t1) * inv;
  acc[2] += norm[1] * (x * t1) * inv;
  acc[3] += norm[2] * (z * (y2 - x2)) * inv;
  acc[4] += norm[2] * (2.0 * x * y * z) * inv;
  acc[5] += norm[3] * (y * (y2 - 3.0 * x2)) * inv;
  acc[6] += norm[3] * (x * (3.0 * y2 - x2)) * inv;
}

#ifdef SEAMS_BOP_X86
// Four bonds side by side. Hardware fma is a per-lane instruction; the
// channels are then added in bond order so the sum matches the scalar loop.
__attribute__((target("avx2,fma"))) void realYlm3Accumulate(const double *x, const double *y,
                                                           const double *z, int nBond, double *acc) {
  const auto &norm = normQ3();
  int bond = 0;
  const __m256d zero = _mm256_setzero_pd();
  const __m256d one = _mm256_set1_pd(1.0);
  const __m256d five = _mm256_set1_pd(5.0);
  const __m256d n0 = _mm256_set1_pd(norm[0]);
  const __m256d n1 = _mm256_set1_pd(norm[1]);
  const __m256d n2 = _mm256_set1_pd(norm[2]);
  const __m256d n3 = _mm256_set1_pd(norm[3]);
  const __m256d three = _mm256_set1_pd(3.0);
  const __m256d two = _mm256_set1_pd(2.0);
  for (; bond + 4 <= nBond; bond += 4) {
    const __m256d X = _mm256_loadu_pd(x + bond);
    const __m256d Y = _mm256_loadu_pd(y + bond);
    const __m256d Z = _mm256_loadu_pd(z + bond);
    const __m256d x2 = _mm256_mul_pd(X, X);
    const __m256d y2 = _mm256_mul_pd(Y, Y);
    const __m256d z2 = _mm256_mul_pd(Z, Z);
    const __m256d r2 = _mm256_add_pd(x2, _mm256_add_pd(y2, z2));
    const __m256d live = _mm256_cmp_pd(r2, zero, _CMP_GT_OQ);
    const __m256d inv = _mm256_blendv_pd(
        zero, _mm256_div_pd(one, _mm256_mul_pd(r2, _mm256_sqrt_pd(r2))), live);
    const __m256d t1 = _mm256_fmsub_pd(five, z2, r2);
    const __m256d q0 = _mm256_mul_pd(
        n0, _mm256_mul_pd(Z, _mm256_fmsub_pd(five, z2, _mm256_mul_pd(three, r2))));
    const __m256d channel[kRealL3] = {
        _mm256_blendv_pd(_mm256_mul_pd(two, n0), _mm256_mul_pd(q0, inv), live),
        _mm256_blendv_pd(zero, _mm256_mul_pd(n1, _mm256_mul_pd(Y, _mm256_mul_pd(t1, inv))), live),
        _mm256_blendv_pd(zero, _mm256_mul_pd(n1, _mm256_mul_pd(X, _mm256_mul_pd(t1, inv))), live),
        _mm256_blendv_pd(zero, _mm256_mul_pd(n2, _mm256_mul_pd(Z, _mm256_mul_pd(_mm256_sub_pd(y2, x2), inv))), live),
        _mm256_blendv_pd(zero, _mm256_mul_pd(n2, _mm256_mul_pd(two, _mm256_mul_pd(X, _mm256_mul_pd(Y, _mm256_mul_pd(Z, inv))))), live),
        _mm256_blendv_pd(zero, _mm256_mul_pd(n3, _mm256_mul_pd(Y, _mm256_mul_pd(_mm256_fmsub_pd(y2, one, _mm256_mul_pd(three, x2)), inv))), live),
        _mm256_blendv_pd(zero, _mm256_mul_pd(n3, _mm256_mul_pd(X, _mm256_mul_pd(_mm256_fmsub_pd(_mm256_mul_pd(three, y2), one, x2), inv))), live)};
    for (int k = 0; k < kRealL3; ++k) {
      alignas(32) double lane[4];
      _mm256_store_pd(lane, channel[k]);
      acc[k] += lane[0];
      acc[k] += lane[1];
      acc[k] += lane[2];
      acc[k] += lane[3];
    }
  }
  for (; bond < nBond; ++bond) {
    addRealYlm3Scalar(x[bond], y[bond], z[bond], norm, acc);
  }
}

__attribute__((target("fma"))) double realL3Dot(const double *a, const double *b) {
  double dot = a[0] * b[0];
  for (int channel = 1; channel < kRealL3; ++channel) {
    dot = __builtin_fma(2.0 * a[channel], b[channel], dot);
  }
  return dot;
}

bool hostHasAvx2Fma() {
  return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
}
#else
// Without x86 SIMD the dispatch below always takes the scalar branch;
// these keep the call sites well-formed.
double realL3DotScalar(const double *a, const double *b);

void realYlm3Accumulate(const double *x, const double *y, const double *z, int nBond,
                        double *acc) {
  const auto &norm = normQ3();
  for (int bond = 0; bond < nBond; ++bond) {
    addRealYlm3Scalar(x[bond], y[bond], z[bond], norm, acc);
  }
}

double realL3Dot(const double *a, const double *b) { return realL3DotScalar(a, b); }

bool hostHasAvx2Fma() { return false; }
#endif

double realL3DotScalar(const double *a, const double *b) {
  double dot = a[0] * b[0];
  for (int channel = 1; channel < kRealL3; ++channel) {
    dot += (2.0 * a[channel]) * b[channel];
  }
  return dot;
}

/**
 * @details Selects the neighbours CHILL and CHILL+ are defined over.
 *
 *  Both schemes fix the coordination at @f$N_b(i) = 4@f$: Moore, de la Llave,
 *  Welke, Scherlis and Molinero (Phys. Chem. Chem. Phys. 12, 4124, 2010) and
 *  Nguyen and Molinero (J. Phys. Chem. B 119, 9369, 2015) both average the
 *  spherical harmonics over the four nearest neighbours, not over whatever
 *  falls inside a distance cutoff. On an ordered lattice the two coincide; in
 *  disordered or interfacial water they do not, which is where the
 *  classification is actually used.
 *
 *  Where fewer than the requested neighbours are available the shorter list
 *  is returned, since there is no further neighbour to choose. A
 *  non-positive @a coordination keeps the whole row, which gives each atom
 *  its own neighbour count for systems without a fixed coordination.
 * @param[in] yCloud The input molSys::PointCloud.
 * @param[in] nList Row-ordered neighbour list, by atom ID.
 * @param[in] iatomIndex Index of the particle whose neighbours are wanted.
 * @param[in] coordination How many nearest neighbours to keep; non-positive
 *  keeps every neighbour in the row.
 * @return Cloud indices of the nearest neighbours, nearest first.
 */
std::vector<int> nearestNeighbourIndices(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, int iatomIndex,
    int coordination) {
  std::vector<std::pair<double, int>> candidates;
  candidates.reserve(nList[iatomIndex].size());
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  for (size_t j = 1; j < nList[iatomIndex].size(); j++) {
    const auto it = yCloud.idIndexMap.find(nList[iatomIndex][j]);
    if (it == yCloud.idIndexMap.end()) {
      continue;
    }
    candidates.emplace_back(gen::periodicDistSq(frame, yCloud, iatomIndex, it->second),
                            it->second);
  }

  const size_t keep = (coordination > 0)
                          ? std::min(static_cast<size_t>(coordination),
                                     candidates.size())
                          : candidates.size();
  std::partial_sort(candidates.begin(), candidates.begin() + keep,
                    candidates.end());

  std::vector<int> nearest;
  nearest.reserve(keep);
  for (size_t k = 0; k < keep; k++) {
    nearest.push_back(candidates[k].second);
  }
  return nearest;
}

void associatedP(int L, double x, std::vector<double> &p) {
  p.assign(static_cast<std::size_t>(L) + 1, 0.0);
  const double s = std::sqrt(std::max(0.0, 1.0 - x * x));
  double pmm = 1.0;
  for (int m = 0; m <= L; m++) {
    if (m > 0) {
      pmm *= -(2.0 * m - 1.0) * s;
    }
    if (m == L) {
      p[static_cast<std::size_t>(m)] = pmm;
    } else {
    double pll2 = pmm;
    double pll1 = x * (2.0 * m + 1.0) * pmm;
    if (L == m + 1) {
      p[static_cast<std::size_t>(m)] = pll1;
      continue;
    }
    double pll = 0.0;
    for (int ll = m + 2; ll <= L; ll++) {
      pll = (x * (2.0 * ll - 1.0) * pll1 - (ll + m - 1.0) * pll2) /
            static_cast<double>(ll - m);
      pll2 = pll1;
      pll1 = pll;
    }
    p[static_cast<std::size_t>(m)] = pll;
    }
  }
}

double nlmNorm(int l, int m) {
  constexpr double pi = 3.14159265358979323846;
  double f = (2.0 * l + 1.0) / (4.0 * pi);
  for (int k = 1; k <= 2 * m; k++) {
    f /= static_cast<double>(l - m + k);
  }
  return std::sqrt(f);
}

void spheriHarmoGeneral(int orderL, std::array<double, 2> angles,
                        std::vector<std::complex<double>> &out) {
  out.assign(static_cast<std::size_t>(2 * orderL + 1), {0.0, 0.0});
  std::vector<double> p;
  associatedP(orderL, std::cos(angles[1]), p);
  const double phi = angles[0];
  for (int m = 0; m <= orderL; m++) {
    const double amp = nlmNorm(orderL, m) * p[static_cast<std::size_t>(m)];
    const std::complex<double> yp =
        amp * std::polar(1.0, static_cast<double>(m) * phi);
    out[static_cast<std::size_t>(orderL + m)] = yp;
    if (m == 0) {
      continue;
    }
    const double sign = (m % 2 == 0) ? 1.0 : -1.0;
    out[static_cast<std::size_t>(orderL - m)] = sign * std::conj(yp);
  }
}

void harmonicPairInto(int orderL, std::array<double, 2> angles,
                      std::vector<std::complex<double>> &out) {
  const AngularTerms terms(angles[1], angles[0]);
  out.resize(static_cast<std::size_t>(2 * orderL + 1));
  for (int m = 0; m <= orderL; m++) {
    const auto [negative, positive] = harmonicPair(orderL, m, terms);
    out[static_cast<std::size_t>(orderL - m)] = negative;
    out[static_cast<std::size_t>(orderL + m)] = positive;
  }
}

} // namespace

void sph::spheriHarmoInto(int orderL, std::array<double, 2> radialCoord,
                          std::vector<std::complex<double>> &out) {
  if (orderL < 0) {
    out.clear();
    return;
  }
  if (orderL == 3 || orderL == 4 || orderL == 6 || orderL == 8) {
    harmonicPairInto(orderL, radialCoord, out);
    return;
  }
  spheriHarmoGeneral(orderL, radialCoord, out);
}

std::vector<std::complex<double>>
sph::spheriHarmo(int orderL, std::array<double, 2> radialCoord) {
  std::vector<std::complex<double>> out;
  spheriHarmoInto(orderL, radialCoord, out);
  return out;
}


/**
 * @details Calculates @f$Q_8@f$ from the shared trigonometric terms. Used for
 *  distinguishing crystal phases whose @f$q_4@f$-@f$q_6@f$ signatures overlap.
 * @param[in] angles The azimuth and polar angles for a particular particle
 * @return a complex vector of length @f$17@f$, ordered @f$m = -8 \ldots 8@f$
 */
std::vector<std::complex<double>>
sph::lookupTableQ8Vec(std::array<double, 2> angles) {
  const AngularTerms terms(angles[1], angles[0]);
  std::vector<std::complex<double>> result(17);
  for (int m = 0; m <= 8; m++) {
    const auto [negative, positive] = harmonicPair(8, m, terms);
    result[8 - m] = negative;
    result[8 + m] = positive;
  }
  return result;
}

/**
 * @details Single order of @f$Q_8@f$.
 * @param[in] m Table index, @f$0 \ldots 16@f$, mapping to @f$m = -8 \ldots 8@f$
 * @param[in] angles The azimuth and polar angles for a particular particle
 * @return The complex value of @f$Y_{8m}@f$
 */
std::complex<double> sph::lookupTableQ8(int m, std::array<double, 2> angles) {
  if (m < 0 || m > 16) {
    return {0.0, 0.0};
  }
  const AngularTerms terms(angles[1], angles[0]);
  const int order = m - 8;
  const auto [negative, positive] = harmonicPair(8, std::abs(order), terms);
  return (order < 0) ? negative : positive;
}

/**
 * @details Function for the azimuth and polar angles, given the Cartesian
 *  coordinates. Pure trigonometry, no external dependencies.
 *  @param[in] cartCoord The Cartesian coordinates of a particular point
 *  @return a double array, holding the azimuth and polar angles
 */
std::array<double, 2> sph::radialCoord(std::array<double, 3> cartCoord) {
  std::array<double, 2> result;
  double r = std::sqrt(cartCoord[0] * cartCoord[0] +
                       cartCoord[1] * cartCoord[1] +
                       cartCoord[2] * cartCoord[2]);
  result[0] = std::atan2(cartCoord[0], cartCoord[1]); // azimuth (phi)
  result[1] = (r > 0.0) ? std::acos(cartCoord[2] / r) : 0.0; // polar (theta)
  return result;
}

/**
 * @details Calculates @f$Q_3@f$ using hard-coded look-up values.
 * @deprecated It is recommended to use the Boost version of this function,
 * sph::spheriHarmo, instead.
 * @param[in] angles The azimuth and polar angles of a particular point
 * @return a complex vector, of length @f$7@f$, calculated using spherical
 * harmonics
 */
std::vector<std::complex<double>>
sph::lookupTableQ3Vec(std::array<double, 2> angles) {
  const AngularTerms terms(angles[1], angles[0]);
  std::vector<std::complex<double>> result(7);

  for (int m = 0; m <= 3; m++) {
    const auto [negative, positive] = harmonicPair(3, m, terms);
    result[3 - m] = negative;
    result[3 + m] = positive;
  }

  return result;
}

/**
 *  @details Look-up hard-coded values for @f$Q_3@f$
 *
 * It is recommended to use the Boost version of this function,
 * sph::spheriHarmo, instead.
 *
 *  @param[in] m An int such that @f$-3<=m<=3@f$
 *  @param[in] angles The azimuth and polar angles for a particular particle
 *  @return a complex vector, of length @f$7@f$, calculated using hard-coded
 *   values
 */
std::complex<double> sph::lookupTableQ3(int m, std::array<double, 2> angles) {
  if (m < 0 || m > 6) {
    return {0.0, 0.0};
  }

  const AngularTerms terms(angles[1], angles[0]);
  const int order = m - 3; // table index 0..6 maps to m = -3..3
  const auto [negative, positive] = harmonicPair(3, std::abs(order), terms);

  return (order < 0) ? negative : positive;
}

/**
 * @details Calculates @f$Q_4@f$ from the shared trigonometric terms.
 *  @f$q_4@f$ paired with @f$q_6@f$ is the standard discriminator between
 *  close-packed and body-centred environments.
 * @param[in] angles The azimuth and polar angles for a particular particle
 * @return a complex vector of length @f$9@f$, ordered @f$m = -4 \ldots 4@f$
 */
std::vector<std::complex<double>>
sph::lookupTableQ4Vec(std::array<double, 2> angles) {
  const AngularTerms terms(angles[1], angles[0]);
  std::vector<std::complex<double>> result(9);

  for (int m = 0; m <= 4; m++) {
    const auto [negative, positive] = harmonicPair(4, m, terms);
    result[4 - m] = negative;
    result[4 + m] = positive;
  }

  return result;
}

/**
 * @details Single order of @f$Q_4@f$.
 * @param[in] m Table index, @f$0 \ldots 8@f$, mapping to @f$m = -4 \ldots 4@f$
 * @param[in] angles The azimuth and polar angles for a particular particle
 * @return The complex value of @f$Y_{4m}@f$
 */
std::complex<double> sph::lookupTableQ4(int m, std::array<double, 2> angles) {
  if (m < 0 || m > 8) {
    return {0.0, 0.0};
  }

  const AngularTerms terms(angles[1], angles[0]);
  const int order = m - 4; // table index 0..8 maps to m = -4..4
  const auto [negative, positive] = harmonicPair(4, std::abs(order), terms);

  return (order < 0) ? negative : positive;
}

/**
 * @details Calculates @f$Q_6@f$ using hard-coded values.
 *
 * It is recommended to use the Boost version of this function,
 * sph::spheriHarmo, instead.
 *
 *  @param[in] angles The azimuth and polar angles for a particular particle
 *  @return a complex vector, of length @f$13@f$, calculated using hard-coded
 *   values
 */
std::vector<std::complex<double>>
sph::lookupTableQ6Vec(std::array<double, 2> angles) {
  const AngularTerms terms(angles[1], angles[0]);
  std::vector<std::complex<double>> result(13);

  for (int m = 0; m <= 6; m++) {
    const auto [negative, positive] = harmonicPair(6, m, terms);
    result[6 - m] = negative;
    result[6 + m] = positive;
  }

  return result;
}

/**
 * @details Hard-coded calculations for determining @f$Q_6@f$.
 *
 * It is recommended to use the general Boost version of this function,
 * sph::spheriHarmo, instead.
 *
 *  @param[in] m An int such that @f$-6<=m<=6@f$
 *  @param[in] angles The azimuth and polar angles for a particular particle
 *  @return a complex vector, of length @f$13@f$, calculated using hard-coded
 *   values
 */
std::complex<double> sph::lookupTableQ6(int m, std::array<double, 2> angles) {
  if (m < 0 || m > 12) {
    return {0.0, 0.0};
  }

  const AngularTerms terms(angles[1], angles[0]);
  const int order = m - 6; // table index 0..12 maps to m = -6..6
  const auto [negative, positive] = harmonicPair(6, std::abs(order), terms);

  return (order < 0) ? negative : positive;
}

/**
 * @details The bond-correlation engine shared by every rule set. Averages the
 *  l=3 spherical harmonics over the rule's nearest-neighbour count, forms the
 *  normalised correlation @f$c_{ij}@f$ for the same bonds, and classifies
 *  each bond against the rule's thresholds (boundaries inclusive). getCorrel
 *  and getCorrelPlus call this with the CHILL and CHILL+ water rules; other
 *  materials pass their own registered rules.
 */
static void classifyBondsLeftFp64(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, const chill::BondClassifier &rule,
    bool isSlice, const std::vector<int> *csrOff = nullptr,
    const std::vector<int> *csrCols = nullptr, bool closed = true) {
  seams::prepareCompute();
  [[maybe_unused]] const seams::TeamWidth team(yCloud.nop);
  constexpr int l = 3;
  static_assert(l == 3, "the closed form below is the l = 3 polynomial");
  const int nop = yCloud.nop;
  const int keepMax = rule.coordinationNumber;
  const gen::FracBox frame = gen::makeFracBox(yCloud);

  int maxId = 0;
  for (int i = 0; i < nop; i++) {
    maxId = std::max(maxId, yCloud.pts[static_cast<std::size_t>(i)].atomID);
  }
  const bool denseIds = maxId > 0 && maxId < nop * 4;
  std::vector<int> idToIndex;
  if (denseIds) {
    idToIndex.assign(static_cast<std::size_t>(maxId) + 1, -1);
    for (int i = 0; i < nop; i++) {
      const int id = yCloud.pts[static_cast<std::size_t>(i)].atomID;
      if (id >= 0) {
        idToIndex[static_cast<std::size_t>(id)] = i;
      }
    }
  }
  const auto resolve = [&](int id) {
    if (denseIds) {
      return (id >= 0 && id <= maxId) ? idToIndex[static_cast<std::size_t>(id)] : -1;
    }
    const auto it = yCloud.idIndexMap.find(id);
    return it == yCloud.idIndexMap.end() ? -1 : it->second;
  };

  // Real l = 3 coefficients, unnormalised by the neighbour count: that
  // factor cancels in c_ij. Coordinates are copied first so the parallel
  // loop does not read a Point while another thread clears that Point's
  // bond list.
  struct Xyz {
    double x, y, z;
  };
  std::vector<Xyz> xyz(static_cast<std::size_t>(nop));
  for (int i = 0; i < nop; i++) {
    auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(i)] = {p.x, p.y, p.z};
    p.c_ij.clear();
  }
  std::vector<double> qAtom(static_cast<std::size_t>(nop) * kRealL3, 0.0);
  std::vector<double> qNorm(static_cast<std::size_t>(nop), 0.0);
  const bool fixedKeep = keepMax > 0;
  // One thread walks atoms in order, so a whole cutoff shell is one
  // buffer instead of a vector per atom.
  const bool serialShell =
      !fixedKeep &&
#ifdef SEAMS_HAS_OPENMP
      omp_get_max_threads() <= 1
#else
      true
#endif
      ;
  std::vector<int> chosen;
  std::vector<std::vector<int>> ragged;
  std::vector<int> flat;
  std::vector<int> flatAt(static_cast<std::size_t>(nop) + 1, 0);
  std::vector<int> nkeep(static_cast<std::size_t>(nop), 0);
  if (fixedKeep) {
    chosen.assign(static_cast<std::size_t>(nop) * static_cast<std::size_t>(keepMax), -1);
  } else if (!serialShell) {
    ragged.resize(static_cast<std::size_t>(nop));
  } else {
    flat.reserve(static_cast<std::size_t>(nop) * 8);
  }

  const bool useAvx = hostHasAvx2Fma();
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nop >= 16384 && omp_get_max_threads() > 1)
#endif
  for (int iatom = 0; iatom < nop; iatom++) {
    if (csrOff == nullptr && static_cast<std::size_t>(iatom) >= nList.size()) {
      continue;
    }
    // The distance and the minimum-image delta are one evaluation. A short
    // shell stays on the stack; a wider cutoff row spills, in the same order.
    struct Near {
      double d2;
      int jatom;
      double x, y, z;
      bool operator<(const Near &o) const {
        return d2 < o.d2 || (d2 == o.d2 && jatom < o.jatom);
      }
    };
    Near stackNear[32];
    std::vector<Near> heapNear;
    Near *candidates = stackNear;
    int nCand = 0;
    const Xyz pi = xyz[static_cast<std::size_t>(iatom)];
    const auto consider = [&](int jatom) {
      if (jatom < 0 || jatom >= nop || jatom == iatom) {
        return;
      }
      const Xyz pj = xyz[static_cast<std::size_t>(jatom)];
      const auto dr = gen::fracDelta(frame, pi.x, pi.y, pi.z, pj.x, pj.y, pj.z);
      if (nCand == 32 && candidates == stackNear) {
        heapNear.assign(stackNear, stackNear + 32);
        candidates = heapNear.data();
      }
      if (candidates == stackNear) {
        candidates[nCand] = {dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2], jatom, dr[0],
                             dr[1], dr[2]};
      } else {
        heapNear.push_back(
            {dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2], jatom, dr[0], dr[1], dr[2]});
        candidates = heapNear.data();
      }
      ++nCand;
    };
    if (csrOff != nullptr && csrCols != nullptr &&
        static_cast<int>(csrOff->size()) > iatom) {
      const int a = (*csrOff)[static_cast<std::size_t>(iatom)];
      const int b = (*csrOff)[static_cast<std::size_t>(iatom) + 1];
      for (int k = a; k < b; ++k) {
        consider((*csrCols)[static_cast<std::size_t>(k)]);
      }
    } else if (static_cast<std::size_t>(iatom) < nList.size()) {
      const auto &row = nList[static_cast<std::size_t>(iatom)];
      for (std::size_t j = 1; j < row.size(); ++j) {
        consider(resolve(static_cast<int>(row[j])));
      }
    }
    const int keep =
        (keepMax > 0) ? std::min(keepMax, nCand) : nCand;
    if (keep == 0) {
      continue;
    }
    std::partial_sort(candidates, candidates + keep, candidates + nCand);
    nkeep[static_cast<std::size_t>(iatom)] = keep;
    double acc[kRealL3] = {};
    for (int j = 0; j < keep;) {
      double xs[4];
      double ys[4];
      double zs[4];
      int batch = 0;
      for (; batch < 4 && j < keep; ++j, ++batch) {
        const int jatom = candidates[j].jatom;
        if (fixedKeep) {
          chosen[static_cast<std::size_t>(iatom) * static_cast<std::size_t>(keepMax) +
                 static_cast<std::size_t>(j)] = jatom;
        } else if (serialShell) {
          flat.push_back(jatom);
        } else {
          ragged[static_cast<std::size_t>(iatom)].push_back(jatom);
        }
        xs[batch] = candidates[j].x;
        ys[batch] = candidates[j].y;
        zs[batch] = candidates[j].z;
      }
      if (useAvx) {
        realYlm3Accumulate(xs, ys, zs, batch, acc);
      } else {
        for (int item = 0; item < batch; ++item) {
          addRealYlm3Scalar(xs[item], ys[item], zs[item], normQ3(), acc);
        }
      }
    }
    double *dst = qAtom.data() + static_cast<std::size_t>(iatom) * kRealL3;
    for (int channel = 0; channel < kRealL3; ++channel) {
      dst[channel] = acc[channel];
    }
    qNorm[static_cast<std::size_t>(iatom)] = useAvx ? realL3Dot(acc, acc) : realL3DotScalar(acc, acc);
  }
  if (serialShell) {
    for (int i = 0; i < nop; i++) {
      flatAt[static_cast<std::size_t>(i + 1)] =
          flatAt[static_cast<std::size_t>(i)] + nkeep[static_cast<std::size_t>(i)];
    }
  }
  const auto neighbourAt = [&](int atom, int slot) {
    if (fixedKeep) {
      return chosen[static_cast<std::size_t>(atom) * static_cast<std::size_t>(keepMax) +
                    static_cast<std::size_t>(slot)];
    }
    if (serialShell) {
      return flat[static_cast<std::size_t>(flatAt[static_cast<std::size_t>(atom)] + slot)];
    }
    return ragged[static_cast<std::size_t>(atom)][static_cast<std::size_t>(slot)];
  };

  // c_ij from the real l = 3 inner product. Norms were formed once per atom.
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nop >= 16384 && omp_get_max_threads() > 1)
#endif
  for (int iatom = 0; iatom < nop; iatom++) {
    if (isSlice && !yCloud.pts[static_cast<std::size_t>(iatom)].inSlice) {
      continue;
    }
    const int keep = nkeep[static_cast<std::size_t>(iatom)];
    if (keep <= 0) {
      continue;
    }
    const double *qI = qAtom.data() + static_cast<std::size_t>(iatom) * kRealL3;
    const double normI = qNorm[static_cast<std::size_t>(iatom)];
    auto &bonds = yCloud.pts[static_cast<std::size_t>(iatom)].c_ij;
    bonds.reserve(static_cast<std::size_t>(keep));
    for (int slot = 0; slot < keep; ++slot) {
      const int jatom = neighbourAt(iatom, slot);
      if (jatom < 0) {
        continue;
      }
      const double *qJ = qAtom.data() + static_cast<std::size_t>(jatom) * kRealL3;
      const double dot = useAvx ? realL3Dot(qI, qJ) : realL3DotScalar(qI, qJ);
      const double cij = correlationOf(dot, normI, qNorm[static_cast<std::size_t>(jatom)]);
      molSys::Result temp_cij;
      temp_cij.c_value = cij;
      temp_cij.classifier = chill::classifyBond(cij, rule, closed);
      bonds.push_back(temp_cij);
    }
  }
}

void chill::classifyBonds(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, const chill::BondClassifier &rule,
    bool isSlice) {
  chill::classifyBonds(yCloud, nList, rule, isSlice, true);
}

void chill::classifyBonds(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, const chill::BondClassifier &rule,
    bool isSlice, bool closed) {
  // Left-to-right fp64 is the path CHILL and CHILL+ run on. A different
  // product reduction or a non-fp64 dot product keeps the pairwise sum.
  if (seams::product::reduction() == seams::product::Reduction::Left &&
      seams::product::mode() == seams::product::Mode::Fp64) {
    classifyBondsLeftFp64(yCloud, nList, rule, isSlice, nullptr, nullptr, closed);
    return;
  }
  constexpr int l = 3;
  chill::QlmAtom QlmTotal; // Qlm for each iatom
  std::vector<std::complex<double>>
      yl; // temp q_lm for each pair of iatom and jatom

  QlmTotal.ptq.resize(yCloud.nop);
  for (auto &point : yCloud.pts) {
    point.c_ij.clear();
  }

  // Loop through the neighbour list
  for (int iatom = 0; iatom < nList.size(); iatom++) {
    // CHILL and CHILL+ average over the four nearest neighbours, not over
    // every neighbour that falls inside the cutoff; other rule sets bring
    // their own coordination
    const std::vector<int> nearest = nearestNeighbourIndices(
        yCloud, nList, iatom, rule.coordinationNumber);
    const int nnumNeighbours = static_cast<int>(nearest.size());
    const bool pairwiseQlm =
        seams::product::reduction() != seams::product::Reduction::Left;
    std::vector<double> qlmRe;
    std::vector<double> qlmIm;
    if (pairwiseQlm) {
      qlmRe.assign(static_cast<std::size_t>(nnumNeighbours) * (2 * l + 1), 0.0);
      qlmIm.assign(static_cast<std::size_t>(nnumNeighbours) * (2 * l + 1), 0.0);
    }
    for (int j = 0; j < nnumNeighbours; j++) {
      const int jatomIndex = nearest[j];

      // Get the relative distance now that the index values are known
      const std::array<double, 3> delta =
          gen::relDist(yCloud, iatom, jatomIndex);
      // radialCoord carries the r > 0 guard, so coincident atoms yield a
      // defined angle instead of acos(0/0)
      const std::array<double, 2> angles = sph::radialCoord(delta);

      // Now add over all nearest neighbours
      sph::spheriHarmoInto(l, angles, yl);
      if (pairwiseQlm) {
        for (int m = 0; m < 2 * l + 1; m++) {
          qlmRe[static_cast<std::size_t>(m) * nnumNeighbours + j] = yl[m].real();
          qlmIm[static_cast<std::size_t>(m) * nnumNeighbours + j] = yl[m].imag();
        }
        continue;
      }
      if (j == 0) {
        QlmTotal.ptq[iatom].ylm = yl;
        continue;
      }
      for (int m = 0; m < 2 * l + 1; m++) {
        QlmTotal.ptq[iatom].ylm[m] += yl[m];
      }
    } // End of loop over nearest neighbours
    if (pairwiseQlm) {
      QlmTotal.ptq[iatom].ylm.assign(static_cast<std::size_t>(2 * l + 1), {0.0, 0.0});
      for (int m = 0; m < 2 * l + 1; m++) {
        const double re = seams::product::reduce(
            qlmRe.data() + static_cast<std::size_t>(m) * nnumNeighbours,
            nnumNeighbours, seams::product::reduction());
        const double im = seams::product::reduce(
            qlmIm.data() + static_cast<std::size_t>(m) * nnumNeighbours,
            nnumNeighbours, seams::product::reduction());
        QlmTotal.ptq[iatom].ylm[static_cast<std::size_t>(m)] = {re, im};
      }
    }

    // Average over the bond count
    QlmTotal.ptq[iatom].ylm =
        gen::avgVector(QlmTotal.ptq[iatom].ylm, l, nnumNeighbours);
  } // End of looping over all iatom

  // ------------------------------------------------
  // Now that you have all qlm for the particles, find c_ij
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // if this is a slice and the particle is not in the slice then skip
    if (isSlice && !yCloud.pts[iatom].inSlice) {
      continue;
    }
    // The same neighbours the q_lm were averaged over, so that the bond
    // count the classification tables are written against matches
    const std::vector<int> nearestBonds = nearestNeighbourIndices(
        yCloud, nList, iatom, rule.coordinationNumber);
    yCloud.pts[iatom].c_ij.reserve(nearestBonds.size());
    const seams::product::Mode productMode = seams::product::mode();
    for (const int jatomIndex : nearestBonds) {
      std::complex<double> dot_product = {0, 0};
      std::complex<double> Inorm = {0, 0};
      std::complex<double> Jnorm = {0, 0};
      if (productMode == seams::product::Mode::Fp64) {
        for (int m = 0; m < 2 * l + 1; m++) {
          const std::complex<double> qI = QlmTotal.ptq[iatom].ylm[m];
          const std::complex<double> qJ = QlmTotal.ptq[jatomIndex].ylm[m];
          dot_product += qI * std::conj(qJ); // unnormalized
          Inorm += qI * std::conj(qI);
          Jnorm += qJ * std::conj(qJ);
        } // end loop over m components
      } else {
        const auto got = seams::product::complexDot(
            QlmTotal.ptq[iatom].ylm.data(),
            QlmTotal.ptq[jatomIndex].ylm.data(), 2 * l + 1, productMode);
        dot_product = {got.re, got.im};
        Inorm = {got.normA, 0.0};
        Jnorm = {got.normB, 0.0};
      }
      const double cij_real =
          correlationOf(dot_product.real(), Inorm.real(), Jnorm.real());
      molSys::Result temp_cij;
      temp_cij.c_value = cij_real;
      temp_cij.classifier = chill::classifyBond(cij_real, rule, closed);
      yCloud.pts[iatom].c_ij.push_back(temp_cij);
    } // end loop over nearest neighbours
  }
}

chill::BondClassifier chill::chillRule() { return {-0.8, -0.2, -0.05, 4}; }

chill::BondClassifier chill::chillPlusRule() { return {-0.8, -0.35, 0.25, 4}; }

namespace {
std::unordered_map<std::string, chill::BondClassifier> g_custom_rules;
std::unordered_map<std::string, chill::BondClassifier> &customRules() {
  return g_custom_rules;
}
} // namespace

chill::BondClassifier chill::bondClassifier(const std::string &name) {
  const auto &custom = customRules();
  const auto found = custom.find(name);
  if (found != custom.end()) {
    return found->second;
  }
  if (name == "CHILL") {
    return chillRule();
  }
  if (name == "CHILL+") {
    return chillPlusRule();
  }
  throw std::out_of_range("bond classifier " + name);
}

void chill::registerBondClassifier(const std::string &name,
                                   const chill::BondClassifier &rule) {
  customRules()[name] = rule;
}

std::vector<std::string> chill::bondClassifierNames() {
  std::vector<std::string> names;
  const auto &custom = customRules();
  names.reserve(custom.size() + 2);
  names.push_back("CHILL");
  names.push_back("CHILL+");
  for (const auto &entry : custom) {
    names.push_back(entry.first);
  }
  std::sort(names.begin(), names.end());
  names.erase(std::unique(names.begin(), names.end()), names.end());
  return names;
}

//! Gets c_ij according to the CHILL rule set
void
chill::getCorrel(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                 const std::vector<std::vector<int>> &nList, bool isSlice,
                 int coordinationNumber) {
  chill::BondClassifier rule = chill::chillRule();
  rule.coordinationNumber = coordinationNumber;
  chill::classifyBonds(yCloud, nList, rule, isSlice);
}

//! Classifies each atom according to the CHILL algorithm without printing
void chill::getIceTypeNoPrint(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, bool isSlice) {
  int ih, ic, water, interIce, unknown, total; // No. of particles of each type
  ih = ic = water = unknown = interIce = total = 0;
  int num_staggrd, num_eclipsd, na;
  molSys::bond_type bondType;
  int nnumNeighbours; // Number of nearest neighbours

  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // if(yCloud.pts[iatom].type!=typeO){continue;}
    // if this is a slice and the particle is not in the slice
    // then skip
    if (isSlice) {
      if (yCloud.pts[iatom].inSlice == false) {
        continue;
      }
    }
    total++; // Update the total number of atoms considered. Change this to
    // check for slices
    num_staggrd = num_eclipsd = na =
        0; // init to zero before loop through neighbours

    // c_ij is the bond list the correlation pass actually produced; the
    // neighbour list is not a proxy for its length
    nnumNeighbours = yCloud.pts[iatom].c_ij.size();
    // Loop through the bond cij and get the number of staggered, eclipsed bonds
    for (int j = 0; j < nnumNeighbours; j++) {
      bondType = yCloud.pts[iatom].c_ij[j].classifier;
      if (bondType == molSys::bond_type::eclipsed) {
        num_eclipsd++;
      } else if (bondType == molSys::bond_type::staggered) {
        num_staggrd++;
      } else {
        na++;
      }
    } // End of loop through neighbours

    yCloud.pts[iatom].iceType = molSys::atom_state_type::unclassified; // default
    if (nnumNeighbours == 4) {
      if (num_staggrd >= 4) {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::cubic;
        ic++;
      } else if (num_eclipsd == 1 && num_staggrd == 3) {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::hexagonal;
        ih++;
      } else if (chill::isInterfacial(yCloud, nList, iatom, num_staggrd, num_eclipsd)) {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::interfacial;
        interIce++;
      } else {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::water;
        water++;
      }
    } else {
      yCloud.pts[iatom].iceType = molSys::atom_state_type::water;
      water++;
    }

  } // End of loop through every iatom

  return;
}

//! Classifies each atom according to the CHILL algorithm
void
chill::getIceType(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                  const std::vector<std::vector<int>> &nList, std::string path,
                  int firstFrame, bool isSlice, std::string outputFileName) {
  int ih, ic, water, interIce, unknown, total; // No. of particles of each type
  ih = ic = water = unknown = interIce = total = 0;
  int num_staggrd, num_eclipsd, na;
  molSys::bond_type bondType;
  int nnumNeighbours; // Number of nearest neighbours

  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // if(yCloud.pts[iatom].type!=typeO){continue;}
    // if this is a slice and the particle is not in the slice
    // then skip
    if (isSlice) {
      if (yCloud.pts[iatom].inSlice == false) {
        continue;
      }
    }
    total++; // Update the total number of atoms considered. Change this to
    // check for slices
    num_staggrd = num_eclipsd = na =
        0; // init to zero before loop through neighbours

    // c_ij is the bond list the correlation pass actually produced; the
    // neighbour list is not a proxy for its length
    nnumNeighbours = yCloud.pts[iatom].c_ij.size();
    // Loop through the bond cij and get the number of staggered, eclipsed bonds
    for (int j = 0; j < nnumNeighbours; j++) {
      bondType = yCloud.pts[iatom].c_ij[j].classifier;
      if (bondType == molSys::bond_type::eclipsed) {
        num_eclipsd++;
      } else if (bondType == molSys::bond_type::staggered) {
        num_staggrd++;
      } else {
        na++;
      }
    } // End of loop through neighbours

    yCloud.pts[iatom].iceType = molSys::atom_state_type::unclassified; // default
    if (nnumNeighbours == 4) {
      if (num_staggrd >= 4) {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::cubic;
        ic++;
      } else if (num_eclipsd == 1 && num_staggrd == 3) {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::hexagonal;
        ih++;
      } else if (chill::isInterfacial(yCloud, nList, iatom, num_staggrd, num_eclipsd)) {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::interfacial;
        interIce++;
      } else {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::water;
        water++;
      }
    } else {
      yCloud.pts[iatom].iceType = molSys::atom_state_type::water;
      water++;
    }

  } // End of loop through every iatom

  // water = total - ic -ih;

  // --------------------
  // Create the directories if needed
  sout::makePath(path);
  std::string outputDirName = path + "bop";
  sout::makePath(outputDirName);
  // --------------------

  // Print to file
  std::ofstream outputFile;
  outputFile.open(path + "bop/" + outputFileName, std::ios_base::app);
  // --------------------
  // Write out the comment line for the first frame
  if (yCloud.currentFrame == firstFrame) {
    outputFile << "Frame Ic Ih Interfacial Water Total \n";
  }
  // --------------------
  outputFile << yCloud.currentFrame << " " << ic << " " << ih << " "
             << interIce << " " << water << " " << total << "\n";
  outputFile.close();

  return;
}

// cmean^4 above this cut, with cmean < 0 and each correlation in [-1, 1],
// forces every one of the four correlations to be <= -0.84395. That is
// staggered under the published edge -0.8. The cut is 1e-6 above the
// highest non-cubic score in discovery/pysr_prefilter.py, including the
// held-out seeds 91001, 92002 and 93003. A shell of four then labels cubic.
constexpr double kSymbolicCubicCut = 0.8528497963226561 + 1e-6;

bool orthoChill(molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
               int typeI, bool certify = false) {
  const int n = yCloud.nop;
  const double ox = yCloud.boxLow.size() > 0 ? yCloud.boxLow[0] : 0.0;
  const double oy = yCloud.boxLow.size() > 1 ? yCloud.boxLow[1] : 0.0;
  const double oz = yCloud.boxLow.size() > 2 ? yCloud.boxLow[2] : 0.0;
  const bool tiltBox = yCloud.box.size() >= 6;
  const gen::FracBox fracBox = gen::makeFracBox(yCloud);
  const double Lx = tiltBox ? fracBox.lx : yCloud.box[0];
  const double Ly = tiltBox ? fracBox.ly : yCloud.box[1];
  const double Lz = tiltBox ? fracBox.lz : yCloud.box[2];
  if (!(Lx > 0.0 && Ly > 0.0 && Lz > 0.0)) {
    return false;
  }
  // This walk is OpenMP only. adoptThreads pauses the team so a Rayon
  // search can have the cores, and that pause is what split the
  // eight-thread samples into a fast group near 8 ms and a slow group
  // near 12 ms.
  seams::prepareCompute();
  // The scratch is shared with the OpenMP team. A second frame running
  // in another worker waits, instead of writing that team's rows.
  static std::mutex scratchMu;
  std::lock_guard<std::mutex> scratchGuard(scratchMu);
  struct Scratch {
    std::vector<double> xyz;
    std::vector<int> head, next;
    std::vector<int> offsets, cols;
    std::vector<int> nkeep, chosen;
    std::vector<double> q, qn;
    int n = -1;
    int ncell = -1;
  };
  static Scratch scratch;
  const int nx = std::max(1, static_cast<int>(std::ceil(Lx / cutoff)));
  const int ny = std::max(1, static_cast<int>(std::ceil(Ly / cutoff)));
  const int nz = std::max(1, static_cast<int>(std::ceil(Lz / cutoff)));
  const int ncell = nx * ny * nz;
  if (scratch.n != n) {
    scratch.xyz.resize(static_cast<std::size_t>(n) * 3);
    scratch.next.assign(static_cast<std::size_t>(n), -1);
    scratch.nkeep.assign(static_cast<std::size_t>(n), 0);
    scratch.chosen.assign(static_cast<std::size_t>(n) * 4, -1);
    scratch.q.assign(static_cast<std::size_t>(n) * kRealL3, 0.0);
    scratch.qn.assign(static_cast<std::size_t>(n), 0.0);
    scratch.n = n;
  }
  if (scratch.ncell != ncell) {
    scratch.head.assign(static_cast<std::size_t>(ncell), -1);
    scratch.ncell = ncell;
  } else {
    std::fill(scratch.head.begin(), scratch.head.end(), -1);
  }
  bool uniform = true;
  for (int i = 0; i < n; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    auto wrap = [](double v, double origin, double L) {
      double t = (v - origin) / L;
      t -= std::floor(t);
      if (!(t < 1.0)) {
        t = 0.0;
      }
      return t * L;
    };
    if (tiltBox) {
      scratch.xyz[static_cast<std::size_t>(3 * i)] = p.x;
      scratch.xyz[static_cast<std::size_t>(3 * i) + 1] = p.y;
      scratch.xyz[static_cast<std::size_t>(3 * i) + 2] = p.z;
    } else {
      scratch.xyz[static_cast<std::size_t>(3 * i)] = wrap(p.x, ox, Lx);
      scratch.xyz[static_cast<std::size_t>(3 * i) + 1] = wrap(p.y, oy, Ly);
      scratch.xyz[static_cast<std::size_t>(3 * i) + 2] = wrap(p.z, oz, Lz);
    }
    if (p.type != typeI) {
      uniform = false;
    }
    if (!certify) {
      yCloud.pts[static_cast<std::size_t>(i)].c_ij.clear();
    }
  }
  if (certify && !uniform) {
    return false;
  }
#ifdef SEAMS_HAS_OPENMP
  // Caps the team and, when that is narrower than the process, binds
  // those workers to that many cores of the allowed set.
  [[maybe_unused]] const seams::TeamWidth teamWidth(n);
#endif
  std::vector<int> subset;
  const std::vector<int> *which = nullptr;
  if (!uniform) {
    subset.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
      if (yCloud.pts[static_cast<std::size_t>(i)].type == typeI) {
        subset.push_back(i);
      }
    }
    which = &subset;
  }
  const bool timing = std::getenv("SEAMS_CHILL_TIME") != nullptr;
  const auto tCsr0 = std::chrono::steady_clock::now();
  const bool csrOk =
      nneigh::orthoFullCsr(yCloud, which, cutoff, scratch.offsets, scratch.cols);
  const auto tCsr1 = std::chrono::steady_clock::now();
  if (tiltBox && !csrOk) {
    return false;
  }
  const double invx = static_cast<double>(nx) / Lx;
  const double invy = static_cast<double>(ny) / Ly;
  const double invz = static_cast<double>(nz) / Lz;
  if (!csrOk) {
  for (int i = 0; i < n; i++) {
    if (!uniform && yCloud.pts[static_cast<std::size_t>(i)].type != typeI) {
      scratch.next[static_cast<std::size_t>(i)] = -1;
      continue;
    }
    const int ix = std::min(nx - 1, static_cast<int>(scratch.xyz[static_cast<std::size_t>(3 * i)] * invx));
    const int iy = std::min(ny - 1, static_cast<int>(scratch.xyz[static_cast<std::size_t>(3 * i) + 1] * invy));
    const int iz = std::min(nz - 1, static_cast<int>(scratch.xyz[static_cast<std::size_t>(3 * i) + 2] * invz));
    const int c = (iz * ny + iy) * nx + ix;
    scratch.next[static_cast<std::size_t>(i)] = scratch.head[static_cast<std::size_t>(c)];
    scratch.head[static_cast<std::size_t>(c)] = i;
  }
  }
  const double cut2 = cutoff * cutoff;
  const double hx = 0.5 * Lx;
  const double hy = 0.5 * Ly;
  const double hz = 0.5 * Lz;
  const auto &norm = normQ3();
  const bool useAvx = hostHasAvx2Fma();
  const chill::BondClassifier rule = chill::chillPlusRule();
  const auto tWalk = std::chrono::steady_clock::now();
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
  for (int i = 0; i < n; i++) {
    struct Near {
      double d2;
      int j;
      double x, y, z;
      bool operator<(const Near &o) const { return d2 < o.d2 || (d2 == o.d2 && j < o.j); }
    };
    Near best[4];
    int keep = 0;
    int shell = 0;
    const double xi = scratch.xyz[static_cast<std::size_t>(3 * i)];
    const double yi = scratch.xyz[static_cast<std::size_t>(3 * i) + 1];
    const double zi = scratch.xyz[static_cast<std::size_t>(3 * i) + 2];
    const bool live = uniform || yCloud.pts[static_cast<std::size_t>(i)].type == typeI;
    if (csrOk && live) {
      const int a = scratch.offsets[static_cast<std::size_t>(i)];
      const int b = scratch.offsets[static_cast<std::size_t>(i) + 1];
      for (int p = a; p < b; p++) {
        const int j = scratch.cols[static_cast<std::size_t>(p)];
        if (j == i) {
          continue;
        }
        double rx, ry, rz;
        if (tiltBox) {
          const auto image = gen::fracDelta(
              fracBox, xi, yi, zi, scratch.xyz[static_cast<std::size_t>(3 * j)],
              scratch.xyz[static_cast<std::size_t>(3 * j) + 1],
              scratch.xyz[static_cast<std::size_t>(3 * j) + 2]);
          rx = image[0];
          ry = image[1];
          rz = image[2];
        } else {
          rx = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
          ry = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
          rz = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
          if (rx > hx) rx -= Lx; else if (rx < -hx) rx += Lx;
          if (ry > hy) ry -= Ly; else if (ry < -hy) ry += Ly;
          if (rz > hz) rz -= Lz; else if (rz < -hz) rz += Lz;
        }
        const double d2 = rx * rx + ry * ry + rz * rz;
        if (!(d2 > 0.0)) {
          continue;
        }
        Near item{d2, j, rx, ry, rz};
        ++shell;
        if (keep < 4) {
          best[keep++] = item;
          if (keep == 4) {
            std::sort(best, best + 4);
          }
        } else if (item < best[3]) {
          best[3] = item;
          std::sort(best, best + 4);
        }
      }
    } else if (live) {
      const int ix = std::min(nx - 1, static_cast<int>(xi * invx));
      const int iy = std::min(ny - 1, static_cast<int>(yi * invy));
      const int iz = std::min(nz - 1, static_cast<int>(zi * invz));
      for (int dz = -1; dz <= 1; dz++) {
        for (int dy = -1; dy <= 1; dy++) {
          for (int dx = -1; dx <= 1; dx++) {
            int jx = ix + dx;
            int jy = iy + dy;
            int jz = iz + dz;
            if (jx < 0) jx += nx; else if (jx >= nx) jx -= nx;
            if (jy < 0) jy += ny; else if (jy >= ny) jy -= ny;
            if (jz < 0) jz += nz; else if (jz >= nz) jz -= nz;
            for (int j = scratch.head[static_cast<std::size_t>((jz * ny + jy) * nx + jx)]; j >= 0;
                 j = scratch.next[static_cast<std::size_t>(j)]) {
              if (j == i) {
                continue;
              }
              double rx = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
              double ry = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
              double rz = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
              if (rx > hx) rx -= Lx; else if (rx < -hx) rx += Lx;
              if (ry > hy) ry -= Ly; else if (ry < -hy) ry += Ly;
              if (rz > hz) rz -= Lz; else if (rz < -hz) rz += Lz;
              const double d2 = rx * rx + ry * ry + rz * rz;
              if (!(d2 < cut2 && d2 > 0.0)) {
                continue;
              }
              Near item{d2, j, rx, ry, rz};
              ++shell;
              if (keep < 4) {
                best[keep++] = item;
                if (keep == 4) {
                  std::sort(best, best + 4);
                }
              } else if (item < best[3]) {
                best[3] = item;
                std::sort(best, best + 4);
              }
            }
          }
        }
      }
    }
    if (keep > 1 && keep < 4) {
      std::sort(best, best + keep);
    }
    scratch.nkeep[static_cast<std::size_t>(i)] = keep + (shell << 8);
    double xs[4], ys[4], zs[4];
    for (int s = 0; s < keep; s++) {
      scratch.chosen[static_cast<std::size_t>(i) * 4 + s] = best[s].j;
      xs[s] = best[s].x;
      ys[s] = best[s].y;
      zs[s] = best[s].z;
    }
    double acc[kRealL3] = {};
    if (keep > 0) {
      if (useAvx) {
        realYlm3Accumulate(xs, ys, zs, keep, acc);
      } else {
        for (int s = 0; s < keep; s++) {
          addRealYlm3Scalar(xs[s], ys[s], zs[s], norm, acc);
        }
      }
    }
    double *dst = scratch.q.data() + static_cast<std::size_t>(i) * kRealL3;
    for (int c = 0; c < kRealL3; c++) {
      dst[c] = acc[c];
    }
    const double hnorm = useAvx ? realL3Dot(acc, acc) : realL3DotScalar(acc, acc);
    scratch.qn[static_cast<std::size_t>(i)] = hnorm;
    scratch.nkeep[static_cast<std::size_t>(i)] = keep + (shell << 8);
  }
  const auto tYlm = std::chrono::steady_clock::now();
  std::vector<std::vector<molSys::Result>> saved;
  if (certify) {
    saved.resize(static_cast<std::size_t>(n));
  }
  std::vector<int> shellCount(static_cast<std::size_t>(n));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
  for (int i = 0; i < n; i++) {
    const int packed = scratch.nkeep[static_cast<std::size_t>(i)];
    const int keep = packed & 255;
    const int shell = packed >> 8;
    shellCount[static_cast<std::size_t>(i)] = shell;
    scratch.nkeep[static_cast<std::size_t>(i)] = keep;
    auto &bonds = certify ? saved[static_cast<std::size_t>(i)]
                          : yCloud.pts[static_cast<std::size_t>(i)].c_ij;
    bonds.clear();
    bonds.reserve(static_cast<std::size_t>(keep));
    const double *qi = scratch.q.data() + static_cast<std::size_t>(i) * kRealL3;
    const double ni = scratch.qn[static_cast<std::size_t>(i)];
    for (int s = 0; s < keep; s++) {
      const int j = scratch.chosen[static_cast<std::size_t>(i) * 4 + s];
      const double *qj = scratch.q.data() + static_cast<std::size_t>(j) * kRealL3;
      const double dot = useAvx ? realL3Dot(qi, qj) : realL3DotScalar(qi, qj);
      molSys::Result bond;
      bond.c_value = correlationOf(dot, ni, scratch.qn[static_cast<std::size_t>(j)]);
      bond.classifier = chill::classifyBond(bond.c_value, rule);
      bonds.push_back(bond);
    }
  }
  if (certify) {
    for (int i = 0; i < n; i++) {
      const int keep = scratch.nkeep[static_cast<std::size_t>(i)];
      const int shell = shellCount[static_cast<std::size_t>(i)];
      const auto &row = saved[static_cast<std::size_t>(i)];
      if (keep != 4 || shell != 4 || row.size() != 4) {
        return false;
      }
      double sum = 0.0;
      for (const auto &bond : row) {
        sum += bond.c_value;
      }
      const double cmean = 0.25 * sum;
      const double score = (cmean * cmean) * (cmean * cmean);
      if (!(cmean < 0.0 && score > kSymbolicCubicCut)) {
        return false;
      }
    }
    for (int i = 0; i < n; i++) {
      auto &point = yCloud.pts[static_cast<std::size_t>(i)];
      point.c_ij = std::move(saved[static_cast<std::size_t>(i)]);
      point.iceType = molSys::atom_state_type::cubic;
    }
    return true;
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
  for (int i = 0; i < n; i++) {
    auto &point = yCloud.pts[static_cast<std::size_t>(i)];
    const int keep = scratch.nkeep[static_cast<std::size_t>(i)];
    const int shell = shellCount[static_cast<std::size_t>(i)];
    int nstag = 0;
    int neclip = 0;
    for (const auto &bond : point.c_ij) {
      if (bond.classifier == molSys::bond_type::staggered) ++nstag;
      else if (bond.classifier == molSys::bond_type::eclipsed) ++neclip;
    }
    const bool four = shell == 4 && keep == 4;
    auto staggeredOf = [&](int j) {
      int count = 0;
      for (const auto &bond : yCloud.pts[static_cast<std::size_t>(j)].c_ij) {
        if (bond.classifier == molSys::bond_type::staggered) ++count;
      }
      return count;
    };
    bool inter = false;
    if (four && nstag == 2) {
      for (int s = 0; s < 4; s++) {
        if (staggeredOf(scratch.chosen[static_cast<std::size_t>(i) * 4 + s]) > 2) inter = true;
      }
    }
    if (four && nstag == 3 && neclip == 0) {
      for (int s = 0; s < 4; s++) {
        if (staggeredOf(scratch.chosen[static_cast<std::size_t>(i) * 4 + s]) > 1) inter = true;
      }
    }
    if (!four) point.iceType = molSys::atom_state_type::water;
    else if (neclip == 0 && nstag == 4) point.iceType = molSys::atom_state_type::cubic;
    else if (neclip == 1 && nstag == 3) point.iceType = molSys::atom_state_type::hexagonal;
    else if (inter) point.iceType = molSys::atom_state_type::interfacial;
    else if (neclip == 4 && nstag == 0) point.iceType = molSys::atom_state_type::clathrate;
    else if (neclip == 3) point.iceType = molSys::atom_state_type::interClathrate;
    else point.iceType = molSys::atom_state_type::water;
  }
  if (timing) {
    const auto tEnd = std::chrono::steady_clock::now();
    const auto ms = [](auto a, auto b) {
      return std::chrono::duration<double, std::milli>(b - a).count();
    };
    std::cerr << "chill-stages csr_ms=" << ms(tCsr0, tCsr1)
              << " walk_ylm_ms=" << ms(tWalk, tYlm)
              << " label_ms=" << ms(tYlm, tEnd) << "\n";
  }
  return true;
}

bool chill::q4q6ChillPlus(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
    int typeI, std::vector<double> &q4, std::vector<double> &q6,
    std::vector<double> &q4bar, std::vector<double> &q6bar) {
  const int n = yCloud.nop;
  if (n <= 0 || !(cutoff > 0.0) ||
      !(yCloud.box.size() == 3 || yCloud.box.size() >= 6)) {
    return false;
  }
  for (int i = 0; i < n; ++i) {
    if (yCloud.pts[static_cast<std::size_t>(i)].type != typeI) {
      return false;
    }
  }
  q4.assign(static_cast<std::size_t>(n), 0.0);
  q6.assign(static_cast<std::size_t>(n), 0.0);
  q4bar.assign(static_cast<std::size_t>(n), 0.0);
  q6bar.assign(static_cast<std::size_t>(n), 0.0);
  seams::prepareCompute();
  const bool tilt = yCloud.box.size() >= 6;
  const gen::FracBox frac = gen::makeFracBox(yCloud);
  const double Lx = tilt ? frac.lx : yCloud.box[0];
  const double Ly = tilt ? frac.ly : yCloud.box[1];
  const double Lz = tilt ? frac.lz : yCloud.box[2];
  const double ox = yCloud.boxLow.size() > 0 ? yCloud.boxLow[0] : 0.0;
  const double oy = yCloud.boxLow.size() > 1 ? yCloud.boxLow[1] : 0.0;
  const double oz = yCloud.boxLow.size() > 2 ? yCloud.boxLow[2] : 0.0;
#ifdef SEAMS_HAS_OPENMP
  const int savedThreads = std::max(1, omp_get_max_threads());
  omp_set_num_threads(seams::teamFor(n));
  struct RestoreThreads {
    int saved;
    ~RestoreThreads() { omp_set_num_threads(saved); }
  } restoreThreads{savedThreads};
  const bool parallel = seams::teamFor(n) > 1;
#else
  const bool parallel = false;
#endif
  std::vector<int> offsets, cols;
  if (!nneigh::orthoFullCsr(yCloud, nullptr, cutoff, offsets, cols)) {
    return false;
  }
  std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
  std::vector<double> qlm4(static_cast<std::size_t>(n) * 18, 0.0);
  std::vector<double> qlm6(static_cast<std::size_t>(n) * 26, 0.0);
  std::vector<double> q3(static_cast<std::size_t>(n) * kRealL3, 0.0);
  std::vector<double> qn(static_cast<std::size_t>(n), 0.0);
  std::vector<int> chosen(static_cast<std::size_t>(n) * 4, -1);
  std::vector<int> nkeep(static_cast<std::size_t>(n), 0);
  std::vector<int> shellCount(static_cast<std::size_t>(n), 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
  for (int i = 0; i < n; ++i) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    if (tilt) {
      double sz = (p.z - frac.oz) / frac.lz;
      sz -= std::floor(sz);
      if (!(sz < 1.0)) sz = 0.0;
      double sy = (p.y - frac.oy - frac.yz * sz) / frac.ly;
      sy -= std::floor(sy);
      if (!(sy < 1.0)) sy = 0.0;
      double sx = (p.x - frac.ox - frac.xy * sy - frac.xz * sz) / frac.lx;
      sx -= std::floor(sx);
      if (!(sx < 1.0)) sx = 0.0;
      xyz[static_cast<std::size_t>(3 * i)] = sx * frac.lx;
      xyz[static_cast<std::size_t>(3 * i) + 1] = sy * frac.ly;
      xyz[static_cast<std::size_t>(3 * i) + 2] = sz * frac.lz;
    } else {
      auto wrap = [](double v, double origin, double L) {
        double t = (v - origin) / L;
        t -= std::floor(t);
        if (!(t < 1.0)) t = 0.0;
        return t * L;
      };
      xyz[static_cast<std::size_t>(3 * i)] = wrap(p.x, ox, Lx);
      xyz[static_cast<std::size_t>(3 * i) + 1] = wrap(p.y, oy, Ly);
      xyz[static_cast<std::size_t>(3 * i) + 2] = wrap(p.z, oz, Lz);
    }
  }
  const double hx = 0.5 * Lx;
  const double hy = 0.5 * Ly;
  const double hz = 0.5 * Lz;
  const double axy = tilt ? frac.xy / frac.ly : 0.0;
  const double axz = tilt ? frac.xz / frac.lz : 0.0;
  const double ayz = tilt ? frac.yz / frac.lz : 0.0;
  constexpr double pi = 3.14159265358979323846;
  const double pre4 = 4.0 * pi / 9.0;
  const double pre6 = 4.0 * pi / 13.0;
  const bool useAvx = hostHasAvx2Fma();
  const auto &norm = normQ3();
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
  for (int i = 0; i < n; ++i) {
    struct Near {
      double d2;
      int j;
      double x, y, z;
      bool operator<(const Near &o) const {
        return d2 < o.d2 || (d2 == o.d2 && j < o.j);
      }
    };
    Near best[4];
    int keep = 0;
    int shell = 0;
    int nUsed = 0;
    double acc4[18] = {};
    double acc6[26] = {};
    double xs[8], ys[8], zs[8];
    double re4[8 * 9], im4[8 * 9], re6[8 * 13], im6[8 * 13];
    int filled = 0;
    const auto flush = [&]() {
      if (filled == 0) return;
      seams::ylm::ylmBatch(4, xs, ys, zs, filled, re4, im4);
      seams::ylm::ylmBatch(6, xs, ys, zs, filled, re6, im6);
      for (int bond = 0; bond < filled; ++bond) {
        const double *r4 = re4 + bond * 9;
        const double *i4 = im4 + bond * 9;
        const double *r6 = re6 + bond * 13;
        const double *i6 = im6 + bond * 13;
        for (int m = 0; m < 9; ++m) {
          acc4[2 * m] += r4[m];
          acc4[2 * m + 1] += i4[m];
        }
        for (int m = 0; m < 13; ++m) {
          acc6[2 * m] += r6[m];
          acc6[2 * m + 1] += i6[m];
        }
      }
      filled = 0;
    };
    const double xi = xyz[static_cast<std::size_t>(3 * i)];
    const double yi = xyz[static_cast<std::size_t>(3 * i) + 1];
    const double zi = xyz[static_cast<std::size_t>(3 * i) + 2];
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    for (int p = a; p < b; ++p) {
      const int j = cols[static_cast<std::size_t>(p)];
      double dx, dy, dz;
      if (tilt) {
        double dX = xi - xyz[static_cast<std::size_t>(3 * j)];
        double dY = yi - xyz[static_cast<std::size_t>(3 * j) + 1];
        double dZ = zi - xyz[static_cast<std::size_t>(3 * j) + 2];
        dX -= Lx * std::round(dX / Lx);
        dY -= Ly * std::round(dY / Ly);
        dZ -= Lz * std::round(dZ / Lz);
        dx = dX + axy * dY + axz * dZ;
        dy = dY + ayz * dZ;
        dz = dZ;
      } else {
        dx = xi - xyz[static_cast<std::size_t>(3 * j)];
        dy = yi - xyz[static_cast<std::size_t>(3 * j) + 1];
        dz = zi - xyz[static_cast<std::size_t>(3 * j) + 2];
        if (dx > hx) dx -= Lx;
        else if (dx < -hx) dx += Lx;
        if (dy > hy) dy -= Ly;
        else if (dy < -hy) dy += Ly;
        if (dz > hz) dz -= Lz;
        else if (dz < -hz) dz += Lz;
      }
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (!(d2 > 0.0)) continue;
      xs[filled] = dx;
      ys[filled] = dy;
      zs[filled] = dz;
      ++nUsed;
      if (++filled == 8) flush();
      Near item{d2, j, dx, dy, dz};
      ++shell;
      if (keep < 4) {
        best[keep++] = item;
        if (keep == 4) std::sort(best, best + 4);
      } else if (item < best[3]) {
        best[3] = item;
        std::sort(best, best + 4);
      }
    }
    flush();
    if (keep > 1 && keep < 4) std::sort(best, best + keep);
    shellCount[static_cast<std::size_t>(i)] = shell;
    nkeep[static_cast<std::size_t>(i)] = keep;
    if (nUsed > 0) {
      const double inv = 1.0 / static_cast<double>(nUsed);
      double sum4 = 0.0;
      double sum6 = 0.0;
      for (int m = 0; m < 9; ++m) {
        acc4[2 * m] *= inv;
        acc4[2 * m + 1] *= inv;
        sum4 += acc4[2 * m] * acc4[2 * m] + acc4[2 * m + 1] * acc4[2 * m + 1];
      }
      for (int m = 0; m < 13; ++m) {
        acc6[2 * m] *= inv;
        acc6[2 * m + 1] *= inv;
        sum6 += acc6[2 * m] * acc6[2 * m] + acc6[2 * m + 1] * acc6[2 * m + 1];
      }
      q4[static_cast<std::size_t>(i)] = std::sqrt(pre4 * sum4);
      q6[static_cast<std::size_t>(i)] = std::sqrt(pre6 * sum6);
      std::copy_n(acc4, 18, qlm4.data() + static_cast<std::size_t>(i) * 18);
      std::copy_n(acc6, 26, qlm6.data() + static_cast<std::size_t>(i) * 26);
    }
    double bx[4], by[4], bz[4];
    for (int s = 0; s < keep; ++s) {
      chosen[static_cast<std::size_t>(i) * 4 + s] = best[s].j;
      bx[s] = best[s].x;
      by[s] = best[s].y;
      bz[s] = best[s].z;
    }
    double acc[kRealL3] = {};
    if (keep > 0) {
      if (useAvx) realYlm3Accumulate(bx, by, bz, keep, acc);
      else {
        for (int s = 0; s < keep; ++s) addRealYlm3Scalar(bx[s], by[s], bz[s], norm, acc);
      }
    }
    double *dst = q3.data() + static_cast<std::size_t>(i) * kRealL3;
    for (int c = 0; c < kRealL3; ++c) dst[c] = acc[c];
    qn[static_cast<std::size_t>(i)] = useAvx ? realL3Dot(acc, acc) : realL3DotScalar(acc, acc);
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
  for (int i = 0; i < n; ++i) {
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    const double inv = 1.0 / static_cast<double>(b - a + 1);
    auto average = [&](const std::vector<double> &qlm, int width, double pre,
                       std::vector<double> &out) {
      double acc[26];
      const double *self = qlm.data() + static_cast<std::size_t>(i) * width;
      for (int c = 0; c < width; ++c) acc[c] = self[c];
      for (int p = a; p < b; ++p) {
        const int j = cols[static_cast<std::size_t>(p)];
        const double *other = qlm.data() + static_cast<std::size_t>(j) * width;
        for (int c = 0; c < width; ++c) acc[c] += other[c];
      }
      double sum = 0.0;
      for (int c = 0; c < width; ++c) {
        acc[c] *= inv;
        sum += acc[c] * acc[c];
      }
      out[static_cast<std::size_t>(i)] = std::sqrt(pre * sum);
    };
    average(qlm4, 18, pre4, q4bar);
    average(qlm6, 26, pre6, q6bar);
  }
  const chill::BondClassifier rule = chill::chillPlusRule();
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
  for (int i = 0; i < n; ++i) {
    const int keep = nkeep[static_cast<std::size_t>(i)];
    auto &bonds = yCloud.pts[static_cast<std::size_t>(i)].c_ij;
    bonds.clear();
    bonds.reserve(static_cast<std::size_t>(keep));
    const double *qi = q3.data() + static_cast<std::size_t>(i) * kRealL3;
    const double ni = qn[static_cast<std::size_t>(i)];
    for (int s = 0; s < keep; ++s) {
      const int j = chosen[static_cast<std::size_t>(i) * 4 + s];
      const double *qj = q3.data() + static_cast<std::size_t>(j) * kRealL3;
      const double dot = useAvx ? realL3Dot(qi, qj) : realL3DotScalar(qi, qj);
      molSys::Result bond;
      bond.c_value = correlationOf(dot, ni, qn[static_cast<std::size_t>(j)]);
      bond.classifier = chill::classifyBond(bond.c_value, rule);
      bonds.push_back(bond);
    }
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
  for (int i = 0; i < n; ++i) {
    auto &point = yCloud.pts[static_cast<std::size_t>(i)];
    const int keep = nkeep[static_cast<std::size_t>(i)];
    const int shell = shellCount[static_cast<std::size_t>(i)];
    int nstag = 0;
    int neclip = 0;
    for (const auto &bond : point.c_ij) {
      if (bond.classifier == molSys::bond_type::staggered) ++nstag;
      else if (bond.classifier == molSys::bond_type::eclipsed) ++neclip;
    }
    const bool four = shell == 4 && keep == 4;
    auto staggeredOf = [&](int j) {
      int count = 0;
      for (const auto &bond : yCloud.pts[static_cast<std::size_t>(j)].c_ij) {
        if (bond.classifier == molSys::bond_type::staggered) ++count;
      }
      return count;
    };
    bool inter = false;
    if (four && nstag == 2) {
      for (int s = 0; s < 4; ++s) {
        if (staggeredOf(chosen[static_cast<std::size_t>(i) * 4 + s]) > 2) inter = true;
      }
    }
    if (four && nstag == 3 && neclip == 0) {
      for (int s = 0; s < 4; ++s) {
        if (staggeredOf(chosen[static_cast<std::size_t>(i) * 4 + s]) > 1) inter = true;
      }
    }
    if (!four) point.iceType = molSys::atom_state_type::water;
    else if (neclip == 0 && nstag == 4) point.iceType = molSys::atom_state_type::cubic;
    else if (neclip == 1 && nstag == 3) point.iceType = molSys::atom_state_type::hexagonal;
    else if (inter) point.iceType = molSys::atom_state_type::interfacial;
    else if (neclip == 4 && nstag == 0) point.iceType = molSys::atom_state_type::clathrate;
    else if (neclip == 3) point.iceType = molSys::atom_state_type::interClathrate;
    else point.iceType = molSys::atom_state_type::water;
  }
  return true;
}

bool chill::chillPlusSymbolicPrefilter(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
    int typeI) {
  if (!(yCloud.nop > 0 && cutoff > 0.0 && yCloud.box.size() >= 3 &&
        yCloud.box.size() < 6)) {
    return false;
  }
  return orthoChill(yCloud, cutoff, typeI, true);
}

bool chill::chillPlusPublishedCutoff(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
    int typeI) {
  if (yCloud.nop > 0 && cutoff > 0.0 &&
      (yCloud.box.size() == 3 || yCloud.box.size() >= 6)) {
    return orthoChill(yCloud, cutoff, typeI, false);
  }
#ifndef SEAMS_HAS_LINKCELL
  (void)yCloud;
  (void)cutoff;
  (void)typeI;
  return false;
#else
  const int n = yCloud.nop;
  if (n <= 0 || !(cutoff > 0.0) || yCloud.box.size() < 3) {
    return false;
  }
  seams::prepareCompute();
  [[maybe_unused]] const seams::TeamWidth team(n);
#ifdef SEAMS_HAS_LINKCELL
  lc_set_threads(seams::teamFor(n));
#endif
  const bool timing = std::getenv("SEAMS_CHILL_TIME") != nullptr;
  const auto t0 = std::chrono::steady_clock::now();
  struct Scratch {
    std::vector<double> xyz;
    std::vector<int> mask;
    std::vector<int> ii, jj, shift;
    std::vector<double> d2;
    std::vector<int> degree, offsets, cols;
    std::vector<int> shell, nkeep, chosen;
    std::vector<double> q, qn;
    int n = -1;
  };
  static Scratch scratch;
  if (scratch.n != n) {
    scratch.xyz.resize(static_cast<std::size_t>(n) * 3);
    scratch.degree.assign(static_cast<std::size_t>(n), 0);
    scratch.shell.assign(static_cast<std::size_t>(n), 0);
    scratch.nkeep.assign(static_cast<std::size_t>(n), 0);
    scratch.chosen.assign(static_cast<std::size_t>(n) * 4, -1);
    scratch.q.assign(static_cast<std::size_t>(n) * kRealL3, 0.0);
    scratch.qn.assign(static_cast<std::size_t>(n), 0.0);
    scratch.n = n;
  } else {
    std::fill(scratch.degree.begin(), scratch.degree.end(), 0);
  }
  bool uniform = true;
  for (int i = 0; i < n; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    scratch.xyz[static_cast<std::size_t>(3 * i)] = p.x;
    scratch.xyz[static_cast<std::size_t>(3 * i) + 1] = p.y;
    scratch.xyz[static_cast<std::size_t>(3 * i) + 2] = p.z;
    if (p.type != typeI) {
      uniform = false;
    }
    yCloud.pts[static_cast<std::size_t>(i)].c_ij.clear();
  }
  const int *maskPtr = nullptr;
  if (!uniform) {
    scratch.mask.assign(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; i++) {
      if (yCloud.pts[static_cast<std::size_t>(i)].type == typeI) {
        scratch.mask[static_cast<std::size_t>(i)] = 1;
      }
    }
    maskPtr = scratch.mask.data();
  }
  const auto tPairs = std::chrono::steady_clock::now();
  // The in-call cell list is slower than the existing pair search on this
  // crystal. Tilted boxes and the ortho cutoff both use that search.
  const bool ortho = false;
  if (ortho) {
    const double Lx = yCloud.box[0];
    const double Ly = yCloud.box[1];
    const double Lz = yCloud.box[2];
    const int nx = std::max(1, static_cast<int>(std::ceil(Lx / cutoff)));
    const int ny = std::max(1, static_cast<int>(std::ceil(Ly / cutoff)));
    const int nz = std::max(1, static_cast<int>(std::ceil(Lz / cutoff)));
    const int ncell = nx * ny * nz;
    std::vector<int> head(static_cast<std::size_t>(ncell), -1);
    std::vector<int> next(static_cast<std::size_t>(n), -1);
    const double invx = static_cast<double>(nx) / Lx;
    const double invy = static_cast<double>(ny) / Ly;
    const double invz = static_cast<double>(nz) / Lz;
    auto wrapped = [](double v, double L) {
      double t = v / L;
      t -= std::floor(t);
      return t * L;
    };
    for (int i = 0; i < n; i++) {
      if (maskPtr != nullptr && maskPtr[i] == 0) {
        continue;
      }
      const double x = wrapped(scratch.xyz[static_cast<std::size_t>(3 * i)], Lx);
      const double y = wrapped(scratch.xyz[static_cast<std::size_t>(3 * i) + 1], Ly);
      const double z = wrapped(scratch.xyz[static_cast<std::size_t>(3 * i) + 2], Lz);
      const int ix = std::min(nx - 1, static_cast<int>(x * invx));
      const int iy = std::min(ny - 1, static_cast<int>(y * invy));
      const int iz = std::min(nz - 1, static_cast<int>(z * invz));
      const int c = (iz * ny + iy) * nx + ix;
      next[static_cast<std::size_t>(i)] = head[static_cast<std::size_t>(c)];
      head[static_cast<std::size_t>(c)] = i;
    }
    const double cut2 = cutoff * cutoff;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
    for (int i = 0; i < n; i++) {
      if (maskPtr != nullptr && maskPtr[i] == 0) {
        continue;
      }
      const double xi = scratch.xyz[static_cast<std::size_t>(3 * i)];
      const double yi = scratch.xyz[static_cast<std::size_t>(3 * i) + 1];
      const double zi = scratch.xyz[static_cast<std::size_t>(3 * i) + 2];
      const double x = wrapped(xi, Lx);
      const double y = wrapped(yi, Ly);
      const double z = wrapped(zi, Lz);
      const int ix = std::min(nx - 1, static_cast<int>(x * invx));
      const int iy = std::min(ny - 1, static_cast<int>(y * invy));
      const int iz = std::min(nz - 1, static_cast<int>(z * invz));
      int count = 0;
      for (int dz = -1; dz <= 1; dz++) {
        for (int dy = -1; dy <= 1; dy++) {
          for (int dx = -1; dx <= 1; dx++) {
            int jx = ix + dx;
            int jy = iy + dy;
            int jz = iz + dz;
            jx %= nx;
            if (jx < 0) {
              jx += nx;
            }
            jy %= ny;
            if (jy < 0) {
              jy += ny;
            }
            jz %= nz;
            if (jz < 0) {
              jz += nz;
            }
            for (int j = head[static_cast<std::size_t>((jz * ny + jy) * nx + jx)]; j >= 0;
                 j = next[static_cast<std::size_t>(j)]) {
              if (j == i) {
                continue;
              }
              double rx = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
              double ry = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
              double rz = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
              rx -= Lx * std::floor(rx / Lx + 0.5);
              ry -= Ly * std::floor(ry / Ly + 0.5);
              rz -= Lz * std::floor(rz / Lz + 0.5);
              if (rx * rx + ry * ry + rz * rz < cut2) {
                ++count;
              }
            }
          }
        }
      }
      scratch.degree[static_cast<std::size_t>(i)] = count;
    }
  } else {
  const linkcell::Cell cell = nneigh::lammpsBoxToLinkcell(yCloud.box, yCloud.boxLow);
  const lc_cell raw = cell.raw();
  std::size_t count = 0;
  if (lc_pairs_within(scratch.xyz.data(), static_cast<std::size_t>(n), &raw, cutoff,
                      maskPtr, cutoff, 0, nullptr, nullptr, nullptr, nullptr, 0,
                      &count) != 0) {
    return false;
  }
  scratch.ii.resize(count);
  scratch.jj.resize(count);
  scratch.shift.resize(count * 3);
  scratch.d2.resize(count);
  std::size_t got = 0;
  if (lc_pairs_within(scratch.xyz.data(), static_cast<std::size_t>(n), &raw, cutoff,
                      maskPtr, cutoff, 0, scratch.ii.data(), scratch.jj.data(),
                      scratch.shift.data(), scratch.d2.data(), count, &got) != 0) {
    return false;
  }
  for (std::size_t t = 0; t < got; t++) {
    const int i = scratch.ii[t];
    const int j = scratch.jj[t];
    if (i >= 0 && j >= 0 && i < n && j < n && i != j) {
      scratch.degree[static_cast<std::size_t>(i)]++;
    }
  }
  }
  scratch.offsets.assign(static_cast<std::size_t>(n) + 1, 0);
  for (int i = 0; i < n; i++) {
    scratch.offsets[static_cast<std::size_t>(i) + 1] =
        scratch.offsets[static_cast<std::size_t>(i)] +
        scratch.degree[static_cast<std::size_t>(i)];
  }
  const int nnz = scratch.offsets[static_cast<std::size_t>(n)];
  scratch.cols.resize(static_cast<std::size_t>(nnz));
  if (ortho) {
    const double Lx = yCloud.box[0];
    const double Ly = yCloud.box[1];
    const double Lz = yCloud.box[2];
    const int nx = std::max(1, static_cast<int>(std::ceil(Lx / cutoff)));
    const int ny = std::max(1, static_cast<int>(std::ceil(Ly / cutoff)));
    const int nz = std::max(1, static_cast<int>(std::ceil(Lz / cutoff)));
    const int ncell = nx * ny * nz;
    std::vector<int> head(static_cast<std::size_t>(ncell), -1);
    std::vector<int> next(static_cast<std::size_t>(n), -1);
    const double invx = static_cast<double>(nx) / Lx;
    const double invy = static_cast<double>(ny) / Ly;
    const double invz = static_cast<double>(nz) / Lz;
    const double cut2 = cutoff * cutoff;
    auto wrapped = [](double v, double L) {
      double t = v / L;
      t -= std::floor(t);
      return t * L;
    };
    for (int i = 0; i < n; i++) {
      if (maskPtr != nullptr && maskPtr[i] == 0) {
        continue;
      }
      const double x = wrapped(scratch.xyz[static_cast<std::size_t>(3 * i)], Lx);
      const double y = wrapped(scratch.xyz[static_cast<std::size_t>(3 * i) + 1], Ly);
      const double z = wrapped(scratch.xyz[static_cast<std::size_t>(3 * i) + 2], Lz);
      const int ix = std::min(nx - 1, static_cast<int>(x * invx));
      const int iy = std::min(ny - 1, static_cast<int>(y * invy));
      const int iz = std::min(nz - 1, static_cast<int>(z * invz));
      const int c = (iz * ny + iy) * nx + ix;
      next[static_cast<std::size_t>(i)] = head[static_cast<std::size_t>(c)];
      head[static_cast<std::size_t>(c)] = i;
    }
    std::vector<int> cursor = scratch.offsets;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
    for (int i = 0; i < n; i++) {
      if (maskPtr != nullptr && maskPtr[i] == 0) {
        continue;
      }
      const double xi = scratch.xyz[static_cast<std::size_t>(3 * i)];
      const double yi = scratch.xyz[static_cast<std::size_t>(3 * i) + 1];
      const double zi = scratch.xyz[static_cast<std::size_t>(3 * i) + 2];
      const double x = wrapped(xi, Lx);
      const double y = wrapped(yi, Ly);
      const double z = wrapped(zi, Lz);
      const int ix = std::min(nx - 1, static_cast<int>(x * invx));
      const int iy = std::min(ny - 1, static_cast<int>(y * invy));
      const int iz = std::min(nz - 1, static_cast<int>(z * invz));
      int slot = cursor[static_cast<std::size_t>(i)];
      for (int dz = -1; dz <= 1; dz++) {
        for (int dy = -1; dy <= 1; dy++) {
          for (int dx = -1; dx <= 1; dx++) {
            int jx = ix + dx;
            int jy = iy + dy;
            int jz = iz + dz;
            jx %= nx;
            if (jx < 0) {
              jx += nx;
            }
            jy %= ny;
            if (jy < 0) {
              jy += ny;
            }
            jz %= nz;
            if (jz < 0) {
              jz += nz;
            }
            for (int j = head[static_cast<std::size_t>((jz * ny + jy) * nx + jx)]; j >= 0;
                 j = next[static_cast<std::size_t>(j)]) {
              if (j == i) {
                continue;
              }
              double rx = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
              double ry = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
              double rz = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
              rx -= Lx * std::floor(rx / Lx + 0.5);
              ry -= Ly * std::floor(ry / Ly + 0.5);
              rz -= Lz * std::floor(rz / Lz + 0.5);
              if (rx * rx + ry * ry + rz * rz < cut2) {
                scratch.cols[static_cast<std::size_t>(slot++)] = j;
              }
            }
          }
        }
      }
    }
  } else {
  std::vector<int> cursor = scratch.offsets;
  for (std::size_t t = 0; t < scratch.ii.size(); t++) {
    const int i = scratch.ii[t];
    const int j = scratch.jj[t];
    if (i >= 0 && j >= 0 && i < n && j < n && i != j) {
      scratch.cols[static_cast<std::size_t>(cursor[static_cast<std::size_t>(i)]++)] = j;
    }
  }
  }
  const auto tCsr = std::chrono::steady_clock::now();
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  const auto &norm = normQ3();
  const bool useAvx = hostHasAvx2Fma();
  const chill::BondClassifier rule = chill::chillPlusRule();
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
  for (int i = 0; i < n; i++) {
    const int begin = scratch.offsets[static_cast<std::size_t>(i)];
    const int end = scratch.offsets[static_cast<std::size_t>(i) + 1];
    struct Near {
      double d2;
      int j;
      double x, y, z;
      bool operator<(const Near &o) const {
        return d2 < o.d2 || (d2 == o.d2 && j < o.j);
      }
    };
    Near stackNear[32];
    std::vector<Near> heapNear;
    Near *near = stackNear;
    int nNear = 0;
    const double xi = scratch.xyz[static_cast<std::size_t>(3 * i)];
    const double yi = scratch.xyz[static_cast<std::size_t>(3 * i) + 1];
    const double zi = scratch.xyz[static_cast<std::size_t>(3 * i) + 2];
    for (int t = begin; t < end; t++) {
      const int j = scratch.cols[static_cast<std::size_t>(t)];
      bool seen = false;
      for (int s = 0; s < nNear; s++) {
        if (near[s].j == j) {
          seen = true;
          break;
        }
      }
      if (seen) {
        continue;
      }
      const auto dr = gen::fracDelta(frame, xi, yi, zi,
                                     scratch.xyz[static_cast<std::size_t>(3 * j)],
                                     scratch.xyz[static_cast<std::size_t>(3 * j) + 1],
                                     scratch.xyz[static_cast<std::size_t>(3 * j) + 2]);
      if (nNear == 32 && near == stackNear) {
        heapNear.assign(stackNear, stackNear + 32);
        near = heapNear.data();
      }
      Near item{dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2], j, dr[0], dr[1], dr[2]};
      if (near == stackNear) {
        near[nNear] = item;
      } else {
        heapNear.push_back(item);
        near = heapNear.data();
      }
      ++nNear;
    }
    const int shell = nNear;
    scratch.shell[static_cast<std::size_t>(i)] = shell;
    const int keep = std::min(4, shell);
    if (keep > 0) {
      std::partial_sort(near, near + keep, near + shell);
    }
    scratch.nkeep[static_cast<std::size_t>(i)] = keep;
    double xs[4], ys[4], zs[4];
    for (int s = 0; s < keep; s++) {
      scratch.chosen[static_cast<std::size_t>(i) * 4 + static_cast<std::size_t>(s)] =
          near[static_cast<std::size_t>(s)].j;
      xs[s] = near[static_cast<std::size_t>(s)].x;
      ys[s] = near[static_cast<std::size_t>(s)].y;
      zs[s] = near[static_cast<std::size_t>(s)].z;
    }
    double acc[kRealL3] = {};
    if (keep > 0) {
      if (useAvx) {
        realYlm3Accumulate(xs, ys, zs, keep, acc);
      } else {
        for (int s = 0; s < keep; s++) {
          addRealYlm3Scalar(xs[s], ys[s], zs[s], norm, acc);
        }
      }
    }
    double *dst = scratch.q.data() + static_cast<std::size_t>(i) * kRealL3;
    for (int c = 0; c < kRealL3; c++) {
      dst[c] = acc[c];
    }
    scratch.qn[static_cast<std::size_t>(i)] =
        useAvx ? realL3Dot(acc, acc) : realL3DotScalar(acc, acc);
  }
  const auto tY = std::chrono::steady_clock::now();
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
  for (int i = 0; i < n; i++) {
    const int keep = scratch.nkeep[static_cast<std::size_t>(i)];
    auto &bonds = yCloud.pts[static_cast<std::size_t>(i)].c_ij;
    bonds.clear();
    bonds.reserve(static_cast<std::size_t>(keep));
    const double *qi = scratch.q.data() + static_cast<std::size_t>(i) * kRealL3;
    const double ni = scratch.qn[static_cast<std::size_t>(i)];
    for (int s = 0; s < keep; s++) {
      const int j =
          scratch.chosen[static_cast<std::size_t>(i) * 4 + static_cast<std::size_t>(s)];
      const double *qj = scratch.q.data() + static_cast<std::size_t>(j) * kRealL3;
      const double dot = useAvx ? realL3Dot(qi, qj) : realL3DotScalar(qi, qj);
      molSys::Result bond;
      bond.c_value = correlationOf(dot, ni, scratch.qn[static_cast<std::size_t>(j)]);
      bond.classifier = chill::classifyBond(bond.c_value, rule);
      bonds.push_back(bond);
    }
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
  for (int i = 0; i < n; i++) {
    auto &point = yCloud.pts[static_cast<std::size_t>(i)];
    const int keep = scratch.nkeep[static_cast<std::size_t>(i)];
    const int shell = scratch.shell[static_cast<std::size_t>(i)];
    int nstag = 0;
    int neclip = 0;
    for (const auto &bond : point.c_ij) {
      if (bond.classifier == molSys::bond_type::staggered) {
        ++nstag;
      } else if (bond.classifier == molSys::bond_type::eclipsed) {
        ++neclip;
      }
    }
    const bool four = shell == 4 && keep == 4;
    auto staggeredOf = [&](int j) {
      int count = 0;
      for (const auto &bond : yCloud.pts[static_cast<std::size_t>(j)].c_ij) {
        if (bond.classifier == molSys::bond_type::staggered) {
          ++count;
        }
      }
      return count;
    };
    bool inter = false;
    if (four && nstag == 2) {
      for (int s = 0; s < 4; s++) {
        const int j =
            scratch.chosen[static_cast<std::size_t>(i) * 4 + static_cast<std::size_t>(s)];
        if (staggeredOf(j) > 2) {
          inter = true;
        }
      }
    }
    if (four && nstag == 3 && neclip == 0) {
      for (int s = 0; s < 4; s++) {
        const int j =
            scratch.chosen[static_cast<std::size_t>(i) * 4 + static_cast<std::size_t>(s)];
        if (staggeredOf(j) > 1) {
          inter = true;
        }
      }
    }
    if (!four) {
      point.iceType = molSys::atom_state_type::water;
    } else if (neclip == 0 && nstag == 4) {
      point.iceType = molSys::atom_state_type::cubic;
    } else if (neclip == 1 && nstag == 3) {
      point.iceType = molSys::atom_state_type::hexagonal;
    } else if (inter) {
      point.iceType = molSys::atom_state_type::interfacial;
    } else if (neclip == 4 && nstag == 0) {
      point.iceType = molSys::atom_state_type::clathrate;
    } else if (neclip == 3) {
      point.iceType = molSys::atom_state_type::interClathrate;
    } else {
      point.iceType = molSys::atom_state_type::water;
    }
  }
  if (timing) {
    const auto tEnd = std::chrono::steady_clock::now();
    const auto ms = [](auto a, auto b) {
      return std::chrono::duration<double, std::milli>(b - a).count();
    };
    std::cerr << "chill-stages n=" << n << " pairs_ms=" << ms(t0, tPairs)
              << " csr_ms=" << ms(tPairs, tCsr) << " ylm_ms=" << ms(tCsr, tY)
              << " label_ms=" << ms(tY, tEnd) << "\n";
  }
  return true;
#endif
}

/**
 *  @details Function for getting the bond order correlations @f$c_{ij}@f$
 * (alternatively
 *   @f$a_{ij}@f$ in certain texts) using the CHILL+ algorithm
 *  @param[in,out] yCloud The output molSys::PointCloud
 *  @param[in] nList Row-ordered neighbour list by atom ID
 *  @param[in] isSlice This decides whether there is a slice or not
 */
void chill::classifyBondsCsr(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &offsets, const std::vector<int> &cols,
    const chill::BondClassifier &rule, bool isSlice) {
  const std::vector<std::vector<int>> none;
  if (seams::product::reduction() == seams::product::Reduction::Left &&
      seams::product::mode() == seams::product::Mode::Fp64) {
    classifyBondsLeftFp64(yCloud, none, rule, isSlice, &offsets, &cols);
    return;
  }
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(yCloud.nop));
  for (int i = 0; i < yCloud.nop && i + 1 < static_cast<int>(offsets.size()); ++i) {
    rows[static_cast<std::size_t>(i)].push_back(
        yCloud.pts[static_cast<std::size_t>(i)].atomID);
    for (int k = offsets[static_cast<std::size_t>(i)];
         k < offsets[static_cast<std::size_t>(i) + 1]; ++k) {
      const int j = cols[static_cast<std::size_t>(k)];
      if (j >= 0 && j < yCloud.nop) {
        rows[static_cast<std::size_t>(i)].push_back(
            yCloud.pts[static_cast<std::size_t>(j)].atomID);
      }
    }
  }
  classifyBonds(yCloud, rows, rule, isSlice);
}

int chill::chillPlusFromCsr(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                            const std::vector<int> &offsets, const std::vector<int> &cols) {
  const int n = yCloud.nop;
  if (n <= 0 || static_cast<int>(offsets.size()) < n + 1) {
    return 0;
  }
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  const bool useAvx = hostHasAvx2Fma();
  const auto &norm = normQ3();
  const chill::BondClassifier rule = chill::chillPlusRule();
  std::vector<double> q(static_cast<std::size_t>(n) * kRealL3, 0.0);
  std::vector<double> qn(static_cast<std::size_t>(n), 0.0);
  std::vector<int> chosen(static_cast<std::size_t>(n) * 4, -1);
  std::vector<int> nkeep(static_cast<std::size_t>(n), 0);
  std::vector<int> shell(static_cast<std::size_t>(n), 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
  for (int i = 0; i < n; ++i) {
    struct Near {
      double d2;
      int j;
      double x, y, z;
      bool operator<(const Near &o) const { return d2 < o.d2 || (d2 == o.d2 && j < o.j); }
    };
    Near best[4];
    int keep = 0;
    int nShell = 0;
    const auto &pi = yCloud.pts[static_cast<std::size_t>(i)];
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    for (int k = a; k < b; ++k) {
      const int j = cols[static_cast<std::size_t>(k)];
      if (j < 0 || j >= n || j == i) {
        continue;
      }
      const auto &pj = yCloud.pts[static_cast<std::size_t>(j)];
      const auto dr = gen::fracDelta(frame, pi.x, pi.y, pi.z, pj.x, pj.y, pj.z);
      const double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
      if (!(d2 > 0.0)) {
        continue;
      }
      ++nShell;
      Near item{d2, j, dr[0], dr[1], dr[2]};
      if (keep < 4) {
        best[keep++] = item;
        if (keep == 4) {
          std::sort(best, best + 4);
        }
      } else if (item < best[3]) {
        best[3] = item;
        std::sort(best, best + 4);
      }
    }
    if (keep > 1 && keep < 4) {
      std::sort(best, best + keep);
    }
    shell[static_cast<std::size_t>(i)] = nShell;
    nkeep[static_cast<std::size_t>(i)] = keep;
    double xs[4], ys[4], zs[4];
    for (int s = 0; s < keep; ++s) {
      chosen[static_cast<std::size_t>(i) * 4 + static_cast<std::size_t>(s)] = best[s].j;
      xs[s] = best[s].x;
      ys[s] = best[s].y;
      zs[s] = best[s].z;
    }
    double acc[kRealL3] = {};
    if (keep > 0) {
      if (useAvx) {
        realYlm3Accumulate(xs, ys, zs, keep, acc);
      } else {
        for (int s = 0; s < keep; ++s) {
          addRealYlm3Scalar(xs[s], ys[s], zs[s], norm, acc);
        }
      }
    }
    double *dst = q.data() + static_cast<std::size_t>(i) * kRealL3;
    for (int c = 0; c < kRealL3; ++c) {
      dst[c] = acc[c];
    }
    qn[static_cast<std::size_t>(i)] = useAvx ? realL3Dot(acc, acc) : realL3DotScalar(acc, acc);
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
  for (int i = 0; i < n; ++i) {
    const int keep = nkeep[static_cast<std::size_t>(i)];
    auto &bonds = yCloud.pts[static_cast<std::size_t>(i)].c_ij;
    bonds.clear();
    bonds.reserve(static_cast<std::size_t>(keep));
    const double *qi = q.data() + static_cast<std::size_t>(i) * kRealL3;
    const double ni = qn[static_cast<std::size_t>(i)];
    for (int s = 0; s < keep; ++s) {
      const int j = chosen[static_cast<std::size_t>(i) * 4 + static_cast<std::size_t>(s)];
      const double *qj = q.data() + static_cast<std::size_t>(j) * kRealL3;
      molSys::Result bond;
      bond.c_value = correlationOf(useAvx ? realL3Dot(qi, qj) : realL3DotScalar(qi, qj), ni,
                                   qn[static_cast<std::size_t>(j)]);
      bond.classifier = chill::classifyBond(bond.c_value, rule);
      bonds.push_back(bond);
    }
  }
  // Bonds are complete before the neighbour test reads them. The published
  // interfacial row looks at a neighbour's staggered count; labelling in
  // the bond loop would see a neighbour that had not been written yet.
  int cubic = 0;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for reduction(+ : cubic) schedule(static) if (n >= 256)
#endif
  for (int i = 0; i < n; ++i) {
    const int keep = nkeep[static_cast<std::size_t>(i)];
    const auto &bonds = yCloud.pts[static_cast<std::size_t>(i)].c_ij;
    int nstag = 0;
    int neclip = 0;
    for (const auto &bond : bonds) {
      if (bond.classifier == molSys::bond_type::staggered) {
        ++nstag;
      } else if (bond.classifier == molSys::bond_type::eclipsed) {
        ++neclip;
      }
    }
    const bool four = shell[static_cast<std::size_t>(i)] == 4 && keep == 4;
    auto staggeredOf = [&](int j) {
      int count = 0;
      for (const auto &bond : yCloud.pts[static_cast<std::size_t>(j)].c_ij) {
        if (bond.classifier == molSys::bond_type::staggered) {
          ++count;
        }
      }
      return count;
    };
    bool inter = false;
    if (four && nstag == 2) {
      for (int s = 0; s < 4; ++s) {
        if (staggeredOf(chosen[static_cast<std::size_t>(i) * 4 + static_cast<std::size_t>(s)]) >
            2) {
          inter = true;
        }
      }
    }
    if (four && nstag == 3 && neclip == 0) {
      for (int s = 0; s < 4; ++s) {
        if (staggeredOf(chosen[static_cast<std::size_t>(i) * 4 + static_cast<std::size_t>(s)]) >
            1) {
          inter = true;
        }
      }
    }
    auto &ice = yCloud.pts[static_cast<std::size_t>(i)].iceType;
    if (!four) {
      ice = molSys::atom_state_type::water;
    } else if (neclip == 0 && nstag == 4) {
      ice = molSys::atom_state_type::cubic;
      ++cubic;
    } else if (neclip == 1 && nstag == 3) {
      ice = molSys::atom_state_type::hexagonal;
    } else if (inter) {
      ice = molSys::atom_state_type::interfacial;
    } else if (neclip == 4 && nstag == 0) {
      ice = molSys::atom_state_type::clathrate;
    } else if (neclip == 3) {
      ice = molSys::atom_state_type::interClathrate;
    } else {
      ice = molSys::atom_state_type::water;
    }
  }
  return cubic;
}

void
chill::getCorrelPlus(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     const std::vector<std::vector<int>> &nList, bool isSlice,
                     int coordinationNumber) {
  chill::getCorrelPlus(yCloud, nList, isSlice, coordinationNumber, true);
}

void
chill::getCorrelPlus(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     const std::vector<std::vector<int>> &nList, bool isSlice,
                     int coordinationNumber, bool closedBounds) {
  chill::BondClassifier rule = chill::chillPlusRule();
  rule.coordinationNumber = coordinationNumber;
  chill::classifyBonds(yCloud, nList, rule, isSlice, closedBounds);
}

/**
 *  @details Function that classifies the #molSys::atom_state_type ice type of
 * each particle, according to the CHILL+ algorithm.
 *  @param[in,out] yCloud The output molSys::PointCloud
 *  @param[in] nList Row-ordered neighbour list by atom ID
 *  @param[in] path Path to the output directory to which ice types are written
 *   out to
 *  @param[in] firstFrame The first frame to be analyzed
 *  @param[in] isSlice This decides whether there is a slice or not
 *  @param[in] outputFileName Name of the output file, to which the ice types
 *   will be written out.
 *   The default file name is "chillPlus.txt"
 */
namespace {
void assignIceTypePlus(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, bool isSlice, bool neighbourTest,
    int &ih, int &ic, int &interIce, int &water, int &clath, int &interClath,
    int &total) {
  ih = ic = water = interIce = total = 0;
  clath = interClath = 0;
  int num_staggrd, num_eclipsd, na;
  molSys::bond_type bondType;
  int nnumNeighbours;

#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1) reduction(+ : ih, ic, water, interIce, total, clath, interClath) private(num_staggrd, num_eclipsd, na, bondType, nnumNeighbours)
#endif
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // if(yCloud.pts[iatom].type!=typeO){continue;}
    // if this is a slice and the particle is not in the slice
    // then skip
    if (isSlice) {
      if (yCloud.pts[iatom].inSlice == false) {
        continue;
      }
    }
    total++; // Update the total number of atoms considered. Change this to a
             // check for slices
    // c_ij is the bond list after the nearest-neighbour cut. The crystal
    // test below counts the cutoff row itself, before that cut.
    nnumNeighbours = yCloud.pts[iatom].c_ij.size();
    num_staggrd = num_eclipsd = na =
        0; // init to zero before loop through neighbours
    // Loop through the bond cij and get the number of staggered, eclipsed bonds
    for (int j = 0; j < nnumNeighbours; j++) {
      bondType = yCloud.pts[iatom].c_ij[j].classifier;
      if (bondType == molSys::bond_type::eclipsed) {
        num_eclipsd++;
      } else if (bondType == molSys::bond_type::staggered) {
        num_staggrd++;
      } else {
        na++;
      }
    } // End of loop through neighbours

    // Add more tests later
    yCloud.pts[iatom].iceType = molSys::atom_state_type::unclassified; // default
    // Table 1 footnote a: a crystal row needs exactly four neighbours in
    // the cutoff shell. Four nearest bonds are not enough when a fifth
    // oxygen also sits inside the shell.
    const bool fourShell =
        chill::cutoffShellCount(yCloud, nList, iatom) == 4 &&
        nnumNeighbours == 4;
    if (fourShell) {
      // Cubic ice
      if (num_eclipsd == 0 && num_staggrd == 4) {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::cubic;
        ic++;
      }
      // Hexagonal
      else if (num_eclipsd == 1 && num_staggrd == 3) {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::hexagonal;
        ih++;
      }
      // Interfacial. The published rule applies Table 1's neighbour
      // test. ovito-chill+ labels two staggered bonds, or three
      // staggered and no eclipsed bond, with no look at the neighbour.
      else if (neighbourTest
                   ? chill::isInterfacial(yCloud, nList, iatom, num_staggrd,
                                          num_eclipsd, true)
                   : (num_staggrd == 2 ||
                      (num_staggrd == 3 && num_eclipsd == 0))) {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::interfacial;
        interIce++;
      }
      // Clathrate
      else if (num_eclipsd == 4 && num_staggrd == 0) {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::clathrate;
        clath++;
      }
      // Interfacial clathrate
      else if (num_eclipsd == 3) {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::interClathrate;
        interClath++;
      }
      // Water
      else {
        yCloud.pts[iatom].iceType = molSys::atom_state_type::water;
        water++;
      }
    } else {
      yCloud.pts[iatom].iceType = molSys::atom_state_type::water;
      water++;
    }

  } // End of loop through every iatom
}
} // namespace

void
chill::getIceTypePlusNoPrint(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, bool isSlice) {
  chill::getIceTypePlusNoPrint(yCloud, nList, isSlice, true);
}

void
chill::getIceTypePlusNoPrint(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, bool isSlice,
    bool neighbourTest) {
  int ih, ic, interIce, water, clath, interClath, total;
  assignIceTypePlus(yCloud, nList, isSlice, neighbourTest, ih, ic, interIce,
                    water, clath, interClath, total);
}

void chill::clearIceTypesExcept(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &types) {
  if (types.empty()) {
    return;
  }
  for (int i = 0; i < yCloud.nop; i++) {
    const int t = yCloud.pts[static_cast<std::size_t>(i)].type;
    if (std::find(types.begin(), types.end(), t) == types.end()) {
      yCloud.pts[static_cast<std::size_t>(i)].iceType =
          molSys::atom_state_type::unclassified;
    }
  }
}

void
chill::getIceTypePlus(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                      const std::vector<std::vector<int>> &nList, std::string path,
                      int firstFrame, bool isSlice, std::string outputFileName) {
  chill::getIceTypePlus(yCloud, nList, std::move(path), firstFrame, isSlice,
                        std::move(outputFileName), true);
}

void
chill::getIceTypePlus(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                      const std::vector<std::vector<int>> &nList, std::string path,
                      int firstFrame, bool isSlice,
                      std::string outputFileName, bool neighbourTest) {
  int ih, ic, interIce, water, clath, interClath, total;
  assignIceTypePlus(yCloud, nList, isSlice, neighbourTest, ih, ic, interIce,
                    water, clath, interClath, total);

  sout::makePath(path);
  std::string outputDirName = path + "bop";
  sout::makePath(outputDirName);

  std::ofstream outputFile;
  outputFile.open(path + "bop/" + outputFileName, std::ios_base::app);
  if (yCloud.currentFrame == firstFrame) {
    outputFile << "Frame Ic Ih Interfacial Clath InterClath Water Total\n";
  }
  outputFile << yCloud.currentFrame << " " << ic << " " << ih << " "
             << interIce << " " << clath << " " << interClath << " " << water
             << " " << total << "\n";
  outputFile.close();
}

// TODO: Add code for slices!
/**
 * @details Function for getting the averaged @f$q_6@f$ parameter.
 *  @param[in,out] yCloud The output molSys::PointCloud
 *  @param[in] nList Row-ordered neighbour list by atom ID
 *  @param[in] isSlice This decides whether there is a slice or not
 *  @return a double vector of the averaged @f$q_6@f$ values.
 */
std::vector<double>
chill::getq6(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
             const std::vector<std::vector<int>> &nList, bool isSlice) {
  //
  int l = 6;
  std::array<double, 3> delta;
  std::array<double, 2> angles;
  chill::QlmAtom QlmTotal; // Qlm for each iatom
  std::vector<std::complex<double>>
      yl; // temp q_lm for each pair of iatom and jatom
  std::complex<double> dot_product = {0, 0};
  std::complex<double> qI = {0, 0};
  std::complex<double> qJ = {0, 0};
  std::complex<double> Inorm = {0, 0};
  std::complex<double> Jnorm = {0, 0};
  std::complex<double> complexDenominator = {0, 0};
  std::complex<double> complexCij = {0, 0};
  double cij_real;
  std::vector<double> resultQ;
  double q_value = 0.0;        // Averaged q value per neighbour pair

  QlmTotal.ptq.resize(yCloud.nop);
  resultQ.resize(yCloud.nop);

  // Loop through every index in yCloud
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // if(yCloud.pts[iatom].type!=typeO){continue;}

    const auto nearest = nearestNeighbourIndices(yCloud, nList, iatom, 4);
    bool first = true;
    for (int jatomIndex : nearest) {
      delta = gen::relDist(yCloud, iatom, jatomIndex);
      double r = std::sqrt(std::pow(delta[0], 2.0) + std::pow(delta[1], 2.0) +
                           std::pow(delta[2], 2.0));
      if (r == 0.0) {
        continue;
      }
      angles[1] = acos(delta[2] / r);
      angles[0] = atan2(delta[0], delta[1]);

      if (first) {
        sph::spheriHarmoInto(l, angles, QlmTotal.ptq[iatom].ylm);
        first = false;
        continue;
      }
      sph::spheriHarmoInto(l, angles, yl);
      for (int m = 0; m < 2 * l + 1; m++) {
        QlmTotal.ptq[iatom].ylm[m] += yl[m];
      }
    }

    if (!nearest.empty()) {
      QlmTotal.ptq[iatom].ylm = gen::avgVector(
          QlmTotal.ptq[iatom].ylm, l, static_cast<int>(nearest.size()));
    }
  } // End of looping over all iatom

  // ------------------------------------------------
  // Now that you have all qlm for the particles,
  // find c_ij
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // if(yCloud.pts[iatom].type!=typeO){continue;}
    // if this is a slice and the particle is not in the slice
    // then skip TODO: UNCOMMENT
    // if(isSlice){
    // 	if(yCloud.pts[iatom].inSlice==false){continue;}
    // }

    const auto nearest = nearestNeighbourIndices(yCloud, nList, iatom, 4);
    q_value = 0.0;
    for (int jatomIndex : nearest) {
      dot_product = {0, 0};
      Inorm = {0, 0};
      Jnorm = {0, 0};

      for (int m = 0; m < 2 * l + 1; m++) {
        qI = QlmTotal.ptq[iatom].ylm[m];
        qJ = QlmTotal.ptq[jatomIndex].ylm[m];
        dot_product = dot_product + (qI * std::conj(qJ));
        Inorm = Inorm + (qI * std::conj(qI));
        Jnorm = Jnorm + (qJ * std::conj(qJ));
      }
      complexDenominator = std::sqrt(Inorm * Jnorm);
      complexCij = dot_product / complexDenominator;
      cij_real = complexCij.real();

      q_value += cij_real;
    }

    if (!nearest.empty()) {
      q_value /= static_cast<double>(nearest.size());
    }
    resultQ[iatom] = q_value;
  }

  // ------------------------------------------------

  return resultQ;
}

/**
 * @details Reclassifies atoms which may have been mis-classified
 *  as water using the averaged @f$q_6@f$ and @f$q_3@f$ parameters.
 *  This function can be called after both averaged @f$q_6@f$ and bond order
 *  correlation function @f$c_{ij}@f$ have been [calculated as described
 *  here](https://pubs.rsc.org/en/content/articlehtml/2011/cp/c1cp22167a).
 *
 *  @param[in,out] yCloud The output molSys::PointCloud
 *  @param[in] q6 Vector containing the previously calculated averaged @f$q_6@f$
 *   values (using chill::getq6)
 */
void chill::reclassifyWater(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    std::vector<double> &q6) {
  // If averaged q6 > 0.5, then consider it to be ice
  // If averaged q3 < -0.75 then it is ih or ic. If q3 < -0.85 then it is cubic,
  // otherwise it is hexagonal
  double avgQ3 = 0.0;
  int nnumNeighbours;

  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // Check if it has been classified as water
    if (yCloud.pts[iatom].iceType == molSys::atom_state_type::water) {
      if (q6[iatom] > 0.5) {
        avgQ3 = 0.0; // init to zero
        // Loop through all c_ij
        nnumNeighbours = yCloud.pts[iatom].c_ij.size();
        for (int j = 0; j < nnumNeighbours; j++) {
          avgQ3 += yCloud.pts[iatom].c_ij[j].c_value;
        }
        avgQ3 /= static_cast<double>(nnumNeighbours);

        // If averaged q3 < -0.75, then reclassify
        if (avgQ3 <= -0.75) {
          if (avgQ3 < -0.85) {
            yCloud.pts[iatom].iceType = molSys::atom_state_type::reCubic;
          } // molSys::cubic
          else {
            yCloud.pts[iatom].iceType = molSys::atom_state_type::reHex;
          } // molSys::hexagonal
        }   // end of reclassification
      }     // check for solid atom!
    }       // end of check for water
  }         // End loop through every iatom

  return;
}

/**
 * @details Prints out the molSys::atom_state_type per-particle ice type, for a
 *  particular frame, to a file.
 *  @param[in] yCloud The input molSys::PointCloud for the current frame
 *  @param[in] path Path to the output directory to which ice types are written
 *   out to
 *  @param[in] firstFrame First frame to be analyzed
 *  @param[in] isSlice Determines whether there is a slice or not
 *  @param[in] outputFileName File name of the output file, to which the
 *   per-particle ice types will be written out. The default file name is
 *   "superChill.txt"
 */
int chill::printIceType(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, std::string path,
    int firstFrame, bool isSlice, std::string outputFileName) {
  int ih, ic, interIce, water, unknown, clath, interClath,
      total; // No. of particles of each type
  ih = ic = water = unknown = interIce = total = 0;
  clath = interClath = 0;

  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // if(yCloud.pts[iatom].type != typeO){continue;}
    // check for slice
    if (isSlice) {
      if (yCloud.pts[iatom].inSlice == false) {
        continue;
      }
    }
    total++;
    if (yCloud.pts[iatom].iceType == molSys::atom_state_type::cubic) {
      ic++;
    } else if (yCloud.pts[iatom].iceType == molSys::atom_state_type::hexagonal) {
      ih++;
    } else if (yCloud.pts[iatom].iceType == molSys::atom_state_type::water) {
      water++;
    } else if (yCloud.pts[iatom].iceType == molSys::atom_state_type::interfacial) {
      interIce++;
    } else if (yCloud.pts[iatom].iceType == molSys::atom_state_type::clathrate) {
      clath++;
    } else if (yCloud.pts[iatom].iceType == molSys::atom_state_type::interClathrate) {
      interClath++;
    } else {
      unknown++;
    }
  }

  // --------------------
  // Create the directories if needed
  sout::makePath(path);
  std::string outputDirName = path + "bop";
  sout::makePath(outputDirName);
  // --------------------
  // Print to file
  std::ofstream outputFile;
  outputFile.open(path + "bop/" + outputFileName, std::ios_base::app);
  // --------------------
  // Write out the comment line
  if (yCloud.currentFrame == firstFrame) {
    outputFile << "Frame Ic Ih Interfacial Clath InterClath Water Total\n";
  }
  // --------------------
  outputFile << yCloud.currentFrame << " " << ic << " " << ih << " "
             << interIce << " " << clath << " " << interClath << " " << water
             << " " << total << "\n";
  outputFile.close();

  return 0;
}

/**
 *  @details Function that checks if the particle with the given atom index
 *   is interfacial or not.
 *  @param[in] yCloud The input molSys::PointCloud
 *  @param[in] nList Row-ordered neighbour list by atom ID
 *  @param[in] iatom The vector index of the current particle
 *  @param[in] num_staggrd The number of staggered bonds that the current
 *   particle participates in
 *  @param[in] num_eclipsd The number of eclipsed bonds that the current
 *   particle participates in
 *  @return a bool; true if the particle is interfacial and otherwise false
 */
bool chill::isInterfacial(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, int iatom, int num_staggrd,
    int num_eclipsd, bool chillPlus) {
  const std::vector<int> nearest =
      nearestNeighbourIndices(yCloud, nList, iatom, 4);
  int neighStaggered =
      0; // number of staggered bonds in the neighbours of iatom

  // INTERFACIAL
  // Condition 1 : only two staggered bonds and at least
  // one neighbor with more than two staggered bonds
  if (num_staggrd == 2) {
    // Loop over the four nearest neighbours, the same star as c_ij
    for (int jatomIndex : nearest) {
      neighStaggered = chill::numStaggered(yCloud, nList, jatomIndex);
      if (neighStaggered > 2) {
        return true;
      }
    } // End loop over nearest neighbours
  }   // end condition 1
  // Condition 2 : three staggered bonds, no eclipsed bond,
  // and at least one neighbor with two staggered bonds
  if (num_staggrd == 3 && num_eclipsd == 0) {
    for (int jatomIndex : nearest) {
      neighStaggered = chill::numStaggered(yCloud, nList, jatomIndex);
      if (chillPlus) {
        if (neighStaggered > 1) {
          return true;
        }
      } else if (neighStaggered == 2) {
        return true;
      }
    }
  }
  // not interfacial
  return false;
}

/**
 *  @details Calculates the number of staggered bonds of an atom
 *   with the given index.
 *  @param[in] yCloud The input molSys::PointCloud
 *  @param[in] nList Row-ordered neighbour list by atom ID
 *  @param[in] jatom The vector index of the current particle
 *  @return an int value, holding the number of staggered bonds of the given
 *   particle
 */
int chill::numStaggered(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, int jatom) {
  int num_staggrd = 0;        // Number of staggered bonds
  molSys::bond_type bondType; // Bond type
  int num_bonds;              // No. of bonds of the jatom
  // Bound on the bonds recorded for jatom, not on its neighbour-list row
  int nnumNeighbours = yCloud.pts[jatom].c_ij.size();

  // Loop over all bonds
  for (int i = 0; i < nnumNeighbours; i++) {
    bondType = yCloud.pts[jatom].c_ij[i].classifier;
    // If the bond is staggered increment the number of staggered bonds
    if (bondType == molSys::bond_type::staggered) {
      num_staggrd++;
    }
  } // end of loop over c_ij

  return num_staggrd;
}

namespace {

int resolveId(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
              const std::vector<int> &idToIndex, int maxId, bool dense, int id) {
  if (dense) {
    return (id >= 0 && id <= maxId) ? idToIndex[static_cast<std::size_t>(id)] : -1;
  }
  const auto it = yCloud.idIndexMap.find(id);
  return it == yCloud.idIndexMap.end() ? -1 : it->second;
}

} // namespace

void chill::bondOrder(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                      const std::vector<std::vector<int>> &nList,
                      const BondClassifier &rule, BondOrder &order,
                      const std::vector<char> *which, bool fresh) {
  seams::prepareCompute();
  [[maybe_unused]] const seams::TeamWidth team(yCloud.nop);
  const int nop = yCloud.nop;
  const int keepMax = rule.coordinationNumber;
  order.keepMax = keepMax;
  const int qn = std::max(nop, 0) * kRealL3;
  const bool sized = static_cast<int>(order.q.size()) == qn &&
                     static_cast<int>(order.norm.size()) == nop &&
                     static_cast<int>(order.nkeep.size()) == nop &&
                     (keepMax > 0
                          ? static_cast<int>(order.chosen.size()) == nop * keepMax
                          : static_cast<int>(order.rows.size()) == std::max(nop, 0));
  if (!sized) {
    order.q.assign(static_cast<std::size_t>(qn), 0.0);
    order.norm.assign(static_cast<std::size_t>(std::max(nop, 0)), 0.0);
    order.nkeep.assign(static_cast<std::size_t>(std::max(nop, 0)), 0);
    order.chosen.clear();
    order.rows.clear();
    if (nop > 0) {
      if (keepMax > 0) {
        order.chosen.assign(static_cast<std::size_t>(nop) * static_cast<std::size_t>(keepMax), -1);
      } else {
        order.rows.assign(static_cast<std::size_t>(nop), {});
      }
    }
  } else if (fresh && nop > 0) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
    for (int i = 0; i < qn; i++) {
      order.q[static_cast<std::size_t>(i)] = 0.0;
    }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
    for (int i = 0; i < nop; i++) {
      order.norm[static_cast<std::size_t>(i)] = 0.0;
      order.nkeep[static_cast<std::size_t>(i)] = 0;
      if (keepMax > 0) {
        for (int s = 0; s < keepMax; s++) {
          order.chosen[static_cast<std::size_t>(i) * static_cast<std::size_t>(keepMax) +
                       static_cast<std::size_t>(s)] = -1;
        }
      } else {
        order.rows[static_cast<std::size_t>(i)].clear();
      }
    }
  }
  if (nop <= 0) {
    return;
  }
  int maxId = 0;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) reduction(max : maxId) if (omp_get_max_threads() > 1)
#endif
  for (int i = 0; i < nop; i++) {
    maxId = std::max(maxId, yCloud.pts[static_cast<std::size_t>(i)].atomID);
  }
  // Atom IDs are frame-wide. A rank's cloud is only its block, so the old
  // nop*4 test left every eight-rank lookup on the hash map.
  const bool dense = maxId > 0 && maxId < 4000000;
  std::vector<int> idToIndex;
  if (dense) {
    idToIndex.assign(static_cast<std::size_t>(maxId) + 1, -1);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
    for (int i = 0; i < nop; i++) {
      const int id = yCloud.pts[static_cast<std::size_t>(i)].atomID;
      if (id >= 0) {
        idToIndex[static_cast<std::size_t>(id)] = i;
      }
    }
  }
  struct Xyz {
    double x, y, z;
  };
  std::vector<Xyz> xyz(static_cast<std::size_t>(nop));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
  for (int i = 0; i < nop; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(i)] = {p.x, p.y, p.z};
  }
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  const bool useAvx = hostHasAvx2Fma();
  const bool fixedKeep = keepMax > 0;
#ifdef SEAMS_HAS_OPENMP
  const double orderT0 = omp_get_wtime();
  const int orderThreads = omp_get_max_threads();
#pragma omp parallel for schedule(static) if (nop >= 16384 && omp_get_max_threads() > 1)
#endif
  for (int iatom = 0; iatom < nop; iatom++) {
    if (which != nullptr &&
        (static_cast<int>(which->size()) <= iatom || (*which)[static_cast<std::size_t>(iatom)] == 0)) {
      continue;
    }
    if (static_cast<std::size_t>(iatom) >= nList.size()) {
      continue;
    }
    struct Near {
      double d2;
      int jatom;
      double x, y, z;
      bool operator<(const Near &o) const {
        return d2 < o.d2 || (d2 == o.d2 && jatom < o.jatom);
      }
    };
    Near stackNear[32];
    std::vector<Near> heapNear;
    Near *candidates = stackNear;
    int nCand = 0;
    const Xyz pi = xyz[static_cast<std::size_t>(iatom)];
    const auto &row = nList[static_cast<std::size_t>(iatom)];
    for (std::size_t j = 1; j < row.size(); j++) {
      const int jatom = resolveId(yCloud, idToIndex, maxId, dense, row[j]);
      if (jatom < 0) {
        continue;
      }
      const Xyz pj = xyz[static_cast<std::size_t>(jatom)];
      const auto dr = gen::fracDelta(frame, pi.x, pi.y, pi.z, pj.x, pj.y, pj.z);
      if (nCand == 32 && candidates == stackNear) {
        heapNear.assign(stackNear, stackNear + 32);
        candidates = heapNear.data();
      }
      const Near item{dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2], jatom, dr[0], dr[1],
                      dr[2]};
      if (candidates == stackNear) {
        candidates[nCand] = item;
      } else {
        heapNear.push_back(item);
        candidates = heapNear.data();
      }
      ++nCand;
    }
    const int keep = fixedKeep ? std::min(keepMax, nCand) : nCand;
    if (keep == 0) {
      continue;
    }
    std::partial_sort(candidates, candidates + keep, candidates + nCand);
    order.nkeep[static_cast<std::size_t>(iatom)] = keep;
    double acc[kRealL3] = {};
    for (int j = 0; j < keep;) {
      double xs[4], ys[4], zs[4];
      int batch = 0;
      for (; batch < 4 && j < keep; ++j, ++batch) {
        if (fixedKeep) {
          order.chosen[static_cast<std::size_t>(iatom) * static_cast<std::size_t>(keepMax) +
                       static_cast<std::size_t>(j)] = candidates[j].jatom;
        } else {
          order.rows[static_cast<std::size_t>(iatom)].push_back(candidates[j].jatom);
        }
        xs[batch] = candidates[j].x;
        ys[batch] = candidates[j].y;
        zs[batch] = candidates[j].z;
      }
      if (useAvx) {
        realYlm3Accumulate(xs, ys, zs, batch, acc);
      } else {
        for (int item = 0; item < batch; ++item) {
          addRealYlm3Scalar(xs[item], ys[item], zs[item], normQ3(), acc);
        }
      }
    }
    double *dst = order.q.data() + static_cast<std::size_t>(iatom) * kRealL3;
    for (int channel = 0; channel < kRealL3; ++channel) {
      dst[channel] = acc[channel];
    }
    order.norm[static_cast<std::size_t>(iatom)] =
        useAvx ? realL3Dot(acc, acc) : realL3DotScalar(acc, acc);
  }
#ifdef SEAMS_HAS_OPENMP
  if (const char *split = std::getenv("SEAMS_SPLIT")) {
    if (split[0] == '1') {
      std::fprintf(stderr, "order threads %d nop %d ms %.3f\n", orderThreads, nop,
                   (omp_get_wtime() - orderT0) * 1.0e3);
    }
  }
#endif
}

void chill::correlFromOrder(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                            const BondOrder &order, const BondClassifier &rule, bool isSlice,
                            const std::vector<char> *which, bool closed) {
  const int nop = yCloud.nop;
  const int keepMax = order.keepMax;
  const bool useAvx = hostHasAvx2Fma();
  const auto neighbourAt = [&](int atom, int slot) {
    if (keepMax > 0) {
      return order.chosen[static_cast<std::size_t>(atom) * static_cast<std::size_t>(keepMax) +
                          static_cast<std::size_t>(slot)];
    }
    return order.rows[static_cast<std::size_t>(atom)][static_cast<std::size_t>(slot)];
  };
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nop >= 16384 && omp_get_max_threads() > 1)
#endif
  for (int iatom = 0; iatom < nop; iatom++) {
    if (which != nullptr &&
        (static_cast<int>(which->size()) <= iatom || (*which)[static_cast<std::size_t>(iatom)] == 0)) {
      continue;
    }
    if (isSlice && !yCloud.pts[static_cast<std::size_t>(iatom)].inSlice) {
      continue;
    }
    const int keep = order.nkeep[static_cast<std::size_t>(iatom)];
    auto &bonds = yCloud.pts[static_cast<std::size_t>(iatom)].c_ij;
    bonds.clear();
    if (keep <= 0) {
      continue;
    }
    const double *qI = order.q.data() + static_cast<std::size_t>(iatom) * kRealL3;
    const double normI = order.norm[static_cast<std::size_t>(iatom)];
    bonds.reserve(static_cast<std::size_t>(keep));
    for (int slot = 0; slot < keep; ++slot) {
      const int jatom = neighbourAt(iatom, slot);
      if (jatom < 0) {
        continue;
      }
      const double *qJ = order.q.data() + static_cast<std::size_t>(jatom) * kRealL3;
      const double dot = useAvx ? realL3Dot(qI, qJ) : realL3DotScalar(qI, qJ);
      const double cij =
          correlationOf(dot, normI, order.norm[static_cast<std::size_t>(jatom)]);
      molSys::Result temp;
      temp.c_value = cij;
      temp.classifier = chill::classifyBond(cij, rule, closed);
      bonds.push_back(temp);
    }
  }
}

#ifdef SEAMS_HAS_MPI
chill::AnalysisPhases gChillPhases;

namespace {

const std::vector<std::vector<int>> &shellRows(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud, double cutoff,
    const std::vector<char> &owned) {
  bool single = true;
  for (const char flag : owned) {
    if (!flag) {
      single = false;
      break;
    }
  }
  static std::vector<std::vector<int>> window;
  if (nneigh::blockWindowRows(cloud, cutoff, single ? nullptr : &owned, window)) {
    return window;
  }
  if (single || cloud.nop < 40000) {
    if (!nneigh::refillNeighListO(cutoff, cloud, 1, window)) {
      window = nneigh::neighListO(cutoff, cloud, 1);
    }
  } else {
    window = nneigh::reusedShellList(cutoff, cloud, &owned);
  }
  return window;
}

struct ChillDuring {
  double cutoff = 3.5;
  chill::BondClassifier rule = chill::chillPlusRule();
  chill::BondOrder order;
  std::vector<char> interior;
  std::vector<std::vector<int>> rows;
  bool ran = false;
};

ChillDuring gDuring;

void duringChillHalo(void *) {
  // Ghost coordinates were copied from the frame this rank already holds,
  // so the owned bond-order kernel can run while the index exchange is in
  // flight. Rows stay in gDuring; copying them allocated one vector per atom.
  const auto t0 = std::chrono::steady_clock::now();
  auto &work = seams::domain::gatheredWork();
  bool single = true;
  for (const char flag : work.share.owned) {
    if (!flag) {
      single = false;
      break;
    }
  }
  if (!nneigh::blockWindowRows(work.cloud, gDuring.cutoff, single ? nullptr : &work.share.owned,
                               gDuring.rows)) {
    gDuring.rows = shellRows(work.cloud, gDuring.cutoff, work.share.owned);
  }
  gChillPhases.rebuildMs =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  chill::bondOrder(work.cloud, gDuring.rows, gDuring.rule, gDuring.order, &work.share.owned,
                   true);
  gDuring.ran = true;
  gChillPhases.computeMs =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() -
      gChillPhases.rebuildMs;
}

void writeStaggered(molSys::PointCloud<molSys::Point<double>, double> &cloud, int atom,
                    int count) {
  auto &bonds = cloud.pts[static_cast<std::size_t>(atom)].c_ij;
  bonds.clear();
  molSys::Result bond;
  bond.c_value = -1.0;
  bond.classifier = molSys::bond_type::staggered;
  bonds.assign(static_cast<std::size_t>(std::max(count, 0)), bond);
}

} // namespace

void chill::distributedChillPlus(molSys::PointCloud<molSys::Point<double>, double> &full,
                                 double cutoff, BondClassifier rule, bool neighbourTest,
                                 MPI_Comm comm, bool closed) {
  int rank = 0;
  int nRanks = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nRanks);
  if (nRanks <= 1) {
    const auto t0 = std::chrono::steady_clock::now();
    thread_local std::vector<std::vector<int>> held;
    if (!nneigh::refillNeighListO(cutoff, full, 1, held)) {
      held = nneigh::neighListO(cutoff, full, 1);
    }
    const auto &rows = held;
    const auto t1 = std::chrono::steady_clock::now();
    classifyBonds(full, rows, rule, false, closed);
    getIceTypePlusNoPrint(full, rows, false, neighbourTest);
    const auto t2 = std::chrono::steady_clock::now();
    gChillPhases.rebuildMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    gChillPhases.computeMs = std::chrono::duration<double, std::milli>(t2 - t1).count();
    gChillPhases.exchangeMs = 0;
    gChillPhases.reduceMs = 0;
    return;
  }
  gChillPhases = {};
  gDuring.cutoff = cutoff;
  gDuring.rule = rule;
  gDuring.ran = false;
  seams::domain::setHaloHook(duringChillHalo, nullptr);
#ifdef SEAMS_HAS_OPENMP
  const double gatherT0 = omp_get_wtime();
#endif
  auto &work = seams::domain::gatherWork(full, cutoff, comm, seams::domain::kWorkChill);
  seams::domain::setHaloHook(nullptr, nullptr);
#ifdef SEAMS_HAS_OPENMP
  const double gatherT1 = omp_get_wtime();
  if (const char *split = std::getenv("SEAMS_SPLIT")) {
    if (split[0] == '1') {
      std::fprintf(stderr, "gather ms %.3f local %zu\n", (gatherT1 - gatherT0) * 1.0e3,
                   work.share.local.size());
    }
  }
#endif
  BondOrder &order = gDuring.order;
  const std::vector<std::vector<int>> &rows = gDuring.ran
                                                   ? gDuring.rows
                                                   : shellRows(work.cloud, cutoff, work.share.owned);
  if (!gDuring.ran) {
    bondOrder(work.cloud, rows, rule, order, &work.share.owned, true);
  }
  const int nop = work.cloud.nop;
  constexpr int stride = kRealL3 + 1;
  std::vector<double> packed(static_cast<std::size_t>(nop) * stride, 0.0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
  for (int i = 0; i < nop; i++) {
    if (!work.share.owned[static_cast<std::size_t>(i)]) {
      continue;
    }
    double *dst = packed.data() + static_cast<std::size_t>(i) * stride;
    const double *src = order.q.data() + static_cast<std::size_t>(i) * kRealL3;
    std::copy_n(src, kRealL3, dst);
    dst[kRealL3] = order.norm[static_cast<std::size_t>(i)];
  }
  seams::domain::postBoundaryValues(packed.data(), stride, comm);
  std::vector<char> interior(static_cast<std::size_t>(nop), 0);
  std::vector<char> boundary(static_cast<std::size_t>(nop), 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
  for (int i = 0; i < nop; i++) {
    if (!work.share.owned[static_cast<std::size_t>(i)]) {
      continue;
    }
    const int keep = order.nkeep[static_cast<std::size_t>(i)];
    bool inside = true;
    for (int slot = 0; slot < keep; ++slot) {
      const int j = rule.coordinationNumber > 0
                        ? order.chosen[static_cast<std::size_t>(i) *
                                           static_cast<std::size_t>(rule.coordinationNumber) +
                                       static_cast<std::size_t>(slot)]
                        : order.rows[static_cast<std::size_t>(i)][static_cast<std::size_t>(slot)];
      if (j < 0 || !work.share.owned[static_cast<std::size_t>(j)]) {
        inside = false;
        break;
      }
    }
    (inside ? interior : boundary)[static_cast<std::size_t>(i)] = 1;
  }
#ifdef SEAMS_HAS_OPENMP
  const double corrT0 = omp_get_wtime();
#endif
  const auto corr0 = std::chrono::steady_clock::now();
  correlFromOrder(work.cloud, order, rule, false, &interior, closed);
  gChillPhases.computeMs +=
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - corr0).count();
  seams::domain::waitBoundaryValues();
  gChillPhases.exchangeMs = seams::domain::takeValueCommMs();
  for (int i = 0; i < nop; i++) {
    if (work.share.owned[static_cast<std::size_t>(i)]) {
      continue;
    }
    const double *src = packed.data() + static_cast<std::size_t>(i) * stride;
    std::copy_n(src, kRealL3, order.q.data() + static_cast<std::size_t>(i) * kRealL3);
    order.norm[static_cast<std::size_t>(i)] = src[kRealL3];
  }
  const auto corr1 = std::chrono::steady_clock::now();
  correlFromOrder(work.cloud, order, rule, false, &boundary, closed);
  gChillPhases.computeMs +=
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - corr1).count();
  // isInterfacial reads the four nearest neighbours' staggered counts.
  // An owned neighbour already has its bonds. A ghost count is needed
  // only when a boundary atom, this frame, is in one of those two rows.
  // Every rank must agree or the exchange deadlocks.
  int localNeed = 0;
  if (neighbourTest) {
    for (int i = 0; i < nop; i++) {
      if (!work.share.owned[static_cast<std::size_t>(i)] ||
          !boundary[static_cast<std::size_t>(i)]) {
        continue;
      }
      int stag = 0;
      int ecl = 0;
      for (const auto &bond : work.cloud.pts[static_cast<std::size_t>(i)].c_ij) {
        if (bond.classifier == molSys::bond_type::staggered) {
          ++stag;
        } else if (bond.classifier == molSys::bond_type::eclipsed) {
          ++ecl;
        }
      }
      if (stag == 2 || (stag == 3 && ecl == 0)) {
        localNeed = 1;
        break;
      }
    }
  }
  int anyNeed = localNeed;
  MPI_Allreduce(&localNeed, &anyNeed, 1, MPI_INT, MPI_MAX, comm);
  if (anyNeed) {
    std::vector<double> stag(static_cast<std::size_t>(nop), 0.0);
    for (int i = 0; i < nop; i++) {
      if (!work.share.owned[static_cast<std::size_t>(i)]) {
        continue;
      }
      int n = 0;
      for (const auto &bond : work.cloud.pts[static_cast<std::size_t>(i)].c_ij) {
        if (bond.classifier == molSys::bond_type::staggered) {
          ++n;
        }
      }
      stag[static_cast<std::size_t>(i)] = static_cast<double>(n);
    }
    seams::domain::postBoundaryValues(stag.data(), 1, comm);
    seams::domain::waitBoundaryValues();
    gChillPhases.exchangeMs += seams::domain::takeValueCommMs();
    for (int i = 0; i < nop; i++) {
      if (work.share.owned[static_cast<std::size_t>(i)]) {
        continue;
      }
      writeStaggered(work.cloud, i, static_cast<int>(std::lround(stag[static_cast<std::size_t>(i)])));
    }
  }
#ifdef SEAMS_HAS_OPENMP
  const double typeT0 = omp_get_wtime();
#endif
  const auto type0 = std::chrono::steady_clock::now();
  getIceTypePlusNoPrint(work.cloud, rows, false, neighbourTest);
  gChillPhases.computeMs +=
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - type0).count();
#ifdef SEAMS_HAS_OPENMP
  if (const char *split = std::getenv("SEAMS_SPLIT")) {
    if (split[0] == '1') {
      const double now = omp_get_wtime();
      std::fprintf(stderr, "chill rest corr %.3f type %.3f\n", (typeT0 - corrT0) * 1.0e3,
                   (now - typeT0) * 1.0e3);
    }
  }
#endif
  if (!rows.empty()) {
    std::vector<double> effort(rows.size(), 0.0);
    for (std::size_t k = 0; k < rows.size() && k < work.share.owned.size(); k++) {
      if (work.share.owned[k]) {
        effort[k] = static_cast<double>(rows[k].size());
      }
    }
    seams::domain::captureLocalWork(seams::domain::kWorkChill, effort.data(),
                                    static_cast<int>(effort.size()), comm);
  }
  const auto red0 = std::chrono::steady_clock::now();
  for (std::size_t k = 0; k < work.share.local.size(); k++) {
    if (!work.share.owned[k]) {
      continue;
    }
    full.pts[static_cast<std::size_t>(work.share.local[k])].iceType =
        work.cloud.pts[k].iceType;
  }
  gChillPhases.reduceMs =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - red0).count();
}

chill::AnalysisPhases chill::lastChillPhases() { return gChillPhases; }
#endif
