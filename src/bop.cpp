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

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <ozaki.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
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
    bool isSlice) {
  constexpr int l = 3;
  constexpr int nComp = 2 * l + 1;
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

  // Real and imaginary parts of q_lm, m = -l..l. Dividing by the neighbour
  // count cancels in c_ij, so the sum is stored unnormalised.
  // Coordinates are copied first so the parallel loop does not read a Point
  // while another thread clears that Point's bond list.
  struct Xyz {
    double x, y, z;
  };
  std::vector<Xyz> xyz(static_cast<std::size_t>(nop));
  for (int i = 0; i < nop; i++) {
    auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(i)] = {p.x, p.y, p.z};
    p.c_ij.clear();
  }
  std::vector<double> qReal(static_cast<std::size_t>(nop) * nComp, 0.0);
  std::vector<double> qImag(static_cast<std::size_t>(nop) * nComp, 0.0);
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

#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
  for (int iatom = 0; iatom < nop; iatom++) {
    if (static_cast<std::size_t>(iatom) >= nList.size()) {
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
    const auto &row = nList[static_cast<std::size_t>(iatom)];
    for (std::size_t j = 1; j < row.size(); j++) {
      const int jatom = resolve(row[j]);
      if (jatom < 0) {
        continue;
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
    }
    const int keep =
        (keepMax > 0) ? std::min(keepMax, nCand) : nCand;
    if (keep == 0) {
      continue;
    }
    std::partial_sort(candidates, candidates + keep, candidates + nCand);
    nkeep[static_cast<std::size_t>(iatom)] = keep;
    std::array<double, nComp> accRe{};
    std::array<double, nComp> accIm{};
    for (int j = 0; j < keep; j++) {
      const int jatom = candidates[j].jatom;
      if (fixedKeep) {
        chosen[static_cast<std::size_t>(iatom) * static_cast<std::size_t>(keepMax) +
               static_cast<std::size_t>(j)] = jatom;
      } else if (serialShell) {
        flat.push_back(jatom);
      } else {
        ragged[static_cast<std::size_t>(iatom)].push_back(jatom);
      }
      const std::array<double, 3> dr = {candidates[j].x, candidates[j].y, candidates[j].z};
      const std::array<double, 2> angles = sph::radialCoord(dr);
      const AngularTerms terms(angles[1], angles[0], l);
      for (int m = 0; m <= l; m++) {
        const auto harm = harmonicPair(l, m, terms);
        accRe[l - m] += harm.first.real();
        accIm[l - m] += harm.first.imag();
        if (m == 0) {
          continue;
        }
        accRe[l + m] += harm.second.real();
        accIm[l + m] += harm.second.imag();
      }
    }
    double *re = qReal.data() + static_cast<std::size_t>(iatom) * nComp;
    double *im = qImag.data() + static_cast<std::size_t>(iatom) * nComp;
    for (int m = 0; m < nComp; m++) {
      re[m] = accRe[m];
      im[m] = accIm[m];
    }
  }
  if (serialShell) {
    for (int i = 0; i < nop; i++) {
      flatAt[static_cast<std::size_t>(i + 1)] =
          flatAt[static_cast<std::size_t>(i)] + nkeep[static_cast<std::size_t>(i)];
    }
  }

  // c_ij from the real part of the complex inner product. The neighbour
  // count cancels between the averaged harmonics and the norms.
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
  for (int iatom = 0; iatom < nop; iatom++) {
    if (isSlice && !yCloud.pts[static_cast<std::size_t>(iatom)].inSlice) {
      continue;
    }
    const int keep = nkeep[static_cast<std::size_t>(iatom)];
    if (keep <= 0) {
      continue;
    }
    const double *reI = qReal.data() + static_cast<std::size_t>(iatom) * nComp;
    const double *imI = qImag.data() + static_cast<std::size_t>(iatom) * nComp;
    double normI = 0.0;
    for (int m = 0; m < nComp; m++) {
      normI += reI[m] * reI[m] + imI[m] * imI[m];
    }
    auto &bonds = yCloud.pts[static_cast<std::size_t>(iatom)].c_ij;
    bonds.reserve(static_cast<std::size_t>(keep));
    for (int j = 0; j < keep; j++) {
      const int jatom =
          fixedKeep ? chosen[static_cast<std::size_t>(iatom) * static_cast<std::size_t>(keepMax) +
                             static_cast<std::size_t>(j)]
          : serialShell ? flat[static_cast<std::size_t>(flatAt[static_cast<std::size_t>(iatom)] + j)]
                        : ragged[static_cast<std::size_t>(iatom)][static_cast<std::size_t>(j)];
      if (jatom < 0) {
        continue;
      }
      const double *reJ = qReal.data() + static_cast<std::size_t>(jatom) * nComp;
      const double *imJ = qImag.data() + static_cast<std::size_t>(jatom) * nComp;
      double dot = 0.0;
      double normJ = 0.0;
      for (int m = 0; m < nComp; m++) {
        dot += reI[m] * reJ[m] + imI[m] * imJ[m];
        normJ += reJ[m] * reJ[m] + imJ[m] * imJ[m];
      }
      const double cij_real = dot / std::sqrt(normI * normJ);
      molSys::Result temp_cij;
      temp_cij.c_value = cij_real;
      temp_cij.classifier = chill::classifyBond(cij_real, rule);
      bonds.push_back(temp_cij);
    }
  }
}

void chill::classifyBonds(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, const chill::BondClassifier &rule,
    bool isSlice) {
  // Left-to-right fp64 is the path CHILL and CHILL+ run on. A different
  // product reduction or a non-fp64 dot product keeps the pairwise sum.
  if (seams::product::reduction() == seams::product::Reduction::Left &&
      seams::product::mode() == seams::product::Mode::Fp64) {
    classifyBondsLeftFp64(yCloud, nList, rule, isSlice);
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
      const std::complex<double> complexCij =
          dot_product / std::sqrt(Inorm * Jnorm);
      molSys::Result temp_cij;
      const double cij_real = complexCij.real();
      temp_cij.c_value = cij_real;
      temp_cij.classifier = chill::classifyBond(cij_real, rule);
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

/**
 *  @details Function for getting the bond order correlations @f$c_{ij}@f$
 * (alternatively
 *   @f$a_{ij}@f$ in certain texts) using the CHILL+ algorithm
 *  @param[in,out] yCloud The output molSys::PointCloud
 *  @param[in] nList Row-ordered neighbour list by atom ID
 *  @param[in] isSlice This decides whether there is a slice or not
 */
void
chill::getCorrelPlus(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     const std::vector<std::vector<int>> &nList, bool isSlice,
                     int coordinationNumber, bool closedBounds) {
  chill::BondClassifier rule = chill::chillPlusRule();
  rule.coordinationNumber = coordinationNumber;
  rule.closed = closedBounds;
  chill::classifyBonds(yCloud, nList, rule, isSlice);
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
