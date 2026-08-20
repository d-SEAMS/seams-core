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
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <unordered_map>
#include <utility>

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

  AngularTerms(double theta, double phi) {
    const double sinTheta = std::sin(theta);
    const double cosTheta = std::cos(theta);
    const std::complex<double> unitPhase(std::cos(phi), std::sin(phi));

    sinPow[0] = 1.0;
    cosPow[0] = 1.0;
    phase[0] = {1.0, 0.0};
    for (std::size_t k = 1; k < sinPow.size(); k++) {
      sinPow[k] = sinPow[k - 1] * sinTheta;
      cosPow[k] = cosPow[k - 1] * cosTheta;
      phase[k] = phase[k - 1] * unitPhase;
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
    switch (absM) {
    case 0:
      return 5.0 * c[3] - 3.0 * c[1];
    case 1:
      return s[1] * (5.0 * c[2] - 1.0);
    case 2:
      return s[2] * c[1];
    case 3:
      return s[3];
    default:
      return 0.0;
    }
  }

  if (orderL == 4) {
    switch (absM) {
    case 0:
      return 35.0 * c[4] - 30.0 * c[2] + 3.0;
    case 1:
      return s[1] * (7.0 * c[3] - 3.0 * c[1]);
    case 2:
      return s[2] * (7.0 * c[2] - 1.0);
    case 3:
      return s[3] * c[1];
    case 4:
      return s[4];
    default:
      return 0.0;
    }
  }

  if (orderL == 6) {
    switch (absM) {
    case 0:
      return 231.0 * c[6] - 315.0 * c[4] + 105.0 * c[2] - 5.0;
    case 1:
      return s[1] * (33.0 * c[5] - 30.0 * c[3] + 5.0 * c[1]);
    case 2:
      return s[2] * (33.0 * c[4] - 18.0 * c[2] + 1.0);
    case 3:
      return s[3] * (11.0 * c[3] - 3.0 * c[1]);
    case 4:
      return s[4] * (11.0 * c[2] - 1.0);
    case 5:
      return s[5] * c[1];
    case 6:
      return s[6];
    default:
      return 0.0;
    }
  }

  if (orderL == 8) {
    switch (absM) {
    case 0:
      return (35.0 + c[2] * (-1260.0 + c[2] * (6930.0 + c[2] * (-12012.0 + c[2] * 6435.0)))) / 128.0;
    case 1:
      return s[1] * c[1] * (-315.0 + c[2] * (3465.0 + c[2] * (-9009.0 + c[2] * 6435.0))) / 16.0;
    case 2:
      return s[2] * (-315.0 + c[2] * (10395.0 + c[2] * (-45045.0 + c[2] * 45045.0))) / 16.0;
    case 3:
      return s[3] * c[1] * (10395.0 + c[2] * (-90090.0 + c[2] * 135135.0)) / 8.0;
    case 4:
      return s[4] * (10395.0 + c[2] * (-270270.0 + c[2] * 675675.0)) / 8.0;
    case 5:
      return s[5] * c[1] * (-135135.0 + c[2] * 675675.0) / 2.0;
    case 6:
      return s[6] * (-135135.0 + c[2] * 2027025.0) / 2.0;
    case 7:
      return s[7] * c[1] * 2027025.0;
    case 8:
      return s[8] * 2027025.0;
    default:
      return 0.0;
    }
  }

  return 0.0;
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
  for (size_t j = 1; j < nList[iatomIndex].size(); j++) {
    const auto it = yCloud.idIndexMap.find(nList[iatomIndex][j]);
    if (it == yCloud.idIndexMap.end()) {
      continue;
    }
    candidates.emplace_back(gen::periodicDistSq(yCloud, iatomIndex, it->second),
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

} // namespace

/**
 * @details Function for calculating spherical harmonics. Dispatches to
 *  the closed forms for l=3, l=4 and l=6.
 *
 *  @param[in] orderL The int value of l (must be 3, 4 or 6)
 *  @param[in] radialCoord Array containing the polar and azimuth angles
 *  @return a complex vector of length 2l+1
 */
std::vector<std::complex<double>>
sph::spheriHarmo(int orderL, std::array<double, 2> radialCoord) {
  if (orderL == 3) {
    return sph::lookupTableQ3Vec(radialCoord);
  } else if (orderL == 4) {
    return sph::lookupTableQ4Vec(radialCoord);
  } else if (orderL == 6) {
    return sph::lookupTableQ6Vec(radialCoord);
  } else if (orderL == 8) {
    return sph::lookupTableQ8Vec(radialCoord);
  }
  // Fallback: return zeros for unsupported l values
  return std::vector<std::complex<double>>(2 * orderL + 1, {0.0, 0.0});
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
void chill::classifyBonds(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, const chill::BondClassifier &rule,
    bool isSlice) {
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
    for (int j = 0; j < nnumNeighbours; j++) {
      const int jatomIndex = nearest[j];

      // Get the relative distance now that the index values are known
      const std::array<double, 3> delta =
          gen::relDist(yCloud, iatom, jatomIndex);
      // radialCoord carries the r > 0 guard, so coincident atoms yield a
      // defined angle instead of acos(0/0)
      const std::array<double, 2> angles = sph::radialCoord(delta);

      // Now add over all nearest neighbours
      if (j == 0) {
        QlmTotal.ptq[iatom].ylm = sph::spheriHarmo(l, angles);
        continue;
      }
      yl = sph::spheriHarmo(l, angles);
      for (int m = 0; m < 2 * l + 1; m++) {
        QlmTotal.ptq[iatom].ylm[m] += yl[m];
      }
    } // End of loop over nearest neighbours

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
    for (const int jatomIndex : nearestBonds) {
      std::complex<double> dot_product = {0, 0};
      std::complex<double> Inorm = {0, 0};
      std::complex<double> Jnorm = {0, 0};
      for (int m = 0; m < 2 * l + 1; m++) {
        const std::complex<double> qI = QlmTotal.ptq[iatom].ylm[m];
        const std::complex<double> qJ = QlmTotal.ptq[jatomIndex].ylm[m];
        dot_product += qI * std::conj(qJ); // unnormalized
        Inorm += qI * std::conj(qI);
        Jnorm += qJ * std::conj(qJ);
      } // end loop over m components
      const std::complex<double> complexCij =
          dot_product / std::sqrt(Inorm * Jnorm);
      molSys::Result temp_cij;
      const double cij_real = complexCij.real();
      temp_cij.c_value = cij_real;
      if (cij_real <= rule.staggeredMax) {
        temp_cij.classifier = molSys::bond_type::staggered;
      } else if (cij_real >= rule.eclipsedMin &&
                 cij_real <= rule.eclipsedMax) {
        temp_cij.classifier = molSys::bond_type::eclipsed;
      } else {
        temp_cij.classifier = molSys::bond_type::out_of_range;
      }
      yCloud.pts[iatom].c_ij.push_back(temp_cij);
    } // end loop over nearest neighbours
  }
}

chill::BondClassifier chill::chillRule() { return {-0.8, -0.2, -0.05, 4}; }

chill::BondClassifier chill::chillPlusRule() { return {-0.8, -0.35, 0.25, 4}; }

namespace {
std::unordered_map<std::string, chill::BondClassifier> &bondRuleRegistry() {
  static std::unordered_map<std::string, chill::BondClassifier> registry = {
      {"CHILL", chill::chillRule()}, {"CHILL+", chill::chillPlusRule()}};
  return registry;
}
} // namespace

chill::BondClassifier chill::bondClassifier(const std::string &name) {
  return bondRuleRegistry().at(name);
}

void chill::registerBondClassifier(const std::string &name,
                                   const chill::BondClassifier &rule) {
  bondRuleRegistry()[name] = rule;
}

std::vector<std::string> chill::bondClassifierNames() {
  std::vector<std::string> names;
  names.reserve(bondRuleRegistry().size());
  for (const auto &entry : bondRuleRegistry()) {
    names.push_back(entry.first);
  }
  std::sort(names.begin(), names.end());
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
                     int coordinationNumber) {
  chill::BondClassifier rule = chill::chillPlusRule();
  rule.coordinationNumber = coordinationNumber;
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
    const std::vector<std::vector<int>> &nList, bool isSlice, int &ih, int &ic,
    int &interIce, int &water, int &clath, int &interClath, int &total) {
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
    // Bound on the bonds that were recorded, not on the neighbour list
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
    if (nnumNeighbours == 4) {
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
      // Interfacial
      else if (chill::isInterfacial(yCloud, nList, iatom, num_staggrd, num_eclipsd,
                             true)) {
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
  int ih, ic, interIce, water, clath, interClath, total;
  assignIceTypePlus(yCloud, nList, isSlice, ih, ic, interIce, water, clath,
                    interClath, total);
}

void
chill::getIceTypePlus(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                      const std::vector<std::vector<int>> &nList, std::string path,
                      int firstFrame, bool isSlice,
                      std::string outputFileName) {
  int ih, ic, interIce, water, clath, interClath, total;
  assignIceTypePlus(yCloud, nList, isSlice, ih, ic, interIce, water, clath,
                    interClath, total);

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
        QlmTotal.ptq[iatom].ylm = sph::lookupTableQ6Vec(angles);
        first = false;
        continue;
      }
      yl = sph::spheriHarmo(l, angles);
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
      if (jatomIndex < 0 || jatomIndex >= yCloud.nop) {
        continue;
      }
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
      if (jatomIndex < 0 || jatomIndex >= yCloud.nop) {
        continue;
      }
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
