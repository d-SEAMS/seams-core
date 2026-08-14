//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#ifndef SEAMS_STRUCTURE_DESC_H_
#define SEAMS_STRUCTURE_DESC_H_

#include <mol_sys.hpp>

#include <string>
#include <vector>

/** @file structure_desc.hpp
 *  @brief Voronoi-weighted order, IRA template overlay, SOAP, and a
 *   linear prototype classifier.
 */

namespace chill {

enum class CrystalKind { other = 0, sc, fcc, hcp, bcc };

struct TemplateHit {
  CrystalKind kind = CrystalKind::other;
  double rmsd = 1e300;
  const char *name = "other";
};

//! Overlay the k nearest neighbours of each particle onto FCC, HCP, BCC and
//! SC shells. Uses IRA when linked; otherwise Horn on a distance-sorted
//! correspondence (high-symmetry shells only).
[[nodiscard]] std::vector<TemplateHit> classifyTemplates(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, int kNeigh = 12);

//! Bartok SOAP power spectrum for one particle (nMax radial Gaussians,
//! spherical harmonics through lMax). Length is nMax*nMax*(lMax+1).
[[nodiscard]] std::vector<double> soapSpectrum(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int iatom,
    const std::vector<std::vector<int>> &nList, int nMax, int lMax,
    double rcut);

//! Bartok SOAP power spectrum for every particle. Length is nop; each
//! row is nMax*nMax*(lMax+1).
[[nodiscard]] std::vector<std::vector<double>> soapSpectrumAll(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, int nMax, int lMax,
    double rcut);

//! One-versus-rest ridge classifier on a dense feature matrix.
struct LinearClassifier {
  int nClasses = 0;
  int nFeat = 0;
  double ridge = 1e-6;
  std::vector<double> weights; // nClasses * nFeat, row-major
  std::vector<std::string> labels;

  void fit(const std::vector<std::vector<double>> &X,
           const std::vector<int> &y);
  [[nodiscard]] int predict(const std::vector<double> &x) const;
};

//! [q4, q6, q8] from the Voronoi-weighted Steinhardt path (Mickel).
[[nodiscard]] std::vector<double> voronoiFeature(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int iatom,
    double candidateCutoff);

//! [q4, q6, q8] for every particle. One tessellation per order (l = 4, 6, 8).
[[nodiscard]] std::vector<std::vector<double>> voronoiFeatures(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double candidateCutoff);

} // namespace chill

#endif // SEAMS_STRUCTURE_DESC_H_
