//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <structure_desc.hpp>

#include <absOrientation.hpp>
#include <bop.hpp>
#include <generic.hpp>
#include <ira_sofi.hpp>
#include <voronoi_qlm.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <numbers>
#include <utility>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
using Vec3 = std::array<double, 3>;

Eigen::MatrixXd shellOf(const Cloud &yCloud,
                        const std::vector<std::vector<int>> &nList, int iatom,
                        int kNeigh) {
  std::vector<std::pair<double, Vec3>> cand;
  if (iatom < 0 || iatom >= static_cast<int>(nList.size())) {
    return {};
  }
  for (size_t j = 1; j < nList[iatom].size(); j++) {
    const int id = nList[iatom][j];
    const auto it = yCloud.idIndexMap.find(id);
    if (it == yCloud.idIndexMap.end()) {
      continue;
    }
    const auto d = gen::relDist(yCloud, iatom, it->second);
    const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    if (r2 > 0.0) {
      cand.push_back({r2, {d[0], d[1], d[2]}});
    }
  }
  const int keep = std::min(kNeigh, static_cast<int>(cand.size()));
  if (keep == 0) {
    return {};
  }
  std::partial_sort(cand.begin(), cand.begin() + keep, cand.end(),
                    [](const auto &a, const auto &b) { return a.first < b.first; });
  Eigen::MatrixXd m(keep, 3);
  for (int k = 0; k < keep; k++) {
    m(k, 0) = cand[static_cast<size_t>(k)].second[0];
    m(k, 1) = cand[static_cast<size_t>(k)].second[1];
    m(k, 2) = cand[static_cast<size_t>(k)].second[2];
  }
  return m;
}

void scaleToUnit(Eigen::MatrixXd &m) {
  double mean = 0.0;
  for (int i = 0; i < m.rows(); i++) {
    mean += m.row(i).norm();
  }
  mean /= static_cast<double>(m.rows());
  if (mean > 0.0) {
    m /= mean;
  }
}

void sortRows(Eigen::MatrixXd &m) {
  std::vector<int> order(static_cast<size_t>(m.rows()));
  for (int i = 0; i < m.rows(); i++) {
    order[static_cast<size_t>(i)] = i;
  }
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    if (std::fabs(m(a, 2) - m(b, 2)) > 1e-9) {
      return m(a, 2) < m(b, 2);
    }
    if (std::fabs(m(a, 1) - m(b, 1)) > 1e-9) {
      return m(a, 1) < m(b, 1);
    }
    return m(a, 0) < m(b, 0);
  });
  Eigen::MatrixXd out(m.rows(), 3);
  for (int i = 0; i < m.rows(); i++) {
    out.row(i) = m.row(order[static_cast<size_t>(i)]);
  }
  m.swap(out);
}

Eigen::MatrixXd fcc12() {
  Eigen::MatrixXd m(12, 3);
  const double s = 1.0 / std::sqrt(2.0);
  int r = 0;
  for (int a = -1; a <= 1; a += 2) {
    for (int b = -1; b <= 1; b += 2) {
      m.row(r++) << a * s, b * s, 0.0;
      m.row(r++) << a * s, 0.0, b * s;
      m.row(r++) << 0.0, a * s, b * s;
    }
  }
  scaleToUnit(m);
  return m;
}

Eigen::MatrixXd hcp12() {
  Eigen::MatrixXd m(12, 3);
  const double a = 1.0;
  const double c = std::sqrt(8.0 / 3.0);
  int r = 0;
  for (int k = 0; k < 6; k++) {
    const double t = k * std::numbers::pi / 3.0;
    m.row(r++) << a * std::cos(t), a * std::sin(t), 0.0;
  }
  for (int sign : { -1, 1 }) {
    for (int k = 0; k < 3; k++) {
      const double t =
          k * 2.0 * std::numbers::pi / 3.0 + std::numbers::pi / 6.0;
      m.row(r++) << (a / std::sqrt(3.0)) * std::cos(t),
          (a / std::sqrt(3.0)) * std::sin(t), sign * 0.5 * c;
    }
  }
  scaleToUnit(m);
  return m;
}

Eigen::MatrixXd bcc8() {
  Eigen::MatrixXd m(8, 3);
  int r = 0;
  for (int x : { -1, 1 }) {
    for (int y : { -1, 1 }) {
      for (int z : { -1, 1 }) {
        m.row(r++) << x, y, z;
      }
    }
  }
  scaleToUnit(m);
  return m;
}

Eigen::MatrixXd sc6() {
  Eigen::MatrixXd m(6, 3);
  m << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
  return m;
}

double overlayRmsd(Eigen::MatrixXd ref, Eigen::MatrixXd tgt) {
  if (ref.rows() != tgt.rows() || ref.rows() == 0) {
    return 1e300;
  }
  scaleToUnit(ref);
  scaleToUnit(tgt);
  if (ira::available()) {
    ira::Match hit;
    if (ira::match(ref, tgt, hit) == 0) {
      return hit.rmsd;
    }
  }
  sortRows(ref);
  sortRows(tgt);
  std::vector<double> quat, perAtom;
  double rmsd = 0.0;
  double scale = 1.0;
  if (absor::hornAbsOrientation(ref, tgt, quat, rmsd, perAtom, scale) != 0) {
    return 1e300;
  }
  return rmsd;
}

} // namespace

std::vector<chill::TemplateHit> chill::classifyTemplates(
    const Cloud &yCloud, const std::vector<std::vector<int>> &nList,
    int kNeigh) {
  std::vector<TemplateHit> out(yCloud.nop);
  const struct {
    CrystalKind kind;
    const char *name;
    Eigen::MatrixXd (*make)();
    int k;
  } templates[] = {{CrystalKind::fcc, "fcc", fcc12, 12},
                   {CrystalKind::hcp, "hcp", hcp12, 12},
                   {CrystalKind::bcc, "bcc", bcc8, 8},
                   {CrystalKind::sc, "sc", sc6, 6}};

  for (int i = 0; i < yCloud.nop; i++) {
    TemplateHit best;
    for (const auto &tmpl : templates) {
      if (kNeigh > 0 && kNeigh < tmpl.k) {
        continue;
      }
      auto shell = shellOf(yCloud, nList, i, tmpl.k);
      if (shell.rows() < tmpl.k) {
        continue;
      }
      if (shell.rows() > tmpl.k) {
        shell.conservativeResize(tmpl.k, 3);
      }
      const double rmsd = overlayRmsd(tmpl.make(), shell);
      if (rmsd < best.rmsd) {
        best.kind = tmpl.kind;
        best.name = tmpl.name;
        best.rmsd = rmsd;
      }
    }
    // A first-neighbour-unit RMSD above 0.35 is not a lattice shell
    if (best.rmsd > 0.35) {
      best.kind = CrystalKind::other;
      best.name = "other";
    }
    out[static_cast<size_t>(i)] = best;
  }
  return out;
}

std::vector<double> chill::soapSpectrum(
    const Cloud &yCloud, int iatom, const std::vector<std::vector<int>> &nList,
    int nMax, int lMax, double rcut) {
  const int nComp = (lMax + 1) * (lMax + 1);
  std::vector<std::complex<double>> coeff(
      static_cast<size_t>(nMax) * static_cast<size_t>(nComp), {0.0, 0.0});
  if (iatom < 0 || iatom >= static_cast<int>(nList.size()) || nMax < 1 ||
      lMax < 0 || rcut <= 0.0) {
    return std::vector<double>(static_cast<size_t>(nMax * nMax * (lMax + 1)),
                               0.0);
  }
  const double sigma = rcut / static_cast<double>(nMax);
  auto addNeighbour = [&](const Vec3 &d) {
    const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (r <= 0.0 || r >= rcut) {
      return;
    }
    const std::array<double, 2> angles = {std::atan2(d[0], d[1]),
                                          std::acos(d[2] / r)};
    for (int n = 0; n < nMax; n++) {
      const double rn = (n + 0.5) * rcut / static_cast<double>(nMax);
      const double g = std::exp(-((r - rn) / sigma) * ((r - rn) / sigma));
      for (int l = 0; l <= lMax; l++) {
        const int base = l * l;
        if (l == 0) {
          coeff[static_cast<size_t>(n) * nComp + base] +=
              g * (0.5 / std::sqrt(std::numbers::pi));
          continue;
        }
        if (l != 3 && l != 4 && l != 6 && l != 8) {
          continue;
        }
        const auto ylm = sph::spheriHarmo(l, angles);
        for (int m = 0; m < 2 * l + 1; m++) {
          coeff[static_cast<size_t>(n) * nComp + base + m] += g * ylm[m];
        }
      }
    }
  };
  for (size_t j = 1; j < nList[iatom].size(); j++) {
    const int id = nList[iatom][j];
    const auto it = yCloud.idIndexMap.find(id);
    if (it == yCloud.idIndexMap.end()) {
      continue;
    }
    const auto d = gen::relDist(yCloud, iatom, it->second);
    addNeighbour({d[0], d[1], d[2]});
  }

  std::vector<double> spec(static_cast<size_t>(nMax * nMax * (lMax + 1)), 0.0);
  int slot = 0;
  for (int n = 0; n < nMax; n++) {
    for (int np = 0; np < nMax; np++) {
      for (int l = 0; l <= lMax; l++) {
        const int base = l * l;
        std::complex<double> acc = 0.0;
        for (int m = 0; m < 2 * l + 1; m++) {
          acc += coeff[static_cast<size_t>(n) * nComp + base + m] *
                 std::conj(coeff[static_cast<size_t>(np) * nComp + base + m]);
        }
        spec[static_cast<size_t>(slot++)] = acc.real();
      }
    }
  }
  return spec;
}

void chill::LinearClassifier::fit(const std::vector<std::vector<double>> &X,
                                  const std::vector<int> &y) {
  nFeat = X.empty() ? 0 : static_cast<int>(X[0].size());
  nClasses = 0;
  for (int lab : y) {
    nClasses = std::max(nClasses, lab + 1);
  }
  weights.assign(static_cast<size_t>(nClasses * nFeat), 0.0);
  if (nFeat == 0 || nClasses == 0 || X.size() != y.size()) {
    return;
  }
  const int n = static_cast<int>(X.size());
  Eigen::MatrixXd A(n, nFeat);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < nFeat; j++) {
      A(i, j) = X[static_cast<size_t>(i)][static_cast<size_t>(j)];
    }
  }
  const Eigen::MatrixXd gram =
      A.transpose() * A +
      ridge * Eigen::MatrixXd::Identity(nFeat, nFeat);
  const Eigen::LDLT<Eigen::MatrixXd> dec(gram);
  for (int c = 0; c < nClasses; c++) {
    Eigen::VectorXd t(n);
    for (int i = 0; i < n; i++) {
      t(i) = (y[static_cast<size_t>(i)] == c) ? 1.0 : 0.0;
    }
    const Eigen::VectorXd w = dec.solve(A.transpose() * t);
    for (int j = 0; j < nFeat; j++) {
      weights[static_cast<size_t>(c * nFeat + j)] = w(j);
    }
  }
}

int chill::LinearClassifier::predict(const std::vector<double> &x) const {
  if (nClasses == 0 || static_cast<int>(x.size()) != nFeat) {
    return 0;
  }
  int best = 0;
  double bestScore = -1e300;
  for (int c = 0; c < nClasses; c++) {
    double s = 0.0;
    for (int j = 0; j < nFeat; j++) {
      s += weights[static_cast<size_t>(c * nFeat + j)] *
           x[static_cast<size_t>(j)];
    }
    if (s > bestScore) {
      bestScore = s;
      best = c;
    }
  }
  return best;
}

std::vector<double> chill::voronoiFeature(const Cloud &yCloud, int iatom,
                                          double candidateCutoff) {
  std::vector<double> feat(3, 0.0);
  if (iatom < 0 || iatom >= yCloud.nop) {
    return feat;
  }
  const int orders[3] = {4, 6, 8};
  for (int k = 0; k < 3; k++) {
    const auto q =
        chill::steinhardtQlVoronoi(yCloud, candidateCutoff, orders[k]);
    feat[static_cast<size_t>(k)] = q.ql[static_cast<size_t>(iatom)];
  }
  return feat;
}
