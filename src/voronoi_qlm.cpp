//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <voronoi_qlm.hpp>

#include <generic.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

namespace {

using Vec3 = std::array<double, 3>;

/**
 * @details One facet of a Voronoi cell, computed in the plane that carries
 *  it. The facet toward neighbour j lies on the perpendicular bisector of the
 *  displacement to j; it is the bisector's large in-plane square clipped
 *  against the half-space of every other candidate's bisector. Working facet
 *  by facet keeps the geometry two-dimensional: each clip is a
 *  Sutherland-Hodgman pass of a convex polygon against a half-plane, which
 *  avoids maintaining full polyhedron topology and its degeneracies.
 * @param[in] target Index into @a disp of the neighbour whose facet is wanted.
 * @param[in] disp Minimum-image displacements to every candidate.
 * @param[in] halfExtent Half-width of the starting square on the plane.
 * @return The facet area; zero when the candidates close the facet off.
 */
double facetArea(size_t target, const std::vector<Vec3> &disp,
                 double halfExtent) {
  const Vec3 &d = disp[target];
  const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  if (r <= 0.0) {
    return 0.0;
  }
  const Vec3 u = {d[0] / r, d[1] / r, d[2] / r};

  // Orthonormal basis in the bisector plane
  Vec3 e1 = (std::fabs(u[0]) < 0.9) ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
  const double proj = e1[0] * u[0] + e1[1] * u[1] + e1[2] * u[2];
  for (int k = 0; k < 3; k++) {
    e1[k] -= proj * u[k];
  }
  const double n1 =
      std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  for (int k = 0; k < 3; k++) {
    e1[k] /= n1;
  }
  const Vec3 e2 = {u[1] * e1[2] - u[2] * e1[1], u[2] * e1[0] - u[0] * e1[2],
                   u[0] * e1[1] - u[1] * e1[0]};

  // Starting polygon: a large square about the plane's foot point (r/2) u
  std::vector<std::array<double, 2>> poly = {{-halfExtent, -halfExtent},
                                             {halfExtent, -halfExtent},
                                             {halfExtent, halfExtent},
                                             {-halfExtent, halfExtent}};
  const double h = 0.5 * r;

  std::vector<std::array<double, 2>> next;
  for (size_t k = 0; k < disp.size(); k++) {
    if (k == target || poly.empty()) {
      continue;
    }
    const Vec3 &dk = disp[k];
    const double rk =
        std::sqrt(dk[0] * dk[0] + dk[1] * dk[1] + dk[2] * dk[2]);
    if (rk <= 0.0) {
      continue;
    }
    // In-plane constraint: (h u + a e1 + b e2) . dk <= rk^2 / 2
    const double c0 = h * (u[0] * dk[0] + u[1] * dk[1] + u[2] * dk[2]);
    const double ca = e1[0] * dk[0] + e1[1] * dk[1] + e1[2] * dk[2];
    const double cb = e2[0] * dk[0] + e2[1] * dk[1] + e2[2] * dk[2];
    const double rhs = 0.5 * rk * rk - c0;

    next.clear();
    const size_t nv = poly.size();
    for (size_t v = 0; v < nv; v++) {
      const auto &p = poly[v];
      const auto &q = poly[(v + 1) % nv];
      const double fp = ca * p[0] + cb * p[1] - rhs;
      const double fq = ca * q[0] + cb * q[1] - rhs;
      const bool inP = fp <= 1e-12;
      const bool inQ = fq <= 1e-12;
      if (inP) {
        next.push_back(p);
      }
      if (inP != inQ) {
        const double t = fp / (fp - fq);
        next.push_back({p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])});
      }
    }
    poly.swap(next);
  }

  if (poly.size() < 3) {
    return 0.0;
  }
  double area = 0.0;
  for (size_t v = 0; v < poly.size(); v++) {
    const auto &p = poly[v];
    const auto &q = poly[(v + 1) % poly.size()];
    area += p[0] * q[1] - q[0] * p[1];
  }
  return 0.5 * std::fabs(area);
}

} // namespace

std::vector<chill::VoronoiWeights> chill::voronoiFacetWeights(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double candidateCutoff) {
  std::vector<chill::VoronoiWeights> result(yCloud.nop);
  if (yCloud.box.size() < 3 || candidateCutoff <= 0.0) {
    return result;
  }
  const double cutoffSq = candidateCutoff * candidateCutoff;

  for (int i = 0; i < yCloud.nop; i++) {
    std::vector<Vec3> disp;
    std::vector<int> who;
    for (int j = 0; j < yCloud.nop; j++) {
      if (j == i) {
        continue;
      }
      const auto d = gen::relDist(yCloud, i, j);
      const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
      if (r2 > 0.0 && r2 <= cutoffSq) {
        disp.push_back({d[0], d[1], d[2]});
        who.push_back(j);
      }
    }

    double total = 0.0;
    std::vector<double> areas(disp.size(), 0.0);
    for (size_t t = 0; t < disp.size(); t++) {
      areas[t] = facetArea(t, disp, candidateCutoff);
      total += areas[t];
    }
    if (total <= 0.0) {
      continue;
    }
    for (size_t t = 0; t < disp.size(); t++) {
      // Facets thinner than one part in 1e9 of the cell surface are clipping
      // residue, not neighbours
      if (areas[t] > 1e-9 * total) {
        result[i].neighbours.push_back(who[t]);
        result[i].weights.push_back(areas[t] / total);
      }
    }
    // Renormalise after dropping residue
    double kept = 0.0;
    for (const double w : result[i].weights) {
      kept += w;
    }
    for (double &w : result[i].weights) {
      w /= kept;
    }
  }
  return result;
}

chill::SteinhardtQl chill::steinhardtQlVoronoi(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double candidateCutoff, int orderL) {
  chill::SteinhardtQl result;
  result.ql.assign(yCloud.nop, 0.0);
  result.qlBar.assign(yCloud.nop, 0.0);
  if (orderL != 3 && orderL != 4 && orderL != 6 && orderL != 8) {
    return result;
  }
  const auto cells = chill::voronoiFacetWeights(yCloud, candidateCutoff);
  const int nComp = 2 * orderL + 1;
  const double prefactor = 4.0 * std::numbers::pi / (2.0 * orderL + 1.0);

  // Weighted q_lm(i) = sum_j w_ij Y_lm(u_ij), the facet weights replacing the
  // 1/N_b of the unweighted parameter
  std::vector<std::complex<double>> qlm(
      static_cast<size_t>(yCloud.nop) * nComp, {0.0, 0.0});
  for (int i = 0; i < yCloud.nop; i++) {
    const auto &cell = cells[i];
    for (size_t k = 0; k < cell.neighbours.size(); k++) {
      const int j = cell.neighbours[k];
      const auto d = gen::relDist(yCloud, i, j);
      const double r =
          std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
      const std::array<double, 2> angles = {std::atan2(d[0], d[1]),
                                            std::acos(d[2] / r)};
      const auto ylm = sph::spheriHarmo(orderL, angles);
      for (int m = 0; m < nComp; m++) {
        qlm[static_cast<size_t>(i) * nComp + m] += cell.weights[k] * ylm[m];
      }
    }
  }

  for (int i = 0; i < yCloud.nop; i++) {
    const size_t row = static_cast<size_t>(i) * nComp;
    double sumLocal = 0.0;
    for (int m = 0; m < nComp; m++) {
      sumLocal += std::norm(qlm[row + m]);
    }
    result.ql[i] = std::sqrt(prefactor * sumLocal);

    // Lechner-Dellago average of the weighted q_lm over the particle and its
    // facet neighbours
    std::vector<std::complex<double>> bar(qlm.begin() + row,
                                          qlm.begin() + row + nComp);
    int contributing = 1;
    for (const int j : cells[i].neighbours) {
      const size_t jRow = static_cast<size_t>(j) * nComp;
      for (int m = 0; m < nComp; m++) {
        bar[m] += qlm[jRow + m];
      }
      contributing++;
    }
    double sumBar = 0.0;
    for (int m = 0; m < nComp; m++) {
      bar[m] /= static_cast<double>(contributing);
      sumBar += std::norm(bar[m]);
    }
    result.qlBar[i] = std::sqrt(prefactor * sumBar);
  }
  return result;
}
