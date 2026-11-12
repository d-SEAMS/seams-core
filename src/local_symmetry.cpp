//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <local_symmetry.hpp>

#include <ira_sofi.hpp>
#include <seams_c_api.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <new>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_SPGLIB
#include <spglib.h>
#endif

namespace {

constexpr double kBond = 2.75;
constexpr double kPi = 3.14159265358979323846;

Eigen::Vector3d tet(int i) {
  const double s = kBond / std::sqrt(3.0);
  const double sgn[4][3] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
  return {s * sgn[i][0], s * sgn[i][1], s * sgn[i][2]};
}

Eigen::Vector3d rotateAround(const Eigen::Vector3d &axis, const Eigen::Vector3d &v,
                             double ang) {
  const Eigen::Vector3d u = axis.normalized();
  const double c = std::cos(ang);
  const double s = std::sin(ang);
  return v * c + u.cross(v) * s + u * u.dot(v) * (1.0 - c);
}

Eigen::MatrixXd fromVectors(const std::vector<Eigen::Vector3d> &pts) {
  Eigen::MatrixXd m(static_cast<Eigen::Index>(pts.size()), 3);
  for (size_t i = 0; i < pts.size(); i++) {
    m.row(static_cast<Eigen::Index>(i)) = pts[i].transpose();
  }
  return m;
}

Eigen::MatrixXd iceShellKind(int nEclipsed) {
  std::vector<Eigen::Vector3d> pts;
  pts.emplace_back(0, 0, 0);
  Eigen::Vector3d v[4];
  for (int i = 0; i < 4; i++) {
    v[i] = tet(i);
    pts.push_back(v[i]);
  }
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (i == j) {
        continue;
      }
      Eigen::Vector3d second = v[i] - v[j];
      if (i < nEclipsed) {
        const Eigen::Vector3d out = rotateAround(v[i], second - v[i], kPi / 3.0);
        second = v[i] + out;
      }
      pts.push_back(second);
    }
  }
  return fromVectors(pts);
}

Eigen::MatrixXd regularRing(int n) {
  const double radius = kBond / (2.0 * std::sin(kPi / static_cast<double>(n)));
  Eigen::MatrixXd m(n, 3);
  for (int k = 0; k < n; k++) {
    const double ang = 2.0 * kPi * static_cast<double>(k) / static_cast<double>(n);
    m(k, 0) = radius * std::cos(ang);
    m(k, 1) = radius * std::sin(ang);
    m(k, 2) = 0.0;
  }
  return m;
}

Eigen::MatrixXd dodecahedron() {
  const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
  const double scale = kBond / (2.0 / phi);
  std::vector<Eigen::Vector3d> pts;
  for (int sx = -1; sx <= 1; sx += 2) {
    for (int sy = -1; sy <= 1; sy += 2) {
      for (int sz = -1; sz <= 1; sz += 2) {
        pts.emplace_back(sx, sy, sz);
      }
    }
  }
  for (int a = -1; a <= 1; a += 2) {
    for (int b = -1; b <= 1; b += 2) {
      const double v[3] = {0.0, static_cast<double>(a) / phi, static_cast<double>(b) * phi};
      for (int r = 0; r < 3; r++) {
        pts.emplace_back(v[r] * scale, v[(r + 1) % 3] * scale, v[(r + 2) % 3] * scale);
      }
    }
  }
  for (size_t i = 0; i < 8; i++) {
    pts[i] *= scale;
  }
  return fromVectors(pts);
}

Eigen::MatrixXd cage51262() {
  // Equal-edge D6d tetrakaidecahedron. The belt radius is the root of the
  // roundness condition derived in the construction (edge = bond).
  auto gap = [](double x) {
    const double u2 = 2.0 * x - x * x;
    const double v2 = (1.0 - (2.0 - std::sqrt(3.0)) * x * x) / 4.0;
    const double u = std::sqrt(u2);
    const double v = std::sqrt(v2);
    return u * v - (x * x - x - 0.5);
  };
  double lo = 1.01;
  double hi = 1.9;
  for (int it = 0; it < 60; it++) {
    const double mid = 0.5 * (lo + hi);
    if (gap(mid) > 0.0) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  const double x = 0.5 * (lo + hi);
  const double u = std::sqrt(2.0 * x - x * x);
  const double v = 0.5 * std::sqrt(1.0 - (2.0 - std::sqrt(3.0)) * x * x);
  const double R = kBond;
  const double P = x * kBond;
  const double H = (u + v) * kBond;
  const double h = v * kBond;
  std::vector<Eigen::Vector3d> pts;
  pts.reserve(24);
  for (int k = 0; k < 6; k++) {
    const double top = k * kPi / 3.0;
    pts.emplace_back(R * std::cos(top), R * std::sin(top), H);
    pts.emplace_back(P * std::cos(top), P * std::sin(top), h);
    const double bot = top + kPi / 6.0;
    pts.emplace_back(P * std::cos(bot), P * std::sin(bot), -h);
    pts.emplace_back(R * std::cos(bot), R * std::sin(bot), -H);
  }
  return fromVectors(pts);
}

std::vector<Eigen::Matrix3d> tetrahedralGroup() {
  std::vector<Eigen::Matrix3d> ops;
  const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  for (int p = 0; p < 6; p++) {
    for (int s0 = -1; s0 <= 1; s0 += 2) {
      for (int s1 = -1; s1 <= 1; s1 += 2) {
        for (int s2 = -1; s2 <= 1; s2 += 2) {
          if (s0 * s1 * s2 != 1) {
            continue;
          }
          Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
          const int sgn[3] = {s0, s1, s2};
          for (int i = 0; i < 3; i++) {
            M(i, perm[p][i]) = static_cast<double>(sgn[i]);
          }
          ops.push_back(M);
        }
      }
    }
  }
  return ops;
}

Eigen::MatrixXd cage51264() {
  // One seed from each Td orbit of the 5^12 6^4 graph: the apex on a
  // threefold axis, a hexagon corner, and an edge site. The other vertices
  // are the images under the tetrahedral group. Seeds are the centred
  // equal-edge embedding of that graph.
  const Eigen::Vector3d seeds[3] = {
      {2.74568415, 2.74568415, 2.74568415},
      {0.16007494, -3.26552686, -3.26552686},
      {-0.95785744, -0.95785744, -4.36797099}};
  std::vector<Eigen::Vector3d> pts;
  for (const auto &M : tetrahedralGroup()) {
    for (const auto &s : seeds) {
      const Eigen::Vector3d v = M * s;
      bool seen = false;
      for (const auto &p : pts) {
        if ((p - v).norm() < 1e-6) {
          seen = true;
          break;
        }
      }
      if (!seen) {
        pts.push_back(v);
      }
    }
  }
  return fromVectors(pts);
}

Eigen::Vector3d centroidOf(const Eigen::MatrixXd &m) {
  Eigen::Vector3d c = Eigen::Vector3d::Zero();
  for (Eigen::Index i = 0; i < m.rows(); i++) {
    c += m.row(i).transpose();
  }
  return c / static_cast<double>(m.rows());
}

Eigen::MatrixXd centered(const Eigen::MatrixXd &m, const Eigen::Vector3d &c) {
  Eigen::MatrixXd out = m;
  for (Eigen::Index i = 0; i < out.rows(); i++) {
    out.row(i) -= c.transpose();
  }
  return out;
}

std::vector<int> hungarian(const std::vector<std::vector<double>> &a) {
  const int n = static_cast<int>(a.size());
  const double inf = 1e100;
  std::vector<double> u(static_cast<size_t>(n + 1), 0.0);
  std::vector<double> v(static_cast<size_t>(n + 1), 0.0);
  std::vector<int> p(static_cast<size_t>(n + 1), 0);
  std::vector<int> way(static_cast<size_t>(n + 1), 0);
  for (int i = 1; i <= n; i++) {
    p[0] = i;
    int j0 = 0;
    std::vector<double> minv(static_cast<size_t>(n + 1), inf);
    std::vector<char> used(static_cast<size_t>(n + 1), 0);
    do {
      used[static_cast<size_t>(j0)] = 1;
      const int i0 = p[static_cast<size_t>(j0)];
      int j1 = 0;
      double delta = inf;
      for (int j = 1; j <= n; j++) {
        if (used[static_cast<size_t>(j)]) {
          continue;
        }
        const double cur = a[static_cast<size_t>(i0 - 1)][static_cast<size_t>(j - 1)] -
                           u[static_cast<size_t>(i0)] - v[static_cast<size_t>(j)];
        if (cur < minv[static_cast<size_t>(j)]) {
          minv[static_cast<size_t>(j)] = cur;
          way[static_cast<size_t>(j)] = j0;
        }
        if (minv[static_cast<size_t>(j)] < delta) {
          delta = minv[static_cast<size_t>(j)];
          j1 = j;
        }
      }
      for (int j = 0; j <= n; j++) {
        if (used[static_cast<size_t>(j)]) {
          u[static_cast<size_t>(p[static_cast<size_t>(j)])] += delta;
          v[static_cast<size_t>(j)] -= delta;
        } else {
          minv[static_cast<size_t>(j)] -= delta;
        }
      }
      j0 = j1;
    } while (p[static_cast<size_t>(j0)] != 0);
    do {
      const int j1 = way[static_cast<size_t>(j0)];
      p[static_cast<size_t>(j0)] = p[static_cast<size_t>(j1)];
      j0 = j1;
    } while (j0 != 0);
  }
  std::vector<int> assign(static_cast<size_t>(n), -1);
  for (int j = 1; j <= n; j++) {
    if (p[static_cast<size_t>(j)] != 0) {
      assign[static_cast<size_t>(p[static_cast<size_t>(j)] - 1)] = j - 1;
    }
  }
  return assign;
}

Eigen::Matrix3d kabsch(const Eigen::MatrixXd &ref, const Eigen::MatrixXd &target,
                       const std::vector<int> &assign) {
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  const int n = static_cast<int>(assign.size());
  for (int i = 0; i < n; i++) {
    const int j = assign[static_cast<size_t>(i)];
    H += target.row(j).transpose() * ref.row(i);
  }
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d U = svd.matrixU();
  Eigen::Matrix3d V = svd.matrixV();
  if ((V * U.transpose()).determinant() < 0.0) {
    V.col(2) *= -1.0;
  }
  return V * U.transpose();
}

double residual(const Eigen::MatrixXd &ref, const Eigen::MatrixXd &moved,
                const std::vector<int> &assign, double &hausdorff) {
  double sum = 0.0;
  hausdorff = 0.0;
  const int n = static_cast<int>(assign.size());
  for (int i = 0; i < n; i++) {
    const int j = assign[static_cast<size_t>(i)];
    const double d = (ref.row(i) - moved.row(j)).norm();
    sum += d * d;
    hausdorff = std::max(hausdorff, d);
  }
  return n > 0 ? std::sqrt(sum / static_cast<double>(n)) : 0.0;
}

std::vector<Eigen::Matrix3d> octahedralRotations() {
  std::vector<Eigen::Matrix3d> ops;
  const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  for (int p = 0; p < 6; p++) {
    for (int s0 = -1; s0 <= 1; s0 += 2) {
      for (int s1 = -1; s1 <= 1; s1 += 2) {
        for (int s2 = -1; s2 <= 1; s2 += 2) {
          Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
          const int sgn[3] = {s0, s1, s2};
          for (int i = 0; i < 3; i++) {
            M(i, perm[p][i]) = static_cast<double>(sgn[i]);
          }
          if (M.determinant() > 0.5) {
            ops.push_back(M);
          }
        }
      }
    }
  }
  return ops;
}

symm::Overlay kabschMatch(const Eigen::MatrixXd &ref, const Eigen::MatrixXd &target) {
  symm::Overlay out;
  const int n = static_cast<int>(ref.rows());
  const Eigen::Vector3d cr = centroidOf(ref);
  const Eigen::Vector3d ct = centroidOf(target);
  const Eigen::MatrixXd A = centered(ref, cr);
  const Eigen::MatrixXd B = centered(target, ct);
  double bestRms = std::numeric_limits<double>::infinity();
  for (const auto &seed : octahedralRotations()) {
    Eigen::Matrix3d R = seed;
    std::vector<int> assign(static_cast<size_t>(n), 0);
    double haus = 0.0;
    double rms = 0.0;
    for (int iter = 0; iter < 6; iter++) {
      Eigen::MatrixXd moved(n, 3);
      for (int j = 0; j < n; j++) {
        moved.row(j) = (R * B.row(j).transpose()).transpose();
      }
      std::vector<std::vector<double>> cost(static_cast<size_t>(n),
                                            std::vector<double>(static_cast<size_t>(n), 0.0));
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          cost[static_cast<size_t>(i)][static_cast<size_t>(j)] =
              (A.row(i) - moved.row(j)).squaredNorm();
        }
      }
      assign = hungarian(cost);
      R = kabsch(A, B, assign);
      for (int j = 0; j < n; j++) {
        moved.row(j) = (R * B.row(j).transpose()).transpose();
      }
      rms = residual(A, moved, assign, haus);
    }
    if (rms < bestRms) {
      bestRms = rms;
      out.rmsd = rms;
      out.hausdorff = haus;
      out.rotation = R;
      out.translation = cr - R * ct;
      out.permutation = assign;
      out.status = 0;
      out.backend = "kabsch";
    }
  }
  return out;
}

int zeroShift(const std::vector<int> &raw, int n) {
  bool hasZero = false;
  bool hasN = false;
  for (int p : raw) {
    if (p == 0) {
      hasZero = true;
    }
    if (p == n) {
      hasN = true;
    }
  }
  return (hasN && !hasZero) ? 1 : 0;
}

symm::Overlay fromIra(const Eigen::MatrixXd &ref, const Eigen::MatrixXd &target) {
  symm::Overlay out;
  ira::Match hit;
  if (ira::match(ref, target, hit) != 0) {
    return out;
  }
  const int n = static_cast<int>(std::min(ref.rows(), target.rows()));
  const int shift = zeroShift(hit.assignment, n);
  out.permutation.assign(static_cast<size_t>(n), -1);
  for (int i = 0; i < n && i < static_cast<int>(hit.assignment.size()); i++) {
    out.permutation[static_cast<size_t>(i)] = hit.assignment[static_cast<size_t>(i)] - shift;
  }
  out.rotation = hit.rotation;
  out.translation = hit.translation;
  out.rmsd = hit.rmsd;
  out.hausdorff = hit.hausdorff;
  out.status = 0;
  out.backend = "ira";
  return out;
}

std::string trimSymbol(const char *text) {
  std::string s = text == nullptr ? "" : text;
  while (!s.empty() && (s.back() == ' ' || s.back() == '\0')) {
    s.pop_back();
  }
  return s;
}

Eigen::Vector3d fracOf(const Eigen::Matrix3d &vectors, const Eigen::Vector3d &cart) {
  Eigen::Vector3d f = vectors.transpose().inverse() * cart;
  for (int k = 0; k < 3; k++) {
    f[k] -= std::floor(f[k]);
  }
  return f;
}

void writeFracCrystal(symm::Lattice &cell, const Eigen::MatrixXd &frac) {
  cell.cartesian.resize(frac.rows(), 3);
  for (Eigen::Index i = 0; i < frac.rows(); i++) {
    Eigen::Vector3d f = frac.row(i).transpose();
    for (int k = 0; k < 3; k++) {
      f[k] -= std::floor(f[k]);
    }
    cell.cartesian.row(i) = (cell.vectors.transpose() * f).transpose();
  }
}

bool wantIra(const char *backend) {
  if (backend != nullptr && std::strcmp(backend, "kabsch") == 0) {
    return false;
  }
  if (backend != nullptr && std::strcmp(backend, "ira") == 0) {
    return true;
  }
  return ira::available();
}

} // namespace

namespace symm {

const char *name(Kind kind) {
  switch (kind) {
  case Kind::iceIh:
    return "ice-Ih";
  case Kind::iceIc:
    return "ice-Ic";
  case Kind::iceIsd:
    return "ice-Isd";
  case Kind::cage512:
    return "5^12";
  case Kind::cage51262:
    return "5^12 6^2";
  case Kind::cage51264:
    return "5^12 6^4";
  case Kind::ring5:
    return "ring-5";
  case Kind::ring6:
    return "ring-6";
  }
  return "unknown";
}

Eigen::MatrixXd points(Kind kind) {
  switch (kind) {
  case Kind::iceIh:
    return iceShellKind(1);
  case Kind::iceIc:
    return iceShellKind(0);
  case Kind::iceIsd:
    return iceShellKind(2);
  case Kind::cage512:
    return dodecahedron();
  case Kind::cage51262:
    return cage51262();
  case Kind::cage51264:
    return cage51264();
  case Kind::ring5:
    return regularRing(5);
  case Kind::ring6:
    return regularRing(6);
  }
  return {};
}

std::vector<Overlay> score(const Eigen::MatrixXd &environment, const char *backend) {
  std::vector<Overlay> hits;
  if (environment.cols() != 3 || environment.rows() < 1) {
    return hits;
  }
  const bool useIra = wantIra(backend);
  if (backend != nullptr && std::strcmp(backend, "ira") == 0 && !ira::available()) {
    return hits;
  }
  for (int id = 0; id < kTemplateCount; id++) {
    const auto kind = static_cast<Kind>(id);
    const Eigen::MatrixXd ref = points(kind);
    if (ref.rows() != environment.rows()) {
      continue;
    }
    Overlay hit = useIra ? fromIra(ref, environment) : kabschMatch(ref, environment);
    if (!useIra && backend != nullptr && std::strcmp(backend, "ira") == 0) {
      continue;
    }
    if (useIra && hit.status != 0) {
      hit = kabschMatch(ref, environment);
    }
    hit.kind = kind;
    hit.name = name(kind);
    if (hit.status == 0) {
      hits.push_back(std::move(hit));
    }
  }
  std::sort(hits.begin(), hits.end(), [](const Overlay &a, const Overlay &b) {
    if (a.rmsd != b.rmsd) {
      return a.rmsd < b.rmsd;
    }
    return static_cast<int>(a.kind) < static_cast<int>(b.kind);
  });
  return hits;
}

Overlay best(const Eigen::MatrixXd &environment, const char *backend) {
  if (backend == nullptr || (std::strcmp(backend, "ira") != 0 &&
                             std::strcmp(backend, "kabsch") != 0)) {
    return matchKernel(environment);
  }
  const auto hits = score(environment, backend);
  if (hits.empty()) {
    return {};
  }
  return hits.front();
}

Group pointGroup(const Eigen::MatrixXd &xyz, double symThr) {
  Group out;
  if (!ira::available()) {
    out.status = 2;
    return out;
  }
  ira::PointGroup pg;
  out.status = ira::pointGroup(xyz, pg, symThr);
  if (out.status == 0) {
    out.symbol = trimSymbol(pg.symbol.c_str());
    out.nOperations = pg.nOperations;
  }
  return out;
}

Space spaceGroup(const Eigen::Matrix3d &vectors, const Eigen::MatrixXd &cartesian,
                 double symprec) {
  Space out;
#ifndef SEAMS_HAS_SPGLIB
  (void)vectors;
  (void)cartesian;
  (void)symprec;
  out.status = 2;
  return out;
#else
  if (cartesian.cols() != 3 || cartesian.rows() < 1 || vectors.determinant() == 0.0) {
    return out;
  }
  const int n = static_cast<int>(cartesian.rows());
  double lattice[3][3];
  for (int comp = 0; comp < 3; comp++) {
    for (int vec = 0; vec < 3; vec++) {
      lattice[comp][vec] = vectors(vec, comp);
    }
  }
  std::vector<double> raw(static_cast<size_t>(n) * 3, 0.0);
  std::vector<int> types(static_cast<size_t>(n), 1);
  for (int i = 0; i < n; i++) {
    const Eigen::Vector3d f =
        fracOf(vectors, cartesian.row(i).transpose());
    raw[static_cast<size_t>(3 * i)] = f[0];
    raw[static_cast<size_t>(3 * i + 1)] = f[1];
    raw[static_cast<size_t>(3 * i + 2)] = f[2];
  }
  // spglib wants position[][3]. The buffer is contiguous in that order.
  SpglibDataset *data = spg_get_dataset(
      lattice, reinterpret_cast<double(*)[3]>(raw.data()), types.data(), n, symprec);
  if (data == nullptr) {
    out.status = 1;
    return out;
  }
  out.status = 0;
  out.number = data->spacegroup_number;
  out.symbol = trimSymbol(data->international_symbol);
  spg_free_dataset(data);
  return out;
#endif
}

bool conventionalCrystalSpace(const Eigen::Matrix3d &vectors,
                              const Eigen::MatrixXd &cartesian,
                              const std::vector<int> &kind, double symprec,
                              Space &out) {
  out = {};
  out.status = 1;
  const int n = static_cast<int>(cartesian.rows());
  if (!spglibLinked() || n < 4 || static_cast<int>(kind.size()) != n ||
      cartesian.cols() != 3) {
    return false;
  }
  const int which = kind[0];
  if (which != static_cast<int>(Kind::iceIh) && which != static_cast<int>(Kind::iceIc)) {
    return false;
  }
  for (int k : kind) {
    if (k != which) {
      return false;
    }
  }
  const double skew = std::abs(vectors(0, 1)) + std::abs(vectors(0, 2)) +
                      std::abs(vectors(1, 0)) + std::abs(vectors(1, 2)) +
                      std::abs(vectors(2, 0)) + std::abs(vectors(2, 1));
  if (skew > 1e-6) {
    return false;
  }
  const double Lx = vectors(0, 0);
  const double Ly = vectors(1, 1);
  const double Lz = vectors(2, 2);
  if (!(Lx > 0.0 && Ly > 0.0 && Lz > 0.0)) {
    return false;
  }
  if (std::abs(Lx - Ly) > 1e-4 * Lx || std::abs(Ly - Lz) > 1e-4 * Lx) {
    return false;
  }
  // Ice Ic has eight molecules in the conventional cube. Ice Ih in a
  // cubic metric is not this cell; eight is tried only so a doubled
  // hexagonal setting that was placed in a cube can still match.
  const int bases[2] = {which == static_cast<int>(Kind::iceIc) ? 8 : 4, 8};
  for (int attempt = 0; attempt < 2; ++attempt) {
    const int per = bases[attempt];
    if (attempt == 1 && per == bases[0]) {
      continue;
    }
    if (per <= 0 || n % per != 0) {
      continue;
    }
    const int m3 = n / per;
    const int m = static_cast<int>(std::lround(std::cbrt(static_cast<double>(m3))));
    if (m < 1 || m * m * m != m3) {
      continue;
    }
    const double a = Lx / static_cast<double>(m);
    Eigen::MatrixXd small(per, 3);
    int got = 0;
    bool spill = false;
    for (int i = 0; i < n; ++i) {
      auto wrap = [](double v, double L) {
        double t = v / L;
        t -= std::floor(t);
        if (!(t < 1.0)) {
          t = 0.0;
        }
        return t * L;
      };
      const double x = wrap(cartesian(i, 0), Lx);
      const double y = wrap(cartesian(i, 1), Ly);
      const double z = wrap(cartesian(i, 2), Lz);
      if (!(x < a && y < a && z < a)) {
        continue;
      }
      if (got >= per) {
        spill = true;
        break;
      }
      small(got, 0) = x;
      small(got, 1) = y;
      small(got, 2) = z;
      ++got;
    }
    if (spill || got != per) {
      continue;
    }
    Eigen::Matrix3d cell = Eigen::Matrix3d::Zero();
    cell(0, 0) = a;
    cell(1, 1) = a;
    cell(2, 2) = a;
    const Space found = spaceGroup(cell, small, symprec);
    if (found.status != 0) {
      continue;
    }
    out = found;
    return true;
  }
  return false;
}

Lattice crystal(Crystal which) {
  Lattice cell;
  const double d = kBond;
  if (which == Crystal::iceIc) {
    const double a = 4.0 * d / std::sqrt(3.0);
    cell.vectors = Eigen::Matrix3d::Identity() * a;
    Eigen::MatrixXd frac(8, 3);
    frac << 0, 0, 0, 0.5, 0.5, 0, 0.5, 0, 0.5, 0, 0.5, 0.5, 0.25, 0.25, 0.25, 0.75, 0.75,
        0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.75;
    writeFracCrystal(cell, frac);
    return cell;
  }
  const double z = 1.0 / 16.0;
  const double c = d / 0.375;
  const double a = c / std::sqrt(8.0 / 3.0);
  cell.vectors << a, 0, 0, -0.5 * a, a * std::sqrt(3.0) / 2.0, 0, 0, 0, c;
  Eigen::MatrixXd frac(4, 3);
  frac << 1.0 / 3.0, 2.0 / 3.0, z, 2.0 / 3.0, 1.0 / 3.0, 0.5 + z, 2.0 / 3.0, 1.0 / 3.0,
      1.0 - z, 1.0 / 3.0, 2.0 / 3.0, 0.5 - z;
  if (which == Crystal::fault) {
    cell.vectors.row(2) *= 2.0;
    Eigen::MatrixXd stacked(8, 3);
    for (int i = 0; i < 4; i++) {
      stacked(i, 0) = frac(i, 0);
      stacked(i, 1) = frac(i, 1);
      stacked(i, 2) = frac(i, 2) * 0.5;
      stacked(i + 4, 0) = frac(i, 0) + 1.0 / 3.0;
      stacked(i + 4, 1) = frac(i, 1);
      stacked(i + 4, 2) = frac(i, 2) * 0.5 + 0.5;
    }
    writeFracCrystal(cell, stacked);
    return cell;
  }
  writeFracCrystal(cell, frac);
  return cell;
}

// Sixteen neighbours inside 5 A. A large orthorhombic cell does not walk
// every atom through a 5x5x5 image stencil: that is n*125 distances and
// it does not use the other cores. Cells of about 2.5 A and a halo that
// covers 5 A hold the same neighbours.
Eigen::MatrixXd iceShellOrtho(const Eigen::MatrixXd &cartesian, double Lx, double Ly, double Lz,
                              int atom) {
  const int n = static_cast<int>(cartesian.rows());
  const double cut = 5.0;
  const double cut2 = cut * cut;
  auto cellsOf = [](double L) {
    int ncell = std::max(1, static_cast<int>(std::floor(L / 2.5)));
    if (L / static_cast<double>(ncell) > 5.0) {
      ncell = std::max(1, static_cast<int>(std::ceil(L / 5.0)));
    }
    return std::max(3, ncell);
  };
  const int nx = cellsOf(Lx);
  const int ny = cellsOf(Ly);
  const int nz = cellsOf(Lz);
  const int ncell = nx * ny * nz;
  std::vector<int> head(static_cast<std::size_t>(ncell), -1);
  std::vector<int> next(static_cast<std::size_t>(n), -1);
  std::vector<double> wrapped(static_cast<std::size_t>(n) * 3);
  auto wrap = [](double v, double L) {
    double t = v / L;
    t -= std::floor(t);
    if (!(t < 1.0)) {
      t = 0.0;
    }
    return t * L;
  };
  for (int i = 0; i < n; i++) {
    const double x = wrap(cartesian(i, 0), Lx);
    const double y = wrap(cartesian(i, 1), Ly);
    const double z = wrap(cartesian(i, 2), Lz);
    wrapped[static_cast<std::size_t>(3 * i)] = x;
    wrapped[static_cast<std::size_t>(3 * i) + 1] = y;
    wrapped[static_cast<std::size_t>(3 * i) + 2] = z;
    const int ix = std::min(nx - 1, static_cast<int>(x / Lx * nx));
    const int iy = std::min(ny - 1, static_cast<int>(y / Ly * ny));
    const int iz = std::min(nz - 1, static_cast<int>(z / Lz * nz));
    const int c = (iz * ny + iy) * nx + ix;
    next[static_cast<std::size_t>(i)] = head[static_cast<std::size_t>(c)];
    head[static_cast<std::size_t>(c)] = i;
  }
  const double ox = wrapped[static_cast<std::size_t>(3 * atom)];
  const double oy = wrapped[static_cast<std::size_t>(3 * atom) + 1];
  const double oz = wrapped[static_cast<std::size_t>(3 * atom) + 2];
  const int ix0 = std::min(nx - 1, static_cast<int>(ox / Lx * nx));
  const int iy0 = std::min(ny - 1, static_cast<int>(oy / Ly * ny));
  const int iz0 = std::min(nz - 1, static_cast<int>(oz / Lz * nz));
  const int rx = std::max(1, static_cast<int>(std::ceil(cut / (Lx / nx))));
  const int ry = std::max(1, static_cast<int>(std::ceil(cut / (Ly / ny))));
  const int rz = std::max(1, static_cast<int>(std::ceil(cut / (Lz / nz))));
  struct Near {
    double d2;
    double x, y, z;
  };
  std::vector<Near> neigh;
  neigh.reserve(64);
  for (int dz = -rz; dz <= rz; dz++) {
    int iz = iz0 + dz;
    int sz = 0;
    while (iz < 0) {
      iz += nz;
      sz--;
    }
    while (iz >= nz) {
      iz -= nz;
      sz++;
    }
    for (int dy = -ry; dy <= ry; dy++) {
      int iy = iy0 + dy;
      int sy = 0;
      while (iy < 0) {
        iy += ny;
        sy--;
      }
      while (iy >= ny) {
        iy -= ny;
        sy++;
      }
      for (int dx = -rx; dx <= rx; dx++) {
        int ix = ix0 + dx;
        int sx = 0;
        while (ix < 0) {
          ix += nx;
          sx--;
        }
        while (ix >= nx) {
          ix -= nx;
          sx++;
        }
        const double shx = sx * Lx;
        const double shy = sy * Ly;
        const double shz = sz * Lz;
        for (int j = head[static_cast<std::size_t>((iz * ny + iy) * nx + ix)]; j >= 0;
             j = next[static_cast<std::size_t>(j)]) {
          if (j == atom && sx == 0 && sy == 0 && sz == 0) {
            continue;
          }
          const double x = wrapped[static_cast<std::size_t>(3 * j)] + shx - ox;
          const double y = wrapped[static_cast<std::size_t>(3 * j) + 1] + shy - oy;
          const double z = wrapped[static_cast<std::size_t>(3 * j) + 2] + shz - oz;
          const double d2 = x * x + y * y + z * z;
          if (d2 <= cut2 && d2 > 1e-16) {
            neigh.push_back({d2, cartesian(atom, 0) + x, cartesian(atom, 1) + y,
                             cartesian(atom, 2) + z});
          }
        }
      }
    }
  }
  if (static_cast<int>(neigh.size()) < 16) {
    return {};
  }
  std::nth_element(neigh.begin(), neigh.begin() + 16, neigh.end(),
                   [](const Near &a, const Near &b) { return a.d2 < b.d2; });
  std::sort(neigh.begin(), neigh.begin() + 16,
            [](const Near &a, const Near &b) { return a.d2 < b.d2; });
  if (neigh[15].d2 > cut2) {
    return {};
  }
  Eigen::MatrixXd shell(17, 3);
  shell(0, 0) = cartesian(atom, 0);
  shell(0, 1) = cartesian(atom, 1);
  shell(0, 2) = cartesian(atom, 2);
  for (int k = 0; k < 16; k++) {
    shell(k + 1, 0) = neigh[static_cast<std::size_t>(k)].x;
    shell(k + 1, 1) = neigh[static_cast<std::size_t>(k)].y;
    shell(k + 1, 2) = neigh[static_cast<std::size_t>(k)].z;
  }
  return shell;
}

Eigen::MatrixXd iceShell(const Eigen::MatrixXd &cartesian, const Eigen::Matrix3d &vectors,
                         int atom) {
  const int n = static_cast<int>(cartesian.rows());
  if (atom < 0 || atom >= n || cartesian.cols() != 3) {
    return {};
  }
  const double skew = std::abs(vectors(0, 1)) + std::abs(vectors(0, 2)) +
                      std::abs(vectors(1, 0)) + std::abs(vectors(1, 2)) +
                      std::abs(vectors(2, 0)) + std::abs(vectors(2, 1));
  const double Lx = vectors(0, 0);
  const double Ly = vectors(1, 1);
  const double Lz = vectors(2, 2);
  if (n >= 64 && skew < 1e-8 && Lx > 12.0 && Ly > 12.0 && Lz > 12.0) {
    return iceShellOrtho(cartesian, Lx, Ly, Lz, atom);
  }
  const Eigen::Vector3d origin = cartesian.row(atom).transpose();
  struct Neigh {
    double dist;
    Eigen::Vector3d at;
  };
  std::vector<Neigh> neigh;
  neigh.reserve(static_cast<size_t>(n) * 125);
  for (int j = 0; j < n; j++) {
    for (int s0 = -2; s0 <= 2; s0++) {
      for (int s1 = -2; s1 <= 2; s1++) {
        for (int s2 = -2; s2 <= 2; s2++) {
          if (j == atom && s0 == 0 && s1 == 0 && s2 == 0) {
            continue;
          }
          const Eigen::Vector3d shift = static_cast<double>(s0) * vectors.row(0).transpose() +
                                        static_cast<double>(s1) * vectors.row(1).transpose() +
                                        static_cast<double>(s2) * vectors.row(2).transpose();
          const Eigen::Vector3d at = cartesian.row(j).transpose() + shift;
          const double dist = (at - origin).norm();
          if (dist > 1e-8) {
            neigh.push_back({dist, at});
          }
        }
      }
    }
  }
  if (static_cast<int>(neigh.size()) < 16) {
    return {};
  }
  std::nth_element(neigh.begin(), neigh.begin() + 16, neigh.end(),
                   [](const Neigh &a, const Neigh &b) { return a.dist < b.dist; });
  std::sort(neigh.begin(), neigh.begin() + 16,
            [](const Neigh &a, const Neigh &b) { return a.dist < b.dist; });
  if (neigh[15].dist > 5.0) {
    return {};
  }
  Eigen::MatrixXd shell(17, 3);
  shell.row(0) = origin.transpose();
  for (int k = 0; k < 16; k++) {
    shell.row(k + 1) = neigh[static_cast<size_t>(k)].at.transpose();
  }
  return shell;
}

bool iraLinked() { return ira::available(); }

bool spglibLinked() {
#ifdef SEAMS_HAS_SPGLIB
  return true;
#else
  return false;
#endif
}

const char *matcher() { return ira::available() ? "ira" : "kabsch"; }

} // namespace symm

namespace {

void copySymbol(const std::string &text, char *symbol, int symbol_len) {
  if (symbol == nullptr || symbol_len <= 0) {
    return;
  }
  std::snprintf(symbol, static_cast<size_t>(symbol_len), "%s", text.c_str());
}

} // namespace

int seams_symmetry_match(const double *points, int n, int *template_id, double *rmsd,
                         double *hausdorff, double *rotation9, double *translation3,
                         int *permutation) {
  if (points == nullptr || template_id == nullptr || rmsd == nullptr || n < 1) {
    return 1;
  }
  Eigen::MatrixXd env(n, 3);
  for (int i = 0; i < n; i++) {
    env(i, 0) = points[3 * i];
    env(i, 1) = points[3 * i + 1];
    env(i, 2) = points[3 * i + 2];
  }
  const symm::Overlay hit = symm::best(env);
  if (hit.status != 0) {
    return 1;
  }
  *template_id = static_cast<int>(hit.kind);
  *rmsd = hit.rmsd;
  if (hausdorff != nullptr) {
    *hausdorff = hit.hausdorff;
  }
  if (rotation9 != nullptr) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        rotation9[3 * i + j] = hit.rotation(i, j);
      }
    }
  }
  if (translation3 != nullptr) {
    translation3[0] = hit.translation[0];
    translation3[1] = hit.translation[1];
    translation3[2] = hit.translation[2];
  }
  if (permutation != nullptr) {
    for (int i = 0; i < n && i < static_cast<int>(hit.permutation.size()); i++) {
      permutation[i] = hit.permutation[static_cast<size_t>(i)];
    }
  }
  return 0;
}

int seams_symmetry_point_group(const double *points, int n, double sym_thr, char *symbol,
                               int symbol_len, int *n_operations) {
  if (points == nullptr || n < 1) {
    return 1;
  }
  Eigen::MatrixXd env(n, 3);
  for (int i = 0; i < n; i++) {
    env(i, 0) = points[3 * i];
    env(i, 1) = points[3 * i + 1];
    env(i, 2) = points[3 * i + 2];
  }
  const symm::Group group = symm::pointGroup(env, sym_thr);
  if (group.status != 0) {
    return group.status;
  }
  copySymbol(group.symbol, symbol, symbol_len);
  if (n_operations != nullptr) {
    *n_operations = group.nOperations;
  }
  return 0;
}

int seams_symmetry_space_group(const double *lattice9, const double *cartesian, int n,
                               double symprec, int *number, char *symbol, int symbol_len) {
  if (lattice9 == nullptr || cartesian == nullptr || number == nullptr || n < 1) {
    return 1;
  }
  Eigen::Matrix3d vectors;
  for (int vec = 0; vec < 3; vec++) {
    for (int comp = 0; comp < 3; comp++) {
      vectors(vec, comp) = lattice9[3 * vec + comp];
    }
  }
  Eigen::MatrixXd xyz(n, 3);
  for (int i = 0; i < n; i++) {
    xyz(i, 0) = cartesian[3 * i];
    xyz(i, 1) = cartesian[3 * i + 1];
    xyz(i, 2) = cartesian[3 * i + 2];
  }
  const symm::Space space = symm::spaceGroup(vectors, xyz, symprec);
  if (space.status != 0) {
    return space.status;
  }
  *number = space.number;
  copySymbol(space.symbol, symbol, symbol_len);
  return 0;
}

int seams_symmetry_template(int id, double *points, int cap, int *n_out) {
  if (id < 0 || id >= symm::kTemplateCount || n_out == nullptr) {
    return 1;
  }
  const Eigen::MatrixXd m = symm::points(static_cast<symm::Kind>(id));
  *n_out = static_cast<int>(m.rows());
  if (points == nullptr || cap < *n_out) {
    return points == nullptr ? 0 : 1;
  }
  for (int i = 0; i < *n_out; i++) {
    points[3 * i] = m(i, 0);
    points[3 * i + 1] = m(i, 1);
    points[3 * i + 2] = m(i, 2);
  }
  return 0;
}

int seams_symmetry_crystal(int which, double *lattice9, double *xyz, int cap, int *n_out) {
  if (which < 0 || which > 2 || n_out == nullptr) {
    return 1;
  }
  const symm::Crystal kinds[3] = {symm::Crystal::iceIh, symm::Crystal::iceIc,
                                  symm::Crystal::fault};
  const symm::Lattice cell = symm::crystal(kinds[which]);
  *n_out = static_cast<int>(cell.cartesian.rows());
  if (lattice9 != nullptr) {
    for (int vec = 0; vec < 3; vec++) {
      for (int comp = 0; comp < 3; comp++) {
        lattice9[3 * vec + comp] = cell.vectors(vec, comp);
      }
    }
  }
  if (xyz == nullptr) {
    return 0;
  }
  if (cap < *n_out) {
    return 1;
  }
  for (int i = 0; i < *n_out; i++) {
    xyz[3 * i] = cell.cartesian(i, 0);
    xyz[3 * i + 1] = cell.cartesian(i, 1);
    xyz[3 * i + 2] = cell.cartesian(i, 2);
  }
  return 0;
}

int seams_symmetry_backend(char *name, int name_len) {
  copySymbol(symm::matcher(), name, name_len);
  return symm::iraLinked() ? 0 : 0;
}

struct SoaCols {
  size_t cap = 0;
  double *x = nullptr;
  double *y = nullptr;
  double *z = nullptr;

  ~SoaCols() {
    ::operator delete(x);
    ::operator delete(y);
    ::operator delete(z);
  }

  void fit(size_t n) {
    if (n <= cap) {
      return;
    }
    ::operator delete(x);
    ::operator delete(y);
    ::operator delete(z);
    x = static_cast<double *>(::operator new(n * sizeof(double)));
    y = static_cast<double *>(::operator new(n * sizeof(double)));
    z = static_cast<double *>(::operator new(n * sizeof(double)));
    cap = n;
  }
};

SoaCols &soaCols() {
  thread_local SoaCols cols;
  return cols;
}

int seams_symmetry_classify_tol(const double *xyz, int n, const double *lattice9,
                                double shell_rmsd, int *labels) {
  if (xyz == nullptr || lattice9 == nullptr || labels == nullptr || n < 1) {
    return 1;
  }
  double lat[3][3];
  for (int vec = 0; vec < 3; vec++) {
    for (int comp = 0; comp < 3; comp++) {
      lat[vec][comp] = lattice9[3 * vec + comp];
    }
  }
  SoaCols &cols = soaCols();
  cols.fit(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    cols.x[i] = xyz[3 * i];
    cols.y[i] = xyz[3 * i + 1];
    cols.z[i] = xyz[3 * i + 2];
  }
  symm::classifyIceArrays(cols.x, cols.y, cols.z, n, lat, labels, shell_rmsd);
  return 0;
}

int seams_symmetry_classify(const double *xyz, int n, const double *lattice9,
                            int *labels) {
  return seams_symmetry_classify_tol(xyz, n, lattice9, symm::kShellRmsd, labels);
}
