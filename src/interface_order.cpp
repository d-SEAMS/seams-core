//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <interface_order.hpp>

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <local_symmetry.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_one_dim.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <utility>

namespace iface {
namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

constexpr int kLabels = 9;
constexpr double kCnaScale = 0.5 * (1.0 + 1.41421356237);

int labelOf(const Cloud &cloud, int i) {
  return static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType);
}

bool crystalLabel(int label) {
  return label == static_cast<int>(molSys::atom_state_type::cubic) ||
         label == static_cast<int>(molSys::atom_state_type::hexagonal);
}

bool faceLabel(int label) {
  return label == static_cast<int>(molSys::atom_state_type::interfacial);
}

Cloud classified(const Cloud &cloud, const std::vector<std::vector<int>> &rows,
                 int coordination, bool closed, bool neighbourTest) {
  Cloud out = cloud;
  chill::getCorrelPlus(out, rows, false, coordination, closed);
  chill::getIceTypePlusNoPrint(out, rows, false, neighbourTest);
  return out;
}

Cause causeOf(int published, int noTest, int shellClosed, int ovito) {
  if (published == ovito) {
    return Cause::agree;
  }
  const bool neighbour = published != noTest;
  const bool shell = noTest != shellClosed;
  const bool edge = shellClosed != ovito;
  const int flips = static_cast<int>(neighbour) + static_cast<int>(shell) +
                    static_cast<int>(edge);
  if (flips == 1 && neighbour) {
    return Cause::neighbour;
  }
  if (flips == 1 && shell) {
    return Cause::q3Shell;
  }
  if (flips == 1 && edge) {
    return Cause::openEdge;
  }
  return Cause::mixed;
}

double axisCoord(const molSys::Point<double> &p, int axis) {
  if (axis == 0) {
    return p.x;
  }
  if (axis == 1) {
    return p.y;
  }
  return p.z;
}

double lengthOf(const Cloud &cloud, int axis) {
  return cloud.box.size() > static_cast<std::size_t>(axis) ? cloud.box[static_cast<std::size_t>(axis)]
                                                           : 0.0;
}

double originOf(const Cloud &cloud, int axis) {
  return cloud.boxLow.size() > static_cast<std::size_t>(axis)
             ? cloud.boxLow[static_cast<std::size_t>(axis)]
             : 0.0;
}

struct SheetSplit {
  int thin = 0;
  int thick = 0;
};

SheetSplit splitSheets(const Cloud &cloud, int axis) {
  SheetSplit split;
  if (cloud.nop <= 0) {
    return split;
  }
  std::vector<double> q(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    q[static_cast<std::size_t>(i)] = axisCoord(cloud.pts[static_cast<std::size_t>(i)], axis);
  }
  std::sort(q.begin(), q.end());
  int i = 0;
  while (i < cloud.nop) {
    int j = i;
    while (j + 1 < cloud.nop &&
           q[static_cast<std::size_t>(j + 1)] - q[static_cast<std::size_t>(j)] <= 1.0) {
      ++j;
    }
    if (q[static_cast<std::size_t>(j)] - q[static_cast<std::size_t>(i)] < 1.0) {
      ++split.thin;
    } else {
      ++split.thick;
    }
    i = j + 1;
  }
  return split;
}

struct FrontFit {
  Profile profile;
};

FrontFit fitPhi(const std::vector<double> &phi, const std::vector<double> &mid, double span,
                bool radial) {
  FrontFit fit;
  const int n = static_cast<int>(phi.size());
  if (n == 0 || !(span > 0.0)) {
    return fit;
  }
  auto crossing = [](double s0, double p0, double s1, double p1, double level) {
    const double t = (level - p0) / (p1 - p0);
    return s0 + t * (s1 - s0);
  };
  if (radial) {
    int imax = 0;
    for (int i = 1; i < n; ++i) {
      if (phi[static_cast<std::size_t>(i)] > phi[static_cast<std::size_t>(imax)]) {
        imax = i;
      }
    }
    double prevS = mid[static_cast<std::size_t>(imax)];
    double prevP = phi[static_cast<std::size_t>(imax)];
    bool have90 = false;
    bool have50 = false;
    bool have10 = false;
    double s90 = 0;
    double s50 = 0;
    double s10 = 0;
    if (prevP < 0.5) {
      return fit;
    }
    for (int i = imax + 1; i < n; ++i) {
      const double s = mid[static_cast<std::size_t>(i)];
      const double p = phi[static_cast<std::size_t>(i)];
      if (!have90 && prevP >= 0.9 && p < 0.9 && p != prevP) {
        s90 = crossing(prevS, prevP, s, p, 0.9);
        have90 = true;
      }
      if (!have50 && prevP >= 0.5 && p < 0.5 && p != prevP) {
        s50 = crossing(prevS, prevP, s, p, 0.5);
        have50 = true;
      }
      if (!have10 && prevP >= 0.1 && p < 0.1 && p != prevP) {
        s10 = crossing(prevS, prevP, s, p, 0.1);
        have10 = true;
      }
      prevS = s;
      prevP = p;
    }
    if (!have50) {
      return fit;
    }
    Front front;
    front.position = s50;
    front.defined = true;
    if (have90 && have10) {
      front.width = std::abs(s10 - s90);
    }
    fit.profile.fronts.push_back(front);
    fit.profile.thickness = 2.0 * s50;
    fit.profile.defined = true;
    return fit;
  }

  int imax = 0;
  for (int i = 1; i < n; ++i) {
    if (phi[static_cast<std::size_t>(i)] > phi[static_cast<std::size_t>(imax)]) {
      imax = i;
    }
  }
  if (phi[static_cast<std::size_t>(imax)] < 0.5) {
    return fit;
  }
  double left50 = 0;
  double right50 = 0;
  bool haveL = false;
  bool haveR = false;
  for (int dir = -1; dir <= 1; dir += 2) {
    double prevS = 0.0;
    double prevP = phi[static_cast<std::size_t>(imax)];
    bool have90 = false;
    bool have50 = false;
    bool have10 = false;
    double s90 = 0;
    double s50 = 0;
    double s10 = 0;
    for (int step = 1; step < n; ++step) {
      const int b = (imax + dir * step + n * 8) % n;
      const double s = step * (span / static_cast<double>(n));
      const double p = phi[static_cast<std::size_t>(b)];
      if (!have90 && prevP >= 0.9 && p < 0.9 && p != prevP) {
        s90 = crossing(prevS, prevP, s, p, 0.9);
        have90 = true;
      }
      if (!have50 && prevP >= 0.5 && p < 0.5 && p != prevP) {
        s50 = crossing(prevS, prevP, s, p, 0.5);
        have50 = true;
      }
      if (!have10 && prevP >= 0.1 && p < 0.1 && p != prevP) {
        s10 = crossing(prevS, prevP, s, p, 0.1);
        have10 = true;
      }
      prevS = s;
      prevP = p;
      if (have50 && have10) {
        break;
      }
    }
    if (!have50) {
      continue;
    }
    Front front;
    const double z0 = mid[static_cast<std::size_t>(imax)];
    double z = z0 + static_cast<double>(dir) * s50;
    z = std::fmod(z, span);
    if (z < 0.0) {
      z += span;
    }
    front.position = z;
    front.defined = true;
    if (have90 && have10) {
      front.width = std::abs(s10 - s90);
    }
    fit.profile.fronts.push_back(front);
    if (dir < 0) {
      left50 = s50;
      haveL = true;
    } else {
      right50 = s50;
      haveR = true;
    }
  }
  if (haveL && haveR) {
    fit.profile.thickness = left50 + right50;
    fit.profile.defined = true;
  }
  return fit;
}

std::vector<double> phiOf(const std::vector<int> &numer, const std::vector<int> &den) {
  std::vector<double> phi(numer.size(), 0.0);
  double last = 0.0;
  for (std::size_t i = 0; i < numer.size(); ++i) {
    if (den[i] > 0) {
      last = static_cast<double>(numer[i]) / static_cast<double>(den[i]);
    }
    phi[i] = last;
  }
  return phi;
}

void fillFractions(Profile &profile, int n, int crystal, int face) {
  if (n <= 0) {
    return;
  }
  profile.crystal = static_cast<double>(crystal) / static_cast<double>(n);
  profile.interfacial = static_cast<double>(face) / static_cast<double>(n);
  profile.solid = profile.crystal + profile.interfacial;
}

int cnaSignature(const double *sx, const double *sy, const double *sz, int nn, double cut2,
                 int j) {
  int common[14];
  int nc = 0;
  for (int k = 0; k < nn; ++k) {
    if (k == j) {
      continue;
    }
    const double dx = sx[j] - sx[k];
    const double dy = sy[j] - sy[k];
    const double dz = sz[j] - sz[k];
    if (dx * dx + dy * dy + dz * dz < cut2 && nc < 14) {
      common[nc++] = k;
    }
  }
  int nBonds = 0;
  int deg[14] = {};
  for (int a = 0; a < nc; ++a) {
    for (int b = a + 1; b < nc; ++b) {
      const int ia = common[a];
      const int ib = common[b];
      const double dx = sx[ia] - sx[ib];
      const double dy = sy[ia] - sy[ib];
      const double dz = sz[ia] - sz[ib];
      if (dx * dx + dy * dy + dz * dz < cut2) {
        ++nBonds;
        ++deg[a];
        ++deg[b];
      }
    }
  }
  int longest = 0;
  for (int a = 0; a < nc; ++a) {
    longest = std::max(longest, deg[a]);
  }
  // Stukowski's third index is the longest bonded chain. A disjoint pair of
  // bonds has degree 1. A path of two bonds has degree 2. A 5-cycle has
  // degree 2 as well, so the icosahedral 555 row is recognised by the bond
  // count together with every common neighbour having two bonds.
  if (nc == 4 && nBonds == 2 && longest == 1) {
    return 421;
  }
  if (nc == 4 && nBonds == 2 && longest >= 2) {
    return 422;
  }
  if (nc == 5 && nBonds == 5) {
    return 555;
  }
  if (nc == 4 && nBonds == 4) {
    return 444;
  }
  if (nc == 6 && nBonds == 6) {
    return 666;
  }
  return 0;
}

int adaptiveCnaLabel(const Cloud &cloud, int iatom, const std::vector<int> &offsets,
                     const std::vector<int> &cols) {
  struct Near {
    double d2;
    double x, y, z;
  };
  Near near[16];
  int nb = 0;
  const int a = offsets[static_cast<std::size_t>(iatom)];
  const int b = offsets[static_cast<std::size_t>(iatom + 1)];
  const gen::FracBox frame = gen::makeFracBox(cloud);
  const auto &pi = cloud.pts[static_cast<std::size_t>(iatom)];
  for (int k = a; k < b; ++k) {
    const int j = cols[static_cast<std::size_t>(k)];
    if (j < 0 || j >= cloud.nop || j == iatom) {
      continue;
    }
    const auto &pj = cloud.pts[static_cast<std::size_t>(j)];
    const auto dr = gen::fracDelta(frame, pi.x, pi.y, pi.z, pj.x, pj.y, pj.z);
    const double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    if (nb < 16) {
      near[nb++] = {d2, dr[0], dr[1], dr[2]};
    } else {
      int worst = 0;
      for (int t = 1; t < nb; ++t) {
        if (near[t].d2 > near[worst].d2) {
          worst = t;
        }
      }
      if (d2 < near[worst].d2) {
        near[worst] = {d2, dr[0], dr[1], dr[2]};
      }
    }
  }
  if (nb < 12) {
    return 0;
  }
  std::sort(near, near + nb, [](const Near &p, const Near &q) { return p.d2 < q.d2; });
  auto score = [&](int nn) {
    double mean = 0.0;
    double sx[14], sy[14], sz[14];
    for (int k = 0; k < nn; ++k) {
      mean += std::sqrt(near[k].d2);
      sx[k] = near[k].x;
      sy[k] = near[k].y;
      sz[k] = near[k].z;
    }
    mean /= static_cast<double>(nn);
    const double cut = kCnaScale * mean;
    const double cut2 = cut * cut;
    int n421 = 0;
    int n422 = 0;
    int n555 = 0;
    int n444 = 0;
    int n666 = 0;
    for (int j = 0; j < nn; ++j) {
      const int sig = cnaSignature(sx, sy, sz, nn, cut2, j);
      n421 += sig == 421;
      n422 += sig == 422;
      n555 += sig == 555;
      n444 += sig == 444;
      n666 += sig == 666;
    }
    if (nn == 12 && n421 == 12) {
      return 1;
    }
    if (nn == 12 && n421 == 6 && n422 == 6) {
      return 2;
    }
    if (nn == 12 && n555 == 12) {
      return 3;
    }
    if (nn == 14 && n444 == 6 && n666 == 8) {
      return 4;
    }
    return 0;
  };
  const int twelve = score(12);
  if (twelve != 0) {
    return twelve;
  }
  if (nb >= 14) {
    return score(14);
  }
  return 0;
}

double tetraRmsd(const double u[4][3]) {
  const double inv = 1.0 / std::sqrt(3.0);
  const double tet[4][3] = {{inv, inv, inv},
                            {inv, -inv, -inv},
                            {-inv, inv, -inv},
                            {-inv, -inv, inv}};
  int order[4] = {0, 1, 2, 3};
  double best = 1.0e9;
  do {
    Eigen::Matrix<double, 3, 4> A;
    Eigen::Matrix<double, 3, 4> B;
    for (int k = 0; k < 4; ++k) {
      A(0, k) = u[k][0];
      A(1, k) = u[k][1];
      A(2, k) = u[k][2];
      const int t = order[k];
      B(0, k) = tet[t][0];
      B(1, k) = tet[t][1];
      B(2, k) = tet[t][2];
    }
    const Eigen::Vector3d ca = A.rowwise().mean();
    const Eigen::Vector3d cb = B.rowwise().mean();
    const Eigen::Matrix<double, 3, 4> Ac = A.colwise() - ca;
    const Eigen::Matrix<double, 3, 4> Bc = B.colwise() - cb;
    const Eigen::Matrix3d H = Ac * Bc.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d V = svd.matrixV();
    const Eigen::Matrix3d U = svd.matrixU();
    if ((V * U.transpose()).determinant() < 0.0) {
      V.col(2) *= -1.0;
    }
    const Eigen::Matrix3d R = V * U.transpose();
    const Eigen::Matrix<double, 3, 4> D = R * Ac - Bc;
    best = std::min(best, std::sqrt(D.squaredNorm() / 4.0));
  } while (std::next_permutation(order, order + 4));
  return best;
}

bool ptmPass(const Cloud &cloud, int iatom, const std::vector<std::vector<int>> &rows) {
  if (iatom < 0 || static_cast<std::size_t>(iatom) >= rows.size()) {
    return false;
  }
  struct Near {
    double d2;
    int j;
  };
  std::vector<Near> near;
  const auto &row = rows[static_cast<std::size_t>(iatom)];
  for (std::size_t k = 1; k < row.size(); ++k) {
    const auto found = cloud.idIndexMap.find(row[k]);
    if (found == cloud.idIndexMap.end()) {
      continue;
    }
    const auto dr = gen::relDist(cloud, iatom, found->second);
    const double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    near.push_back({d2, found->second});
  }
  if (static_cast<int>(near.size()) < 4) {
    return false;
  }
  std::sort(near.begin(), near.end(),
            [](const Near &p, const Near &q) { return p.d2 < q.d2; });
  double u[4][3];
  for (int k = 0; k < 4; ++k) {
    const auto dr = gen::relDist(cloud, iatom, near[static_cast<std::size_t>(k)].j);
    const double nrm = std::sqrt(near[static_cast<std::size_t>(k)].d2);
    if (!(nrm > 1.0e-8)) {
      return false;
    }
    u[k][0] = -dr[0] / nrm;
    u[k][1] = -dr[1] / nrm;
    u[k][2] = -dr[2] / nrm;
  }
  // OVITO Polyhedral Template Matching accepts a diamond environment when
  // the RMSD of the neighbour star is below 0.1.
  return tetraRmsd(u) < 0.1;
}

int countFrames(const std::string &path) {
  std::ifstream in(path);
  int n = 0;
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("ITEM: TIMESTEP", 0) == 0) {
      ++n;
    }
  }
  return n;
}

void addPoint(Cloud &cloud, double x, double y, double z) {
  molSys::Point<double> p;
  p.x = x;
  p.y = y;
  p.z = z;
  p.type = 1;
  p.atomID = cloud.nop + 1;
  p.molID = p.atomID;
  cloud.idIndexMap[p.atomID] = cloud.nop;
  cloud.pts.push_back(p);
  cloud.nop += 1;
}

double wrap(double q, double L) {
  if (!(L > 0.0)) {
    return q;
  }
  q = std::fmod(q, L);
  if (q < 0.0) {
    q += L;
  }
  return q;
}

} // namespace

Study studyCloud(Cloud cloud, double cutoff, double binWidth, bool census) {
  Study study;
  study.n = cloud.nop;
  if (cloud.nop <= 0 || cloud.box.size() < 3) {
    return study;
  }
  const auto rows = nneigh::neighListO(cutoff, cloud, 1);
  const Cloud published = classified(cloud, rows, 4, true, true);
  const Cloud noTest = classified(cloud, rows, 4, true, false);
  const Cloud shellClosed = classified(cloud, rows, 0, true, false);
  const Cloud ovito = classified(cloud, rows, 0, false, false);

  std::vector<char> disagree(static_cast<std::size_t>(cloud.nop), 0);
  std::vector<char> isCrystal(static_cast<std::size_t>(cloud.nop), 0);
  int crystalPubN = 0;
  int facePubN = 0;
  int crystalOviN = 0;
  int faceOviN = 0;
  for (int i = 0; i < cloud.nop; ++i) {
    const int lp = labelOf(published, i);
    const int lo = labelOf(ovito, i);
    const Cause cause = causeOf(lp, labelOf(noTest, i), labelOf(shellClosed, i), lo);
    ++study.causeCount[static_cast<int>(cause)];
    if (lp >= 0 && lp < kLabels) {
      ++study.countPub[lp];
    }
    if (lo >= 0 && lo < kLabels) {
      ++study.countOvi[lo];
    }
    if (lp >= 0 && lp < kLabels && lo >= 0 && lo < kLabels) {
      ++study.transition[lp][lo];
    }
    const int shell = std::min(chill::cutoffShellCount(cloud, rows, i), 7);
    if (cause == Cause::agree) {
      ++study.shellAgree[shell];
    } else {
      ++study.nDisagree;
      ++study.shellDisagree[shell];
      disagree[static_cast<std::size_t>(i)] = 1;
    }
    if (crystalLabel(lp)) {
      ++crystalPubN;
      isCrystal[static_cast<std::size_t>(i)] = 1;
    }
    facePubN += faceLabel(lp) ? 1 : 0;
    crystalOviN += crystalLabel(lo) ? 1 : 0;
    faceOviN += faceLabel(lo) ? 1 : 0;
  }

  if (census) {
  const auto indexRows = nneigh::neighbourListByIndex(cloud, rows);
  const auto rings = primitive::ringNetwork(indexRows, 6);
  std::vector<char> ringHit[3] = {std::vector<char>(static_cast<std::size_t>(cloud.nop), 0),
                                  std::vector<char>(static_cast<std::size_t>(cloud.nop), 0),
                                  std::vector<char>(static_cast<std::size_t>(cloud.nop), 0)};
  std::vector<char> prismAtom(static_cast<std::size_t>(cloud.nop), 0);
  std::vector<double> rmsd(static_cast<std::size_t>(cloud.nop), 0.0);
  for (int size = 4; size <= 6; ++size) {
    const auto one = ring::getSingleRingSize(rings, size);
    if (one.empty()) {
      continue;
    }
    if (one.size() > 400) {
      continue;
    }
    std::vector<ring::strucType> types(one.size(), ring::strucType::unclassified);
    int perfect = 0;
    int imperfect = 0;
    const auto listed =
        ring::findPrisms(one, types, perfect, imperfect, indexRows, cloud, rmsd, false);
    study.prismBlocks += perfect;
    for (int ringIndex : listed) {
      if (ringIndex < 0 || static_cast<std::size_t>(ringIndex) >= one.size()) {
        continue;
      }
      for (int atom : one[static_cast<std::size_t>(ringIndex)]) {
        if (atom >= 0 && atom < cloud.nop) {
          prismAtom[static_cast<std::size_t>(atom)] = 1;
        }
      }
    }
  }
  for (const auto &ring : rings) {
    const int size = static_cast<int>(ring.size());
    if (size < 4 || size > 6) {
      continue;
    }
    for (int atom : ring) {
      if (atom >= 0 && atom < cloud.nop) {
        ringHit[size - 4][static_cast<std::size_t>(atom)] = 1;
      }
    }
  }
  for (int i = 0; i < cloud.nop; ++i) {
    const bool bad = disagree[static_cast<std::size_t>(i)] != 0;
    for (int s = 0; s < 3; ++s) {
      if (ringHit[s][static_cast<std::size_t>(i)] == 0) {
        continue;
      }
      if (bad) {
        ++study.ringDisagree[s];
      } else {
        ++study.ringAgree[s];
      }
    }
    study.prismAtoms += prismAtom[static_cast<std::size_t>(i)] ? 1 : 0;
  }

  int sheetAxis = 2;
  int sheetThin = 0;
  int mostThin = 0;
  for (int axis = 0; axis < 3; ++axis) {
    const SheetSplit split = splitSheets(cloud, axis);
    if (split.thin > mostThin) {
      mostThin = split.thin;
    }
    // One or two thin slabs and no thick slab is a monolayer or a bilayer.
    // A nanotube has a thin ring at every step along its axis, so it is
    // not either of those.
    if (split.thick == 0 && (split.thin == 1 || split.thin == 2) && split.thin > sheetThin) {
      sheetThin = split.thin;
      sheetAxis = axis;
    }
  }
  study.thinSheets = sheetThin > 0 ? sheetThin : mostThin;
  std::vector<char> sheetAtom(static_cast<std::size_t>(cloud.nop), 0);
  if (sheetThin == 1 || sheetThin == 2) {
    for (const auto &ring : rings) {
      if (ring.size() != 6) {
        continue;
      }
      bool prismRing = false;
      int members = 0;
      for (int atom : ring) {
        if (atom < 0 || atom >= cloud.nop) {
          continue;
        }
        prismRing = prismRing || prismAtom[static_cast<std::size_t>(atom)] != 0;
        ++members;
      }
      if (prismRing || members < 6) {
        continue;
      }
      const auto &p0 = cloud.pts[static_cast<std::size_t>(ring[0])];
      const auto &p1 = cloud.pts[static_cast<std::size_t>(ring[1])];
      const auto &p2 = cloud.pts[static_cast<std::size_t>(ring[2])];
      const double e1[3] = {p1.x - p0.x, p1.y - p0.y, p1.z - p0.z};
      const double e2[3] = {p2.x - p0.x, p2.y - p0.y, p2.z - p0.z};
      const double nrm[3] = {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                             e1[0] * e2[1] - e1[1] * e2[0]};
      const double nn = std::sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
      if (!(nn > 1.0e-8)) {
        continue;
      }
      if (std::abs(nrm[sheetAxis] / nn) < 0.85) {
        continue;
      }
      double meanQ = 0.0;
      for (int atom : ring) {
        if (atom < 0 || atom >= cloud.nop) {
          continue;
        }
        meanQ += axisCoord(cloud.pts[static_cast<std::size_t>(atom)], sheetAxis);
      }
      meanQ /= static_cast<double>(members);
      double spread = 0.0;
      for (int atom : ring) {
        if (atom < 0 || atom >= cloud.nop) {
          continue;
        }
        spread = std::max(
            spread,
            std::abs(axisCoord(cloud.pts[static_cast<std::size_t>(atom)], sheetAxis) - meanQ));
      }
      if (spread > 0.6) {
        continue;
      }
      for (int atom : ring) {
        sheetAtom[static_cast<std::size_t>(atom)] = 1;
      }
    }
  }
  for (int i = 0; i < cloud.nop; ++i) {
    if (sheetAtom[static_cast<std::size_t>(i)] == 0) {
      continue;
    }
    if (sheetThin == 2 && prismAtom[static_cast<std::size_t>(i)] == 0) {
      ++study.bilayerAtoms;
    } else if (sheetThin == 1) {
      ++study.polygonAtoms;
    }
  }

  std::vector<int> offsets;
  std::vector<int> cols;
  nneigh::orthoFullCsr(cloud, nullptr, cutoff, offsets, cols, false);
  std::vector<int> diamond;
  symm::classifyIceCloud(cloud, diamond);
  for (int kind : diamond) {
    if (kind == 1) {
      ++study.diamondCubic;
    } else if (kind == 0) {
      ++study.diamondHex;
    } else {
      ++study.diamondOther;
    }
  }

  std::vector<int> wideOff;
  std::vector<int> wideCol;
  nneigh::orthoFullCsr(cloud, nullptr, 8.0, wideOff, wideCol, false);
  const auto q6 = chill::steinhardtQl(cloud, rows, 6);
  double qDisagree = 0.0;
  double qAgree = 0.0;
  int nDisagreeQ = 0;
  int nAgreeQ = 0;
  for (int i = 0; i < cloud.nop; ++i) {
    const int cna = adaptiveCnaLabel(cloud, i, wideOff, wideCol);
    if (cna == 0) {
      ++study.cnaOther;
    } else {
      ++study.cnaCrystal;
    }
    const bool pass = ptmPass(cloud, i, rows);
    study.ptmPass += pass ? 1 : 0;
    if (disagree[static_cast<std::size_t>(i)] != 0) {
      study.ptmPassDisagree += pass ? 1 : 0;
      if (i < static_cast<int>(q6.ql.size()) && std::isfinite(q6.ql[static_cast<std::size_t>(i)])) {
        qDisagree += q6.ql[static_cast<std::size_t>(i)];
        ++nDisagreeQ;
      }
    } else if (i < static_cast<int>(q6.ql.size()) &&
               std::isfinite(q6.ql[static_cast<std::size_t>(i)])) {
      qAgree += q6.ql[static_cast<std::size_t>(i)];
      ++nAgreeQ;
    }
  }
  if (nDisagreeQ > 0) {
    study.q6Disagree = qDisagree / static_cast<double>(nDisagreeQ);
  }
  if (nAgreeQ > 0) {
    study.q6Agree = qAgree / static_cast<double>(nAgreeQ);
  }
  }

  const double lx = lengthOf(cloud, 0);
  const double ly = lengthOf(cloud, 1);
  const double lz = lengthOf(cloud, 2);
  int longAxis = 0;
  if (ly > lengthOf(cloud, longAxis)) {
    longAxis = 1;
  }
  if (lz > lengthOf(cloud, longAxis)) {
    longAxis = 2;
  }
  const double longL = lengthOf(cloud, longAxis);
  const double midL = longAxis == 0 ? std::max(ly, lz) : longAxis == 1 ? std::max(lx, lz)
                                                                       : std::max(lx, ly);
  study.kind = (longL > 1.3 * midL) ? Kind::slab : Kind::nucleus;
  study.axis = study.kind == Kind::slab ? longAxis : -1;

  const double dr = binWidth > 0.0 ? binWidth : 1.0;
  std::vector<int> den;
  std::vector<int> cPub;
  std::vector<int> sPub;
  std::vector<int> cOvi;
  std::vector<int> sOvi;
  std::vector<int> dcount;
  std::vector<double> mid;
  auto pushAtom = [&](int bin, int i, int lp, int lo) {
    if (bin < 0) {
      return;
    }
    if (static_cast<int>(den.size()) <= bin) {
      den.resize(static_cast<std::size_t>(bin + 1), 0);
      cPub.resize(den.size(), 0);
      sPub.resize(den.size(), 0);
      cOvi.resize(den.size(), 0);
      sOvi.resize(den.size(), 0);
      dcount.resize(den.size(), 0);
    }
    ++den[static_cast<std::size_t>(bin)];
    cPub[static_cast<std::size_t>(bin)] += crystalLabel(lp) ? 1 : 0;
    sPub[static_cast<std::size_t>(bin)] += (crystalLabel(lp) || faceLabel(lp)) ? 1 : 0;
    cOvi[static_cast<std::size_t>(bin)] += crystalLabel(lo) ? 1 : 0;
    sOvi[static_cast<std::size_t>(bin)] += (crystalLabel(lo) || faceLabel(lo)) ? 1 : 0;
    dcount[static_cast<std::size_t>(bin)] += disagree[static_cast<std::size_t>(i)] ? 1 : 0;
  };

  std::array<double, 3> com{0.0, 0.0, 0.0};
  int anchor = -1;
  for (int i = 0; i < cloud.nop; ++i) {
    if (isCrystal[static_cast<std::size_t>(i)] != 0) {
      anchor = i;
      break;
    }
  }
  if (study.kind == Kind::nucleus && anchor >= 0) {
    std::array<double, 3> acc{0.0, 0.0, 0.0};
    int nC = 0;
    for (int i = 0; i < cloud.nop; ++i) {
      if (isCrystal[static_cast<std::size_t>(i)] == 0) {
        continue;
      }
      const auto drv = gen::relDistFromPoint(cloud, anchor, cloud.pts[static_cast<std::size_t>(i)].x,
                                             cloud.pts[static_cast<std::size_t>(i)].y,
                                             cloud.pts[static_cast<std::size_t>(i)].z);
      acc[0] += drv[0];
      acc[1] += drv[1];
      acc[2] += drv[2];
      ++nC;
    }
    if (nC > 0) {
      com[0] = cloud.pts[static_cast<std::size_t>(anchor)].x - acc[0] / static_cast<double>(nC);
      com[1] = cloud.pts[static_cast<std::size_t>(anchor)].y - acc[1] / static_cast<double>(nC);
      com[2] = cloud.pts[static_cast<std::size_t>(anchor)].z - acc[2] / static_cast<double>(nC);
    }
  }

  if (study.kind == Kind::slab) {
    const int nb = std::max(8, static_cast<int>(std::lround(longL / dr)));
    const double width = longL / static_cast<double>(nb);
    den.assign(static_cast<std::size_t>(nb), 0);
    cPub.assign(den.size(), 0);
    sPub.assign(den.size(), 0);
    cOvi.assign(den.size(), 0);
    sOvi.assign(den.size(), 0);
    dcount.assign(den.size(), 0);
    mid.resize(static_cast<std::size_t>(nb));
    for (int b = 0; b < nb; ++b) {
      mid[static_cast<std::size_t>(b)] = (static_cast<double>(b) + 0.5) * width;
    }
    const double origin = originOf(cloud, longAxis);
    for (int i = 0; i < cloud.nop; ++i) {
      double q = axisCoord(cloud.pts[static_cast<std::size_t>(i)], longAxis) - origin;
      q = std::fmod(q, longL);
      if (q < 0.0) {
        q += longL;
      }
      int bin = static_cast<int>(q / width);
      if (bin >= nb) {
        bin = nb - 1;
      }
      pushAtom(bin, i, labelOf(published, i), labelOf(ovito, i));
    }
  } else {
    const double rMax = 0.5 * std::min(lx, std::min(ly, lz));
    const int nb = std::max(4, static_cast<int>(std::ceil(rMax / dr)));
    den.assign(static_cast<std::size_t>(nb), 0);
    cPub.assign(den.size(), 0);
    sPub.assign(den.size(), 0);
    cOvi.assign(den.size(), 0);
    sOvi.assign(den.size(), 0);
    dcount.assign(den.size(), 0);
    mid.resize(static_cast<std::size_t>(nb));
    for (int b = 0; b < nb; ++b) {
      mid[static_cast<std::size_t>(b)] = (static_cast<double>(b) + 0.5) * dr;
    }
    for (int i = 0; i < cloud.nop; ++i) {
      const double r = gen::unWrappedDistFromPoint(
          cloud, i, std::vector<double>{com[0], com[1], com[2]});
      int bin = static_cast<int>(r / dr);
      if (bin >= nb) {
        bin = nb - 1;
      }
      pushAtom(bin, i, labelOf(published, i), labelOf(ovito, i));
    }
  }

  study.bins.resize(den.size());
  for (std::size_t b = 0; b < den.size(); ++b) {
    study.bins[b].mid = mid[b];
    study.bins[b].n = den[b];
    study.bins[b].disagree = dcount[b];
    if (den[b] > 0) {
      const double inv = 1.0 / static_cast<double>(den[b]);
      study.bins[b].crystalPub = cPub[b] * inv;
      study.bins[b].solidPub = sPub[b] * inv;
      study.bins[b].crystalOvi = cOvi[b] * inv;
      study.bins[b].solidOvi = sOvi[b] * inv;
    }
  }
  const double span = study.kind == Kind::slab ? longL : (mid.empty() ? 0.0 : mid.back() + dr);
  const bool radial = study.kind == Kind::nucleus;
  auto applyFit = [&](Profile &profile, const std::vector<int> &numer, int crystalN, int faceN) {
    fillFractions(profile, cloud.nop, crystalN, faceN);
    const auto phi = phiOf(numer, den);
    auto fitted = fitPhi(phi, mid, span, radial);
    profile.fronts = std::move(fitted.profile.fronts);
    profile.thickness = fitted.profile.thickness;
    profile.defined = fitted.profile.defined;
  };
  applyFit(study.solidPub, sPub, crystalPubN, facePubN);
  applyFit(study.crystalPub, cPub, crystalPubN, facePubN);
  applyFit(study.solidOvi, sOvi, crystalOviN, faceOviN);
  applyFit(study.crystalOvi, cOvi, crystalOviN, faceOviN);

  auto depthTo = [&](const Profile &profile, int i) {
    if (study.kind == Kind::nucleus) {
      if (profile.fronts.empty() || !profile.fronts[0].defined) {
        return -1.0;
      }
      const double r = gen::unWrappedDistFromPoint(
          cloud, i, std::vector<double>{com[0], com[1], com[2]});
      return std::abs(r - profile.fronts[0].position);
    }
    double best = -1.0;
    const double origin = originOf(cloud, longAxis);
    double q = axisCoord(cloud.pts[static_cast<std::size_t>(i)], longAxis) - origin;
    q = std::fmod(q, longL);
    if (q < 0.0) {
      q += longL;
    }
    for (const auto &front : profile.fronts) {
      if (!front.defined) {
        continue;
      }
      double d = std::abs(q - front.position);
      d = std::min(d, longL - d);
      if (best < 0.0 || d < best) {
        best = d;
      }
    }
    return best;
  };
  double sumBad = 0.0;
  double sumGood = 0.0;
  for (int i = 0; i < cloud.nop; ++i) {
    const double d = depthTo(study.solidPub, i);
    if (!(d >= 0.0)) {
      continue;
    }
    if (disagree[static_cast<std::size_t>(i)] != 0) {
      sumBad += d;
      ++study.depthDisagreeN;
    } else {
      sumGood += d;
      ++study.depthAgreeN;
    }
  }
  if (study.depthDisagreeN > 0) {
    study.depthDisagree = sumBad / static_cast<double>(study.depthDisagreeN);
  }
  if (study.depthAgreeN > 0) {
    study.depthAgree = sumGood / static_cast<double>(study.depthAgreeN);
  }
  return study;
}

Study studyDump(const std::string &path, int frame, double cutoff, double binWidth, bool census) {
  Cloud cloud;
  cloud = sinp::readLammpsTrjO(path, frame, cloud, 1);
  return studyCloud(std::move(cloud), cutoff, binWidth, census);
}

namespace {

const char *labelName(int label) {
  switch (label) {
  case 0:
    return "cubic";
  case 1:
    return "hexagonal";
  case 2:
    return "water";
  case 3:
    return "interfacial";
  case 4:
    return "clathrate";
  case 5:
    return "interClathrate";
  case 6:
    return "unclassified";
  case 7:
    return "reCubic";
  case 8:
    return "reHex";
  default:
    return "unknown";
  }
}

void writeProfile(std::ostream &os, const char *name, const Profile &profile) {
  os << name << "_crystal " << profile.crystal << "\n";
  os << name << "_interfacial " << profile.interfacial << "\n";
  os << name << "_solid " << profile.solid << "\n";
  os << name << "_defined " << (profile.defined ? 1 : 0) << "\n";
  os << name << "_thickness " << profile.thickness << "\n";
  os << name << "_nfront " << profile.fronts.size() << "\n";
  for (std::size_t i = 0; i < profile.fronts.size(); ++i) {
    os << name << "_front " << i << " " << profile.fronts[i].position << " "
       << profile.fronts[i].width << " " << (profile.fronts[i].defined ? 1 : 0) << "\n";
  }
}

} // namespace

std::string formatStudy(const Study &study, const std::string &tag) {
  std::ostringstream os;
  os.setf(std::ios::fixed);
  os.precision(6);
  os << "# iface study\n";
  os << "tag " << tag << "\n";
  os << "n " << study.n << "\n";
  os << "disagree " << study.nDisagree << "\n";
  os << "cause neighbour_test " << study.causeCount[static_cast<int>(Cause::neighbour)] << "\n";
  os << "cause q3_shell " << study.causeCount[static_cast<int>(Cause::q3Shell)] << "\n";
  os << "cause open_edge " << study.causeCount[static_cast<int>(Cause::openEdge)] << "\n";
  os << "cause mixed " << study.causeCount[static_cast<int>(Cause::mixed)] << "\n";
  os << "kind " << (study.kind == Kind::slab ? "slab" : "nucleus") << "\n";
  os << "axis " << study.axis << "\n";
  for (int i = 0; i < kLabels; ++i) {
    os << "count_pub " << labelName(i) << " " << study.countPub[i] << "\n";
    os << "count_ovi " << labelName(i) << " " << study.countOvi[i] << "\n";
  }
  for (int a = 0; a < kLabels; ++a) {
    for (int b = 0; b < kLabels; ++b) {
      if (study.transition[a][b] == 0 || a == b) {
        continue;
      }
      os << "transition " << labelName(a) << " " << labelName(b) << " "
         << study.transition[a][b] << "\n";
    }
  }
  os << "shell_disagree";
  for (int s : study.shellDisagree) {
    os << " " << s;
  }
  os << "\nshell_agree";
  for (int s : study.shellAgree) {
    os << " " << s;
  }
  os << "\nring_disagree " << study.ringDisagree[0] << " " << study.ringDisagree[1] << " "
     << study.ringDisagree[2] << "\n";
  os << "ring_agree " << study.ringAgree[0] << " " << study.ringAgree[1] << " "
     << study.ringAgree[2] << "\n";
  os << "prism_atoms " << study.prismAtoms << "\n";
  os << "prism_blocks " << study.prismBlocks << "\n";
  os << "bilayer_atoms " << study.bilayerAtoms << "\n";
  os << "polygon_atoms " << study.polygonAtoms << "\n";
  os << "thin_sheets " << study.thinSheets << "\n";
  os << "diamond_cubic " << study.diamondCubic << "\n";
  os << "diamond_hex " << study.diamondHex << "\n";
  os << "diamond_other " << study.diamondOther << "\n";
  os << "cna_other " << study.cnaOther << "\n";
  os << "cna_crystal " << study.cnaCrystal << "\n";
  os << "ptm_pass " << study.ptmPass << "\n";
  os << "ptm_pass_disagree " << study.ptmPassDisagree << "\n";
  os << "q6_disagree " << study.q6Disagree << "\n";
  os << "q6_agree " << study.q6Agree << "\n";
  os << "depth_disagree " << study.depthDisagree << " " << study.depthDisagreeN << "\n";
  os << "depth_agree " << study.depthAgree << " " << study.depthAgreeN << "\n";
  writeProfile(os, "solid_pub", study.solidPub);
  writeProfile(os, "crystal_pub", study.crystalPub);
  writeProfile(os, "solid_ovi", study.solidOvi);
  writeProfile(os, "crystal_ovi", study.crystalOvi);
  os << "bin mid n disagree crystal_pub solid_pub crystal_ovi solid_ovi\n";
  for (const auto &bin : study.bins) {
    os << "bin " << bin.mid << " " << bin.n << " " << bin.disagree << " " << bin.crystalPub
       << " " << bin.solidPub << " " << bin.crystalOvi << " " << bin.solidOvi << "\n";
  }
  return os.str();
}

double slope(const std::vector<double> &x, const std::vector<double> &y) {
  const int n = static_cast<int>(std::min(x.size(), y.size()));
  if (n < 2) {
    return 0.0;
  }
  double mx = 0.0;
  double my = 0.0;
  for (int i = 0; i < n; ++i) {
    mx += x[static_cast<std::size_t>(i)];
    my += y[static_cast<std::size_t>(i)];
  }
  mx /= static_cast<double>(n);
  my /= static_cast<double>(n);
  double num = 0.0;
  double den = 0.0;
  for (int i = 0; i < n; ++i) {
    const double dx = x[static_cast<std::size_t>(i)] - mx;
    num += dx * (y[static_cast<std::size_t>(i)] - my);
    den += dx * dx;
  }
  if (!(den > 0.0)) {
    return 0.0;
  }
  return num / den;
}

Growth studyTrajectory(const std::string &path, double dumpEveryPs, double cutoff) {
  Growth growth;
  const int nFrames = countFrames(path);
  std::vector<double> time;
  std::vector<double> solidPub;
  std::vector<double> solidOvi;
  std::vector<double> crystalPub;
  std::vector<double> crystalOvi;
  for (int frame = 1; frame <= nFrames; ++frame) {
    GrowthFrame row;
    row.frame = frame;
    row.timePs = static_cast<double>(frame - 1) * dumpEveryPs;
    row.study = studyDump(path, frame, cutoff, 1.0, false);
    if (row.study.n <= 0) {
      break;
    }
    if (row.study.solidPub.defined) {
      time.push_back(row.timePs);
      solidPub.push_back(row.study.solidPub.thickness);
      solidOvi.push_back(row.study.solidOvi.thickness);
      crystalPub.push_back(row.study.crystalPub.thickness);
      crystalOvi.push_back(row.study.crystalOvi.thickness);
    }
    growth.frames.push_back(std::move(row));
  }
  growth.rateSolidPub = slope(time, solidPub);
  growth.rateSolidOvi = slope(time, solidOvi);
  growth.rateCrystalPub = slope(time, crystalPub);
  growth.rateCrystalOvi = slope(time, crystalOvi);
  return growth;
}

std::string formatGrowth(const Growth &growth) {
  std::ostringstream os;
  os.setf(std::ios::fixed);
  os.precision(6);
  os << "# iface growth\n";
  os << "nframes " << growth.frames.size() << "\n";
  os << "rate_solid_pub " << growth.rateSolidPub << "\n";
  os << "rate_solid_ovi " << growth.rateSolidOvi << "\n";
  os << "rate_crystal_pub " << growth.rateCrystalPub << "\n";
  os << "rate_crystal_ovi " << growth.rateCrystalOvi << "\n";
  os << "frame t_ps n disagree face_pub face_ovi solid_pub solid_ovi "
        "thick_solid_pub thick_solid_ovi thick_crystal_pub thick_crystal_ovi "
        "width_solid_pub width_solid_ovi pos0_pub pos1_pub pos0_ovi pos1_ovi\n";
  auto widthMean = [](const Profile &profile) {
    double sum = 0.0;
    int n = 0;
    for (const auto &front : profile.fronts) {
      if (front.defined && front.width > 0.0) {
        sum += front.width;
        ++n;
      }
    }
    return n > 0 ? sum / static_cast<double>(n) : 0.0;
  };
  auto posAt = [](const Profile &profile, std::size_t i) {
    return i < profile.fronts.size() ? profile.fronts[i].position : 0.0;
  };
  for (const auto &row : growth.frames) {
    const auto &s = row.study;
    os << "frame " << row.frame << " " << row.timePs << " " << s.n << " " << s.nDisagree << " "
       << s.solidPub.interfacial << " " << s.solidOvi.interfacial << " " << s.solidPub.solid
       << " " << s.solidOvi.solid << " " << s.solidPub.thickness << " " << s.solidOvi.thickness
       << " " << s.crystalPub.thickness << " " << s.crystalOvi.thickness << " "
       << widthMean(s.solidPub) << " " << widthMean(s.solidOvi) << " " << posAt(s.solidPub, 0)
       << " " << posAt(s.solidPub, 1) << " " << posAt(s.solidOvi, 0) << " "
       << posAt(s.solidOvi, 1) << "\n";
  }
  return os.str();
}

std::string formatQ6Growth(const std::string &path, double dumpEveryPs, double cutoff) {
  const int nFrames = countFrames(path);
  const int nBins = 50;
  std::vector<double> thresholds;
  for (int i = 0; i <= 40; ++i) {
    thresholds.push_back(0.20 + 0.01 * static_cast<double>(i));
  }
  std::ostringstream os;
  os.setf(std::ios::fixed);
  os.precision(6);
  os << "# iface q6\n";
  os << "cutoff " << cutoff << "\n";
  os << "nframes " << nFrames << "\n";
  os << "thresholds";
  for (double t : thresholds) {
    os << " " << t;
  }
  os << "\n";
  os << "# frame t_ps n axis length n_ice n_liq bar <counts> local <counts>\n";
  os << "# hist_ice and hist_liq are frame 1, bins of width 0.02 from 0 to 1,\n";
  os << "# ice core is the central fifth of the long axis; liquid core is the outer tenths.\n";
  bool wroteHist = false;
  for (int frame = 1; frame <= nFrames; ++frame) {
    Cloud cloud;
    cloud = sinp::readLammpsTrjO(path, frame, cloud, 1);
    if (cloud.nop <= 0) {
      break;
    }
    const auto q6 = chill::steinhardtQlCutoff(cloud, cutoff, 6, 1, true);
    int axis = 0;
    double length = lengthOf(cloud, 0);
    for (int a = 1; a < 3; ++a) {
      const double edge = lengthOf(cloud, a);
      if (edge > length) {
        length = edge;
        axis = a;
      }
    }
    const double origin = originOf(cloud, axis);
    std::vector<int> iceHist(static_cast<std::size_t>(nBins), 0);
    std::vector<int> liqHist(static_cast<std::size_t>(nBins), 0);
    int nIce = 0;
    int nLiq = 0;
    std::vector<int> barCount(thresholds.size(), 0);
    std::vector<int> localCount(thresholds.size(), 0);
    for (int i = 0; i < cloud.nop; ++i) {
      double q = axisCoord(cloud.pts[static_cast<std::size_t>(i)], axis) - origin;
      q = std::fmod(q, length);
      if (q < 0.0) {
        q += length;
      }
      const double frac = length > 0.0 ? q / length : 0.0;
      const bool ice = frac >= 0.40 && frac < 0.60;
      const bool liq = frac < 0.10 || frac >= 0.90;
      const double bar = i < static_cast<int>(q6.qlBar.size()) ? q6.qlBar[static_cast<std::size_t>(i)] : 0.0;
      const double local = i < static_cast<int>(q6.ql.size()) ? q6.ql[static_cast<std::size_t>(i)] : 0.0;
      if (ice) {
        ++nIce;
      }
      if (liq) {
        ++nLiq;
      }
      if (frame == 1 && (ice || liq) && std::isfinite(bar)) {
        int bin = static_cast<int>(bar / 0.02);
        if (bin < 0) {
          bin = 0;
        }
        if (bin >= nBins) {
          bin = nBins - 1;
        }
        if (ice) {
          ++iceHist[static_cast<std::size_t>(bin)];
        }
        if (liq) {
          ++liqHist[static_cast<std::size_t>(bin)];
        }
      }
      for (std::size_t k = 0; k < thresholds.size(); ++k) {
        if (std::isfinite(bar) && bar >= thresholds[k]) {
          ++barCount[k];
        }
        if (std::isfinite(local) && local >= thresholds[k]) {
          ++localCount[k];
        }
      }
    }
    os << "frame " << frame << " " << static_cast<double>(frame - 1) * dumpEveryPs << " " << cloud.nop
       << " " << axis << " " << length << " " << nIce << " " << nLiq << " bar";
    for (int c : barCount) {
      os << " " << c;
    }
    os << " local";
    for (int c : localCount) {
      os << " " << c;
    }
    os << "\n";
    if (frame == 1) {
      os << "hist_ice";
      for (int c : iceHist) {
        os << " " << c;
      }
      os << "\nhist_liq";
      for (int c : liqHist) {
        os << " " << c;
      }
      os << "\n";
      wroteHist = true;
    }
  }
  if (!wroteHist) {
    os << "hist_ice";
    for (int i = 0; i < nBins; ++i) {
      os << " 0";
    }
    os << "\nhist_liq";
    for (int i = 0; i < nBins; ++i) {
      os << " 0";
    }
    os << "\n";
  }
  return os.str();
}

Cloud nanotube(int ring, int stacks, double bond) {
  Cloud cloud;
  const double pi = std::acos(-1.0);
  const double radius = bond / (2.0 * std::sin(pi / static_cast<double>(ring)));
  const double length = static_cast<double>(stacks) * bond;
  cloud.box = {12.0, 12.0, length};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  for (int s = 0; s < stacks; ++s) {
    for (int k = 0; k < ring; ++k) {
      const double ang = 2.0 * pi * static_cast<double>(k) / static_cast<double>(ring);
      addPoint(cloud, 6.0 + radius * std::cos(ang), 6.0 + radius * std::sin(ang),
               static_cast<double>(s) * bond);
    }
  }
  return cloud;
}

Cloud honeycombSheet(int nx, int ny, int layers, double bond) {
  Cloud cloud;
  const double a = bond * std::sqrt(3.0);
  const double height = a * std::sqrt(3.0);
  cloud.box = {static_cast<double>(nx) * a, static_cast<double>(ny) * height, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  const double lx = cloud.box[0];
  const double ly = cloud.box[1];
  for (int layer = 0; layer < layers; ++layer) {
    const double z = 5.0 + static_cast<double>(layer) * bond;
    const double shift = layer == 0 ? 0.0 : bond;
    for (int ix = 0; ix < nx; ++ix) {
      for (int iy = 0; iy < ny; ++iy) {
        const double x0 = static_cast<double>(ix) * a;
        const double y0 = static_cast<double>(iy) * height;
        const double site[4][2] = {
            {0.0, 0.0},
            {0.0, bond},
            {0.5 * a, 0.5 * height},
            {0.5 * a, 0.5 * height + bond},
        };
        for (const auto &s : site) {
          addPoint(cloud, wrap(x0 + s[0], lx), wrap(y0 + s[1] + shift, ly), z);
        }
      }
    }
  }
  return cloud;
}

} // namespace iface
