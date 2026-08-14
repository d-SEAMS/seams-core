//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <ira_sofi.hpp>

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#ifdef SEAMS_HAS_IRA
#include <iralib_interf.h>
#endif

namespace {

void flatten(const Eigen::MatrixXd &pts, std::vector<double> &coords,
             std::vector<int> &typ) {
  const int n = static_cast<int>(pts.rows());
  coords.resize(static_cast<size_t>(3 * n));
  typ.assign(static_cast<size_t>(n), 1);
  for (int i = 0; i < n; i++) {
    coords[static_cast<size_t>(3 * i)] = pts(i, 0);
    coords[static_cast<size_t>(3 * i + 1)] = pts(i, 1);
    coords[static_cast<size_t>(3 * i + 2)] = pts(i, 2);
  }
}

#ifdef SEAMS_HAS_IRA
// ira_mod.py fills rotation[i][j] from the flat buffer in C row-major order.
Eigen::Matrix3d rotationFromC(const double *rmat) {
  Eigen::Matrix3d R;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      R(i, j) = rmat[3 * i + j];
    }
  }
  return R;
}

void fillQuat(const Eigen::Matrix3d &R, std::vector<double> &quat) {
  const Eigen::Quaterniond q(R);
  quat = {q.w(), q.x(), q.y(), q.z()};
}
#endif

} // namespace

namespace ira {

bool available() {
#ifdef SEAMS_HAS_IRA
  return true;
#else
  return false;
#endif
}

int match(const Eigen::MatrixXd &ref, const Eigen::MatrixXd &target, Match &out,
          double kmaxFactor) {
#ifndef SEAMS_HAS_IRA
  (void)ref;
  (void)target;
  (void)out;
  (void)kmaxFactor;
  return 1;
#else
  if (ref.cols() != 3 || target.cols() != 3 || ref.rows() < 1 ||
      target.rows() < 1) {
    return 1;
  }
  const int nat1 = static_cast<int>(ref.rows());
  const int nat2 = static_cast<int>(target.rows());
  std::vector<double> coords1;
  std::vector<double> coords2;
  std::vector<int> typ1;
  std::vector<int> typ2;
  flatten(ref, coords1, typ1);
  flatten(target, coords2, typ2);
  std::vector<int> cand1(static_cast<size_t>(nat1), 0);
  std::vector<int> cand2(static_cast<size_t>(nat2), 0);
  if (nat1 == nat2) {
    cand1[0] = -1;
    cand2[0] = -1;
  } else {
    cand1[0] = 1;
    for (int i = 0; i < nat2; i++) {
      cand2[static_cast<size_t>(i)] = i + 1;
    }
  }
  std::vector<double> rmat(9, 0.0);
  std::vector<double> tr(3, 0.0);
  std::vector<int> perm(static_cast<size_t>(nat2), 0);
  double *rmatP = rmat.data();
  double *trP = tr.data();
  int *permP = perm.data();
  double hd = 0.0;
  int err = 0;
  libira_match(nat1, typ1.data(), coords1.data(), cand1.data(), nat2,
               typ2.data(), coords2.data(), cand2.data(), kmaxFactor, &rmatP,
               &trP, &permP, &hd, &err);
  if (err != 0) {
    return err;
  }
  out.rotation = rotationFromC(rmatP);
  out.translation = Eigen::Vector3d(trP[0], trP[1], trP[2]);
  out.hausdorff = hd;
  out.assignment.assign(permP, permP + nat2);
  fillQuat(out.rotation, out.quat);
  const int nPair = std::min(nat1, nat2);
  bool hasZero = false;
  bool hasN = false;
  for (int p : out.assignment) {
    if (p == 0) {
      hasZero = true;
    }
    if (p == nat2) {
      hasN = true;
    }
  }
  const int shift = (hasN && !hasZero) ? 1 : 0;
  double sum = 0.0;
  int used = 0;
  for (int i = 0; i < nPair; i++) {
    const int j = out.assignment[static_cast<size_t>(i)] - shift;
    if (j < 0 || j >= nat2) {
      continue;
    }
    const Eigen::Vector3d a(ref(i, 0), ref(i, 1), ref(i, 2));
    const Eigen::Vector3d b(target(j, 0), target(j, 1), target(j, 2));
    const Eigen::Vector3d d = out.rotation * a + out.translation - b;
    sum += d.squaredNorm();
    used++;
  }
  out.rmsd = used > 0 ? std::sqrt(sum / static_cast<double>(used)) : 0.0;
  return 0;
#endif
}

int pointGroup(const Eigen::MatrixXd &coords, PointGroup &out, double symThr) {
#ifndef SEAMS_HAS_IRA
  (void)coords;
  (void)out;
  (void)symThr;
  return 1;
#else
  if (coords.cols() != 3 || coords.rows() < 1) {
    return 1;
  }
  const int nat = static_cast<int>(coords.rows());
  std::vector<double> xyz;
  std::vector<int> typ;
  flatten(coords, xyz, typ);
  Eigen::Vector3d gc = Eigen::Vector3d::Zero();
  for (int i = 0; i < nat; i++) {
    gc += Eigen::Vector3d(xyz[static_cast<size_t>(3 * i)],
                          xyz[static_cast<size_t>(3 * i + 1)],
                          xyz[static_cast<size_t>(3 * i + 2)]);
  }
  gc /= static_cast<double>(nat);
  for (int i = 0; i < nat; i++) {
    xyz[static_cast<size_t>(3 * i)] -= gc[0];
    xyz[static_cast<size_t>(3 * i + 1)] -= gc[1];
    xyz[static_cast<size_t>(3 * i + 2)] -= gc[2];
  }
  constexpr int kMax = 400;
  std::vector<double> mat(static_cast<size_t>(9 * kMax), 0.0);
  std::vector<int> perm(static_cast<size_t>(nat * kMax), 0);
  std::vector<char> op(static_cast<size_t>(kMax + 1), '\0');
  std::vector<int> nData(kMax, 0);
  std::vector<int> pData(kMax, 0);
  std::vector<double> ax(static_cast<size_t>(3 * kMax), 0.0);
  std::vector<double> angle(kMax, 0.0);
  std::vector<double> dh(kMax, 0.0);
  std::vector<char> pg(11, '\0');
  std::vector<double> prin(static_cast<size_t>(3 * kMax), 0.0);
  int nmat = 0;
  int nPrin = 0;
  int err = 0;
  double *matP = mat.data();
  int *permP = perm.data();
  char *opP = op.data();
  int *nP = nData.data();
  int *pP = pData.data();
  double *axP = ax.data();
  double *angP = angle.data();
  double *dhP = dh.data();
  char *pgP = pg.data();
  double *prinP = prin.data();
  libira_compute_all(nat, typ.data(), xyz.data(), symThr, 0, &nmat, &matP,
                     &permP, &opP, &nP, &pP, &axP, &angP, &dhP, &pgP, &nPrin,
                     &prinP, &err);
  if (err < 0) {
    return err;
  }
  out.symbol = pgP;
  out.nOperations = nmat;
  out.operations.clear();
  out.operations.reserve(static_cast<size_t>(nmat));
  for (int k = 0; k < nmat; k++) {
    out.operations.push_back(rotationFromC(matP + 9 * k));
  }
  return 0;
#endif
}

bool orient(const Eigen::MatrixXd &ref, const Eigen::MatrixXd &target,
            std::vector<double> &quat, double &rmsd) {
  Match m;
  if (match(ref, target, m) != 0) {
    return false;
  }
  quat = m.quat;
  rmsd = m.rmsd;
  return true;
}

} // namespace ira
