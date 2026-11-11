//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Serial timings of the symmetry wrappers against the libraries they call.
// Release build, one thread. Prints median and interquartile range.

#include <local_symmetry.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_IRA
#include <iralib_interf.h>
#endif
#ifdef SEAMS_HAS_SPGLIB
#include <spglib.h>
#endif

namespace {

double percentile(std::vector<double> v, double p) {
  if (v.empty()) {
    return 0;
  }
  std::sort(v.begin(), v.end());
  const double idx = (p / 100.0) * static_cast<double>(v.size() - 1);
  const size_t lo = static_cast<size_t>(std::floor(idx));
  const size_t hi = static_cast<size_t>(std::ceil(idx));
  const double w = idx - static_cast<double>(lo);
  return v[lo] * (1.0 - w) + v[hi] * w;
}

template <typename Fn> std::vector<double> repeat(int n, Fn &&fn) {
  fn();
  std::vector<double> samples;
  samples.reserve(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    const auto t0 = std::chrono::steady_clock::now();
    fn();
    const auto t1 = std::chrono::steady_clock::now();
    samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
  }
  return samples;
}

void report(const char *label, const std::vector<double> &samples) {
  const double med = percentile(samples, 50);
  const double iqr = percentile(samples, 75) - percentile(samples, 25);
  std::cout << label << " n=" << samples.size() << " median_us=" << med
            << " iqr_us=" << iqr << "\n";
}

Eigen::MatrixXd shifted(const Eigen::MatrixXd &src, int k) {
  Eigen::MatrixXd out = src;
  for (Eigen::Index i = 0; i < out.rows(); i++) {
    out(i, 0) += 0.01 * std::sin(static_cast<double>(i + k));
    out(i, 1) += 0.01 * std::cos(static_cast<double>(i + k));
  }
  return out;
}

#ifdef SEAMS_HAS_IRA
void directIra(const Eigen::MatrixXd &ref, const Eigen::MatrixXd &target) {
  const int n = static_cast<int>(ref.rows());
  std::vector<double> c1(static_cast<size_t>(3 * n));
  std::vector<double> c2(static_cast<size_t>(3 * n));
  std::vector<int> typ(static_cast<size_t>(n), 1);
  for (int i = 0; i < n; i++) {
    c1[static_cast<size_t>(3 * i)] = ref(i, 0);
    c1[static_cast<size_t>(3 * i + 1)] = ref(i, 1);
    c1[static_cast<size_t>(3 * i + 2)] = ref(i, 2);
    c2[static_cast<size_t>(3 * i)] = target(i, 0);
    c2[static_cast<size_t>(3 * i + 1)] = target(i, 1);
    c2[static_cast<size_t>(3 * i + 2)] = target(i, 2);
  }
  std::vector<int> cand(static_cast<size_t>(n), 0);
  cand[0] = -1;
  std::vector<double> rmat(9, 0.0);
  std::vector<double> tr(3, 0.0);
  std::vector<int> perm(static_cast<size_t>(n), 0);
  double *rmatP = rmat.data();
  double *trP = tr.data();
  int *permP = perm.data();
  double hd = 0;
  int err = 0;
  libira_match(n, typ.data(), c1.data(), cand.data(), n, typ.data(), c2.data(), cand.data(),
               1.8, &rmatP, &trP, &permP, &hd, &err);
  if (err != 0) {
    std::cerr << "direct ira err " << err << "\n";
  }
}

void directSofi(const Eigen::MatrixXd &xyz) {
  const int n = static_cast<int>(xyz.rows());
  std::vector<double> c(static_cast<size_t>(3 * n));
  std::vector<int> typ(static_cast<size_t>(n), 1);
  for (int i = 0; i < n; i++) {
    c[static_cast<size_t>(3 * i)] = xyz(i, 0);
    c[static_cast<size_t>(3 * i + 1)] = xyz(i, 1);
    c[static_cast<size_t>(3 * i + 2)] = xyz(i, 2);
  }
  constexpr int kMax = 400;
  std::vector<double> mat(static_cast<size_t>(9 * kMax));
  std::vector<int> perm(static_cast<size_t>(n * kMax));
  std::vector<char> op(static_cast<size_t>(kMax + 1), 0);
  std::vector<int> nd(static_cast<size_t>(kMax));
  std::vector<int> pd(static_cast<size_t>(kMax));
  std::vector<double> ax(static_cast<size_t>(3 * kMax));
  std::vector<double> ang(static_cast<size_t>(kMax));
  std::vector<double> dh(static_cast<size_t>(kMax));
  std::vector<double> prin(static_cast<size_t>(3 * kMax));
  std::vector<char> pg(16, 0);
  int nmat = 0;
  int nprin = 0;
  int err = 0;
  double *matP = mat.data();
  int *permP = perm.data();
  char *opP = op.data();
  int *nP = nd.data();
  int *pP = pd.data();
  double *axP = ax.data();
  double *angP = ang.data();
  double *dhP = dh.data();
  char *pgP = pg.data();
  double *prinP = prin.data();
  libira_compute_all(n, typ.data(), c.data(), 0.08, 0, &nmat, &matP, &permP, &opP, &nP, &pP,
                     &axP, &angP, &dhP, &pgP, &nprin, &prinP, &err);
  if (err < 0) {
    std::cerr << "direct sofi err " << err << "\n";
  }
}
#endif

#ifdef SEAMS_HAS_SPGLIB
void directSpglib(const symm::Lattice &cell) {
  const int n = static_cast<int>(cell.cartesian.rows());
  double lattice[3][3];
  for (int comp = 0; comp < 3; comp++) {
    for (int vec = 0; vec < 3; vec++) {
      lattice[comp][vec] = cell.vectors(vec, comp);
    }
  }
  std::vector<double> raw(static_cast<size_t>(3 * n));
  std::vector<int> types(static_cast<size_t>(n), 1);
  const Eigen::Matrix3d inv = cell.vectors.transpose().inverse();
  for (int i = 0; i < n; i++) {
    Eigen::Vector3d f = inv * cell.cartesian.row(i).transpose();
    for (int k = 0; k < 3; k++) {
      f[k] -= std::floor(f[k]);
      raw[static_cast<size_t>(3 * i + k)] = f[k];
    }
  }
  SpglibDataset *data = spg_get_dataset(
      lattice, reinterpret_cast<double (*)[3]>(raw.data()), types.data(), n, 1e-3);
  if (data == nullptr || data->spacegroup_number <= 0) {
    std::cerr << "direct spglib failed\n";
  }
  if (data != nullptr) {
    spg_free_dataset(data);
  }
}
#endif

} // namespace

int main() {
  constexpr int kRepeats = 8;
  constexpr int kMatch = 40;
  std::cout << "matcher " << symm::matcher() << " ira " << symm::iraLinked() << " spglib "
            << symm::spglibLinked() << "\n";
  const auto ice = symm::points(symm::Kind::iceIh);
  const auto cage = symm::points(symm::Kind::cage512);
  const auto cell = symm::crystal(symm::Crystal::iceIh);
  std::vector<Eigen::MatrixXd> copies;
  copies.reserve(kMatch);
  for (int k = 0; k < kMatch; k++) {
    copies.push_back(shifted(k % 2 == 0 ? ice : cage, k));
  }

  report("overlay", repeat(kRepeats, [&] {
           for (const auto &env : copies) {
             (void)symm::best(env);
           }
         }));
  std::vector<Eigen::MatrixXd> cages;
  cages.reserve(kMatch);
  for (int k = 0; k < kMatch; k++) {
    cages.push_back(shifted(cage, k));
  }
  report("overlay-one-template", repeat(kRepeats, [&] {
           for (const auto &env : cages) {
             (void)symm::best(env);
           }
         }));
#ifdef SEAMS_HAS_IRA
  report("direct-ira", repeat(kRepeats, [&] {
           for (int k = 0; k < kMatch; k++) {
             const auto &ref = k % 2 == 0 ? ice : cage;
             directIra(ref, copies[static_cast<size_t>(k)]);
           }
         }));
  report("direct-ira-one-template", repeat(kRepeats, [&] {
           for (const auto &env : cages) {
             directIra(cage, env);
           }
         }));
  report("point-group", repeat(kRepeats, [&] {
           for (int k = 0; k < 8; k++) {
             (void)symm::pointGroup(cage, 0.08);
           }
         }));
  report("direct-sofi", repeat(kRepeats, [&] {
           for (int k = 0; k < 8; k++) {
             directSofi(cage);
           }
         }));
#else
  std::cout << "direct-ira off\npoint-group off\ndirect-sofi off\n";
#endif
#ifdef SEAMS_HAS_SPGLIB
  report("space-group", repeat(kRepeats, [&] {
           for (int k = 0; k < 20; k++) {
             (void)symm::spaceGroup(cell.vectors, cell.cartesian, 1e-3);
           }
         }));
  report("direct-spglib", repeat(kRepeats, [&] {
           for (int k = 0; k < 20; k++) {
             directSpglib(cell);
           }
         }));
  {
    const auto ic = symm::crystal(symm::Crystal::iceIc);
    const int nrep = 3;
    const int base = static_cast<int>(ic.cartesian.rows());
    Eigen::MatrixXd big(base * nrep * nrep * nrep, 3);
    int row = 0;
    for (int ix = 0; ix < nrep; ix++) {
      for (int iy = 0; iy < nrep; iy++) {
        for (int iz = 0; iz < nrep; iz++) {
          const Eigen::Vector3d shift = static_cast<double>(ix) * ic.vectors.row(0).transpose() +
                                        static_cast<double>(iy) * ic.vectors.row(1).transpose() +
                                        static_cast<double>(iz) * ic.vectors.row(2).transpose();
          for (int i = 0; i < base; i++) {
            big.row(row++) = ic.cartesian.row(i) + shift.transpose();
          }
        }
      }
    }
    Eigen::Matrix3d box = ic.vectors * static_cast<double>(nrep);
    int labelled = 0;
    std::string name;
    report("ice-ic-shells-216", repeat(kRepeats, [&] {
             labelled = 0;
             for (int i = 0; i < big.rows(); i++) {
               const Eigen::MatrixXd shell = symm::iceShell(big, box, i);
               if (shell.rows() != 17) {
                 continue;
               }
               const symm::Overlay hit = symm::best(shell);
               if (hit.status == 0) {
                 labelled++;
                 name = hit.name;
               }
             }
           }));
    std::cout << "ice-ic-shells labelled " << labelled << " as " << name << "\n";
  }
#else
  std::cout << "space-group off\ndirect-spglib off\n";
#endif
  return 0;
}
