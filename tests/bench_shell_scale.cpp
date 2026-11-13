//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Raw samples for the in-tree ice kernel. Release -O3. One warm-up, then
// eight calls, and at least six. Threads come from argv (default 1).

#include <local_symmetry.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif
#ifdef SEAMS_HAS_IRA
#include <iralib_interf.h>
#endif

namespace {

constexpr double kBond = 2.75;

Eigen::MatrixXd iceIc(int repetitions, Eigen::Matrix3d &box) {
  const double lattice = 4.0 * kBond / std::sqrt(3.0);
  const int n = 8 * repetitions * repetitions * repetitions;
  Eigen::MatrixXd xyz(n, 3);
  const double basis[8][3] = {{0, 0, 0},         {0.5, 0.5, 0}, {0.5, 0, 0.5},
                              {0, 0.5, 0.5},     {0.25, 0.25, 0.25},
                              {0.75, 0.75, 0.25},{0.75, 0.25, 0.75},
                              {0.25, 0.75, 0.75}};
  int row = 0;
  for (int ix = 0; ix < repetitions; ix++) {
    for (int iy = 0; iy < repetitions; iy++) {
      for (int iz = 0; iz < repetitions; iz++) {
        for (int b = 0; b < 8; b++) {
          xyz(row, 0) = (ix + basis[b][0]) * lattice;
          xyz(row, 1) = (iy + basis[b][1]) * lattice;
          xyz(row, 2) = (iz + basis[b][2]) * lattice;
          row++;
        }
      }
    }
  }
  box = Eigen::Matrix3d::Identity() * (repetitions * lattice);
  return xyz;
}

double msSince(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

void quantile(const std::vector<double> &sorted, double q, double &out) {
  if (sorted.size() == 1) {
    out = sorted[0];
    return;
  }
  const double pos = q * static_cast<double>(sorted.size() - 1);
  const size_t lo = static_cast<size_t>(pos);
  const double frac = pos - static_cast<double>(lo);
  const size_t hi = std::min(lo + 1, sorted.size() - 1);
  out = sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

void report(const char *label, int molecules, int threads, std::vector<double> samples) {
  std::sort(samples.begin(), samples.end());
  double q1 = 0;
  double q2 = 0;
  double q3 = 0;
  quantile(samples, 0.25, q1);
  quantile(samples, 0.50, q2);
  quantile(samples, 0.75, q3);
  const double mid = samples.size() % 2 == 0
                         ? 0.5 * (samples[samples.size() / 2 - 1] + samples[samples.size() / 2])
                         : samples[samples.size() / 2];
  std::cout << label << " molecules=" << molecules << " threads=" << threads
            << " n=" << samples.size() << " median_ms=" << mid << " q1_ms=" << q1
            << " q3_ms=" << q3 << " iqr_ms=" << (q3 - q1) << " per_atom_us="
            << (mid * 1000.0 / static_cast<double>(molecules)) << " raw";
  for (double v : samples) {
    std::cout << " " << v;
  }
  std::cout << "\n";
  (void)q2;
}

#ifdef SEAMS_HAS_IRA
double timeIraAtoms(const Eigen::MatrixXd &xyz, const Eigen::Matrix3d &box, int sample) {
  const int n = std::min(sample, static_cast<int>(xyz.rows()));
  const auto t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < n; i++) {
    const auto shell = symm::iceShell(xyz, box, i);
    if (shell.rows() != 17) {
      continue;
    }
    (void)symm::best(shell, "ira");
  }
  return msSince(t0);
}
#endif

} // namespace

int main(int argc, char **argv) {
  int threads = 1;
  if (argc > 1) {
    threads = std::atoi(argv[1]);
  }
  if (threads < 1) {
    threads = 1;
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(threads);
#endif
  const int reps[] = {3, 11, 24};
  constexpr int kWarm = 1;
  constexpr int kN = 8;
  std::cout << "build release -O3 kernel threads_requested=" << threads << "\n";
  for (int rep : reps) {
    Eigen::Matrix3d box;
    const Eigen::MatrixXd xyz = iceIc(rep, box);
    const int molecules = static_cast<int>(xyz.rows());
    std::vector<double> samples;
    int labelled = 0;
    for (int i = 0; i < kWarm + kN; i++) {
      std::vector<int> kind;
      const auto t0 = std::chrono::steady_clock::now();
      symm::classifyIce(xyz, box, kind);
      const double dt = msSince(t0);
      if (i >= kWarm) {
        samples.push_back(dt);
      }
      labelled = 0;
      for (int k : kind) {
        if (k == static_cast<int>(symm::Kind::iceIc)) {
          labelled++;
        }
      }
    }
    std::cout << "kernel_ice_ic labelled " << labelled << " of " << molecules << "\n";
    report("kernel", molecules, threads, samples);

#ifdef SEAMS_HAS_IRA
    if (molecules == 216 && threads == 1) {
      std::vector<double> iraSamples;
      for (int i = 0; i < kWarm + kN; i++) {
        const double dt = timeIraAtoms(xyz, box, molecules);
        if (i >= kWarm) {
          iraSamples.push_back(dt);
        }
      }
      report("ira_full_cell", molecules, 1, iraSamples);
    }
#else
    std::cout << "ira off\n";
#endif
  }
  if (threads == 1) {
    const auto shell = symm::points(symm::Kind::iceIc);
    std::vector<double> iraShell;
    std::vector<double> cshm;
    constexpr int kCalls = 40;
    for (int i = 0; i < kWarm + kN; i++) {
      const auto t0 = std::chrono::steady_clock::now();
      for (int c = 0; c < kCalls; c++) {
        (void)symm::best(shell, "ira");
      }
      const double dt = msSince(t0) * static_cast<double>(216) / kCalls;
      if (i >= kWarm) {
        iraShell.push_back(dt);
      }
    }
    std::cout << "ira_note " << kCalls
              << " library matches of one 17-atom shell, scaled to 216 atoms, serial\n";
    report("ira_shell_scaled216", 216, 1, iraShell);
    for (int i = 0; i < kWarm + kN; i++) {
      const auto t0 = std::chrono::steady_clock::now();
      for (int c = 0; c < kCalls; c++) {
        (void)symm::shapeMeasure(shell);
      }
      const double dt = msSince(t0) * static_cast<double>(216) / kCalls;
      if (i >= kWarm) {
        cshm.push_back(dt);
      }
    }
    std::cout << "cshm_note in-tree Alvarez measure on one shell, scaled to 216, serial\n";
    report("cshm_shell_scaled216", 216, 1, cshm);
    const auto cell = symm::crystal(symm::Crystal::iceIc);
    std::vector<double> spg;
    for (int i = 0; i < kWarm + kN; i++) {
      const auto t0 = std::chrono::steady_clock::now();
      const auto space = symm::spaceGroup(cell.vectors, cell.cartesian, 1e-3);
      const double dt = msSince(t0);
      if (i == kWarm) {
        std::cout << "spglib " << space.number << " " << space.symbol << "\n";
      }
      if (i >= kWarm) {
        spg.push_back(dt);
      }
    }
    report("spglib_ice_ic_cell", static_cast<int>(cell.cartesian.rows()), 1, spg);
  }
  return 0;
}
