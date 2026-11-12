//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Times the SVD-and-permutation namer on ice shells, cages and rings.
// Release -O3. One warm-up, then six calls. Threads come from argv.

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

namespace {

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

void report(const char *label, int nshell, int threads, std::vector<double> samples) {
  std::sort(samples.begin(), samples.end());
  double q1 = 0;
  double q3 = 0;
  quantile(samples, 0.25, q1);
  quantile(samples, 0.75, q3);
  const double mid = samples.size() % 2 == 0
                         ? 0.5 * (samples[samples.size() / 2 - 1] + samples[samples.size() / 2])
                         : samples[samples.size() / 2];
  std::cout << label << " shells=" << nshell << " threads=" << threads
            << " n=" << samples.size() << " median_ms=" << mid << " q1_ms=" << q1
            << " q3_ms=" << q3 << " iqr_ms=" << (q3 - q1) << " per_shell_us="
            << (mid * 1000.0 / static_cast<double>(nshell)) << " raw";
  for (double v : samples) {
    std::cout << " " << v;
  }
  std::cout << "\n";
}

Eigen::Matrix3d rotation(int seed) {
  const double a = 0.37 * static_cast<double>(seed + 1);
  const double c = std::cos(a);
  const double s = std::sin(a);
  Eigen::Vector3d axis(1.0, 0.4, 0.2);
  axis.normalize();
  Eigen::Matrix3d K;
  K << 0, -axis[2], axis[1], axis[2], 0, -axis[0], -axis[1], axis[0], 0;
  return Eigen::Matrix3d::Identity() + s * K + (1.0 - c) * (K * K);
}

std::vector<Eigen::MatrixXd> spun(symm::Kind kind, int count) {
  const Eigen::MatrixXd src = symm::points(kind);
  std::vector<Eigen::MatrixXd> out;
  out.reserve(static_cast<size_t>(count));
  for (int i = 0; i < count; i++) {
    const Eigen::Matrix3d R = rotation(i);
    Eigen::MatrixXd env(src.rows(), 3);
    for (Eigen::Index row = 0; row < src.rows(); row++) {
      env.row(row) = (R * src.row(row).transpose()).transpose();
      env(row, 0) += 0.04 * std::sin(0.2 * static_cast<double>(i + row));
    }
    out.push_back(std::move(env));
  }
  return out;
}

void timeKernel(const char *label, symm::Kind kind, int count, int threads) {
  const std::vector<Eigen::MatrixXd> envs = spun(kind, count);
  std::vector<symm::Overlay> named;
  symm::matchAll(envs, named);
  std::vector<double> samples;
  for (int n = 0; n < 6; n++) {
    const auto t0 = std::chrono::steady_clock::now();
    symm::matchAll(envs, named);
    samples.push_back(msSince(t0));
  }
  int namedOk = 0;
  for (const symm::Overlay &hit : named) {
    if (hit.status == 0 && hit.kind == kind) {
      namedOk++;
    }
  }
  std::cout << label << " recognised " << namedOk << " of " << count << "\n";
  report(label, count, threads, samples);
}

void timeIra(const char *label, symm::Kind kind, int count) {
  const std::vector<Eigen::MatrixXd> envs = spun(kind, count);
  (void)symm::best(envs[0], "ira");
  std::vector<double> samples;
  for (int n = 0; n < 6; n++) {
    const auto t0 = std::chrono::steady_clock::now();
    for (const Eigen::MatrixXd &env : envs) {
      (void)symm::best(env, "ira");
    }
    samples.push_back(msSince(t0));
  }
  report(label, count, 1, samples);
}

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
  std::cout << "naming threads=" << threads << "\n";
  timeKernel("kernel-ice-Ic", symm::Kind::iceIc, 216, threads);
  timeKernel("kernel-ice-Ic", symm::Kind::iceIc, 2048, threads);
  timeKernel("kernel-512", symm::Kind::cage512, 216, threads);
  timeKernel("kernel-51262", symm::Kind::cage51262, 216, threads);
  timeKernel("kernel-ring6", symm::Kind::ring6, 216, threads);
  timeKernel("kernel-ring6", symm::Kind::ring6, 2048, threads);
  {
    const std::vector<Eigen::MatrixXd> envs = spun(symm::Kind::iceIc, 216);
    (void)symm::shapeMeasure(envs[0]);
    std::vector<double> samples;
    for (int n = 0; n < 6; n++) {
      const auto t0 = std::chrono::steady_clock::now();
      double sink = 0;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) reduction(+ : sink)
#endif
      for (int i = 0; i < static_cast<int>(envs.size()); i++) {
        sink += symm::shapeMeasure(envs[static_cast<size_t>(i)]);
      }
      samples.push_back(msSince(t0));
      if (sink == 1e300) {
        std::cout << "unused\n";
      }
    }
    report("kernel-cshm", 216, threads, samples);
  }
  if (threads == 1 && symm::iraLinked()) {
    timeIra("ira-ice", symm::Kind::iceIc, 216);
    timeIra("ira-cage", symm::Kind::cage512, 64);
    timeIra("ira-ring", symm::Kind::ring6, 216);
  }
  std::cout << "ptm has no cage or ring name\n";
  return 0;
}
