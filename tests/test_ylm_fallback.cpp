// Copyright (c) d-SEAMS developers
#include <catch2/catch_test_macros.hpp>

#include <bop.hpp>
#include <cpu_features.hpp>
#include <local_symmetry.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <shell_hot.hpp>
#include <ylm_batch.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>

#if defined(SEAMS_HAS_OPENMP)
#include <omp.h>
#endif

// Every call into the AVX-512 kernel archive goes through
// seams::ylm::avx512Kernels, every call into the AVX2 shell archive through
// symm::avx2Kernels, and the AVX2 cutoff stash scan checks the same CPU
// bits. With the kernels turned off the generic path runs on any CPU. Its
// values must match the kernel's, and its labels and neighbour lists must
// be the same.

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

struct ScalarOnly {
  ScalarOnly() { seams::cpu::allowSimdKernels(false); }
  ~ScalarOnly() { seams::cpu::allowSimdKernels(true); }
  ScalarOnly(const ScalarOnly &) = delete;
  ScalarOnly &operator=(const ScalarOnly &) = delete;
};

Cloud iceIc(int cells) {
  Cloud cloud;
  const double bond = 2.75;
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  const double side = cells * lattice;
  cloud.box = {side, side, side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  int id = 1;
  for (int i = 0; i < cells; i++) {
    for (int j = 0; j < cells; j++) {
      for (int k = 0; k < cells; k++) {
        for (int b = 0; b < 8; b++) {
          molSys::Point<double> p;
          const int f = b % 4;
          const double shift = b < 4 ? 0.0 : 0.25;
          p.x = std::fmod((i + fcc[f][0] + shift) * lattice, side);
          p.y = std::fmod((j + fcc[f][1] + shift) * lattice, side);
          p.z = std::fmod((k + fcc[f][2] + shift) * lattice, side);
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          cloud.idIndexMap[id] = cloud.nop;
          cloud.pts.push_back(p);
          ++cloud.nop;
          ++id;
        }
      }
    }
  }
  return cloud;
}

// A jittered crystal: four bonds on most atoms, so the eight-atom row
// kernel runs on some tiles and declines others.
Cloud jittered(int cells, double amplitude, unsigned seed) {
  Cloud cloud = iceIc(cells);
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(-amplitude, amplitude);
  for (auto &p : cloud.pts) {
    p.x = std::fmod(p.x + u(rng) + cloud.box[0], cloud.box[0]);
    p.y = std::fmod(p.y + u(rng) + cloud.box[1], cloud.box[1]);
    p.z = std::fmod(p.z + u(rng) + cloud.box[2], cloud.box[2]);
  }
  return cloud;
}

Cloud liquid(int n, double side, unsigned seed) {
  Cloud cloud;
  cloud.box = {side, side, side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(0.0, side);
  for (int i = 0; i < n; i++) {
    molSys::Point<double> p;
    p.x = u(rng);
    p.y = u(rng);
    p.z = u(rng);
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    cloud.idIndexMap[p.atomID] = cloud.nop;
    cloud.pts.push_back(p);
    ++cloud.nop;
  }
  return cloud;
}

// Cartesian points from fractional ones in a cell with every tilt set, so
// the full cutoff list takes the triclinic stash scan.
Cloud tiltedLiquid(int n, double side, unsigned seed) {
  const double xy = 1.5;
  const double xz = 0.8;
  const double yz = 1.0;
  Cloud cloud;
  cloud.box = {side + xy + xz, side + yz, side, xy, xz, yz};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int i = 0; i < n; i++) {
    const double fa = u(rng);
    const double fb = u(rng);
    const double fc = u(rng);
    molSys::Point<double> p;
    p.x = fa * side + fb * xy + fc * xz;
    p.y = fb * side + fc * yz;
    p.z = fc * side;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    cloud.idIndexMap[p.atomID] = cloud.nop;
    cloud.pts.push_back(p);
    ++cloud.nop;
  }
  return cloud;
}

struct OneThread {
#if defined(SEAMS_HAS_OPENMP)
  int saved = omp_get_max_threads();
  OneThread() { omp_set_num_threads(1); }
  ~OneThread() { omp_set_num_threads(saved); }
#else
  OneThread() = default;
#endif
  OneThread(const OneThread &) = delete;
  OneThread &operator=(const OneThread &) = delete;
};

void requireClose(const std::vector<double> &kernel, const std::vector<double> &scalar) {
  REQUIRE(kernel.size() == scalar.size());
  for (std::size_t i = 0; i < kernel.size(); i++) {
    const double scale = std::max(1.0, std::fabs(scalar[i]));
    REQUIRE(std::fabs(kernel[i] - scalar[i]) <= 1e-12 * scale);
  }
}

void requireCloseQl(const chill::SteinhardtQl &kernel, const chill::SteinhardtQl &scalar) {
  requireClose(kernel.ql, scalar.ql);
  requireClose(kernel.qlBar, scalar.qlBar);
}

} // namespace

TEST_CASE("the AVX-512 kernels need F, DQ and FMA together", "[ylm][avx512]") {
  const auto &cpu = seams::cpu::features();
#if defined(SEAMS_YLM_AVX512)
  const bool all = cpu.avx512f && cpu.avx512dq && cpu.fma;
  REQUIRE(seams::ylm::avx512Kernels() == all);
  if (!all) {
    WARN("this CPU lacks AVX-512F, AVX-512DQ or FMA; both sides of the comparisons are scalar");
  }
#else
  (void)cpu;
  REQUIRE_FALSE(seams::ylm::avx512Kernels());
#endif
  {
    ScalarOnly scalar;
    REQUIRE_FALSE(seams::ylm::avx512Kernels());
  }
#if defined(SEAMS_YLM_AVX512)
  REQUIRE(seams::ylm::avx512Kernels() == all);
#endif
}

TEST_CASE("ylmBatch and ylmLocalQl match on the scalar path", "[ylm][avx512]") {
  std::mt19937 rng(41017);
  std::uniform_real_distribution<double> u(-3.0, 3.0);
  constexpr int n = 37;
  std::vector<double> dx(n * 8), dy(n * 8), dz(n * 8);
  for (std::size_t i = 0; i < dx.size(); i++) {
    dx[i] = u(rng);
    dy[i] = u(rng);
    dz[i] = u(rng);
  }
  for (int orderL : {3, 4, 6}) {
    const std::size_t len = static_cast<std::size_t>(n) * static_cast<std::size_t>(2 * orderL + 1);
    std::vector<double> re(len), im(len), reScalar(len), imScalar(len);
    seams::ylm::ylmBatch(orderL, dx.data(), dy.data(), dz.data(), n, re.data(), im.data());
    {
      ScalarOnly scalar;
      seams::ylm::ylmBatch(orderL, dx.data(), dy.data(), dz.data(), n, reScalar.data(),
                           imScalar.data());
    }
    requireClose(re, reScalar);
    requireClose(im, imScalar);
  }
  const int nUsed[8] = {0, 1, 4, 7, 12, 25, 36, 37};
  for (int orderL : {4, 6}) {
    const double prefactor = 4.0 * M_PI / (2.0 * orderL + 1.0);
    std::vector<double> ql(8), qlScalar(8);
    REQUIRE(seams::ylm::ylmLocalQl(orderL, dx.data(), dy.data(), dz.data(), n, nUsed,
                                   prefactor, ql.data()));
    {
      ScalarOnly scalar;
      REQUIRE(seams::ylm::ylmLocalQl(orderL, dx.data(), dy.data(), dz.data(), n, nUsed,
                                     prefactor, qlScalar.data()));
    }
    requireClose(ql, qlScalar);
  }
}

TEST_CASE("cutoff q4, q6 and q12 match on the scalar path", "[ylm][avx512]") {
  const Cloud clouds[] = {iceIc(3), jittered(3, 0.15, 52021), liquid(240, 16.0, 91001)};
  for (const Cloud &cloud : clouds) {
    for (int orderL : {4, 6, 12}) {
      for (bool average : {false, true}) {
        const auto kernel = chill::steinhardtQlCutoff(cloud, 3.5, orderL, 1, average);
        chill::SteinhardtQl scalar;
        {
          ScalarOnly off;
          scalar = chill::steinhardtQlCutoff(cloud, 3.5, orderL, 1, average);
        }
        requireCloseQl(kernel, scalar);
      }
    }
    chill::SteinhardtQl k4, k6, s4, s6;
    REQUIRE(chill::steinhardtQ4Q6Cutoff(cloud, 3.5, 1, true, k4, k6));
    {
      ScalarOnly off;
      REQUIRE(chill::steinhardtQ4Q6Cutoff(cloud, 3.5, 1, true, s4, s6));
    }
    requireCloseQl(k4, s4);
    requireCloseQl(k6, s6);
  }
}

TEST_CASE("neighbour-list q12 matches on the scalar path", "[ylm][avx512]") {
  // Twelve and more bonds per atom fill the eight-bond l = 12 kernel and
  // leave a tail.
  auto cloud = liquid(240, 16.0, 77003);
  const auto rows = nneigh::neighListO(5.0, cloud, 1);
  const auto kernel = chill::steinhardtQl(cloud, rows, 12);
  chill::SteinhardtQl scalar;
  {
    ScalarOnly off;
    scalar = chill::steinhardtQl(cloud, rows, 12);
  }
  requireCloseQl(kernel, scalar);
}

TEST_CASE("the AVX2 shell kernels need AVX2 and FMA together", "[symmetry][avx2]") {
  const auto &cpu = seams::cpu::features();
#if defined(SEAMS_SHELL_AVX2)
  const bool both = cpu.avx2 && cpu.fma;
  REQUIRE(symm::avx2Kernels() == both);
  if (!both) {
    WARN("this CPU lacks AVX2 or FMA; both sides of the comparisons are generic");
  }
#else
  (void)cpu;
  REQUIRE_FALSE(symm::avx2Kernels());
#endif
  {
    ScalarOnly generic;
    REQUIRE_FALSE(symm::avx2Kernels());
  }
#if defined(SEAMS_SHELL_AVX2)
  REQUIRE(symm::avx2Kernels() == both);
#endif
}

TEST_CASE("ice labels match on the generic shell path", "[symmetry][avx2]") {
  const Cloud clouds[] = {iceIc(3), jittered(3, 0.15, 52021), jittered(3, 0.35, 52023),
                          liquid(240, 16.0, 91001)};
  for (const Cloud &cloud : clouds) {
    for (bool openFace : {false, true}) {
      std::vector<int> kernel;
      std::vector<int> generic;
      symm::classifyIceCloud(cloud, kernel, openFace);
      {
        ScalarOnly off;
        symm::classifyIceCloud(cloud, generic, openFace);
      }
      REQUIRE(kernel == generic);
    }
  }

  // The thermal TIP4P/2005 slab of test_local_symmetry, under the default
  // gates and the earlier 0.12, 0.94, 1.15.
  std::ifstream in;
  for (const char *path : {"tests/data/tip4p252-oxygens.txt",
                           "../tests/data/tip4p252-oxygens.txt"}) {
    in.open(path);
    if (in) {
      break;
    }
    in.clear();
  }
  REQUIRE(in);
  int n = 0;
  double lx = 0, ly = 0, lz = 0;
  in >> n >> lx >> ly >> lz;
  REQUIRE(n == 4096);
  Eigen::MatrixXd xyz(n, 3);
  for (int i = 0; i < n; i++) {
    in >> xyz(i, 0) >> xyz(i, 1) >> xyz(i, 2);
  }
  REQUIRE(in);
  Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
  box(0, 0) = lx;
  box(1, 1) = ly;
  box(2, 2) = lz;
  struct Gates {
    double rmsd, q, ratio;
    bool openFace;
  };
  const Gates gates[] = {{symm::kShellRmsd, symm::kThermalQ, symm::kThermalBondRatio, false},
                         {symm::kShellRmsd, symm::kThermalQ, symm::kThermalBondRatio, true},
                         {0.12, 0.94, 1.15, false}};
  for (const Gates &g : gates) {
    std::vector<int> kernel;
    std::vector<int> generic;
    symm::classifyIce(xyz, box, kernel, g.rmsd, g.q, g.ratio, g.openFace);
    {
      ScalarOnly off;
      symm::classifyIce(xyz, box, generic, g.rmsd, g.q, g.ratio, g.openFace);
    }
    REQUIRE(kernel == generic);
  }
}

TEST_CASE("diamond signatures match on the generic shell path", "[symmetry][avx2]") {
  // FCC and HCP second shells, nearest-neighbour length sqrt(2), with
  // growing noise so every return of diamondCna comes up.
  const double fcc[12][3] = {
      {-1, -1, 0}, {-1, 1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 0, -1}, {-1, 0, 1},
      {1, 0, -1},  {1, 0, 1},  {0, -1, -1}, {0, -1, 1}, {0, 1, -1}, {0, 1, 1},
  };
  const double r3 = std::sqrt(3.0);
  const double h = 2.0 / std::sqrt(3.0);
  const double hcp[12][3] = {
      {std::sqrt(2.0), 0, 0},       {-std::sqrt(2.0), 0, 0},
      {std::sqrt(0.5), r3 / std::sqrt(2.0), 0},  {-std::sqrt(0.5), r3 / std::sqrt(2.0), 0},
      {std::sqrt(0.5), -r3 / std::sqrt(2.0), 0}, {-std::sqrt(0.5), -r3 / std::sqrt(2.0), 0},
      {0, std::sqrt(2.0 / 3.0), h},  {std::sqrt(0.5), -std::sqrt(1.0 / 6.0), h},
      {-std::sqrt(0.5), -std::sqrt(1.0 / 6.0), h}, {0, std::sqrt(2.0 / 3.0), -h},
      {std::sqrt(0.5), -std::sqrt(1.0 / 6.0), -h}, {-std::sqrt(0.5), -std::sqrt(1.0 / 6.0), -h},
  };
  std::mt19937 rng(60311);
  for (const auto *shell : {fcc, hcp}) {
    for (double amplitude : {0.0, 0.05, 0.12, 0.2, 0.3}) {
      for (int trial = 0; trial < 64; trial++) {
        std::normal_distribution<double> noise(0.0, amplitude > 0.0 ? amplitude : 1.0);
        double sx[12], sy[12], sz[12];
        for (int i = 0; i < 12; i++) {
          sx[i] = shell[i][0] + (amplitude > 0.0 ? noise(rng) : 0.0);
          sy[i] = shell[i][1] + (amplitude > 0.0 ? noise(rng) : 0.0);
          sz[i] = shell[i][2] + (amplitude > 0.0 ? noise(rng) : 0.0);
        }
        for (int ambiguous : {9, 12}) {
          int k421 = -1, k422 = -1, kGood = -1;
          int g421 = -1, g422 = -1, gGood = -1;
          const int kernel = symm::diamondCna(sx, sy, sz, ambiguous, &k421, &k422, &kGood);
          int generic = 0;
          {
            ScalarOnly off;
            generic = symm::diamondCna(sx, sy, sz, ambiguous, &g421, &g422, &gGood);
          }
          REQUIRE(kernel == generic);
          REQUIRE(k421 == g421);
          REQUIRE(k422 == g422);
          REQUIRE(kGood == gGood);
        }
      }
    }
  }
}

TEST_CASE("triclinic cutoff lists match on the generic stash scan", "[neighbours][avx2]") {
  OneThread one;
  for (unsigned seed : {17011u, 17013u}) {
    for (int nAtoms : {240, 900}) {
      const Cloud cloud = tiltedLiquid(nAtoms, 16.0, seed);
      for (double rc : {3.2, 3.5}) {
        for (bool branchy : {false, true}) {
          std::vector<int> kOffsets, kCols, gOffsets, gCols;
          REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, rc, kOffsets, kCols, false, nullptr,
                                       nullptr, nullptr, branchy));
          {
            ScalarOnly off;
            REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, rc, gOffsets, gCols, false, nullptr,
                                         nullptr, nullptr, branchy));
          }
          REQUIRE_FALSE(kCols.empty());
          REQUIRE(kOffsets == gOffsets);
          REQUIRE(kCols == gCols);
        }
      }
    }
  }
}
