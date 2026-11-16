// Copyright (c) d-SEAMS developers
#include <catch2/catch_test_macros.hpp>

#include <bop.hpp>
#include <cpu_features.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ylm_batch.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

// Every call into the AVX-512 kernel archive goes through
// seams::ylm::avx512Kernels. With the kernels turned off the scalar path
// runs on any CPU, and each value it gives must match the kernel's.

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

struct ScalarOnly {
  ScalarOnly() { seams::ylm::allowAvx512Kernels(false); }
  ~ScalarOnly() { seams::ylm::allowAvx512Kernels(true); }
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
