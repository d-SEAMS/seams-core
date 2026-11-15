// The cutoff neighbour list, local Steinhardt q_l and the ring census at
// one, two, four and eight OpenMP threads. Every width must return the
// one-thread arrays exactly: the same CSR entries in the same order, and
// the same q_l bits.

#include <catch2/catch_test_macros.hpp>

#include <bop.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_c_api.h>

#include "ice_ic.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

constexpr int kWidths[] = {1, 2, 4, 8};

struct Width {
  int saved = 1;
  explicit Width(int n) {
#ifdef SEAMS_HAS_OPENMP
    saved = omp_get_max_threads();
    omp_set_dynamic(0);
    omp_set_num_threads(n);
#else
    (void)n;
#endif
  }
  ~Width() {
#ifdef SEAMS_HAS_OPENMP
    omp_set_num_threads(saved);
#endif
  }
  Width(const Width &) = delete;
  Width &operator=(const Width &) = delete;
};

struct Lcg {
  std::uint64_t s;
  explicit Lcg(std::uint64_t seed) : s(seed) {}
  double next() {
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(s >> 11) * (1.0 / 9007199254740992.0);
  }
};

// Uniform points in an orthorhombic box. A tenth of them sit one box
// length outside, so the list has to wrap them.
Cloud randomCloud(int n, double lx, double ly, double lz, std::uint64_t seed) {
  Cloud cloud;
  Lcg rng(seed);
  cloud.nop = n;
  cloud.currentFrame = 1;
  cloud.box = {lx, ly, lz};
  cloud.boxLow = {-1.5, 0.25, 0.0};
  cloud.pts.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    auto &p = cloud.pts[static_cast<std::size_t>(i)];
    p.x = cloud.boxLow[0] + rng.next() * lx;
    p.y = cloud.boxLow[1] + rng.next() * ly;
    p.z = cloud.boxLow[2] + rng.next() * lz;
    if (i % 10 == 3) {
      p.x += lx;
      p.z -= lz;
    }
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
  }
  return cloud;
}

Cloud jittered(Cloud cloud, double amplitude, std::uint64_t seed) {
  Lcg rng(seed);
  for (auto &p : cloud.pts) {
    p.x += (rng.next() - 0.5) * amplitude;
    p.y += (rng.next() - 0.5) * amplitude;
    p.z += (rng.next() - 0.5) * amplitude;
  }
  return cloud;
}

struct Csr {
  bool ok = false;
  std::vector<int> offsets;
  std::vector<int> cols;
};

Csr csrAt(int width, const Cloud &cloud, const std::vector<int> *subset, double rc,
          bool half = false, bool branchy = false) {
  const Width team(width);
  Csr out;
  out.ok = nneigh::orthoFullCsr(cloud, subset, rc, out.offsets, out.cols, half, nullptr,
                                nullptr, nullptr, branchy);
  return out;
}

// Every pair in the nearest image, rows ascending. Valid while the box is
// at least twice the cutoff on every axis.
Csr bruteCsr(const Cloud &cloud, const std::vector<int> *subset, double rc, bool half) {
  const int n = cloud.nop;
  std::vector<char> in(static_cast<std::size_t>(n), subset == nullptr ? 1 : 0);
  if (subset != nullptr) {
    for (const int id : *subset) {
      in[static_cast<std::size_t>(id)] = 1;
    }
  }
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
  const double l[3] = {cloud.box[0], cloud.box[1], cloud.box[2]};
  for (int i = 0; i < n; ++i) {
    if (!in[static_cast<std::size_t>(i)]) {
      continue;
    }
    const auto &a = cloud.pts[static_cast<std::size_t>(i)];
    for (int j = 0; j < n; ++j) {
      if (j == i || !in[static_cast<std::size_t>(j)] || (half && j < i)) {
        continue;
      }
      const auto &b = cloud.pts[static_cast<std::size_t>(j)];
      double d[3] = {a.x - b.x, a.y - b.y, a.z - b.z};
      double d2 = 0.0;
      for (int k = 0; k < 3; ++k) {
        d[k] -= l[k] * std::round(d[k] / l[k]);
        d2 += d[k] * d[k];
      }
      if (d2 < rc * rc) {
        rows[static_cast<std::size_t>(i)].push_back(j);
      }
    }
  }
  Csr out;
  out.ok = true;
  out.offsets.assign(static_cast<std::size_t>(n) + 1, 0);
  for (int i = 0; i < n; ++i) {
    out.offsets[static_cast<std::size_t>(i) + 1] =
        out.offsets[static_cast<std::size_t>(i)] +
        static_cast<int>(rows[static_cast<std::size_t>(i)].size());
    out.cols.insert(out.cols.end(), rows[static_cast<std::size_t>(i)].begin(),
                    rows[static_cast<std::size_t>(i)].end());
  }
  return out;
}

void requireSameCsr(const Csr &got, const Csr &want) {
  REQUIRE(got.ok == want.ok);
  REQUIRE(got.offsets == want.offsets);
  REQUIRE(got.cols == want.cols);
}

// Every width against one thread, full and half, with and without the
// merged z runs.
void requireWidthInvariant(const Cloud &cloud, const std::vector<int> *subset, double rc) {
  for (const bool half : {false, true}) {
    for (const bool branchy : {false, true}) {
      const Csr serial = csrAt(1, cloud, subset, rc, half, branchy);
      REQUIRE(serial.ok);
      for (const int width : kWidths) {
        INFO("width " << width << " half " << half << " branchy " << branchy);
        requireSameCsr(csrAt(width, cloud, subset, rc, half, branchy), serial);
      }
    }
  }
}

std::vector<int> everyThird(int n) {
  std::vector<int> subset;
  for (int i = 1; i < n; i += 3) {
    subset.push_back(i);
  }
  return subset;
}

chill::SteinhardtQl qlAt(int width, const Cloud &cloud, double rc, int l, int typeI,
                         bool average) {
  const Width team(width);
  return chill::steinhardtQlCutoff(cloud, rc, l, typeI, average);
}

void requireSameBits(const std::vector<double> &got, const std::vector<double> &want) {
  REQUIRE(got.size() == want.size());
  for (std::size_t i = 0; i < want.size(); ++i) {
    INFO("atom " << i);
    REQUIRE(std::memcmp(&got[i], &want[i], sizeof(double)) == 0);
  }
}

std::vector<double> flatXyz(const Cloud &cloud) {
  std::vector<double> xyz;
  xyz.reserve(static_cast<std::size_t>(cloud.nop) * 3);
  for (const auto &p : cloud.pts) {
    xyz.push_back(p.x);
    xyz.push_back(p.y);
    xyz.push_back(p.z);
  }
  return xyz;
}

std::vector<double> box9Of(const Cloud &cloud) {
  return {cloud.boxLow[0], cloud.boxLow[1], cloud.boxLow[2], cloud.box[0], cloud.box[1],
          cloud.box[2], 0.0,             0.0,             0.0};
}

} // namespace

TEST_CASE("the cutoff CSR of ice Ic is the same at every width", "[thread_invariance]") {
  const Cloud ice = iceIc(13824);
  REQUIRE(ice.nop == 13824);
  requireWidthInvariant(ice, nullptr, 3.5);
  const Cloud warm = jittered(ice, 0.6, 7);
  requireWidthInvariant(warm, nullptr, 3.5);
  const auto subset = everyThird(warm.nop);
  requireWidthInvariant(warm, &subset, 3.5);
}

TEST_CASE("the cutoff CSR of a random liquid is the same at every width and matches "
          "brute force",
          "[thread_invariance]") {
  const Cloud liquid = randomCloud(6000, 56.0, 57.0, 58.0, 11);
  requireWidthInvariant(liquid, nullptr, 3.5);
  const auto subset = everyThird(liquid.nop);
  requireWidthInvariant(liquid, &subset, 3.5);
  for (const bool half : {false, true}) {
    requireSameCsr(csrAt(8, liquid, nullptr, 3.5, half), bruteCsr(liquid, nullptr, 3.5, half));
    requireSameCsr(csrAt(8, liquid, &subset, 3.5, half), bruteCsr(liquid, &subset, 3.5, half));
  }
}

TEST_CASE("two cells along one axis give the same CSR at every width", "[thread_invariance]") {
  // Two cells reach the same neighbour column through both images.
  const Cloud flatZ = randomCloud(4096, 60.0, 61.0, 8.0, 21);
  const Cloud flatY = randomCloud(4096, 60.0, 8.0, 61.0, 22);
  const Cloud flatX = randomCloud(4096, 8.0, 60.0, 61.0, 23);
  for (const Cloud *cloud : {&flatZ, &flatY, &flatX}) {
    requireWidthInvariant(*cloud, nullptr, 3.5);
    requireSameCsr(csrAt(8, *cloud, nullptr, 3.5), bruteCsr(*cloud, nullptr, 3.5, false));
  }
}

TEST_CASE("rows longer than the slot stash give the same CSR at every width",
          "[thread_invariance]") {
  // About 64 neighbours inside 7.9 A of ice Ic.
  const Cloud ice = jittered(iceIc(13824), 0.4, 5);
  requireWidthInvariant(ice, nullptr, 7.9);
  const Cloud dense = randomCloud(3000, 20.0, 21.0, 22.0, 31);
  requireWidthInvariant(dense, nullptr, 6.5);
  requireSameCsr(csrAt(8, dense, nullptr, 6.5), bruteCsr(dense, nullptr, 6.5, false));
}

TEST_CASE("a small frame and a tilted box give the same CSR at every width",
          "[thread_invariance]") {
  const Cloud small = randomCloud(300, 12.0, 13.0, 14.0, 41);
  requireWidthInvariant(small, nullptr, 3.5);
  requireSameCsr(csrAt(4, small, nullptr, 3.5), bruteCsr(small, nullptr, 3.5, false));
  Cloud tilted = jittered(iceIc(4096), 0.3, 9);
  tilted.box = {tilted.box[0], tilted.box[1], tilted.box[2], 1.5, 0.0, 0.0};
  const Csr serial = csrAt(1, tilted, nullptr, 3.5);
  REQUIRE(serial.ok);
  for (const int width : kWidths) {
    INFO("width " << width);
    requireSameCsr(csrAt(width, tilted, nullptr, 3.5), serial);
  }
}

TEST_CASE("local q4 and q6 keep every bit at every width", "[thread_invariance]") {
  const Cloud warm = jittered(iceIc(13824), 0.6, 3);
  const Cloud liquid = randomCloud(6000, 56.0, 57.0, 58.0, 13);
  for (const Cloud *cloud : {&warm, &liquid}) {
    for (const int l : {3, 4, 6}) {
      const auto serial = qlAt(1, *cloud, 3.5, l, 1, false);
      REQUIRE(serial.ql.size() == static_cast<std::size_t>(cloud->nop));
      for (const int width : kWidths) {
        INFO("l " << l << " width " << width);
        const auto got = qlAt(width, *cloud, 3.5, l, 1, false);
        requireSameBits(got.ql, serial.ql);
        requireSameBits(got.qlBar, serial.qlBar);
      }
    }
  }
}

TEST_CASE("averaged q6 and a mixed-type q4 keep every bit at every width",
          "[thread_invariance]") {
  Cloud mixed = jittered(iceIc(13824), 0.6, 4);
  for (int i = 0; i < mixed.nop; i += 5) {
    mixed.pts[static_cast<std::size_t>(i)].type = 2;
  }
  const auto serial6 = qlAt(1, mixed, 3.5, 6, 1, true);
  const auto serial4 = qlAt(1, mixed, 3.5, 4, 1, false);
  for (const int width : kWidths) {
    INFO("width " << width);
    const auto got6 = qlAt(width, mixed, 3.5, 6, 1, true);
    requireSameBits(got6.ql, serial6.ql);
    requireSameBits(got6.qlBar, serial6.qlBar);
    const auto got4 = qlAt(width, mixed, 3.5, 4, 1, false);
    requireSameBits(got4.ql, serial4.ql);
  }
}

TEST_CASE("the C API CSR, q_l and ring census are the same at every width",
          "[thread_invariance][c_api]") {
  const Cloud warm = jittered(iceIc(13824), 0.5, 8);
  const auto xyz = flatXyz(warm);
  const auto box9 = box9Of(warm);
  const int n = warm.nop;

  const Csr direct = csrAt(1, warm, nullptr, 3.5);
  std::vector<int> wantCounts(7, 0);
  for (const auto &ring : primitive::ringNetwork(direct.offsets, direct.cols, 6)) {
    if (ring.size() < wantCounts.size()) {
      wantCounts[ring.size()] += 1;
    }
  }
  REQUIRE(wantCounts[6] > 0);

  std::vector<double> serialQ4;
  std::vector<double> serialQ6;
  for (const int width : kWidths) {
    INFO("width " << width);
    const Width team(width);
    std::vector<int> offsets(static_cast<std::size_t>(n) + 1, -1);
    std::vector<int> cols(direct.cols.size() + 16, -1);
    int nCols = -1;
    REQUIRE(seams_cutoff_csr(xyz.data(), n, box9.data(), 3.5, offsets.data(), cols.data(),
                             static_cast<int>(cols.size()), &nCols) == 0);
    REQUIRE(nCols == static_cast<int>(direct.cols.size()));
    cols.resize(static_cast<std::size_t>(nCols));
    REQUIRE(offsets == direct.offsets);
    REQUIRE(cols == direct.cols);

    std::vector<double> q4(static_cast<std::size_t>(n), -1.0);
    std::vector<double> q6(static_cast<std::size_t>(n), -1.0);
    REQUIRE(seams_ql_cutoff(xyz.data(), n, box9.data(), 3.5, 4, q4.data()) == 0);
    REQUIRE(seams_ql_cutoff(xyz.data(), n, box9.data(), 3.5, 6, q6.data()) == 0);
    if (width == 1) {
      serialQ4 = q4;
      serialQ6 = q6;
      requireSameBits(q4, qlAt(1, warm, 3.5, 4, 1, false).ql);
    }
    requireSameBits(q4, serialQ4);
    requireSameBits(q6, serialQ6);

    std::vector<int> counts(7, -1);
    REQUIRE(seams_primitive_rings(xyz.data(), n, box9.data(), 3.5, 6, counts.data(),
                                  static_cast<int>(counts.size())) == 0);
    REQUIRE(counts == wantCounts);
    std::vector<int> sized(7, -1);
    primitive::ringSizeCounts(direct.offsets, direct.cols, 6, sized.data(),
                              static_cast<int>(sized.size()));
    REQUIRE(sized == wantCounts);
  }
}

TEST_CASE("the C API reuses its frame buffers across frames of different sizes",
          "[thread_invariance][c_api]") {
  const Cloud big = jittered(iceIc(13824), 0.5, 12);
  const Cloud small = randomCloud(2500, 30.0, 31.0, 32.0, 17);
  const Cloud other = jittered(iceIc(4096), 0.5, 19);
  for (const Cloud *cloud : {&big, &small, &other, &big}) {
    const auto xyz = flatXyz(*cloud);
    const auto box9 = box9Of(*cloud);
    const int n = cloud->nop;
    const Cloud &fresh = *cloud;
    const Csr want = csrAt(8, fresh, nullptr, 3.5);
    std::vector<int> offsets(static_cast<std::size_t>(n) + 1, -1);
    int nCols = -1;
    REQUIRE(seams_cutoff_csr(xyz.data(), n, box9.data(), 3.5, offsets.data(), nullptr, 0,
                             &nCols) == 0);
    REQUIRE(nCols == static_cast<int>(want.cols.size()));
    std::vector<int> cols(static_cast<std::size_t>(nCols), -1);
    REQUIRE(seams_cutoff_csr(xyz.data(), n, box9.data(), 3.5, offsets.data(), cols.data(),
                             nCols, &nCols) == 0);
    REQUIRE(offsets == want.offsets);
    REQUIRE(cols == want.cols);
    std::vector<double> q6(static_cast<std::size_t>(n), -1.0);
    REQUIRE(seams_ql_cutoff(xyz.data(), n, box9.data(), 3.5, 6, q6.data()) == 0);
    requireSameBits(q6, qlAt(8, fresh, 3.5, 6, 1, false).ql);
  }
}
