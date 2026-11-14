// The cutoff neighbour list, local Steinhardt q_l, the ring census, CHILL
// and CHILL+, the shell match, rings with their symmetry codes, cages,
// prisms, the monolayer census and the 2D RDF at one, two, four and eight
// OpenMP threads. Every width must return the one-thread arrays exactly:
// the same CSR entries in the same order, the same labels, rings and
// counts, and the same floating-point bits.

#include <catch2/catch_test_macros.hpp>

#include <bop.hpp>
#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <rdf2d.hpp>
#include <seams.h>
#include <seams_c_api.h>
#include <topo_bulk.hpp>
#include <topo_one_dim.hpp>

#include "ice_ic.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

constexpr int kWidths[] = {1, 2, 4, 8};

// The library narrows a call's team when the machine looks contended, and
// meson runs these tests beside others. Every width below must be the
// width asked for.
const bool kFixedWidth = setenv("SEAMS_FAIR_WIDTH", "0", 1) == 0;

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

// Honeycomb sheets stacked straight up, 2.75 A bonds in and between the
// sheets. Each hexagon and each square between two sheets is a prism face
// or base, and no other pair is inside 3.5 A.
Cloud stackedHoneycomb(int nx, int ny, int nz) {
  const double bond = 2.75;
  const double a = bond * std::sqrt(3.0);
  const double basis[4][2] = {
      {0.0, 0.0}, {0.5 * a, 0.5 * bond}, {0.5 * a, 1.5 * bond}, {0.0, 2.0 * bond}};
  Cloud cloud;
  cloud.currentFrame = 1;
  cloud.box = {nx * a, ny * 3.0 * bond, nz * bond};
  cloud.boxLow = {0.0, 0.0, 0.0};
  for (int iz = 0; iz < nz; ++iz) {
    for (int ix = 0; ix < nx; ++ix) {
      for (int iy = 0; iy < ny; ++iy) {
        for (const auto &site : basis) {
          molSys::Point<double> point;
          point.x = ix * a + site[0];
          point.y = iy * 3.0 * bond + site[1];
          point.z = iz * bond;
          point.type = 1;
          point.atomID = static_cast<int>(cloud.pts.size()) + 1;
          point.molID = point.atomID;
          cloud.idIndexMap[point.atomID] = static_cast<int>(cloud.pts.size());
          cloud.pts.push_back(point);
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  return cloud;
}

std::vector<std::vector<int>> indexRows(const Cloud &cloud, double rc) {
  return nneigh::neighbourListByIndex(cloud, nneigh::neighListO(rc, cloud, 1));
}

std::vector<std::vector<int>> ofSize(const std::vector<std::vector<int>> &rings,
                                     std::size_t size) {
  std::vector<std::vector<int>> out;
  for (const auto &ring : rings) {
    if (ring.size() == size) {
      out.push_back(ring);
    }
  }
  return out;
}

// One frame handle of the stable C API.
struct System {
  seams_system *handle = seams_system_create();
  explicit System(const Cloud &cloud, const std::vector<int> *types = nullptr) {
    set(cloud, types);
  }
  ~System() { seams_system_destroy(handle); }
  System(const System &) = delete;
  System &operator=(const System &) = delete;
  void set(const Cloud &cloud, const std::vector<int> *types = nullptr) {
    const auto xyz = flatXyz(cloud);
    REQUIRE(seams_system_set_ortho(handle, xyz.data(), cloud.nop, cloud.boxLow[0],
                                   cloud.boxLow[1], cloud.boxLow[2], cloud.box[0],
                                   cloud.box[1], cloud.box[2],
                                   types == nullptr ? nullptr : types->data(),
                                   nullptr) == SEAMS_OK);
  }
};

struct Labels {
  std::vector<int> labels;
  std::vector<int> counts;
  bool operator==(const Labels &) const = default;
};

enum class Chill { published, ovito, plusPublished, plusOvito };

Labels chillOf(seams_system *sys, int n, Chill kind, int type) {
  Labels out;
  out.labels.assign(static_cast<std::size_t>(n), -7);
  out.counts.assign(9, -1);
  const int rule = (kind == Chill::ovito || kind == Chill::plusOvito) ? SEAMS_CHILL_OVITO
                                                                      : SEAMS_CHILL_PUBLISHED;
  if (kind == Chill::plusPublished || kind == Chill::plusOvito) {
    REQUIRE(seams_chill_plus(sys, 3.5, type, SEAMS_GRAPH_CUTOFF, rule, out.labels.data(),
                             out.counts.data()) == SEAMS_OK);
  } else {
    REQUIRE(seams_chill(sys, 3.5, type, SEAMS_GRAPH_CUTOFF, rule, out.labels.data(),
                        out.labels.size(), out.counts.data()) == SEAMS_OK);
  }
  return out;
}

struct RingSet {
  std::vector<int> sizes;
  std::vector<int> symmetry;
  std::vector<int> vertices;
  std::vector<int> histogram;
  bool operator==(const RingSet &) const = default;
};

RingSet ringsOf(seams_system *sys, int type) {
  seams_rings *rings = nullptr;
  REQUIRE(seams_rings_find(sys, 3.5, type, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK);
  const int count = seams_rings_count(rings);
  RingSet out;
  out.sizes.assign(static_cast<std::size_t>(count), -1);
  out.symmetry.assign(static_cast<std::size_t>(count), -9);
  REQUIRE(seams_rings_labels(rings, out.sizes.data(), out.symmetry.data()) == SEAMS_OK);
  out.histogram.assign(7, -1);
  REQUIRE(seams_rings_histogram(rings, out.histogram.data(), 7) == SEAMS_OK);
  int atoms[16];
  for (int i = 0; i < count; ++i) {
    std::size_t got = 0;
    REQUIRE(seams_rings_vertices(rings, i, atoms, 16, &got) == SEAMS_OK);
    out.vertices.insert(out.vertices.end(), atoms, atoms + got);
  }
  seams_rings_destroy(rings);
  return out;
}

struct CageSet {
  std::vector<int> counts;
  std::vector<int> perAtom;
  bool operator==(const CageSet &) const = default;
};

CageSet cagesOf(seams_system *sys, int n, int type) {
  seams_cage_set *cages = nullptr;
  REQUIRE(seams_cages_find(sys, 3.5, type, SEAMS_GRAPH_CUTOFF, 4, 0, &cages) == SEAMS_OK);
  CageSet out;
  out.counts.assign(4, -1);
  REQUIRE(seams_cages_counts(cages, &out.counts[0], &out.counts[1], &out.counts[2],
                             &out.counts[3]) == SEAMS_OK);
  out.perAtom.assign(static_cast<std::size_t>(n), -1);
  REQUIRE(seams_cages_labels(cages, out.perAtom.data(), out.perAtom.size()) == SEAMS_OK);
  seams_cages_destroy(cages);
  return out;
}

std::vector<int> prismsOf(seams_system *sys, int type) {
  std::vector<int> perfect(7, -1);
  std::vector<int> imperfect(7, -1);
  REQUIRE(seams_prisms(sys, 3.5, type, perfect.data(), imperfect.data(), 7) == SEAMS_OK);
  perfect.insert(perfect.end(), imperfect.begin(), imperfect.end());
  return perfect;
}

std::vector<int> monolayerOf(seams_system *sys, const Cloud &cloud, int type) {
  std::vector<int> counts(7, -1);
  REQUIRE(seams_monolayer(sys, 3.5, type, 6, cloud.box[0] * cloud.box[1], counts.data(),
                          counts.size()) == SEAMS_OK);
  return counts;
}

int distinctLabels(const std::vector<int> &labels) {
  std::vector<int> seen(labels);
  std::sort(seen.begin(), seen.end());
  return static_cast<int>(std::unique(seen.begin(), seen.end()) - seen.begin());
}

// Every type-2 atom is a fifth of the frame, so a species filter drops
// atoms from the middle of the arrays.
std::vector<int> everyFifthType2(int n) {
  std::vector<int> types(static_cast<std::size_t>(n), 1);
  for (int i = 2; i < n; i += 5) {
    types[static_cast<std::size_t>(i)] = 2;
  }
  return types;
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

TEST_CASE("CHILL and CHILL+ labels and counts are the same at every width",
          "[thread_invariance][c_api]") {
  REQUIRE(kFixedWidth);
  const Cloud warm = jittered(iceIc(32768), 0.8, 21);
  const auto types = everyFifthType2(warm.nop);
  System plain(warm);
  System mixed(warm, &types);
  for (const Chill kind : {Chill::published, Chill::ovito, Chill::plusPublished,
                           Chill::plusOvito}) {
    for (const auto &[sys, type] : {std::pair{plain.handle, 0}, std::pair{mixed.handle, 1}}) {
      Labels serial;
      {
        const Width team(1);
        serial = chillOf(sys, warm.nop, kind, type);
      }
      REQUIRE(distinctLabels(serial.labels) >= 2);
      for (const int width : kWidths) {
        INFO("kind " << static_cast<int>(kind) << " type " << type << " width " << width);
        const Width team(width);
        REQUIRE(chillOf(sys, warm.nop, kind, type) == serial);
      }
    }
  }
}

TEST_CASE("the shell match keeps its labels at every width", "[thread_invariance][c_api]") {
  REQUIRE(kFixedWidth);
  const Cloud warm = jittered(iceIc(32768), 0.8, 22);
  const auto xyz = flatXyz(warm);
  const double lattice9[9] = {warm.box[0], 0.0, 0.0, 0.0, warm.box[1], 0.0,
                              0.0,         0.0, 0.0 + warm.box[2]};
  System sys(warm);
  std::vector<int> serialApi;
  std::vector<int> serialFlat;
  for (const int width : kWidths) {
    INFO("width " << width);
    const Width team(width);
    std::vector<int> api(static_cast<std::size_t>(warm.nop), -9);
    REQUIRE(seams_system_symmetry(sys.handle, 0.0, api.data(), api.size()) == SEAMS_OK);
    std::vector<int> flat(static_cast<std::size_t>(warm.nop), -9);
    REQUIRE(seams_symmetry_classify(xyz.data(), warm.nop, lattice9, flat.data()) == 0);
    if (width == 1) {
      serialApi = api;
      serialFlat = flat;
      REQUIRE(distinctLabels(api) >= 2);
    }
    REQUIRE(api == serialApi);
    REQUIRE(flat == serialFlat);
  }
}

TEST_CASE("rings, their symmetry codes and the cage planes are the same at every width",
          "[thread_invariance][c_api]") {
  REQUIRE(kFixedWidth);
  const Cloud warm = jittered(iceIc(32768), 0.8, 23);
  System sys(warm);
  RingSet serialSet;
  {
    const Width team(1);
    serialSet = ringsOf(sys.handle, 0);
  }
  REQUIRE(serialSet.histogram[6] > 0);
  REQUIRE(distinctLabels(serialSet.symmetry) >= 2);

  const auto idx = indexRows(warm, 3.5);
  std::vector<std::vector<int>> serialRings;
  {
    const Width team(1);
    serialRings = primitive::ringNetwork(idx, 6);
  }
  const auto six = ofSize(serialRings, 6);
  REQUIRE(six.size() > 2048);
  const auto wantAff = ring::cageAffiliation(six, idx);
  const auto wantPlanes = ring::stackingPlanes(six, idx);
  const auto contractible = ring::contractibleRings(serialRings, warm);
  REQUIRE(ring::contractibleRings(std::vector<std::vector<int>>(serialRings), warm) ==
          contractible);
  for (const int width : kWidths) {
    INFO("width " << width);
    const Width team(width);
    REQUIRE(ringsOf(sys.handle, 0) == serialSet);
    REQUIRE(primitive::ringNetwork(idx, 6) == serialRings);
    const auto both = ring::cagePlanes(six, idx);
    REQUIRE(both.affiliation.hc == wantAff.hc);
    REQUIRE(both.affiliation.ddc == wantAff.ddc);
    REQUIRE(both.planes.basal == wantPlanes.basal);
    REQUIRE(both.planes.equatorial == wantPlanes.equatorial);
    REQUIRE(ring::cageAffiliation(six, idx).hc == wantAff.hc);
    REQUIRE(ring::stackingPlanes(six, idx).equatorial == wantPlanes.equatorial);
  }
}

TEST_CASE("cages and their per-atom codes are the same at every width",
          "[thread_invariance][c_api]") {
  REQUIRE(kFixedWidth);
  const Cloud cold = jittered(iceIc(32768), 0.3, 24);
  const Cloud warm = jittered(iceIc(32768), 0.8, 25);
  for (const Cloud *cloud : {&cold, &warm}) {
    System sys(*cloud);
    CageSet serial;
    {
      const Width team(1);
      serial = cagesOf(sys.handle, cloud->nop, 0);
    }
    REQUIRE(serial.counts[1] + serial.counts[3] > 0);
    for (const int width : kWidths) {
      INFO("width " << width);
      const Width team(width);
      REQUIRE(cagesOf(sys.handle, cloud->nop, 0) == serial);
    }
  }
}

TEST_CASE("prisms are the same at every width", "[thread_invariance][c_api]") {
  REQUIRE(kFixedWidth);
  for (const double amplitude : {0.0, 0.2}) {
    INFO("jitter " << amplitude);
    Cloud stack = jittered(stackedHoneycomb(12, 10, 16), amplitude, 26);
    const auto idx = indexRows(stack, 3.5);
    std::vector<std::vector<int>> rings;
    {
      const Width team(1);
      rings = primitive::ringNetwork(idx, 6);
    }
    System sys(stack);
    std::vector<int> serialApi;
    {
      const Width team(1);
      serialApi = prismsOf(sys.handle, 0);
    }
    for (const std::size_t size : {std::size_t{4}, std::size_t{6}}) {
      const auto group = ofSize(rings, size);
      REQUIRE(group.size() >= 2048);
      std::vector<ring::strucType> serialKind;
      std::vector<int> serialList;
      int serialPerfect = -1;
      for (const int width : kWidths) {
        INFO("size " << size << " width " << width);
        const Width team(width);
        std::vector<ring::strucType> kind(group.size());
        std::vector<double> rmsd(static_cast<std::size_t>(stack.nop), -1.0);
        int perfect = 0;
        int imperfect = 0;
        auto list =
            ring::findPrisms(group, kind, perfect, imperfect, idx, stack, rmsd, false);
        if (width == 1) {
          serialKind = kind;
          serialList = list;
          serialPerfect = perfect;
          // Squares stacked on squares share a bond, so only the
          // hexagons pair into prisms.
          REQUIRE((size == 4 || perfect > 0));
        }
        REQUIRE(perfect == serialPerfect);
        REQUIRE(list == serialList);
        REQUIRE(kind == serialKind);
      }
    }
    for (const int width : kWidths) {
      INFO("width " << width);
      const Width team(width);
      REQUIRE(prismsOf(sys.handle, 0) == serialApi);
    }
  }
}

TEST_CASE("the monolayer census, the 2D RDF and k-nearest rows are the same at every width",
          "[thread_invariance][c_api]") {
  REQUIRE(kFixedWidth);
  const Cloud warm = jittered(iceIc(32768), 0.8, 27);
  System sys(warm);
  std::vector<int> serialCensus;
  std::vector<int> serialHistogram;
  std::vector<double> serialLengths;
  std::vector<std::vector<int>> serialMutual;
  std::vector<std::vector<int>> serialAll;
  for (const int width : kWidths) {
    INFO("width " << width);
    const Width team(width);
    const auto census = monolayerOf(sys.handle, warm, 0);
    const auto histogram = rdf2::sampleRDF_AA(warm, 8.0, 0.05, 160);
    const auto lengths = rdf2::getSystemLengths(warm);
    const auto mutual = nneigh::kNearestNeighbourList(warm, 4, 5.0, 1, true);
    const auto all = nneigh::kNearestNeighbourList(warm, 4, 5.0, 1, false);
    if (width == 1) {
      serialCensus = census;
      serialHistogram = histogram;
      serialLengths = lengths;
      serialMutual = mutual;
      serialAll = all;
      REQUIRE(census[6] > 0);
      REQUIRE(std::count_if(histogram.begin(), histogram.end(), [](int v) { return v > 0; }) >
              10);
    }
    REQUIRE(census == serialCensus);
    REQUIRE(histogram == serialHistogram);
    requireSameBits(lengths, serialLengths);
    REQUIRE(mutual == serialMutual);
    REQUIRE(all == serialAll);
  }
}

TEST_CASE("one frame handle gives each kernel the frame a fresh handle would",
          "[thread_invariance][c_api]") {
  // The handle keeps one filtered copy of the frame between calls. A new
  // species filter, a new frame, or a kernel that does not use the copy
  // in between must not leak into the next call.
  REQUIRE(kFixedWidth);
  const Cloud warm = jittered(iceIc(32768), 0.8, 28);
  const Cloud small = jittered(iceIc(13824), 0.8, 29);
  const auto types = everyFifthType2(warm.nop);
  const auto smallTypes = everyFifthType2(small.nop);
  const Width team(8);
  System reused(warm, &types);
  const auto fresh = [](const Cloud &cloud, const std::vector<int> &kinds) {
    return std::make_unique<System>(cloud, &kinds);
  };
  REQUIRE(chillOf(reused.handle, warm.nop, Chill::plusPublished, 1) ==
          chillOf(fresh(warm, types)->handle, warm.nop, Chill::plusPublished, 1));
  REQUIRE(cagesOf(reused.handle, warm.nop, 0) == cagesOf(fresh(warm, types)->handle, warm.nop, 0));
  REQUIRE(chillOf(reused.handle, warm.nop, Chill::published, 2) ==
          chillOf(fresh(warm, types)->handle, warm.nop, Chill::published, 2));
  std::vector<int> labels(static_cast<std::size_t>(warm.nop), -9);
  REQUIRE(seams_system_symmetry(reused.handle, 0.0, labels.data(), labels.size()) == SEAMS_OK);
  REQUIRE(ringsOf(reused.handle, 1) == ringsOf(fresh(warm, types)->handle, 1));
  REQUIRE(prismsOf(reused.handle, 1) == prismsOf(fresh(warm, types)->handle, 1));
  REQUIRE(monolayerOf(reused.handle, warm, 1) == monolayerOf(fresh(warm, types)->handle, warm, 1));
  reused.set(small, &smallTypes);
  REQUIRE(chillOf(reused.handle, small.nop, Chill::plusOvito, 1) ==
          chillOf(fresh(small, smallTypes)->handle, small.nop, Chill::plusOvito, 1));
  REQUIRE(cagesOf(reused.handle, small.nop, 1) ==
          cagesOf(fresh(small, smallTypes)->handle, small.nop, 1));
}
