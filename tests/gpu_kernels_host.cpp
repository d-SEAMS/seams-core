// The CUDA kernels of src/gpu_batch.cpp compiled for the host: each
// kernel runs its grid one thread at a time, under the per-frame decisions
// of gpu::analyzeResident (grow an overflowed buffer and rerun, or label a
// frame with a row past kRowCap on the host). Every frame must match
// tum::hostCageCounts per atom, and ddc_affil must pass the repeated-witness
// fixtures of test_cage_affiliation.cpp. Exits 1 on any mismatch.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <random>
#include <vector>

#define __global__
#define __device__
struct Dim3 {
  int x = 0;
};
static Dim3 blockIdx, blockDim{1}, threadIdx;
inline int atomicAdd(int *p, int v) {
  const int o = *p;
  *p += v;
  return o;
}
inline int atomicMax(int *p, int v) {
  const int o = *p;
  if (v > o)
    *p = v;
  return o;
}
namespace k {
#include "gpu_batch_kernels.inc"
}
#undef __global__
#undef __device__

#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <tum_offload.hpp>

using Graph = std::vector<std::vector<int>>;
namespace {

// Ideal hexagonal ice: the lonsdaleite O network at 2.76 A, nx x ny x nz
// orthogonal cells of eight molecules in a periodic box
molSys::PointCloud<molSys::Point<double>, double> iceIh(int nx, int ny,
                                                        int nz) {
  const double a = 2.76 * std::sqrt(8.0 / 3.0);
  const double c = a * std::sqrt(8.0 / 3.0);
  const double b = a * std::sqrt(3.0);
  const double basis[4][3] = {{1.0 / 3, 2.0 / 3, 0.0},
                              {2.0 / 3, 1.0 / 3, 0.5},
                              {1.0 / 3, 2.0 / 3, 3.0 / 8},
                              {2.0 / 3, 1.0 / 3, 7.0 / 8}};
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud.box = {nx * a, ny * b, nz * c};
  yCloud.boxLow = {0.0, 0.0, 0.0};
  yCloud.currentFrame = 1;
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        for (const double shift : {0.0, 0.5}) {
          for (const auto &f : basis) {
            molSys::Point<double> p;
            p.type = 1;
            p.atomID = yCloud.nop;
            p.x =
                std::fmod(shift * a + f[0] * a - f[1] * a / 2 + ix * a, nx * a);
            p.y = std::fmod(shift * b + f[1] * b / 2 + iy * b, ny * b);
            p.z = (f[2] + iz) * c;
            yCloud.idIndexMap[p.atomID] = yCloud.nop;
            yCloud.pts.push_back(p);
            yCloud.nop++;
          }
        }
      }
    }
  }
  return yCloud;
}

} // namespace

template <class F> void grid(int n, F f) {
  for (int t = 0; t < n; ++t) {
    blockIdx.x = t;
    blockDim.x = 1;
    threadIdx.x = 0;
    f();
  }
}

Graph levi(int q) {
  std::vector<std::array<int, 3>> v;
  for (int a = 0; a < q; ++a)
    for (int b = 0; b < q; ++b)
      for (int c = 0; c < q; ++c) {
        const int lead = a != 0 ? a : (b != 0 ? b : c);
        if (lead == 1)
          v.push_back({a, b, c});
      }
  const int m = static_cast<int>(v.size());
  Graph g(2 * m);
  for (int i = 0; i < 2 * m; ++i)
    g[i].push_back(i);
  for (int p = 0; p < m; ++p)
    for (int l = 0; l < m; ++l)
      if ((v[p][0] * v[l][0] + v[p][1] * v[l][1] + v[p][2] * v[l][2]) % q ==
          0) {
        g[p].push_back(m + l);
        g[m + l].push_back(p);
      }
  return g;
}

Graph pad(Graph g, int n) {
  for (int i = static_cast<int>(g.size()); i < n; ++i)
    g.push_back({i});
  return g;
}

Graph relabel(const Graph &g, unsigned seed) {
  std::vector<int> perm(g.size());
  std::iota(perm.begin(), perm.end(), 0);
  std::mt19937 rng(seed);
  std::shuffle(perm.begin(), perm.end(), rng);
  Graph out(g.size());
  for (std::size_t a = 0; a < g.size(); ++a) {
    out[perm[a]].push_back(perm[a]);
    for (std::size_t t = 1; t < g[a].size(); ++t)
      out[perm[a]].push_back(perm[g[a][t]]);
  }
  return out;
}

struct Out {
  std::vector<int> six, hc, ddc, nRings;
  int retries = 0;
  std::vector<int> hostFrames;
};

// One batch of frames of equal size; graphs have degree at most 4.
Out batch(const std::vector<Graph> &frames, int maxPer, int maxRings,
          int maxPairs) {
  const int kRowCap = 32;
  const int nF = static_cast<int>(frames.size());
  const int n = static_cast<int>(frames[0].size());
  std::vector<int> mdeg(nF * n), mcols(nF * n * 4, -1);
  for (int f = 0; f < nF; ++f)
    for (int i = 0; i < n; ++i) {
      const auto &row = frames[f][i];
      mdeg[f * n + i] = static_cast<int>(row.size()) - 1;
      for (std::size_t t = 1; t < row.size(); ++t)
        mcols[(f * n + i) * 4 + t - 1] = row[t];
    }
  Out out;
  out.six.assign(nF * n, 0);
  out.hc.assign(nF * n, 0);
  out.ddc.assign(nF * n, 0);
  out.nRings.assign(nF, 0);
  std::vector<char> done(nF, 0), onHost(nF, 0);
  for (;;) {
    std::vector<int> nRings(nF, 0), rowMax(nF, 0), nPairs(nF, 0);
    std::vector<int> ringAtoms(static_cast<std::size_t>(nF) * maxRings * 6, -1);
    std::vector<int> tc(nF * n, 0),
        through(static_cast<std::size_t>(nF) * n * maxPer, -1);
    std::vector<int> hc(static_cast<std::size_t>(nF) * maxRings, 0),
        ddc(hc.size(), 0);
    std::vector<int> pairs(static_cast<std::size_t>(nF) * maxPairs * 2, -1);
    std::vector<int> aHc(nF * n, 0), aDdc(nF * n, 0), six(nF * n, 0);
    grid(nF * n, [&] {
      k::enum_six(mdeg.data(), mcols.data(), n, nF, maxRings, nRings.data(),
                  ringAtoms.data());
    });
    grid(nF * maxRings, [&] {
      k::invert_rings(nRings.data(), ringAtoms.data(), n, nF, maxRings, maxPer,
                      tc.data(), through.data());
    });
    grid(nF * n, [&] {
      k::sort_through(tc.data(), through.data(), n, nF, maxPer, rowMax.data());
    });
    grid(nF * maxRings, [&] {
      k::emit_basal(nRings.data(), ringAtoms.data(), mdeg.data(), mcols.data(),
                    tc.data(), through.data(), n, nF, maxRings, maxPer,
                    maxPairs, nPairs.data(), pairs.data());
    });
    grid(nF * maxPairs, [&] {
      k::apply_hc(nPairs.data(), pairs.data(), ringAtoms.data(), tc.data(),
                  through.data(), n, nF, maxRings, maxPer, maxPairs, hc.data());
    });
    grid(nF * maxRings, [&] {
      k::ddc_affil(nRings.data(), ringAtoms.data(), tc.data(), through.data(),
                   hc.data(), n, nF, maxRings, maxPer, ddc.data());
    });
    grid(nF * maxRings, [&] {
      k::atom_ice(nRings.data(), ringAtoms.data(), hc.data(), ddc.data(), n, nF,
                  maxRings, aHc.data(), aDdc.data(), six.data());
    });
    int needPer = maxPer, needRings = maxRings, needPairs = maxPairs;
    bool pending = false;
    for (int f = 0; f < nF; ++f) {
      if (done[f] || onHost[f])
        continue;
      if (rowMax[f] > kRowCap) {
        onHost[f] = 1;
        continue;
      }
      if (nRings[f] <= maxRings && rowMax[f] <= maxPer &&
          nPairs[f] <= maxPairs) {
        std::copy(six.begin() + f * n, six.begin() + (f + 1) * n,
                  out.six.begin() + f * n);
        std::copy(aHc.begin() + f * n, aHc.begin() + (f + 1) * n,
                  out.hc.begin() + f * n);
        std::copy(aDdc.begin() + f * n, aDdc.begin() + (f + 1) * n,
                  out.ddc.begin() + f * n);
        out.nRings[f] = nRings[f];
        done[f] = 1;
        continue;
      }
      needPer = std::max(needPer, rowMax[f]);
      needRings = std::max(needRings, nRings[f]);
      needPairs = std::max(needPairs, nPairs[f]);
      pending = true;
    }
    if (!pending)
      break;
    maxPer = needPer;
    maxRings = needRings;
    maxPairs = needPairs;
    ++out.retries;
  }
  for (int f = 0; f < nF; ++f) {
    if (!onHost[f])
      continue;
    const auto h = tum::hostCageCounts(frames[f]);
    std::copy(h.atomSix.begin(), h.atomSix.end(), out.six.begin() + f * n);
    std::copy(h.atomHc.begin(), h.atomHc.end(), out.hc.begin() + f * n);
    std::copy(h.atomDdc.begin(), h.atomDdc.end(), out.ddc.begin() + f * n);
    out.nRings[f] = h.nSix;
    out.hostFrames.push_back(f);
  }
  return out;
}

bool run(const char *name, const std::vector<Graph> &frames, int per, int rings,
         int pairs) {
  const Out o = batch(frames, per, rings, pairs);
  const int n = static_cast<int>(frames[0].size());
  bool ok = true;
  std::printf("%s: frames %zu atoms %d start caps row %d rings %d pairs %d, "
              "retries %d, host frames",
              name, frames.size(), n, per, rings, pairs, o.retries);
  for (int f : o.hostFrames)
    std::printf(" %d", f);
  std::printf("\n");
  for (std::size_t f = 0; f < frames.size(); ++f) {
    const auto h = tum::hostCageCounts(frames[f]);
    const bool same =
        o.nRings[f] == h.nSix &&
        std::equal(h.atomSix.begin(), h.atomSix.end(), o.six.begin() + f * n) &&
        std::equal(h.atomHc.begin(), h.atomHc.end(), o.hc.begin() + f * n) &&
        std::equal(h.atomDdc.begin(), h.atomDdc.end(), o.ddc.begin() + f * n);
    std::printf("  frame %zu nSix %d HC rings %d DDC rings %d HC atoms %d DDC "
                "atoms %d: kernels==host %s\n",
                f, h.nSix, h.nHcRings, h.nDdcRings, h.nHcAtoms, h.nDdcAtoms,
                same ? "yes" : "NO");
    ok = ok && same;
  }
  return ok;
}

// The repeated-witness fixtures of test_cage_affiliation.cpp: ring 0 has
// T_k entries that pass every intersection only with a P_k repeated, so it
// and its peripherals are DDC only once the extra ring supplies a sixth
// distinct peripheral; ring 6 is a filler and never is
bool ddcFixtures() {
  struct Fixture {
    const char *name;
    Graph rings;
    std::vector<int> extra;
  };
  const std::vector<Fixture> fixtures = {
      {"P_1 = P_0",
       {{0, 1, 2, 3, 4, 5},
        {0, 1, 2, 3, 6, 7},
        {2, 3, 4, 6, 8, 9},
        {4, 5, 0, 6, 7, 10},
        {3, 4, 5, 7, 6, 11},
        {5, 0, 1, 7, 12, 13},
        {1, 2, 14, 15, 16, 17}},
       {1, 2, 3, 7, 6, 12}},
      {"P_2 = P_0",
       {{0, 1, 2, 3, 4, 5},
        {0, 1, 2, 3, 4, 6},
        {4, 5, 0, 6, 8, 9},
        {1, 2, 3, 7, 10, 11},
        {3, 4, 5, 7, 10, 12},
        {5, 0, 1, 7, 11, 13},
        {2, 14, 15, 16, 17, 18}},
       {2, 3, 4, 6, 19, 20}},
  };
  const int maxPer = 32;
  bool ok = true;
  for (const auto &fx : fixtures) {
    for (const bool withExtra : {false, true}) {
      Graph rings = fx.rings;
      if (withExtra)
        rings.push_back(fx.extra);
      const int nR = static_cast<int>(rings.size());
      int n = 0;
      std::vector<int> ringAtoms;
      for (const auto &r : rings) {
        ringAtoms.insert(ringAtoms.end(), r.begin(), r.end());
        n = std::max(n, *std::max_element(r.begin(), r.end()) + 1);
      }
      std::vector<int> tc(n, 0),
          through(static_cast<std::size_t>(n) * maxPer, -1);
      for (int r = 0; r < nR; ++r)
        for (const int a : rings[r])
          through[a * maxPer + tc[a]++] = r;
      std::vector<int> hc(nR, 0), ddc(nR, 0);
      grid(nR, [&] {
        k::ddc_affil(&nR, ringAtoms.data(), tc.data(), through.data(),
                     hc.data(), n, 1, nR, maxPer, ddc.data());
      });
      bool same = true;
      for (int r = 0; r < nR; ++r)
        same = same && (ddc[r] != 0) == (withExtra && r != 6);
      std::printf("DDC fixture %s, extra ring %s: DDC rings %d: %s\n", fx.name,
                  withExtra ? "yes" : "no",
                  static_cast<int>(std::count_if(ddc.begin(), ddc.end(),
                                                 [](int v) { return v != 0; })),
                  same ? "as expected" : "NO");
      ok = ok && same;
    }
  }
  return ok;
}

int main() {
  const Graph heawood = pad(levi(2), 26), pg3 = levi(3);
  auto y = iceIh(4, 3, 3);
  auto nl = nneigh::neighListO(3.5, y, 1);
  const Graph ih = nneigh::neighbourListByIndex(y, nl);
  bool ok = true;
  // Planned caps of analyzeResident: 16 per row, 16 rings per atom, 8 pairs per
  // ring
  ok &= run("Heawood + PG(2,3), planned caps", {heawood, pg3}, 16, 26 * 16,
            26 * 16 * 8);
  ok &= run("Heawood + PG(2,3), caps 1", {heawood, pg3, relabel(heawood, 7)}, 1,
            1, 1);
  ok &= run("ice Ih + relabelled, planned caps", {ih, relabel(ih, 3)}, 16,
            288 * 16, 288 * 16 * 8);
  ok &= run("ice Ih + relabelled, caps 1", {ih, relabel(ih, 3)}, 1, 1, 1);
  ok &= ddcFixtures();
  std::printf("%s\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
