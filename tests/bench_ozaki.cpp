// Timings and label checks for the Ozaki / error-free product paths.
// One process, one OpenMP team (the size is OMP_NUM_THREADS).
// Writes one JSON object to the path in argv[1].

#include <bop.hpp>
#include <neighbours.hpp>
#include <ozaki.hpp>
#include <rdf.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Clock = std::chrono::steady_clock;
using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

struct Stats {
  double median = 0;
  double lo = 0;
  double hi = 0;
};

Stats summarise(std::vector<double> v) {
  Stats s;
  if (v.empty()) {
    return s;
  }
  std::sort(v.begin(), v.end());
  s.lo = v.front();
  s.hi = v.back();
  s.median = v[v.size() / 2];
  return s;
}

template <typename Fn> Stats time11(Fn &&fn) {
  fn();
  std::vector<double> samples;
  samples.reserve(11);
  for (int i = 0; i < 11; i++) {
    const auto t0 = Clock::now();
    fn();
    const auto t1 = Clock::now();
    samples.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
  }
  return summarise(samples);
}

Cloud loadLikeCli(const std::string &path, int frame) {
  Cloud cloud;
  cloud = sinp::readLammpsTrjO(path, frame, cloud, 2);
  if (cloud.nop > 0) {
    return cloud;
  }
  Cloud again;
  return sinp::readLammpsTrjO(path, frame, again, 1);
}

std::uint64_t hashLabels(const Cloud &cloud) {
  std::uint64_t h = 14695981039346656037ull;
  for (const auto &pt : cloud.pts) {
    const auto v = static_cast<unsigned>(pt.iceType);
    h ^= v;
    h *= 1099511628211ull;
  }
  return h;
}

std::string histOf(const Cloud &cloud) {
  std::map<int, int> h;
  for (const auto &pt : cloud.pts) {
    h[static_cast<int>(pt.iceType)]++;
  }
  std::ostringstream os;
  bool first = true;
  for (const auto &kv : h) {
    if (!first) {
      os << ",";
    }
    first = false;
    os << kv.first << ":" << kv.second;
  }
  return os.str();
}

const char *jsonNum(double x, char *buf, std::size_t n) {
  if (!std::isfinite(x)) {
    std::snprintf(buf, n, "null");
  } else {
    std::snprintf(buf, n, "%.9g", x);
  }
  return buf;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: bench_ozaki out.json\n");
    return 2;
  }
  const std::string outPath = argv[1];
  int threads = 1;
#ifdef SEAMS_HAS_OPENMP
  threads = omp_get_max_threads();
#endif
  const std::filesystem::path trajDir = std::filesystem::path(argv[0]).parent_path();
  // The binary lives in the build tree. Trajectories are passed by the
  // working directory: the harness runs this from the repository root.
  const std::filesystem::path root = "input/traj";
  std::vector<std::string> files;
  for (const auto &ent : std::filesystem::directory_iterator(root)) {
    if (ent.path().extension() == ".lammpstrj") {
      files.push_back(ent.path().string());
    }
  }
  std::sort(files.begin(), files.end());

  const seams::product::Mode modes[] = {
      seams::product::Mode::Fp64, seams::product::Mode::Dot2,
      seams::product::Mode::OzakiI, seams::product::Mode::OzakiII};

  std::FILE *out = std::fopen(outPath.c_str(), "w");
  if (out == nullptr) {
    std::perror(outPath.c_str());
    return 1;
  }
  std::fprintf(out, "{\n");
  std::fprintf(out, "  \"threads\": %d,\n", threads);
  std::fprintf(out, "  \"reps\": 11,\n");
  std::fprintf(out, "  \"warmup_discarded\": 1,\n");
  std::fprintf(out, "  \"amx\": %s,\n", seams::product::amxAvailable() ? "true" : "false");
  std::fprintf(out, "  \"avx512_vnni\": %s,\n",
               seams::product::vnniAvailable() ? "true" : "false");
  std::fprintf(out, "  \"timing\": \"steady_clock milliseconds, one untimed call then 11\",\n");
  std::fprintf(out, "  \"labels\": [\n");

  bool firstFile = true;
  for (const auto &file : files) {
    const int nframes = sinp::nLammpsFrames(file);
    for (int frame = 1; frame <= nframes; frame++) {
      Cloud base = loadLikeCli(file, frame);
      if (base.nop == 0) {
        continue;
      }
      const double cand = 3.5 + 1.5;
      const int typ = base.pts[0].type;
      auto nList = nneigh::kNearestNeighbourList(base, 4, cand, typ, true);
      if (!firstFile) {
        std::fprintf(out, ",\n");
      }
      firstFile = false;
      std::fprintf(out, "    {\"file\": \"%s\", \"frame\": %d, \"nop\": %d",
                   file.c_str(), frame, base.nop);
      std::uint64_t fp64Hash = 0;
      std::string fp64Hist;
      std::vector<double> fp64C;
      double minMargin = 1e300;
      for (auto mode : modes) {
        seams::product::setMode(mode);
        seams::product::setReduction(seams::product::Reduction::Left);
        Cloud cloud = base;
        chill::getCorrelPlus(cloud, nList, false);
        chill::getIceTypePlusNoPrint(cloud, nList, false);
        const auto h = hashLabels(cloud);
        const auto hist = histOf(cloud);
        if (mode == seams::product::Mode::Fp64) {
          fp64Hash = h;
          fp64Hist = hist;
          for (const auto &pt : cloud.pts) {
            for (const auto &c : pt.c_ij) {
              fp64C.push_back(c.c_value);
              for (double thr : {-0.8, -0.35, 0.25}) {
                minMargin = std::min(minMargin, std::fabs(c.c_value - thr));
              }
            }
          }
        }
        int mismatches = 0;
        double maxAbs = 0.0;
        std::size_t k = 0;
        for (const auto &pt : cloud.pts) {
          for (const auto &c : pt.c_ij) {
            if (k < fp64C.size()) {
              maxAbs = std::max(maxAbs, std::fabs(c.c_value - fp64C[k]));
            }
            ++k;
          }
        }
        if (h != fp64Hash) {
          mismatches = 1;
        }
        char mb[64];
        std::fprintf(out,
                     ", \"%s\": {\"hash\": \"%llu\", \"hist\": \"%s\", \"same_labels\": %s, "
                     "\"max_abs_dc\": %s}",
                     seams::product::modeName(mode),
                     static_cast<unsigned long long>(h), hist.c_str(),
                     mismatches ? "false" : "true", jsonNum(maxAbs, mb, 64));
      }
      char margin[64];
      std::fprintf(out, ", \"min_threshold_margin\": %s}", jsonNum(minMargin, margin, 64));
    }
  }
  std::fprintf(out, "\n  ],\n");

  std::fprintf(out, "  \"kernels\": [\n");
  bool firstK = true;
  auto emitKernel = [&](const std::string &file, int frame, const char *name, const Stats &s) {
    if (!firstK) {
      std::fprintf(out, ",\n");
    }
    firstK = false;
    char a[64], b[64], c[64];
    std::fprintf(out,
                 "    {\"file\": \"%s\", \"frame\": %d, \"kernel\": \"%s\", "
                 "\"median_ms\": %s, \"min_ms\": %s, \"max_ms\": %s}",
                 file.c_str(), frame, name, jsonNum(s.median, a, 64), jsonNum(s.lo, b, 64),
                 jsonNum(s.hi, c, 64));
  };

  for (const auto &file : files) {
    Cloud base = loadLikeCli(file, 1);
    if (base.nop == 0) {
      continue;
    }
    const double cand = 5.0;
    const int typ = base.pts[0].type;
    Cloud cloud = base;
    auto nList = nneigh::kNearestNeighbourList(cloud, 4, cand, typ, true);
    emitKernel(file, 1, "knn4", time11([&] { nList = nneigh::kNearestNeighbourList(base, 4, cand, typ, true); }));
    for (auto mode : modes) {
      std::string key = std::string("chill_plus_") + seams::product::modeName(mode);
      emitKernel(file, 1, key.c_str(), time11([&] {
                   seams::product::setMode(mode);
                   seams::product::setReduction(seams::product::Reduction::Left);
                   cloud = base;
                   chill::getCorrelPlus(cloud, nList, false);
                   chill::getIceTypePlusNoPrint(cloud, nList, false);
                 }));
    }
    seams::product::setMode(seams::product::Mode::Fp64);
    seams::product::setReduction(seams::product::Reduction::Left);
    emitKernel(file, 1, "steinhardt_q3_q12_fp64", time11([&] {
                 chill::steinhardtQl(base, nList, 3);
                 chill::steinhardtQl(base, nList, 12);
               }));
    seams::product::setReduction(seams::product::Reduction::Pairwise);
    emitKernel(file, 1, "steinhardt_q3_q12_pairwise", time11([&] {
                 chill::steinhardtQl(base, nList, 3);
                 chill::steinhardtQl(base, nList, 12);
               }));
    seams::product::setReduction(seams::product::Reduction::Left);
    seams::product::setMode(seams::product::Mode::OzakiI);
    emitKernel(file, 1, "steinhardt_q3_q12_ozaki1", time11([&] {
                 chill::steinhardtQl(base, nList, 3);
                 chill::steinhardtQl(base, nList, 12);
               }));
    seams::product::setMode(seams::product::Mode::Fp64);
    emitKernel(file, 1, "rdf", time11([&] { (void)rdf::partialRdf(base, typ, typ, 8.0, 80); }));

    // q_l bit identity across the reduction, one frame.
    seams::product::setMode(seams::product::Mode::Fp64);
    seams::product::setReduction(seams::product::Reduction::Left);
    const auto qLeft = chill::steinhardtQl(base, nList, 6);
    seams::product::setReduction(seams::product::Reduction::Pairwise);
    const auto qPair = chill::steinhardtQl(base, nList, 6);
    double maxDq = 0.0;
    int nDiff = 0;
    for (std::size_t i = 0; i < qLeft.ql.size(); i++) {
      const double d = std::fabs(qLeft.ql[i] - qPair.ql[i]);
      maxDq = std::max(maxDq, d);
      if (qLeft.ql[i] != qPair.ql[i]) {
        ++nDiff;
      }
    }
    if (!firstK) {
      std::fprintf(out, ",\n");
    }
    firstK = false;
    char md[64];
    std::fprintf(out,
                 "    {\"file\": \"%s\", \"frame\": 1, \"kernel\": \"q6_left_vs_pairwise\", "
                 "\"n_bits_differ\": %d, \"max_abs\": %s}",
                 file.c_str(), nDiff, jsonNum(maxDq, md, 64));
    seams::product::setReduction(seams::product::Reduction::Left);
  }

  // Synthetic GEMMs. These shapes are not d-SEAMS kernels. They locate the
  // size at which Ozaki-II INT8 can beat an AVX-512 FP64 GEMM on this CPU.
  std::fprintf(out, "\n  ],\n  \"synthetic_gemm\": [\n");
  bool firstG = true;
  for (int n : {32, 64, 128, 256, 512}) {
    const auto nn = static_cast<std::size_t>(n);
    std::vector<double> A(nn * nn);
    std::vector<double> B(nn * nn);
    std::vector<double> C(nn * nn);
    std::srand(1);
    for (std::size_t i = 0; i < A.size(); i++) {
      A[i] = (std::rand() % 2000) / 1000.0 - 1.0;
      B[i] = (std::rand() % 2000) / 1000.0 - 1.0;
    }
    const auto fp = time11([&] {
      seams::product::gemm(n, n, n, A.data(), n, B.data(), n, C.data(), n,
                           seams::product::Mode::Fp64, nullptr);
    });
    const auto oz = time11([&] {
      seams::product::gemm(n, n, n, A.data(), n, B.data(), n, C.data(), n,
                           seams::product::Mode::OzakiII, nullptr,
                           seams::product::Int8Backend::Auto);
    });
    // One accuracy sample against the FP64 result, not a certificate.
    std::vector<double> Cfp(C.size());
    std::vector<double> Coz(C.size());
    seams::product::gemm(n, n, n, A.data(), n, B.data(), n, Cfp.data(), n,
                         seams::product::Mode::Fp64, nullptr);
    const auto info = seams::product::gemm(n, n, n, A.data(), n, B.data(), n, Coz.data(), n,
                                           seams::product::Mode::OzakiII, nullptr,
                                           seams::product::Int8Backend::Auto);
    double maxAbs = 0.0;
    for (std::size_t i = 0; i < Cfp.size(); i++) {
      maxAbs = std::max(maxAbs, std::fabs(Cfp[i] - Coz[i]));
    }
    if (!firstG) {
      std::fprintf(out, ",\n");
    }
    firstG = false;
    char a[64], b[64], c[64], d[64], e[64], f[64], g[64], h[64], i[64];
    const double flops = 2.0 * static_cast<double>(n) * n * n;
    std::fprintf(out,
                 "    {\"n\": %d, \"backend\": \"%s\", \"moduli\": %d, "
                 "\"fp64_median_ms\": %s, \"fp64_min_ms\": %s, \"fp64_max_ms\": %s, "
                 "\"ozaki2_median_ms\": %s, \"ozaki2_min_ms\": %s, \"ozaki2_max_ms\": %s, "
                 "\"max_abs_vs_fp64\": %s, \"fp64_gflops\": %s, \"ozaki2_gflops\": %s}",
                 n, info.backend, info.moduli, jsonNum(fp.median, a, 64), jsonNum(fp.lo, b, 64),
                 jsonNum(fp.hi, c, 64), jsonNum(oz.median, d, 64), jsonNum(oz.lo, e, 64),
                 jsonNum(oz.hi, f, 64), jsonNum(maxAbs, g, 64),
                 jsonNum(fp.median > 0 ? flops / (fp.median * 1e6) : 0, h, 64),
                 jsonNum(oz.median > 0 ? flops / (oz.median * 1e6) : 0, i, 64));
  }
  // Short dot, the CHILL shape: length 7 complex is 14 real multiplies.
  {
    std::vector<double> x(14), y(14);
    for (int i = 0; i < 14; i++) {
      x[static_cast<std::size_t>(i)] = 0.1 * (i - 7);
      y[static_cast<std::size_t>(i)] = 0.05 * (i + 1);
    }
    auto batch = [&](seams::product::Mode m) {
      double acc = 0;
      for (int i = 0; i < 20000; i++) {
        acc += seams::product::realDot(x.data(), y.data(), 14, m).value;
      }
      if (acc == 1e300) {
        std::fprintf(stderr, "unreachable\n");
      }
    };
    const auto fp = time11([&] { batch(seams::product::Mode::Fp64); });
    const auto d2 = time11([&] { batch(seams::product::Mode::Dot2); });
    const auto o1 = time11([&] { batch(seams::product::Mode::OzakiI); });
    const auto o2 = time11([&] { batch(seams::product::Mode::OzakiII); });
    std::fprintf(out, ",\n");
    char a[64], b[64], c[64], d[64];
    std::fprintf(out,
                 "    {\"n\": 14, \"note\": \"20000 real dots of length 14 per sample\", "
                 "\"fp64_median_ms\": %s, \"dot2_median_ms\": %s, \"ozaki1_median_ms\": %s, "
                 "\"ozaki2_median_ms\": %s}",
                 jsonNum(fp.median, a, 64), jsonNum(d2.median, b, 64), jsonNum(o1.median, c, 64),
                 jsonNum(o2.median, d, 64));
  }
  std::fprintf(out, "\n  ]\n}\n");
  std::fclose(out);
  (void)trajDir;
  seams::product::setMode(seams::product::Mode::Fp64);
  seams::product::setReduction(seams::product::Reduction::Left);
  return 0;
}
