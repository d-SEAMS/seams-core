/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Carried primitive-ring update against a full ringNetwork on a
** production-stride dump. The graphs are the cutoff list and the seeded
** pair used by `seams cages` (mutual and union k-nearest, k = 4, candidate
** cutoff 5 A, rings through 6). A transition stays inside the ball when
** the update enumerates fewer sources than the molecule count. It falls
** back when that count is the whole vertex set.
**
**   bench_hop_stride time DUMP [--graph cutoff|seeded|both]
**                    [--reps 8] [--warmup 1] [--stride N] [--last N]
**   bench_hop_stride crossover DUMP [--frame-a 1] [--frame-b 2]
**                    [--graph cutoff|seeded|both] [--reps 8] [--warmup 1]
*/

#include <edge_set.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Clock = std::chrono::steady_clock;
using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
using Rows = std::vector<std::vector<int>>;

constexpr double kCutoff = 3.5;
constexpr double kCandidate = 5.0;
constexpr int kNeighbours = 4;
constexpr int kDepth = 6;

double ms(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

double quantile7(std::vector<double> values, double p) {
  if (values.empty()) {
    return 0.0;
  }
  std::sort(values.begin(), values.end());
  if (values.size() == 1) {
    return values[0];
  }
  const double h = (static_cast<double>(values.size()) - 1.0) * p;
  const auto lo = static_cast<std::size_t>(std::floor(h));
  const auto hi = static_cast<std::size_t>(std::ceil(h));
  const double w = h - static_cast<double>(lo);
  return values[lo] * (1.0 - w) + values[hi] * w;
}

int cpuCount() {
  const long n = sysconf(_SC_NPROCESSORS_ONLN);
  return n > 0 ? static_cast<int>(n) : 1;
}

int threadCount() {
#ifdef SEAMS_HAS_OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

Rows byIndex(Cloud &cloud, const std::vector<std::vector<int>> &ids) {
  return nneigh::neighbourListByIndex(cloud, ids);
}

struct Graphs {
  Rows cutoff;
  Rows strict;
  Rows permissive;
  bool seededSame = false;
};

Graphs graphsOf(Cloud &cloud, int atomType) {
  Graphs g;
  g.cutoff = byIndex(cloud, nneigh::neighListO(kCutoff, cloud, atomType));
  const auto pair =
      nneigh::kNearestNeighbourPair(cloud, kNeighbours, kCandidate, atomType);
  g.strict = byIndex(cloud, pair.first);
  g.permissive = byIndex(cloud, pair.second);
  g.seededSame = seams::provenEqual(seams::edgeSetOf(g.strict),
                                    seams::edgeSetOf(g.permissive));
  return g;
}

struct Series {
  std::string name;
  std::vector<Rows> frames;
  std::vector<int> changed;
  std::vector<int> recomputed;
  std::vector<char> same;
  std::vector<std::vector<double>> carry;
  std::vector<std::vector<double>> rebuild;
};

void classify(Series &series) {
  const int nFrames = static_cast<int>(series.frames.size());
  series.changed.assign(static_cast<std::size_t>(nFrames), 0);
  series.recomputed.assign(static_cast<std::size_t>(nFrames), 0);
  series.same.assign(static_cast<std::size_t>(nFrames), 1);
  primitive::RingUpdater updater(kDepth);
  for (int frame = 0; frame < nFrames; ++frame) {
    const auto &rows = series.frames[static_cast<std::size_t>(frame)];
    const auto full = primitive::ringNetwork(rows, kDepth);
    const auto &carried = updater.update(rows);
    series.same[static_cast<std::size_t>(frame)] = carried == full ? 1 : 0;
    series.changed[static_cast<std::size_t>(frame)] = updater.lastChangedBonds();
    series.recomputed[static_cast<std::size_t>(frame)] =
        updater.lastRecomputedSources();
  }
}

void timeSeries(Series &series, int warmup, int reps) {
  const int nFrames = static_cast<int>(series.frames.size());
  series.carry.assign(static_cast<std::size_t>(nFrames), {});
  series.rebuild.assign(static_cast<std::size_t>(nFrames), {});
  for (int pass = 0; pass < warmup + reps; ++pass) {
    primitive::RingUpdater updater(kDepth);
    for (int frame = 0; frame < nFrames; ++frame) {
      const auto t0 = Clock::now();
      (void)updater.update(series.frames[static_cast<std::size_t>(frame)]);
      const double dt = ms(t0, Clock::now());
      if (pass >= warmup) {
        series.carry[static_cast<std::size_t>(frame)].push_back(dt);
      }
    }
  }
  for (int frame = 0; frame < nFrames; ++frame) {
    const auto &rows = series.frames[static_cast<std::size_t>(frame)];
    for (int pass = 0; pass < warmup; ++pass) {
      (void)primitive::ringNetwork(rows, kDepth);
    }
    for (int rep = 0; rep < reps; ++rep) {
      const auto t0 = Clock::now();
      (void)primitive::ringNetwork(rows, kDepth);
      series.rebuild[static_cast<std::size_t>(frame)].push_back(ms(t0, Clock::now()));
    }
  }
}

bool fallbackAt(const Series &series, int frame) {
  const int nop = static_cast<int>(series.frames[static_cast<std::size_t>(frame)].size());
  return series.recomputed[static_cast<std::size_t>(frame)] == nop && nop > 0;
}

void printSamples(const Series &series, int stride) {
  for (int frame = 0; frame < static_cast<int>(series.frames.size()); ++frame) {
    const auto &carry = series.carry[static_cast<std::size_t>(frame)];
    const auto &rebuild = series.rebuild[static_cast<std::size_t>(frame)];
    const int n = static_cast<int>(std::min(carry.size(), rebuild.size()));
    for (int rep = 0; rep < n; ++rep) {
      std::printf("sample stride %d graph %s frame %d rep %d carry_ms %.6f rebuild_ms %.6f\n",
                  stride, series.name.c_str(), frame + 1, rep + 1, carry[static_cast<std::size_t>(rep)],
                  rebuild[static_cast<std::size_t>(rep)]);
    }
  }
}

void printSummary(const Series &series, int stride, bool transitionsOnly) {
  std::vector<double> carryMed;
  std::vector<double> rebuildMed;
  std::vector<double> bonds;
  int inside = 0;
  int fell = 0;
  int identical = 0;
  int transitions = 0;
  double carrySum = 0.0;
  double rebuildSum = 0.0;
  const int nFrames = static_cast<int>(series.frames.size());
  const int first = transitionsOnly ? 1 : 0;
  for (int frame = first; frame < nFrames; ++frame) {
    const auto &c = series.carry[static_cast<std::size_t>(frame)];
    const auto &r = series.rebuild[static_cast<std::size_t>(frame)];
    if (c.empty() || r.empty()) {
      continue;
    }
    const double cm = quantile7(c, 0.5);
    const double rm = quantile7(r, 0.5);
    carryMed.push_back(cm);
    rebuildMed.push_back(rm);
    carrySum += cm;
    rebuildSum += rm;
    bonds.push_back(static_cast<double>(series.changed[static_cast<std::size_t>(frame)]));
    const bool fellBack = fallbackAt(series, frame);
    inside += fellBack ? 0 : 1;
    fell += fellBack ? 1 : 0;
    identical += series.same[static_cast<std::size_t>(frame)] ? 1 : 0;
    ++transitions;
    std::printf(
        "frame stride %d graph %s index %d changed_bonds %d recomputed %d nop %zu "
        "path %s carry_ms %.6f rebuild_ms %.6f identical %s\n",
        stride, series.name.c_str(), frame + 1,
        series.changed[static_cast<std::size_t>(frame)],
        series.recomputed[static_cast<std::size_t>(frame)],
        series.frames[static_cast<std::size_t>(frame)].size(),
        fellBack ? "fallback" : "ball", cm, rm,
        series.same[static_cast<std::size_t>(frame)] ? "yes" : "no");
  }
  const double speedup = carrySum > 0.0 ? rebuildSum / carrySum : 0.0;
  const bool allSame = transitions > 0 && identical == transitions;
  std::printf(
      "summary stride %d graph %s transitions %d inside_ball %d fell_back %d "
      "median_changed_bonds %.6f bonds_q1 %.6f bonds_q3 %.6f "
      "carry_ms %.6f carry_q1 %.6f carry_q3 %.6f "
      "rebuild_ms %.6f rebuild_q1 %.6f rebuild_q3 %.6f "
      "speedup %.6f identical %s\n",
      stride, series.name.c_str(), transitions, inside, fell, quantile7(bonds, 0.5),
      quantile7(bonds, 0.25), quantile7(bonds, 0.75), quantile7(carryMed, 0.5),
      quantile7(carryMed, 0.25), quantile7(carryMed, 0.75), quantile7(rebuildMed, 0.5),
      quantile7(rebuildMed, 0.25), quantile7(rebuildMed, 0.75), speedup,
      allSame ? "yes" : "no");
}

Series combineSeeded(const Series &strict, const Series &permissive, const std::vector<char> &same) {
  Series out;
  out.name = "seeded";
  const int nFrames = static_cast<int>(strict.frames.size());
  out.frames.resize(static_cast<std::size_t>(nFrames));
  out.changed = permissive.changed;
  out.recomputed.assign(static_cast<std::size_t>(nFrames), 0);
  out.same.assign(static_cast<std::size_t>(nFrames), 1);
  out.carry.assign(static_cast<std::size_t>(nFrames), {});
  out.rebuild.assign(static_cast<std::size_t>(nFrames), {});
  for (int frame = 0; frame < nFrames; ++frame) {
    const bool one = same[static_cast<std::size_t>(frame)] != 0;
    const int nop = static_cast<int>(strict.frames[static_cast<std::size_t>(frame)].size());
    // A shared graph is one ring search, matching seams cages.
    const int strictN = strict.recomputed[static_cast<std::size_t>(frame)];
    const int permN = one ? strictN : permissive.recomputed[static_cast<std::size_t>(frame)];
    const bool fell = strictN == nop || permN == nop;
    out.recomputed[static_cast<std::size_t>(frame)] = fell ? nop : std::max(strictN, permN);
    out.same[static_cast<std::size_t>(frame)] =
        strict.same[static_cast<std::size_t>(frame)] &&
                (one || permissive.same[static_cast<std::size_t>(frame)])
            ? 1
            : 0;
    out.changed[static_cast<std::size_t>(frame)] =
        one ? strict.changed[static_cast<std::size_t>(frame)]
            : std::max(strict.changed[static_cast<std::size_t>(frame)],
                       permissive.changed[static_cast<std::size_t>(frame)]);
    const auto &sc = strict.carry[static_cast<std::size_t>(frame)];
    const auto &pc = permissive.carry[static_cast<std::size_t>(frame)];
    const auto &sr = strict.rebuild[static_cast<std::size_t>(frame)];
    const auto &pr = permissive.rebuild[static_cast<std::size_t>(frame)];
    const int n = static_cast<int>(sc.size());
    for (int rep = 0; rep < n && rep < static_cast<int>(sr.size()); ++rep) {
      const double extraC = one || rep >= static_cast<int>(pc.size())
                                ? 0.0
                                : pc[static_cast<std::size_t>(rep)];
      const double extraR = one || rep >= static_cast<int>(pr.size())
                                ? 0.0
                                : pr[static_cast<std::size_t>(rep)];
      out.carry[static_cast<std::size_t>(frame)].push_back(sc[static_cast<std::size_t>(rep)] + extraC);
      out.rebuild[static_cast<std::size_t>(frame)].push_back(
          sr[static_cast<std::size_t>(rep)] + extraR);
    }
    // The combined row has no neighbour rows of its own. nop is recovered
    // from recomputed on a fallback and from the strict frame otherwise.
    if (!fell) {
      out.frames[static_cast<std::size_t>(frame)].resize(static_cast<std::size_t>(nop));
    } else {
      out.frames[static_cast<std::size_t>(frame)].resize(static_cast<std::size_t>(nop));
    }
  }
  return out;
}

bool want(const std::string &graph, const char *name) {
  return graph == "both" || graph == name;
}

int runTime(const std::string &dump, const std::string &graph, int reps, int warmup,
            int stride, int last, int atomType) {
  const int nFile = sinp::nLammpsFrames(dump);
  const int nFrames = last > 0 ? std::min(last, nFile) : nFile;
  if (nFrames < 2) {
    std::fprintf(stderr, "need at least two frames in %s\n", dump.c_str());
    return 2;
  }
  Series cutoff;
  Series strict;
  Series permissive;
  cutoff.name = "cutoff";
  strict.name = "seeded-strict";
  permissive.name = "seeded-union";
  std::vector<char> seededSame(static_cast<std::size_t>(nFrames), 0);
  const bool doCutoff = want(graph, "cutoff");
  const bool doSeeded = graph == "both" || graph == "seeded";
  for (int frame = 1; frame <= nFrames; ++frame) {
    Cloud cloud;
    cloud = sinp::readLammpsTrjO(dump, frame, cloud, atomType);
    if (cloud.nop == 0) {
      std::fprintf(stderr, "empty frame %d\n", frame);
      return 2;
    }
    Graphs g = graphsOf(cloud, atomType);
    seededSame[static_cast<std::size_t>(frame - 1)] = g.seededSame ? 1 : 0;
    if (doCutoff) {
      cutoff.frames.push_back(std::move(g.cutoff));
    }
    if (doSeeded) {
      strict.frames.push_back(std::move(g.strict));
      if (!g.seededSame) {
        permissive.frames.push_back(std::move(g.permissive));
      } else {
        permissive.frames.push_back(strict.frames.back());
      }
    }
  }
  std::printf("# cores %d\n# omp_threads %d\n# reps %d\n# warmup %d\n# quantile type-7\n",
              cpuCount(), threadCount(), reps, warmup);
  std::printf("# depth %d\n# cutoff %.1f\n# k %d\n# candidate %.1f\n", kDepth, kCutoff,
              kNeighbours, kCandidate);
  bool identical = true;
  auto finish = [&](Series &series) {
    classify(series);
    timeSeries(series, warmup, reps);
    for (const char ok : series.same) {
      identical = identical && ok;
    }
    printSamples(series, stride);
    printSummary(series, stride, true);
  };
  if (doCutoff) {
    finish(cutoff);
  }
  if (doSeeded) {
    finish(strict);
    finish(permissive);
    Series seeded = combineSeeded(strict, permissive, seededSame);
    // classify() was already run on the parts. The combined row must not
    // recompute rings; its flags were copied. printSummary reads them.
    printSamples(seeded, stride);
    printSummary(seeded, stride, true);
    for (const char ok : seeded.same) {
      identical = identical && ok;
    }
  }
  if (!identical) {
    std::fprintf(stderr, "carried rings differ from a full search\n");
    return 1;
  }
  return 0;
}

struct SweepPoint {
  double radius = 0;
  int bonds = 0;
  int recomputed = 0;
  int nop = 0;
  double carry = 0;
  double rebuild = 0;
  bool identical = false;
};

SweepPoint timeSplice(const Rows &base, const Rows &moved, int warmup, int reps) {
  SweepPoint point;
  point.nop = static_cast<int>(moved.size());
  const auto full = primitive::ringNetwork(moved, kDepth);
  primitive::RingUpdater updater(kDepth);
  (void)updater.update(base);
  const auto &carried = updater.update(moved);
  point.identical = carried == full;
  point.bonds = updater.lastChangedBonds();
  point.recomputed = updater.lastRecomputedSources();
  for (int pass = 0; pass < warmup; ++pass) {
    primitive::RingUpdater warm(kDepth);
    (void)warm.update(base);
    (void)warm.update(moved);
    (void)primitive::ringNetwork(moved, kDepth);
  }
  std::vector<double> carry;
  std::vector<double> rebuild;
  for (int rep = 0; rep < reps; ++rep) {
    primitive::RingUpdater timed(kDepth);
    (void)timed.update(base);
    const auto t0 = Clock::now();
    (void)timed.update(moved);
    carry.push_back(ms(t0, Clock::now()));
    const auto t1 = Clock::now();
    (void)primitive::ringNetwork(moved, kDepth);
    rebuild.push_back(ms(t1, Clock::now()));
  }
  point.carry = quantile7(carry, 0.5);
  point.rebuild = quantile7(rebuild, 0.5);
  return point;
}

Cloud splice(const Cloud &a, const Cloud &b, double radius) {
  Cloud c = a;
  std::vector<int> index(static_cast<std::size_t>(b.nop + 1), -1);
  for (int i = 0; i < b.nop; ++i) {
    const int id = b.pts[static_cast<std::size_t>(i)].atomID;
    if (id >= 0 && id < static_cast<int>(index.size())) {
      index[static_cast<std::size_t>(id)] = i;
    }
  }
  for (int i = 0; i < c.nop; ++i) {
    double d2 = 0.0;
    const double p[3] = {c.pts[static_cast<std::size_t>(i)].x,
                         c.pts[static_cast<std::size_t>(i)].y,
                         c.pts[static_cast<std::size_t>(i)].z};
    for (int axis = 0; axis < 3; ++axis) {
      double d = p[axis] - (c.boxLow[static_cast<std::size_t>(axis)] +
                            0.5 * c.box[static_cast<std::size_t>(axis)]);
      const double length = c.box[static_cast<std::size_t>(axis)];
      if (length > 0.0) {
        d -= length * std::round(d / length);
      }
      d2 += d * d;
    }
    const int id = c.pts[static_cast<std::size_t>(i)].atomID;
    if (d2 >= radius * radius || id < 0 || id >= static_cast<int>(index.size())) {
      continue;
    }
    const int j = index[static_cast<std::size_t>(id)];
    if (j < 0) {
      continue;
    }
    c.pts[static_cast<std::size_t>(i)].x = b.pts[static_cast<std::size_t>(j)].x;
    c.pts[static_cast<std::size_t>(i)].y = b.pts[static_cast<std::size_t>(j)].y;
    c.pts[static_cast<std::size_t>(i)].z = b.pts[static_cast<std::size_t>(j)].z;
  }
  return c;
}

void reportCrossover(const std::string &name, std::vector<SweepPoint> points) {
  std::sort(points.begin(), points.end(),
            [](const SweepPoint &a, const SweepPoint &b) { return a.bonds < b.bonds; });
  int still = -1;
  int loses = -1;
  for (const auto &point : points) {
    std::printf(
        "sweep graph %s radius %.3f changed_bonds %d recomputed %d nop %d "
        "carry_ms %.6f rebuild_ms %.6f speedup %.6f path %s identical %s\n",
        name.c_str(), point.radius, point.bonds, point.recomputed, point.nop, point.carry,
        point.rebuild, point.carry > 0.0 ? point.rebuild / point.carry : 0.0,
        point.recomputed == point.nop ? "fallback" : "ball",
        point.identical ? "yes" : "no");
  }
  for (std::size_t i = 0; i < points.size(); ++i) {
    bool restLoses = true;
    for (std::size_t j = i; j < points.size(); ++j) {
      if (points[j].carry < points[j].rebuild) {
        restLoses = false;
        break;
      }
    }
    if (restLoses) {
      loses = points[i].bonds;
      break;
    }
    still = points[i].bonds;
  }
  std::printf("crossover graph %s still_wins %d loses %d\n", name.c_str(), still, loses);
}

int runCrossover(const std::string &dump, const std::string &graph, int frameA, int frameB,
                 int reps, int warmup, int atomType) {
  Cloud a;
  Cloud b;
  a = sinp::readLammpsTrjO(dump, frameA, a, atomType);
  b = sinp::readLammpsTrjO(dump, frameB, b, atomType);
  if (a.nop == 0 || b.nop == 0) {
    std::fprintf(stderr, "empty crossover frame\n");
    return 2;
  }
  const Graphs base = graphsOf(a, atomType);
  const double radii[] = {0, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20,
                          24, 28, 34, 42, 1e9};
  std::vector<SweepPoint> cutoff;
  std::vector<SweepPoint> strict;
  std::vector<SweepPoint> permissive;
  bool identical = true;
  for (const double radius : radii) {
    Cloud moved = splice(a, b, radius);
    Graphs g = graphsOf(moved, atomType);
    if (want(graph, "cutoff") || graph == "both") {
      SweepPoint point = timeSplice(base.cutoff, g.cutoff, warmup, reps);
      point.radius = radius >= 1e8 ? -1.0 : radius;
      identical = identical && point.identical;
      cutoff.push_back(point);
    }
    if (graph == "both" || graph == "seeded") {
      SweepPoint s = timeSplice(base.strict, g.strict, warmup, reps);
      s.radius = radius >= 1e8 ? -1.0 : radius;
      identical = identical && s.identical;
      strict.push_back(s);
      SweepPoint u = timeSplice(base.permissive, g.permissive, warmup, reps);
      u.radius = radius >= 1e8 ? -1.0 : radius;
      identical = identical && u.identical;
      permissive.push_back(u);
    }
  }
  std::printf("# cores %d\n# omp_threads %d\n# reps %d\n# warmup %d\n# quantile type-7\n",
              cpuCount(), threadCount(), reps, warmup);
  std::printf("# crossover frames %d %d\n", frameA, frameB);
  if (!cutoff.empty()) {
    reportCrossover("cutoff", cutoff);
  }
  if (!strict.empty()) {
    reportCrossover("seeded-strict", strict);
    reportCrossover("seeded-union", permissive);
  }
  if (!identical) {
    std::fprintf(stderr, "crossover rings differ from a full search\n");
    return 1;
  }
  return 0;
}

const char *flag(int argc, char **argv, const char *name, const char *fallback) {
  for (int i = 0; i < argc - 1; ++i) {
    if (std::strcmp(argv[i], name) == 0) {
      return argv[i + 1];
    }
  }
  return fallback;
}

int numberFlag(int argc, char **argv, const char *name, int fallback) {
  const char *value = flag(argc, argv, name, nullptr);
  return value == nullptr ? fallback : std::atoi(value);
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::fprintf(stderr,
                 "usage: bench_hop_stride time|crossover DUMP [--graph cutoff|seeded|both] "
                 "[--reps 8] [--warmup 1] [--stride N] [--last N] "
                 "[--frame-a 1] [--frame-b 2] [--type 1]\n");
    return 2;
  }
  const std::string mode = argv[1];
  const std::string dump = argv[2];
  const std::string graph = flag(argc, argv, "--graph", "both");
  const int reps = numberFlag(argc, argv, "--reps", 8);
  const int warmup = numberFlag(argc, argv, "--warmup", 1);
  const int atomType = numberFlag(argc, argv, "--type", 1);
  if (reps < 1 || warmup < 0) {
    std::fprintf(stderr, "reps and warmup must be non-negative, reps at least 1\n");
    return 2;
  }
  if (graph != "cutoff" && graph != "seeded" && graph != "both") {
    std::fprintf(stderr, "graph must be cutoff, seeded or both\n");
    return 2;
  }
  if (mode == "time") {
    return runTime(dump, graph, reps, warmup, numberFlag(argc, argv, "--stride", 0),
                   numberFlag(argc, argv, "--last", 0), atomType);
  }
  if (mode == "crossover") {
    return runCrossover(dump, graph, numberFlag(argc, argv, "--frame-a", 1),
                        numberFlag(argc, argv, "--frame-b", 2), reps, warmup, atomType);
  }
  std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
  return 2;
}
