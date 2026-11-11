#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <thread_team.hpp>

#include <cmath>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

molSys::PointCloud<molSys::Point<double>, double> jitteredIce(int cells) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double bond = 2.75;
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  const double side = cells * lattice;
  const double basis[8][3] = {{0, 0, 0},       {0.5, 0.5, 0},   {0.5, 0, 0.5},
                              {0, 0.5, 0.5},   {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  unsigned s = 1;
  int id = 1;
  for (int i = 0; i < cells; i++) {
    for (int j = 0; j < cells; j++) {
      for (int k = 0; k < cells; k++) {
        for (const auto &b : basis) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          const double xyz[3] = {i + b[0], j + b[1], k + b[2]};
          double *out[3] = {&p.x, &p.y, &p.z};
          for (int c = 0; c < 3; c++) {
            s = s * 1664525u + 1013904223u;
            const double jitter = (static_cast<int>(s % 1000) / 1000.0 - 0.5) * 0.3;
            *out[c] = std::fmod((xyz[c] * lattice + jitter), side);
            if (*out[c] < 0.0) {
              *out[c] += side;
            }
          }
          cloud.pts.push_back(p);
          cloud.idIndexMap[id] = id - 1;
          id++;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.currentFrame = 1;
  cloud.box = {side, side, side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

struct Snapshot {
  std::vector<int> published;
  std::vector<int> ovito;
  std::vector<double> q6;
  std::vector<double> q6bar;
  std::vector<int> rings;
  std::vector<std::vector<int>> ringList;
};

std::vector<std::vector<int>> rowsByIndex(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const std::vector<std::vector<int>> &byId) {
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; i++) {
    rows[static_cast<std::size_t>(i)].push_back(i);
    if (static_cast<int>(byId.size()) <= i) {
      continue;
    }
    for (std::size_t k = 1; k < byId[static_cast<std::size_t>(i)].size(); k++) {
      const auto it = cloud.idIndexMap.find(byId[static_cast<std::size_t>(i)][k]);
      if (it != cloud.idIndexMap.end()) {
        rows[static_cast<std::size_t>(i)].push_back(it->second);
      }
    }
  }
  return rows;
}

std::vector<int> census(const std::vector<std::vector<int>> &rings) {
  std::vector<int> counts(7, 0);
  for (const auto &ring : rings) {
    if (ring.size() < counts.size()) {
      counts[ring.size()]++;
    }
  }
  return counts;
}

Snapshot runOnce(const molSys::PointCloud<molSys::Point<double>, double> &base, int threads) {
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(threads);
#else
  (void)threads;
#endif
  auto cloud = base;
  const auto nList = nneigh::neighListO(3.5, cloud, 1);
  Snapshot snap;
  snap.published.resize(static_cast<std::size_t>(cloud.nop));
  snap.ovito.resize(static_cast<std::size_t>(cloud.nop));
  {
    auto labeled = cloud;
    chill::applyChillPlus(labeled, nList, chill::ChillPlusRule::Published);
    for (int i = 0; i < cloud.nop; i++) {
      snap.published[static_cast<std::size_t>(i)] =
          static_cast<int>(labeled.pts[static_cast<std::size_t>(i)].iceType);
    }
  }
  {
    auto labeled = cloud;
    chill::applyChillPlus(labeled, nList, chill::ChillPlusRule::Ovito);
    for (int i = 0; i < cloud.nop; i++) {
      snap.ovito[static_cast<std::size_t>(i)] =
          static_cast<int>(labeled.pts[static_cast<std::size_t>(i)].iceType);
    }
  }
  const auto q = chill::steinhardtQl(cloud, nList, 6);
  snap.q6 = q.ql;
  snap.q6bar = q.qlBar;
  snap.ringList = primitive::ringNetwork(rowsByIndex(cloud, nList), 6);
  snap.rings = census(snap.ringList);
  return snap;
}

} // namespace

TEST_CASE("one, two and four threads match the serial CHILL+, q6 and rings",
          "[thread]") {
  const auto base = jitteredIce(4);
  const Snapshot serial = runOnce(base, 1);
  REQUIRE(lc_num_threads() == 1);
  for (const int threads : {2, 4}) {
    const Snapshot got = runOnce(base, threads);
    REQUIRE(lc_num_threads() == threads);
    REQUIRE(got.published == serial.published);
    REQUIRE(got.ovito == serial.ovito);
    REQUIRE(got.rings == serial.rings);
    REQUIRE(got.ringList == serial.ringList);
    REQUIRE(got.q6.size() == serial.q6.size());
    for (std::size_t i = 0; i < serial.q6.size(); i++) {
      REQUIRE_THAT(got.q6[i], Catch::Matchers::WithinAbs(serial.q6[i], 1e-15));
      REQUIRE_THAT(got.q6bar[i], Catch::Matchers::WithinAbs(serial.q6bar[i], 1e-15));
    }
  }
}
