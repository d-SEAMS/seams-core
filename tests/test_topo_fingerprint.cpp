#include <catch2/catch_test_macros.hpp>

#include <cage_canon.hpp>
#include <topo_fingerprint.hpp>

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

namespace {

// Neighbour rows of a periodic cubic-diamond lattice with `reps` cells per
// side (8 atoms per cell, four neighbours each), by index, rows leading
// with the atom itself.
topo::Rows diamondRows(int reps) {
  const int n = 8 * reps * reps * reps;
  // basis: fcc sites and their (1/4,1/4,1/4) partners, in quarter-cell units
  const int fcc[4][3] = {{0, 0, 0}, {2, 2, 0}, {2, 0, 2}, {0, 2, 2}};
  auto index = [&](int x, int y, int z) {
    const int L = 4 * reps;
    x = ((x % L) + L) % L;
    y = ((y % L) + L) % L;
    z = ((z % L) + L) % L;
    const int cx = x / 4, cy = y / 4, cz = z / 4;
    const int fx = x % 4, fy = y % 4, fz = z % 4;
    int b = -1;
    for (int k = 0; k < 4; k++) {
      if (fx == fcc[k][0] && fy == fcc[k][1] && fz == fcc[k][2]) {
        b = k;
      }
      if (fx == fcc[k][0] + 1 && fy == fcc[k][1] + 1 && fz == fcc[k][2] + 1) {
        b = 4 + k;
      }
    }
    REQUIRE(b >= 0);
    return ((cx * reps + cy) * reps + cz) * 8 + b;
  };
  topo::Rows rows(static_cast<std::size_t>(n));
  for (int cx = 0; cx < reps; cx++) {
    for (int cy = 0; cy < reps; cy++) {
      for (int cz = 0; cz < reps; cz++) {
        for (int k = 0; k < 4; k++) {
          const int x = 4 * cx + fcc[k][0], y = 4 * cy + fcc[k][1], z = 4 * cz + fcc[k][2];
          const int a = index(x, y, z);
          const int b = index(x + 1, y + 1, z + 1);
          // a bonds to b and to the three other partners at (+1,-1,-1) etc.
          const int partners[4][3] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
          for (const auto &p : partners) {
            const int q = index(x + p[0], y + p[1], z + p[2]);
            rows[static_cast<std::size_t>(a)].push_back(q);
            rows[static_cast<std::size_t>(q)].push_back(a);
          }
          (void)b;
        }
      }
    }
  }
  for (int i = 0; i < n; i++) {
    auto &r = rows[static_cast<std::size_t>(i)];
    std::sort(r.begin(), r.end());
    r.erase(std::unique(r.begin(), r.end()), r.end());
    r.insert(r.begin(), i);
    REQUIRE(r.size() == 5);
  }
  return rows;
}

topo::Rows permute(const topo::Rows &rows, const std::vector<int> &perm) {
  const int n = static_cast<int>(rows.size());
  topo::Rows out(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    std::vector<int> row;
    for (int v : rows[static_cast<std::size_t>(i)]) {
      row.push_back(perm[static_cast<std::size_t>(v)]);
    }
    std::sort(row.begin() + 1, row.end());
    out[static_cast<std::size_t>(perm[static_cast<std::size_t>(i)])] = row;
  }
  return out;
}

} // namespace

TEST_CASE("every atom of a diamond lattice carries the same local key", "[topo]") {
  const auto rows = diamondRows(3);
  const auto fp = topo::fingerprint(rows, 2, 7);
  REQUIRE(fp.classes.size() == 1);
  REQUIRE(fp.atomKeys.size() == rows.size());
  REQUIRE(fp.ringCensus[6] == 2 * static_cast<int>(rows.size()));
  REQUIRE(fp.ringCensus[3] == 0);
  REQUIRE(fp.ringCensus[4] == 0);
  REQUIRE(fp.ringCensus[5] == 0);
  const auto lk = topo::localKey(rows, 0, 2);
  REQUIRE(lk.vertices == 1 + 4 + 12);
  REQUIRE(lk.edges == 4 + 12);
  REQUIRE(lk.method == (cage::nautyAvailable() ? "nauty" : "wl"));
}

TEST_CASE("the frame key is independent of atom numbering", "[topo]") {
  const auto rows = diamondRows(3);
  std::vector<int> perm(rows.size());
  std::iota(perm.begin(), perm.end(), 0);
  std::mt19937 rng(7);
  std::shuffle(perm.begin(), perm.end(), rng);
  const auto a = topo::fingerprint(rows, 2, 7);
  const auto b = topo::fingerprint(permute(rows, perm), 2, 7);
  REQUIRE(a.key == b.key);
  REQUIRE(a.classes == b.classes);
  REQUIRE(a.ringCensus == b.ringCensus);
}

TEST_CASE("one broken bond splits the classes and changes the frame key", "[topo]") {
  auto rows = diamondRows(3);
  const auto whole = topo::fingerprint(rows, 2, 7);
  // remove the bond between atom 0 and its first neighbour
  const int u = 0;
  const int v = rows[0][1];
  rows[0].erase(rows[0].begin() + 1);
  auto &rv = rows[static_cast<std::size_t>(v)];
  rv.erase(std::find(rv.begin() + 1, rv.end(), u));
  const auto broken = topo::fingerprint(rows, 2, 7);
  REQUIRE(broken.key != whole.key);
  REQUIRE(broken.classes.size() > 1);
  // the two endpoints share an environment (three neighbours, a vacancy
  // opposite), and the perfect atoms keep theirs
  REQUIRE(broken.atomKeys[0] == broken.atomKeys[static_cast<std::size_t>(v)]);
  REQUIRE(broken.atomKeys[0] != whole.atomKeys[0]);
  REQUIRE(broken.ringCensus[6] < whole.ringCensus[6]);
}

TEST_CASE("the refinement hash separates a triangle from a path", "[topo]") {
  const std::vector<std::vector<int>> triangle = {{1, 2}, {0, 2}, {0, 1}};
  const std::vector<std::vector<int>> path = {{1}, {0, 2}, {1}};
  REQUIRE(topo::wlHash(triangle, 0, 2) != topo::wlHash(path, 0, 2));
  REQUIRE(topo::wlHash(path, 0, 2) != topo::wlHash(path, 1, 2));
  REQUIRE(topo::wlHash(triangle, 0, 2) == topo::wlHash(triangle, 1, 2));
}

TEST_CASE("hop neighbourhoods grow by one shell per hop", "[topo]") {
  const auto rows = diamondRows(3);
  REQUIRE(topo::hopNeighbourhood(rows, 5, 0).size() == 1);
  REQUIRE(topo::hopNeighbourhood(rows, 5, 1).size() == 5);
  REQUIRE(topo::hopNeighbourhood(rows, 5, 2).size() == 17);
  REQUIRE(topo::hopNeighbourhood(rows, 5, 2)[0] == 5);
}

TEST_CASE("vertex colours split the classes and travel with a permutation", "[topo]") {
  const auto rows = diamondRows(3);
  const int n = static_cast<int>(rows.size());
  // the two sublattices of the diamond cell are the two species of zincblende
  std::vector<int> colours(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    colours[static_cast<std::size_t>(i)] = (i % 8 < 4) ? 1 : 2;
  }
  const auto fp = topo::fingerprint(rows, 2, 7, colours);
  REQUIRE(fp.classes.size() == 2);
  for (const auto &kv : fp.classes) {
    REQUIRE(kv.second == n / 2);
  }
  REQUIRE(fp.key != topo::fingerprint(rows, 2, 7).key);
  std::vector<int> perm(rows.size());
  std::iota(perm.begin(), perm.end(), 0);
  std::mt19937 rng(11);
  std::shuffle(perm.begin(), perm.end(), rng);
  std::vector<int> permColours(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    permColours[static_cast<std::size_t>(perm[static_cast<std::size_t>(i)])] =
        colours[static_cast<std::size_t>(i)];
  }
  const auto fq = topo::fingerprint(permute(rows, perm), 2, 7, permColours);
  REQUIRE(fq.key == fp.key);
  REQUIRE(fq.classes == fp.classes);
  // the local key of an atom and of its image agree
  REQUIRE(fq.atomKeys[static_cast<std::size_t>(perm[3])] == fp.atomKeys[3]);
}
