#include <catch2/catch_test_macros.hpp>

#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <map>
#include <vector>

namespace {

std::vector<std::vector<int>>
sixMembered(const std::vector<std::vector<int>> &rings) {
  std::vector<std::vector<int>> six;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  return six;
}

std::vector<int> canonicalKeyOf(const std::vector<int> &ringAtoms) {
  const size_t k = ringAtoms.size();
  std::vector<int> best, candidate(k);
  for (int dir = 0; dir < 2; dir++) {
    for (size_t start = 0; start < k; start++) {
      for (size_t m = 0; m < k; m++) {
        const size_t idx = (dir == 0) ? (start + m) % k : (start + k - m) % k;
        candidate[m] = ringAtoms[idx];
      }
      if (best.empty() || candidate < best) {
        best = candidate;
      }
    }
  }
  return best;
}

// The twelve-atom hexagonal cage used across the topo_bulk tests
molSys::PointCloud<molSys::Point<double>, double> hcCloud() {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud.box = {53.9690018, 54.5289994, 51.257};
  yCloud.boxLow = {0.0, 0.0, 0.0};
  yCloud.currentFrame = 1;
  molSys::Point<double> iPoint;
  iPoint.type = 1;
  double coords[12][3] = {
      {8.995, 10.3859997, 15.0939999},     {6.7459998, 14.2810001, 15.0939999},
      {4.4970002, 10.3859997, 15.0939999}, {8.995, 12.9829998, 14.1949997},
      {6.7459998, 9.0880003, 14.1949997},  {4.4970002, 12.9829998, 14.1949997},
      {8.995, 12.9829998, 11.4329996},     {6.7459998, 9.0880003, 11.4329996},
      {4.4970002, 12.9829998, 11.4329996}, {8.995, 10.3859997, 10.5340004},
      {6.7459998, 14.2810001, 10.5340004}, {4.4970002, 10.3859997, 10.5340004}};
  for (int i = 0; i < 12; i++) {
    iPoint.atomID = i;
    iPoint.x = coords[i][0];
    iPoint.y = coords[i][1];
    iPoint.z = coords[i][2];
    yCloud.pts.push_back(iPoint);
    yCloud.idIndexMap[i] = i;
  }
  yCloud.nop = 12;
  return yCloud;
}

} // namespace

TEST_CASE("cageAffiliation matches the greedy classification on mW",
          "[cage_affiliation]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 7));
  REQUIRE(six.size() == 8192);

  const auto affiliation = ring::cageAffiliation(six, idx);

  std::vector<ring::strucType> rt(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cl;
  auto hc = ring::findHC(six, rt, idx, cl);
  auto ddc = ring::findDDC(six, rt, hc, cl);

  // In cubic ice every six-ring is DDC-classified and none HC-classified;
  // the order-free predicates agree ring for ring
  REQUIRE(ddc.size() == six.size());
  REQUIRE(hc.empty());
  for (size_t i = 0; i < six.size(); i++) {
    REQUIRE(affiliation.ddc[i]);
    REQUIRE_FALSE(affiliation.hc[i]);
  }
}

TEST_CASE("cageAffiliation marks the rings of an isolated hexagonal cage",
          "[cage_affiliation]") {
  auto yCloud = hcCloud();
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 7));
  REQUIRE_FALSE(six.empty());

  const auto affiliation = ring::cageAffiliation(six, idx);

  std::vector<ring::strucType> rt(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cl;
  auto listHC = ring::findHC(six, rt, idx, cl);
  REQUIRE(cl.size() == 1);
  REQUIRE(cl[0].rings.size() == 5);

  // The affiliation covers every ring of the greedy cage
  for (const int r : cl[0].rings) {
    REQUIRE(affiliation.hc[r]);
  }
  // And no ring outside the greedy HC list is HC-affiliated here
  for (size_t i = 0; i < six.size(); i++) {
    const bool inGreedy =
        std::find(listHC.begin(), listHC.end(), static_cast<int>(i)) !=
        listHC.end();
    REQUIRE(affiliation.hc[i] == inGreedy);
  }
}

TEST_CASE("cageAffiliation is invariant under ring permutation",
          "[cage_affiliation]") {
  auto yCloud = hcCloud();
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 7));

  const auto reference = ring::cageAffiliation(six, idx);
  std::map<std::vector<int>, std::pair<bool, bool>> byKey;
  for (size_t i = 0; i < six.size(); i++) {
    byKey[canonicalKeyOf(six[i])] = {reference.hc[i], reference.ddc[i]};
  }

  auto shuffled = six;
  std::sort(shuffled.rbegin(), shuffled.rend());
  const auto permuted = ring::cageAffiliation(shuffled, idx);
  for (size_t i = 0; i < shuffled.size(); i++) {
    const auto expected = byKey.at(canonicalKeyOf(shuffled[i]));
    REQUIRE(permuted.hc[i] == expected.first);
    REQUIRE(permuted.ddc[i] == expected.second);
  }
}

TEST_CASE("AffiliationUpdater equals batch across the mW trajectory",
          "[cage_affiliation]") {
  ring::AffiliationUpdater updater;
  for (int frame = 1; frame <= 11; frame++) {
    molSys::PointCloud<molSys::Point<double>, double> yCloud;
    yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", frame, yCloud, 1);
    auto nList = nneigh::neighListO(3.5, yCloud, 1);
    auto idx = nneigh::neighbourListByIndex(yCloud, nList);
    auto six = sixMembered(primitive::ringNetwork(idx, 7));

    const auto &incremental = updater.update(six, idx);
    const auto batch = ring::cageAffiliation(six, idx);
    REQUIRE(incremental.hc == batch.hc);
    REQUIRE(incremental.ddc == batch.ddc);
    if (frame == 1) {
      REQUIRE(updater.lastReclassified() == static_cast<int>(six.size()));
    } else {
      // Thermal jitter never crosses the bond cutoff in this trajectory
      REQUIRE(updater.lastReclassified() == 0);
    }
  }
}

TEST_CASE("AffiliationUpdater tracks synthetic topology churn exactly",
          "[cage_affiliation]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);

  auto removeBond = [](std::vector<std::vector<int>> &rows, int a, int b) {
    auto drop = [](std::vector<int> &row, int value) {
      row.erase(std::remove(row.begin() + 1, row.end(), value), row.end());
    };
    drop(rows[a], b);
    drop(rows[b], a);
  };

  ring::AffiliationUpdater updater;
  uint64_t rng = 88172645463325252ULL;
  auto rows = idx;
  for (int step = 0; step < 4; step++) {
    if (step > 0) {
      // Remove three pseudo-random bonds per step, forcing real ring churn
      for (int cut = 0; cut < 3; cut++) {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        const int a = static_cast<int>(rng % rows.size());
        if (rows[a].size() > 2) {
          const int b = rows[a][1 + (rng >> 32) % (rows[a].size() - 1)];
          removeBond(rows, a, b);
        }
      }
    }
    auto six = sixMembered(primitive::ringNetwork(rows, 7));
    const auto &incremental = updater.update(six, rows);
    const auto batch = ring::cageAffiliation(six, rows);
    REQUIRE(incremental.hc == batch.hc);
    REQUIRE(incremental.ddc == batch.ddc);
    if (step > 0) {
      // The dirty closure stays local: far fewer rings than the network
      REQUIRE(updater.lastReclassified() < static_cast<int>(six.size()) / 2);
    }
  }
}

TEST_CASE("AffiliationUpdater falls back to batch on whole-graph churn",
          "[cage_affiliation]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 7));

  ring::AffiliationUpdater updater;
  (void)updater.update(six, idx);

  auto emptyRows = idx;
  for (auto &row : emptyRows) {
    if (!row.empty()) {
      row.resize(1);
    }
  }
  const auto &incremental = updater.update(six, emptyRows);
  const auto batch = ring::cageAffiliation(six, emptyRows);
  REQUIRE(incremental.hc == batch.hc);
  REQUIRE(incremental.ddc == batch.ddc);
  REQUIRE(updater.lastReclassified() == static_cast<int>(six.size()));
}

TEST_CASE("cage membership is not a function of the four-neighbour star",
          "[cage_affiliation]") {
  // Two configurations that agree exactly on the closed star of a vertex --
  // same induced graph on {v} and its four neighbours, same coordinates for
  // those five atoms -- but differ in the vertex's cage label. Any per-atom
  // score that reads only the star must therefore assign both the same
  // label, and cannot equal the topological classification. The second
  // configuration displaces second-shell atoms only, which changes no
  // coordinate and no edge inside the star.
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);

  const int v = 137; // an arbitrary bulk vertex
  REQUIRE(idx[v].size() == 5);
  std::vector<int> star(idx[v].begin(), idx[v].end()); // v plus neighbours
  std::vector<bool> inStar(yCloud.nop, false);
  for (const int a : star) {
    inStar[a] = true;
  }

  // Second shell: neighbours of v's neighbours, outside the star
  std::vector<int> secondShell;
  for (size_t k = 1; k < idx[v].size(); k++) {
    for (size_t m = 1; m < idx[idx[v][k]].size(); m++) {
      const int w = idx[idx[v][k]][m];
      if (!inStar[w]) {
        secondShell.push_back(w);
      }
    }
  }
  std::sort(secondShell.begin(), secondShell.end());
  secondShell.erase(std::unique(secondShell.begin(), secondShell.end()),
                    secondShell.end());
  REQUIRE_FALSE(secondShell.empty());

  auto labelOf = [](const molSys::PointCloud<molSys::Point<double>, double>
                        &cloud,
                    int vertex) {
    auto nl = nneigh::neighListO(3.5, cloud, 1);
    auto ix = nneigh::neighbourListByIndex(cloud, nl);
    auto rings = primitive::ringNetwork(ix, 7);
    std::vector<std::vector<int>> six;
    for (const auto &r : rings) {
      if (r.size() == 6) {
        six.push_back(r);
      }
    }
    bool ddc = false;
    if (!six.empty()) {
      const auto affiliation = ring::cageAffiliation(six, ix);
      for (size_t i = 0; i < six.size(); i++) {
        if (!affiliation.ddc[i]) {
          continue;
        }
        if (std::find(six[i].begin(), six[i].end(), vertex) != six[i].end()) {
          ddc = true;
          break;
        }
      }
    }
    return ddc;
  };

  REQUIRE(labelOf(yCloud, v)); // cubic in the pristine crystal

  // Displace every second-shell atom far away: outside the star, so the
  // star's coordinates and induced edges are untouched
  auto broken = yCloud;
  for (const int w : secondShell) {
    broken.pts[w].x += 20.0;
  }

  // The closed star is bit-identical between the two configurations
  auto brokenIdx = nneigh::neighbourListByIndex(
      broken, nneigh::neighListO(3.5, broken, 1));
  for (const int a : star) {
    REQUIRE(yCloud.pts[a].x == broken.pts[a].x);
    REQUIRE(yCloud.pts[a].y == broken.pts[a].y);
    REQUIRE(yCloud.pts[a].z == broken.pts[a].z);
  }
  // Induced star edges: neighbour rows restricted to star members agree
  for (const int a : star) {
    std::vector<int> before, after;
    for (size_t m = 1; m < idx[a].size(); m++) {
      if (inStar[idx[a][m]]) {
        before.push_back(idx[a][m]);
      }
    }
    for (size_t m = 1; m < brokenIdx[a].size(); m++) {
      if (inStar[brokenIdx[a][m]]) {
        after.push_back(brokenIdx[a][m]);
      }
    }
    std::sort(before.begin(), before.end());
    std::sort(after.begin(), after.end());
    REQUIRE(before == after);
  }

  // Yet the cage label of v flips
  REQUIRE_FALSE(labelOf(broken, v));
}

TEST_CASE("seeded affiliation: structural zero, consistency, recovery",
          "[cage_affiliation]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 7));
  REQUIRE(six.size() == 8192);

  // Consistency: identical graphs reproduce the plain affiliation atoms
  const auto both = ring::seededCageAffiliation(six, idx, six, idx);
  const auto plain = ring::cageAffiliation(six, idx);
  std::vector<bool> plainDdcAtoms(yCloud.nop, false);
  for (size_t i = 0; i < six.size(); i++) {
    if (plain.ddc[i]) {
      for (const int a : six[i]) {
        plainDdcAtoms[a] = true;
      }
    }
  }
  REQUIRE(both.ddc == plainDdcAtoms);

  // Structural zero: an empty strict graph accepts nothing, however much
  // the permissive graph affiliates
  std::vector<std::vector<int>> selfOnly(idx.size());
  for (size_t i = 0; i < idx.size(); i++) {
    selfOnly[i].push_back(idx[i].empty() ? static_cast<int>(i) : idx[i][0]);
  }
  std::vector<std::vector<int>> noRings;
  const auto zero = ring::seededCageAffiliation(noRings, selfOnly, six, idx);
  REQUIRE(std::none_of(zero.hc.begin(), zero.hc.end(),
                       [](bool b) { return b; }));
  REQUIRE(std::none_of(zero.ddc.begin(), zero.ddc.end(),
                       [](bool b) { return b; }));

  // Recovery: isolating atoms in the strict graph loses them (affiliation
  // survives single cuts, since each atom sits in many rings); seeding the
  // intact permissive graph recovers every one of them
  auto cut = idx;
  for (const int victim : {17, 911, 2048, 3333, 4000}) {
    for (size_t m = 1; m < cut[victim].size(); m++) {
      auto &row = cut[cut[victim][m]];
      row.erase(std::remove(row.begin() + 1, row.end(), victim), row.end());
    }
    cut[victim].resize(1);
  }
  auto sixCut = sixMembered(primitive::ringNetwork(cut, 7));
  const auto strictOnly = ring::cageAffiliation(sixCut, cut);
  std::vector<bool> strictAtoms(yCloud.nop, false);
  for (size_t i = 0; i < sixCut.size(); i++) {
    if (strictOnly.ddc[i]) {
      for (const int a : sixCut[i]) {
        strictAtoms[a] = true;
      }
    }
  }
  const int lost = static_cast<int>(
      std::count(strictAtoms.begin(), strictAtoms.end(), false));
  REQUIRE(lost > 0);
  const auto seeded =
      ring::seededCageAffiliation(sixCut, cut, six, idx);
  REQUIRE(seeded.ddc == plainDdcAtoms);
}
