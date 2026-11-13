#include <catch2/catch_test_macros.hpp>

#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <order_parameter.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>
#include <tum_device.hpp>
#include <tum_offload.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

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
            p.x = std::fmod(shift * a + f[0] * a - f[1] * a / 2 + ix * a,
                            nx * a);
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

TEST_CASE("seeded affiliation on an empty framed cloud is empty",
          "[cage_affiliation]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = 0;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  auto knn = nneigh::kNearestNeighbourList(cloud, 4, 5.5, 1, true);
  auto uni = nneigh::kNearestNeighbourList(cloud, 4, 5.5, 1, false);
  auto idxS = nneigh::neighbourListByIndex(cloud, knn);
  auto idxU = nneigh::neighbourListByIndex(cloud, uni);
  auto sixS = sixMembered(primitive::ringNetwork(idxS, 6));
  auto sixU = sixMembered(primitive::ringNetwork(idxU, 6));
  const auto seeded = ring::seededCageAffiliation(sixS, idxS, sixU, idxU);
  REQUIRE(knn.empty());
  REQUIRE(sixS.empty());
  REQUIRE(seeded.hc.empty());
  REQUIRE(seeded.ddc.empty());
}

TEST_CASE("cageAffiliation matches the greedy classification on mW",
          "[cage_affiliation]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 7));
  REQUIRE(six.size() == 8192);

  std::vector<ring::strucType> rt(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cl;
  auto hc = ring::findHC(six, rt, idx, cl);
  auto ddc = ring::findDDC(six, rt, hc, cl);

  // In cubic ice every six-ring is DDC-classified and none HC-classified;
  // the order-free predicates agree ring for ring
  REQUIRE(ddc.size() == six.size());
  REQUIRE(hc.empty());
#ifdef SEAMS_HAS_OPENMP
  const int maxThreads = omp_get_max_threads();
  for (const int threads : {1, 4}) {
    omp_set_num_threads(threads);
#endif
    const auto affiliation = ring::cageAffiliation(six, idx);
    for (size_t i = 0; i < six.size(); i++) {
      REQUIRE(affiliation.ddc[i]);
      REQUIRE_FALSE(affiliation.hc[i]);
    }
#ifdef SEAMS_HAS_OPENMP
  }
  omp_set_num_threads(maxThreads);
#endif
}

TEST_CASE("stackingPlanes marks HC basals and keeps them off DDC equatorials",
          "[cage_affiliation]") {
  auto yCloud = hcCloud();
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 7));
  const auto planes = ring::stackingPlanes(six, idx);
  REQUIRE(planes.basal.size() == six.size());
  REQUIRE(planes.equatorial.size() == six.size());
  int nBasal = 0;
  int nEq = 0;
  int nBoth = 0;
  for (size_t i = 0; i < six.size(); i++) {
    if (planes.basal[i]) {
      ++nBasal;
    }
    if (planes.equatorial[i]) {
      ++nEq;
    }
    if (planes.basal[i] && planes.equatorial[i]) {
      ++nBoth;
    }
  }
  REQUIRE(nBasal == 2);
  REQUIRE(nEq == 0);
  REQUIRE(nBoth == 0);
  const auto aff = ring::cageAffiliation(six, idx);
  for (size_t i = 0; i < six.size(); i++) {
    if (planes.basal[i]) {
      REQUIRE(aff.hc[i]);
    }
    if (planes.equatorial[i]) {
      REQUIRE(aff.ddc[i]);
    }
  }
  const auto tum =
      topoparam::tumLayerStack(yCloud, six, planes.basal, planes.equatorial, 2,
                               3.7);
  REQUIRE(std::count(tum.sequence.begin(), tum.sequence.end(), 'H') == 2);
  REQUIRE(tum.sequence.find('C') == std::string::npos);
  REQUIRE(tum.phiC == 0.0);
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

TEST_CASE("cageAffiliation marks hexagonal ice HC on every thread count",
          "[cage_affiliation]") {
  // Enough rings for the threaded sweeps
  auto yCloud = iceIh(8, 5, 5);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 6));
  REQUIRE(six.size() >= 1024);

  std::vector<ring::strucType> rt(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cl;
  auto hc = ring::findHC(six, rt, idx, cl);
  REQUIRE(ring::findDDC(six, rt, hc, cl).empty());
  REQUIRE(hc.size() == six.size());

#ifdef SEAMS_HAS_OPENMP
  const int maxThreads = omp_get_max_threads();
  for (const int threads : {1, 4}) {
    omp_set_num_threads(threads);
#endif
    const auto affiliation = ring::cageAffiliation(six, idx);
    for (size_t i = 0; i < six.size(); i++) {
      REQUIRE(affiliation.hc[i]);
      REQUIRE_FALSE(affiliation.ddc[i]);
    }
#ifdef SEAMS_HAS_OPENMP
  }
  omp_set_num_threads(maxThreads);
#endif
}

TEST_CASE("parallel cage affiliation matches the one-thread sweep",
          "[cage_affiliation]") {
#ifdef SEAMS_HAS_OPENMP
  auto yCloud = iceIh(8, 5, 5);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 6));
  REQUIRE(six.size() >= 1024);

  struct RestoreThreads {
    int saved = omp_get_max_threads();
    ~RestoreThreads() { omp_set_num_threads(saved); }
  } restore;
  const int parallelThreads = std::max(2, std::min(restore.saved, 4));
  auto at = [&](int threads, const std::vector<std::vector<int>> &rows) {
    omp_set_num_threads(threads);
    return ring::cageAffiliation(six, rows);
  };

  const auto serial = at(1, idx);
  const auto parallel = at(parallelThreads, idx);
  REQUIRE(serial.hc == parallel.hc);
  REQUIRE(serial.ddc == parallel.ddc);
  REQUIRE(std::all_of(serial.hc.begin(), serial.hc.end(),
                      [](bool on) { return on; }));
  REQUIRE(std::none_of(serial.ddc.begin(), serial.ddc.end(),
                       [](bool on) { return on; }));

  // Drop every bond that touches every fifth atom. The flags are no longer
  // uniform, so a sweep that marks the wrong union cannot hide in an
  // all-hexagonal crystal.
  auto damaged = idx;
  for (std::size_t i = 0; i < damaged.size(); i += 5) {
    if (damaged[i].empty()) {
      continue;
    }
    const int self = damaged[i][0];
    for (auto &row : damaged) {
      if (row.size() <= 1) {
        continue;
      }
      if (row[0] == self) {
        row.resize(1);
        continue;
      }
      row.erase(std::remove(row.begin() + 1, row.end(), self), row.end());
    }
  }
  const auto serialCut = at(1, damaged);
  const auto parallelCut = at(parallelThreads, damaged);
  REQUIRE(serialCut.hc == parallelCut.hc);
  REQUIRE(serialCut.ddc == parallelCut.ddc);
  const bool hcMixed =
      std::find(serialCut.hc.begin(), serialCut.hc.end(), false) !=
          serialCut.hc.end() &&
      std::find(serialCut.hc.begin(), serialCut.hc.end(), true) !=
          serialCut.hc.end();
  const bool ddcMixed =
      std::find(serialCut.ddc.begin(), serialCut.ddc.end(), false) !=
          serialCut.ddc.end() &&
      std::find(serialCut.ddc.begin(), serialCut.ddc.end(), true) !=
          serialCut.ddc.end();
  REQUIRE((hcMixed || ddcMixed));
  REQUIRE(serialCut.hc != serial.hc);
#else
  SUCCEED("the sweep is the one-thread path without OpenMP");
#endif
}

TEST_CASE("hop-bound path grows its ring, row and pair buffers on ice Ih",
          "[cage_affiliation][offload]") {
  auto yCloud = iceIh(4, 3, 3);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  const auto host = tum::hostCageCounts(idx);
  REQUIRE(host.nSix > 0);
  REQUIRE(host.nHcRings == host.nSix);
  tum::SpecializedCaps caps;
  caps.rowCap = 1;
  caps.ringCap = 1;
  caps.pairCap = 1;
  for (const bool target : {false, true}) {
    INFO("target " << target);
    caps.forceTarget = target;
    const auto dev = tum::specializedCageCounts(idx, caps);
    // The ring table overflows, then the rows, then the basal pairs
    REQUIRE(dev.retries == 3);
    REQUIRE_FALSE(dev.hostFallback);
    REQUIRE(dev.nSix == host.nSix);
    REQUIRE(dev.nHcRings == host.nHcRings);
    REQUIRE(dev.nDdcRings == host.nDdcRings);
    REQUIRE(dev.atomSix == host.atomSix);
    REQUIRE(dev.atomHc == host.atomHc);
    REQUIRE(dev.atomDdc == host.atomDdc);
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
  nneigh::SkinNeighborList skin(3.5, 2.0, 1);
  for (int frame = 1; frame <= 11; frame++) {
    molSys::PointCloud<molSys::Point<double>, double> yCloud;
    yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", frame, yCloud, 1);
    auto nList = skin.update(yCloud);
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
  // Exercise the restricted update on every step, not the batch fallback
  updater.setBatchFraction(1.0);
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

TEST_CASE("AffiliationUpdater matches batch when every neighbour row changes",
          "[cage_affiliation]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 7));

  ring::AffiliationUpdater updater;
  const auto primed = updater.update(six, idx);
  REQUIRE(primed.hc == ring::cageAffiliation(six, idx).hc);
  REQUIRE(primed.ddc == ring::cageAffiliation(six, idx).ddc);

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
  // The frame that crossed into the whole-frame sweep leaves the next
  // frame on that sweep, with the same row count.
  const auto &again = updater.update(six, emptyRows);
  const auto againBatch = ring::cageAffiliation(six, emptyRows);
  REQUIRE(again.hc == againBatch.hc);
  REQUIRE(again.ddc == againBatch.ddc);
}

TEST_CASE("AffiliationUpdater matches an unsorted ring list and an empty frame",
          "[cage_affiliation]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 7));

  ring::AffiliationUpdater updater;
  updater.setBatchFraction(1.0);
  const auto &first = updater.update(six, idx);
  const auto batch = ring::cageAffiliation(six, idx);
  REQUIRE(first.hc == batch.hc);
  REQUIRE(first.ddc == batch.ddc);

  auto flipped = six;
  std::reverse(flipped.begin(), flipped.end());
  const auto &second = updater.update(flipped, idx);
  const auto batchFlip = ring::cageAffiliation(flipped, idx);
  REQUIRE(second.hc == batchFlip.hc);
  REQUIRE(second.ddc == batchFlip.ddc);

  std::vector<std::vector<int>> none;
  const auto &cleared = updater.update(none, idx);
  REQUIRE(cleared.hc.empty());
  REQUIRE(cleared.ddc.empty());

  // A tiny batch fraction sends a real edit through the backoff, then the
  // following frame skips the diff.
  ring::AffiliationUpdater backoff;
  backoff.setBatchFraction(0.001);
  const auto primedBackoff = backoff.update(six, idx);
  REQUIRE(primedBackoff.hc == ring::cageAffiliation(six, idx).hc);
  REQUIRE(primedBackoff.ddc == ring::cageAffiliation(six, idx).ddc);
  auto cut = idx;
  if (cut.size() > 2 && cut[0].size() > 1) {
    const int other = cut[0][1];
    cut[0].resize(1);
    cut[static_cast<std::size_t>(other)].erase(
        std::remove(cut[static_cast<std::size_t>(other)].begin(),
                    cut[static_cast<std::size_t>(other)].end(), 0),
        cut[static_cast<std::size_t>(other)].end());
  }
  auto sixCut = sixMembered(primitive::ringNetwork(cut, 7));
  const auto &churn = backoff.update(sixCut, cut);
  const auto churnBatch = ring::cageAffiliation(sixCut, cut);
  REQUIRE(churn.hc == churnBatch.hc);
  REQUIRE(churn.ddc == churnBatch.ddc);
  const auto &skipped = backoff.update(six, idx);
  const auto skipBatch = ring::cageAffiliation(six, idx);
  REQUIRE(skipped.hc == skipBatch.hc);
  REQUIRE(skipped.ddc == skipBatch.ddc);
  // The skipped frame stored no pair answers. The same graph again
  // rebuilds them instead of carrying a diff.
  const auto &rebuilt = backoff.update(six, idx);
  const auto rebuiltBatch = ring::cageAffiliation(six, idx);
  REQUIRE(rebuilt.hc == rebuiltBatch.hc);
  REQUIRE(rebuilt.ddc == rebuiltBatch.ddc);
}

TEST_CASE("AffiliationUpdater falls back when the ring set churns",
          "[cage_affiliation]") {
  // The neighbour rows stay put, so the changed-row check does not return.
  // Every new ring is dirty. Past the batch fraction that has to call
  // fellBack before returning, or the following frames never store a diff.
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 7));
  REQUIRE(six.size() >= 8);
  REQUIRE(idx.size() >= 6);

  ring::AffiliationUpdater churn;
  churn.setBatchFraction(0.05);
  const auto &primed = churn.update(six, idx);
  REQUIRE(primed.hc == ring::cageAffiliation(six, idx).hc);

  const int nAtoms = static_cast<int>(idx.size());
  std::vector<std::vector<int>> fresh;
  fresh.reserve(six.size());
  for (int i = 0; i < static_cast<int>(six.size()); ++i) {
    std::vector<int> ring(6);
    for (int k = 0; k < 6; ++k) {
      ring[static_cast<std::size_t>(k)] = (i * 6 + k) % nAtoms;
    }
    fresh.push_back(std::move(ring));
  }
  const auto expect = ring::cageAffiliation(fresh, idx);
  const auto &swept = churn.update(fresh, idx);
  REQUIRE(swept.hc == expect.hc);
  REQUIRE(swept.ddc == expect.ddc);
  REQUIRE(churn.lastReclassified() == static_cast<int>(fresh.size()));

  // The sweep does not store pairs. The next frame rebuilds them, and the
  // frame after that carries. Later frames are not forced onto the sweep.
  churn.update(fresh, idx);
  REQUIRE(churn.lastReclassified() == static_cast<int>(fresh.size()));
  churn.update(fresh, idx);
  REQUIRE(churn.lastReclassified() == 0);
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

TEST_CASE("ring-adjacent completion: identity on a crystal, empty on nothing",
          "[cage_affiliation]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 7));
  REQUIRE(six.size() == 8192);

  // Every atom of the perfect crystal is already accepted; completion is
  // the identity there
  const auto plain = ring::seededCageAffiliation(six, idx, six, idx, false);
  const auto completed = ring::seededCageAffiliation(six, idx, six, idx, true);
  REQUIRE(plain.hc == completed.hc);
  REQUIRE(plain.ddc == completed.ddc);
  REQUIRE(std::count(completed.ddc.begin(), completed.ddc.end(), true) ==
          yCloud.nop);

  // Remove one atom's labels by hand: every six-ring through it now has
  // five labelled vertices, and the completion restores the sixth
  ring::SeededAtomLabels holed = plain;
  holed.ddc[137] = false;
  const auto filled = ring::ringAdjacentCompletion(holed, six);
  REQUIRE(filled.ddc[137]);
  REQUIRE(filled.ddc == plain.ddc);

  // With no labelled vertex no ring has five: structural zero
  ring::SeededAtomLabels empty;
  empty.hc.assign(yCloud.nop, false);
  empty.ddc.assign(yCloud.nop, false);
  const auto still = ring::ringAdjacentCompletion(empty, six);
  REQUIRE(std::count(still.hc.begin(), still.hc.end(), true) == 0);
  REQUIRE(std::count(still.ddc.begin(), still.ddc.end(), true) == 0);
}

namespace {

// Ring E = {0..5} and peripherals on its triplets T_k = (k, k+1, k+2).
// Without the extra ring the only choice of P_0..P_5 that meets every
// equatorial condition puts one ring in two slots: Q holds T_0 and T_1
// in the first set, Q2 holds T_0, T_1 and T_2 in the second. The extra
// ring is a distinct witness for the repeated slot.
struct DdcFixture {
  const char *name;
  std::vector<std::vector<int>> rings;
  std::vector<int> extra;
};

std::vector<DdcFixture> repeatedWitnessFixtures() {
  return {
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
}

// tum::device::ddcFrom over every ring of an explicit ring list
std::vector<int> deviceDdc(const std::vector<std::vector<int>> &rings,
                           int nAtoms) {
  const int nRings = static_cast<int>(rings.size());
  const int maxPer = tum::device::kRowCap;
  std::vector<int> ringAtoms;
  for (const auto &r : rings) {
    ringAtoms.insert(ringAtoms.end(), r.begin(), r.end());
  }
  std::vector<int> count(nAtoms, 0);
  std::vector<int> through(static_cast<std::size_t>(nAtoms) * maxPer, -1);
  for (int r = 0; r < nRings; ++r) {
    for (const int a : rings[r]) {
      through[a * maxPer + count[a]++] = r;
    }
  }
  std::vector<int> hc(nRings, 0);
  std::vector<int> ddc(nRings, 0);
  for (int r = 0; r < nRings; ++r) {
    tum::device::ddcFrom(r, &nRings, ringAtoms.data(), count.data(),
                         through.data(), hc.data(), nAtoms, maxPer,
                         ddc.data());
  }
  return ddc;
}

} // namespace

TEST_CASE("DDC needs six distinct peripheral rings", "[cage_affiliation]") {
  for (const auto &fx : repeatedWitnessFixtures()) {
    INFO(fx.name);
    for (const bool withExtra : {false, true}) {
      INFO("extra ring " << withExtra);
      auto rings = fx.rings;
      if (withExtra) {
        rings.push_back(fx.extra);
      }
      int nAtoms = 0;
      for (const auto &r : rings) {
        nAtoms = std::max(nAtoms, *std::max_element(r.begin(), r.end()) + 1);
      }
      // No bonds: no basal pair, so no HC ring
      std::vector<std::vector<int>> nList(nAtoms);
      for (int a = 0; a < nAtoms; ++a) {
        nList[a] = {a};
      }
      const auto index = ring::buildRingSearchIndex(rings, nAtoms);
      // Series 05 keeps a peripheral set only when the seven rings span
      // 14 molecules and each peripheral shares exactly three equatorial
      // atoms. These witnesses share four, so the host rejects both.
      REQUIRE_FALSE(ring::equatorialConditions(rings, 0, index, nullptr, nullptr,
                                              nullptr));

      const auto aff = ring::cageAffiliation(rings, nList);
      REQUIRE(aff.hc.size() == rings.size());
      REQUIRE(aff.ddc.size() == rings.size());
      for (std::size_t r = 0; r < rings.size(); ++r) {
        INFO("ring " << r);
        REQUIRE_FALSE(aff.hc[r]);
        REQUIRE_FALSE(aff.ddc[r]);
      }

      std::vector<ring::strucType> types(rings.size(),
                                         ring::strucType::unclassified);
      std::vector<cage::Cage> cages;
      const auto listDDC = ring::findDDC(rings, types, {}, cages);
      REQUIRE(std::count(listDDC.begin(), listDDC.end(), 0) == 0);
      REQUIRE(cages.empty());
    }
  }
}

TEST_CASE("AffiliationUpdater matches a forced HC/DDC edge flip with fallback off",
          "[cage_affiliation]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/nucleation_f1.lammpstrj", 1, yCloud, 1);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto rows = nneigh::neighbourListByIndex(yCloud, nList);
  auto six = sixMembered(primitive::ringNetwork(rows, 6));
  REQUIRE_FALSE(six.empty());

  ring::AffiliationUpdater updater;
  updater.setBatchFraction(1.0);
  const auto &first = updater.update(six, rows);
  const auto batch0 = ring::cageAffiliation(six, rows);
  REQUIRE(first.hc == batch0.hc);
  REQUIRE(first.ddc == batch0.ddc);
  const int nHc =
      static_cast<int>(std::count(batch0.hc.begin(), batch0.hc.end(), true));
  const int nDdc =
      static_cast<int>(std::count(batch0.ddc.begin(), batch0.ddc.end(), true));
  REQUIRE(nHc > 0);
  REQUIRE(nDdc > 0);

  auto dropBond = [](std::vector<std::vector<int>> &graph, int a, int b) {
    auto drop = [](std::vector<int> &row, int value) {
      row.erase(std::remove(row.begin() + 1, row.end(), value), row.end());
    };
    drop(graph[a], b);
    drop(graph[b], a);
  };
  auto addBond = [](std::vector<std::vector<int>> &graph, int a, int b) {
    if (std::find(graph[a].begin(), graph[a].end(), b) == graph[a].end()) {
      graph[a].push_back(b);
    }
    if (std::find(graph[b].begin(), graph[b].end(), a) == graph[b].end()) {
      graph[b].push_back(a);
    }
  };

  bool sawSwitch = false;
  const auto baseRows = rows;
  const auto baseHc = batch0.hc;
  const auto baseDdc = batch0.ddc;
  const int nAtoms = static_cast<int>(rows.size());
  for (int a = 0; a < nAtoms && !sawSwitch; ++a) {
    if (rows[a].size() < 3) {
      continue;
    }
    const int b = rows[a][1];
    int c = -1;
    for (int cand = 0; cand < nAtoms; ++cand) {
      if (cand == a || cand == b) {
        continue;
      }
      if (std::find(rows[a].begin(), rows[a].end(), cand) == rows[a].end()) {
        c = cand;
        break;
      }
    }
    if (c < 0) {
      continue;
    }
    auto flipped = baseRows;
    dropBond(flipped, a, b);
    addBond(flipped, a, c);
    auto sixFlip = sixMembered(primitive::ringNetwork(flipped, 6));
    ring::AffiliationUpdater step;
    step.setBatchFraction(1.0);
    const auto primedStep = step.update(six, baseRows);
    REQUIRE(primedStep.hc == ring::cageAffiliation(six, baseRows).hc);
    REQUIRE(primedStep.ddc == ring::cageAffiliation(six, baseRows).ddc);
    const auto &inc = step.update(sixFlip, flipped);
    const auto batch = ring::cageAffiliation(sixFlip, flipped);
    REQUIRE(inc.hc == batch.hc);
    REQUIRE(inc.ddc == batch.ddc);
    REQUIRE(step.lastReclassified() > 0);
    REQUIRE(step.lastReclassified() < static_cast<int>(sixFlip.size()));
    const bool hcChanged = inc.hc != baseHc;
    const bool ddcChanged = inc.ddc != baseDdc;
    if (hcChanged || ddcChanged) {
      sawSwitch = true;
    }
  }
  REQUIRE(sawSwitch);
}

TEST_CASE("AffiliationUpdater covers a moved owner, a short batch, and a dropped ring",
          "[cage_affiliation]") {
  ring::AffiliationUpdater donor;
  ring::AffiliationUpdater moved(std::move(donor));
  ring::AffiliationUpdater assigned;
  assigned = std::move(moved);

  std::vector<std::vector<int>> rows(8);
  for (int a = 0; a < 8; ++a) {
    rows[a] = {a};
  }
  std::vector<std::vector<int>> rings = {
      {0, 1, 2, 3, 4, 5},
      {0, 1, 2, 3, 4, 6},
      {0, 1, 2, 3, 4, 7},
      {1, 2, 3, 4, 5, 6},
  };
  ring::AffiliationUpdater tight;
  tight.setBatchFraction(0.4);
  const auto &primed = tight.update(rings, rows);
  const auto primedBatch = ring::cageAffiliation(rings, rows);
  REQUIRE(primed.hc == primedBatch.hc);
  rows[0].push_back(1);
  const auto &over = tight.update(rings, rows);
  const auto overBatch = ring::cageAffiliation(rings, rows);
  REQUIRE(over.hc == overBatch.hc);
  const auto &hot = tight.update(rings, rows);
  const auto hotBatch = ring::cageAffiliation(rings, rows);
  REQUIRE(hot.hc == hotBatch.hc);

  // Equal spans, so the first sampled atom stays the pivot and the
  // negative index is an arm. The span and the arm stamp both see it.
  std::vector<std::vector<int>> tied(6);
  for (int a = 0; a < 6; ++a) {
    tied[a] = {a};
  }
  const auto tiedAff =
      ring::cageAffiliation(std::vector<std::vector<int>>{{0, 1, -1, 2, 3, 4}}, tied);
  REQUIRE(tiedAff.hc.size() == 1);

  // Every basal row names atom 6, so the pivot walks the ring that holds
  // a negative index. A five-ring in the same list takes the other walk.
  std::vector<std::vector<int>> wide(13);
  for (int a = 0; a < 13; ++a) {
    wide[a] = {a, 6, 11};
  }
  std::vector<std::vector<int>> mixed = {
      {0, 1, 2, 3, 4, 5},
      {-1, 6, 7, 8, 9, 10},
  };
  const auto mixedAff = ring::cageAffiliation(mixed, wide);
  REQUIRE(mixedAff.hc.size() == mixed.size());
  std::vector<std::vector<int>> withFive = mixed;
  withFive.push_back({0, 1, 2, 100, 11});
  const auto fiveAff = ring::cageAffiliation(withFive, wide);
  REQUIRE(fiveAff.hc.size() == withFive.size());
  const auto fivePlanes = ring::stackingPlanes(withFive, wide);
  REQUIRE(fivePlanes.basal.size() == withFive.size());
  const auto seeded =
      ring::seededCageAffiliation(withFive, wide, withFive, wide, false);
  REQUIRE(seeded.hc.size() == wide.size());
  const auto noRings = ring::stackingPlanes({}, wide);
  REQUIRE(noRings.basal.empty());
  ring::SeededAtomLabels labels;
  labels.hc.assign(3, false);
  labels.ddc.assign(3, false);
  const auto completed =
      ring::ringAdjacentCompletion(labels, std::vector<std::vector<int>>{{0, 1, -1, 100}});
  REQUIRE(completed.hc.size() == 3);

  ring::AffiliationUpdater negative;
  negative.setBatchFraction(1.0);
  const auto primedNeg = negative.update(mixed, wide);
  REQUIRE(primedNeg.hc == ring::cageAffiliation(mixed, wide).hc);
  REQUIRE(primedNeg.ddc == ring::cageAffiliation(mixed, wide).ddc);
  auto withNegative = mixed;
  withNegative[0].push_back(-1);
  withNegative[0].pop_back();
  withNegative[0][0] = -1;
  const auto &marked = negative.update(withNegative, wide);
  REQUIRE(marked.hc.size() == withNegative.size());

  std::vector<std::vector<int>> fpRows(9);
  for (int a = 0; a < 9; ++a) {
    fpRows[a] = {a};
  }
  std::vector<std::vector<int>> firstRings = {
      {1, 2, 3, 4, 5, 0},
      {},
      {0, 1, 2, 3, 4, 6},
  };
  std::vector<std::vector<int>> secondRings = {
      {2, 3, 4, 5, 7, 8},
      {},
      {0, 1, 2, 3, 4, 5},
  };
  ring::AffiliationUpdater finger;
  finger.setBatchFraction(1.0);
  const auto primedFinger = finger.update(firstRings, fpRows);
  REQUIRE(primedFinger.hc == ring::cageAffiliation(firstRings, fpRows).hc);
  REQUIRE(primedFinger.ddc == ring::cageAffiliation(firstRings, fpRows).ddc);
  const auto &fingered = finger.update(secondRings, fpRows);
  const auto fingerBatch = ring::cageAffiliation(secondRings, fpRows);
  REQUIRE(fingered.hc == fingerBatch.hc);

  auto cloud = iceIh(4, 3, 3);
  auto nList = nneigh::neighListO(3.5, cloud, 1);
  auto idx = nneigh::neighbourListByIndex(cloud, nList);
  auto six = sixMembered(primitive::ringNetwork(idx, 6));
  REQUIRE(six.size() > 2);
  ring::AffiliationUpdater drop;
  drop.setBatchFraction(1.0);
  const auto primedDrop = drop.update(six, idx);
  REQUIRE(primedDrop.hc == ring::cageAffiliation(six, idx).hc);
  REQUIRE(primedDrop.ddc == ring::cageAffiliation(six, idx).ddc);
  auto shorter = six;
  const auto planes = ring::stackingPlanes(six, idx);
  std::size_t dropAt = shorter.size();
  for (std::size_t i = 0; i < planes.basal.size(); ++i) {
    if (planes.basal[i]) {
      dropAt = i;
    }
  }
  REQUIRE(dropAt < shorter.size());
  shorter.erase(shorter.begin() + static_cast<std::ptrdiff_t>(dropAt));
  const auto &kept = drop.update(shorter, idx);
  const auto keptBatch = ring::cageAffiliation(shorter, idx);
  REQUIRE(kept.hc == keptBatch.hc);
  REQUIRE(kept.ddc == keptBatch.ddc);

  std::vector<std::size_t> basals;
  for (std::size_t i = 0; i < planes.basal.size(); ++i) {
    if (planes.basal[i]) {
      basals.push_back(i);
    }
  }
  for (const std::size_t at : {basals.front(), basals[basals.size() / 2]}) {
    ring::AffiliationUpdater again;
    again.setBatchFraction(1.0);
    const auto primedAgain = again.update(six, idx);
    REQUIRE(primedAgain.hc == ring::cageAffiliation(six, idx).hc);
    REQUIRE(primedAgain.ddc == ring::cageAffiliation(six, idx).ddc);
    auto cut = six;
    cut.erase(cut.begin() + static_cast<std::ptrdiff_t>(at));
    const auto &part = again.update(cut, idx);
    const auto partBatch = ring::cageAffiliation(cut, idx);
    REQUIRE(part.hc == partBatch.hc);
    REQUIRE(part.ddc == partBatch.ddc);
  }
}

TEST_CASE("Affiliation edges the sweep used to leave closed",
          "[cage_affiliation]") {
  // The lesser neighbour of the least atom sits behind it, so the cycle
  // hash walks backward. The updater must still match the forward writing.
  std::vector<std::vector<int>> rows(6);
  for (int a = 0; a < 6; ++a) {
    rows[a] = {a};
  }
  const std::vector<std::vector<int>> backward = {{1, 0, 3, 4, 5, 2}};
  const std::vector<std::vector<int>> forward = {{0, 2, 5, 4, 3, 1}};
  ring::AffiliationUpdater walk;
  walk.setBatchFraction(1.0);
  const auto &backAff = walk.update(backward, rows);
  const auto backBatch = ring::cageAffiliation(backward, rows);
  REQUIRE(backAff.hc == backBatch.hc);
  REQUIRE(backAff.ddc == backBatch.ddc);
  const auto &fwdAff = walk.update(forward, rows);
  const auto fwdBatch = ring::cageAffiliation(forward, rows);
  REQUIRE(fwdAff.hc == fwdBatch.hc);
  REQUIRE(fwdAff.hc == backAff.hc);

  // These two triples are not the same cycle, and they share a hash.
  // The matcher has to reject the cross pair.
  const std::vector<std::vector<int>> collided = {
      {1, 1277, 8208187},
      {0, 1, 2},
  };
  std::vector<std::vector<int>> wide(8208188);
  for (int a : {0, 1, 2, 1277, 8208187}) {
    wide[static_cast<std::size_t>(a)] = {a};
  }
  ring::AffiliationUpdater collide;
  collide.setBatchFraction(1.0);
  const auto &first = collide.update(collided, wide);
  const auto firstBatch = ring::cageAffiliation(collided, wide);
  REQUIRE(first.hc.size() == 2);
  REQUIRE(first.hc == firstBatch.hc);
  const std::vector<std::vector<int>> collidedAgain = {
      {0, 1, 2},
      {1, 1277, 8208187},
  };
  const auto &second = collide.update(collidedAgain, wide);
  const auto secondBatch = ring::cageAffiliation(collidedAgain, wide);
  REQUIRE(second.hc == secondBatch.hc);
  REQUIRE(second.ddc == secondBatch.ddc);

  // Enough generations that the 8-bit stamp counter wraps. The rings
  // are disjoint, so the answer is empty and a second pass matches it.
  constexpr int nWrap = 10000;
  std::vector<std::vector<int>> wrapRows(static_cast<std::size_t>(nWrap * 6));
  std::vector<std::vector<int>> many(static_cast<std::size_t>(nWrap));
  for (int i = 0; i < nWrap; ++i) {
    std::vector<int> ring;
    for (int k = 0; k < 6; ++k) {
      const int atom = i * 6 + k;
      ring.push_back(atom);
      wrapRows[static_cast<std::size_t>(atom)] = {atom};
    }
    many[static_cast<std::size_t>(i)] = std::move(ring);
  }
  const auto wrapped = ring::cageAffiliation(many, wrapRows);
  REQUIRE(wrapped.hc.size() == many.size());
  REQUIRE(std::count(wrapped.hc.begin(), wrapped.hc.end(), true) == 0);
  const auto wrappedAgain = ring::cageAffiliation(many, wrapRows);
  REQUIRE(wrappedAgain.hc == wrapped.hc);
  REQUIRE(wrappedAgain.ddc == wrapped.ddc);
  const auto small = ring::cageAffiliation(forward, rows);
  REQUIRE(small.hc == fwdBatch.hc);

  // First use of a ring fills its basal list. The planes and the
  // affiliation describe the same rings.
  const auto planes = ring::stackingPlanes(forward, rows);
  REQUIRE(planes.basal.size() == forward.size());
  for (std::size_t i = 0; i < planes.basal.size(); ++i) {
    if (planes.basal[i]) {
      REQUIRE(small.hc[i]);
    }
  }

  // A neighbour id outside the cloud is dropped by the flood and does
  // not change which atoms keep the seed.
  ring::CageAffiliation strict;
  strict.hc = {true};
  strict.ddc = {false};
  ring::CageAffiliation permissive = strict;
  std::vector<std::vector<int>> nList(3);
  for (int a = 0; a < 3; ++a) {
    nList[static_cast<std::size_t>(a)] = {a};
  }
  nList[0].push_back(-1);
  nList[0].push_back(40);
  nList[0].push_back(1);
  const std::vector<std::vector<int>> one = {{0, 1, 2}};
  const auto flooded =
      ring::seededCageAffiliation(one, strict, one, permissive, nList, false);
  std::vector<std::vector<int>> clean = nList;
  clean[0] = {0, 1};
  const auto plain =
      ring::seededCageAffiliation(one, strict, one, permissive, clean, false);
  REQUIRE(flooded.hc == plain.hc);
  REQUIRE(flooded.ddc == plain.ddc);
  REQUIRE(flooded.hc[0]);
  REQUIRE(flooded.hc[1]);
}

TEST_CASE("Affiliation reads a backward cycle and an index outside the cloud",
          "[cage_affiliation]") {
  // Two rings, out of order, so the matcher hashes. The least atom of the
  // first ring has its lesser neighbour behind it. The next frame writes
  // that same cycle toward the lesser neighbour.
  std::vector<std::vector<int>> rows(8);
  for (int a = 0; a < 8; ++a) {
    rows[static_cast<std::size_t>(a)] = {a};
  }
  const std::vector<std::vector<int>> backward = {
      {1, 0, 3, 4, 5, 2},
      {0, 1, 2, 3, 6, 7},
  };
  const std::vector<std::vector<int>> forward = {
      {0, 1, 2, 5, 4, 3},
      {0, 1, 2, 3, 6, 7},
  };
  ring::AffiliationUpdater walk;
  walk.setBatchFraction(1.0);
  const auto &back = walk.update(backward, rows);
  REQUIRE(back.hc == ring::cageAffiliation(backward, rows).hc);
  REQUIRE(back.ddc == ring::cageAffiliation(backward, rows).ddc);
  const auto &fwd = walk.update(forward, rows);
  REQUIRE(fwd.hc == ring::cageAffiliation(forward, rows).hc);
  REQUIRE(fwd.ddc == ring::cageAffiliation(forward, rows).ddc);

  // A six-ring atom past the cloud, a negative ring atom, and neighbour
  // ids that are negative or past the cloud. Atom 5 is in the cloud and
  // past every atom the rings actually contain, except the rejected one.
  std::vector<std::vector<int>> wide(6);
  for (int a = 0; a < 6; ++a) {
    wide[static_cast<std::size_t>(a)] = {a, -1, 50, 5};
  }
  const std::vector<std::vector<int>> odd = {
      {0, 1, 2, 3, 4, 100},
      {-4, 0, 1, 2, 3, 4},
  };
  const auto oddAff = ring::cageAffiliation(odd, wide);
  REQUIRE(oddAff.hc.size() == odd.size());
  REQUIRE(oddAff.ddc.size() == odd.size());
  REQUIRE(ring::stackingPlanes(odd, wide).basal.size() == odd.size());
  const auto oddAgain = ring::cageAffiliation(odd, wide);
  REQUIRE(oddAgain.hc == oddAff.hc);
  REQUIRE(oddAgain.ddc == oddAff.ddc);
  // Neighbour 5 is in this cloud and past every atom the ring lists.
  std::vector<std::vector<int>> pastRows(6);
  for (int a = 0; a < 6; ++a) {
    pastRows[static_cast<std::size_t>(a)] = {a, 5};
  }
  const auto pastRing =
      ring::cageAffiliation(std::vector<std::vector<int>>{{0, 1, 2, 3, 4, 0}}, pastRows);
  REQUIRE(pastRing.hc.size() == 1);

  // Drop a ring that shares no atom with the ring that stays. The rows
  // do not change, and they name a negative id and an id past the cloud.
  std::vector<std::vector<int>> cloud(12);
  for (int a = 0; a < 12; ++a) {
    cloud[static_cast<std::size_t>(a)] = {a, -2, 80};
  }
  const std::vector<std::vector<int>> both = {
      {0, 1, 2, 3, 4, 5},
      {6, 7, 8, 9, 10, 100},
  };
  ring::AffiliationUpdater drop;
  drop.setBatchFraction(1.0);
  REQUIRE(drop.update(both, cloud).hc == ring::cageAffiliation(both, cloud).hc);
  const std::vector<std::vector<int>> left = {{0, 1, 2, 3, 4, 5}};
  const auto &kept = drop.update(left, cloud);
  REQUIRE(kept.hc == ring::cageAffiliation(left, cloud).hc);
  REQUIRE(kept.ddc == ring::cageAffiliation(left, cloud).ddc);

  // Any changed row crosses a zero fraction and the next cloud is longer.
  std::vector<std::vector<int>> hotRows(6);
  for (int a = 0; a < 6; ++a) {
    hotRows[static_cast<std::size_t>(a)] = {a};
  }
  const std::vector<std::vector<int>> hex = {{0, 1, 2, 3, 4, 5}};
  ring::AffiliationUpdater hot;
  hot.setBatchFraction(0.0);
  REQUIRE(hot.update(hex, hotRows).hc == ring::cageAffiliation(hex, hotRows).hc);
  auto changed = hotRows;
  changed[0].push_back(1);
  REQUIRE(hot.update(hex, changed).hc == ring::cageAffiliation(hex, changed).hc);
  changed.push_back({0});
  const auto &resized = hot.update(hex, changed);
  REQUIRE(resized.hc == ring::cageAffiliation(hex, changed).hc);
  REQUIRE(resized.ddc == ring::cageAffiliation(hex, changed).ddc);

  // Strict and permissive clouds label different atoms. An empty label
  // vector returns before the walk, and a ring atom past the cloud does not.
  ring::CageAffiliation strictOnly;
  strictOnly.hc = {true};
  strictOnly.ddc = {false};
  ring::CageAffiliation permOnly;
  permOnly.hc = {false};
  permOnly.ddc = {true};
  std::vector<std::vector<int>> bond(2);
  bond[0] = {0, -1, 40, 1};
  bond[1] = {1};
  const std::vector<std::vector<int>> pair = {{0, 1}};
  const auto mixed =
      ring::seededCageAffiliation(pair, strictOnly, pair, permOnly, bond, true);
  REQUIRE(mixed.hc.size() == bond.size());
  REQUIRE(mixed.hc[0]);
  REQUIRE(mixed.ddc.size() == bond.size());
  ring::SeededAtomLabels none;
  REQUIRE(ring::ringAdjacentCompletion(none, pair).hc.empty());
  ring::SeededAtomLabels oneAtom;
  oneAtom.hc = {false};
  oneAtom.ddc = {false};
  const auto noRing = ring::ringAdjacentCompletion(oneAtom, {});
  REQUIRE(noRing.hc.size() == 1);
  REQUIRE(!noRing.hc[0]);
  ring::SeededAtomLabels two;
  two.hc = {true, false};
  two.ddc = {false, false};
  const auto past = ring::ringAdjacentCompletion(two, {{0, 9, -3}});
  REQUIRE(past.hc.size() == 2);
  REQUIRE(past.hc[0]);
  REQUIRE(!past.hc[1]);
}

TEST_CASE("a sized warm start matches a fresh affiliation", "[cage_affiliation]") {
  std::vector<std::vector<int>> hotRows(6);
  for (int a = 0; a < 6; ++a) {
    hotRows[static_cast<std::size_t>(a)] = {a};
  }
  const std::vector<std::vector<int>> hex = {{0, 1, 2, 3, 4, 5}};
  ring::AffiliationUpdater sized;
  sized.setBatchFraction(1.0);
  const auto primedSized = sized.update(hex, hotRows);
  REQUIRE(primedSized.hc == ring::cageAffiliation(hex, hotRows).hc);
  REQUIRE(primedSized.ddc == ring::cageAffiliation(hex, hotRows).ddc);
  auto longer = hotRows;
  longer.push_back({1});
  REQUIRE(sized.update(hex, longer).hc == ring::cageAffiliation(hex, longer).hc);

#ifdef SEAMS_HAS_OPENMP
  std::vector<std::vector<int>> manyRows(1024);
  for (int a = 0; a < 1024; ++a) {
    manyRows[static_cast<std::size_t>(a)] = {a % 6};
  }
  ring::AffiliationUpdater inside;
  inside.setBatchFraction(1.0);
  const auto primedInside = inside.update(hex, manyRows);
  REQUIRE(primedInside.hc.size() == 1u);
  REQUIRE(primedInside.ddc.size() == 1u);
  int inner = 0;
#pragma omp parallel
  {
    if (omp_get_thread_num() == 0) {
      inner = static_cast<int>(inside.update(hex, manyRows).hc.size());
    }
  }
  REQUIRE(inner == 1);
#endif
}

TEST_CASE("Affiliation reads a pivot past the cloud and a one-sided seed",
          "[cage_affiliation]") {
  // The first sampled atom is past the cloud, so its span is empty and
  // it stays the pivot. The affiliation is still defined.
  std::vector<std::vector<int>> rows(6);
  for (int a = 0; a < 6; ++a) {
    rows[static_cast<std::size_t>(a)] = {a, -3, (a + 1) % 6};
  }
  const std::vector<std::vector<int>> pivot = {{100, 1, 2, 3, 4, 5}};
  const auto pivotAff = ring::cageAffiliation(pivot, rows);
  REQUIRE(pivotAff.hc.size() == 1);
  REQUIRE(pivotAff.ddc.size() == 1);
  REQUIRE(ring::stackingPlanes(pivot, rows).basal.size() == 1);

  // The ring that stays names an atom past the cloud. Dropping the other
  // ring, with the rows unchanged, must match a fresh sweep.
  std::vector<std::vector<int>> keptRows(6);
  for (int a = 0; a < 6; ++a) {
    keptRows[static_cast<std::size_t>(a)] = {a};
  }
  const std::vector<std::vector<int>> both = {
      {0, 1, 2, 3, 4, 100},
      {0, 1, 2, 3, 5, 4},
  };
  ring::AffiliationUpdater keep;
  keep.setBatchFraction(10.0);
  REQUIRE(keep.update(both, keptRows).hc == ring::cageAffiliation(both, keptRows).hc);
  const std::vector<std::vector<int>> stay = {{100, 0, 1, 2, 3, 4}};
  const auto &kept = keep.update(stay, keptRows);
  REQUIRE(kept.hc == ring::cageAffiliation(stay, keptRows).hc);
  REQUIRE(kept.ddc == ring::cageAffiliation(stay, keptRows).ddc);

  // {1,0,2} walks backward from 0. The next frame writes the same cycle
  // forward, and the lists are out of order so the matcher hashes them.
  std::vector<std::vector<int>> hashRows(9);
  for (int a = 0; a < 9; ++a) {
    hashRows[static_cast<std::size_t>(a)] = {a};
  }
  const std::vector<std::vector<int>> back = {
      {8, 7, 6, 5, 4, 3},
      {1, 0, 2},
  };
  const std::vector<std::vector<int>> fwd = {
      {3, 8, 7, 6, 5, 4},
      {0, 1, 2},
  };
  ring::AffiliationUpdater walk;
  walk.setBatchFraction(1.0);
  const auto &backAff = walk.update(back, hashRows);
  REQUIRE(backAff.hc == ring::cageAffiliation(back, hashRows).hc);
  const auto &fwdAff = walk.update(fwd, hashRows);
  REQUIRE(fwdAff.hc == ring::cageAffiliation(fwd, hashRows).hc);
  REQUIRE(fwdAff.ddc == ring::cageAffiliation(fwd, hashRows).ddc);

  // Same hash, same order, still unsorted. The second ring meets an
  // old ring the first ring already claimed.
  const std::vector<std::vector<int>> collided = {
      {1, 1277, 8208187},
      {0, 1, 2},
  };
  std::vector<std::vector<int>> hashWide(8208188);
  for (int a : {0, 1, 2, 1277, 8208187}) {
    hashWide[static_cast<std::size_t>(a)] = {a};
  }
  ring::AffiliationUpdater collide;
  collide.setBatchFraction(1.0);
  REQUIRE(collide.update(collided, hashWide).hc ==
          ring::cageAffiliation(collided, hashWide).hc);
  const auto &again = collide.update(collided, hashWide);
  REQUIRE(again.hc == ring::cageAffiliation(collided, hashWide).hc);
  REQUIRE(again.ddc == ring::cageAffiliation(collided, hashWide).ddc);

  // Atom 0 is a strict seed. Atom 1 is only in the permissive cloud.
  // The bond carries the seed onto atom 1.
  ring::CageAffiliation strictSeed;
  strictSeed.hc = {true};
  strictSeed.ddc = {true};
  ring::CageAffiliation permSeed;
  permSeed.hc = {true};
  permSeed.ddc = {true};
  std::vector<std::vector<int>> bond(3);
  bond[0] = {0, 1};
  bond[1] = {1, 0};
  bond[2] = {2};
  const auto sided = ring::seededCageAffiliation(
      {{0}}, strictSeed, {{1}}, permSeed, bond, false);
  REQUIRE(sided.hc.size() == 3);
  REQUIRE(sided.hc[0]);
  REQUIRE(sided.hc[1]);
  REQUIRE(!sided.hc[2]);
  REQUIRE(sided.ddc[0]);
  REQUIRE(sided.ddc[1]);
  REQUIRE(!sided.ddc[2]);

#ifdef SEAMS_HAS_OPENMP
  std::vector<std::vector<int>> many(1024);
  for (int a = 0; a < 1024; ++a) {
    many[static_cast<std::size_t>(a)] = {a % 6};
  }
  const std::vector<std::vector<int>> hex = {{0, 1, 2, 3, 4, 5}};
  int inner = 0;
  std::vector<std::vector<int>> wideRows(static_cast<std::size_t>(1024 * 6));
  std::vector<std::vector<int>> manyRings(1024);
  for (int i = 0; i < 1024; ++i) {
    std::vector<int> ring;
    for (int k = 0; k < 6; ++k) {
      const int atom = i * 6 + k;
      ring.push_back(atom);
      wideRows[static_cast<std::size_t>(atom)] = {atom};
    }
    manyRings[static_cast<std::size_t>(i)] = std::move(ring);
  }
  int wideCount = 0;
#pragma omp parallel num_threads(2)
  {
    if (omp_get_thread_num() == 0) {
      ring::AffiliationUpdater nested;
      nested.setBatchFraction(1.0);
      inner = static_cast<int>(nested.update(hex, many).hc.size());
      wideCount = static_cast<int>(ring::cageAffiliation(manyRings, wideRows).hc.size());
    }
  }
  REQUIRE(inner == 1);
  REQUIRE(wideCount == 1024);
  REQUIRE(ring::cageAffiliation(hex, many).hc.size() == 1);
#endif
}
