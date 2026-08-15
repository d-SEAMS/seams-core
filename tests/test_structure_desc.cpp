#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <ira_sofi.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <structure_desc.hpp>
#include <voronoi_qlm.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

Cloud lattice(const std::vector<std::array<double, 3>> &basis, int reps,
              double a) {
  Cloud cloud;
  int id = 1;
  for (int i = 0; i < reps; i++) {
    for (int j = 0; j < reps; j++) {
      for (int k = 0; k < reps; k++) {
        for (const auto &b : basis) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          p.x = (i + b[0]) * a;
          p.y = (j + b[1]) * a;
          p.z = (k + b[2]) * a;
          cloud.pts.push_back(p);
          cloud.idIndexMap[id] = id - 1;
          id++;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.currentFrame = 1;
  const double L = reps * a;
  cloud.box = {L, L, L};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

Cloud fcc() {
  return lattice(
      {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}}, 3,
      4.0);
}

Cloud bcc() {
  return lattice({{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, 3, 4.0);
}

Cloud hcp() {
  Cloud cloud;
  const double a = 4.0;
  const double c = a * std::sqrt(8.0 / 3.0);
  const double ly = a * std::sqrt(3.0);
  const std::array<std::array<double, 3>, 4> basis = {{
      {{0.0, 0.0, 0.0}},
      {{0.5 * a, 0.5 * ly, 0.0}},
      {{0.5 * a, ly / 6.0, 0.5 * c}},
      {{0.0, 2.0 * ly / 3.0, 0.5 * c}},
  }};
  int id = 1;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (const auto &b : basis) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          p.x = i * a + b[0];
          p.y = j * ly + b[1];
          p.z = k * c + b[2];
          cloud.pts.push_back(p);
          cloud.idIndexMap[id] = id - 1;
          id++;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.currentFrame = 1;
  cloud.box = {3.0 * a, 3.0 * ly, 3.0 * c};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

} // namespace

TEST_CASE("IRA/Horn templates assign FCC, HCP and BCC lattices",
          "[structure_desc]") {
  auto fccCloud = fcc();
  auto bccCloud = bcc();
  auto hcpCloud = hcp();
  auto fccN = nneigh::neighListO(3.2, fccCloud, 1);
  auto bccN = nneigh::neighListO(4.0, bccCloud, 1);
  auto hcpN = nneigh::neighListO(1.2 * 4.0, hcpCloud, 1);

  auto fccHit = chill::classifyTemplates(fccCloud, fccN, 12);
  auto bccHit = chill::classifyTemplates(bccCloud, bccN, 8);
  auto hcpHit = chill::classifyTemplates(hcpCloud, hcpN, 12);

  int fccOk = 0;
  int bccOk = 0;
  int hcpOk = 0;
  int hcpClose = 0;
  for (const auto &h : fccHit) {
    if (h.kind == chill::CrystalKind::fcc && h.rmsd < 0.2) {
      fccOk++;
    }
  }
  for (const auto &h : bccHit) {
    if (h.kind == chill::CrystalKind::bcc && h.rmsd < 0.2) {
      bccOk++;
    }
  }
  for (const auto &h : hcpHit) {
    if (h.kind == chill::CrystalKind::hcp && h.rmsd < 0.2) {
      hcpOk++;
    }
    if ((h.kind == chill::CrystalKind::hcp ||
         h.kind == chill::CrystalKind::fcc) &&
        h.rmsd < 0.2) {
      hcpClose++;
    }
  }
  REQUIRE(fccOk > fccCloud.nop / 2);
  REQUIRE(bccOk > bccCloud.nop / 2);
  if (ira::available()) {
    REQUIRE(hcpOk > hcpCloud.nop / 2);
  } else {
    REQUIRE(hcpClose > hcpCloud.nop / 2);
  }
}

TEST_CASE("SOAP spectrum is finite and rotationally stable on FCC",
          "[structure_desc]") {
  auto cloud = fcc();
  auto nList = nneigh::neighListO(3.2, cloud, 1);
  auto a = chill::soapSpectrum(cloud, 0, nList, 3, 6, 3.2);
  auto b = chill::soapSpectrum(cloud, 1, nList, 3, 6, 3.2);
  REQUIRE(a.size() == 3 * 3 * 7);
  double na = 0.0;
  double nb = 0.0;
  double dot = 0.0;
  for (size_t i = 0; i < a.size(); i++) {
    REQUIRE(std::isfinite(a[i]));
    na += a[i] * a[i];
    nb += b[i] * b[i];
    dot += a[i] * b[i];
  }
  REQUIRE(na > 0.0);
  REQUIRE(dot / std::sqrt(na * nb) > 0.99);
}

TEST_CASE("linear classifier separates FCC from BCC on Voronoi features",
          "[structure_desc]") {
  auto fccCloud = fcc();
  auto bccCloud = bcc();
  auto pack = [](const Cloud &cloud, double cut) {
    const auto q4 = chill::steinhardtQlVoronoi(cloud, cut, 4);
    const auto q6 = chill::steinhardtQlVoronoi(cloud, cut, 6);
    const auto q8 = chill::steinhardtQlVoronoi(cloud, cut, 8);
    std::vector<std::vector<double>> rows;
    for (int i = 0; i < cloud.nop; i++) {
      rows.push_back({q4.ql[static_cast<size_t>(i)],
                      q6.ql[static_cast<size_t>(i)],
                      q8.ql[static_cast<size_t>(i)]});
    }
    return rows;
  };
  auto fccX = pack(fccCloud, 4.8);
  auto bccX = pack(bccCloud, 5.6);
  std::vector<std::vector<double>> X = fccX;
  std::vector<int> y(fccX.size(), 0);
  X.insert(X.end(), bccX.begin(), bccX.end());
  y.insert(y.end(), bccX.size(), 1);
  chill::LinearClassifier clf;
  clf.labels = {"fcc", "bcc"};
  clf.fit(X, y);

  int fccRight = 0;
  int bccRight = 0;
  for (const auto &row : fccX) {
    if (clf.predict(row) == 0) {
      fccRight++;
    }
  }
  for (const auto &row : bccX) {
    if (clf.predict(row) == 1) {
      bccRight++;
    }
  }
  REQUIRE(fccRight == fccCloud.nop);
  REQUIRE(bccRight == bccCloud.nop);
}

TEST_CASE("soapSpectrumAll matches soapSpectrum for atom 0",
          "[structure_desc]") {
  auto cloud = fcc();
  auto nList = nneigh::neighListO(3.2, cloud, 1);
  auto one = chill::soapSpectrum(cloud, 0, nList, 3, 6, 3.2);
  auto all = chill::soapSpectrumAll(cloud, nList, 3, 6, 3.2);
  REQUIRE(all.size() == static_cast<size_t>(cloud.nop));
  REQUIRE(all[0].size() == one.size());
  for (size_t i = 0; i < one.size(); i++) {
    REQUIRE_THAT(all[0][i], Catch::Matchers::WithinAbs(one[i], 1e-12));
  }
}

TEST_CASE("voronoiFeatures matches voronoiFeature for atom 0",
          "[structure_desc]") {
  auto cloud = fcc();
  auto all = chill::voronoiFeatures(cloud, 4.8);
  auto one = chill::voronoiFeature(cloud, 0, 4.8);
  REQUIRE(all.size() == static_cast<size_t>(cloud.nop));
  REQUIRE(all[0].size() == 3);
  REQUIRE(one.size() == 3);
  for (size_t i = 0; i < 3; i++) {
    REQUIRE_THAT(all[0][i], Catch::Matchers::WithinAbs(one[i], 1e-12));
  }
}
