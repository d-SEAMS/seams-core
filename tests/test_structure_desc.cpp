#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
#include <generic.hpp>
#include <ira_sofi.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <structure_desc.hpp>
#include <voronoi_qlm.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <numbers>
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

namespace {

// n atoms at water density, uniform from a portable xorshift, in an
// orthorhombic or tilted box; atom i has ID 3 i + 7
Cloud xorshiftFrame(int n, double tilt) {
  Cloud cloud;
  const double L = std::cbrt(n / 0.0332);
  const double xy = tilt;
  const double xz = -0.5 * tilt;
  const double yz = 0.3 * tilt;
  if (tilt == 0.0) {
    cloud.box = {L, L, L};
    cloud.boxLow = {0.0, 0.0, 0.0};
  } else {
    const double xmin = std::min({0.0, xy, xz, xy + xz});
    const double xmax = std::max({0.0, xy, xz, xy + xz});
    cloud.box = {L + xmax - xmin, L + std::max(0.0, yz) - std::min(0.0, yz), L,
                 xy, xz, yz};
    cloud.boxLow = {xmin, std::min(0.0, yz), 0.0};
  }
  unsigned long long state = 0x2545F4914F6CDD1DULL;
  auto unit = [&state]() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return static_cast<double>(state >> 11) / 9007199254740992.0;
  };
  for (int i = 0; i < n; i++) {
    const double sx = unit();
    const double sy = unit();
    const double sz = unit();
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = 3 * i + 7;
    p.molID = p.atomID;
    p.x = L * sx + xy * sy + xz * sz;
    p.y = L * sy + yz * sz;
    p.z = L * sz;
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = i;
  }
  cloud.nop = n;
  cloud.currentFrame = 1;
  return cloud;
}

// The power spectrum term by term: every Y_lm evaluated afresh for each
// neighbour, radial function and degree
std::vector<double> soapByTerms(const Cloud &cloud,
                                const std::vector<std::vector<int>> &nList,
                                int iatom, int nMax, int lMax, double rcut) {
  const int nComp = (lMax + 1) * (lMax + 1);
  const double sigma = rcut / static_cast<double>(nMax);
  std::vector<std::complex<double>> coeff(
      static_cast<size_t>(nMax) * static_cast<size_t>(nComp), {0.0, 0.0});
  for (size_t j = 1; j < nList[static_cast<size_t>(iatom)].size(); j++) {
    const int jatom = cloud.idIndexMap.at(nList[static_cast<size_t>(iatom)][j]);
    const auto d = gen::relDist(cloud, iatom, jatom);
    const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (r <= 0.0 || r >= rcut) {
      continue;
    }
    const std::array<double, 2> angles = {std::atan2(d[0], d[1]),
                                          std::acos(d[2] / r)};
    for (int n = 0; n < nMax; n++) {
      const double rn = (n + 0.5) * rcut / static_cast<double>(nMax);
      const double g = std::exp(-((r - rn) / sigma) * ((r - rn) / sigma));
      for (int l = 0; l <= lMax; l++) {
        if (l == 0) {
          coeff[static_cast<size_t>(n) * nComp] +=
              g * (0.5 / std::sqrt(std::numbers::pi));
          continue;
        }
        const auto yl = sph::spheriHarmo(l, angles);
        for (int m = 0; m < 2 * l + 1; m++) {
          coeff[static_cast<size_t>(n) * nComp + l * l + m] +=
              g * yl[static_cast<size_t>(m)];
        }
      }
    }
  }
  std::vector<double> spec;
  for (int n = 0; n < nMax; n++) {
    for (int np = 0; np < nMax; np++) {
      for (int l = 0; l <= lMax; l++) {
        std::complex<double> acc = 0.0;
        for (int m = 0; m < 2 * l + 1; m++) {
          acc += coeff[static_cast<size_t>(n) * nComp + l * l + m] *
                 std::conj(coeff[static_cast<size_t>(np) * nComp + l * l + m]);
        }
        spec.push_back(acc.real());
      }
    }
  }
  return spec;
}

} // namespace

TEST_CASE("SOAP equals its term-by-term expansion on disordered frames",
          "[structure_desc]") {
  for (const double tilt : {0.0, 3.0}) {
    INFO("tilt " << tilt);
    const auto cloud = xorshiftFrame(400, tilt);
    const auto nList = nneigh::neighListO(5.0, cloud, 1);
    const auto all = chill::soapSpectrumAll(cloud, nList, 4, 6, 5.0);
    REQUIRE(all.size() == 400);
    for (int i = 0; i < cloud.nop; i++) {
      REQUIRE(all[static_cast<size_t>(i)] ==
              soapByTerms(cloud, nList, i, 4, 6, 5.0));
    }
    REQUIRE(chill::soapSpectrum(cloud, 7, nList, 4, 6, 5.0) == all[7]);
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

TEST_CASE("spheriHarmo l=6 matches the Q6 lookup table", "[structure_desc]") {
  const std::array<double, 2> angles = {0.4, 1.1};
  const auto a = sph::spheriHarmo(6, angles);
  const auto b = sph::lookupTableQ6Vec(angles);
  REQUIRE(a.size() == 13);
  REQUIRE(b.size() == 13);
  for (size_t i = 0; i < a.size(); i++) {
    REQUIRE_THAT(a[i].real(), Catch::Matchers::WithinAbs(b[i].real(), 1e-12));
    REQUIRE_THAT(a[i].imag(), Catch::Matchers::WithinAbs(b[i].imag(), 1e-12));
  }
}

TEST_CASE("soapSpectrum includes l=1 and l=2", "[structure_desc]") {
  auto cloud = fcc();
  auto nList = nneigh::neighListO(3.2, cloud, 1);
  auto spec = chill::soapSpectrum(cloud, 0, nList, 2, 2, 3.2);
  REQUIRE(spec.size() == 2 * 2 * 3);
  double l1 = 0.0;
  double l2 = 0.0;
  // slots: for n,np in 0..1, l in 0..2
  int slot = 0;
  for (int n = 0; n < 2; n++) {
    for (int np = 0; np < 2; np++) {
      for (int l = 0; l <= 2; l++) {
        if (l == 1) {
          l1 += spec[static_cast<size_t>(slot)] * spec[static_cast<size_t>(slot)];
        }
        if (l == 2) {
          l2 += spec[static_cast<size_t>(slot)] * spec[static_cast<size_t>(slot)];
        }
        slot++;
      }
    }
  }
  REQUIRE(l1 > 0.0);
  REQUIRE(l2 > 0.0);
}
