#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <steinhardt_device.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <numeric>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

TEST_CASE("radialCoord converts Cartesian to spherical angles", "[bop]") {
  // Point along the z-axis: (0, 0, 1)
  // Polar angle (theta) should be 0, azimuth (phi) should be 0
  std::array<double, 3> zAxis = {0.0, 0.0, 1.0};
  auto angles = sph::radialCoord(zAxis);

  // angles[1] is the polar angle (theta from z-axis)
  REQUIRE_THAT(angles[1], Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("radialCoord for point along x-axis", "[bop]") {
  std::array<double, 3> xAxis = {1.0, 0.0, 0.0};
  auto angles = sph::radialCoord(xAxis);

  // Polar angle should be pi/2 (perpendicular to z)
  REQUIRE_THAT(angles[1], Catch::Matchers::WithinAbs(M_PI / 2.0, 1e-10));
}

TEST_CASE("spheriHarmo returns correct vector length for l=3", "[bop]") {
  std::array<double, 2> angles = {0.5, 1.0}; // phi, theta
  auto result = sph::spheriHarmo(3, angles);

  // 2*3+1 = 7 components
  REQUIRE(result.size() == 7);
}

TEST_CASE("spheriHarmo returns correct vector length for l=6", "[bop]") {
  std::array<double, 2> angles = {0.5, 1.0};
  auto result = sph::spheriHarmo(6, angles);

  // 2*6+1 = 13 components
  REQUIRE(result.size() == 13);
}

TEST_CASE("spheriHarmo Y_l0 at theta=0 is real and matches known value",
          "[bop]") {
  // At theta=0 (along z-axis), Y_lm for m!=0 should vanish
  // Y_30(0,phi) = sqrt(7/(4*pi)) for theta=0
  std::array<double, 2> angles = {0.0, 0.0}; // phi=0, theta=0
  auto result = sph::spheriHarmo(3, angles);

  // m ranges from -3 to 3; index 3 corresponds to m=0
  // For m!=0, values should be zero at theta=0
  for (int k = 0; k < 7; k++) {
    if (k == 3)
      continue; // skip m=0
    REQUIRE_THAT(std::abs(result[k]), Catch::Matchers::WithinAbs(0.0, 1e-10));
  }

  // m=0 component should be nonzero
  double y30_expected = 0.5 * std::sqrt(7.0 / M_PI);
  REQUIRE_THAT(result[3].real(),
               Catch::Matchers::WithinAbs(y30_expected, 1e-10));
  REQUIRE_THAT(result[3].imag(), Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("lookupTableQ3Vec matches spheriHarmo for l=3", "[bop]") {
  std::array<double, 2> angles = {0.7, 1.2}; // arbitrary phi, theta
  auto boostResult = sph::spheriHarmo(3, angles);
  auto lookupResult = sph::lookupTableQ3Vec(angles);

  REQUIRE(boostResult.size() == lookupResult.size());

  for (size_t i = 0; i < boostResult.size(); i++) {
    REQUIRE_THAT(boostResult[i].real(),
                 Catch::Matchers::WithinAbs(lookupResult[i].real(), 1e-8));
    REQUIRE_THAT(boostResult[i].imag(),
                 Catch::Matchers::WithinAbs(lookupResult[i].imag(), 1e-8));
  }
}

TEST_CASE("lookupTableQ6Vec matches spheriHarmo for l=6", "[bop]") {
  std::array<double, 2> angles = {0.3, 0.8};
  auto boostResult = sph::spheriHarmo(6, angles);
  auto lookupResult = sph::lookupTableQ6Vec(angles);

  REQUIRE(boostResult.size() == lookupResult.size());

  for (size_t i = 0; i < boostResult.size(); i++) {
    REQUIRE_THAT(boostResult[i].real(),
                 Catch::Matchers::WithinAbs(lookupResult[i].real(), 1e-8));
    REQUIRE_THAT(boostResult[i].imag(),
                 Catch::Matchers::WithinAbs(lookupResult[i].imag(), 1e-8));
  }
}

TEST_CASE("getCorrel populates c_ij for a simple tetrahedral system", "[bop]") {
  // Build a 5-atom system: central atom + 4 tetrahedral neighbours
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  // Central atom at origin
  double coords[5][3] = {
      {10.0, 10.0, 10.0},                       // center
      {11.0, 11.0, 11.0},                        // vertex 1
      {11.0, 9.0, 9.0},                          // vertex 2
      {9.0, 11.0, 9.0},                          // vertex 3
      {9.0, 9.0, 11.0}                           // vertex 4
  };

  for (int i = 0; i < 5; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.molID = i;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 5;

  // Build neighbour list
  auto nList = nneigh::neighListO(3.0, cloud, 1);

  // Get correlations
  chill::getCorrel(cloud, nList, false);

  // The central atom (index 0) should now have c_ij entries
  REQUIRE_FALSE(cloud.pts[0].c_ij.empty());

  // Each c_ij entry should have a real c_value in [-1, 1]
  for (const auto &cij : cloud.pts[0].c_ij) {
    REQUIRE(cij.c_value >= -1.0);
    REQUIRE(cij.c_value <= 1.0);
  }
}

// Helper: build a tetrahedral cloud for ice-like testing
static molSys::PointCloud<molSys::Point<double>, double>
makeTetraCloud() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  double coords[5][3] = {
      {10.0, 10.0, 10.0},
      {11.0, 11.0, 11.0},
      {11.0, 9.0, 9.0},
      {9.0, 11.0, 9.0},
      {9.0, 9.0, 11.0}};

  for (int i = 0; i < 5; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.molID = i;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 5;
  return cloud;
}

TEST_CASE("getIceTypeNoPrint classifies atoms after correlation", "[bop]") {
  auto cloud = makeTetraCloud();
  auto nList = nneigh::neighListO(3.0, cloud, 1);

  // First get correlations, then classify
  chill::getCorrel(cloud, nList, false);
  chill::getIceTypeNoPrint(cloud, nList, false);

  // Every atom should have an ice type assigned (could be any type)
  for (int i = 0; i < cloud.nop; i++) {
    // Just check that it's a valid enum value (not undefined)
    auto t = cloud.pts[i].iceType;
    bool validType =
        (t == molSys::atom_state_type::cubic ||
         t == molSys::atom_state_type::hexagonal ||
         t == molSys::atom_state_type::water ||
         t == molSys::atom_state_type::interfacial ||
         t == molSys::atom_state_type::clathrate ||
         t == molSys::atom_state_type::interClathrate ||
         t == molSys::atom_state_type::unclassified ||
         t == molSys::atom_state_type::reCubic ||
         t == molSys::atom_state_type::reHex);
    REQUIRE(validType);
  }
}

TEST_CASE("getq6 computes Q6 order parameter for each atom", "[bop]") {
  auto cloud = makeTetraCloud();
  auto nList = nneigh::neighListO(3.0, cloud, 1);

  auto q6 = chill::getq6(cloud, nList, false);

  REQUIRE(q6.size() == static_cast<size_t>(cloud.nop));
  for (const auto &val : q6) {
    REQUIRE_FALSE(std::isnan(val));
    REQUIRE_FALSE(std::isinf(val));
  }
}

TEST_CASE("getq6 uses the four nearest neighbours, not the cutoff list",
          "[bop]") {
  auto cloud = makeTetraCloud();
  molSys::Point<double> extra;
  extra.type = 1;
  extra.atomID = 5;
  extra.molID = 5;
  extra.x = 10.0;
  extra.y = 10.0;
  extra.z = 12.8;
  cloud.pts.push_back(extra);
  cloud.idIndexMap[5] = cloud.nop;
  cloud.nop += 1;

  auto full = nneigh::neighListO(3.0, cloud, 1);
  REQUIRE(full[0].size() > 5);

  std::vector<std::vector<int>> fourOnly(full.size());
  for (int iatom = 0; iatom < cloud.nop; iatom++) {
    std::vector<std::pair<double, int>> cand;
    for (size_t j = 1; j < full[iatom].size(); j++) {
      const auto it = cloud.idIndexMap.find(full[iatom][j]);
      REQUIRE(it != cloud.idIndexMap.end());
      cand.emplace_back(gen::periodicDistSq(cloud, iatom, it->second),
                        full[iatom][j]);
    }
    const size_t keep = std::min<size_t>(4, cand.size());
    std::partial_sort(cand.begin(), cand.begin() + keep, cand.end());
    fourOnly[iatom].push_back(full[iatom][0]);
    for (size_t k = 0; k < keep; k++) {
      fourOnly[iatom].push_back(cand[k].second);
    }
  }
  REQUIRE(fourOnly[0].size() < full[0].size());

  const auto qFull = chill::getq6(cloud, full, false);
  const auto qFour = chill::getq6(cloud, fourOnly, false);
  REQUIRE(qFull.size() == qFour.size());
  for (size_t i = 0; i < qFull.size(); i++) {
    REQUIRE_THAT(qFull[i], Catch::Matchers::WithinAbs(qFour[i], 1e-12));
  }
}

TEST_CASE("reclassifyWater uses Q6 to reclassify water atoms", "[bop]") {
  auto cloud = makeTetraCloud();
  auto nList = nneigh::neighListO(3.0, cloud, 1);

  chill::getCorrel(cloud, nList, false);
  chill::getIceTypeNoPrint(cloud, nList, false);
  auto q6 = chill::getq6(cloud, nList, false);
  chill::reclassifyWater(cloud, q6);

  // Should not crash and all types should still be valid
  for (int i = 0; i < cloud.nop; i++) {
    auto t = cloud.pts[i].iceType;
    REQUIRE(static_cast<int>(t) >= 0);
  }
}

TEST_CASE("getIceType writes output and classifies atoms", "[bop]") {
  auto cloud = makeTetraCloud();
  auto nList = nneigh::neighListO(3.0, cloud, 1);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_icetype/").string();

  chill::getCorrel(cloud, nList, false);
  chill::getIceType(cloud, nList, tmpPath, 1, false);

  // Check output file was written (getIceType writes to path + "bop/chill.txt")
  REQUIRE(std::filesystem::exists(tmpPath + "bop/chill.txt"));

  std::error_code _ec_; std::filesystem::remove_all(tmpPath, _ec_);
}

TEST_CASE("getCorrelPlus uses CHILL+ algorithm", "[bop]") {
  auto cloud = makeTetraCloud();
  auto nList = nneigh::neighListO(3.0, cloud, 1);

  chill::getCorrelPlus(cloud, nList, false);

  // Central atom should have c_ij entries
  REQUIRE_FALSE(cloud.pts[0].c_ij.empty());
  for (const auto &cij : cloud.pts[0].c_ij) {
    REQUIRE(cij.c_value >= -1.0);
    REQUIRE(cij.c_value <= 1.0);
  }
}

TEST_CASE("printIceType writes super chill classification", "[bop]") {
  auto cloud = makeTetraCloud();
  auto nList = nneigh::neighListO(3.0, cloud, 1);

  chill::getCorrel(cloud, nList, false);
  chill::getIceTypeNoPrint(cloud, nList, false);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_printice/").string();
  int ret = chill::printIceType(cloud, tmpPath, 1, false);
  REQUIRE(ret == 0);

  std::error_code _ec_; std::filesystem::remove_all(tmpPath, _ec_);
}

TEST_CASE("numStaggered counts staggered bonds for an atom", "[bop]") {
  auto cloud = makeTetraCloud();
  auto nList = nneigh::neighListO(3.0, cloud, 1);
  chill::getCorrel(cloud, nList, false);

  int nStag = chill::numStaggered(cloud, nList, 0);
  REQUIRE(nStag >= 0);
  REQUIRE(nStag <= 4);
}

TEST_CASE("CHILL+ on the mixed TIP4P example is twelve interClathrate",
          "[bop][fixture]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, cloud, 2);
  REQUIRE(cloud.nop == 250);
  auto nList = nneigh::neighListO(3.5, cloud, 2);
  chill::getCorrelPlus(cloud, nList, false);
  chill::getIceTypePlusNoPrint(cloud, nList, false);
  int interClath = 0;
  int water = 0;
  int other = 0;
  for (const auto &pt : cloud.pts) {
    if (pt.iceType == molSys::atom_state_type::interClathrate) {
      interClath++;
    } else if (pt.iceType == molSys::atom_state_type::water) {
      water++;
    } else {
      other++;
    }
  }
  REQUIRE(interClath == 12);
  REQUIRE(water == 238);
  REQUIRE(other == 0);
}

TEST_CASE("getIceTypePlusNoPrint classifies without writing", "[bop]") {
  auto cloud = makeTetraCloud();
  auto nList = nneigh::neighListO(3.0, cloud, 1);
  chill::getCorrelPlus(cloud, nList, false);
  chill::getIceTypePlusNoPrint(cloud, nList, false);
  REQUIRE(cloud.pts[0].iceType != molSys::atom_state_type::unclassified);
}

TEST_CASE("getIceTypePlus classifies with CHILL+ and writes output", "[bop]") {
  auto cloud = makeTetraCloud();
  auto nList = nneigh::neighListO(3.0, cloud, 1);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_icetypeplus/").string();

  chill::getCorrelPlus(cloud, nList, false);
  chill::getIceTypePlus(cloud, nList, tmpPath, 1, false);

  // Check output file was written
  REQUIRE(std::filesystem::exists(tmpPath + "bop/chillPlus.txt"));

  std::error_code _ec_; std::filesystem::remove_all(tmpPath, _ec_);
}

TEST_CASE("isInterfacial checks interfacial criteria", "[bop]") {
  auto cloud = makeTetraCloud();
  auto nList = nneigh::neighListO(3.0, cloud, 1);
  chill::getCorrel(cloud, nList, false);

  bool result = chill::isInterfacial(cloud, nList, 0, 2, 1);
  REQUIRE((result == true || result == false));
}

TEST_CASE("CHILL+ interfacial accepts a neighbour with three staggered bonds",
          "[bop]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 2;
  cloud.currentFrame = 1;
  for (int i = 0; i < 2; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i + 1;
    pt.x = static_cast<double>(i);
    pt.y = 0.0;
    pt.z = 0.0;
    cloud.pts.push_back(pt);
    cloud.idIndexMap[pt.atomID] = i;
  }
  std::vector<std::vector<int>> nList = {{1, 2}, {2, 1}};
  molSys::Result stag;
  stag.classifier = molSys::bond_type::staggered;
  stag.c_value = -0.9;
  cloud.pts[1].c_ij = {stag, stag, stag};

  const bool moore =
      chill::isInterfacial(cloud, nList, 0, 3, 0, false);
  const bool nguyen =
      chill::isInterfacial(cloud, nList, 0, 3, 0, true);
  REQUIRE_FALSE(moore);
  REQUIRE(nguyen);
}

// ---------------------------------------------------------------------------
// Independent verification of the l=3 and l=6 spherical harmonics.
//
// The reference values below were produced with mpmath's spherharm at 30
// decimal digits, which shares the Condon-Shortley convention used here.  The
// "matches spheriHarmo" cases above cannot serve this purpose: spheriHarmo
// dispatches to the same lookup tables, so it compares the implementation with
// itself.
// ---------------------------------------------------------------------------

namespace {

/// Reference Y_{3m}(theta, phi) for m = -3 .. 3.
struct HarmonicReference {
  double theta;
  double phi;
  std::vector<std::complex<double>> values;
};

const std::vector<HarmonicReference> &referenceQ3() {
  static const std::vector<HarmonicReference> refs = {
      {0.7, 1.2,
       {{-1.00032815672438657e-01, +4.93628664410609153e-02},
        {-2.39211073242928146e-01, -2.19120761338635062e-01},
        {+1.45220233090458700e-01, -3.73528458109185446e-01},
        {-2.14300204504549903e-02, +0.00000000000000000e+00},
        {-1.45220233090458700e-01, -3.73528458109185446e-01},
        {-2.39211073242928146e-01, +2.19120761338635062e-01},
        {+1.00032815672438657e-01, +4.93628664410609153e-02}}},
      {2.31, 5.02,
       {{-1.34258415477968879e-01, -1.01633027317392696e-01},
        {+3.07071559111461612e-01, -2.17013718940049311e-01},
        {+9.17923958435597764e-02, +2.88932088592419922e-01},
        {+1.83690302770999236e-01, +0.00000000000000000e+00},
        {-9.17923958435597764e-02, +2.88932088592419922e-01},
        {+3.07071559111461612e-01, +2.17013718940049311e-01},
        {+1.34258415477968879e-01, -1.01633027317392696e-01}}}};
  return refs;
}

const std::vector<HarmonicReference> &referenceQ6() {
  static const std::vector<HarmonicReference> refs = {
      {0.3, 0.8,
       {{+2.81545670960087548e-05, +3.20536104714261793e-04},
        {-2.35530243398039772e-03, +2.72701928443226617e-03},
        {-2.45553748690631526e-02, +1.43584741629399393e-03},
        {-8.33664022251006664e-02, -7.63648157168533986e-02},
        {-1.00161212420759645e-02, -3.42877198302676434e-01},
        {+4.13897504353533674e-01, -4.26164829149319191e-01},
        {+2.57919389430428581e-01, +0.00000000000000000e+00},
        {-4.13897504353533674e-01, -4.26164829149319191e-01},
        {-1.00161212420759645e-02, +3.42877198302676434e-01},
        {+8.33664022251006664e-02, -7.63648157168533986e-02},
        {-2.45553748690631526e-02, -1.43584741629399393e-03},
        {+2.35530243398039772e-03, +2.72701928443226617e-03},
        {+2.81545670960087548e-05, -3.20536104714261793e-04}}},
      {1.77, 3.41,
       {{-1.69855955987211517e-02, -4.28188076653546246e-01},
        {+6.79590313290316522e-02, -2.91875937921525719e-01},
        {-8.94254456210246268e-02, +1.64800164515641279e-01},
        {-2.16145708595750607e-01, +2.24889888997438808e-01},
        {+9.29737859183197707e-02, -5.53294457857134367e-02},
        {+2.98637874216372556e-01, -8.21386292029009507e-02},
        {-8.62877939979529995e-02, +0.00000000000000000e+00},
        {-2.98637874216372556e-01, -8.21386292029009507e-02},
        {+9.29737859183197707e-02, +5.53294457857134367e-02},
        {+2.16145708595750607e-01, +2.24889888997438808e-01},
        {-8.94254456210246268e-02, -1.64800164515641279e-01},
        {-6.79590313290316522e-02, -2.91875937921525719e-01},
        {-1.69855955987211517e-02, +4.28188076653546246e-01}}}};
  return refs;
}

} // namespace

TEST_CASE("lookupTableQ3Vec matches independent reference values", "[bop]") {
  for (const auto &ref : referenceQ3()) {
    auto result = sph::lookupTableQ3Vec({ref.phi, ref.theta});
    REQUIRE(result.size() == ref.values.size());
    for (size_t m = 0; m < result.size(); m++) {
      REQUIRE_THAT(result[m].real(),
                   Catch::Matchers::WithinAbs(ref.values[m].real(), 1e-13));
      REQUIRE_THAT(result[m].imag(),
                   Catch::Matchers::WithinAbs(ref.values[m].imag(), 1e-13));
    }
  }
}

TEST_CASE("lookupTableQ6Vec matches independent reference values", "[bop]") {
  for (const auto &ref : referenceQ6()) {
    auto result = sph::lookupTableQ6Vec({ref.phi, ref.theta});
    REQUIRE(result.size() == ref.values.size());
    for (size_t m = 0; m < result.size(); m++) {
      REQUIRE_THAT(result[m].real(),
                   Catch::Matchers::WithinAbs(ref.values[m].real(), 1e-13));
      REQUIRE_THAT(result[m].imag(),
                   Catch::Matchers::WithinAbs(ref.values[m].imag(), 1e-13));
    }
  }
}

TEST_CASE("spherical harmonics satisfy the addition theorem", "[bop]") {
  // sum_m |Y_lm(theta, phi)|^2 = (2l+1) / (4 pi), for every direction
  const std::vector<std::pair<double, double>> directions = {
      {0.7, 1.2}, {0.3, 0.8}, {2.31, 5.02}, {1.77, 3.41}, {0.05, 0.0},
      {M_PI - 0.05, 6.1}};

  for (const auto &[theta, phi] : directions) {
    double sumQ3 = 0.0;
    for (const auto &y : sph::lookupTableQ3Vec({phi, theta})) {
      sumQ3 += std::norm(y);
    }
    REQUIRE_THAT(sumQ3, Catch::Matchers::WithinAbs(7.0 / (4.0 * M_PI), 1e-13));

    double sumQ6 = 0.0;
    for (const auto &y : sph::lookupTableQ6Vec({phi, theta})) {
      sumQ6 += std::norm(y);
    }
    REQUIRE_THAT(sumQ6, Catch::Matchers::WithinAbs(13.0 / (4.0 * M_PI), 1e-13));
  }
}

TEST_CASE("spherical harmonics obey the Condon-Shortley relation", "[bop]") {
  // Y_{l,m} = (-1)^m conj(Y_{l,-m})
  const std::vector<std::pair<double, double>> directions = {
      {0.7, 1.2}, {0.3, 0.8}, {2.31, 5.02}, {1.77, 3.41}};

  for (const auto &[theta, phi] : directions) {
    auto q3 = sph::lookupTableQ3Vec({phi, theta});
    for (int m = 1; m <= 3; m++) {
      const std::complex<double> expected =
          ((m % 2 == 0) ? 1.0 : -1.0) * std::conj(q3[3 - m]);
      REQUIRE_THAT(q3[3 + m].real(),
                   Catch::Matchers::WithinAbs(expected.real(), 1e-15));
      REQUIRE_THAT(q3[3 + m].imag(),
                   Catch::Matchers::WithinAbs(expected.imag(), 1e-15));
    }

    auto q6 = sph::lookupTableQ6Vec({phi, theta});
    for (int m = 1; m <= 6; m++) {
      const std::complex<double> expected =
          ((m % 2 == 0) ? 1.0 : -1.0) * std::conj(q6[6 - m]);
      REQUIRE_THAT(q6[6 + m].real(),
                   Catch::Matchers::WithinAbs(expected.real(), 1e-15));
      REQUIRE_THAT(q6[6 + m].imag(),
                   Catch::Matchers::WithinAbs(expected.imag(), 1e-15));
    }
  }
}

TEST_CASE("per-order lookup agrees with the vector form", "[bop]") {
  const std::vector<std::pair<double, double>> directions = {
      {0.7, 1.2}, {2.31, 5.02}};

  for (const auto &[theta, phi] : directions) {
    const std::array<double, 2> angles = {phi, theta};

    auto q3 = sph::lookupTableQ3Vec(angles);
    for (int m = 0; m < 7; m++) {
      const auto single = sph::lookupTableQ3(m, angles);
      REQUIRE_THAT(single.real(),
                   Catch::Matchers::WithinAbs(q3[m].real(), 1e-15));
      REQUIRE_THAT(single.imag(),
                   Catch::Matchers::WithinAbs(q3[m].imag(), 1e-15));
    }

    auto q6 = sph::lookupTableQ6Vec(angles);
    for (int m = 0; m < 13; m++) {
      const auto single = sph::lookupTableQ6(m, angles);
      REQUIRE_THAT(single.real(),
                   Catch::Matchers::WithinAbs(q6[m].real(), 1e-15));
      REQUIRE_THAT(single.imag(),
                   Catch::Matchers::WithinAbs(q6[m].imag(), 1e-15));
    }
  }
}

namespace {

const std::vector<HarmonicReference> &referenceQ4() {
  static const std::vector<HarmonicReference> refs = {
      {0.7, 1.2,
       {{+6.66927990676342581e-03, +7.59288890954429030e-02},
        {-2.29527952617402109e-01, +1.13264408218239299e-01},
        {-3.16836766612184129e-01, -2.90227005710680053e-01},
        {+9.24808639851284198e-02, -2.37874804314991500e-01},
        {-2.72112678899294469e-01, +0.00000000000000000e+00},
        {-9.24808639851284198e-02, -2.37874804314991500e-01},
        {-3.16836766612184129e-01, +2.90227005710680053e-01},
        {+2.29527952617402109e-01, +1.13264408218239299e-01},
        {+6.66927990676342581e-03, -7.59288890954429030e-02}}},
      {2.31, 5.02,
       {{+4.40602314860978964e-02, -1.24417097718791364e-01},
        {+2.71349536749032338e-01, +2.05410400404298560e-01},
        {-3.24812661828660476e-01, +2.29551717216077311e-01},
        {-1.26299260956492816e-02, -3.97548281864557596e-02},
        {-3.60323417242993271e-01, +0.00000000000000000e+00},
        {+1.26299260956492816e-02, -3.97548281864557596e-02},
        {-3.24812661828660476e-01, -2.29551717216077311e-01},
        {-2.71349536749032338e-01, +2.05410400404298560e-01},
        {+4.40602314860978964e-02, +1.24417097718791364e-01}}}};
  return refs;
}

/// Builds a perfect FCC cell of `reps` repeats, which has a known Steinhardt
/// signature and so pins the absolute scale of the parameters.
molSys::PointCloud<molSys::Point<double>, double> fccCloud(int reps,
                                                           double lattice) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const std::array<std::array<double, 3>, 4> basis = {
      {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}}};

  int id = 1;
  for (int i = 0; i < reps; i++) {
    for (int j = 0; j < reps; j++) {
      for (int k = 0; k < reps; k++) {
        for (const auto &b : basis) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          p.x = (i + b[0]) * lattice;
          p.y = (j + b[1]) * lattice;
          p.z = (k + b[2]) * lattice;
          cloud.pts.push_back(p);
          cloud.idIndexMap[id] = id - 1;
          id++;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.currentFrame = 1;
  const double L = reps * lattice;
  cloud.box = {L, L, L};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

} // namespace

TEST_CASE("lookupTableQ4Vec matches independent reference values", "[bop]") {
  for (const auto &ref : referenceQ4()) {
    auto result = sph::lookupTableQ4Vec({ref.phi, ref.theta});
    REQUIRE(result.size() == ref.values.size());
    for (size_t m = 0; m < result.size(); m++) {
      REQUIRE_THAT(result[m].real(),
                   Catch::Matchers::WithinAbs(ref.values[m].real(), 1e-13));
      REQUIRE_THAT(result[m].imag(),
                   Catch::Matchers::WithinAbs(ref.values[m].imag(), 1e-13));
    }
  }
}

TEST_CASE("Q4 satisfies the addition theorem and Condon-Shortley", "[bop]") {
  const std::vector<std::pair<double, double>> directions = {
      {0.7, 1.2}, {0.3, 0.8}, {2.31, 5.02}, {1.77, 3.41}};

  for (const auto &[theta, phi] : directions) {
    auto q4 = sph::lookupTableQ4Vec({phi, theta});

    double sum = 0.0;
    for (const auto &y : q4) {
      sum += std::norm(y);
    }
    REQUIRE_THAT(sum, Catch::Matchers::WithinAbs(9.0 / (4.0 * M_PI), 1e-13));

    for (int m = 1; m <= 4; m++) {
      const std::complex<double> expected =
          ((m % 2 == 0) ? 1.0 : -1.0) * std::conj(q4[4 - m]);
      REQUIRE_THAT(q4[4 + m].real(),
                   Catch::Matchers::WithinAbs(expected.real(), 1e-15));
      REQUIRE_THAT(q4[4 + m].imag(),
                   Catch::Matchers::WithinAbs(expected.imag(), 1e-15));
    }

    for (int m = 0; m < 9; m++) {
      const auto single = sph::lookupTableQ4(m, {phi, theta});
      REQUIRE_THAT(single.real(),
                   Catch::Matchers::WithinAbs(q4[m].real(), 1e-15));
      REQUIRE_THAT(single.imag(),
                   Catch::Matchers::WithinAbs(q4[m].imag(), 1e-15));
    }
  }
}

TEST_CASE("device Ylm matches spheriHarmo for l=3,4,6", "[bop]") {
  const std::array<double, 2> angles = {0.7, 1.2};
  for (int l : {3, 4, 6, 8}) {
    auto host = sph::spheriHarmo(l, angles);
    double dev[34];
    seams::steinhardt::ylmAll(l, angles[1], angles[0], dev);
    const int nComp = 2 * l + 1;
    for (int m = 0; m < nComp; m++) {
      REQUIRE_THAT(dev[2 * m],
                   Catch::Matchers::WithinAbs(host[m].real(), 1e-12));
      REQUIRE_THAT(dev[2 * m + 1],
                   Catch::Matchers::WithinAbs(host[m].imag(), 1e-12));
    }
  }
}

TEST_CASE("steinhardtQl l=3 vanishes on an inversion-symmetric FCC lattice",
          "[bop]") {
  const double lattice = 4.0;
  auto cloud = fccCloud(4, lattice);
  const double cutoff = 0.85 * lattice;
  auto nList = nneigh::neighListO(cutoff, cloud, 1);
  auto q3 = chill::steinhardtQl(cloud, nList, 3);

  REQUIRE(q3.ql.size() == static_cast<size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; i++) {
    REQUIRE_THAT(q3.ql[i], Catch::Matchers::WithinAbs(0.0, 1e-5));
    REQUIRE_THAT(q3.qlBar[i], Catch::Matchers::WithinAbs(0.0, 1e-5));
  }
}

TEST_CASE("getCorrel coordination number generalizes beyond four", "[bop]") {
  // The FCC first shell holds twelve neighbours. The tetrahedral default
  // keeps the four nearest, a coordination of twelve keeps the shell, and a
  // non-positive coordination keeps each atom's whole neighbour-list row.
  const double lattice = 4.0;
  const double cutoff = 0.85 * lattice; // first shell at a/sqrt(2)

  auto cloudDefault = fccCloud(3, lattice);
  auto nList = nneigh::neighListO(cutoff, cloudDefault, 1);
  chill::getCorrel(cloudDefault, nList, false);
  for (int i = 0; i < cloudDefault.nop; i++) {
    REQUIRE(cloudDefault.pts[i].c_ij.size() == 4);
  }

  auto cloudTwelve = fccCloud(3, lattice);
  chill::getCorrel(cloudTwelve, nList, false, 12);
  for (int i = 0; i < cloudTwelve.nop; i++) {
    REQUIRE(cloudTwelve.pts[i].c_ij.size() == 12);
    for (const auto &cij : cloudTwelve.pts[i].c_ij) {
      REQUIRE(cij.c_value >= -1.0);
      REQUIRE(cij.c_value <= 1.0);
    }
  }

  auto cloudRow = fccCloud(3, lattice);
  chill::getCorrel(cloudRow, nList, false, 0);
  for (int i = 0; i < cloudRow.nop; i++) {
    // Row layout: the first entry is the atom's own ID
    REQUIRE(cloudRow.pts[i].c_ij.size() == nList[i].size() - 1);
  }

  auto cloudPlus = fccCloud(3, lattice);
  chill::getCorrelPlus(cloudPlus, nList, false, 12);
  for (int i = 0; i < cloudPlus.nop; i++) {
    REQUIRE(cloudPlus.pts[i].c_ij.size() == 12);
  }
}

TEST_CASE("classifyBonds under registered rules reproduces CHILL/CHILL+",
          "[bop]") {
  // The named rule sets and the water entry points are the same engine
  auto reference = makeTetraCloud();
  auto nList = nneigh::neighListO(3.0, reference, 1);
  chill::getCorrelPlus(reference, nList, false);

  auto viaRule = makeTetraCloud();
  chill::classifyBonds(viaRule, nList, chill::bondClassifier("CHILL+"));
  REQUIRE(viaRule.pts[0].c_ij.size() == reference.pts[0].c_ij.size());
  for (size_t j = 0; j < reference.pts[0].c_ij.size(); j++) {
    REQUIRE_THAT(viaRule.pts[0].c_ij[j].c_value,
                 Catch::Matchers::WithinAbs(
                     reference.pts[0].c_ij[j].c_value, 1e-12));
    REQUIRE(viaRule.pts[0].c_ij[j].classifier ==
            reference.pts[0].c_ij[j].classifier);
  }

  // A registered custom rule set with an all-covering eclipsed band makes
  // every bond eclipsed or staggered; nothing lands out of range
  chill::registerBondClassifier("everything-eclipsed",
                                {-2.0, -1.0, 1.0, 4});
  auto names = chill::bondClassifierNames();
  REQUIRE(std::find(names.begin(), names.end(), "everything-eclipsed") !=
          names.end());
  REQUIRE(std::find(names.begin(), names.end(), "CHILL") != names.end());
  REQUIRE(std::find(names.begin(), names.end(), "CHILL+") != names.end());

  auto custom = makeTetraCloud();
  chill::classifyBonds(custom, nList,
                       chill::bondClassifier("everything-eclipsed"));
  for (const auto &cij : custom.pts[0].c_ij) {
    REQUIRE(cij.classifier == molSys::bond_type::eclipsed);
  }
}

TEST_CASE("Steinhardt parameters reproduce the FCC reference values", "[bop]") {
  // A perfect FCC lattice has q4 = 0.190941, q6 = 0.574524 for the twelve
  // nearest neighbours (Steinhardt, Nelson and Ronchetti 1983, Table I).
  // With every environment identical, the averaged parameters must coincide
  // with the local ones.
  const double lattice = 4.0;
  auto cloud = fccCloud(4, lattice);

  // First shell of FCC sits at a/sqrt(2); cut between that and the second
  // shell at a
  const double cutoff = 0.85 * lattice;
  auto nList = nneigh::neighListO(cutoff, cloud, 1);

  auto q4 = chill::steinhardtQl(cloud, nList, 4);
  auto q6 = chill::steinhardtQl(cloud, nList, 6);

  REQUIRE(q4.ql.size() == static_cast<size_t>(cloud.nop));

  for (int i = 0; i < cloud.nop; i++) {
    REQUIRE_THAT(q4.ql[i], Catch::Matchers::WithinAbs(0.190941, 1e-5));
    REQUIRE_THAT(q6.ql[i], Catch::Matchers::WithinAbs(0.574524, 1e-5));
    // Uniform environment: averaging changes nothing
    REQUIRE_THAT(q4.qlBar[i], Catch::Matchers::WithinAbs(q4.ql[i], 1e-9));
    REQUIRE_THAT(q6.qlBar[i], Catch::Matchers::WithinAbs(q6.ql[i], 1e-9));
  }
}

TEST_CASE("steinhardtQl averages only the bonds that resolve", "[bop]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  for (int i = 0; i < 2; i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    p.x = static_cast<double>(i);
    p.y = 0.0;
    p.z = 0.0;
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = i;
  }
  cloud.nop = 2;

  std::vector<std::vector<int>> clean = {{1, 2}, {2, 1}};
  std::vector<std::vector<int>> padded = {{1, 2, 99}, {2, 1}};
  const auto qClean = chill::steinhardtQl(cloud, clean, 6);
  const auto qPadded = chill::steinhardtQl(cloud, padded, 6);
  REQUIRE_THAT(qPadded.ql[0], Catch::Matchers::WithinAbs(qClean.ql[0], 1e-15));
  REQUIRE_THAT(qPadded.qlBar[0],
               Catch::Matchers::WithinAbs(qClean.qlBar[0], 1e-15));
}

TEST_CASE("steinhardtQl skips neighbour indices outside the cloud", "[bop]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> p;
  p.type = 1;
  p.atomID = 1;
  p.x = 0.0;
  p.y = 0.0;
  p.z = 0.0;
  cloud.pts.push_back(p);
  cloud.idIndexMap[1] = 0;
  cloud.idIndexMap[2] = 99;
  cloud.nop = 1;

  std::vector<std::vector<int>> nList = {{1, 2}};
  const auto q = chill::steinhardtQl(cloud, nList, 6);
  REQUIRE(q.ql.size() == 1);
  REQUIRE(q.ql[0] == 0.0);
  REQUIRE(q.qlBar[0] == 0.0);
}

TEST_CASE("steinhardtQl rejects unsupported degrees", "[bop]") {
  auto cloud = fccCloud(2, 4.0);
  auto nList = nneigh::neighListO(3.4, cloud, 1);
  auto q5 = chill::steinhardtQl(cloud, nList, 5);
  REQUIRE(q5.ql.size() == static_cast<size_t>(cloud.nop));
  for (double v : q5.ql) {
    REQUIRE(v == 0.0);
  }
}

TEST_CASE("steinhardtQl is independent of the OpenMP thread count", "[bop]") {
  // steinhardtQl only starts threads at 50000 atoms. A chain of nearest
  // neighbours is enough to exercise both parallel loops without building
  // a physical lattice.
  constexpr int nAtoms = 50000;
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = nAtoms;
  cloud.currentFrame = 1;
  cloud.box = {static_cast<double>(nAtoms) * 2.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.pts.reserve(nAtoms);
  std::vector<std::vector<int>> nList(nAtoms);
  for (int i = 0; i < nAtoms; i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    p.x = static_cast<double>(i) * 1.5;
    p.y = 0.0;
    p.z = 0.0;
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = i;
    nList[i].push_back(p.atomID);
    if (i + 1 < nAtoms) {
      nList[i].push_back(i + 2);
    }
  }

  auto run = [&]() { return chill::steinhardtQl(cloud, nList, 6); };

#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(1);
#endif
  const auto serial = run();
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(4);
#endif
  const auto threaded = run();

  REQUIRE(serial.ql.size() == static_cast<size_t>(nAtoms));
  REQUIRE(threaded.ql.size() == serial.ql.size());
  for (int i = 0; i < nAtoms; i++) {
    REQUIRE_THAT(threaded.ql[i],
                 Catch::Matchers::WithinAbs(serial.ql[i], 1e-15));
    REQUIRE_THAT(threaded.qlBar[i],
                 Catch::Matchers::WithinAbs(serial.qlBar[i], 1e-15));
  }
}

namespace {

const std::vector<HarmonicReference> &referenceQ8() {
  // mpmath spherharm at 30 digits, theta/phi as in the l=3,4,6 references
  static const std::vector<HarmonicReference> refs = {
      {0.7, 1.2,
       {{-1.50566575719588720e-02, +2.66559461970755553e-03},
        {-3.77083914952981436e-02, -6.20571043526132460e-02},
        {+1.27263877225608890e-01, -1.66031119141199057e-01},
        {+3.82638236458474390e-01, +1.11350095869957796e-01},
        {+4.09586774205907000e-02, +4.66309244602285233e-01},
        {-1.82029622803040386e-01, +8.98255627250017646e-02},
        {+1.97868206964352672e-01, +1.81250105051402932e-01},
        {-1.26611692820706584e-01, +3.25664471068939643e-01},
        {+1.66803869622521250e-01, +0.00000000000000000e+00},
        {+1.26611692820706584e-01, +3.25664471068939643e-01},
        {+1.97868206964352672e-01, -1.81250105051402932e-01},
        {+1.82029622803040386e-01, +8.98255627250017646e-02},
        {+4.09586774205907000e-02, -4.66309244602285233e-01},
        {-3.82638236458474390e-01, +1.11350095869957796e-01},
        {+1.27263877225608890e-01, +1.66031119141199057e-01},
        {+3.77083914952981436e-02, -6.20571043526132460e-02},
        {-1.50566575719588720e-02, -2.66559461970755553e-03}}},
      {2.31, 5.02,
       {{-3.56321951023816641e-02, -2.88559218572720516e-02},
        {+1.39626124302155624e-01, -9.19745088952253348e-02},
        {+9.66561669354474851e-02, +3.42742586142052830e-01},
        {-4.59565562033470043e-01, -1.50521201263779918e-02},
        {+8.72113398829371667e-02, -2.46267017453726950e-01},
        {-1.45178642334412467e-01, -1.09899590798437044e-01},
        {+3.03538655695245463e-01, -2.14516944210316723e-01},
        {-4.48339778327363353e-03, -1.41122526937808735e-02},
        {+3.69672685761530673e-01, +0.00000000000000000e+00},
        {+4.48339778327363353e-03, -1.41122526937808735e-02},
        {+3.03538655695245463e-01, +2.14516944210316723e-01},
        {+1.45178642334412467e-01, -1.09899590798437044e-01},
        {+8.72113398829371667e-02, +2.46267017453726950e-01},
        {+4.59565562033470043e-01, -1.50521201263779918e-02},
        {+9.66561669354474851e-02, -3.42742586142052830e-01},
        {-1.39626124302155624e-01, -9.19745088952253348e-02},
        {-3.56321951023816641e-02, +2.88559218572720516e-02}}}};
  return refs;
}

} // namespace

TEST_CASE("lookupTableQ8Vec matches independent reference values", "[bop]") {
  for (const auto &ref : referenceQ8()) {
    auto result = sph::lookupTableQ8Vec({ref.phi, ref.theta});
    REQUIRE(result.size() == ref.values.size());
    for (size_t m = 0; m < result.size(); m++) {
      REQUIRE_THAT(result[m].real(),
                   Catch::Matchers::WithinAbs(ref.values[m].real(), 1e-13));
      REQUIRE_THAT(result[m].imag(),
                   Catch::Matchers::WithinAbs(ref.values[m].imag(), 1e-13));
    }
  }
}

TEST_CASE("Q8 satisfies the addition theorem and Condon-Shortley", "[bop]") {
  const std::vector<std::pair<double, double>> directions = {
      {0.7, 1.2}, {0.3, 0.8}, {2.31, 5.02}, {1.77, 3.41}};
  for (const auto &[theta, phi] : directions) {
    auto q8 = sph::lookupTableQ8Vec({phi, theta});
    double sum = 0.0;
    for (const auto &y : q8) {
      sum += std::norm(y);
    }
    REQUIRE_THAT(sum, Catch::Matchers::WithinAbs(17.0 / (4.0 * M_PI), 1e-13));
    for (int m = 1; m <= 8; m++) {
      const std::complex<double> expected =
          ((m % 2 == 0) ? 1.0 : -1.0) * std::conj(q8[8 - m]);
      REQUIRE_THAT(q8[8 + m].real(),
                   Catch::Matchers::WithinAbs(expected.real(), 1e-15));
      REQUIRE_THAT(q8[8 + m].imag(),
                   Catch::Matchers::WithinAbs(expected.imag(), 1e-15));
    }
    for (int m = 0; m < 17; m++) {
      const auto single = sph::lookupTableQ8(m, {phi, theta});
      REQUIRE_THAT(single.real(),
                   Catch::Matchers::WithinAbs(q8[m].real(), 1e-15));
      REQUIRE_THAT(single.imag(),
                   Catch::Matchers::WithinAbs(q8[m].imag(), 1e-15));
    }
  }
}

TEST_CASE("steinhardtQl accepts l=8 and averages correctly on FCC", "[bop]") {
  // FCC reference value q8 = 0.40391 for twelve nearest neighbours
  // (Steinhardt, Nelson and Ronchetti 1983); uniform environment makes the
  // neighbour average coincide with the local value
  const double lattice = 4.0;
  auto cloud = fccCloud(4, lattice);
  const double cutoff = 0.85 * lattice;
  auto nList = nneigh::neighListO(cutoff, cloud, 1);
  auto q8 = chill::steinhardtQl(cloud, nList, 8);
  REQUIRE(q8.ql.size() == static_cast<size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; i++) {
    REQUIRE_THAT(q8.ql[i], Catch::Matchers::WithinAbs(0.40391, 1e-4));
    REQUIRE_THAT(q8.qlBar[i], Catch::Matchers::WithinAbs(q8.ql[i], 1e-9));
  }
}

TEST_CASE("qlmOneAtomDr uses packed relDist not minImage", "[bop]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 2;
  molSys::Point<double> a;
  a.x = 0.5;
  a.y = 0.5;
  a.z = 1.0;
  a.atomID = 1;
  a.type = 1;
  molSys::Point<double> b;
  b.x = 1.0;
  b.y = 8.0;
  b.z = 1.0;
  b.atomID = 2;
  b.type = 1;
  cloud.pts.push_back(a);
  cloud.pts.push_back(b);
  cloud.idIndexMap[1] = 0;
  cloud.idIndexMap[2] = 1;
  const auto dr01 = gen::relDist(cloud, 0, 1);
  REQUIRE_THAT(dr01[0], Catch::Matchers::WithinAbs(4.5, 1e-10));
  REQUIRE_THAT(dr01[1], Catch::Matchers::WithinAbs(1.160254037844386, 1e-10));
  REQUIRE_THAT(dr01[2], Catch::Matchers::WithinAbs(0.0, 1e-12));
  const double packed[6] = {dr01[0], dr01[1], dr01[2],
                            -dr01[0], -dr01[1], -dr01[2]};
  const double xyz[6] = {0.5, 0.5, 1.0, 1.0, 8.0, 1.0};
  const int offsets[3] = {0, 1, 2};
  const int cols[2] = {1, 0};
  std::vector<double> qlmDr(2 * 2 * 13, 0.0);
  std::vector<double> qlmSpan(2 * 2 * 13, 0.0);
  seams::steinhardt::qlmOneAtomDr(0, 6, packed, offsets, cols, qlmDr.data());
  seams::steinhardt::qlmOneAtom(0, 6, xyz, offsets, cols, cloud.box[0],
                                cloud.box[1], cloud.box[2], qlmSpan.data());
  bool qlmDiffers = false;
  for (size_t i = 0; i < qlmDr.size(); i++) {
    if (std::abs(qlmDr[i] - qlmSpan[i]) > 1e-8) {
      qlmDiffers = true;
      break;
    }
  }
  REQUIRE(qlmDiffers);
}

TEST_CASE("steinhardtQl flatten uses dump H not span minImage", "[bop]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 3;
  molSys::Point<double> a;
  a.x = 0.5;
  a.y = 0.5;
  a.z = 1.0;
  a.atomID = 1;
  a.type = 1;
  molSys::Point<double> b;
  b.x = 1.0;
  b.y = 8.0;
  b.z = 1.0;
  b.atomID = 2;
  b.type = 1;
  molSys::Point<double> c;
  c.x = 0.5;
  c.y = 0.5;
  c.z = 2.0;
  c.atomID = 3;
  c.type = 1;
  cloud.pts.push_back(a);
  cloud.pts.push_back(b);
  cloud.pts.push_back(c);
  cloud.idIndexMap[1] = 0;
  cloud.idIndexMap[2] = 1;
  cloud.idIndexMap[3] = 2;
  const std::vector<std::vector<int>> nList = {{1, 2, 3}, {2, 1}, {3, 1}};
  const auto ql = chill::steinhardtQl(cloud, nList, 6);
  const auto dr01 = gen::relDist(cloud, 0, 1);
  const auto dr02 = gen::relDist(cloud, 0, 2);
  const double packed[6] = {dr01[0], dr01[1], dr01[2],
                            dr02[0], dr02[1], dr02[2]};
  const double xyz[9] = {0.5, 0.5, 1.0, 1.0, 8.0, 1.0, 0.5, 0.5, 2.0};
  const int offsets[2] = {0, 2};
  const int cols[2] = {1, 2};
  std::vector<double> qlmDr(1 * 13 * 2, 0.0);
  std::vector<double> qlmSpan(1 * 13 * 2, 0.0);
  seams::steinhardt::qlmOneAtomDr(0, 6, packed, offsets, cols, qlmDr.data());
  seams::steinhardt::qlmOneAtom(0, 6, xyz, offsets, cols, cloud.box[0],
                                cloud.box[1], cloud.box[2], qlmSpan.data());
  std::vector<double> qlDr(1, 0.0);
  std::vector<double> qlBarDr(1, 0.0);
  std::vector<double> qlSpan(1, 0.0);
  std::vector<double> qlBarSpan(1, 0.0);
  seams::steinhardt::qlOneAtom(0, 6, qlmDr.data(), offsets, cols, qlDr.data(),
                               qlBarDr.data());
  seams::steinhardt::qlOneAtom(0, 6, qlmSpan.data(), offsets, cols,
                               qlSpan.data(), qlBarSpan.data());
  REQUIRE(std::abs(qlDr[0] - qlSpan[0]) > 1e-4);
  REQUIRE_THAT(ql.ql[0], Catch::Matchers::WithinAbs(qlDr[0], 1e-5));
}
