#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <cmath>
#include <complex>
#include <filesystem>
#include <numeric>

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
  // Q6 values should be finite
  for (const auto &val : q6) {
    REQUIRE_FALSE(std::isnan(val));
    REQUIRE_FALSE(std::isinf(val));
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
  // Just check it doesn't crash; result depends on geometry
  REQUIRE((result == true || result == false));
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
