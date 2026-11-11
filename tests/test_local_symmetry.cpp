#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <local_symmetry.hpp>
#include <seams_c_api.h>

#include <Eigen/Geometry>

#include <cmath>
#include <string>
#include <vector>

namespace {

Eigen::MatrixXd rotatePermute(const Eigen::MatrixXd &src, double degrees, bool swap) {
  const double rad = degrees * 3.14159265358979323846 / 180.0;
  Eigen::Matrix3d R;
  R << std::cos(rad), -std::sin(rad), 0, std::sin(rad), std::cos(rad), 0, 0, 0, 1;
  Eigen::MatrixXd out(src.rows(), 3);
  for (Eigen::Index i = 0; i < src.rows(); i++) {
    const Eigen::Index j = swap && i < 2 ? 1 - i : i;
    out.row(j) = (R * src.row(i).transpose()).transpose();
    out(j, 0) += 0.02 * std::sin(static_cast<double>(i) + 1.0);
    out(j, 1) += 0.02 * std::cos(static_cast<double>(i) + 0.3);
    out(j, 2) += 0.015 * std::sin(0.5 * static_cast<double>(i));
  }
  return out;
}

double reconstructed(const Eigen::MatrixXd &ref, const Eigen::MatrixXd &env,
                     const symm::Overlay &hit) {
  double sum = 0.0;
  int used = 0;
  for (int i = 0; i < static_cast<int>(hit.permutation.size()); i++) {
    const int j = hit.permutation[static_cast<size_t>(i)];
    if (j < 0 || j >= env.rows()) {
      continue;
    }
    const Eigen::Vector3d pred =
        hit.rotation * env.row(j).transpose() + hit.translation;
    sum += (pred - ref.row(i).transpose()).squaredNorm();
    used++;
  }
  return used > 0 ? std::sqrt(sum / static_cast<double>(used)) : 1e9;
}

int neighboursNear(const Eigen::MatrixXd &m, int i, double lo, double hi) {
  int n = 0;
  for (Eigen::Index j = 0; j < m.rows(); j++) {
    if (j == i) {
      continue;
    }
    const double d = (m.row(i) - m.row(j)).norm();
    if (d >= lo && d <= hi) {
      n++;
    }
  }
  return n;
}

std::vector<int> chillLabels(const Eigen::MatrixXd &m) {
  const int n = static_cast<int>(m.rows());
  std::vector<double> xyz(static_cast<size_t>(3 * n));
  for (int i = 0; i < n; i++) {
    xyz[static_cast<size_t>(3 * i)] = m(i, 0) + 20.0;
    xyz[static_cast<size_t>(3 * i + 1)] = m(i, 1) + 20.0;
    xyz[static_cast<size_t>(3 * i + 2)] = m(i, 2) + 20.0;
  }
  const double box[6] = {0, 0, 0, 50, 50, 50};
  std::vector<int> labels(static_cast<size_t>(n), -1);
  REQUIRE(seams_chill_plus(xyz.data(), n, box, labels.data()) == 0);
  return labels;
}

} // namespace

TEST_CASE("ice shells and cages have the coordination the templates claim",
          "[symmetry]") {
  for (symm::Kind kind : {symm::Kind::iceIh, symm::Kind::iceIc, symm::Kind::iceIsd}) {
    const auto m = symm::points(kind);
    REQUIRE(m.rows() == 17);
    REQUIRE(neighboursNear(m, 0, 2.7, 2.8) == 4);
    REQUIRE(neighboursNear(m, 0, 4.4, 4.6) == 12);
  }
  const auto dod = symm::points(symm::Kind::cage512);
  const auto sI = symm::points(symm::Kind::cage51262);
  const auto sII = symm::points(symm::Kind::cage51264);
  REQUIRE(dod.rows() == 20);
  REQUIRE(sI.rows() == 24);
  REQUIRE(sII.rows() == 28);
  for (int i = 0; i < dod.rows(); i++) {
    REQUIRE(neighboursNear(dod, i, 2.7, 2.8) == 3);
  }
  for (int i = 0; i < sI.rows(); i++) {
    REQUIRE(neighboursNear(sI, i, 2.7, 2.85) == 3);
  }
  for (int i = 0; i < sII.rows(); i++) {
    REQUIRE(neighboursNear(sII, i, 2.6, 2.9) == 3);
  }
  REQUIRE(symm::points(symm::Kind::ring5).rows() == 5);
  REQUIRE(symm::points(symm::Kind::ring6).rows() == 6);
}

TEST_CASE("perturbed ice shells and cages match their own template", "[symmetry]") {
  const symm::Kind kinds[] = {symm::Kind::iceIh,   symm::Kind::iceIc,
                              symm::Kind::iceIsd,  symm::Kind::cage512,
                              symm::Kind::cage51262, symm::Kind::cage51264,
                              symm::Kind::ring5,   symm::Kind::ring6};
  for (symm::Kind kind : kinds) {
    const auto ref = symm::points(kind);
    const auto env = rotatePermute(ref, 0.0, true);
    const symm::Overlay hit = symm::best(env);
    INFO(symm::name(kind));
    REQUIRE(hit.status == 0);
    REQUIRE(hit.name == symm::name(kind));
    REQUIRE(hit.rmsd < 0.08);
    REQUIRE(reconstructed(ref, env, hit) < 0.08);
    REQUIRE(hit.permutation.size() == static_cast<size_t>(ref.rows()));
  }
}

TEST_CASE("IRA recovers a rotated cage and reports its point group", "[symmetry]") {
  if (!symm::iraLinked()) {
    SKIP("libira is not linked");
  }
  REQUIRE(std::string(symm::matcher()) == "ira");
  const auto ref = symm::points(symm::Kind::cage512);
  const auto env = rotatePermute(ref, 37.0, true);
  const symm::Overlay hit = symm::best(env, "ira");
  REQUIRE(hit.status == 0);
  REQUIRE(hit.name == "5^12");
  REQUIRE(hit.backend == "ira");
  REQUIRE(hit.rmsd < 0.08);
  REQUIRE(reconstructed(ref, env, hit) < 0.08);

  const symm::Group dod = symm::pointGroup(ref, 0.08);
  REQUIRE(dod.status == 0);
  REQUIRE(dod.symbol == "Ih");
  REQUIRE(dod.nOperations == 120);
  const symm::Group hex = symm::pointGroup(symm::points(symm::Kind::ring6), 0.05);
  REQUIRE(hex.symbol == "D6h");
  REQUIRE(hex.nOperations == 24);
  const symm::Group pen = symm::pointGroup(symm::points(symm::Kind::ring5), 0.05);
  REQUIRE(pen.symbol == "D5h");
  const symm::Group large = symm::pointGroup(symm::points(symm::Kind::cage51264), 0.08);
  REQUIRE(large.status == 0);
  REQUIRE(large.symbol == "Td");
  REQUIRE(large.nOperations == 24);
  const symm::Group sI = symm::pointGroup(symm::points(symm::Kind::cage51262), 0.08);
  REQUIRE(sI.symbol == "D6d");
  REQUIRE(sI.nOperations == 24);
  REQUIRE(symm::pointGroup(symm::points(symm::Kind::iceIc), 0.08).symbol == "Td");
  REQUIRE(symm::pointGroup(symm::points(symm::Kind::iceIh), 0.08).symbol == "C3v");
  REQUIRE(symm::pointGroup(symm::points(symm::Kind::iceIsd), 0.08).symbol == "C2v");
}

TEST_CASE("space groups separate ice Ih, ice Ic and a stacking fault", "[symmetry]") {
  if (!symm::spglibLinked()) {
    SKIP("spglib is not linked");
  }
  const auto ih = symm::crystal(symm::Crystal::iceIh);
  const auto ic = symm::crystal(symm::Crystal::iceIc);
  const auto fault = symm::crystal(symm::Crystal::fault);
  const symm::Space gIh = symm::spaceGroup(ih.vectors, ih.cartesian, 1e-3);
  const symm::Space gIc = symm::spaceGroup(ic.vectors, ic.cartesian, 1e-3);
  const symm::Space gFault = symm::spaceGroup(fault.vectors, fault.cartesian, 1e-3);
  REQUIRE(gIh.number == 194);
  REQUIRE(gIh.symbol.find("P6_3/mmc") != std::string::npos);
  REQUIRE(gIc.number == 227);
  REQUIRE(gIc.symbol.find("Fd-3m") != std::string::npos);
  REQUIRE(gFault.number != 194);
  REQUIRE(gFault.number != 227);
  REQUIRE(gFault.number == 15);
}

TEST_CASE("template names add what CHILL+ does not say", "[symmetry]") {
  const auto dod = symm::points(symm::Kind::cage512);
  const auto sI = symm::points(symm::Kind::cage51262);
  const auto sII = symm::points(symm::Kind::cage51264);
  REQUIRE(symm::best(dod).name == "5^12");
  REQUIRE(symm::best(sI).name == "5^12 6^2");
  REQUIRE(symm::best(sII).name == "5^12 6^4");
  const auto labDod = chillLabels(dod);
  const auto labSI = chillLabels(sI);
  const auto labSII = chillLabels(sII);
  for (int label : labDod) {
    REQUIRE(label == 2);
  }
  for (int label : labSI) {
    REQUIRE(label == 2);
  }
  for (int label : labSII) {
    REQUIRE(label == 2);
  }

  const auto ih = symm::points(symm::Kind::iceIh);
  const auto isd = symm::points(symm::Kind::iceIsd);
  const int centreIh = chillLabels(ih)[0];
  const int centreIsd = chillLabels(isd)[0];
  REQUIRE(symm::best(ih).name == "ice-Ih");
  REQUIRE(symm::best(isd).name == "ice-Isd");
  // CHILL+ has cubic, hexagonal, water, interfacial, clathrate. It has no
  // stacking-fault id. The fault shell is not labelled hexagonal.
  REQUIRE(centreIh == 1);
  REQUIRE(centreIsd != 1);
  REQUIRE(centreIsd != 0);
}

TEST_CASE("a periodic ice shell matches the same template as the ideal", "[symmetry]") {
  const auto ic = symm::crystal(symm::Crystal::iceIc);
  const auto ih = symm::crystal(symm::Crystal::iceIh);
  const Eigen::MatrixXd icShell = symm::iceShell(ic.cartesian, ic.vectors, 0);
  const Eigen::MatrixXd ihShell = symm::iceShell(ih.cartesian, ih.vectors, 0);
  REQUIRE(icShell.rows() == 17);
  REQUIRE(ihShell.rows() == 17);
  REQUIRE(symm::best(icShell).name == "ice-Ic");
  REQUIRE(symm::best(ihShell).name == "ice-Ih");
  REQUIRE(symm::best(icShell).rmsd < 1e-4);
  REQUIRE(symm::best(ihShell).rmsd < 1e-4);
}

TEST_CASE("the C API returns the template, the rotation and the space group",
          "[symmetry]") {
  const auto ref = symm::points(symm::Kind::ring6);
  const int n = static_cast<int>(ref.rows());
  std::vector<double> raw(static_cast<size_t>(3 * n));
  for (int i = 0; i < n; i++) {
    raw[static_cast<size_t>(3 * i)] = ref(i, 1);
    raw[static_cast<size_t>(3 * i + 1)] = ref(i, 0);
    raw[static_cast<size_t>(3 * i + 2)] = ref(i, 2);
  }
  int id = -1;
  double rmsd = -1;
  double haus = -1;
  double rotation[9];
  double translation[3];
  std::vector<int> perm(static_cast<size_t>(n), -1);
  REQUIRE(seams_symmetry_match(raw.data(), n, &id, &rmsd, &haus, rotation, translation,
                               perm.data()) == 0);
  REQUIRE(id == static_cast<int>(symm::Kind::ring6));
  REQUIRE(rmsd < 1e-6);
  char symbol[16] = {};
  int nops = 0;
  const int pg = seams_symmetry_point_group(raw.data(), n, 0.05, symbol, 16, &nops);
  if (symm::iraLinked()) {
    REQUIRE(pg == 0);
    REQUIRE(std::string(symbol) == "D6h");
    REQUIRE(nops == 24);
  } else {
    REQUIRE(pg == 2);
  }
  const auto cell = symm::crystal(symm::Crystal::iceIc);
  double lattice[9];
  for (int vec = 0; vec < 3; vec++) {
    for (int comp = 0; comp < 3; comp++) {
      lattice[3 * vec + comp] = cell.vectors(vec, comp);
    }
  }
  std::vector<double> cart(static_cast<size_t>(3 * cell.cartesian.rows()));
  for (int i = 0; i < cell.cartesian.rows(); i++) {
    cart[static_cast<size_t>(3 * i)] = cell.cartesian(i, 0);
    cart[static_cast<size_t>(3 * i + 1)] = cell.cartesian(i, 1);
    cart[static_cast<size_t>(3 * i + 2)] = cell.cartesian(i, 2);
  }
  int number = 0;
  char sg[16] = {};
  const int sp = seams_symmetry_space_group(lattice, cart.data(),
                                            static_cast<int>(cell.cartesian.rows()),
                                            1e-3, &number, sg, 16);
  if (symm::spglibLinked()) {
    REQUIRE(sp == 0);
    REQUIRE(number == 227);
  } else {
    REQUIRE(sp == 2);
  }
  int count = 0;
  REQUIRE(seams_symmetry_template(3, nullptr, 0, &count) == 0);
  REQUIRE(count == 20);
}
