#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diamond_cna.hpp>
#include <local_symmetry.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_c_api.h>

#include <Eigen/Geometry>

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <random>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <system_error>
#include <unistd.h>
#include <vector>

namespace symm {
int diamondCna(const double *sx, const double *sy, const double *sz, int ambiguous, int *n421 = nullptr,
               int *n422 = nullptr, int *nGood = nullptr);
}

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
  REQUIRE(seams_chill_plus_flat(xyz.data(), n, box, labels.data()) == 0);
  return labels;
}

// One file per call. Parallel meson test must not share a name, and the
// file is removed when the owner goes out of scope, including on failure.
struct ScratchFile {
  std::filesystem::path path;

  explicit ScratchFile(const char *tag) {
    static std::atomic<unsigned> seq{0};
    const unsigned id = seq.fetch_add(1, std::memory_order_relaxed);
    path = std::filesystem::temp_directory_path() /
           ("csm-" + std::to_string(getpid()) + "-" + std::to_string(id) +
            "-" + tag);
  }

  ScratchFile(const ScratchFile &) = delete;
  ScratchFile &operator=(const ScratchFile &) = delete;

  ~ScratchFile() {
    std::error_code ec;
    std::filesystem::remove(path, ec);
  }

  std::string str() const { return path.string(); }
};

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

TEST_CASE("a cubic ice Ic supercell has the conventional space group", "[symmetry]") {
  if (!symm::spglibLinked()) {
    SKIP("spglib is not linked");
  }
  const auto cell = symm::crystal(symm::Crystal::iceIc);
  const int m = 3;
  const int n0 = static_cast<int>(cell.cartesian.rows());
  symm::Lattice big;
  big.vectors = cell.vectors * static_cast<double>(m);
  big.cartesian.resize(static_cast<Eigen::Index>(n0 * m * m * m), 3);
  int w = 0;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      for (int k = 0; k < m; k++) {
        const Eigen::RowVector3d shift =
            static_cast<double>(i) * cell.vectors.row(0) +
            static_cast<double>(j) * cell.vectors.row(1) +
            static_cast<double>(k) * cell.vectors.row(2);
        for (int a = 0; a < n0; a++) {
          big.cartesian.row(w++) = cell.cartesian.row(a) + shift;
        }
      }
    }
  }
  std::vector<int> kind(static_cast<std::size_t>(big.cartesian.rows()),
                        static_cast<int>(symm::Kind::iceIc));
  symm::Space fast;
  REQUIRE(symm::conventionalCrystalSpace(big.vectors, big.cartesian, kind, 1e-3, fast));
  const symm::Space full = symm::spaceGroup(big.vectors, big.cartesian, 1e-3);
  REQUIRE(fast.status == 0);
  REQUIRE(fast.number == full.number);
  REQUIRE(fast.symbol == full.symbol);
  REQUIRE(fast.number == 227);
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

TEST_CASE("the cell list labels every crystal oxygen the way the shell match does",
          "[symmetry]") {
  if (!symm::iraLinked()) {
    SKIP("libira is not linked");
  }
  const symm::Crystal cells[] = {symm::Crystal::iceIh, symm::Crystal::iceIc};
  for (symm::Crystal which : cells) {
    const auto cell = symm::crystal(which);
    std::vector<int> kind;
    symm::classifyIce(cell.cartesian, cell.vectors, kind);
    REQUIRE(static_cast<int>(kind.size()) == cell.cartesian.rows());
    for (int i = 0; i < cell.cartesian.rows(); i++) {
      const auto shell = symm::iceShell(cell.cartesian, cell.vectors, i);
      INFO(symm::name(which == symm::Crystal::fault ? symm::Kind::iceIsd
                                                     : symm::Kind::iceIh)
           << " atom " << i);
      if (shell.rows() != 17) {
        REQUIRE(kind[static_cast<size_t>(i)] == -1);
        continue;
      }
      const auto ref = symm::best(shell, "ira");
      REQUIRE(ref.status == 0);
      REQUIRE(kind[static_cast<size_t>(i)] == static_cast<int>(ref.kind));
    }
  }
}

TEST_CASE("the stacking-fault cell is ice Ic and ice Isd under the reference match",
          "[symmetry]") {
  if (!symm::iraLinked()) {
    SKIP("libira is not linked");
  }
  const auto cell = symm::crystal(symm::Crystal::fault);
  int nIc = 0;
  int nIsd = 0;
  for (int i = 0; i < cell.cartesian.rows(); i++) {
    const auto shell = symm::iceShell(cell.cartesian, cell.vectors, i);
    REQUIRE(shell.rows() == 17);
    const auto ref = symm::best(shell, "ira");
    if (ref.kind == symm::Kind::iceIc) {
      nIc++;
    }
    if (ref.kind == symm::Kind::iceIsd) {
      nIsd++;
    }
  }
  REQUIRE(nIc == 4);
  REQUIRE(nIsd == 4);
}

TEST_CASE("an ice Ic supercell is cubic at every oxygen", "[symmetry]") {
  const auto cell = symm::crystal(symm::Crystal::iceIc);
  const int nrep = 3;
  const int base = static_cast<int>(cell.cartesian.rows());
  Eigen::MatrixXd big(base * nrep * nrep * nrep, 3);
  int row = 0;
  for (int ix = 0; ix < nrep; ix++) {
    for (int iy = 0; iy < nrep; iy++) {
      for (int iz = 0; iz < nrep; iz++) {
        const Eigen::Vector3d shift =
            static_cast<double>(ix) * cell.vectors.row(0).transpose() +
            static_cast<double>(iy) * cell.vectors.row(1).transpose() +
            static_cast<double>(iz) * cell.vectors.row(2).transpose();
        for (int i = 0; i < base; i++) {
          big.row(row++) = cell.cartesian.row(i) + shift.transpose();
        }
      }
    }
  }
  Eigen::Matrix3d box = cell.vectors * static_cast<double>(nrep);
  std::vector<int> kind;
  symm::classifyIce(big, box, kind);
  REQUIRE(static_cast<int>(kind.size()) == big.rows());
  for (int label : kind) {
    REQUIRE(label == static_cast<int>(symm::Kind::iceIc));
  }
}

TEST_CASE("ideal Ic, Ih and a gaussian displacement keep their labels", "[symmetry]") {
  std::mt19937 rng(20261113);
  std::normal_distribution<double> noise(0.0, 0.05);
  const symm::Crystal cells[] = {symm::Crystal::iceIc, symm::Crystal::iceIh};
  const int expect[] = {static_cast<int>(symm::Kind::iceIc),
                        static_cast<int>(symm::Kind::iceIh)};
  for (int which = 0; which < 2; ++which) {
    const auto cell = symm::crystal(cells[which]);
    std::vector<int> kind;
    symm::classifyIce(cell.cartesian, cell.vectors, kind);
    for (int label : kind) {
      REQUIRE(label == expect[which]);
    }
    Eigen::MatrixXd moved = cell.cartesian;
    for (int row = 0; row < moved.rows(); ++row) {
      for (int axis = 0; axis < 3; ++axis) {
        moved(row, axis) += noise(rng);
      }
    }
    symm::classifyIce(moved, cell.vectors, kind);
    for (int label : kind) {
      REQUIRE(label == expect[which]);
    }
  }
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

TEST_CASE("the eclipsed cutoff sits in the gap between staggered and eclipsed",
          "[symmetry]") {
  const auto ic = symm::crystal(symm::Crystal::iceIc);
  const auto ih = symm::crystal(symm::Crystal::iceIh);
  std::vector<int> kind;
  symm::classifyIce(ic.cartesian, ic.vectors, kind);
  REQUIRE(static_cast<int>(kind.size()) == ic.cartesian.rows());
  for (int label : kind) {
    REQUIRE(label == static_cast<int>(symm::Kind::iceIc));
  }
  symm::classifyIce(ic.cartesian, ic.vectors, kind, 0.05);
  for (int label : kind) {
    REQUIRE(label == static_cast<int>(symm::Kind::iceIc));
  }
  symm::classifyIce(ic.cartesian, ic.vectors, kind, symm::kShellRmsd);
  for (int label : kind) {
    REQUIRE(label == static_cast<int>(symm::Kind::iceIc));
  }

  symm::classifyIce(ih.cartesian, ih.vectors, kind);
  for (int label : kind) {
    REQUIRE(label == static_cast<int>(symm::Kind::iceIh));
  }
  symm::classifyIce(ih.cartesian, ih.vectors, kind, 0.05);
  for (int label : kind) {
    REQUIRE(label == static_cast<int>(symm::Kind::iceIh));
  }

  std::vector<double> raw(static_cast<size_t>(3 * ih.cartesian.rows()));
  double lattice[9];
  for (int i = 0; i < ih.cartesian.rows(); i++) {
    raw[static_cast<size_t>(3 * i)] = ih.cartesian(i, 0);
    raw[static_cast<size_t>(3 * i + 1)] = ih.cartesian(i, 1);
    raw[static_cast<size_t>(3 * i + 2)] = ih.cartesian(i, 2);
  }
  for (int vec = 0; vec < 3; vec++) {
    for (int comp = 0; comp < 3; comp++) {
      lattice[3 * vec + comp] = ih.vectors(vec, comp);
    }
  }
  std::vector<int> labels(static_cast<size_t>(ih.cartesian.rows()), 9);
  REQUIRE(seams_symmetry_classify_tol(raw.data(), static_cast<int>(ih.cartesian.rows()),
                                      lattice, 0.05, labels.data()) == 0);
  for (int label : labels) {
    REQUIRE(label == static_cast<int>(symm::Kind::iceIh));
  }
  REQUIRE(seams_symmetry_classify(raw.data(), static_cast<int>(ih.cartesian.rows()), lattice,
                                  labels.data()) == 0);
  for (int label : labels) {
    REQUIRE(label == static_cast<int>(symm::Kind::iceIh));
  }

  const Eigen::MatrixXd icShell = symm::iceShell(ic.cartesian, ic.vectors, 0);
  REQUIRE(symm::iceTemplateRmsd(icShell, symm::Kind::iceIc) < 1e-6);
  REQUIRE(symm::iceTemplateRmsd(icShell, symm::Kind::iceIh) > 0.05);
  REQUIRE(symm::iceTemplateRmsd(icShell, symm::Kind::iceIsd) > 0.05);
  double cosine[4];
  REQUIRE(symm::bondCosines(icShell, cosine) == 4);
  REQUIRE(cosine[0] < 0.6);
  REQUIRE(cosine[3] > 0.4);
  const Eigen::MatrixXd ihShell = symm::iceShell(ih.cartesian, ih.vectors, 0);
  REQUIRE(symm::bondCosines(ihShell, cosine) == 4);
  REQUIRE(cosine[0] > 0.95);
  REQUIRE(cosine[1] < 0.6);
}

TEST_CASE("parallel naming matches the serial shell, cage and ring results", "[symmetry]") {
  const symm::Kind kinds[] = {symm::Kind::iceIc,  symm::Kind::iceIh, symm::Kind::iceIsd,
                              symm::Kind::cage512, symm::Kind::cage51262,
                              symm::Kind::cage51264, symm::Kind::ring5, symm::Kind::ring6};
  std::vector<Eigen::MatrixXd> envs;
  for (int copy = 0; copy < 4; copy++) {
    for (symm::Kind kind : kinds) {
      envs.push_back(symm::points(kind));
    }
  }
  REQUIRE(envs.size() >= 32);
  std::vector<symm::Overlay> named;
  symm::matchAll(envs, named);
  REQUIRE(named.size() == envs.size());
  for (size_t i = 0; i < envs.size(); i++) {
    const symm::Overlay one = symm::matchKernel(envs[i]);
    REQUIRE(named[i].status == 0);
    REQUIRE(one.status == 0);
    REQUIRE(named[i].name == one.name);
    REQUIRE(std::abs(named[i].rmsd - one.rmsd) < 1e-9);
  }
}

void writeXyz(const std::string &path, const Eigen::MatrixXd &xyz) {
  std::ofstream out(path);
  out << xyz.rows() << "\n\n";
  for (Eigen::Index i = 0; i < xyz.rows(); i++) {
    out << "O " << xyz(i, 0) << " " << xyz(i, 1) << " " << xyz(i, 2) << "\n";
  }
}

double shapeProgram(const std::string &bin, const std::string &studied,
                    const std::string &reference) {
  const ScratchFile out("shape-out.txt");
  const std::string cmd = bin + " " + studied + " 0 " + reference + " 1 1 > " + out.str();
  REQUIRE(std::system(cmd.c_str()) == 0);
  std::ifstream in(out.path);
  std::string line;
  double value = -1.0;
  while (std::getline(in, line)) {
    if (line.rfind("SHAPE", 0) == 0) {
      std::istringstream row(line);
      std::string word;
      row >> word >> word >> value;
    }
  }
  REQUIRE(value >= 0.0);
  return value;
}

TEST_CASE("the in-tree shape measure matches the shape program", "[symmetry]") {
  const char *bin = std::getenv("CSM_SHAPE");
  if (bin == nullptr || bin[0] == '\0') {
    SKIP("shape program was not built");
  }
  const Eigen::MatrixXd ref = symm::points(symm::Kind::iceIc);
  Eigen::MatrixXd moved = ref;
  moved(7, 0) += 0.35;
  moved(7, 1) -= 0.12;
  moved(7, 2) += 0.08;
  const symm::Overlay hit = symm::matchKernel(moved);
  REQUIRE(hit.status == 0);
  REQUIRE(hit.kind == symm::Kind::iceIc);
  Eigen::MatrixXd paired(moved.rows(), 3);
  for (Eigen::Index i = 0; i < moved.rows(); i++) {
    paired.row(i) = moved.row(hit.permutation[static_cast<size_t>(i)]);
  }
  const ScratchFile refPath("shell-ref.xyz");
  const ScratchFile movedPath("shell-moved.xyz");
  writeXyz(refPath.str(), ref);
  writeXyz(movedPath.str(), paired);
  const double theirs = shapeProgram(bin, movedPath.str(), refPath.str());
  const double ours = symm::shapeMeasure(moved);
  REQUIRE(std::abs(ours - theirs) < 0.05);
}

// FCC second shell, nearest-neighbour length sqrt(2). The first four sites
// are pushed out by 1.2. Each of those four loses the (4,2,1) pattern, and
// the other eight keep it. The production ambiguous count is 9, so an early
// exit on "cannot reach 9" returns before those eight are counted.
TEST_CASE("eight good diamond signatures stay a partial cubic label", "[symmetry]") {
  const double ideal[12][3] = {
      {-1, -1, 0}, {-1, 1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 0, -1}, {-1, 0, 1},
      {1, 0, -1},  {1, 0, 1},  {0, -1, -1}, {0, -1, 1}, {0, 1, -1}, {0, 1, 1},
  };
  double sx[12], sy[12], sz[12];
  for (int i = 0; i < 12; i++) {
    const double scale = i < 4 ? 1.2 : 1.0;
    sx[i] = ideal[i][0] * scale;
    sy[i] = ideal[i][1] * scale;
    sz[i] = ideal[i][2] * scale;
  }
  int n421 = -1;
  int n422 = -1;
  int nGood = -1;
  const int partial = symm::diamondCna(sx, sy, sz, 9, &n421, &n422, &nGood);
  REQUIRE(nGood == 8);
  REQUIRE(n421 == 8);
  REQUIRE(n422 == 0);
  REQUIRE(partial == 3);

  for (int i = 0; i < 12; i++) {
    sx[i] = ideal[i][0];
    sy[i] = ideal[i][1];
    sz[i] = ideal[i][2];
  }
  n421 = n422 = nGood = -1;
  const int cubic = symm::diamondCna(sx, sy, sz, 9, &n421, &n422, &nGood);
  REQUIRE(nGood == 12);
  REQUIRE(n421 == 12);
  REQUIRE(cubic == 1);
}

TEST_CASE("shell kernel matches diamondFromCsr on ideal cubic ice", "[symmetry]") {
  const double oo = 2.75;
  const double a = oo * 4.0 / std::sqrt(3.0);
  const int nrep = 2;
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double length = nrep * a;
  cloud.box = {length, length, length};
  cloud.boxLow = {0, 0, 0};
  for (int i = 0; i < nrep; ++i) {
    for (int j = 0; j < nrep; ++j) {
      for (int k = 0; k < nrep; ++k) {
        for (const auto &shift : fcc) {
          for (const auto &site : basis) {
            molSys::Point<double> p;
            p.x = (i + shift[0] + site[0]) * a;
            p.y = (j + shift[1] + site[1]) * a;
            p.z = (k + shift[2] + site[2]) * a;
            p.type = 1;
            p.atomID = cloud.nop + 1;
            cloud.idIndexMap[p.atomID] = cloud.nop;
            cloud.pts.push_back(p);
            cloud.nop += 1;
          }
        }
      }
    }
  }
  std::vector<int> offsets;
  std::vector<int> cols;
  REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, 3.5, offsets, cols, false));
  std::vector<int> diamond;
  symm::diamondFromCsr(cloud, offsets, cols, diamond);
  std::vector<int> shell;
  symm::classifyIceCloud(cloud, shell);
  REQUIRE(shell == diamond);
}
