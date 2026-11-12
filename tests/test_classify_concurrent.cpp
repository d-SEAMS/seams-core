#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <local_symmetry.hpp>
#include <seams_c_api.h>

#include <omp.h>

#include <cmath>
#include <cstdio>
#include <thread>
#include <vector>

namespace {

Eigen::MatrixXd iceBlock(int nrep) {
  const auto cell = symm::crystal(symm::Crystal::iceIc);
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
  return big;
}

Eigen::Matrix3d iceBox(int nrep) {
  return static_cast<double>(nrep) * symm::crystal(symm::Crystal::iceIc).vectors;
}

void flatten(const Eigen::MatrixXd &xyz, std::vector<double> &raw, double lattice[9],
             const Eigen::Matrix3d &box) {
  raw.resize(static_cast<size_t>(xyz.rows() * 3));
  for (int i = 0; i < xyz.rows(); i++) {
    raw[static_cast<size_t>(3 * i)] = xyz(i, 0);
    raw[static_cast<size_t>(3 * i + 1)] = xyz(i, 1);
    raw[static_cast<size_t>(3 * i + 2)] = xyz(i, 2);
  }
  for (int vec = 0; vec < 3; vec++) {
    for (int comp = 0; comp < 3; comp++) {
      lattice[3 * vec + comp] = box(vec, comp);
    }
  }
}

// Raw cosine used by idealKind: the second-shell bond against another
// first-shell bond, with no projection. points() stores the centre, the
// four first neighbours, then three seconds of each first neighbour.
void rawSquares(const Eigen::MatrixXd &shell, double out[4]) {
  REQUIRE(shell.rows() == 17);
  double first[4][3];
  double b2[4];
  for (int a = 0; a < 4; a++) {
    for (int c = 0; c < 3; c++) {
      first[a][c] = shell(1 + a, c) - shell(0, c);
    }
    b2[a] = first[a][0] * first[a][0] + first[a][1] * first[a][1] + first[a][2] * first[a][2];
  }
  for (int a = 0; a < 4; a++) {
    double best = 0.0;
    for (int e = 0; e < 3; e++) {
      const int row = 5 + 3 * a + e;
      double vx = shell(row, 0) - shell(0, 0) - first[a][0];
      double vy = shell(row, 1) - shell(0, 1) - first[a][1];
      double vz = shell(row, 2) - shell(0, 2) - first[a][2];
      const double vn2 = vx * vx + vy * vy + vz * vz;
      for (int b = 0; b < 4; b++) {
        if (b == a) {
          continue;
        }
        const double dot = vx * first[b][0] + vy * first[b][1] + vz * first[b][2];
        if (dot <= 0.0) {
          continue;
        }
        const double ratio = (dot * dot) / (vn2 * b2[b]);
        if (ratio > best) {
          best = ratio;
        }
      }
    }
    out[a] = best;
  }
}

double minimg(double d, double L) {
  return d - L * std::round(d / L);
}

} // namespace

TEST_CASE("two classify calls run at once", "[symmetry]") {
  const int saved = omp_get_max_threads();
  omp_set_num_threads(1);

  const Eigen::MatrixXd small = iceBlock(1);
  const Eigen::MatrixXd large = iceBlock(2);
  const Eigen::Matrix3d boxS = iceBox(1);
  const Eigen::Matrix3d boxL = iceBox(2);

  std::vector<int> serialS;
  std::vector<int> serialL;
  symm::classifyIce(small, boxS, serialS);
  symm::classifyIce(large, boxL, serialL);
  REQUIRE(serialS.size() == static_cast<size_t>(small.rows()));
  REQUIRE(serialL.size() == static_cast<size_t>(large.rows()));
  for (int label : serialS) {
    REQUIRE(label == static_cast<int>(symm::Kind::iceIc));
  }
  for (int label : serialL) {
    REQUIRE(label == static_cast<int>(symm::Kind::iceIc));
  }

  std::vector<double> rawS;
  std::vector<double> rawL;
  double latS[9];
  double latL[9];
  flatten(small, rawS, latS, boxS);
  flatten(large, rawL, latL, boxL);

  constexpr int kLoops = 8;
  std::vector<int> gotS(serialS.size());
  std::vector<int> gotL(serialL.size());
  std::thread left([&] {
    for (int n = 0; n < kLoops; n++) {
      std::vector<int> kind;
      symm::classifyIce(small, boxS, kind);
      REQUIRE(kind == serialS);
      std::vector<int> labels(rawL.size() / 3, -9);
      REQUIRE(seams_symmetry_classify_tol(rawL.data(), static_cast<int>(rawL.size() / 3), latL,
                                          symm::kShellRmsd, labels.data()) == 0);
      REQUIRE(labels == serialL);
    }
    gotS = serialS;
  });
  std::thread right([&] {
    for (int n = 0; n < kLoops; n++) {
      std::vector<int> kind;
      symm::classifyIce(large, boxL, kind);
      REQUIRE(kind == serialL);
      std::vector<int> labels(rawS.size() / 3, -9);
      REQUIRE(seams_symmetry_classify_tol(rawS.data(), static_cast<int>(rawS.size() / 3), latS,
                                          symm::kShellRmsd, labels.data()) == 0);
      REQUIRE(labels == serialS);
    }
    gotL = serialL;
  });
  left.join();
  right.join();
  REQUIRE(gotS == serialS);
  REQUIRE(gotL == serialL);
  omp_set_num_threads(saved);
}

TEST_CASE("ideal ice cosines are the raw dot product and the azimuthal projection",
          "[symmetry]") {
  const Eigen::MatrixXd ic = symm::points(symm::Kind::iceIc);
  const Eigen::MatrixXd ih = symm::points(symm::Kind::iceIh);
  double azIc[4];
  double azIh[4];
  REQUIRE(symm::bondCosines(ic, azIc) == 4);
  REQUIRE(symm::bondCosines(ih, azIh) == 4);
  double rawIc[4];
  double rawIh[4];
  rawSquares(ic, rawIc);
  rawSquares(ih, rawIh);

  // bondCosines sorts high to low. idealKind's raw square is per arm, in
  // template order: ice Ih eclipses the first arm.
  std::printf("azimuthal ice-Ic (projection onto the plane normal to the bond,"
              " armMaxCosine and kEclipsedCosine):");
  for (int a = 0; a < 4; a++) {
    std::printf(" %.6f", azIc[a]);
  }
  std::printf("\n");
  std::printf("azimuthal ice-Ih (same projection):");
  for (int a = 0; a < 4; a++) {
    std::printf(" %.6f", azIh[a]);
  }
  std::printf("\n");
  std::printf("raw ice-Ic (idealKind, second-shell bond against another first-shell bond):");
  for (int a = 0; a < 4; a++) {
    std::printf(" cos %.6f cos2 %.6f", std::sqrt(rawIc[a]), rawIc[a]);
  }
  std::printf("\n");
  std::printf("raw ice-Ih (idealKind, same raw dot product):");
  for (int a = 0; a < 4; a++) {
    std::printf(" cos %.6f cos2 %.6f", std::sqrt(rawIh[a]), rawIh[a]);
  }
  std::printf("\n");
  std::fflush(stdout);

  for (int a = 0; a < 4; a++) {
    REQUIRE_THAT(azIc[a], Catch::Matchers::WithinAbs(0.5, 1e-9));
    REQUIRE_THAT(std::sqrt(rawIc[a]), Catch::Matchers::WithinAbs(1.0 / 3.0, 1e-9));
    REQUIRE(rawIc[a] <= 0.25);
    REQUIRE(rawIc[a] < 0.36);
  }
  REQUIRE_THAT(azIh[0], Catch::Matchers::WithinAbs(1.0, 1e-9));
  for (int a = 1; a < 4; a++) {
    REQUIRE_THAT(azIh[a], Catch::Matchers::WithinAbs(0.5, 1e-9));
  }
  REQUIRE_THAT(std::sqrt(rawIh[0]), Catch::Matchers::WithinAbs(7.0 / 9.0, 1e-9));
  REQUIRE(rawIh[0] >= 0.36);
  for (int a = 1; a < 4; a++) {
    REQUIRE_THAT(std::sqrt(rawIh[a]), Catch::Matchers::WithinAbs(1.0 / 3.0, 1e-9));
    REQUIRE(rawIh[a] <= 0.25);
  }
}

TEST_CASE("a dense shell keeps the crystal neighbours", "[symmetry]") {
  const int saved = omp_get_max_threads();
  omp_set_num_threads(1);

  const Eigen::MatrixXd ice = iceBlock(3);
  const int nIce = static_cast<int>(ice.rows());
  constexpr int kSide = 9;
  constexpr double kSpace = 1.2;
  const int nDense = kSide * kSide * kSide;
  const int n = nIce + nDense;
  Eigen::MatrixXd frame(n, 3);
  const Eigen::Vector3d iceShift(40.0, 40.0, 40.0);
  for (int i = 0; i < nIce; i++) {
    frame.row(i) = ice.row(i) + iceShift.transpose();
  }
  int filled = nIce;
  for (int ix = 0; ix < kSide; ix++) {
    for (int iy = 0; iy < kSide; iy++) {
      for (int iz = 0; iz < kSide; iz++) {
        frame(filled, 0) = ix * kSpace;
        frame(filled, 1) = iy * kSpace;
        frame(filled, 2) = iz * kSpace;
        filled++;
      }
    }
  }
  REQUIRE(filled == n);

  const double side = 80.0;
  Eigen::Matrix3d box = Eigen::Matrix3d::Identity() * side;

  std::vector<int> alone;
  Eigen::MatrixXd iceOnly = frame.topRows(nIce);
  symm::classifyIce(iceOnly, box, alone);
  int nCrystal = 0;
  for (int label : alone) {
    if (label == static_cast<int>(symm::Kind::iceIc)) {
      nCrystal++;
    }
  }
  REQUIRE(nCrystal > 0);

  int pairs = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      const double dx = minimg(frame(j, 0) - frame(i, 0), side);
      const double dy = minimg(frame(j, 1) - frame(i, 1), side);
      const double dz = minimg(frame(j, 2) - frame(i, 2), side);
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 >= 2.15 * 2.15 && d2 <= 3.5 * 3.5) {
        pairs++;
      }
    }
  }
  REQUIRE(pairs > n * 16);

  std::vector<int> both;
  symm::classifyIce(frame, box, both);
  REQUIRE(static_cast<int>(both.size()) == n);
  for (int i = 0; i < nIce; i++) {
    if (alone[static_cast<size_t>(i)] == static_cast<int>(symm::Kind::iceIc)) {
      REQUIRE(both[static_cast<size_t>(i)] == static_cast<int>(symm::Kind::iceIc));
    }
  }
  omp_set_num_threads(saved);
}
