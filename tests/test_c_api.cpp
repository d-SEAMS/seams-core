#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <diamond_cna.hpp>
#include <neighbours.hpp>
#include <seams_c_api.h>
#include <seams_input.hpp>

#include <cmath>
#include <vector>

TEST_CASE("seams_chill_plus labels a four-atom tetrahedron without crashing",
          "[c_api]") {
  const double xyz[12] = {0, 0, 0, 2.7, 0, 0, 1.35, 2.34, 0, 1.35, 0.78, 2.21};
  const double box[6] = {0, 0, 0, 20, 20, 20};
  int labels[4] = {-1, -1, -1, -1};
  REQUIRE(seams_chill_plus(xyz, 4, box, labels) == 0);
  for (int i = 0; i < 4; i++) {
    REQUIRE(labels[i] >= 0);
    REQUIRE(labels[i] <= 8);
  }
}

TEST_CASE("seams_chill_plus is seams_chill_plus_graph on the 3.5 A cutoff list",
          "[c_api]") {
  const double xyz[12] = {0, 0, 0, 2.7, 0, 0, 1.35, 2.34, 0, 1.35, 0.78, 2.21};
  const double box[6] = {0, 0, 0, 20, 20, 20};
  int byDefault[4] = {-1, -1, -1, -1};
  int byCutoff[4] = {-2, -2, -2, -2};
  int byNull[4] = {-3, -3, -3, -3};
  int byMutual[4] = {-4, -4, -4, -4};
  REQUIRE(seams_chill_plus(xyz, 4, box, byDefault) == 0);
  REQUIRE(seams_chill_plus_graph(xyz, 4, box, 3.5, "cutoff", byCutoff) == 0);
  REQUIRE(seams_chill_plus_graph(xyz, 4, box, 3.5, nullptr, byNull) == 0);
  REQUIRE(seams_chill_plus_graph(xyz, 4, box, 3.5, "mutual", byMutual) == 0);
  for (int i = 0; i < 4; i++) {
    REQUIRE(byDefault[i] == byCutoff[i]);
    REQUIRE(byNull[i] == byCutoff[i]);
    REQUIRE(byMutual[i] >= 0);
    REQUIRE(byMutual[i] <= 8);
  }
  int unused[4] = {};
  REQUIRE(seams_chill_plus_graph(xyz, 4, box, 3.5, "knn-union", unused) == 1);
  REQUIRE(seams_chill_plus_graph(xyz, 4, box, 0.0, "cutoff", unused) == 1);
  const double zero[6] = {0, 0, 0, 0, 0, 0};
  REQUIRE(seams_chill_plus_graph(xyz, 4, zero, 3.5, "cutoff", unused) == 1);
}

TEST_CASE("seams_chill_plus_graph_lammps matches the ortho call and accepts tilt",
          "[c_api]") {
  const double xyz[6] = {0, 0, 0, 2.0, 0, 0};
  const double ortho[6] = {0, 0, 0, 20, 20, 20};
  const double box9[9] = {0, 0, 0, 20, 20, 20, 0, 0, 0};
  int a[2] = {-1, -1};
  int b[2] = {-1, -1};
  REQUIRE(seams_chill_plus_graph(xyz, 2, ortho, 3.5, "cutoff", a) == 0);
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 2, box9, 3.5, "cutoff", b) == 0);
  REQUIRE(a[0] == b[0]);
  REQUIRE(a[1] == b[1]);
  const double tilted[9] = {0, 0, 0, 10, 8.660254037844386, 10, 5, 0, 0};
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 2, tilted, 3.5, "cutoff", b) == 0);
}

TEST_CASE("seams_q4_q6 returns finite cutoff values and rejects a bad box",
          "[c_api]") {
  const double xyz[12] = {0, 0, 0, 2.7, 0, 0, 1.35, 2.34, 0, 1.35, 0.78, 2.21};
  const double box[6] = {0, 0, 0, 20, 20, 20};
  const double box9[9] = {0, 0, 0, 20, 20, 20, 0, 0, 0};
  double q4[4] = {};
  double q6[4] = {};
  double q4bar[4] = {};
  double q6bar[4] = {};
  double q4b[4] = {};
  double q6b[4] = {};
  REQUIRE(seams_q4_q6(xyz, 4, box, 3.5, q4, q6, q4bar, q6bar) == 0);
  REQUIRE(seams_q4_q6_lammps(xyz, 4, box9, 3.5, q4b, q6b, nullptr, nullptr) == 0);
  for (int i = 0; i < 4; i++) {
    REQUIRE(std::isfinite(q4[i]));
    REQUIRE(std::isfinite(q6[i]));
    REQUIRE(std::isfinite(q4bar[i]));
    REQUIRE(std::isfinite(q6bar[i]));
    REQUIRE(q4[i] == q4b[i]);
    REQUIRE(q6[i] == q6b[i]);
  }
  REQUIRE(seams_q4_q6(nullptr, 4, box, 3.5, q4, q6, nullptr, nullptr) == 1);
  REQUIRE(seams_q4_q6(xyz, 4, box, 0.0, q4, q6, nullptr, nullptr) == 1);
  const double zero[6] = {0, 0, 0, 0, 0, 0};
  REQUIRE(seams_q4_q6(xyz, 4, zero, 3.5, q4, q6, nullptr, nullptr) == 1);
  REQUIRE(seams_q4_q6_lammps(xyz, 4, nullptr, 3.5, q4, q6, nullptr, nullptr) == 1);
}

namespace {

void iceIc(int nrep, std::vector<double> &xyz, double box9[9]) {
  const double bond = 2.75;
  const double a = 4.0 * bond / std::sqrt(3.0);
  const double basis[8][3] = {
      {0, 0, 0},         {0.5, 0.5, 0},     {0.5, 0, 0.5},     {0, 0.5, 0.5},
      {0.25, 0.25, 0.25},{0.75, 0.75, 0.25},{0.75, 0.25, 0.75},{0.25, 0.75, 0.75}};
  const double side = nrep * a;
  xyz.clear();
  for (int i = 0; i < nrep; i++) {
    for (int j = 0; j < nrep; j++) {
      for (int k = 0; k < nrep; k++) {
        for (int b = 0; b < 8; b++) {
          xyz.push_back((i + basis[b][0]) * a);
          xyz.push_back((j + basis[b][1]) * a);
          xyz.push_back((k + basis[b][2]) * a);
        }
      }
    }
  }
  box9[0] = box9[1] = box9[2] = 0;
  box9[3] = box9[4] = box9[5] = side;
  box9[6] = box9[7] = box9[8] = 0;
}

void iceIh(int nx, int ny, int nz, std::vector<double> &xyz, double box9[9]) {
  const double d = 2.75;
  const double z = 1.0 / 16.0;
  const double c = d / 0.375;
  const double a = c / std::sqrt(8.0 / 3.0);
  const double basis[4][3] = {{1.0 / 3.0, 2.0 / 3.0, z},
                              {2.0 / 3.0, 1.0 / 3.0, 0.5 + z},
                              {2.0 / 3.0, 1.0 / 3.0, 1.0 - z},
                              {1.0 / 3.0, 2.0 / 3.0, 0.5 - z}};
  const double lx = nx * a;
  const double ly = ny * a * std::sqrt(3.0) / 2.0;
  const double lz = nz * c;
  const double xy = ny * (-0.5 * a);
  box9[0] = box9[1] = box9[2] = 0;
  box9[3] = lx;
  box9[4] = ly;
  box9[5] = lz;
  box9[6] = xy;
  box9[7] = box9[8] = 0;
  xyz.clear();
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        for (int b = 0; b < 4; b++) {
          const double sx = (i + basis[b][0]) / nx;
          const double sy = (j + basis[b][1]) / ny;
          const double sz = (k + basis[b][2]) / nz;
          xyz.push_back(sx * lx + sy * xy);
          xyz.push_back(sy * ly);
          xyz.push_back(sz * lz);
        }
      }
    }
  }
}

} // namespace

TEST_CASE("cutoff list, q_l, diamond and primitive rings on ideal ice Ic",
          "[c_api]") {
  std::vector<double> xyz;
  double box9[9];
  iceIc(2, xyz, box9);
  const int n = static_cast<int>(xyz.size() / 3);
  REQUIRE(n == 64);
  std::vector<int> offsets(static_cast<std::size_t>(n + 1));
  int nCols = 0;
  REQUIRE(seams_cutoff_csr(xyz.data(), n, box9, 3.5, offsets.data(), nullptr, 0,
                           &nCols) == 0);
  REQUIRE(nCols == 4 * n);
  std::vector<int> cols(static_cast<std::size_t>(nCols));
  REQUIRE(seams_cutoff_csr(xyz.data(), n, box9, 3.5, offsets.data(), cols.data(),
                           nCols, &nCols) == 0);
  for (int i = 0; i < n; i++) {
    REQUIRE(offsets[static_cast<std::size_t>(i + 1)] -
                offsets[static_cast<std::size_t>(i)] ==
            4);
  }
  std::vector<double> q4(static_cast<std::size_t>(n));
  std::vector<double> q6(static_cast<std::size_t>(n));
  REQUIRE(seams_ql_cutoff(xyz.data(), n, box9, 3.5, 4, q4.data()) == 0);
  REQUIRE(seams_ql_cutoff(xyz.data(), n, box9, 3.5, 6, q6.data()) == 0);
  for (int i = 0; i < n; i++) {
    REQUIRE(q4[static_cast<std::size_t>(i)] == Catch::Approx(0.509175077).margin(1e-6));
    REQUIRE(q6[static_cast<std::size_t>(i)] == Catch::Approx(0.628539361).margin(1e-6));
  }
  std::vector<int> labels(static_cast<std::size_t>(n), -3);
  REQUIRE(seams_diamond_cna(xyz.data(), n, box9, 3.5, labels.data()) == 0);
  for (int label : labels) {
    REQUIRE(label == 1);
  }
  int counts[7] = {};
  REQUIRE(seams_primitive_rings(xyz.data(), n, box9, 3.5, 6, counts, 7) == 0);
  REQUIRE(counts[3] == 0);
  REQUIRE(counts[4] == 0);
  REQUIRE(counts[5] == 0);
  REQUIRE(counts[6] == 2 * n);
  REQUIRE(seams_ql_cutoff(xyz.data(), n, box9, 3.5, 5, q4.data()) == 1);
  REQUIRE(seams_cutoff_csr(nullptr, n, box9, 3.5, offsets.data(), nullptr, 0,
                           &nCols) == 1);
}

TEST_CASE("diamond common-neighbour test labels ideal ice Ih hexagonal",
          "[c_api]") {
  std::vector<double> xyz;
  double box9[9];
  iceIh(4, 4, 2, xyz, box9);
  const int n = static_cast<int>(xyz.size() / 3);
  REQUIRE(n == 128);
  std::vector<int> offsets(static_cast<std::size_t>(n + 1));
  int nCols = 0;
  REQUIRE(seams_cutoff_csr(xyz.data(), n, box9, 3.5, offsets.data(), nullptr, 0,
                           &nCols) == 0);
  REQUIRE(nCols == 4 * n);
  std::vector<int> labels(static_cast<std::size_t>(n), -3);
  REQUIRE(seams_diamond_cna(xyz.data(), n, box9, 3.5, labels.data()) == 0);
  for (int label : labels) {
    REQUIRE(label == 0);
  }
  int counts[7] = {};
  REQUIRE(seams_primitive_rings(xyz.data(), n, box9, 3.5, 6, counts, 7) == 0);
  REQUIRE(counts[3] == 0);
  REQUIRE(counts[4] == 0);
  REQUIRE(counts[5] == 0);
  REQUIRE(counts[6] == 2 * n);
}

TEST_CASE("crowded diamond labels match the cutoff list", "[c_api]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO("traj/dump-l.lammpstrj", 1, cloud, 1);
  REQUIRE(cloud.nop == 512);
  REQUIRE(cloud.box.size() >= 3);
  std::vector<double> xyz(static_cast<std::size_t>(cloud.nop) * 3);
  for (int i = 0; i < cloud.nop; i++) {
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(3 * i)] = p.x;
    xyz[static_cast<std::size_t>(3 * i + 1)] = p.y;
    xyz[static_cast<std::size_t>(3 * i + 2)] = p.z;
  }
  double box9[9] = {cloud.boxLow[0], cloud.boxLow[1], cloud.boxLow[2],
                    cloud.box[0],    cloud.box[1],    cloud.box[2],
                    0,               0,               0};
  std::vector<int> fast(static_cast<std::size_t>(cloud.nop), 9);
  REQUIRE(seams_diamond_cna(xyz.data(), cloud.nop, box9, 3.5, fast.data()) == 0);
  std::vector<int> off;
  std::vector<int> col;
  REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, 3.5, off, col, false));
  std::vector<int> full;
  symm::diamondFromCsr(cloud, off, col, full);
  REQUIRE(full.size() == fast.size());
  for (int i = 0; i < cloud.nop; i++) {
    REQUIRE(fast[static_cast<std::size_t>(i)] == full[static_cast<std::size_t>(i)]);
  }
}
