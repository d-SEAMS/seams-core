#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <seams.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

struct IceFrame {
  std::vector<double> xyz;
  double box[9] = {};
  int n = 0;
};

// Eight molecules per orthorhombic cell. ncell = 8 is 4096 molecules.
IceFrame iceIh(int ncell) {
  const double bond = 2.76;
  const double a = bond * std::sqrt(8.0 / 3.0);
  const double c = a * std::sqrt(8.0 / 3.0);
  const double b = a * std::sqrt(3.0);
  const double basis[4][3] = {{1.0 / 3.0, 2.0 / 3.0, 0.0},
                              {2.0 / 3.0, 1.0 / 3.0, 0.5},
                              {1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0},
                              {2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0}};
  IceFrame ice;
  const double lx = ncell * a;
  const double ly = ncell * b;
  const double lz = ncell * c;
  ice.box[3] = lx;
  ice.box[4] = ly;
  ice.box[5] = lz;
  ice.n = ncell * ncell * ncell * 8;
  ice.xyz.resize(static_cast<std::size_t>(ice.n) * 3);
  int atom = 0;
  for (int ix = 0; ix < ncell; ++ix) {
    for (int iy = 0; iy < ncell; ++iy) {
      for (int iz = 0; iz < ncell; ++iz) {
        for (const double shift : {0.0, 0.5}) {
          for (const auto &site : basis) {
            double x = std::fmod(shift * a + site[0] * a - site[1] * a / 2.0 + ix * a, lx);
            double y = std::fmod(shift * b + site[1] * b / 2.0 + iy * b, ly);
            const double z = (site[2] + iz) * c;
            if (x < 0.0) {
              x += lx;
            }
            if (y < 0.0) {
              y += ly;
            }
            ice.xyz[static_cast<std::size_t>(3 * atom)] = x;
            ice.xyz[static_cast<std::size_t>(3 * atom + 1)] = y;
            ice.xyz[static_cast<std::size_t>(3 * atom + 2)] = z;
            ++atom;
          }
        }
      }
    }
  }
  return ice;
}

struct ThreadRestore {
  int saved;
  ThreadRestore() : saved(seams_get_num_threads()) {}
  ~ThreadRestore() { seams_set_num_threads(saved); }
};

}  // namespace

TEST_CASE("one thread and eight threads agree on knn, q6 and CHILL", "[omp][c_api]") {
  const ThreadRestore restore;
  const IceFrame ice = iceIh(8);
  seams_system *sys = seams_system_create();
  REQUIRE(sys != nullptr);
  std::vector<int> types(static_cast<std::size_t>(ice.n), 1);
  std::vector<int> ids(static_cast<std::size_t>(ice.n));
  for (int i = 0; i < ice.n; ++i) {
    ids[static_cast<std::size_t>(i)] = i + 1;
  }
  REQUIRE(seams_system_set_frame(sys, ice.xyz.data(), ice.n, ice.box, types.data(), ids.data()) ==
          SEAMS_OK);

  REQUIRE(seams_set_num_threads(1) == SEAMS_OK);
  seams_neighbours *serialCut = nullptr;
  REQUIRE(seams_neighbours_build(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, &serialCut) == SEAMS_OK);
  std::size_t cutNoff = 0;
  std::size_t cutNnz = 0;
  REQUIRE(seams_neighbours_copy(serialCut, nullptr, 0, &cutNoff, nullptr, 0, &cutNnz) == SEAMS_OK);
  std::vector<int> cutOff1(cutNoff);
  std::vector<int> cutCol1(cutNnz);
  REQUIRE(seams_neighbours_copy(serialCut, cutOff1.data(), cutOff1.size(), &cutNoff, cutCol1.data(),
                                cutCol1.size(), &cutNnz) == SEAMS_OK);
  REQUIRE(cutNnz == static_cast<std::size_t>(ice.n) * 4);
  for (int i = 0; i < ice.n; ++i) {
    REQUIRE(cutOff1[static_cast<std::size_t>(i) + 1] - cutOff1[static_cast<std::size_t>(i)] == 4);
  }

  seams_neighbours *serialNb = nullptr;
  REQUIRE(seams_neighbours_knn(sys, 4, 3.5, 1, 1, &serialNb) == SEAMS_OK);
  std::size_t serialNoff = 0;
  std::size_t serialNnz = 0;
  REQUIRE(seams_neighbours_copy(serialNb, nullptr, 0, &serialNoff, nullptr, 0, &serialNnz) ==
          SEAMS_OK);
  std::vector<int> off1(serialNoff);
  std::vector<int> col1(serialNnz);
  REQUIRE(seams_neighbours_copy(serialNb, off1.data(), off1.size(), &serialNoff, col1.data(),
                                col1.size(), &serialNnz) == SEAMS_OK);

  std::vector<double> q1(static_cast<std::size_t>(ice.n));
  std::vector<double> bar1(static_cast<std::size_t>(ice.n));
  REQUIRE(seams_steinhardt(sys, 3.5, 1, 6, q1.data(), q1.size(), bar1.data(), bar1.size()) ==
          SEAMS_OK);
  std::vector<int> lab1(static_cast<std::size_t>(ice.n));
  int counts1[9] = {};
  REQUIRE(seams_chill(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, lab1.data(),
                      lab1.size(), counts1) == SEAMS_OK);
  REQUIRE(counts1[SEAMS_LABEL_HEXAGONAL] == ice.n);
  std::vector<int> plus1(static_cast<std::size_t>(ice.n));
  int plusCounts1[9] = {};
  REQUIRE(seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_MUTUAL, SEAMS_CHILL_PUBLISHED, plus1.data(),
                           plusCounts1) == SEAMS_OK);
  REQUIRE(plusCounts1[SEAMS_LABEL_HEXAGONAL] == ice.n);

  seams_neighbours *serialOneWay = nullptr;
  REQUIRE(seams_neighbours_knn(sys, 4, 3.5, 1, 0, &serialOneWay) == SEAMS_OK);
  std::size_t oneNoff = 0;
  std::size_t oneNnz = 0;
  REQUIRE(seams_neighbours_copy(serialOneWay, nullptr, 0, &oneNoff, nullptr, 0, &oneNnz) ==
          SEAMS_OK);
  std::vector<int> oneOff(oneNoff);
  std::vector<int> oneCol(oneNnz);
  REQUIRE(seams_neighbours_copy(serialOneWay, oneOff.data(), oneOff.size(), &oneNoff, oneCol.data(),
                                oneCol.size(), &oneNnz) == SEAMS_OK);
  for (int i = 0; i < ice.n; ++i) {
    const int a = cutOff1[static_cast<std::size_t>(i)];
    const int b = cutOff1[static_cast<std::size_t>(i) + 1];
    const int c = oneOff[static_cast<std::size_t>(i)];
    const int d = oneOff[static_cast<std::size_t>(i) + 1];
    REQUIRE(b - a == d - c);
    std::vector<int> fromCut(cutCol1.begin() + a, cutCol1.begin() + b);
    std::vector<int> fromKnn(oneCol.begin() + c, oneCol.begin() + d);
    std::sort(fromCut.begin(), fromCut.end());
    std::sort(fromKnn.begin(), fromKnn.end());
    REQUIRE(fromCut == fromKnn);
  }

  REQUIRE(seams_set_num_threads(8) == SEAMS_OK);
  seams_neighbours *wideCut = nullptr;
  REQUIRE(seams_neighbours_build(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, &wideCut) == SEAMS_OK);
  std::size_t wideCutNoff = 0;
  std::size_t wideCutNnz = 0;
  REQUIRE(seams_neighbours_copy(wideCut, nullptr, 0, &wideCutNoff, nullptr, 0, &wideCutNnz) ==
          SEAMS_OK);
  REQUIRE(wideCutNoff == cutNoff);
  REQUIRE(wideCutNnz == cutNnz);
  std::vector<int> cutOff8(wideCutNoff);
  std::vector<int> cutCol8(wideCutNnz);
  REQUIRE(seams_neighbours_copy(wideCut, cutOff8.data(), cutOff8.size(), &wideCutNoff,
                                cutCol8.data(), cutCol8.size(), &wideCutNnz) == SEAMS_OK);
  REQUIRE(cutOff8 == cutOff1);
  REQUIRE(cutCol8 == cutCol1);
  seams_neighbours *wideNb = nullptr;
  REQUIRE(seams_neighbours_knn(sys, 4, 3.5, 1, 1, &wideNb) == SEAMS_OK);
  std::size_t wideNoff = 0;
  std::size_t wideNnz = 0;
  REQUIRE(seams_neighbours_copy(wideNb, nullptr, 0, &wideNoff, nullptr, 0, &wideNnz) == SEAMS_OK);
  REQUIRE(wideNoff == serialNoff);
  REQUIRE(wideNnz == serialNnz);
  std::vector<int> off8(wideNoff);
  std::vector<int> col8(wideNnz);
  REQUIRE(seams_neighbours_copy(wideNb, off8.data(), off8.size(), &wideNoff, col8.data(),
                                col8.size(), &wideNnz) == SEAMS_OK);
  REQUIRE(off8 == off1);
  REQUIRE(col8 == col1);

  std::vector<double> q8(static_cast<std::size_t>(ice.n));
  std::vector<double> bar8(static_cast<std::size_t>(ice.n));
  REQUIRE(seams_steinhardt(sys, 3.5, 1, 6, q8.data(), q8.size(), bar8.data(), bar8.size()) ==
          SEAMS_OK);
  for (int i = 0; i < ice.n; ++i) {
    REQUIRE(q8[static_cast<std::size_t>(i)] ==
            Catch::Approx(q1[static_cast<std::size_t>(i)]).margin(1e-12));
    REQUIRE(bar8[static_cast<std::size_t>(i)] ==
            Catch::Approx(bar1[static_cast<std::size_t>(i)]).margin(1e-12));
  }
  std::vector<int> lab8(static_cast<std::size_t>(ice.n));
  int counts8[9] = {};
  REQUIRE(seams_chill(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, lab8.data(),
                      lab8.size(), counts8) == SEAMS_OK);
  for (int i = 0; i < ice.n; ++i) {
    REQUIRE(lab8[static_cast<std::size_t>(i)] == lab1[static_cast<std::size_t>(i)]);
  }
  for (int label = 0; label < 9; ++label) {
    REQUIRE(counts8[label] == counts1[label]);
  }
  std::vector<int> plus8(static_cast<std::size_t>(ice.n));
  int plusCounts8[9] = {};
  REQUIRE(seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_MUTUAL, SEAMS_CHILL_PUBLISHED, plus8.data(),
                           plusCounts8) == SEAMS_OK);
  REQUIRE(plus8 == plus1);
  for (int label = 0; label < 9; ++label) {
    REQUIRE(plusCounts8[label] == plusCounts1[label]);
  }

  seams_neighbours *wideOneWay = nullptr;
  REQUIRE(seams_neighbours_knn(sys, 4, 3.5, 1, 0, &wideOneWay) == SEAMS_OK);
  std::size_t wideOneNoff = 0;
  std::size_t wideOneNnz = 0;
  REQUIRE(seams_neighbours_copy(wideOneWay, nullptr, 0, &wideOneNoff, nullptr, 0, &wideOneNnz) ==
          SEAMS_OK);
  REQUIRE(wideOneNoff == oneNoff);
  REQUIRE(wideOneNnz == oneNnz);
  std::vector<int> oneOff8(wideOneNoff);
  std::vector<int> oneCol8(wideOneNnz);
  REQUIRE(seams_neighbours_copy(wideOneWay, oneOff8.data(), oneOff8.size(), &wideOneNoff,
                                oneCol8.data(), oneCol8.size(), &wideOneNnz) == SEAMS_OK);
  REQUIRE(oneOff8 == oneOff);
  REQUIRE(oneCol8 == oneCol);

  seams_neighbours_destroy(serialCut);
  seams_neighbours_destroy(wideCut);
  seams_neighbours_destroy(serialOneWay);
  seams_neighbours_destroy(wideOneWay);
  seams_neighbours_destroy(serialNb);
  seams_neighbours_destroy(wideNb);
  seams_system_destroy(sys);
}

TEST_CASE("CHILL labels match at one and eight threads on the three ice cells",
          "[omp][c_api][chill]") {
  const ThreadRestore restore;
  for (const int ncell : {8, 16, 24}) {
    const IceFrame ice = iceIh(ncell);
    seams_system *sys = seams_system_create();
    REQUIRE(sys != nullptr);
    std::vector<int> types(static_cast<std::size_t>(ice.n), 1);
    std::vector<int> ids(static_cast<std::size_t>(ice.n));
    for (int i = 0; i < ice.n; ++i) {
      ids[static_cast<std::size_t>(i)] = i + 1;
    }
    REQUIRE(seams_system_set_frame(sys, ice.xyz.data(), ice.n, ice.box, types.data(),
                                   ids.data()) == SEAMS_OK);

    std::vector<int> lab1(static_cast<std::size_t>(ice.n));
    std::vector<int> lab8(static_cast<std::size_t>(ice.n));
    int counts1[9] = {};
    int counts8[9] = {};
    REQUIRE(seams_set_num_threads(1) == SEAMS_OK);
    REQUIRE(seams_chill(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, lab1.data(),
                        lab1.size(), counts1) == SEAMS_OK);
    REQUIRE(seams_set_num_threads(8) == SEAMS_OK);
    REQUIRE(seams_chill(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, lab8.data(),
                        lab8.size(), counts8) == SEAMS_OK);
    REQUIRE(lab8 == lab1);
    REQUIRE(counts1[SEAMS_LABEL_HEXAGONAL] == ice.n);
    for (int label = 0; label < 9; ++label) {
      REQUIRE(counts8[label] == counts1[label]);
    }
    for (int i = 0; i < ice.n; ++i) {
      REQUIRE(lab1[static_cast<std::size_t>(i)] == SEAMS_LABEL_HEXAGONAL);
    }

    std::vector<int> plus1(static_cast<std::size_t>(ice.n));
    std::vector<int> plus8(static_cast<std::size_t>(ice.n));
    int plusCounts1[9] = {};
    int plusCounts8[9] = {};
    REQUIRE(seams_set_num_threads(1) == SEAMS_OK);
    REQUIRE(seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                             plus1.data(), plusCounts1) == SEAMS_OK);
    REQUIRE(seams_set_num_threads(8) == SEAMS_OK);
    REQUIRE(seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                             plus8.data(), plusCounts8) == SEAMS_OK);
    REQUIRE(plus8 == plus1);
    REQUIRE(plusCounts1[SEAMS_LABEL_HEXAGONAL] == ice.n);

    std::vector<int> ov1(static_cast<std::size_t>(ice.n));
    std::vector<int> ov8(static_cast<std::size_t>(ice.n));
    int ovCounts1[9] = {};
    int ovCounts8[9] = {};
    REQUIRE(seams_set_num_threads(1) == SEAMS_OK);
    REQUIRE(seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, ov1.data(),
                             ovCounts1) == SEAMS_OK);
    REQUIRE(seams_set_num_threads(8) == SEAMS_OK);
    REQUIRE(seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, ov8.data(),
                             ovCounts8) == SEAMS_OK);
    REQUIRE(ov8 == ov1);
    REQUIRE(ovCounts1[SEAMS_LABEL_HEXAGONAL] == ice.n);
    for (int label = 0; label < 9; ++label) {
      REQUIRE(ovCounts8[label] == ovCounts1[label]);
    }
    seams_system_destroy(sys);
  }
}
