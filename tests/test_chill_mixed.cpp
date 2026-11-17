// Copyright (c) d-SEAMS developers
#include <catch2/catch_test_macros.hpp>
#include <seams.h>

#include <array>
#include <fstream>
#include <thread>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Frame {
  std::vector<double> xyz;
  double box[9] = {};
  int n = 0;
};

Frame readDump(const std::string &path) {
  std::ifstream in(path);
  REQUIRE(in.good());
  Frame frame;
  std::string line;
  while (std::getline(in, line)) {
    if (line == "ITEM: NUMBER OF ATOMS") {
      std::getline(in, line);
      frame.n = std::stoi(line);
    } else if (line.rfind("ITEM: BOX BOUNDS", 0) == 0) {
      for (int axis = 0; axis < 3; ++axis) {
        std::getline(in, line);
        std::istringstream row(line);
        double lo = 0;
        double hi = 0;
        row >> lo >> hi;
        frame.box[3 + axis] = hi - lo;
      }
    } else if (line.rfind("ITEM: ATOMS", 0) == 0) {
      frame.xyz.resize(static_cast<std::size_t>(frame.n) * 3);
      for (int i = 0; i < frame.n; ++i) {
        std::getline(in, line);
        std::istringstream row(line);
        int id = 0;
        int type = 0;
        double x = 0;
        double y = 0;
        double z = 0;
        row >> id >> type >> x >> y >> z;
        frame.xyz[static_cast<std::size_t>(3 * i)] = x;
        frame.xyz[static_cast<std::size_t>(3 * i + 1)] = y;
        frame.xyz[static_cast<std::size_t>(3 * i + 2)] = z;
      }
    }
  }
  REQUIRE(frame.n > 0);
  REQUIRE(static_cast<int>(frame.xyz.size()) == frame.n * 3);
  return frame;
}

std::vector<int> readLabels(const std::string &path) {
  std::ifstream in(path);
  REQUIRE(in.good());
  std::string line;
  std::string body;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (line.find('=') != std::string::npos) {
      continue;
    }
    body = line;
  }
  std::vector<int> labels;
  std::istringstream row(body);
  int value = 0;
  while (row >> value) {
    labels.push_back(value);
  }
  return labels;
}

void requireMatch(const Frame &frame, const std::vector<int> &reference, bool plus) {
  REQUIRE(static_cast<int>(reference.size()) == frame.n);
  seams_system *sys = seams_system_create();
  REQUIRE(sys != nullptr);
  std::vector<int> types(static_cast<std::size_t>(frame.n), 1);
  std::vector<int> ids(static_cast<std::size_t>(frame.n));
  for (int i = 0; i < frame.n; ++i) {
    ids[static_cast<std::size_t>(i)] = i + 1;
  }
  REQUIRE(seams_system_set_frame(sys, frame.xyz.data(), frame.n, frame.box, types.data(),
                                 ids.data()) == SEAMS_OK);
  for (const int threads : {1, 8}) {
    REQUIRE(seams_set_num_threads(threads) == SEAMS_OK);
    std::vector<int> got(static_cast<std::size_t>(frame.n));
    int counts[9] = {};
    if (plus) {
      REQUIRE(seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, got.data(),
                               counts) == SEAMS_OK);
    } else {
      REQUIRE(seams_chill(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, got.data(),
                          got.size(), counts) == SEAMS_OK);
    }
    REQUIRE(got == reference);
    int histogram[9] = {};
    for (int label : got) {
      REQUIRE(label >= 0);
      REQUIRE(label < 9);
      histogram[label] += 1;
    }
    for (int label = 0; label < 9; ++label) {
      REQUIRE(counts[label] == histogram[label]);
    }
  }
  seams_system_destroy(sys);
}

}  // namespace

TEST_CASE("mixed-phase CHILL labels match the pre-change implementation", "[chill][mixed]") {
  const char *names[] = {"interface", "stacking", "melted"};
  for (const char *name : names) {
    const std::string dir = "tests/fixtures/mixed/";
    const Frame frame = readDump(dir + name + ".lammpstrj");
    const std::vector<int> chill = readLabels(dir + name + ".chill.txt");
    const std::vector<int> plus = readLabels(dir + name + ".chill_plus.txt");
    requireMatch(frame, chill, false);
    requireMatch(frame, plus, true);
    seams_system *sys = seams_system_create();
    REQUIRE(sys != nullptr);
    std::vector<int> types(static_cast<std::size_t>(frame.n), 1);
    std::vector<int> ids(static_cast<std::size_t>(frame.n));
    for (int i = 0; i < frame.n; ++i) {
      ids[static_cast<std::size_t>(i)] = i + 1;
    }
    REQUIRE(seams_system_set_frame(sys, frame.xyz.data(), frame.n, frame.box, types.data(),
                                   ids.data()) == SEAMS_OK);
    std::vector<int> ovito1(static_cast<std::size_t>(frame.n));
    std::vector<int> ovito8(static_cast<std::size_t>(frame.n));
    int ovitoCounts[9] = {};
    REQUIRE(seams_set_num_threads(1) == SEAMS_OK);
    REQUIRE(seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, ovito1.data(),
                             ovitoCounts) == SEAMS_OK);
    REQUIRE(seams_set_num_threads(8) == SEAMS_OK);
    REQUIRE(seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, ovito8.data(),
                             ovitoCounts) == SEAMS_OK);
    REQUIRE(ovito8 == ovito1);
    seams_system_destroy(sys);
    if (std::string(name) == "stacking") {
      REQUIRE(ovito1 == plus);
    } else {
      REQUIRE(ovito1 != plus);
    }
    // The three frames are not a single phase. A regression that relabels
    // the interface, the stacking fault, or the melt changes this count.
    int classes = 0;
    bool seen[9] = {};
    for (int label : chill) {
      if (!seen[label]) {
        seen[label] = true;
        classes += 1;
      }
    }
    REQUIRE(classes >= 2);
  }
}

namespace {

std::vector<int> labelsWithTypes(const Frame &frame, const std::vector<int> &types, bool plus) {
  seams_system *sys = seams_system_create();
  REQUIRE(sys != nullptr);
  std::vector<int> ids(static_cast<std::size_t>(frame.n));
  for (int i = 0; i < frame.n; ++i) {
    ids[static_cast<std::size_t>(i)] = i + 1;
  }
  REQUIRE(seams_system_set_frame(sys, frame.xyz.data(), frame.n, frame.box, types.data(),
                                 ids.data()) == SEAMS_OK);
  std::vector<int> got(static_cast<std::size_t>(frame.n));
  int counts[9] = {};
  if (plus) {
    REQUIRE(seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_MUTUAL, SEAMS_CHILL_PUBLISHED, got.data(),
                             counts) == SEAMS_OK);
  } else {
    REQUIRE(seams_chill(sys, 3.5, 0, SEAMS_GRAPH_MUTUAL, SEAMS_CHILL_PUBLISHED, got.data(),
                        got.size(), counts) == SEAMS_OK);
  }
  seams_system_destroy(sys);
  return got;
}

}  // namespace

TEST_CASE("mutual-graph CHILL with type 0 does not depend on the type numbers",
          "[chill][mutual][types]") {
  const Frame frame = readDump("tests/fixtures/mixed/interface.lammpstrj");
  const std::size_t n = static_cast<std::size_t>(frame.n);
  const std::vector<int> ones(n, 1);
  const std::vector<int> twos(n, 2);
  std::vector<int> mixed(n, 1);
  for (std::size_t i = 1; i < n; i += 2) {
    mixed[i] = 2;
  }
  REQUIRE(seams_set_num_threads(1) == SEAMS_OK);
  for (const bool plus : {false, true}) {
    const std::vector<int> reference = labelsWithTypes(frame, ones, plus);
    bool varied = false;
    for (int label : reference) {
      varied = varied || label != reference[0];
    }
    REQUIRE(varied);
    // Every atom of type 2: the same single-species frame.
    REQUIRE(labelsWithTypes(frame, twos, plus) == reference);
    // Mixed numbers with type 0 search every atom, as the copied path does.
    REQUIRE(labelsWithTypes(frame, mixed, plus) == reference);
  }
}

namespace {

// Diamond lattice (cubic ice oxygens): nearest neighbours 2.75 A apart,
// the second shell 4.5 A away. The 3.5 A ring graph is bipartite.
Frame diamondFrame() {
  const double a = 6.35;
  const int cells = 3;
  const std::array<std::array<double, 3>, 8> basis = {{{0, 0, 0},
                                                       {0, 0.5, 0.5},
                                                       {0.5, 0, 0.5},
                                                       {0.5, 0.5, 0},
                                                       {0.25, 0.25, 0.25},
                                                       {0.25, 0.75, 0.75},
                                                       {0.75, 0.25, 0.75},
                                                       {0.75, 0.75, 0.25}}};
  Frame frame;
  for (int ix = 0; ix < cells; ++ix) {
    for (int iy = 0; iy < cells; ++iy) {
      for (int iz = 0; iz < cells; ++iz) {
        for (const auto &b : basis) {
          frame.xyz.push_back((ix + b[0]) * a);
          frame.xyz.push_back((iy + b[1]) * a);
          frame.xyz.push_back((iz + b[2]) * a);
        }
      }
    }
  }
  frame.n = static_cast<int>(frame.xyz.size() / 3);
  frame.box[3] = frame.box[4] = frame.box[5] = cells * a;
  return frame;
}

std::array<int, 13> ringHistogram(seams_system *sys) {
  seams_rings *rings = nullptr;
  REQUIRE(seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK);
  std::array<int, 13> hist = {};
  REQUIRE(seams_rings_histogram(rings, hist.data(), 13) == SEAMS_OK);
  seams_rings_destroy(rings);
  return hist;
}

seams_system *loaded(const Frame &frame) {
  seams_system *sys = seams_system_create();
  REQUIRE(sys != nullptr);
  REQUIRE(seams_system_set_frame(sys, frame.xyz.data(), frame.n, frame.box, nullptr, nullptr) ==
          SEAMS_OK);
  return sys;
}

}  // namespace

TEST_CASE("concurrent ring searches keep their own odd rings", "[rings][threads]") {
  const Frame crystal = diamondFrame();
  const Frame melt = readDump("tests/fixtures/mixed/melted.lammpstrj");
  seams_system *a = loaded(crystal);
  seams_system *b = loaded(melt);
  REQUIRE(seams_set_num_threads(2) == SEAMS_OK);
  const std::array<int, 13> crystalRef = ringHistogram(a);
  const std::array<int, 13> meltRef = ringHistogram(b);
  REQUIRE(crystalRef[6] > 0);
  REQUIRE(crystalRef[5] == 0);
  REQUIRE(meltRef[5] > 0);
  // Distinct handles may be used from distinct threads at once.
  std::vector<std::array<int, 13>> crystalRuns(40);
  std::vector<std::array<int, 13>> meltRuns(40);
  std::thread ta([&] {
    for (auto &h : crystalRuns) {
      seams_rings *rings = nullptr;
      if (seams_rings_find(a, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK) {
        seams_rings_histogram(rings, h.data(), 13);
        seams_rings_destroy(rings);
      }
    }
  });
  std::thread tb([&] {
    for (auto &h : meltRuns) {
      seams_rings *rings = nullptr;
      if (seams_rings_find(b, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK) {
        seams_rings_histogram(rings, h.data(), 13);
        seams_rings_destroy(rings);
      }
    }
  });
  ta.join();
  tb.join();
  for (const auto &h : crystalRuns) {
    REQUIRE(h == crystalRef);
  }
  for (const auto &h : meltRuns) {
    REQUIRE(h == meltRef);
  }
  seams_system_destroy(a);
  seams_system_destroy(b);
}
