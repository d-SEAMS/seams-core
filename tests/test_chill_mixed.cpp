#include <catch2/catch_test_macros.hpp>
#include <seams.h>

#include <fstream>
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
