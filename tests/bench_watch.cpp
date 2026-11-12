// Release-only timer for the related-work watch.
// One warm-up, then the requested number of calls. The frame load is
// outside the timer. CHILL+ is seams_chill_plus (four mutual neighbours).
// Rings are the cutoff graph at 3.5 A, primitive, at most six members.

#include <seams_c_api.h>

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Frame {
  int n = 0;
  double side = 0;
  std::vector<double> xyz;
  double box[6] = {};
};

Frame load(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("cannot read " + path);
  }
  std::uint32_t n = 0;
  double side = 0;
  in.read(reinterpret_cast<char *>(&n), sizeof n);
  in.read(reinterpret_cast<char *>(&side), sizeof side);
  Frame frame;
  frame.n = static_cast<int>(n);
  frame.side = side;
  frame.xyz.resize(static_cast<std::size_t>(n) * 3);
  in.read(reinterpret_cast<char *>(frame.xyz.data()),
          static_cast<std::streamsize>(frame.xyz.size() * sizeof(double)));
  if (!in) {
    throw std::runtime_error("short frame " + path);
  }
  frame.box[3] = side;
  frame.box[4] = side;
  frame.box[5] = side;
  return frame;
}

molSys::PointCloud<molSys::Point<double>, double> cloudOf(const Frame &frame) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = frame.n;
  cloud.currentFrame = 1;
  cloud.box = {frame.side, frame.side, frame.side};
  cloud.boxLow = {0, 0, 0};
  cloud.pts.resize(static_cast<std::size_t>(frame.n));
  for (int i = 0; i < frame.n; i++) {
    auto &p = cloud.pts[static_cast<std::size_t>(i)];
    p.x = frame.xyz[3 * i];
    p.y = frame.xyz[3 * i + 1];
    p.z = frame.xyz[3 * i + 2];
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    cloud.idIndexMap[i + 1] = i;
  }
  return cloud;
}

std::vector<int> chill(const Frame &frame) {
  std::vector<int> labels(static_cast<std::size_t>(frame.n), -1);
  if (seams_chill_plus_flat(frame.xyz.data(), frame.n, frame.box, labels.data()) != 0) {
    throw std::runtime_error("seams_chill_plus failed");
  }
  return labels;
}

std::vector<int> rings(const Frame &frame) {
  auto cloud = cloudOf(frame);
  const auto rows = nneigh::neighListO(3.5, cloud, 1);
  const auto indexed = nneigh::neighbourListByIndex(cloud, rows);
  const auto found = primitive::ringNetwork(indexed, 6);
  std::vector<int> census(4, 0);
  for (const auto &ring : found) {
    const int size = static_cast<int>(ring.size());
    if (size >= 3 && size <= 6) {
      census[static_cast<std::size_t>(size - 3)] += 1;
    }
  }
  return census;
}

template <typename Fn>
std::vector<double> timeCalls(int warmup, int repeats, Fn &&fn) {
  using clock = std::chrono::steady_clock;
  for (int i = 0; i < warmup; i++) {
    fn();
  }
  std::vector<double> samples;
  samples.reserve(static_cast<std::size_t>(repeats));
  for (int i = 0; i < repeats; i++) {
    const auto t0 = clock::now();
    fn();
    const auto dt = std::chrono::duration<double, std::milli>(clock::now() - t0);
    samples.push_back(dt.count());
  }
  return samples;
}

void printSamples(const std::vector<double> &samples) {
  std::cout << "\"samples_ms\":[";
  for (std::size_t i = 0; i < samples.size(); i++) {
    if (i) {
      std::cout << ',';
    }
    std::cout.setf(std::ios::fixed);
    std::cout.precision(6);
    std::cout << samples[i];
  }
  std::cout << ']';
}

} // namespace

int main(int argc, char **argv) {
#ifndef SEAMS_WATCH_RELEASE
  std::cerr << "bench_watch requires meson -Dbuildtype=release\n";
  return 2;
#else
  if (argc == 2 && std::string(argv[1]) == "--build") {
    std::cout << "release\n";
    return 0;
  }
  if (argc != 5) {
    std::cerr << "usage: bench_watch chill|rings frame.bin warmup repeats\n";
    return 2;
  }
  const std::string operation = argv[1];
  const int warmup = std::stoi(argv[3]);
  const int repeats = std::stoi(argv[4]);
  const Frame frame = load(argv[2]);
  if (operation == "chill") {
    const auto samples = timeCalls(warmup, repeats, [&] { chill(frame); });
    const auto labels = chill(frame);
    std::cout << "{\"operation\":\"chill+\",\"build\":\"release\",\"molecules\":"
              << frame.n << ',';
    printSamples(samples);
    std::cout << ",\"labels\":[";
    for (int i = 0; i < frame.n; i++) {
      if (i) {
        std::cout << ',';
      }
      std::cout << labels[static_cast<std::size_t>(i)];
    }
    std::cout << "]}\n";
    return 0;
  }
  if (operation == "rings") {
    const auto samples = timeCalls(warmup, repeats, [&] { rings(frame); });
    const auto census = rings(frame);
    std::cout << "{\"operation\":\"rings\",\"build\":\"release\",\"molecules\":"
              << frame.n << ',';
    printSamples(samples);
    std::cout << ",\"census\":[" << census[0] << ',' << census[1] << ','
              << census[2] << ',' << census[3] << "]}\n";
    return 0;
  }
  std::cerr << "unknown operation " << operation << '\n';
  return 2;
#endif
}
