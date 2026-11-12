// Time q4 and q6 on ideal ice Ic. One warm-up, then the repeats.
// Build: g++ -O3 -std=c++20 discovery/bench_ql.cpp -Ibindings/cxx
//        -L${SEAMS_BEFORE}/bindings -lseams_lang -o bench_ql

#include <seams_lang.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void ice_ic(int cells, std::vector<double> &xyz, double &side) {
  const double bond = 2.75;
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  side = cells * lattice;
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  xyz.clear();
  xyz.reserve(static_cast<std::size_t>(cells) * cells * cells * 8 * 3);
  for (int i = 0; i < cells; i++) {
    for (int j = 0; j < cells; j++) {
      for (int k = 0; k < cells; k++) {
        for (int b = 0; b < 8; b++) {
          const int f = b % 4;
          const double shift = b < 4 ? 0.0 : 0.25;
          const double x = std::fmod((i + fcc[f][0] + shift) * lattice, side);
          const double y = std::fmod((j + fcc[f][1] + shift) * lattice, side);
          const double z = std::fmod((k + fcc[f][2] + shift) * lattice, side);
          xyz.push_back(x);
          xyz.push_back(y);
          xyz.push_back(z);
        }
      }
    }
  }
}

double time_call(seams_frame *frame, int order, int repeats) {
  const int n = seams_frame_nop(frame);
  std::vector<double> ql(static_cast<std::size_t>(n));
  std::vector<double> bar(static_cast<std::size_t>(n));
  if (seams_frame_steinhardt(frame, 3.5, 1, order, ql.data(), bar.data()) != 0) {
    std::cerr << seams_last_error() << "\n";
    std::exit(1);
  }
  double acc = 0;
  for (int r = 0; r < repeats; r++) {
    const auto t0 = std::chrono::steady_clock::now();
    if (seams_frame_steinhardt(frame, 3.5, 1, order, ql.data(), bar.data()) != 0) {
      std::cerr << seams_last_error() << "\n";
      std::exit(1);
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << ms << "\n";
    acc += ms;
  }
  return acc / repeats;
}

} // namespace

int main(int argc, char **argv) {
  const int cells = argc > 1 ? std::atoi(argv[1]) : 24;
  const int repeats = argc > 2 ? std::atoi(argv[2]) : 3;
  const int order = argc > 3 ? std::atoi(argv[3]) : 6;
  std::vector<double> xyz;
  double side = 0;
  ice_ic(cells, xyz, side);
  const int n = static_cast<int>(xyz.size() / 3);
  std::vector<double> box = {0, 0, 0, side, side, side, 0, 0, 0};
  std::vector<int> types(static_cast<std::size_t>(n), 1);
  std::vector<int> ids(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    ids[static_cast<std::size_t>(i)] = i + 1;
  }
  seams_frame *frame = seams_frame_from_arrays(xyz.data(), n, box.data(), types.data(), ids.data());
  if (frame == nullptr) {
    std::cerr << "frame failed\n";
    return 1;
  }
  const double ms = time_call(frame, order, repeats);
  std::cout << "molecules=" << n << " l=" << order << " ms=" << ms << "\n";
  seams_frame_free(frame);
  return 0;
}
