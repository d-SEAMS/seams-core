/*
** Split the k-nearest path: pack, lc_knearest, unpack, symmetrize.
**   profile_knn TRAJ [frame] [type] [reps]
*/
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <linkcell.h>
#include <neighbours.hpp>
#include <seams_input.hpp>

namespace {
using Clock = std::chrono::steady_clock;

double ms(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}
} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: profile_knn TRAJ [frame] [type] [reps]\n");
    return 2;
  }
  const std::string traj = argv[1];
  const int frame = argc > 2 ? std::atoi(argv[2]) : 1;
  const int typeI = argc > 3 ? std::atoi(argv[3]) : 1;
  const int reps = argc > 4 ? std::atoi(argv[4]) : 20;

  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO(traj, frame, cloud, typeI);
  if (cloud.nop == 0) {
    std::fprintf(stderr, "empty frame\n");
    return 1;
  }
  const int n = cloud.nop;
  const int k = 4;
  const double hint = 0.75 * 5.5;

  double tPack = 0, tLc = 0, tUnpack = 0, tSym = 0, tPair = 0;
  std::vector<std::vector<int>> nominated(static_cast<std::size_t>(n));

  for (int r = 0; r < reps; r++) {
    auto a0 = Clock::now();
    std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
    std::vector<int> mask(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; i++) {
      const auto &p = cloud.pts[static_cast<std::size_t>(i)];
      xyz[static_cast<std::size_t>(i) * 3 + 0] = p.x;
      xyz[static_cast<std::size_t>(i) * 3 + 1] = p.y;
      xyz[static_cast<std::size_t>(i) * 3 + 2] = p.z;
      mask[static_cast<std::size_t>(i)] = p.type == typeI ? 1 : 0;
    }
    lc_cell box = nneigh::lammpsBoxToLcCell(cloud.box, cloud.boxLow);
    std::vector<int> out(static_cast<std::size_t>(n) * static_cast<std::size_t>(k),
                         -1);
    auto a1 = Clock::now();
    if (lc_knearest(xyz.data(), n, &box, k, mask.data(), hint, out.data()) !=
        0) {
      std::fprintf(stderr, "%s\n", lc_last_error());
      return 1;
    }
    auto a2 = Clock::now();
    for (int i = 0; i < n; i++) {
      if (mask[static_cast<std::size_t>(i)] == 0) {
        continue;
      }
      auto &row = nominated[static_cast<std::size_t>(i)];
      row.clear();
      row.reserve(static_cast<std::size_t>(k));
      for (int t = 0; t < k; t++) {
        const int j =
            out[static_cast<std::size_t>(i) * static_cast<std::size_t>(k) +
                static_cast<std::size_t>(t)];
        if (j >= 0) {
          row.push_back(j);
        }
      }
    }
    auto a3 = Clock::now();
    auto both = nneigh::kNearestNeighbourPair(cloud, k, 5.5, typeI);
    auto a4 = Clock::now();
    (void)both;
    tPack += ms(a0, a1);
    tLc += ms(a1, a2);
    tUnpack += ms(a2, a3);
    tPair += ms(a3, a4);
  }

  // One extra pass: symmetrize only, after a fresh nominate via the public API.
  {
    auto nom = nominated;
    for (int r = 0; r < reps; r++) {
      auto s0 = Clock::now();
      auto both = nneigh::kNearestNeighbourPair(cloud, k, 5.5, typeI);
      (void)both;
      auto s1 = Clock::now();
      tSym += ms(s0, s1);
    }
  }

  std::printf("# nop %d reps %d\n", n, reps);
  std::printf("pack_ms %.3f\n", tPack / reps);
  std::printf("lc_knearest_ms %.3f\n", tLc / reps);
  std::printf("unpack_ms %.3f\n", tUnpack / reps);
  std::printf("pair_api_ms %.3f\n", tPair / reps);
  std::printf("pair_again_ms %.3f\n", tSym / reps);
  return 0;
}
