/*
** Time the hydrate signature search the way the language binding does:
** a cutoff neighbour list, primitive rings, and one closed-cage search
** for 5^12, then again for 5^12 6^2, then again for 5^12 6^4.
** The replicate is built in memory from the unit-cell dump.
*/
#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

double msSince(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - t0)
      .count();
}

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

Cloud tile(Cloud src, int n) {
  if (n <= 1) {
    return src;
  }
  const double lx = src.box.size() > 0 ? src.box[0] : 0.0;
  const double ly = src.box.size() > 1 ? src.box[1] : 0.0;
  const double lz = src.box.size() > 2 ? src.box[2] : 0.0;
  Cloud out;
  out.box = src.box;
  if (out.box.size() >= 1) {
    out.box[0] = lx * n;
  }
  if (out.box.size() >= 2) {
    out.box[1] = ly * n;
  }
  if (out.box.size() >= 3) {
    out.box[2] = lz * n;
  }
  out.boxLow = src.boxLow;
  const int base = src.nop;
  out.nop = base * n * n * n;
  out.pts.resize(static_cast<std::size_t>(out.nop));
  int at = 0;
  for (int ix = 0; ix < n; ++ix) {
    for (int iy = 0; iy < n; ++iy) {
      for (int iz = 0; iz < n; ++iz) {
        for (int i = 0; i < base; ++i) {
          auto point = src.pts[static_cast<std::size_t>(i)];
          point.x += lx * ix;
          point.y += ly * iy;
          point.z += lz * iz;
          point.atomID = at + 1;
          out.idIndexMap[point.atomID] = at;
          out.pts[static_cast<std::size_t>(at)] = point;
          ++at;
        }
      }
    }
  }
  return out;
}

int signatureCount(const Cloud &cloud, const char *spec) {
  const auto signature = cage::Signature::parse(spec);
  const int depth = std::max(signature.maxRingSize(), 3);
  auto listed = nneigh::neighListO(3.5, cloud, 1);
  auto idx = nneigh::neighbourListByIndex(cloud, listed);
  const auto rings =
      ring::cageRings(primitive::ringNetwork(idx, depth), cloud, 3.5);
  const auto found = cage::findBySignatures(
      rings, idx, std::vector<cage::Signature>{signature}, {});
  return static_cast<int>(found.at(0).size());
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 5) {
    std::cerr << "usage: bench_trace_threads TRAJ tile warmup repeats\n";
    return 2;
  }
  const std::string traj = argv[1];
  const int ntile = std::atoi(argv[2]);
  const int warmup = std::atoi(argv[3]);
  const int repeats = std::atoi(argv[4]);
  if (ntile < 1 || warmup < 0 || repeats < 1) {
    std::cerr << "bad tile, warmup, or repeats\n";
    return 2;
  }
  Cloud raw;
  raw = sinp::readLammpsTrjO(traj, 1, raw, 1);
  if (raw.nop <= 0) {
    std::cerr << "empty " << traj << "\n";
    return 1;
  }
  const Cloud cloud = tile(std::move(raw), ntile);
  const int n512 = signatureCount(cloud, "512");
  const int n62 = signatureCount(cloud, "51262");
  const int n64 = signatureCount(cloud, "51264");
  std::cout << "atoms " << cloud.nop << "\n";
  std::cout << "census " << n512 << " " << n62 << " " << n64 << "\n";
  std::cout.flush();
  for (int rep = 0; rep < warmup + repeats; ++rep) {
    const auto t0 = std::chrono::steady_clock::now();
    const int a = signatureCount(cloud, "512");
    const int b = signatureCount(cloud, "51262");
    const int c = signatureCount(cloud, "51264");
    const double ms = msSince(t0);
    if (a != n512 || b != n62 || c != n64) {
      std::cerr << "census changed during timing\n";
      return 1;
    }
    if (rep >= warmup) {
      std::cout << "ms " << ms << "\n";
      std::cout.flush();
    }
  }
  return 0;
}
