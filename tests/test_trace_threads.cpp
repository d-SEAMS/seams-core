// The closed-cage search and the hexagonal-cage assembly must return the
// same cages at one thread and at a full team. Membership is the sorted
// atom ids of each face, not the ring index, so a reordering of the ring
// list cannot hide a different cage.
#include <cage.hpp>
#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
using Face = std::vector<int>;
using CageMem = std::vector<Face>;

Cloud loadFrame(const std::string &path, int frame) {
  Cloud cloud;
  cloud = sinp::readLammpsTrjO(path, frame, cloud, 1);
  REQUIRE(cloud.nop > 0);
  return cloud;
}

Cloud tile(Cloud src, int n) {
  REQUIRE(n >= 1);
  if (n == 1) {
    return src;
  }
  const double lx = src.box.at(0);
  const double ly = src.box.at(1);
  const double lz = src.box.at(2);
  Cloud out;
  out.box = src.box;
  out.box[0] = lx * n;
  out.box[1] = ly * n;
  out.box[2] = lz * n;
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
          out.pts[static_cast<std::size_t>(at)] = std::move(point);
          ++at;
        }
      }
    }
  }
  return out;
}

Face sortedFace(std::vector<int> ring) {
  std::sort(ring.begin(), ring.end());
  return ring;
}

CageMem facesOf(const std::vector<cage::FoundCage> &found,
                const std::vector<std::vector<int>> &rings) {
  CageMem cages;
  for (const auto &cage : found) {
    std::vector<Face> raw;
    for (const int fi : cage.faces) {
      raw.push_back(sortedFace(rings.at(static_cast<std::size_t>(fi))));
    }
    std::sort(raw.begin(), raw.end());
    Face flat;
    for (const auto &face : raw) {
      flat.push_back(static_cast<int>(face.size()));
      flat.insert(flat.end(), face.begin(), face.end());
    }
    cages.push_back(std::move(flat));
  }
  std::sort(cages.begin(), cages.end());
  return cages;
}

CageMem iceFaces(const Cloud &cloud, const std::vector<std::vector<int>> &idx) {
  const auto rings = primitive::ringNetwork(idx, 6);
  std::vector<std::vector<int>> six;
  for (const auto &ring : rings) {
    if (ring.size() == 6) {
      six.push_back(ring);
    }
  }
  const auto index = ring::buildRingSearchIndex(six, cloud);
  std::vector<ring::strucType> ringType(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cages;
  const auto hc = ring::findHC(six, ringType, idx, cages, index);
  (void)ring::findDDC(six, ringType, hc, cages, index);
  CageMem out;
  for (const auto &cage : cages) {
    std::vector<Face> raw;
    for (const int ri : cage.rings) {
      raw.push_back(sortedFace(six.at(static_cast<std::size_t>(ri))));
    }
    std::sort(raw.begin(), raw.end());
    Face flat;
    flat.push_back(cage.type == cage::cageType::HexC ? 0 : 1);
    for (const auto &face : raw) {
      flat.push_back(static_cast<int>(face.size()));
      flat.insert(flat.end(), face.begin(), face.end());
    }
    out.push_back(std::move(flat));
  }
  std::sort(out.begin(), out.end());
  return out;
}

struct Snapshot {
  CageMem pent;
  CageMem hex2;
  CageMem hex4;
  CageMem ice;
  int n512 = 0;
  int n62 = 0;
  int n64 = 0;
  int nDdc = 0;
};

Snapshot take(const Cloud &cloud) {
  auto listed = nneigh::neighListLinkcell(3.5, cloud, 1);
  if (listed.empty()) {
    listed = nneigh::neighListO(3.5, cloud, 1);
  }
  auto idx = nneigh::neighbourListByIndex(cloud, listed);
  const auto rings6 =
      ring::cageRings(primitive::ringNetwork(idx, 6), cloud, 3.5);
  const auto rings5 =
      ring::cageRings(primitive::ringNetwork(idx, 5), cloud, 3.5);
  const auto sigs = std::vector<cage::Signature>{
      cage::Signature::parse("512"),
      cage::Signature::parse("51262"),
      cage::Signature::parse("51264"),
  };
  // 5^12 only needs pentagons. The binding searches it on a depth-5
  // network, and the two hexagonal signatures share the depth-6 network.
  const auto pent = cage::findBySignatures(
      rings5, idx, std::vector<cage::Signature>{sigs[0]}, {});
  const auto hex = cage::findBySignatures(
      rings6, idx, std::vector<cage::Signature>{sigs[1], sigs[2]}, {});
  Snapshot shot;
  shot.pent = facesOf(pent.at(0), rings5);
  shot.hex2 = facesOf(hex.at(0), rings6);
  shot.hex4 = facesOf(hex.at(1), rings6);
  shot.n512 = static_cast<int>(shot.pent.size());
  shot.n62 = static_cast<int>(shot.hex2.size());
  shot.n64 = static_cast<int>(shot.hex4.size());
  shot.ice = iceFaces(cloud, idx);
  for (const auto &cage : shot.ice) {
    if (!cage.empty() && cage[0] == 1) {
      ++shot.nDdc;
    }
  }
  return shot;
}

void setThreads(int n) {
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(n);
#else
  (void)n;
#endif
}

int teamWidth() {
#ifdef SEAMS_HAS_OPENMP
  const unsigned hw = std::thread::hardware_concurrency();
  const int cap = hw == 0 ? 1 : static_cast<int>(hw);
  return std::max(1, std::min(8, cap));
#else
  return 1;
#endif
}

void requireSame(const Snapshot &serial, const Snapshot &other, const std::string &label) {
  INFO(label);
  REQUIRE(other.n512 == serial.n512);
  REQUIRE(other.n62 == serial.n62);
  REQUIRE(other.n64 == serial.n64);
  REQUIRE(other.nDdc == serial.nDdc);
  REQUIRE(other.pent == serial.pent);
  REQUIRE(other.hex2 == serial.hex2);
  REQUIRE(other.hex4 == serial.hex4);
  REQUIRE(other.ice == serial.ice);
}

void expectMatch(const Cloud &cloud, const std::string &label, int want512, int want62,
                 int want64, int wantDdc) {
  setThreads(1);
  const Snapshot serial = take(cloud);
  if (want512 >= 0) {
    REQUIRE(serial.n512 == want512);
  }
  if (want62 >= 0) {
    REQUIRE(serial.n62 == want62);
  }
  if (want64 >= 0) {
    REQUIRE(serial.n64 == want64);
  }
  if (wantDdc >= 0) {
    REQUIRE(serial.nDdc == wantDdc);
  }
  const int wide = teamWidth();
  for (int threads = 2; threads <= wide; threads *= 2) {
    setThreads(threads);
    requireSame(serial, take(cloud), label + " threads " + std::to_string(threads));
  }
  setThreads(1);
}

} // namespace

TEST_CASE("signature cages match the serial census on GenIce and mW",
          "[trace][threads]") {
  const Cloud sI = loadFrame("traj/genice_sI.lammpstrj", 1);
  expectMatch(sI, "sI", -1, -1, -1, -1);
  expectMatch(tile(sI, 2), "sI 2x2x2", 16, 48, 0, 0);

  const Cloud sII = loadFrame("traj/genice_sII.lammpstrj", 1);
  expectMatch(sII, "sII", 16, 0, 8, 0);
  expectMatch(tile(sII, 2), "sII 2x2x2", 128, 0, 64, 0);

  const Cloud sH = loadFrame("traj/genice_sH.lammpstrj", 1);
  setThreads(1);
  const Snapshot serialH = take(sH);
  const int wide = teamWidth();
  for (int threads = 2; threads <= wide; threads *= 2) {
    setThreads(threads);
    requireSame(serialH, take(sH), "sH threads " + std::to_string(threads));
  }

  const Cloud mw = loadFrame("traj/mW_cubic.lammpstrj", 1);
  REQUIRE(mw.nop == 4096);
  expectMatch(mw, "mW", 0, 0, 0, 2048);
  setThreads(1);
}

TEST_CASE("nucleation frames keep cage membership at every width",
          "[trace][threads]") {
  const int frames = sinp::nLammpsFrames("../data/hop/nucleation.lammpstrj");
  REQUIRE(frames >= 148);
  const int wide = teamWidth();
  for (int frame = 1; frame <= 148; ++frame) {
    const Cloud cloud = loadFrame("../data/hop/nucleation.lammpstrj", frame);
    setThreads(1);
    const Snapshot serial = take(cloud);
    for (int threads = 2; threads <= wide; threads *= 2) {
      setThreads(threads);
      requireSame(serial, take(cloud),
                  "nucleation frame " + std::to_string(frame) + " threads " +
                      std::to_string(threads));
    }
  }
  setThreads(1);
}
