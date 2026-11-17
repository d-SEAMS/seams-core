/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
** Copyright (c) d-SEAMS developers
**
** The DDC search used to build a rotation key for every ring and to split
** the conflict graph by pointer jumping. Ideal ice distinguishes rings by
** their atoms, and the split walks each arc a constant number of times.
** seams_ddc_cycle_evals and seams_ddc_match_steps count that work. A build
** without those counters, or one that still does the eager key and the
** jumping split, fails the bounds below. The cage lines are the membership
** of the same frames.
*/

#include <catch2/catch_test_macros.hpp>

#include <cage.hpp>
#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <dlfcn.h>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

double wrap(double value, double side) {
  double wrapped = std::fmod(value, side);
  if (wrapped < 0.0) {
    wrapped += side;
  }
  return wrapped;
}

void addPoint(Cloud &cloud, double x, double y, double z) {
  molSys::Point<double> point;
  point.type = 1;
  point.atomID = static_cast<int>(cloud.pts.size());
  point.x = x;
  point.y = y;
  point.z = z;
  cloud.idIndexMap[point.atomID] = point.atomID;
  cloud.pts.push_back(point);
}

// Cubic diamond, eight sites per conventional cell. Bond length 2.75 angstrom.
Cloud cubicIce(int cells) {
  const double bond = 2.75;
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  const double side = cells * lattice;
  const double basis[8][3] = {
      {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5},
      {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25}, {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75},
  };
  Cloud cloud;
  cloud.box = {side, side, side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  for (int i = 0; i < cells; i++) {
    for (int j = 0; j < cells; j++) {
      for (int k = 0; k < cells; k++) {
        for (const auto &site : basis) {
          addPoint(cloud, wrap((i + site[0]) * lattice, side),
                   wrap((j + site[1]) * lattice, side),
                   wrap((k + site[2]) * lattice, side));
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  return cloud;
}

// Lonsdaleite, eight sites per cell. Same bond length.
Cloud hexagonalIce(int cells) {
  const double bond = 2.75;
  const double lattice = bond / ((3.0 / 8.0) * std::sqrt(8.0 / 3.0));
  const double height = lattice * std::sqrt(8.0 / 3.0);
  const double xLength = lattice;
  const double yLength = lattice * std::sqrt(3.0);
  const double hex[4][3] = {
      {1.0 / 3.0, 2.0 / 3.0, 0.0},
      {2.0 / 3.0, 1.0 / 3.0, 0.5},
      {1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0},
      {2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0},
  };
  double basis[8][3];
  int nBasis = 0;
  for (const auto &site : hex) {
    const double x = (site[0] - 0.5 * site[1]) * lattice;
    const double y = site[1] * std::sqrt(3.0) * lattice / 2.0;
    const double z = site[2] * height;
    for (int shift = 0; shift < 2; shift++) {
      basis[nBasis][0] = std::fmod(x + (shift ? xLength / 2.0 : 0.0), xLength);
      basis[nBasis][1] = std::fmod(y + (shift ? yLength / 2.0 : 0.0), yLength);
      basis[nBasis][2] = z;
      if (basis[nBasis][0] < 0.0) {
        basis[nBasis][0] += xLength;
      }
      if (basis[nBasis][1] < 0.0) {
        basis[nBasis][1] += yLength;
      }
      nBasis++;
    }
  }
  const double boxX = cells * xLength;
  const double boxY = cells * yLength;
  const double boxZ = cells * height;
  Cloud cloud;
  cloud.box = {boxX, boxY, boxZ};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  for (int i = 0; i < cells; i++) {
    for (int j = 0; j < cells; j++) {
      for (int k = 0; k < cells; k++) {
        for (int b = 0; b < nBasis; b++) {
          addPoint(cloud, wrap(i * xLength + basis[b][0], boxX),
                   wrap(j * yLength + basis[b][1], boxY),
                   wrap(k * height + basis[b][2], boxZ));
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  return cloud;
}

struct CageRun {
  int hcCages = 0;
  int ddcCages = 0;
  int hcRings = 0;
  int ddcRings = 0;
  std::uint64_t membership = 0;
  std::uint64_t cycles = 0;
  std::uint64_t steps = 0;
  std::uint64_t copies = 0;
  std::uint64_t digits = 0;
  std::uint64_t buckets = 0;
};

std::uint64_t fnv1a(const std::string &text) {
  std::uint64_t hash = 14695981039346656037ULL;
  for (const char ch : text) {
    hash ^= static_cast<unsigned char>(ch);
    hash *= 1099511628211ULL;
  }
  return hash;
}

using Counter = std::uint64_t (*)();

Counter counter(const char *name) {
  return reinterpret_cast<Counter>(dlsym(RTLD_DEFAULT, name));
}

CageRun classify(Cloud &cloud) {
  auto neighbours = nneigh::neighListO(3.5, cloud, 1);
  auto indexed = nneigh::neighbourListByIndex(cloud, neighbours);
  auto rings = primitive::ringNetwork(indexed, 7);
  std::vector<std::vector<int>> six;
  for (const auto &ring : rings) {
    if (ring.size() == 6) {
      six.push_back(ring);
    }
  }
  const auto index = ring::buildRingSearchIndex(six, cloud);
  std::vector<ring::strucType> ringType(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cages;
  auto listHC = ring::findHC(six, ringType, indexed, cages, index);
  auto listDDC = ring::findDDC(six, ringType, listHC, cages, index);
  std::vector<std::string> lines;
  for (const auto &cage : cages) {
    std::vector<std::vector<int>> faces;
    for (int id : cage.rings) {
      auto face = six[static_cast<std::size_t>(id)];
      std::sort(face.begin(), face.end());
      faces.push_back(std::move(face));
    }
    std::sort(faces.begin(), faces.end());
    std::string line = std::to_string(static_cast<int>(cage.type));
    for (const auto &face : faces) {
      line.push_back('|');
      for (int atom : face) {
        line += std::to_string(atom);
        line.push_back(',');
      }
    }
    lines.push_back(std::move(line));
  }
  std::sort(lines.begin(), lines.end());
  std::string body;
  for (const auto &line : lines) {
    body += line;
    body.push_back('\n');
  }
  CageRun run;
  for (const auto &cage : cages) {
    if (cage.type == cage::cageType::HexC) {
      run.hcCages++;
    } else if (cage.type == cage::cageType::DoubleDiaC) {
      run.ddcCages++;
    }
  }
  run.hcRings = static_cast<int>(listHC.size());
  run.ddcRings = static_cast<int>(listDDC.size());
  run.membership = fnv1a(body);
  const Counter cycles = counter("seams_ddc_cycle_evals");
  const Counter steps = counter("seams_ddc_match_steps");
  const Counter copies = counter("seams_ddc_row_copies");
  const Counter digits = counter("seams_ddc_key_digits");
  const Counter buckets = counter("seams_ddc_bucket_steps");
  REQUIRE(cycles != nullptr);
  REQUIRE(steps != nullptr);
  REQUIRE(copies != nullptr);
  REQUIRE(digits != nullptr);
  REQUIRE(buckets != nullptr);
  run.cycles = cycles();
  run.steps = steps();
  run.copies = copies();
  run.digits = digits();
  run.buckets = buckets();
  return run;
}

struct Expect {
  int molecules;
  int hcCages;
  int ddcCages;
  int hcRings;
  int ddcRings;
  std::uint64_t membership;
};

} // namespace

TEST_CASE("ideal ice keeps its cages and the DDC search stays linear",
          "[ddc_regression]") {
  const Expect cubic[] = {
      {216, 0, 108, 0, 432, 0x3ab0cafb12d6f17cULL},
      {512, 0, 256, 0, 1024, 0x4fba50541829278eULL},
      {1000, 0, 500, 0, 2000, 0xa5a626ba25b61504ULL},
      {1728, 0, 864, 0, 3456, 0x08c1f4933225e407ULL},
      {2744, 0, 1372, 0, 5488, 0x0e60627a0edf4dc9ULL},
      {4096, 0, 2048, 0, 8192, 0x56a32896faee60dcULL},
  };
  const Expect hexagonal[] = {
      {216, 108, 0, 324, 0, 0x7b2611839cb75d58ULL},
      {512, 256, 0, 1024, 0, 0xe62cc2d673c97e47ULL},
      {1000, 500, 0, 2000, 0, 0xb32c75816c7a3efbULL},
      {1728, 864, 0, 3456, 0, 0xdf16d7ade467e4d7ULL},
      {2744, 1372, 0, 5488, 0, 0x5ab2bfb9e78c6409ULL},
      {4096, 2048, 0, 8192, 0, 0x98ea6260ca8de0ebULL},
  };
  for (const Expect &want : cubic) {
    const int cells = static_cast<int>(std::lround(std::cbrt(want.molecules / 8.0)));
    Cloud cloud = cubicIce(cells);
    REQUIRE(cloud.nop == want.molecules);
    const CageRun run = classify(cloud);
    INFO("cubic " << want.molecules << " membership " << std::hex << run.membership);
    REQUIRE(run.hcCages == want.hcCages);
    REQUIRE(run.ddcCages == want.ddcCages);
    REQUIRE(run.hcRings == want.hcRings);
    REQUIRE(run.ddcRings == want.ddcRings);
    REQUIRE(run.membership == want.membership);
    // No rotation key: the atom ids already order the rings.
    REQUIRE(run.cycles == 0);
    // One walk of each arc at degree 4 and at degree 2. The jumping split
    // visits each arc once per doubling step, which is past this bound.
    REQUIRE(run.steps >= static_cast<std::uint64_t>(2 * want.molecules));
    REQUIRE(run.steps <= static_cast<std::uint64_t>(8 * want.molecules));
    // Ice rows hold six neighbours, so the inline eight never spill.
    // Clique ids are four 16-bit digits. Five digits of the same ranks,
    // or a copy of every conflict row, is past these bounds.
    REQUIRE(run.copies == 0);
    REQUIRE(run.digits <= 16 * static_cast<std::uint64_t>(want.molecules));
    REQUIRE(run.buckets <= 8 * static_cast<std::uint64_t>(want.molecules));
  }
  for (const Expect &want : hexagonal) {
    const int cells = static_cast<int>(std::lround(std::cbrt(want.molecules / 8.0)));
    Cloud cloud = hexagonalIce(cells);
    REQUIRE(cloud.nop == want.molecules);
    const CageRun run = classify(cloud);
    INFO("hexagonal " << want.molecules << " membership " << std::hex << run.membership);
    REQUIRE(run.hcCages == want.hcCages);
    REQUIRE(run.ddcCages == want.ddcCages);
    REQUIRE(run.hcRings == want.hcRings);
    REQUIRE(run.ddcRings == want.ddcRings);
    REQUIRE(run.membership == want.membership);
    REQUIRE(run.cycles == 0);
  }
}

TEST_CASE("structure I keeps two 5^12 cages and six 5^12 6^2 cages",
          "[ddc_regression]") {
  Cloud cloud;
  cloud = sinp::readLammpsTrjO("traj/genice_sI.lammpstrj", 1, cloud, 1);
  REQUIRE(cloud.nop == 46);
  auto neighbours = nneigh::neighListO(3.5, cloud, 1);
  neighbours = nneigh::neighbourListByIndex(cloud, neighbours);
  const auto rings = primitive::ringNetwork(neighbours, 6);
  const auto dodeca =
      cage::findBySignature(rings, neighbours, cage::Signature::parse("512"));
  const auto tetra =
      cage::findBySignature(rings, neighbours, cage::Signature::parse("51262"));
  REQUIRE(dodeca.size() == 2);
  REQUIRE(tetra.size() == 6);
}
