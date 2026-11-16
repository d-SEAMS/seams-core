//-----------------------------------------------------------------------------------
// d-SEAMS
// SPDX-License-Identifier: MIT
//
// Compare the shell classifier with and without open-face second-shell
// padding. CHILL+ is the published rule and is not modified. Ice for both
// scores is cubic or hexagonal: a d-SEAMS label of any non-negative kind
// against a CHILL+ atom that is neither cubic nor hexagonal is a false
// positive.
//-----------------------------------------------------------------------------------

#include <chill_plus_scheme.hpp>
#include <local_symmetry.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sched.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

bool chillIce(molSys::atom_state_type kind) {
  return kind == molSys::atom_state_type::cubic || kind == molSys::atom_state_type::hexagonal;
}

void shiftToOrigin(Cloud &cloud) {
  if (cloud.boxLow.size() < 3) {
    return;
  }
  const double z0[3] = {cloud.boxLow[0], cloud.boxLow[1], cloud.boxLow[2]};
  if (z0[0] == 0.0 && z0[1] == 0.0 && z0[2] == 0.0) {
    return;
  }
  for (auto &point : cloud.pts) {
    point.x -= z0[0];
    point.y -= z0[1];
    point.z -= z0[2];
  }
  cloud.boxLow = {0.0, 0.0, 0.0};
}

struct Score {
  int n = 0;
  int chill = 0;
  int off = 0;
  int on = 0;
  int fpOff = 0;
  int fpOn = 0;
  int changed = 0;
  int agreeOff = 0;
  int agreeOn = 0;
};

Score scoreCloud(Cloud &cloud) {
  Score out;
  out.n = cloud.nop;
  if (cloud.nop <= 0) {
    return out;
  }
  auto rows = nneigh::neighListO(3.5, cloud, 1);
  chill::applyChillPlus(cloud, rows, chill::ChillPlusRule::Published, false);
  std::vector<int> off;
  std::vector<int> on;
  symm::classifyIceCloud(cloud, off, false);
  symm::classifyIceCloud(cloud, on, true);
  for (int i = 0; i < cloud.nop; ++i) {
    const bool chill = chillIce(cloud.pts[static_cast<std::size_t>(i)].iceType);
    const bool a = off[static_cast<std::size_t>(i)] >= 0;
    const bool b = on[static_cast<std::size_t>(i)] >= 0;
    out.chill += chill ? 1 : 0;
    out.off += a ? 1 : 0;
    out.on += b ? 1 : 0;
    out.fpOff += (a && !chill) ? 1 : 0;
    out.fpOn += (b && !chill) ? 1 : 0;
    out.changed += off[static_cast<std::size_t>(i)] != on[static_cast<std::size_t>(i)] ? 1 : 0;
    out.agreeOff += (a == chill) ? 1 : 0;
    out.agreeOn += (b == chill) ? 1 : 0;
  }
  return out;
}

void printScore(const char *record, const std::string &lane, const std::string &temperature, int frame,
                int layer, double zLo, const Score &row) {
  std::cout << record << "," << lane << "," << temperature << "," << frame << "," << layer << ","
            << zLo << "," << row.n << "," << row.chill << "," << row.off << "," << row.on << ","
            << row.fpOff << "," << row.fpOn << "," << row.changed << "," << row.agreeOff << ","
            << row.agreeOn << "\n";
}

Cloud loadDump(const std::string &path, int frame) {
  Cloud cloud;
  cloud = sinp::readLammpsTrjO(path, frame, cloud, 1, false, {0, 0, 0}, {0, 0, 0});
  shiftToOrigin(cloud);
  return cloud;
}

Cloud loadOxygenXyz(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot read xyz");
  }
  std::string countLine;
  std::string comment;
  if (!std::getline(in, countLine) || !std::getline(in, comment)) {
    throw std::runtime_error("xyz header is short");
  }
  Cloud cloud;
  cloud.box = {1, 1, 1};
  cloud.boxLow = {0, 0, 0};
  const auto lattice = comment.find("Lattice=\"");
  if (lattice != std::string::npos) {
    const auto start = lattice + 9;
    const auto end = comment.find('"', start);
    std::istringstream numbers(comment.substr(start, end - start));
    double cell[9] = {};
    for (double &value : cell) {
      if (!(numbers >> value)) {
        throw std::runtime_error("xyz lattice is short");
      }
    }
    cloud.box = {cell[0], cell[4], cell[8], cell[3], cell[6], cell[7]};
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.size() < 2 || line[0] != 'O' || (line[1] != ' ' && line[1] != '\t')) {
      continue;
    }
    std::istringstream row(line);
    std::string species;
    molSys::Point<double> point;
    row >> species >> point.x >> point.y >> point.z;
    if (!row) {
      continue;
    }
    point.type = 1;
    point.atomID = cloud.nop + 1;
    point.molID = point.atomID;
    cloud.idIndexMap[point.atomID] = cloud.nop;
    cloud.pts.push_back(point);
    cloud.nop += 1;
  }
  if (cloud.nop <= 0) {
    throw std::runtime_error("xyz has no oxygen");
  }
  return cloud;
}

void layersOf(Cloud &cloud, const std::string &lane, const std::string &temperature, int frame,
              double bin) {
  const Score whole = scoreCloud(cloud);
  printScore("frame", lane, temperature, frame, -1, 0.0, whole);
  const double span = cloud.box.size() > 2 ? cloud.box[2] : 1.0;
  const int nLayer = std::max(1, static_cast<int>(std::ceil(span / bin)));
  // Layer counts use the whole-frame labels. A slice classified on its own
  // would drop neighbours that sit across the bin edge.
  std::vector<int> off;
  std::vector<int> on;
  std::vector<char> chill(static_cast<std::size_t>(cloud.nop), 0);
  symm::classifyIceCloud(cloud, off, false);
  symm::classifyIceCloud(cloud, on, true);
  for (int i = 0; i < cloud.nop; ++i) {
    chill[static_cast<std::size_t>(i)] =
        chillIce(cloud.pts[static_cast<std::size_t>(i)].iceType) ? 1 : 0;
  }
  std::vector<int> which(static_cast<std::size_t>(cloud.nop), 0);
  for (int i = 0; i < cloud.nop; ++i) {
    int layer = static_cast<int>(std::floor(cloud.pts[static_cast<std::size_t>(i)].z / bin));
    if (layer < 0) {
      layer = 0;
    }
    if (layer >= nLayer) {
      layer = nLayer - 1;
    }
    which[static_cast<std::size_t>(i)] = layer;
  }
  for (int layer = 0; layer < nLayer; ++layer) {
    Score row;
    for (int i = 0; i < cloud.nop; ++i) {
      if (which[static_cast<std::size_t>(i)] != layer) {
        continue;
      }
      const bool c = chill[static_cast<std::size_t>(i)] != 0;
      const bool a = off[static_cast<std::size_t>(i)] >= 0;
      const bool b = on[static_cast<std::size_t>(i)] >= 0;
      row.n += 1;
      row.chill += c ? 1 : 0;
      row.off += a ? 1 : 0;
      row.on += b ? 1 : 0;
      row.fpOff += (a && !c) ? 1 : 0;
      row.fpOn += (b && !c) ? 1 : 0;
      row.changed += off[static_cast<std::size_t>(i)] != on[static_cast<std::size_t>(i)] ? 1 : 0;
      row.agreeOff += (a == c) ? 1 : 0;
      row.agreeOn += (b == c) ? 1 : 0;
    }
    if (row.n == 0) {
      continue;
    }
    printScore("layer", lane, temperature, frame, layer, layer * bin, row);
  }
}

Cloud tile(const Cloud &cloud, int copies) {
  Cloud out;
  const double lx = cloud.box.size() > 0 ? cloud.box[0] : 1.0;
  const double ly = cloud.box.size() > 1 ? cloud.box[1] : lx;
  const double lz = cloud.box.size() > 2 ? cloud.box[2] : lx;
  out.box = {lx * copies, ly * copies, lz * copies};
  out.boxLow = {0, 0, 0};
  out.pts.reserve(cloud.pts.size() * static_cast<std::size_t>(copies * copies * copies));
  for (int ix = 0; ix < copies; ++ix) {
    for (int iy = 0; iy < copies; ++iy) {
      for (int iz = 0; iz < copies; ++iz) {
        for (const auto &point : cloud.pts) {
          molSys::Point<double> kept = point;
          kept.x += ix * lx;
          kept.y += iy * ly;
          kept.z += iz * lz;
          kept.type = 1;
          kept.atomID = out.nop + 1;
          kept.molID = kept.atomID;
          out.idIndexMap[kept.atomID] = out.nop;
          out.pts.push_back(kept);
          out.nop += 1;
        }
      }
    }
  }
  return out;
}

double seconds(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

void timeOne(Cloud &cloud, int reps) {
  std::vector<int> kind;
  symm::classifyIceCloud(cloud, kind, false);
  symm::classifyIceCloud(cloud, kind, true);
  for (int rep = 1; rep <= reps; ++rep) {
    const auto offStart = std::chrono::steady_clock::now();
    symm::classifyIceCloud(cloud, kind, false);
    const double off = seconds(offStart);
    const auto onStart = std::chrono::steady_clock::now();
    symm::classifyIceCloud(cloud, kind, true);
    const double on = seconds(onStart);
    std::cout << "time," << cloud.nop << ",0," << rep << "," << off << "\n";
    std::cout << "time," << cloud.nop << ",1," << rep << "," << on << "\n";
  }
}

void printCores() {
  const long online = sysconf(_SC_NPROCESSORS_ONLN);
  cpu_set_t set;
  CPU_ZERO(&set);
  int affinity = 0;
  if (sched_getaffinity(0, sizeof(set), &set) == 0) {
    affinity = CPU_COUNT(&set);
  }
#ifdef _OPENMP
  const int threads = omp_get_max_threads();
#else
  const int threads = 1;
#endif
  std::cout << "cores," << online << "," << affinity << "," << threads << "\n";
}

} // namespace

int main(int argc, char **argv) {
  std::cout.setf(std::ios::fixed);
  std::cout.precision(6);
  try {
    if (argc < 2) {
      std::cerr << "usage: padding_probe liquid|layers|xyz|time ...\n";
      return 2;
    }
    const std::string mode = argv[1];
    if (mode == "liquid") {
      if (argc < 6) {
        return 2;
      }
      const std::string path = argv[2];
      const std::string temperature = argv[3];
      const int first = std::stoi(argv[4]);
      const int last = std::stoi(argv[5]);
      for (int frame = first; frame <= last; ++frame) {
        Cloud cloud = loadDump(path, frame);
        if (cloud.nop <= 0) {
          std::cerr << "empty frame " << frame << "\n";
          return 2;
        }
        printScore("liquid", "mw-bulk", temperature, frame, -1, 0.0, scoreCloud(cloud));
      }
      return 0;
    }
    if (mode == "layers") {
      if (argc < 9) {
        return 2;
      }
      const std::string path = argv[2];
      const std::string lane = argv[3];
      const std::string temperature = argv[4];
      const int first = std::stoi(argv[5]);
      const int last = std::stoi(argv[6]);
      const int stride = std::stoi(argv[7]);
      const double bin = std::stod(argv[8]);
      for (int frame = first; frame <= last; frame += stride) {
        Cloud cloud = loadDump(path, frame);
        if (cloud.nop <= 0) {
          std::cerr << "empty frame " << frame << "\n";
          return 2;
        }
        layersOf(cloud, lane, temperature, frame, bin);
      }
      return 0;
    }
    if (mode == "xyz") {
      if (argc < 6) {
        return 2;
      }
      Cloud cloud = loadOxygenXyz(argv[2]);
      layersOf(cloud, argv[3], argv[4], 1, std::stod(argv[5]));
      return 0;
    }
    if (mode == "time") {
      if (argc < 5) {
        return 2;
      }
      const int frame = std::stoi(argv[3]);
      const int reps = std::stoi(argv[4]);
      if (reps < 8) {
        std::cerr << "timing asks for at least 8 replicates\n";
        return 2;
      }
      Cloud liquid = loadDump(argv[2], frame);
      if (liquid.nop != 4096) {
        std::cerr << "timing frame has " << liquid.nop << " atoms\n";
        return 2;
      }
      printCores();
      timeOne(liquid, reps);
      Cloud large = tile(liquid, 3);
      if (large.nop != 110592) {
        std::cerr << "tile has " << large.nop << " atoms\n";
        return 2;
      }
      timeOne(large, reps);
      return 0;
    }
  } catch (const std::exception &error) {
    std::cerr << error.what() << "\n";
    return 2;
  }
  std::cerr << "unknown mode\n";
  return 2;
}
