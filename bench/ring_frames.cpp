// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// Per-frame seams_rings_find time on a LAMMPS dump, through the public C API
// as the engine adapters call it. --out writes every ring (vertices, size,
// symmetry code), the cage label of every atom and the prism counts, so two
// builds can be compared byte for byte.
//
// usage: seams_ring_frames DUMP [--type T] [--cutoff R] [--reps N] [--out FILE]

#include "seams.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Frame {
  long step = 0;
  double box9[9] = {};
  std::vector<double> xyz;
  std::vector<int> types;
};

bool readFrame(std::istream &in, Frame &f) {
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("ITEM: TIMESTEP", 0) == 0) break;
  }
  if (!in) return false;
  std::getline(in, line);
  f.step = std::atol(line.c_str());
  std::getline(in, line);
  std::getline(in, line);
  const int n = std::atoi(line.c_str());
  std::getline(in, line);
  double lo[3], hi[3];
  for (int d = 0; d < 3; ++d) {
    std::getline(in, line);
    std::istringstream row(line);
    row >> lo[d] >> hi[d];
  }
  std::getline(in, line);
  std::istringstream head(line.substr(std::strlen("ITEM: ATOMS")));
  std::vector<std::string> names;
  for (std::string name; head >> name;) names.push_back(name);
  auto col = [&](const char *name) {
    const auto at = std::find(names.begin(), names.end(), name);
    return at == names.end() ? -1 : static_cast<int>(at - names.begin());
  };
  const int cid = col("id"), ctype = col("type"), cx = col("x"), cy = col("y"), cz = col("z");
  if (cid < 0 || cx < 0 || cy < 0 || cz < 0) return false;
  std::vector<long> ids(n);
  std::vector<double> xyz(3 * static_cast<std::size_t>(n));
  std::vector<int> types(n, 1);
  std::vector<std::string> parts(names.size());
  for (int i = 0; i < n; ++i) {
    std::getline(in, line);
    std::istringstream row(line);
    for (auto &p : parts) row >> p;
    ids[i] = std::atol(parts[cid].c_str());
    if (ctype >= 0) types[i] = std::atoi(parts[ctype].c_str());
    xyz[3 * i] = std::atof(parts[cx].c_str());
    xyz[3 * i + 1] = std::atof(parts[cy].c_str());
    xyz[3 * i + 2] = std::atof(parts[cz].c_str());
  }
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int a, int b) { return ids[a] < ids[b]; });
  f.xyz.resize(xyz.size());
  f.types.resize(n);
  for (int k = 0; k < n; ++k) {
    for (int d = 0; d < 3; ++d) f.xyz[3 * k + d] = xyz[3 * order[k] + d];
    f.types[k] = types[order[k]];
  }
  for (int d = 0; d < 3; ++d) {
    f.box9[d] = lo[d];
    f.box9[3 + d] = hi[d] - lo[d];
  }
  return true;
}

int check(seams_status s, const char *what) {
  if (s == SEAMS_OK) return 0;
  std::fprintf(stderr, "%s: %s: %s\n", what, seams_status_string(s), seams_last_error());
  return 1;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s DUMP [--type T] [--cutoff R] [--reps N] [--out FILE]\n",
                 argv[0]);
    return 2;
  }
  const char *dump = argv[1];
  int type = 0;
  int reps = 1;
  double cutoff = 3.5;
  const char *outPath = nullptr;
  for (int i = 2; i + 1 < argc; i += 2) {
    if (std::strcmp(argv[i], "--type") == 0) {
      type = std::atoi(argv[i + 1]);
    } else if (std::strcmp(argv[i], "--cutoff") == 0) {
      cutoff = std::atof(argv[i + 1]);
    } else if (std::strcmp(argv[i], "--reps") == 0) {
      reps = std::max(1, std::atoi(argv[i + 1]));
    } else if (std::strcmp(argv[i], "--out") == 0) {
      outPath = argv[i + 1];
    } else {
      std::fprintf(stderr, "unknown option %s\n", argv[i]);
      return 2;
    }
  }
  std::ifstream in(dump);
  if (!in) {
    std::fprintf(stderr, "cannot read %s\n", dump);
    return 2;
  }
  std::vector<Frame> frames;
  for (Frame f; readFrame(in, f); f = Frame()) frames.push_back(std::move(f));
  if (frames.empty()) {
    std::fprintf(stderr, "%s has no frames\n", dump);
    return 2;
  }
  std::FILE *out = outPath != nullptr ? std::fopen(outPath, "w") : nullptr;
  if (outPath != nullptr && out == nullptr) {
    std::fprintf(stderr, "cannot write %s\n", outPath);
    return 2;
  }
  double total = 0.0;
  long calls = 0;
  for (const Frame &f : frames) {
    const int n = static_cast<int>(f.types.size());
    seams_system *sys = seams_system_create();
    if (sys == nullptr ||
        check(seams_system_set_frame(sys, f.xyz.data(), n, f.box9, f.types.data(), nullptr),
              "set_frame") != 0) {
      return 1;
    }
    seams_rings *rings = nullptr;
    for (int r = 0; r < reps; ++r) {
      seams_rings_destroy(rings);
      rings = nullptr;
      const auto t0 = std::chrono::steady_clock::now();
      if (check(seams_rings_find(sys, cutoff, type, SEAMS_GRAPH_CUTOFF, 6, &rings),
                "rings_find") != 0) {
        return 1;
      }
      total += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      ++calls;
    }
    if (out != nullptr) {
      const int count = seams_rings_count(rings);
      std::vector<int> sizes(count), symmetry(count);
      if (count > 0 && check(seams_rings_labels(rings, sizes.data(), symmetry.data()),
                             "rings_labels") != 0) {
        return 1;
      }
      std::fprintf(out, "frame %ld rings %d\n", f.step, count);
      std::vector<int> atoms(16);
      for (int k = 0; k < count; ++k) {
        std::size_t got = 0;
        if (check(seams_rings_vertices(rings, k, atoms.data(), atoms.size(), &got),
                  "rings_vertices") != 0) {
          return 1;
        }
        std::fprintf(out, "%d %d", sizes[k], symmetry[k]);
        for (std::size_t v = 0; v < got; ++v) std::fprintf(out, " %d", atoms[v]);
        std::fputc('\n', out);
      }
      seams_cage_set *cages = nullptr;
      if (check(seams_cages_find(sys, cutoff, type, SEAMS_GRAPH_CUTOFF, 0, 0, &cages),
                "cages_find") != 0) {
        return 1;
      }
      int hex = 0, cubic = 0, water = 0, both = 0;
      std::vector<int> per(static_cast<std::size_t>(n));
      if (check(seams_cages_counts(cages, &hex, &cubic, &water, &both), "cages_counts") != 0 ||
          check(seams_cages_labels(cages, per.data(), per.size()), "cages_labels") != 0) {
        return 1;
      }
      seams_cages_destroy(cages);
      std::fprintf(out, "cages hex %d cubic %d water %d both %d\n", hex, cubic, water, both);
      for (int i = 0; i < n; ++i) std::fprintf(out, "%d%c", per[i], i + 1 == n ? '\n' : ' ');
      int perfect[7] = {}, imperfect[7] = {};
      if (check(seams_prisms(sys, cutoff, type, perfect, imperfect, 7), "prisms") != 0) return 1;
      std::fprintf(out, "prisms");
      for (int s = 3; s <= 6; ++s) std::fprintf(out, " %d:%d/%d", s, perfect[s], imperfect[s]);
      std::fputc('\n', out);
    }
    seams_rings_destroy(rings);
    seams_system_destroy(sys);
  }
  if (out != nullptr) std::fclose(out);
  std::printf("frames %zu atoms %zu calls %ld s_per_frame %.6e\n", frames.size(),
              frames.front().types.size(), calls, total / static_cast<double>(calls));
  return 0;
}
