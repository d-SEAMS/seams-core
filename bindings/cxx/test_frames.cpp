// Copyright (c) d-SEAMS developers
#include "../../tests/fresh_temp.hpp"
#include "seams_lang.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

namespace {

int check(bool ok, const char *message) {
  if (!ok) {
    std::cerr << message << "\n";
    return 1;
  }
  return 0;
}

// A diamond lattice of type 1 with a type 2 atom beside each site.
// Selecting type 1 must give the rows a type 1 graph gives, with the
// type 2 atoms on no cycle.
int chordless_selection() {
  const double a = 2.75 * 4.0 / std::sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  std::vector<double> xyz;
  std::vector<int> types;
  std::vector<int> ids;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        for (const auto &shift : fcc) {
          for (double site : {0.0, 0.25}) {
            const double x = (i + shift[0] + site) * a;
            const double y = (j + shift[1] + site) * a;
            const double z = (k + shift[2] + site) * a;
            xyz.insert(xyz.end(), {x, y, z, x + 0.6, y + 0.6, z});
            types.insert(types.end(), {1, 2});
          }
        }
      }
    }
  }
  const int n = static_cast<int>(types.size());
  for (int i = 0; i < n; ++i) {
    ids.push_back(i + 1);
  }
  const double box[9] = {0, 0, 0, 2 * a, 2 * a, 2 * a, 0, 0, 0};
  seams_frame *typed = seams_frame_from_arrays(xyz.data(), n, box, types.data(), ids.data());
  seams_frame *picked = seams_frame_from_arrays(xyz.data(), n, box, types.data(), ids.data());
  int bad = check(typed != nullptr && picked != nullptr, "chordless frames");
  if (bad != 0) {
    return bad;
  }
  bad += check(seams_frame_select(picked, 1, nullptr, 0) == 0, "select type 1");
  int want_n = 0;
  int got_n = 0;
  std::vector<int> want_per(static_cast<std::size_t>(n), -1);
  std::vector<int> got_per(static_cast<std::size_t>(n), -1);
  bad += check(seams_frame_chordless(typed, 3.5, 1, "cutoff", 6, nullptr, 0, &want_n,
                                     want_per.data()) == 0,
               "typed count");
  bad += check(seams_frame_chordless(picked, 3.5, 0, "cutoff", 6, nullptr, 0, &got_n,
                                     got_per.data()) == 0,
               "selected count");
  bad += check(want_n == 64 * 12 / 6 && got_n == want_n, "twelve six-cycles per site");
  std::vector<int> want(static_cast<std::size_t>(want_n) * 6);
  std::vector<int> got(want.size());
  bad += check(seams_frame_chordless(typed, 3.5, 1, nullptr, 6, want.data(),
                                     static_cast<int>(want.size()), &want_n, nullptr) == 0,
               "typed rows");
  bad += check(seams_frame_chordless(picked, 3.5, 0, "", 6, got.data(),
                                     static_cast<int>(got.size()), &got_n, nullptr) == 0,
               "selected rows");
  bad += check(want == got && want_per == got_per, "selection matches the type 1 graph");
  for (int i = 0; i < n; ++i) {
    bad += check(got_per[static_cast<std::size_t>(i)] == (types[static_cast<std::size_t>(i)] == 1 ? 12 : 0),
                 "per-atom counts");
  }
  bad += check(seams_frame_chordless(typed, 3.5, 1, "cutoff", 6, got.data(),
                                     static_cast<int>(got.size()) - 1, &got_n, nullptr) == 1 &&
                   got_n == want_n,
               "short buffer still counts");
  bad += check(seams_frame_chordless(typed, 3.5, 1, "seeded", 6, nullptr, 0, &got_n, nullptr) == 1,
               "seeded is not a chordless graph");
  bad += check(seams_frame_chordless(typed, 3.5, 1, "cutoff", 13, nullptr, 0, &got_n, nullptr) == 1,
               "size 13");
  bad += check(seams_frame_chordless(typed, 1.0, 1, "mutual", 6, nullptr, 0, &got_n, nullptr) == 0 &&
                   got_n == 0,
               "mutual chordless inside cutoff + 1.5");
  int hex_cages = -1;
  int cubic_cages = -1;
  int water_cages = -1;
  bad += check(seams_frame_cages(typed, 1.0, 1, 4, "seeded", 0, nullptr, &hex_cages,
                                 &cubic_cages, &water_cages, nullptr) == 0 &&
                   hex_cages == 0 && cubic_cages == 0 && water_cages == 64,
               "seeded cages inside cutoff + 1.5");
  std::vector<double> far_q6(static_cast<std::size_t>(n), 1.0);
  bad += check(seams_frame_steinhardt_graph(typed, 1.0, 1, 6, "mutual", far_q6.data(), nullptr) == 0 &&
                   far_q6[0] == 0.0,
               "mutual steinhardt inside cutoff + 1.5");
  seams_frame *empty = seams_frame_read("traj/ice1h.lammpstrj", 1, 9);
  if (empty == nullptr) {
    const auto empty_dump = dseams_test::freshDir() / "empty.lammpstrj";
    {
      std::ofstream out(empty_dump);
      out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
             "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
             "ITEM: ATOMS id type x y z\n";
    }
    empty = seams_frame_read(empty_dump.string().c_str(), 1, 0);
    std::filesystem::remove_all(empty_dump.parent_path());
  }
  double empty_q4 = 0.0;
  double empty_q6 = 0.0;
  int empty_counts[9] = {9, 9, 9, 9, 9, 9, 9, 9, 9};
  bad += check(empty != nullptr &&
                   seams_frame_q4_q6_chill(empty, 3.5, 0, &empty_q4, &empty_q6, nullptr, nullptr,
                                           nullptr, empty_counts) == 0 &&
                   empty_counts[0] == 0 && empty_counts[6] == 0,
               "empty q4 q6 chill zeros counts");
  seams_frame_free(empty);
  seams_frame_free(typed);
  seams_frame_free(picked);
  return bad;
}

int rdf2d_elongated_triclinic() {
  const double lx = 1.0e4;
  const double ly = 12.0;
  const double lz = 12.0;
  const double xy = 3.0;
  const double xz = 3.0;
  const double yz = 0.0;
  const int n = 40;
  unsigned s = 12345;
  auto u = [&]() {
    s = s * 1103515245u + 12345u;
    return static_cast<double>(s >> 8) / static_cast<double>(1 << 24);
  };
  std::vector<double> xyz;
  std::vector<int> types(static_cast<std::size_t>(n), 1);
  std::vector<int> ids;
  for (int i = 0; i < n; ++i) {
    const double fa = u();
    const double fb = u();
    const double fc = u();
    xyz.push_back(fa * lx + fb * xy + fc * xz);
    xyz.push_back(fb * ly + fc * yz);
    xyz.push_back(fc * lz);
    ids.push_back(i + 1);
  }
  const double box[9] = {0.0, 0.0, 0.0, lx, ly, lz, xy, xz, yz};
  seams_frame *frame = seams_frame_from_arrays(xyz.data(), n, box, types.data(), ids.data());
  int bad = check(frame != nullptr, "rdf2d elongated frame");
  if (bad != 0) {
    return bad;
  }
  const double cutoff = 0.5 * lx;
  const double bw = 0.5;
  const int nbin = seams_frame_rdf2d_nbin(cutoff, bw);
  bad += check(nbin > 0, "rdf2d nbin > 0");
  std::vector<double> r(static_cast<std::size_t>(nbin), 0.0);
  std::vector<double> g(static_cast<std::size_t>(nbin), 0.0);
  bad += check(seams_frame_rdf2d(frame, cutoff, bw, nullptr, r.data(), g.data()) == 0,
               "seams_frame_rdf2d elongated");
  seams_frame_free(frame);

  auto exactMin = [&](double dx, double dy, double dz) {
    const double kw = std::round(-dz / lz);
    const double zw = dz + kw * lz;
    const double y0w = dy + kw * yz;
    const double jw = std::round(-y0w / ly);
    const double yw = y0w + jw * ly;
    const double x0w = dx + jw * xy + kw * xz;
    const double xw = x0w + std::round(-x0w / lx) * lx;
    double best = xw * xw + yw * yw + zw * zw;
    const double R = std::sqrt(best);
    const long k0 = static_cast<long>(std::floor((-dz - R) / lz)) - 1;
    const long k1 = static_cast<long>(std::ceil((-dz + R) / lz)) + 1;
    for (long k = k0; k <= k1; ++k) {
      const double z = dz + static_cast<double>(k) * lz;
      if (z * z > best) {
        continue;
      }
      const double rem = std::sqrt(best - z * z);
      const double y0 = dy + static_cast<double>(k) * yz;
      const long j0 = static_cast<long>(std::floor((-y0 - rem) / ly)) - 1;
      const long j1 = static_cast<long>(std::ceil((-y0 + rem) / ly)) + 1;
      for (long j = j0; j <= j1; ++j) {
        const double y = y0 + static_cast<double>(j) * ly;
        const double x0 = dx + static_cast<double>(j) * xy + static_cast<double>(k) * xz;
        const double ic = std::round(-x0 / lx);
        for (double i = ic - 1.0; i <= ic + 1.0; i += 1.0) {
          const double x = x0 + i * lx;
          best = std::min(best, x * x + y * y + z * z);
        }
      }
    }
    return std::sqrt(best);
  };

  std::vector<int> hist(static_cast<std::size_t>(nbin), 0);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double dist = exactMin(
          xyz[3 * static_cast<std::size_t>(j)] - xyz[3 * static_cast<std::size_t>(i)],
          xyz[3 * static_cast<std::size_t>(j) + 1] - xyz[3 * static_cast<std::size_t>(i) + 1],
          xyz[3 * static_cast<std::size_t>(j) + 2] - xyz[3 * static_cast<std::size_t>(i) + 2]);
      if (dist <= cutoff) {
        hist[static_cast<std::size_t>(std::min(static_cast<int>(dist / bw), nbin - 1))] += 2;
      }
    }
  }
  const double volume = lx * ly * lz;
  const double density = static_cast<double>(n) / volume;
  int mismatches = 0;
  for (int b = 0; b < nbin; ++b) {
    const double rb = bw * (static_cast<double>(b) + 0.5);
    const double factor = rb > lz ? lz / (2.0 * rb) : 1.0 - rb / (2.0 * lz);
    const double binVol = 4.0 * std::numbers::pi * std::pow(bw, 3.0) *
                          (std::pow(b + 1, 3.0) - std::pow(b, 3.0)) / 3.0;
    const double refG = static_cast<double>(hist[static_cast<std::size_t>(b)]) /
                        (binVol * static_cast<double>(n) * density * factor);
    if (std::abs(g[static_cast<std::size_t>(b)] - refG) > 1e-10 * (1.0 + std::abs(refG))) {
      ++mismatches;
    }
  }
  bad += check(mismatches == 0, "seams_frame_rdf2d matches exact brute force");
  return bad;
}

} // namespace

int main() {
  const auto root = dseams_test::freshDir();
  auto write = [&](const char *name, const std::string &body) {
    const auto path = root / name;
    std::ofstream out(path);
    out << body;
    return path.string();
  };
  const auto noBox = write("nobox.gro", "title\n1\nonly\n");
  const auto badGro = write("bad.gro", "title\nabc\n");
  const auto shortGro = write("short.gro", "title\n2\none\n");
  const auto blank = write("blank.xyz", "\n1\ncomment\nO 0 0 0\n");
  const auto half = write("half.xyz", "2\ncomment\nO 0 0 0\n");
  const auto one = write("one.xyz", "1\ncomment\nO 0 0 0\n");
  const auto pdb = write("one.pdb", "ATOM bare\nEND\n");
  int bad = 0;
  bad += check(seams_path_nframes(badGro.c_str()) == -1, "bad gro");
  bad += check(seams_path_nframes(shortGro.c_str()) == -1, "short gro");
  bad += check(seams_path_nframes("missing-file.gro") == -1, "missing gro");
  bad += check(seams_path_nframes(one.c_str()) == 1, "one xyz");
  bad += check(seams_path_nframes(pdb.c_str()) >= 0, "one pdb");
  const int noBoxFrames = seams_path_nframes(noBox.c_str());
  if (noBoxFrames == 1) {
    bad += check(seams_path_nframes(blank.c_str()) == 0, "blank xyz");
    bad += check(seams_path_nframes(half.c_str()) == 1, "half xyz");
  } else {
    bad += check(noBoxFrames == -1, "gro without a box");
    bad += check(seams_path_nframes(blank.c_str()) == 1, "blank xyz");
    bad += check(seams_path_nframes(half.c_str()) == -1, "half xyz");
  }
  bad += chordless_selection();
  bad += rdf2d_elongated_triclinic();
  std::filesystem::remove_all(root);
  return bad == 0 ? 0 : 1;
}
