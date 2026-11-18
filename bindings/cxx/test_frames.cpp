// Copyright (c) d-SEAMS developers
#include "../../tests/fresh_temp.hpp"
#include "seams_lang.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
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
  seams_frame_free(typed);
  seams_frame_free(picked);
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
  std::filesystem::remove_all(root);
  return bad == 0 ? 0 : 1;
}
