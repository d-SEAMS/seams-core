// Copyright (c) d-SEAMS developers
#include "../../tests/fresh_temp.hpp"
#include "seams_lang.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {

int check(bool ok, const char *message) {
  if (!ok) {
    std::cerr << message << "\n";
    return 1;
  }
  return 0;
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
  std::filesystem::remove_all(root);
  return bad == 0 ? 0 : 1;
}
