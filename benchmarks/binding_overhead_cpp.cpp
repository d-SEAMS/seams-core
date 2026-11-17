/* Copyright (c) d-SEAMS developers
 * Times seams::System, the header wrapper over the C API.
 * Same cell, cutoff, and orders as the C driver.
 */
#include "seams.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

namespace {

double sink = 0.0;

double now_sec() {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

void retune(int threads) {
  if (threads < 1) {
    throw std::runtime_error("thread count");
  }
  std::string text = std::to_string(threads);
  if (setenv("OMP_NUM_THREADS", text.c_str(), 1) != 0 ||
      setenv("PLUMED_NUM_THREADS", text.c_str(), 1) != 0) {
    throw std::runtime_error("setenv");
  }
  if (seams_set_num_threads(threads) != SEAMS_OK) {
    throw std::runtime_error(seams_last_error());
  }
}

struct Held {
  seams::System sys;
  bool loaded = false;
  void ensure(const std::string &path) {
    if (!loaded) {
      sys.read(path, 1, -1);
      loaded = true;
    }
  }
};

double run(Held &held, const std::string &op, const std::string &path) {
  if (op == "read") {
    double t0 = now_sec();
    held.sys.read(path, 1, -1);
    double dt = now_sec() - t0;
    held.loaded = true;
    sink += held.sys.count();
    return dt;
  }
  held.ensure(path);
  if (op == "knn") {
    double t0 = now_sec();
    {
      auto nb = held.sys.knn(4, 3.5, true);
      sink += static_cast<double>(nb.nnz());
    }
    return now_sec() - t0;
  }
  if (op == "steinhardt") {
    double t0 = now_sec();
    {
      auto got = held.sys.steinhardt(3.5, 6, true);
      sink += got.first.front() + got.second.front();
    }
    return now_sec() - t0;
  }
  if (op == "chill") {
    double t0 = now_sec();
    {
      auto got = held.sys.chill(3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED);
      sink += got.second[SEAMS_LABEL_HEXAGONAL];
    }
    return now_sec() - t0;
  }
  if (op == "rings") {
    double t0 = now_sec();
    {
      auto rings = held.sys.rings(3.5, 6);
      (void)rings;
    }
    return now_sec() - t0;
  }
  throw std::runtime_error(op);
}

}  // namespace

int main() {
  std::ios::sync_with_stdio(false);
  std::cout.setf(std::ios::unitbuf);
  std::unordered_map<std::string, Held> frames;
  std::string line;
  try {
    while (std::getline(std::cin, line)) {
      if (line == "quit" || line.rfind("quit", 0) == 0) {
        break;
      }
      std::istringstream in(line);
      int threads = 0;
      std::string op;
      std::string path;
      if (!(in >> threads >> op >> path)) {
        throw std::runtime_error("command");
      }
      retune(threads);
      std::cout << run(frames[path], op, path) * 1e3 << "\n";
    }
  } catch (const std::exception &ex) {
    std::cerr << "error " << ex.what() << "\n";
    return 1;
  }
  std::cerr << "sink " << sink << "\n";
  return 0;
}
