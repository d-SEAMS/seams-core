// Copyright (c) d-SEAMS developers
#pragma once

#include <atomic>
#include <filesystem>
#include <string>

#include <unistd.h>

namespace dseams_test {

// One directory per call. A fixed name under the shared temporary
// directory is a second process's file, and a mapped dump of that
// file faults when the other process replaces it.
inline std::filesystem::path freshDir() {
  static std::atomic<unsigned> seq{0};
  auto dir = std::filesystem::temp_directory_path() /
             ("dseams-" + std::to_string(::getpid()) + "-" +
              std::to_string(seq.fetch_add(1, std::memory_order_relaxed)));
  std::filesystem::create_directories(dir);
  return dir;
}

inline std::filesystem::path freshFile(const std::string &name) {
  return freshDir() / name;
}

} // namespace dseams_test
