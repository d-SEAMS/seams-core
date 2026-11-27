// Copyright (c) d-SEAMS developers
#pragma once

#include <atomic>
#include <filesystem>
#include <string>
#include <system_error>

#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif

namespace dseams_test {

// One directory per call. A fixed name under the shared temporary
// directory is a second process's file, and a mapped dump of that
// file faults when the other process replaces it.
inline std::filesystem::path freshDir() {
  static std::atomic<unsigned> seq{0};
#if defined(_WIN32)
  const auto pid = _getpid();
#else
  const auto pid = ::getpid();
#endif
  auto dir = std::filesystem::temp_directory_path() /
             ("dseams-" + std::to_string(pid) + "-" +
              std::to_string(seq.fetch_add(1, std::memory_order_relaxed)));
  std::filesystem::create_directories(dir);
  return dir;
}

inline std::filesystem::path freshFile(const std::string &name) {
  return freshDir() / name;
}

// Directory paths in these tests are concatenated with a file name.
inline std::string freshDirSlash() {
  auto text = freshDir().string();
  if (text.empty() || (text.back() != '/' && text.back() != '\\')) {
    text.push_back('/');
  }
  return text;
}

// Writers that use ../output follow the process working directory.
// A fresh child makes that directory private to this process.
class PrivateCwd {
public:
  PrivateCwd()
      : root_(freshDir()), previous_(std::filesystem::current_path()) {
    const auto cwd = root_ / "cwd";
    std::filesystem::create_directories(cwd);
    std::filesystem::current_path(cwd);
  }

  ~PrivateCwd() {
    std::error_code ec;
    std::filesystem::current_path(previous_, ec);
    std::filesystem::remove_all(root_, ec);
  }

  PrivateCwd(const PrivateCwd &) = delete;
  PrivateCwd &operator=(const PrivateCwd &) = delete;

private:
  std::filesystem::path root_;
  std::filesystem::path previous_;
};

} // namespace dseams_test
