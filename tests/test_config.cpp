#include <catch2/catch_test_macros.hpp>

#include <seams_config.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

const char *kKeys[] = {
    "SEAMS_CONFIG",  "SEAMS_FRAME",    "SEAMS_LAST",     "SEAMS_JOBS",
    "SEAMS_TYPE",    "SEAMS_CUTOFF",   "SEAMS_K",        "SEAMS_GRAPH",
    "SEAMS_RESIDENT", "SEAMS_CELL",    "SEAMS_TPP",      "SEAMS_BLOCK",
    "SEAMS_OFFLOAD", "LINKCELL_TPP",   "LINKCELL_BLOCK", "YODA_FENNEL_PATH",
    "YODA_LUA_PATH"};

struct EnvGuard {
  std::vector<std::pair<std::string, std::string>> saved;
  std::vector<std::string> had;

  EnvGuard() {
    for (const char *k : kKeys) {
      if (const char *v = std::getenv(k)) {
        had.emplace_back(k);
        saved.emplace_back(k, v);
      }
      unsetenv(k);
    }
  }

  ~EnvGuard() {
    for (const char *k : kKeys) {
      unsetenv(k);
    }
    for (const auto &[k, v] : saved) {
      setenv(k.c_str(), v.c_str(), 1);
    }
  }
};

std::filesystem::path tmpEnv(const std::string &body) {
  auto p = std::filesystem::temp_directory_path() /
           ("seams-cfg-" + std::to_string(getpid()) + ".env");
  std::ofstream out(p);
  out << body;
  return p;
}

} // namespace

TEST_CASE("defaults when the environment is empty") {
  EnvGuard g;
  const auto cfg = seams::cfg::load();
  REQUIRE(cfg.frame == 1);
  REQUIRE(cfg.last == 0);
  REQUIRE(cfg.jobs == 1);
  REQUIRE(cfg.cutoff == 3.5);
  REQUIRE(cfg.k == 4);
  REQUIRE(cfg.graph == "seeded");
  REQUIRE(cfg.resident == 0.80);
  REQUIRE(cfg.cell == 3.0);
  REQUIRE(cfg.tpp == 0);
  REQUIRE(cfg.block == 0);
  REQUIRE_FALSE(cfg.offload);
  REQUIRE(cfg.source.empty());
}

TEST_CASE("environment fills the runtime table") {
  EnvGuard g;
  setenv("SEAMS_K", "16", 1);
  setenv("SEAMS_CUTOFF", "3.25", 1);
  setenv("SEAMS_RESIDENT", "0.7", 1);
  setenv("LINKCELL_TPP", "4", 1);
  setenv("SEAMS_OFFLOAD", "1", 1);
  const auto cfg = seams::cfg::load();
  REQUIRE(cfg.k == 16);
  REQUIRE(cfg.cutoff == 3.25);
  REQUIRE(cfg.resident == 0.7);
  REQUIRE(cfg.tpp == 4);
  REQUIRE(cfg.offload);
}

TEST_CASE("file sets only unset keys") {
  EnvGuard g;
  const auto path = tmpEnv("SEAMS_K=8\nSEAMS_GRAPH=knn\n# comment\n");
  setenv("SEAMS_K", "16", 1);
  const auto cfg = seams::cfg::load(path.string());
  REQUIRE(cfg.k == 16);
  REQUIRE(cfg.graph == "knn");
  REQUIRE(cfg.source == path.string());
  std::filesystem::remove(path);
}

TEST_CASE("missing SEAMS_CONFIG is an error") {
  EnvGuard g;
  REQUIRE_THROWS_AS(seams::cfg::load("/no/such/seams.env"), std::runtime_error);
}

TEST_CASE("exportEnviron writes CLI values for getenv readers") {
  EnvGuard g;
  auto cfg = seams::cfg::load();
  cfg.tpp = 8;
  cfg.block = 256;
  cfg.cell = 3.0;
  seams::cfg::exportEnviron(cfg);
  REQUIRE(std::string(std::getenv("LINKCELL_TPP")) == "8");
  REQUIRE(std::string(std::getenv("LINKCELL_BLOCK")) == "256");
  REQUIRE(std::atof(std::getenv("SEAMS_CELL")) == 3.0);
}

TEST_CASE("dump lists resolved keys") {
  EnvGuard g;
  setenv("SEAMS_JOBS", "4", 1);
  std::ostringstream os;
  seams::cfg::dump(seams::cfg::load(), os);
  const auto text = os.str();
  REQUIRE(text.find("SEAMS_JOBS=4") != std::string::npos);
  REQUIRE(text.find("SEAMS_GRAPH=seeded") != std::string::npos);
}

TEST_CASE("pathFromArgv reads --config") {
  char a[] = "seams";
  char b[] = "--config";
  char c[] = "/tmp/x.env";
  char *argv[] = {a, b, c};
  REQUIRE(seams::cfg::pathFromArgv(3, argv) == "/tmp/x.env");
  char d[] = "--config=/tmp/y.env";
  char *argv2[] = {a, d};
  REQUIRE(seams::cfg::pathFromArgv(2, argv2) == "/tmp/y.env");
}
