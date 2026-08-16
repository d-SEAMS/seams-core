#include <seams_config.hpp>

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>

namespace seams::cfg {
namespace {

constexpr const char *kConfigEnv = "SEAMS_CONFIG";
constexpr const char *kCwdFile = "seams.env";

bool validKey(std::string_view key) {
  if (key.empty() || (!std::isalpha(static_cast<unsigned char>(key[0])) &&
                      key[0] != '_')) {
    return false;
  }
  for (char c : key) {
    if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_') {
      return false;
    }
  }
  return true;
}

std::string strip(std::string_view in) {
  std::size_t a = 0;
  std::size_t b = in.size();
  while (a < b && std::isspace(static_cast<unsigned char>(in[a]))) {
    ++a;
  }
  while (b > a && std::isspace(static_cast<unsigned char>(in[b - 1]))) {
    --b;
  }
  if (b - a >= 2) {
    const char q = in[a];
    if ((q == '"' || q == '\'') && in[b - 1] == q) {
      ++a;
      --b;
    }
  }
  return std::string(in.substr(a, b - a));
}

const char *env(const char *key) {
  const char *v = std::getenv(key);
  if (v == nullptr || v[0] == '\0') {
    return nullptr;
  }
  return v;
}

int envInt(const char *key, int fallback) {
  const char *v = env(key);
  if (v == nullptr) {
    return fallback;
  }
  return std::atoi(v);
}

double envDouble(const char *key, double fallback) {
  const char *v = env(key);
  if (v == nullptr) {
    return fallback;
  }
  return std::atof(v);
}

bool envBool(const char *key, bool fallback) {
  const char *v = env(key);
  if (v == nullptr) {
    return fallback;
  }
  return !(v[0] == '0' && v[1] == '\0') && v[0] != 'f' && v[0] != 'F' &&
         v[0] != 'n' && v[0] != 'N';
}

std::string envStr(const char *key) {
  const char *v = env(key);
  return v == nullptr ? std::string() : std::string(v);
}

void putEnv(const char *key, const std::string &val, bool overwrite) {
  if (!overwrite && env(key) != nullptr) {
    return;
  }
#if defined(_WIN32)
  _putenv_s(key, val.c_str());
#else
  setenv(key, val.c_str(), overwrite ? 1 : 0);
#endif
}

int firstPositive(int a, int b) {
  if (a > 0) {
    return a;
  }
  return b;
}

} // namespace

} // namespace seams::cfg

namespace site {

const char *familyName(Family f) {
  switch (f) {
  case Family::waterIce:
    return "waterIce";
  case Family::ionicLiquid:
    return "ionicLiquid";
  case Family::moltenSalt:
    return "moltenSalt";
  case Family::des:
    return "des";
  case Family::electrolyte:
    return "electrolyte";
  case Family::confinedIL:
    return "confinedIL";
  case Family::confinedWater:
    return "confinedWater";
  case Family::networkFormer:
    return "networkFormer";
  }
  return "waterIce";
}

Family parseFamily(std::string_view name) {
  if (name == "waterIce" || name == "water-ice" || name == "water_ice") {
    return Family::waterIce;
  }
  if (name == "ionicLiquid" || name == "ionic-liquid" ||
      name == "ionic_liquid" || name == "il") {
    return Family::ionicLiquid;
  }
  if (name == "moltenSalt" || name == "molten-salt" ||
      name == "molten_salt") {
    return Family::moltenSalt;
  }
  if (name == "des") {
    return Family::des;
  }
  if (name == "electrolyte") {
    return Family::electrolyte;
  }
  if (name == "confinedIL" || name == "confined-il" ||
      name == "confined_il") {
    return Family::confinedIL;
  }
  if (name == "confinedWater" || name == "confined-water" ||
      name == "confined_water") {
    return Family::confinedWater;
  }
  if (name == "networkFormer" || name == "network-former" ||
      name == "network_former" || name == "network") {
    return Family::networkFormer;
  }
  throw std::invalid_argument("unknown family '" + std::string(name) +
                              "' (want waterIce, ionicLiquid, moltenSalt, "
                              "des, electrolyte, confinedIL, confinedWater, "
                              "networkFormer)");
}

bool iceScoreAllowed(Family f) { return f == Family::waterIce; }

const char *refuseIceScore(Family f) {
  switch (f) {
  case Family::waterIce:
    return "";
  case Family::ionicLiquid:
    return "CHILL/TUM refused for family ionicLiquid";
  case Family::moltenSalt:
    return "CHILL/TUM refused for family moltenSalt";
  case Family::des:
    return "CHILL/TUM refused for family des";
  case Family::electrolyte:
    return "CHILL/TUM refused for family electrolyte "
           "(name a waterIce subset)";
  case Family::confinedIL:
    return "CHILL/TUM refused for family confinedIL";
  case Family::confinedWater:
    return "CHILL/TUM refused for family confinedWater";
  case Family::networkFormer:
    return "CHILL/TUM refused for family networkFormer";
  }
  return "CHILL/TUM refused";
}

} // namespace site

namespace seams::cfg {

std::string pathFromArgv(int argc, char **argv) {
  for (int i = 1; i < argc; ++i) {
    const std::string_view a(argv[i]);
    if (a == "--config") {
      if (i + 1 < argc) {
        return argv[i + 1];
      }
      return {};
    }
    constexpr std::string_view p = "--config=";
    if (a.size() > p.size() && a.substr(0, p.size()) == p) {
      return std::string(a.substr(p.size()));
    }
  }
  return {};
}

void applyFile(const std::string &path, bool required) {
  std::ifstream in(path);
  if (!in) {
    if (required) {
      throw std::runtime_error("SEAMS_CONFIG file not found: " + path);
    }
    return;
  }
  std::string line;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    auto hash = line.find('#');
    if (hash != std::string::npos) {
      line = line.substr(0, hash);
    }
    line = strip(line);
    if (line.empty()) {
      continue;
    }
    if (line.rfind("export ", 0) == 0) {
      line = strip(line.substr(7));
    }
    const auto eq = line.find('=');
    if (eq == std::string::npos) {
      throw std::runtime_error(path + ":" + std::to_string(lineno) +
                               ": expected KEY=VAL");
    }
    const std::string key = strip(line.substr(0, eq));
    const std::string val = strip(line.substr(eq + 1));
    if (!validKey(key)) {
      throw std::runtime_error(path + ":" + std::to_string(lineno) +
                               ": invalid key");
    }
    putEnv(key.c_str(), val, false);
  }
}

Runtime load(std::string_view explicitFile) {
  std::string chosen(explicitFile);
  if (chosen.empty()) {
    if (const char *fromEnv = env(kConfigEnv)) {
      chosen = fromEnv;
    } else {
      std::ifstream probe(kCwdFile);
      if (probe) {
        chosen = kCwdFile;
      }
    }
  }
  if (!chosen.empty()) {
    applyFile(chosen, true);
  }
  Runtime r;
  r.frame = envInt("SEAMS_FRAME", r.frame);
  r.last = envInt("SEAMS_LAST", r.last);
  r.jobs = envInt("SEAMS_JOBS", r.jobs);
  r.type = envInt("SEAMS_TYPE", r.type);
  r.cutoff = envDouble("SEAMS_CUTOFF", r.cutoff);
  r.k = envInt("SEAMS_K", r.k);
  if (const char *g = env("SEAMS_GRAPH")) {
    r.graph = g;
  }
  if (const char *fam = env("SEAMS_FAMILY")) {
    r.family = site::parseFamily(fam);
  }
  r.resident = envDouble("SEAMS_RESIDENT", r.resident);
  r.cell = envDouble("SEAMS_CELL", r.cell);
  r.tpp = firstPositive(envInt("LINKCELL_TPP", 0), envInt("SEAMS_TPP", 0));
  r.block =
      firstPositive(envInt("LINKCELL_BLOCK", 0), envInt("SEAMS_BLOCK", 0));
  r.offload = envBool("SEAMS_OFFLOAD", r.offload);
  r.fennel = envStr("YODA_FENNEL_PATH");
  r.luaPath = envStr("YODA_LUA_PATH");
  r.source = chosen;
  return r;
}

void exportEnviron(const Runtime &cfg) {
  putEnv("SEAMS_FRAME", std::to_string(cfg.frame), true);
  putEnv("SEAMS_LAST", std::to_string(cfg.last), true);
  putEnv("SEAMS_JOBS", std::to_string(cfg.jobs), true);
  putEnv("SEAMS_TYPE", std::to_string(cfg.type), true);
  putEnv("SEAMS_CUTOFF", std::to_string(cfg.cutoff), true);
  putEnv("SEAMS_K", std::to_string(cfg.k), true);
  putEnv("SEAMS_GRAPH", cfg.graph, true);
  putEnv("SEAMS_FAMILY", site::familyName(cfg.family), true);
  putEnv("SEAMS_RESIDENT", std::to_string(cfg.resident), true);
  putEnv("SEAMS_CELL", std::to_string(cfg.cell), true);
  if (cfg.tpp > 0) {
    putEnv("LINKCELL_TPP", std::to_string(cfg.tpp), true);
    putEnv("SEAMS_TPP", std::to_string(cfg.tpp), true);
  }
  if (cfg.block > 0) {
    putEnv("LINKCELL_BLOCK", std::to_string(cfg.block), true);
    putEnv("SEAMS_BLOCK", std::to_string(cfg.block), true);
  }
  putEnv("SEAMS_OFFLOAD", cfg.offload ? "1" : "0", true);
  if (!cfg.fennel.empty()) {
    putEnv("YODA_FENNEL_PATH", cfg.fennel, true);
  }
  if (!cfg.luaPath.empty()) {
    putEnv("YODA_LUA_PATH", cfg.luaPath, true);
  }
}

void dump(const Runtime &cfg, std::ostream &os) {
  os << "SEAMS_FRAME=" << cfg.frame << "\n";
  os << "SEAMS_LAST=" << cfg.last << "\n";
  os << "SEAMS_JOBS=" << cfg.jobs << "\n";
  os << "SEAMS_TYPE=" << cfg.type << "\n";
  os << "SEAMS_CUTOFF=" << cfg.cutoff << "\n";
  os << "SEAMS_K=" << cfg.k << "\n";
  os << "SEAMS_GRAPH=" << cfg.graph << "\n";
  os << "SEAMS_FAMILY=" << site::familyName(cfg.family) << "\n";
  os << "SEAMS_RESIDENT=" << cfg.resident << "\n";
  os << "SEAMS_CELL=" << cfg.cell << "\n";
  if (cfg.tpp > 0) {
    os << "LINKCELL_TPP=" << cfg.tpp << "\n";
  }
  if (cfg.block > 0) {
    os << "LINKCELL_BLOCK=" << cfg.block << "\n";
  }
  os << "SEAMS_OFFLOAD=" << (cfg.offload ? "1" : "0") << "\n";
  if (!cfg.fennel.empty()) {
    os << "YODA_FENNEL_PATH=" << cfg.fennel << "\n";
  }
  if (!cfg.luaPath.empty()) {
    os << "YODA_LUA_PATH=" << cfg.luaPath << "\n";
  }
  if (!cfg.source.empty()) {
    os << "# source=" << cfg.source << "\n";
  }
}

} // namespace seams::cfg
