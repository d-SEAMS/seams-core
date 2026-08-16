#ifndef SEAMS_CONFIG_HPP_
#define SEAMS_CONFIG_HPP_

#include <ostream>
#include <string>
#include <string_view>

#include <site.hpp>

/** @file seams_config.hpp
 *  @brief Twelve-factor runtime knobs for the engine and the fronts.
 *
 *  Config is what varies between deploys, not the analysis itself.
 *  Defaults live in code. An optional dotenv file (SEAMS_CONFIG, or
 *  ./seams.env) fills unset variables. The process environment wins
 *  over the file. CLI flags win over the environment. Compile-time
 *  meson options stay out of this table.
 */

namespace site {
bool iceScoreAllowed(Family f);
const char *refuseIceScore(Family f);
const char *familyName(Family f);
Family parseFamily(std::string_view name);
} // namespace site

namespace seams::cfg {

struct Runtime {
  int frame = 1;
  int last = 0;
  int jobs = 1;
  int type = 0;
  double cutoff = 3.5;
  int k = 4;
  std::string graph = "seeded";
  site::Family family = site::Family::waterIce;
  double resident = 0.80;
  double cell = 3.0;
  int tpp = 0;
  int block = 0;
  bool offload = false;
  std::string fennel;
  std::string luaPath;
  std::string source;
};

/** Scan argv for --config PATH / --config=PATH. Empty if absent. */
std::string pathFromArgv(int argc, char **argv);

/** Apply a dotenv file: set KEY=VAL only when KEY is not already in
 *  the environment. Throws if required and the file is missing. */
void applyFile(const std::string &path, bool required);

/** Load file (if any) then read the environment into a Runtime. */
Runtime load(std::string_view explicitFile = {});

/** Write Runtime fields into the environment so libraries that only
 *  call getenv (linkcell, OpenMP offload) see the CLI result. */
void exportEnviron(const Runtime &cfg);

void dump(const Runtime &cfg, std::ostream &os);

} // namespace seams::cfg

#endif
