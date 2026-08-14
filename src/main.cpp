//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <opt_parser.h>
#include <seams_yaml.hpp>

#include <lua_api.hpp>
#include <mol_sys.hpp>

#include <rang.hpp>
#include <sol/sol.hpp>

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

//! Where the vendored Fennel compiler lives. The build embeds the source-root
//! path at configure time; YODA_FENNEL_PATH overrides it at runtime (e.g. for
//! installed binaries running away from a checkout).
std::string fennelSourcePath() {
  const char *env = std::getenv("YODA_FENNEL_PATH");
  if (env != nullptr && *env != '\0') {
    return env;
  }
#ifdef YODA_FENNEL_LUA
  return YODA_FENNEL_LUA;
#else
  return "src/include/external/fennel/fennel.lua";
#endif
}

//! Runs a user script in protected mode. Paths ending in .fnl are compiled by
//! the vendored Fennel (a lisp that compiles to Lua) against the same
//! globals; anything else runs as plain Lua. Errors in either language are
//! reported on stderr and returned as a non-zero exit code instead of
//! aborting the process
[[nodiscard]] int runScriptFile(sol::state &lua, const std::string &path) {
  const std::string fnlExt = ".fnl";
  const bool isFennel =
      path.size() >= fnlExt.size() &&
      path.compare(path.size() - fnlExt.size(), fnlExt.size(), fnlExt) == 0;
  if (!isFennel) {
    const sol::protected_function_result result =
        lua.safe_script_file(path, sol::script_pass_on_error);
    if (!result.valid()) {
      const sol::error err = result;
      std::cerr << rang::fg::red << "Lua error in " << path << ":\n"
                << rang::fg::reset << err.what() << "\n";
      return 1;
    }
    return 0;
  }
  sol::table loaded = lua["package"]["loaded"];
  sol::table fennel;
  if (sol::object cached = loaded["fennel"]; cached.is<sol::table>()) {
    fennel = cached.as<sol::table>();
  } else {
    const std::string fennelLua = fennelSourcePath();
    auto mod = lua.safe_script_file(fennelLua, sol::script_pass_on_error);
    if (!mod.valid()) {
      const sol::error err = mod;
      std::cerr << rang::fg::red << "Cannot load the Fennel compiler from "
                << fennelLua << " (set YODA_FENNEL_PATH to override):\n"
                << rang::fg::reset << err.what() << "\n";
      return 1;
    }
    fennel = mod.get<sol::table>();
    loaded["fennel"] = fennel;
  }
  sol::protected_function dofile = fennel["dofile"];
  sol::table opts = lua.create_table();
  opts["allowedGlobals"] = false;
  auto run = dofile(path, opts);
  if (!run.valid()) {
    const sol::error err = run;
    std::cerr << rang::fg::red << "Fennel error in " << path << ":\n"
              << rang::fg::reset << err.what() << "\n";
    return 1;
  }
  return 0;
}

} // namespace

int main(int argc, char *argv[]) {
  auto result = parse(argc, argv);
  const auto cfg = seams::loadLuaConfig(result["c"].as<std::string>());
  const std::string tFile = cfg.trajectory;

  sol::state lua;
  lua.open_libraries();
  luaApi::registerAll(lua);

  if (cfg.topoTwoDim) {
    if (runScriptFile(lua, cfg.variables) != 0) {
      return 1;
    }
    molSys::PointCloud<molSys::Point<double>, double> resCloud;
    std::vector<std::vector<int>> nList, hbnList;
    std::vector<std::vector<int>> rings;
    std::vector<double> rdfValues;
    auto lscript = lua.get<std::string>("functionScript");
    lua["doBOP"] = cfg.bulkUse;
    lua["topoOneDim"] = cfg.topoOneDim;
    lua["topoTwoDim"] = cfg.topoTwoDim;
    lua["topoBulk"] = cfg.bulkUse;
    lua["nList"] = &nList;
    lua["hbnList"] = &hbnList;
    lua["resCloud"] = &resCloud;
    lua["trajectory"] = tFile;
    lua["ringsAllSizes"] = &rings;
    lua["rdf"] = &rdfValues;
    if (runScriptFile(lua, lscript) != 0) {
      return 1;
    }
  }

  if (cfg.topoOneDim) {
    if (runScriptFile(lua, cfg.variables) != 0) {
      return 1;
    }
    molSys::PointCloud<molSys::Point<double>, double> resCloud, oCloud, hCloud;
    std::vector<std::vector<int>> nList, hbnList, rings;
    int atomID = 0;
    auto lscript = lua.get<std::string>("functionScript");
    lua["doBOP"] = cfg.bulkUse;
    lua["topoOneDim"] = cfg.topoOneDim;
    lua["topoTwoDim"] = cfg.topoTwoDim;
    lua["topoBulk"] = cfg.bulkUse;
    lua["nList"] = &nList;
    lua["hbnList"] = &hbnList;
    lua["resCloud"] = &resCloud;
    lua["oCloud"] = &oCloud;
    lua["hCloud"] = &hCloud;
    lua["trajectory"] = tFile;
    lua["ringsAllSizes"] = &rings;
    lua["lowestAtomID"] = &atomID;
    if (runScriptFile(lua, lscript) != 0) {
      return 1;
    }
  }

  if (cfg.bulkUse) {
    if (runScriptFile(lua, cfg.variables) != 0) {
      return 1;
    }
    molSys::PointCloud<molSys::Point<double>, double> resCloud, solCloud, oCloud,
        hCloud;
    std::vector<std::vector<int>> nList, hbnList, iceList, rings;
    std::vector<double> avgQ6;
    auto lscript = lua.get<std::string>("functionScript");
    lua["doBOP"] = cfg.bulkBop;
    lua["topoOneDim"] = cfg.topoOneDim;
    lua["topoTwoDim"] = cfg.topoTwoDim;
    lua["topoBulk"] = cfg.bulkTopo;
    lua["nList"] = &nList;
    lua["hbnList"] = &hbnList;
    lua["iceNeighbourList"] = &iceList;
    lua["resCloud"] = &resCloud;
    lua["oCloud"] = &oCloud;
    lua["hCloud"] = &hCloud;
    lua["clusterCloud"] = &solCloud;
    lua["avgQ6"] = &avgQ6;
    lua["trajectory"] = tFile;
    lua["ringsAllSizes"] = &rings;
    if (runScriptFile(lua, lscript) != 0) {
      return 1;
    }
  }

  if (cfg.structureDesc) {
    if (runScriptFile(lua, cfg.variables) != 0) {
      return 1;
    }
    molSys::PointCloud<molSys::Point<double>, double> resCloud;
    std::vector<std::vector<int>> nList;
    auto lscript = lua.get<std::string>("functionScript");
    lua["resCloud"] = &resCloud;
    lua["nList"] = &nList;
    lua["trajectory"] = tFile;
    if (runScriptFile(lua, lscript) != 0) {
      return 1;
    }
  }

  std::cout << rang::style::bold << "Welcome to the Black Parade.\nYou ran:-\n"
            << rang::style::reset
            << "\nBulk Ice Analysis: " << (cfg.bulkUse ? "true" : "false")
            << "\nQuasi-one-dimensional Ice Analysis: "
            << (cfg.topoOneDim ? "true" : "false")
            << "\nQuasi-two-dimensional Ice Analysis: "
            << (cfg.topoTwoDim ? "true" : "false")
            << "\nStructure descriptors: "
            << (cfg.structureDesc ? "true" : "false") << "\n";
  return 0;
}
