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

#include <iostream>
#include <string>
#include <vector>

namespace {

//! Runs a Lua file in protected mode; a script error is reported on stderr
//! and returned as a non-zero exit code instead of aborting the process
[[nodiscard]] int runScriptFile(sol::state &lua, const std::string &path) {
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
