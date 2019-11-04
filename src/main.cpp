///////////////////////////////////////////////////////////////////////////////////////////
//
//
// d-SEAMS molecular dynamics analysis engine code
// Copyright (C) <2018-present> Amrita Goswami, Rohit Goswami
// amrita16thaug646[at]gmail.com, r95g10[at]gmail.com

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////////////////

// Standard Library
#include <array>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <string>

// Internal Libraries
#include "opt_parser.h"

// Newer pointCloud
#include <bond.hpp>
#include <bop.hpp>
#include <cluster.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>
#include <topo_bulk.hpp>
#include <topo_one_dim.hpp>

// Externally bundled-input libraries
// #include <cxxopts.hpp>
#include <sol.hpp>

// Managed with Conan
#include <fmt/core.h>
#include <yaml-cpp/yaml.h>
#include <rang.hpp>

int main(int argc, char *argv[]) {
  // Parse Things
  auto result = parse(argc, argv);
  auto &arguments = result.arguments();
  // Initialize yaml config
  YAML::Node config = YAML::LoadFile(result["c"].as<std::string>());
  // This is a dummy used to figure out the order of options (cmd > yml)
  std::string script, tFile;
  // Initialize Lua
  sol::state lua;
  // Use all libraries
  lua.open_libraries();

  // --------------------------------------
  // iceType Determination Block
  if (config["iceType"]["use"].as<bool>()) {
    // Do error handling if the traj is not there
    // Get the damn file
    // Determine script location [cmd > iceType Traj > Traj]
    if (config["iceType"]["trajectory"]) {
      tFile = config["iceType"]["trajectory"].as<std::string>();
    } else {
      tFile = config["trajectory"].as<std::string>();
    }
    // Get variables
    std::string vars = config["iceType"]["variables"].as<std::string>();

    // Use the script
    lua.script_file(vars);

    // Get Variables (explicitly)
    // -----------------
    // Bulk
    auto rc = lua.get<double>("cutoffRadius");
    auto oType = lua.get<int>("oxygenAtomType");
    auto tFrame = lua.get<int>("targetFrame");
    auto fFrame = lua.get<int>("finalFrame");
    auto fGap = lua.get<int>("frameGap");
    auto dumpName = lua.get<std::string>("dumpName");
    auto advUse = lua.get<bool>("defineFunctions");
    auto outFileChillPlus = lua.get<std::string>("chillPlus_noMod");
    auto outFileChill = lua.get<std::string>("chill_mod");
    auto outFileSuper = lua.get<std::string>("chillPlus_mod");
    auto outCluster = lua.get<std::string>("largest_ice_cluster_name");
    // -----------------
    // Topological Network Ring lua variables
    auto hType =
        lua.get<int>("hydrogenAtomType");  // If you want to use the hydrogen
                                           // atoms to get the HBN
    auto maxDepth = lua.get<int>("maxDepth");  // If you want to use the
                                               // hydrogen atoms to get the HBN
    // -----------------
    // Variables for output directory writing: TODO: shift to config??
    auto doBOP = lua.get<bool>("doBOP");  // If you want to do BOP analysis
    auto topoOneDim = lua.get<bool>(
        "topoOneDim");  // If you want to do the topological analysis for INTs
    auto topoTwoDim = lua.get<bool>(
        "topoTwoDim");  // If you want to do topological analysis for monolayers
    auto topoBulk = lua.get<bool>(
        "topoBulk");  // If you want topological network analysis for bulk
    auto printCages = lua.get<bool>(
        "printCages");  // If you want topological network analysis for bulk
    // -----------------
    // Bulk/Common Variables defined in C++
    // Variables which must be declared in C++
    //
    // Newer pointCloud (rescloud -> ice structure, solcloud -> largest cluster)
    molSys::PointCloud<molSys::Point<double>, double> resCloud, solCloud;
    // Some neighbor
    std::vector<std::vector<int>> nList, hbnList;
    // For averaged q6
    std::vector<double> avgQ6;
    // For the list of all rings (of all sizes)
    std::vector<std::vector<int>> ringsAllSizes;
    std::vector<std::vector<int>> rings;
    // -----------------
    // Variables defined in C++ specific to confined systems

    if (advUse == true) {
      // This section basically only registers functions and handles the rest in
      // lua Use the functions defined here
      auto lscript = lua.get<std::string>("functionScript");
      // Transfer variables to lua
      lua["nList"] = &nList;
      lua["hbnList"] = &hbnList;
      lua["resCloud"] = &resCloud;
      lua["clusterCloud"] = &solCloud;
      lua["avgQ6"] = &avgQ6;
      lua["trajectory"] = tFile;
      // Confined ice stuff
      lua["ringsAllSizes"] = &rings;
      // Register functions
      //
      // Writing stuff
      lua.set_function("writeDump", sout::writeDump);
      lua.set_function("writeHistogram", sout::writeHisto);
      // Generic requirements
      lua.set_function("readFrame", sinp::readLammpsTrjO);
      lua.set_function("neighborList", nneigh::neighListO);
      // CHILL+ and modifications
      lua.set_function("chillPlus_cij", chill::getCorrelPlus);
      lua.set_function("chillPlus_iceType", chill::getIceTypePlus);
      // CHILL functions
      lua.set_function("chill_cij", chill::getCorrel);
      lua.set_function("chill_iceType", chill::getIceType);
      // Reclassify using q6
      lua.set_function("averageQ6", chill::getq6);
      lua.set_function("modifyChill", chill::reclassifyWater);
      lua.set_function("percentage_Ice", chill::printIceType);
      // Largest ice cluster
      lua.set_function("clusterAnalysis", clump::clusterAnalysis);
      // -----------------
      // Topological Network Methods
      // Generic requirements (read in only inside the slice)
      lua.set_function("readFrameOnlyOne", sinp::readLammpsTrjreduced);
      lua.set_function("getHbondNetwork", bond::populateHbonds);
      lua.set_function("hBondNetworkByIndex", nneigh::neighbourListByIndex);
      // -----------------
      // Primitive rings
      lua.set_function("getPrimitiveRings", primitive::ringNetwork);
      // -----------------
      // Quasi-one-dimensional ice
      lua.set_function("prismAnalysis", ring::prismAnalysis);
      // -----------------
      // Bulk ice, using the topological network criterion
      lua.set_function("bulkTopologicalNetworkCriterion",
                       ring::topoBulkAnalysis);
      // --------------------------
      // Use the script
      lua.script_file(lscript);
      // --------------------------

    }  // If adv=true
  }    // end of ice type determination block
  // --------------------------------------

  std::cout << rang::style::bold
            << fmt::format("Welcome to the Black Parade.\nYou ran:-\n")
            << rang::style::reset
            << fmt::format("RDF 3D Analysis: {}",
                           config["rdf3D"]["use"].as<bool>())
            << fmt::format("\nRDF 2D Analysis: {}",
                           config["rdf2D"]["use"].as<bool>())
            << fmt::format("\nPhase Transition Analysis: {}",
                           config["transition"]["use"].as<bool>())
            << fmt::format("\nIce Structure: {}\n",
                           config["iceType"]["use"].as<bool>());
  return 0;
}
