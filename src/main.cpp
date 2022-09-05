///////////////////////////////////////////////////////////////////////////////////////////
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
#include <bulkTUM.hpp>
#include <cluster.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <rdf2d.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>
#include <topo_bulk.hpp>
#include <topo_one_dim.hpp>
#include <topo_two_dim.hpp>
#include <selection.hpp>

// Externally bundled-input libraries
// #include <cxxopts.hpp>
#include <rang.hpp>
#include <sol.hpp>

#include <fmt/core.h>
#include <yaml-cpp/yaml.h>

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
  // Get the trajectory string
  if (config["trajectory"]) {
    tFile = config["trajectory"].as<std::string>();
  } // end of getting the trajectory
  // Get variable file string
  std::string vars = config["variables"].as<std::string>();
  // --------------------------------------
  // Structure determination block for TWO-DIMENSIONAL ICE
  if (config["topoTwoDim"]["use"].as<bool>()) {
    // Use the variables script
    lua.script_file(vars);
    // -----------------
    // Variables which must be declared in C++
    //
    // Newer pointCloud (rescloud -> ice structure, solcloud -> largest cluster)
    molSys::PointCloud<molSys::Point<double>, double> resCloud;
    // Some neighbor lists
    std::vector<std::vector<int>> nList, hbnList;
    // For the list of all rings (of all sizes)
    std::vector<std::vector<int>> ringsAllSizes;
    std::vector<std::vector<int>> rings;
    // RDF stuff
    std::vector<double> rdfValues; // RDF vector
    // -----------------
    // This section basically only registers functions and handles the rest in
    // lua Use the functions defined here
    auto lscript = lua.get<std::string>("functionScript");
    // Transfer variables to lua
    lua["doBOP"] = config["bulk"]["use"].as<bool>();
    lua["topoOneDim"] = config["topoOneDim"]["use"].as<bool>();
    lua["topoTwoDim"] = config["topoTwoDim"]["use"].as<bool>();
    lua["topoBulk"] = config["bulk"]["use"].as<bool>();
    //
    lua["nList"] = &nList;
    lua["hbnList"] = &hbnList;
    lua["resCloud"] = &resCloud;
    lua["trajectory"] = tFile;
    // Confined ice stuff
    lua["ringsAllSizes"] = &rings;
    // RDF stuff
    lua["rdf"] = &rdfValues;
    // -----------------
    // Register functions
    //
    // Writing stuff
    // Generic requirements
    lua.set_function("readFrameOnlyOne", sinp::readLammpsTrjreduced);
    lua.set_function("readFrameOnlyOneAllAtoms", sinp::readLammpsTrj); // reads in all atoms regardless of type  
    lua.set_function("neighborList", nneigh::neighListO);
    // -----------------
    // Topological Network Method Specific Functions
    // Generic requirements (read in only inside the slice)
    lua.set_function("getHbondNetwork", bond::populateHbonds);
    lua.set_function("bondNetworkByIndex", nneigh::neighbourListByIndex);
    // -----------------
    // Primitive rings
    lua.set_function("getPrimitiveRings", primitive::ringNetwork);
    // -----------------
    // Quasi-two-dimensional ice
    lua.set_function("ringAnalysis", ring::polygonRingAnalysis);
    // --------------------------
    // RDF functions
    lua.set_function("calcRDF", rdf2::rdf2Danalysis_AA);
    // --------------------------
    // Use the script
    lua.script_file(lscript);
    // --------------------------

  } // end of two-dimensional ice block
  // --------------------------------------
  // Structure determination block for ONE-DIMENSIONAL ICE
  if (config["topoOneDim"]["use"].as<bool>()) {
    // Use the script
    lua.script_file(vars);
    // -----------------
    // Variables which must be declared in C++
    //
    // Newer pointCloud (rescloud -> ice structure, solcloud -> largest cluster)
    molSys::PointCloud<molSys::Point<double>, double> resCloud;
    molSys::PointCloud<molSys::Point<double>, double> oCloud; // O atom pointCloud 
    molSys::PointCloud<molSys::Point<double>, double> hCloud; // H atom pointCloud
    // Some neighbor
    std::vector<std::vector<int>> nList, hbnList;
    // For the list of all rings (of all sizes)
    std::vector<std::vector<int>> ringsAllSizes;
    std::vector<std::vector<int>> rings;
    int atomID;
    // -----------------
    // This section basically only registers functions and handles the rest in
    // lua Use the functions defined here
    auto lscript = lua.get<std::string>("functionScript");
    // Transfer variables to lua
    lua["doBOP"] = config["bulk"]["use"].as<bool>();
    lua["topoOneDim"] = config["topoOneDim"]["use"].as<bool>();
    lua["topoTwoDim"] = config["topoTwoDim"]["use"].as<bool>();
    lua["topoBulk"] = config["bulk"]["use"].as<bool>();
    //
    lua["nList"] = &nList;
    lua["hbnList"] = &hbnList;
    lua["resCloud"] = &resCloud;
    lua["oCloud"] = &oCloud;
    lua["hCloud"] = &hCloud;
    lua["trajectory"] = tFile;
    // Confined ice stuff
    lua["ringsAllSizes"] = &rings;
    lua["lowestAtomID"] = &atomID;
    // Register functions
    //
    // Writing stuff
    // Generic requirements
    lua.set_function("readFrameOnlyOne", sinp::readLammpsTrjreduced);
    lua.set_function("readFrameOnlyOneAllAtoms", sinp::readLammpsTrj); // reads in all atoms regardless of type  
    lua.set_function("neighborList", nneigh::neighListO);
    // -----------------
    // Topological Network Method Specific Functions
    // Generic requirements (read in only inside the slice)
    lua.set_function("getHbondNetwork", bond::populateHbonds);
    lua.set_function("bondNetworkByIndex", nneigh::neighbourListByIndex);
    // -----------------
    // Primitive rings
    lua.set_function("getPrimitiveRings", primitive::ringNetwork);
    // -----------------
    // Quasi-one-dimensional ice
    lua.set_function("prismAnalysis", ring::prismAnalysis);
    // --------------------------
    // Use the script
    lua.script_file(lscript);
    // --------------------------

  } // end of one-dimensional ice block
  // --------------------------------------
  // Ice Structure Determination for BULK ICE
  if (config["bulk"]["use"].as<bool>()) {
    // Use the variables script
    lua.script_file(vars);
    // Variables which must be declared in C++
    //
    // Newer pointCloud (rescloud -> ice structure, solcloud -> largest
    molSys::PointCloud<molSys::Point<double>, double> resCloud, solCloud;
    // PointCloud for O atoms and H atoms separately 
    molSys::PointCloud<molSys::Point<double>, double> oCloud, hCloud;
    // Some neighbor
    std::vector<std::vector<int>> nList,
        hbnList; // Neighbour lists (by cutoff and hydrogen-bonded neighbour
                 // lists)
    std::vector<std::vector<int>>
        iceList; // Neighbour list for the largest ice cluster
    // For averaged q6
    std::vector<double> avgQ6;
    // For the list of all rings (of all sizes)
    std::vector<std::vector<int>> ringsAllSizes;
    std::vector<std::vector<int>> rings;
    // -----------------
    // Variables defined in C++ specific to confined systems

    // This section basically only registers functions and handles the rest in
    // lua lua Use the functions defined here
    auto lscript = lua.get<std::string>("functionScript");
    // Transfer variables to lua
    lua["doBOP"] = config["bulk"]["bondOrderParameters"].as<bool>();
    lua["topoOneDim"] = config["topoOneDim"]["use"].as<bool>();
    lua["topoTwoDim"] = config["topoTwoDim"]["use"].as<bool>();
    lua["topoBulk"] = config["bulk"]["topologicalNetworkCriterion"].as<bool>();
    //
    lua["nList"] = &nList;
    lua["hbnList"] = &hbnList;
    lua["iceNeighbourList"] = &iceList;
    lua["resCloud"] = &resCloud;
    lua["oCloud"] = &oCloud;
    lua["hCloud"] = &hCloud;
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
    lua.set_function("recenterCluster", clump::recenterClusterCloud);
    // -----------------
    // Selection Functions
    // Function for getting an output PointCloud of a particular atom type from an existing PointCloud
    lua.set_function("getPointCloudAtomsOfOneAtomType", gen::getPointCloudOneAtomType);
    // -----------------
    // Topological Network Methods
    // Generic requirements (read in only inside the slice)
    lua.set_function("readFrameOnlyOne", sinp::readLammpsTrjreduced);
    lua.set_function("readFrameOnlyOneAllAtoms", sinp::readLammpsTrj); // reads in all atoms regardless of type  
    lua.set_function("getHbondNetwork", bond::populateHbonds);
    lua.set_function("getHbondNetworkFromClouds", bond::populateHbondsWithInputClouds);
    lua.set_function("bondNetworkByIndex", nneigh::neighbourListByIndex);
    // -----------------
    // Primitive rings
    lua.set_function("getPrimitiveRings", primitive::ringNetwork);
    // Function for just getting and writing out the ring numbers
    lua.set_function("bulkRingNumberAnalysis", ring::bulkPolygonRingAnalysis);
    // -----------------
    // Bulk ice, using the topological network criterion
    lua.set_function("bulkTopologicalNetworkCriterion", ring::topoBulkAnalysis);
    // --------------------------
    // Bulk ice, using the TUM (Topological Unit Matching Criterion). No need to
    // use bulkTopologicalNetworkCriterion if you use this function
    lua.set_function("bulkTopoUnitMatching", tum3::topoUnitMatchingBulk);
    // --------------------------
    // Use the script
    lua.script_file(lscript);
    // --------------------------

  } // end of bulk ice structure determination block
  // --------------------------------------

  std::cout << rang::style::bold
            << fmt::format("Welcome to the Black Parade.\nYou ran:-\n")
            << rang::style::reset
            << fmt::format("\nBulk Ice Analysis: {}",
                           config["bulk"]["use"].as<bool>())
            << fmt::format("\nQuasi-one-dimensional Ice Analysis: {}",
                           config["topoOneDim"]["use"].as<bool>())
            << fmt::format("\nQuasi-two-dimensional Ice Analysis: {}",
                           config["topoTwoDim"]["use"].as<bool>())
            << "\n";

  return 0;
}
