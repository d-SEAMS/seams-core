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
#include <bop.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

// External bundled libraries
// #include <cxxopts.hpp>
#include <sol.hpp>

// Managed with Conan
#include <fmt/core.h>
#include <rang.hpp>
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

    // Variables which must be declared in C++
    // Newer pointCloud (rescloud -> ice structure, solcloud -> largest cluster)
    MolSys::PointCloud<MolSys::Point<double>, double> resCloud, solCloud;
    // For averaged q6
    std::vector<double> avgQ6;

    if (advUse == true) {
      // This section basically only registers functions and handles the rest in lua
      // Use the functions defined here
      auto lscript = lua.get<std::string>("functionScript");
      // Transfer variables to lua
      lua["resCloud"] = &resCloud;
      lua["clusterCloud"] = &solCloud;
      lua["avgQ6"] = &avgQ6;
      lua["trajectory"] = tFile;
      // Register functions
      // Writing stuff
      lua.set_function("writeDump", gen::writeDump);
      lua.set_function("writeHistogram", gen::writeHisto);
      // Generic requirements
      lua.set_function("readFrame", MolSys::readLammpsTrjO);
      lua.set_function("neighborList", nneigh::neighListO);
      // CHILL+ and modifications
      lua.set_function("chillPlus_cij", chill::getCorrelPlus);
      lua.set_function("chillPlus_iceType", chill::getIceTypePlus);
      lua.set_function("averageQ6", chill::getq6);
      lua.set_function("modifyChill", chill::reclassifyWater);
      lua.set_function("percentage_Ice", chill::printIceType);
      // Largest ice cluster
      lua.set_function("create_cluster", chill::getIceCloud);
      lua.set_function("largest_cluster", chill::largestIceCluster);
      lua.set_function("writeCluster", gen::writeCluster);
      // Use the script
      lua.script_file(lscript);
      std::cout << "\nTest\n";
    } else {
      // Analyze

      // Delete later
      std::ofstream cijFile;
      cijFile.open("cij.txt");
      cijFile << "Cij\n";
      cijFile.close();
      std::ofstream q3File;
      q3File.open("q3.txt");
      q3File << "Q3\n";
      q3File.close();
      std::ofstream q6File;
      q6File.open("q6.txt");
      q6File << "Q6\n";
      q6File.close();

      // For overwriting old files
      // and for printing the first line of output files
      std::ofstream clusterFile;
      for (int frame = tFrame; frame <= fFrame; frame += fGap) {
        // Read in a frame
        resCloud = MolSys::readLammpsTrj(tFile, frame, &resCloud, oType);
        // // Sort according to atom ID (OPTIONAL)
        // std::sort(resCloud.pts.begin(), resCloud.pts.end(), gen::compareByAtomID);
        // Update the neighbour lists
        resCloud = nneigh::neighListO(rc, &resCloud, oType);

        // ------------------------------
        // If you want to use CHILL+
        // Calculate c_ij
        resCloud = chill::getCorrelPlus(&resCloud, false);
        // Print first line to file
        if (frame == tFrame) {
          std::ofstream chill;
          chill.open(outFileChillPlus);
          chill << "Frame Ic Ih Interfacial Clath InterClath Water Total\n";
          chill.close();
        }
        // Classify according to CHILL+
        resCloud = chill::getIceTypePlus(&resCloud, false, outFileChillPlus);
        // ------------------------------

        // Get the averaged q6 per atom
        // Greater than 0.5 means ice
        // Update the neighbour lists
        // resCloud = nneigh::neighList(3.2, &resCloud, oxyType);
        avgQ6 = chill::getq6(&resCloud, false);

        // Reclassify according to averaged q3 and q6
        resCloud = chill::reclassifyWater(&resCloud, &avgQ6);

        // --------------------
        // Print modified parameter
        // Print first line to file
        if (frame == tFrame) {
          std::ofstream chill;
          chill.open(outFileSuper);
          chill << "Frame Ic Ih Interfacial Clath InterClath Water Total\n";
          chill.close();
        }
        // Print out and calculate the number
        // and percentage of the ice types after reclassification
        chill::printIceType(&resCloud, false, outFileSuper);
        // ---------------------

        // ---------------------
        // Get the largest ice cluster
        MolSys::PointCloud<MolSys::Point<double>, double> solCloud;
        int largestIceCluster;
        solCloud = chill::getIceCloud(&resCloud, &solCloud);
        largestIceCluster =
            chill::largestIceCluster(&solCloud, rc, true, false);
        // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        // Write out the largest ice cluster to a file
        if (frame == tFrame) {
          clusterFile.open("cluster.txt");
          clusterFile << "Frame NumberInCluster\n";
          clusterFile << solCloud.currentFrame << " " << largestIceCluster
                      << "\n";
          clusterFile.close();
        } else {
          gen::writeCluster(&solCloud, outCluster, false, largestIceCluster);
        }
        // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        solCloud = clearPointCloud(&solCloud);
        // ---------------------

        // // Print to file here (cij etc)
        // std::string outName = "chill-"+std::to_string(frame);
        // int val = gen::prettyPrintYoda(&resCloud, outName);
        // Print to file here
        gen::writeDump(&resCloud, dumpName);

        // Write out Cij, Q3, Q6 to files
        gen::writeHisto(&resCloud, avgQ6);
      }
    }
  }
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
