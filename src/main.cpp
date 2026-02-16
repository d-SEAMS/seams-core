///////////////////////////////////////////////////////////////////////////////////////////
//
// d-SEAMS molecular dynamics analysis engine code
// Copyright (C) <2018--present> Amrita Goswami, Rohit Goswami
// amrita16thaug646[at]gmail.com, r95g10[at]gmail.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
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
#include <sol/sol.hpp>

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
    double cutoffRadius = lua["cutoffRadius"];
    int oxygenAtomType = lua["oxygenAtomType"];
    int hydrogenAtomType = lua["hydrogenAtomType"];
    int targetFrame = lua["targetFrame"];
    int finalFrame = lua["finalFrame"];
    int frameGap = lua["frameGap"];
    int maxDepth = lua["maxDepth"];
    bool isSlice = lua["isSlice"];

    sol::table sliceLowerLimitsTable = lua["sliceLowerLimits"];
    std::array<double, 3> sliceLowerLimits = {sliceLowerLimitsTable[1],
                                              sliceLowerLimitsTable[2],
                                              sliceLowerLimitsTable[3]};

    sol::table sliceUpperLimitsTable = lua["sliceUpperLimits"];
    std::array<double, 3> sliceUpperLimits = {sliceUpperLimitsTable[1],
                                              sliceUpperLimitsTable[2],
                                              sliceUpperLimitsTable[3]};

    std::string outDir = lua["outDir"];
    std::string functionScript = lua["functionScript"];

    std::double confiningSheetArea = lua["confiningSheetArea"];
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
    // Script equivalent
        for (int frame = targetFrame; frame <= finalFrame; frame += frameGap) {
      resCloud =
          sinp::readLammpsTrjreduced(tFile, frame, oxygenAtomType, isSlice,
                           sliceLowerLimits, sliceUpperLimits); // Get the frame
      nList = nneigh::neighListO(cutoffRadius, &resCloud,
                           oxygenAtomType); // Calculate the neighborlist by ID
      hbnList =
          bond::populateHbonds(tFile, &resCloud, nList, frame,
                          hydrogenAtomType); // Get the hydrogen-bonded network
                                             // for the current frame
      hbnList = nneigh::neighbourListByIndex(
          &resCloud, hbnList); // Hydrogen-bonded network using indices not IDs
      rings = primitive::ringNetwork(
          hbnList, maxDepth); // Gets every ring (non-primitives included)
      ring::prismAnalysis(
          outDir, rings, hbnList, &resCloud, maxDepth, &atomID, targetFrame,
          frame,
          false); // Does the prism analysis for quasi-one-dimensional ice
    }
    // --------------------------

  } // end of two-dimensional ice block
  // --------------------------------------
  // Structure determination block for ONE-DIMENSIONAL ICE
  if (config["topoOneDim"]["use"].as<bool>()) {
    // Use the script
    lua.script_file(vars);
    double cutoffRadius = lua["cutoffRadius"];
    int oxygenAtomType = lua["oxygenAtomType"];
    int hydrogenAtomType = lua["hydrogenAtomType"];
    int targetFrame = lua["targetFrame"];
    int finalFrame = lua["finalFrame"];
    int frameGap = lua["frameGap"];
    int maxDepth = lua["maxDepth"];
    bool isSlice = lua["isSlice"];

    sol::table sliceLowerLimitsTable = lua["sliceLowerLimits"];
    std::array<double, 3> sliceLowerLimits = {sliceLowerLimitsTable[1],
                                              sliceLowerLimitsTable[2],
                                              sliceLowerLimitsTable[3]};

    sol::table sliceUpperLimitsTable = lua["sliceUpperLimits"];
    std::array<double, 3> sliceUpperLimits = {sliceUpperLimitsTable[1],
                                              sliceUpperLimitsTable[2],
                                              sliceUpperLimitsTable[3]};

    std::string outDir = lua["outDir"];
    std::string functionScript = lua["functionScript"];
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
    // Script equivalent
    for (int frame = targetFrame; frame <= finalFrame; frame += frameGap) {
      resCloud =
          sinp::readLammpsTrjO(tFile, frame, oxygenAtomType, isSlice,
                           sliceLowerLimits, sliceUpperLimits); // Get the frame
      nList = nneigh::neighListO(cutoffRadius, &resCloud,
                           oxygenAtomType); // Calculate the neighborlist by ID
      hbnList =
          bond::populateHbonds(tFile, &resCloud, nList, frame,
                          hydrogenAtomType); // Get the hydrogen-bonded network
                                             // for the current frame
      hbnList = nneigh::neighbourListByIndex(
          &resCloud, hbnList); // Hydrogen-bonded network using indices not IDs
      rings = primitive::ringNetwork(
          hbnList, maxDepth); // Gets every ring (non-primitives included)
      ring::prismAnalysis(
          outDir, rings, hbnList, &resCloud, maxDepth, &atomID, targetFrame,
          frame,
          false); // Does the prism analysis for quasi-one-dimensional ice
    }
    // --------------------------

  } // end of one-dimensional ice block
  // --------------------------------------
  // Ice Structure Determination for BULK ICE
  if (config["bulk"]["use"].as<bool>()) {
    // Use the variables script
    lua.script_file(vars);
    //bulk_use
    double cutoffRadius = lua["cutoffRadius"];
    int oxygenAtomType = lua["oxygenAtomType"];
    int hydrogenAtomType = lua["hydrogenAtomType"];
    int targetFrame = lua["targetFrame"];
    int finalFrame = lua["finalFrame"];
    int frameGap = lua["frameGap"];
    int maxDepth = lua["maxDepth"];
    bool isSlice = lua["isSlice"];

    sol::table sliceLowerLimitsTable = lua["sliceLowerLimits"];
    std::array<double, 3> sliceLowerLimits = {sliceLowerLimitsTable[1],
                                              sliceLowerLimitsTable[2],
                                              sliceLowerLimitsTable[3]};

    sol::table sliceUpperLimitsTable = lua["sliceUpperLimits"];
    std::array<double, 3> sliceUpperLimits = {sliceUpperLimitsTable[1],
                                              sliceUpperLimitsTable[2],
                                              sliceUpperLimitsTable[3]};

    std::string outDir = lua["outDir"];
    std::bool onlyTetrahedral = lua["onlyTetrahedral"];

    //bondOrderParameters
    std::string dumpName =  lua["dumpName"]; //Output file name
    std::string chillPlus_mod=  lua["chillPlus_mod"]; //This the modified file
    std::string chillPlus_noMod=  lua["chillPlus_noMod"]; //This is the standard file name
    std::string chill_noMod=  lua["chill_noMod"]; //This is the standard file name
    std::string largest_ice_cluster_name=  lua["largest_ice_cluster_name"]; //This is the standard file name 
    std::string dumpChillP=  lua["dumpChillP"]; //Output dump file for CHILL+ classification
    std::string dumpSupaaP=  lua["dumpSupaaP"]; // Output dump file for the SUPER CHILL+ classification
    std::string largestClusterDump =  lua["largestClusterDump"]; // Output dump file for the largest ice cluster
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
    lua.set_function("selectInSingleSlice", gen::moleculesInSingleSlice);
    lua.set_function("selectEdgeAtomsInRingsWithinSlice", ring::getEdgeMoleculesInRings);
    lua.set_function("selectAtomsInSliceWithRingEdgeAtoms", ring::printSliceGetEdgeMoleculesInRings);
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
    // Script equivalent_bulk_use
    for (int frame = targetFrame; frame <= finalFrame; frame += frameGap) {
      resCloud =
          sinp::readLammpsTrjreduced(tFile, frame, oxygenAtomType, isSlice,
                           sliceLowerLimits, sliceUpperLimits); // Get the frame
      nList = nneigh::neighListO(cutoffRadius, &resCloud,
                           oxygenAtomType); // Calculate the neighborlist by ID
      hbnList =
          bond::populateHbonds(tFile, &resCloud, nList, frame,
                          hydrogenAtomType); // Get the hydrogen-bonded network
                                             // for the current frame
      hbnList = nneigh::neighbourListByIndex(
          &resCloud, hbnList); // Hydrogen-bonded network using indices not IDs
      rings = primitive::ringNetwork(
          hbnList, maxDepth); // Gets every ring (non-primitives included)
      ring::bulkPolygonRingAnalysis(
          outDir, rings, hbnList, &resCloud, maxDepth,
          targetFrame,
          true); // Does the prism analysis for quasi-one-dimensional ice
    }

    // --------------------------

  if (config["bulk"]["topologicalNetworkCriterion"].as<bool>()) {
    // Script equivalent_bulk_topologicalNetworkCriterion
    for (int frame = targetFrame; frame <= finalFrame; frame += frameGap) {
      resCloud =
          sinp::readLammpsTrjreduced(tFile, frame, &resCloud, oxygenAtomType, isSlice,
                           sliceLowerLimits, sliceUpperLimits); // Get the frame
      nList = nneigh::neighListO(cutoffRadius, &resCloud,
                           oxygenAtomType); // Calculate the neighborlist by ID
      clump::clusterAnalysis(outDir, &solCloud, &resCloud, nList, &iceList,
          cutoffRadius, targetFrame, "q6");
      rings = primitive::ringNetwork(
          iceNeighbourList, maxDepth); // Gets every ring (non-primitives included)
      tum3::topoUnitMatchingBulk(
          outDir, rings, &iceList, &solCloud,
          targetFrame,
          true,
          true); // Finds DDCs and HCs
    }
  }

      // --------------------------
  if (config["bulk"]["bondOrderParameters"].as<bool>()) {
    // Script equivalent_bulk_bondOrderParameters
    for (int frame = targetFrame; frame <= finalFrame; frame += frameGap) {
      resCloud =
          sinp::readLammpsTrjO(tFile, frame, &resCloud, oxygenAtomType, isSlice,
                           sliceLowerLimits, sliceUpperLimits); // Get the frame
      nList = nneigh::neighListO(cutoffRadius, &resCloud,
                           oxygenAtomType); // Calculate the neighborlist by ID
      resCloud = clump::getCorrelPlus(&resCloud, nList, isSlice);
      resCloud = clump::getIceTypePlus(&resCloud, nList, outDir, targetFrame, isSlice, chillPlus_noMod);
      sout::writeDump(&resCloud, outDir, dumpChillP);
      avgQ6 = chill::getq6(&resCloud, nList, isSlice);
      resCloud = chill::reclassifyWater(&resCloud, avgQ6);
      chill::printIceType(&resCloud, outDir, targetFrame, isSlice, chillPlus_mod);
      sout::writeDump(&resCloud, outDir, dumpSupaaP);
      clump::clusterAnalysis(outDir, &solCloud, &resCloud, nList, &iceList, cutoffRadius, targetFrame, "q6");
      clump::recenterClusterCloud(&solCloud, &iceList);
      sout::writeDump(&resCloud, outDir, largestClusterDump);
    }
  }

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
