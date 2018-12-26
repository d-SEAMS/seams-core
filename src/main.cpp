///////////////////////////////////////////////////////////////////////////////////////////
//
//
// Structure Factor code
// MIT License

// Copyright (c) 2018 Amrita Goswami, Rohit Goswami
//
// amrita16thaug646[at]gmail.com, r95g10[at]gmail.com
// This code has been written for the purpose of obtaining the 1-D analogue of the Structure
// factor for a confined system, from a lammps trajectory file.
//
//
///////////////////////////////////////////////////////////////////////////////////////////

// Standard Library
#include <array>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <string>

// Internal Libraries
#include "density.h"
#include "molecular_system.h"
#include "molecule.h"
#include "opt_parser.h"
#include "output.h"
#include "parameter.h"
#include "rdf2D.h"
#include "rdf3D.h"
#include "spherical_harmonics.h"
#include "structure_factor.h"
#include "transition.h"

// Newer pointCloud
#include <bop.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <rdf.hpp>

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

  // The program reads the parameter file inside the input folder The main obejct
  //     is created.It hold all the functions and data used in the
  //         analysis.
  CMolecularSystem *m_MolSys = new CMolecularSystem;
  // The parameterfile is read
  if (result["f"].count() > 0) {
    m_MolSys->parameter->readParameter(result["f"].as<std::string>());
  } else {
    m_MolSys->parameter->readParameter(config["file"].as<std::string>());
  }
  // System is initalized, memory allocated, ...
  m_MolSys->InitializeSystem();

  // Total number of steps in the trajectory
  int traj_steps = m_MolSys->parameter->nsteps;
  std::cout << " The total number of steps in the trajectory is " << traj_steps
            << "\n";

  // Get random step info at a frame number
  int frame = 4;   // 400 2000
  int nsteps = 10; // 50 100
  m_MolSys->readParticleFile(frame);

  // // ----------------------------------------
  // // 3D RDF (single step)
  // // Create object for 3D RDF
  // Rdf3D *rdf1 = new Rdf3D;
  // // Testing 3D rdf function
  // rdf1->initRDF3D(*m_MolSys, 0.01);
  // // Get the 3D RDF for one step
  // rdf1->singleRDF3D(*m_MolSys); // default ID=1
  // // Print the RDF
  // rdf1->printRDF3D();

  // --------------------------------------------
  // // Single frame Rdf2D
  // // Create object for 2D RDF
  // Rdf2D *rdf = new Rdf2D;
  //  // Testing 2D rdf function. RDF calculated is incorrect if the wrong volume is set
  // double volume = (8)*m_MolSys->parameter->boxx*m_MolSys->parameter->boxy;
  // rdf->initRDFxy(*m_MolSys, 0.05, volume);
  // // Get the 2D RDF for one step
  // rdf->singleRDFxy(*m_MolSys, 17.85, 0.8, 2, 2);
  // // Print the RDF
  // rdf->printRDF2D();
  // // Free the memory
  // rdf->deleteRDF2D();

  // ----------------------------------------------
  //Rdf3D over multiple frames
  if (config["rdf3D"]["use"].as<bool>()) {
    // Create object for 3D RDF
    Rdf3D *rdf3D = new Rdf3D;
    // Testing 3D rdf function. RDF calculated is incorrect if the wrong volume is set
    double volume = m_MolSys->parameter->boxx * m_MolSys->parameter->boxy *
                    m_MolSys->parameter->boxz;
    rdf3D->initRDF3D(*m_MolSys, 0.01, volume);
    // Loop through steps
    for (int istep = 1; istep <= nsteps; istep++) {
      // Get the coordinates at a particule step
      m_MolSys->readParticleFile(frame + istep);
      // Get the 3D RDF at this step
      rdf3D->accumulateRDF3D(*m_MolSys, 2, 2);
    }

    // Normalizes the RDF (required for multiple steps. This
    // is called automatically in the single step RDF function)
    rdf3D->normalizeRDF3D();
    // Print the RDF
    rdf3D->printRDF3D();
    // Cleanup
    rdf3D->deleteRDF3D();
  }

  // // ----------------------------------------------
  // //Structure Factor from 3D RDF
  // StructureFactor *s_k = new StructureFactor;
  // s_k->initStrucFactor(*rdf3D, m_MolSys->parameter->boxx, m_MolSys->parameter->boxy, m_MolSys->parameter->boxz);
  // // ----------------------------------------------

  // ----------------------------------------------
  // Rdf2D over multiple frames
  if (config["rdf2D"]["use"].as<bool>()) {
    // Create object for 2D RDF
    Rdf2D *rdf = new Rdf2D;
    // Testing 2D rdf function. RDF calculated is incorrect if the wrong volume is set
    // double volume = (8)*m_MolSys->parameter->boxx*m_MolSys->parameter->boxy;
    rdf->initRDFxy(*m_MolSys, 0.03, 15);
    // Loop through steps
    for (int istep = 1; istep <= nsteps; istep++) {
      // Get the coordinates at a particule step
      m_MolSys->readParticleFile(frame + istep);
      // Get the 2D RDF at this step
      rdf->accumulateRDFxy(*m_MolSys, 17.85, 0.8, 2, 2);
    }

    // Normalizes the RDF (required for multiple steps. This
    // is called automatically in the single step RDF function)
    // The width of the layer is the argument
    rdf->normalizeRDF2D(0.8);
    // Print the RDF
    rdf->printRDF2D();
    // Cleanup
    rdf->deleteRDF2D();
  }

  // // ----------------------------------------------
  // //Structure Factor from RDF
  // StructureFactor *s_k = new StructureFactor;
  // s_k->initStrucFactor(*rdf, m_MolSys->parameter->boxx, m_MolSys->parameter->boxy);
  // // ----------------------------------------------

  // --------------------------------------
  // Transition system block
  if (config["transition"]["use"].as<bool>()) {
    trans::TransitionSystem *t_sys = new trans::TransitionSystem;
    // Bind Function to class instance
    lua.set_function("transition_probability",
                     &trans::TransitionSystem::mightTrans, t_sys);
    // Pass variables to lua
    lua["trajectory_file"] = m_MolSys->parameter->trajFile;
    lua["steps_in_trajectory"] = traj_steps;
    lua["number_of_particles"] = m_MolSys->parameter->nop;
    // Determine script location
    if (result["s"].count() > 0) {
      script = result["s"].as<std::string>();
    } else {
      script = config["transition"]["script"].as<std::string>();
    }
    // Bind Harmonics

    // Run the script
    lua.script_file(script);
    // Free the memory.
    t_sys->cleanUp();
  }
  // --------------------------------------

  // Free the memory.
  m_MolSys->deleteMolecules();

  // --------------------------------------
  // iceType Determination Block
  if (config["iceType"]["use"].as<bool>()) {

    // Newer pointCloud
    molSys::PointCloud<molSys::Point<double>, double> resCloud;

    // Do error handling if the traj is not there
    // Get the damn file
    // Determine script location [cmd > iceType Traj > Traj]
    if (config["iceType"]["trajectory"]) {
      tFile = config["iceType"]["trajectory"].as<std::string>();
    } else {
      tFile = config["trajectory"].as<std::string>();
    }
    // Determine script location
    if (result["s"].count() > 0) {
      script = result["s"].as<std::string>();
    } else {
      script = config["iceType"]["script"].as<std::string>();
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

    // Analyze
    // For averaged q6
    std::vector<double> avgQ6;

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
    // For changing the names, change it here
    // TODO: Fix this
    std::string outFileChillPlus = "chillPlus.txt"; // Default
    std::string outFileChill = "chill.txt";
    std::string outFileSuper = "superChill.txt";
    std::ofstream clusterFile;
    for (int frame = tFrame; frame <= fFrame; frame += fGap) {
      // Read in a frame
      resCloud = molSys::readLammpsTrj(tFile, frame, &resCloud, oType);
      // // Sort according to atom ID (OPTIONAL)
      // std::sort(resCloud.pts.begin(), resCloud.pts.end(), gen::compareByAtomID);
      // Update the neighbour lists
      resCloud = nneigh::neighList(rc, &resCloud, oType);

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
      resCloud = chill::getIceTypePlus(&resCloud, false);
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
      chill::printIceType(&resCloud, false);
      // ---------------------

      // ---------------------
      // Get the largest ice cluster
      molSys::PointCloud<molSys::Point<double>, double> solCloud;
      int largestIceCluster;
      solCloud = chill::getIceCloud(&resCloud, &solCloud);
      largestIceCluster = chill::largestIceCluster(&solCloud, rc, true, false);
      // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      // Write out the largest ice cluster to a file
      if (frame == tFrame) {
        clusterFile.open("cluster.txt");
        clusterFile << "Frame NumberInCluster\n";
        clusterFile << solCloud.currentFrame << " " << largestIceCluster
                    << "\n";
        clusterFile.close();
      } else {
        clusterFile.open("cluster.txt");
        clusterFile << solCloud.currentFrame << " " << largestIceCluster
                    << "\n";
        clusterFile.close();
      }
      // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      solCloud = clearPointCloud(&solCloud);
      // ---------------------

      // // Print to file here (cij etc)
      // std::string outName = "chill-"+std::to_string(frame);
      // int val = gen::prettyPrintYoda(&resCloud, outName);

      lua["resCloud"] = &resCloud;
      // Register functions
      lua.set_function("writeDump", gen::writeDump);
      // Use the script
      lua.script_file(script);
      // Print to file here
      std::string dumpName = "wat.lammpstrj";
      gen::writeDump(&resCloud, dumpName);

      // Write out Cij, Q3, Q6 to files
      gen::writeHisto(&resCloud, avgQ6);
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
