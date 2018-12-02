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
  std::string script;
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

  std::cout << rang::style::bold
            << fmt::format("Welcome to the Black Parade.\nYou ran:-\n")
            << rang::style::reset
            << fmt::format("RDF 3D Analysis: {}",
                           config["rdf3D"]["use"].as<bool>())
            << fmt::format("\nRDF 2D Analysis: {}",
                           config["rdf2D"]["use"].as<bool>())
            << fmt::format("\nPhase Transition Analysis: {}\n",
                           config["transition"]["use"].as<bool>());
  return 0;
}
