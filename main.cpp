///////////////////////////////////////////////////////////////////////////////////////////
//
//
// Structure Factor code
// MIT License

// Copyright (c) 2018 Amrita Goswami
//
// amrita16thaug646@gmail.com
// This code has been written for the purpose of obtaining the 1-D analogue of the Structure
// factor for a confined system, from a lammps trajectory file. 
//
//
///////////////////////////////////////////////////////////////////////////////////////////

#include "molecular_system.h"
#include "molecule.h"
#include "parameter.h"
#include "analysis.h"
#include "output.h"
#include <ctime>
#include <sstream>
#include <string>
#include <cstdlib>

int main()
{
    // The program reads the parameter file inside the input folder
    // The main obejct is created. It hold all the functions and data
    // used in the analysis.
    CMolecularSystem *m_MolSys = new CMolecularSystem;
    // The parameterfile is read
    m_MolSys->parameter->readParameter();
    // System is initalized, memory allocated, ...
    m_MolSys->InitializeSystem();

    //Get random step
    m_MolSys->readParticleFile(150);

    // Create object for 3D RDF 
    CAnalysis *woot = new CAnalysis;
    // Testing 3D rdf function
    woot->initRDF3D(*m_MolSys);
    // Get the 3D RDF for one step
    woot->singleRDF3D(*m_MolSys); 
    
    // Free the memory 
    woot->deleteRDF3D();

    //Free the memory.
    m_MolSys->deleteMolecules();

   
    cout << "Welcome to the Black Parade \n";
    return 0;
}
