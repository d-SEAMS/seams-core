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
#include "rdf3D.h"
#include "rdf2D.h"
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
    m_MolSys->readParticleFile(50);

    // // Create object for 3D RDF 
    // Rdf3D *woot = new Rdf3D;
    // // Testing 3D rdf function
    // woot->initRDF3D(*m_MolSys, 0.05); // 4 
    // // Get the 3D RDF for one step
    // woot->singleRDF3D(*m_MolSys, 2, 1); // default ID=1

    // // Print the RDF 
    // woot->printRDF3D();
    
    // // Free the memory 
    // woot->deleteRDF3D();

    // Create object for 2D RDF
    Rdf2D *rdf = new Rdf2D; 
     // Testing 2D rdf function. RDF calculated is incorrect if the wrong volume is set
    double volume = 8*m_MolSys->parameter->boxx*m_MolSys->parameter->boxy;
    rdf->initRDFxy(*m_MolSys, 0.05, volume); 
    // Get the 2D RDF for one step
    rdf->singleRDFxy(*m_MolSys, 19.5, 22, 2, 1);
    // Print the RDF 
    rdf->printRDF2D();
    // Free the memory 
    rdf->deleteRDF2D();

    //Free the memory.
    m_MolSys->deleteMolecules();
   
    std::cout << "Welcome to the Black Parade \n";
    return 0;
}
