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
#include "structure_factor.h"
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

    // Total number of steps in the trajectory
    int traj_steps = m_MolSys->parameter->nsteps;
    std::cout<<" The total number of steps in the trajectory is "<< traj_steps << "\n";

    //Get random step info at a frame number
    int frame = 400; // 400 2000
    int nsteps = 100; // 50 100
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
    // Create object for 3D RDF
    Rdf3D *rdf3D = new Rdf3D; 
     // Testing 3D rdf function. RDF calculated is incorrect if the wrong volume is set
    double volume = m_MolSys->parameter->boxx*m_MolSys->parameter->boxy*m_MolSys->parameter->boxz;
    rdf3D->initRDF3D(*m_MolSys, 0.01, volume); 
    // Loop through steps
    for (int istep=1; istep<=50; istep++)
    {
        // Get the coordinates at a particule step
        m_MolSys->readParticleFile(frame+istep);
        // Get the 3D RDF at this step
        rdf3D->accumulateRDF3D(*m_MolSys, 2, 2);
    }

    // Normalizes the RDF (required for multiple steps. This
    // is called automatically in the single step RDF function)
    rdf3D->normalizeRDF3D();
    // Print the RDF 
    rdf3D->printRDF3D();

    // ----------------------------------------------
    //Structure Factor from 3D RDF
    StructureFactor *s_k = new StructureFactor;
    s_k->initStrucFactor(*rdf3D, m_MolSys->parameter->boxx, m_MolSys->parameter->boxy, m_MolSys->parameter->boxz);
    // ----------------------------------------------

    
    
    // // ----------------------------------------------
    // //Rdf2D over multiple frames
    // // Create object for 2D RDF
    // Rdf2D *rdf = new Rdf2D; 
    //  // Testing 2D rdf function. RDF calculated is incorrect if the wrong volume is set
    // double volume = (8)*m_MolSys->parameter->boxx*m_MolSys->parameter->boxy;
    // rdf->initRDFxy(*m_MolSys, 0.03, volume); 
    // // Loop through steps
    // for (int istep=1; istep<=nsteps; istep++)
    // {
    //     // Get the coordinates at a particule step
    //     m_MolSys->readParticleFile(frame+istep);
    //     // Get the 2D RDF at this step
    //     rdf->accumulateRDFxy(*m_MolSys, 17.85, 0.8, 2, 2);
    // }

    // // Normalizes the RDF (required for multiple steps. This
    // // is called automatically in the single step RDF function)
    // // The width of the layer is the argument 
    // rdf->normalizeRDF2D(0.8);
    // // Print the RDF 
    // rdf->printRDF2D();
    // // ----------------------------------------------
    // //Structure Factor from RDF
    // StructureFactor *s_k = new StructureFactor;
    // s_k->initStrucFactor(*rdf, m_MolSys->parameter->boxx, m_MolSys->parameter->boxy);
    // // ----------------------------------------------

    //Free the memory. 
    // rdf3D->deleteRDF3D();
    m_MolSys->deleteMolecules();
    // rdf1->deleteRDF3D();
    // rdf->deleteRDF2D();
    // s_k->deleteStrucFactor();
   
    std::cout << "Welcome to the Black Parade \n";
    return 0;
}
