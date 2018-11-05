#ifndef _MOLECULAR_SYSTEM_H
#define _MOLECULAR_SYSTEM_H

#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <vector>
#include "molecule.h"
#include "parameter.h"
#include "rdf3D.h"

/*! \brief Class for information in each frame.
 *         This class creates an object containing an array of CMolecule objects 
 with positions, type etc. 
 *
 The functions in this class read in the lammps trajectory or an xyz file, and 
 save the information in a particular snapshot. This can be passed to other objects for
 post-processing. 

 For processing lammps trajectories, 
 - First read in the entire trajectory to get the number of frames inside the trajectory
 file using readWholeTrj()
 - Now read in information at a particular step number using overloaded function
 readParticleFile(step). The lammps IDs of each atom are stored along with the positions

 For processing an xyz file:
 - Use the readParticleFile() to read in the positions of the particles 
 The type ID is set as 1 by default
 */

//Value used to indicate end of neighbour lists
const int nilvalue = 33333333;
//Pi
const double pi = 3.141592653589793;


class NumberOfParticlesNotDefinedException
{
};

class CMolecularSystem {
    
  public:
    //the main object where all properties of all particles are saved
    CMolecularSystem();
    virtual ~CMolecularSystem();
    //Properties of one single molecule
    CMolecule* molecules;
    CParameter* parameter;
    
    //Init the system
    void initializeMolecules(int);
    void initializeMolecules();
    //and delete the System afterwards
    void deleteMolecules();

    //System can be initialized from a lammps trajectory xyz-File
    void InitializeSystem();
    void readParticleFile();
    // Initialize from the lammps trajectory 
    void readWholeTrj();
    // To read from a lammpstrj file use the overloaded
    // function
    void readParticleFile(int );

    // Reads in the box dimensions from the lammps trajectory
    double getBoxLength(std::string );

};

#endif
