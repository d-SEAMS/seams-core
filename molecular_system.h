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
#include <sstream>
#include "molecule.h"
#include "parameter.h"
#include "analysis.h"

using namespace std;

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

};

#endif
