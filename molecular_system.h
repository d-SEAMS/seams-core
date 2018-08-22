#ifndef _MOLECULAR_SYSTEM_H
#define _MOLECULAR_SYSTEM_H

#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "molecule.h"
#include "parameter.h"

using namespace std;

//Value used to indicate end of neighbour lists
const int nilvalue = 33333333;
//Pi
const double pi = 3.141592653589793;


class NumberOfParticlesNotDefinedException
{
};

class CMolecularSystem {
  private:
    //converst int to string
    void IntToString(int , string& );
    

    //Get all neighboring particles within neighbordistance
    void get_AllNeighbors();

    
    //get distance between two particles with nearest image
    void get_distance(int,int,double&,double&,double&);
    //get distance between two particles with nearest image
    void get_distancePosition(int,double, double,double,double&,double&,double&);
    //get absolute distance between two particles
    double get_absDistance(int,int);
    //Applies periodic boundary conditions on particle i
    void makeperiodic(int);
    
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

    //System can be initialized from a perfect Fcc, Bcc or Hcp
    //structure, reads in a xyz-File
    // void InitializeSystem();
    // void readParticleFile();
    // To read from a lammpstrj file use the overloaded
    // function 
    void readParticleFile(int );

    //**********************************************************
    //Output
    //**********************************************************


};

#endif
