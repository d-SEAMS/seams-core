#ifndef _PARAMETER_H
#define _PARAMETER_H

#include <iostream>
#include <exception>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

const int NUMBER_OF_PARAMETERS = 200;
const double PI = 3.14159265;

class ExceptionBadLineInParameterFile{

};

struct s_rawParameter
{
  string name;
  string value;
};

class CParameter {
  private:
    s_rawParameter rawParameter[NUMBER_OF_PARAMETERS];
  public:
    CParameter();
    virtual ~CParameter();
    void readParameter();
    void checkParameter();  
    //Number of Particles
    int nop;
    //Boxsize x and y
    double boxx, boxy, boxz;
    //XYZFile
    string xyzFile;
    double neighbordistance;
    // LAMMPS trajectory file
    string trajFile;
    // Total number of steps; starting step and ending step number
    int nsteps;

    // Parameters for histogram 
    double binwidth;
    double max_radius;
 };


#endif
