#ifndef _PARAMETER_H
#define _PARAMETER_H

#include <iostream>
#include <exception>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const int NUMBER_OF_PARAMETERS = 200;
const double PI = 3.14159265;

class ExceptionBadLineInParameterFile{

};

struct s_rawParameter
{
  std::string name;
  std::string value;
};

class CParameter {
  private:
    // Moved s_rawParameter to public
  public:
    CParameter();
    virtual ~CParameter();
    void readParameter(); 
    // To read input files
    s_rawParameter rawParameter[NUMBER_OF_PARAMETERS];
    //Number of Particles
    int nop;
    //Boxsize x and y
    double boxx, boxy, boxz;
    //XYZFile
    std::string xyzFile;
    double neighbordistance;
    // LAMMPS trajectory file
    std::string trajFile;
    // Total number of steps; starting step and ending step number
    int nsteps;
 };


#endif
