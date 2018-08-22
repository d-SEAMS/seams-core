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
 
// const int SLIQ = 0;
// const int SFCC = 1;
// const int SBCC = 2;
// const int SHCP = 3;
// const int SUND = -1;
// const int SOBER = 4;


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
    void readParameter(int);
    void checkParameter();
    void readWholeTrj();  
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
 };


#endif
