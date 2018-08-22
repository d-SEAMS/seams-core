#ifndef _MOLECULE_H     
#define _MOLECULE_H

#include <iostream>
using namespace std;

const int MAXNUMBEROFNEIGHBORS = 100;

class CMolecule {
  private:
    //The position of the particle
  public:
    CMolecule();
    virtual ~CMolecule();
    double posx,posy,posz;
    double vx,vy,vz;
    double fx,fy,fz;
    double rfx,rfy,rfz;
    double potential;
    
    void set_position(double,double,double);
    double get_posx();
    double get_posy();
    double get_posz();
    int neighbors[MAXNUMBEROFNEIGHBORS];
    double neighbordist[MAXNUMBEROFNEIGHBORS];
    int n_neighbors;
};

#endif
