#ifndef _MOLECULE_H     
#define _MOLECULE_H

#include <iostream>

class CMolecule {
  private:
    //The position of the particle
    double posx,posy,posz;
  public:
    CMolecule();
    virtual ~CMolecule();
    double vx,vy,vz;
    double fx,fy,fz;
    double rfx,rfy,rfz;
    double potential;
    
    void set_position(double,double,double);
    double get_posx();
    double get_posy();
    double get_posz();
    // Lammps trajectory type ID
    int type;
    // Lammps molecule ID
    int molID;
    // Lammps atom ID
    int atomID;
};

#endif
