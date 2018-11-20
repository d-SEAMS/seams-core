#include "geometry.h"
#include "molecular_system.h"
#include "molecule.h"
#include <bits/stdc++.h>
#include <cmath>

/********************************************/ /**
 *  Constructor
 ***********************************************/
CLayer::CLayer() {}
/********************************************/ /**
 *  Destructor
 ***********************************************/
CLayer::~CLayer() {}

/********************************************/ /**
 *  Checks if the atom is inside the layer with a 
 midpoint of r_layer.
 ***********************************************/
bool CLayer::atomInsideLayer(double r, double r_layer, double dr) {
  if ((r - r_layer) <= 0.5 * dr) {
    return true;
  } else {
    return false;
  }
}

//-------------------------------------------------------------------------------------------------------
// VOLUME CLASS
//-------------------------------------------------------------------------------------------------------

/********************************************/ /**
 *  Constructor
 ***********************************************/
CVolume::CVolume() {
  // Volume defined by the user
  this->xlo = 0.0;
  this->xhi = 0.0;
  this->ylo = 0.0;
  this->yhi = 0.0;
  this->zlo = 0.0;
  this->zhi = 0.0;
  this->notset = true;
  this->volFlag = 0;
  this->n_iatoms = 0;
  this->n_jatoms = 0;
  // Set pointers to null
  this->iIndex = nullptr;
  this->jIndex = nullptr;
}
/********************************************/ /**
 *  Destructor
 ***********************************************/
CVolume::~CVolume() {
  delete[] iIndex;
  delete[] jIndex;
}

/********************************************/ /**
 *  Makes a list of atom IDs of the particles of type I
 within a user-defined volume. If the volume is not defined, the 
 box volume is taken by default
 ***********************************************/
void CVolume::getAtomListI(class CMolecularSystem &molSys, int typeI,
                           double xlo, double xhi, double ylo, double yhi,
                           double zlo, double zhi) {
  int n_iatoms = 0;
  // int n_jatoms=0;
  int ii = 0; // Current index of array iIndex being filled
  // int jj=0; // Current index of array jIndex being filled

  // Set the limits
  this->xlo = xlo;
  this->xhi = xhi;
  this->ylo = ylo;
  this->yhi = yhi;
  this->zlo = zlo;
  this->zhi = zhi;
  if (xlo == xhi && ylo == yhi && zlo == zhi) {
    this->notset = true;
  }

  // Check that the user-defined volume is within the box limits
  this->checkVolume();

  // Create arrays for vectors holding indices for particles
  // of type I and J
  this->iIndex = new int[molSys.parameter->nop];
  // this->jIndex  = new int[molSys.parameter->nop];

  // Loop through all the atoms and make the list of atom IDs
  // that are of type I.
  for (int iatom = 0; iatom < molSys.parameter->nop; iatom++) {
    // Check if the atom type is of type I
    if (molSys.molecules[iatom].type == typeI || typeI == -1) {
      if (atomInsideVol(molSys, iatom, xlo, xhi, ylo, yhi, zlo, zhi) == true) {
        n_iatoms += 1;
        this->iIndex[ii] = iatom; // Put atom ID in iIndex array
        ii += 1;
      }
    }
  }

  // Set the n_iatoms value

  // Check to make sure that the atom number is not zero
  if (n_iatoms == 0) {
    std::cerr << "You have entered an incorrect type ID\n";
    typeI = -1;
    this->n_iatoms = molSys.parameter->nop;
    for (int iatom = 0; iatom < molSys.parameter->nop; iatom++) {
      this->iIndex[ii] = iatom;
    }
    return;
  }

  // Otherwise update the number of iatoms
  this->n_iatoms = n_iatoms;
  return;
}

/********************************************/ /**
 *  Checks if the atom is inside the user-defined volume, irrespective of type
 ***********************************************/
bool CVolume::atomInsideVol(class CMolecularSystem &molSys, int iatom,
                            double xlo, double xhi, double ylo, double yhi,
                            double zlo, double zhi) {
  double x_atom;
  double y_atom;
  double z_atom;

  this->xlo = xlo;
  this->xhi = xhi;
  this->ylo = ylo;
  this->yhi = yhi;
  this->zlo = zlo;
  this->zhi = zhi;
  if (xlo == xhi && ylo == yhi && zlo == zhi) {
    this->notset = true;
    return true;
  }

  // Check that the user-defined volume is within the box limits
  this->checkVolume();
  if (this->notset == true) {
    return true;
  }

  // Now that the volume entered is feasible, check if the atom is within
  // the limits
  x_atom = molSys.molecules[iatom].get_posx();
  y_atom = molSys.molecules[iatom].get_posy();
  z_atom = molSys.molecules[iatom].get_posz();

  this->atomCoordLimits(x_atom, this->xlo, this->xhi);
  this->atomCoordLimits(y_atom, this->ylo, this->yhi);
  this->atomCoordLimits(z_atom, this->zlo, this->zhi);

  if (volFlag == 3) {
    return true;
  } else {
    return false;
  }
}

/********************************************/ /**
 *  Checks if the user-entered volume is correct and
 updates box limits
 ***********************************************/
void CVolume::checkVolume() {
  // Volume lengths
  double x_length = this->xhi - this->xlo;
  double y_length = this->yhi - this->ylo;
  double z_length = this->zhi - this->zlo;

  if (x_length < 0.0 || y_length < 0.0 || z_length < 0.0) {
    std::cerr << "You have entered an incorrect volume.\n";
    this->notset = true;
  } else {
    this->notset = false;
  }
}

/********************************************/ /**
 *  Check if the atom is within a particular dimension range
 ***********************************************/
void CVolume::atomCoordLimits(double r_atom, double r_min, double r_max) {
  if (r_min == 0 && r_max == 0) {
    this->volFlag += 1;
  } else if (r_atom >= r_min && r_atom <= r_max) {
    this->volFlag += 1;
  }
}

/********************************************/ /**
 *  Checks if the atom is inside the user-defined volume, irrespective of type
 This is used when xlo, xhigh etc have been defined separately and do not need to be checked again
 ***********************************************/
bool CVolume::atomInsideVol(class CMolecularSystem &molSys, int iatom) {
  double x_atom;
  double y_atom;
  double z_atom;

  // Check if the atom is within
  // the limits
  x_atom = molSys.molecules[iatom].get_posx();
  y_atom = molSys.molecules[iatom].get_posy();
  z_atom = molSys.molecules[iatom].get_posz();

  this->atomCoordLimits(x_atom, this->xlo, this->xhi);
  this->atomCoordLimits(y_atom, this->ylo, this->yhi);
  this->atomCoordLimits(z_atom, this->zlo, this->zhi);

  if (volFlag == 3) {
    return true;
  } else {
    return false;
  }
}

//-------------------------------------------------------------------------------------------------------
// GENERIC CLASS
//-------------------------------------------------------------------------------------------------------

/********************************************/ /**
 *  Constructor
 ***********************************************/
CGeneric::CGeneric() {}
/********************************************/ /**
 *  Destructor
 ***********************************************/
CGeneric::~CGeneric() {}

// DISTANCE

/********************************************/ /**
 *  Returns the absolute distance between two particles
 with particle indices iatom and jatom (x[iatom] - x[jatom])
 ***********************************************/
double CGeneric::getAbsDistance(int iatom, int jatom,
                                class CMolecularSystem &molSys) {
  double dr[3]; // Relative distance between wrapped coordinates
  double box[3] = {molSys.parameter->boxx, molSys.parameter->boxy,
                   molSys.parameter->boxz};
  double r2 = 0.0; // Squared absolute distance

  // Get the relative distance in the x, y, z dim
  dr[0] =
      molSys.molecules[iatom].get_posx() - molSys.molecules[jatom].get_posx();
  dr[1] =
      molSys.molecules[iatom].get_posy() - molSys.molecules[jatom].get_posy();
  dr[2] =
      molSys.molecules[iatom].get_posz() - molSys.molecules[jatom].get_posz();

  // Get the squared absolute distance
  for (int k = 0; k < 3; k++) {
    // Correct for periodicity
    dr[k] -= box[k] * round(dr[k] / box[k]);

    r2 += pow(dr[k], 2.0);
  }

  return sqrt(r2);
}

// Overload

double CGeneric::getAbsDistance(int iatom, class CMolecularSystem *frameOne,
                                class CMolecularSystem *frameTwo) {
  double dr[3]; // Relative distance between wrapped coordinates
  double box[3] = {frameOne->parameter->boxx, frameOne->parameter->boxy,
                   frameOne->parameter->boxz};
  double r2 = 0.0; // Squared absolute distance

  // Get the relative distance in the x, y, z dim
  dr[0] = fdim(frameOne->molecules[iatom].get_posx(),
               frameTwo->molecules[iatom].get_posx());
  dr[1] = fdim(frameOne->molecules[iatom].get_posy(),
               frameTwo->molecules[iatom].get_posy());
  dr[2] = fdim(frameOne->molecules[iatom].get_posz(),
               frameTwo->molecules[iatom].get_posz());

  // Get the squared absolute distance
  for (int k = 0; k < 3; k++) {
    // Correct for periodicity
    dr[k] -= box[k] * round(dr[k] / box[k]);

    r2 += pow(dr[k], 2.0);
  }

  return sqrt(r2);
}
