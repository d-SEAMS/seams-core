#include "layer.h"
#include "molecular_system.h"
#include "molecule.h"

/********************************************//**
 *  Constructor
 ***********************************************/
CLayer::CLayer()
{
	// Volume defined by the user
	this->xlo = 0.0;
	this->xhi = 0.0;
	this->ylo = 0.0;
	this->yhi = 0.0;
	this->zlo = 0.0;
	this->zhi = 0.0;
	this->notset = true;
	this->volFlag = 0;
}
/********************************************//**
 *  Destructor
 ***********************************************/
CLayer::~CLayer()
{
}

/********************************************//**
 *  Checks if the atom is inside the layer with a 
 midpoint of r_layer.
 ***********************************************/ 
bool CLayer::atomInsideLayer(double r, double r_layer, double dr)
{
	if ((r-r_layer)<=0.5*dr){return true;}
	else {return false;}
}

/********************************************//**
 *  Checks if the atom is inside the user-defined volume
 ***********************************************/ 
bool CLayer::atomInsideVol(class CMolecularSystem& molSys, int iatom, double xlo, double xhi, double ylo, double yhi, double zlo, double zhi)
{
	double x_atom;
	double y_atom;
	double z_atom; 

	this->xlo=xlo; this->xhi=xhi; this->ylo=ylo; this->yhi=yhi; this->zlo=zlo; this->zhi=zhi;
	if (xlo==xhi && ylo==yhi && zlo==zhi){this->notset=true; return true;} 

	// Check that the user-defined volume is within the box limits
	this->checkVolume();
	if (this->notset==true){return true;} 

	// Now that the volume entered is feasible, check if the atom is within 
	// the limits
	x_atom = molSys.molecules[iatom].get_posx();
	y_atom = molSys.molecules[iatom].get_posy();
	z_atom = molSys.molecules[iatom].get_posz();

	this->atomCoordLimits(x_atom, this->xlo, this->xhi);
	this->atomCoordLimits(y_atom, this->ylo, this->yhi);
	this->atomCoordLimits(z_atom, this->zlo, this->zhi);

	if(volFlag == 3){return true;}
	else {return false;}
}

/********************************************//**
 *  Checks if the user-entered volume is correct and
 updates box limits
 ***********************************************/
void CLayer::checkVolume()
{
	// Volume lengths
	double x_length = this->xhi - this->xlo;
	double y_length = this->yhi - this->ylo;
	double z_length = this->zhi - this->zlo;

	if (x_length < 0.0 || y_length < 0.0 || z_length < 0.0)
		{
			std::cerr<<"You have entered an incorrect volume.\n";
			this->notset=true;
		}
	else {this->notset=false;}
}

/********************************************//**
 *  Check if the atom is within a particular dimension range
 ***********************************************/
void CLayer::atomCoordLimits(double r_atom, double r_min, double r_max)
{
	if (r_min == 0 && r_max == 0){this->volFlag += 1;}
	else if (r_atom >=r_min && r_atom <= r_max){this->volFlag += 1;}
}