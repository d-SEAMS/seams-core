#include "molecule.h"

//Constructor
CMolecule::CMolecule()
{
	// Default type ID
	this->type = 1;
}
// Destructor
CMolecule::~CMolecule()
{
}


//Sets the position of the molecule
void CMolecule::set_position(double x, double y, double z)
{
  this->posx = x;
  this->posy = y;
  this->posz = z;
}

//Get the x-Component of the position
double CMolecule::get_posx()
{
  return this->posx;
}

//Get the y-Component of the position
double CMolecule::get_posy()
{
  return this->posy;
}

//Get the z-Component of the position
double CMolecule::get_posz()
{
  return this->posz;
}
