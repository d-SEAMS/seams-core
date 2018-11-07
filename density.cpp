#include "density.h"

/********************************************//**
 *  Constructor
 ***********************************************/
Density::Density()
{
	this->binwidth = -1.0;
	this->typeI = -1;
	this->number = NULL; 
	this->coord = NULL; 
}
/********************************************//**
 *  Destructor
 ***********************************************/
Density::~Density()
{
	delete [] number;
	delete [] coord;
}


/********************************************//**
 *  Initializes the number density. The binning dimension 
 is the z dimension
 ***********************************************/
void Density::initNumberZ(class CMolecularSystem& molSys,double binwidth, int typeI, double xlo, double xhi, double ylo, double yhi, double zlo, double zhi)
{
	double z_length = zhi-zlo; // Length for binning
	// If the user has not entered z dimension lengths, the box length is taken as default
	if (z_length==0.0){z_length=molSys.parameter->boxz;}

	// Assign values
	this->binwidth = binwidth;
	this->typeI = typeI;
	this->xlo=xlo; this->xhi = xhi; 
	this->ylo=ylo; this->yhi=yhi;
	this->zlo=zhi; this->zhi=zhi;


}


/********************************************//**
 *  Gets the number density of a particle of type I
 in a user-defined volume, for a single frame? Same func for all. If the volume is not defined, the entire 
 box volume is taken. Binning is done in the z dimension
 ***********************************************/
void Density::NumberSingleFrameZ(class CMolecularSystem& molSys,double binwidth, int typeI,double xlo,double xhi,double ylo,double yhi,double zlo,double zhi)
{
	//
}


/********************************************//**
 *  Gets the number of bins in the binning dimension.
 Takes maximum length in which binning will be done, and the CMolecularSystem object
 as the arguments
 ***********************************************/
void Density::getBins(class CMolecularSystem& molSys, double max_length)
{
	//
}