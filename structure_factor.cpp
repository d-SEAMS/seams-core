#include "structure_factor.h"
#include "molecular_system.h"
#include "molecule.h"

// Constructor
StructureFactor::StructureFactor()
{
  this->nbin=-1;
}

StructureFactor::~StructureFactor()
{
  delete [] strucFactor;
  delete [] k;
}

/********************************************//**
 *  Initializes the histogram array for structure factor

 It takes the CMolecularSystem object reference and binwidth as arguments.
 Optional arguments include the maximum radius upto which 
 the 2-D RDF will be calculated and the desired volume. If not set, the default
 values are half the simulation box and the volume of the simulation box respectively
 ***********************************************/
void StructureFactor::initStrucFactor(class Rdf2D& rdf, double box_length1, double box_lenth2)
{
    // Get the number of bins for the RDF 
    this->nbin = rdf.binsInRDF();
    // Get the number of bins for the structure factor
    this->getBins(box_length1, box_lenth2);
    // Initialize the array for structure factor
    this->strucFactor = new double[this->sbin];
    // Initialize the structure factor array to zero
    this->initToZero();
    // Get an array for k
    this->k = new double[this->sbin];
    // Get the values of the inverse distance k
    this->getK();
    // Structure factor calculation
    // Print to file
    this->printStrucFactor();
}

/********************************************//**
 *  Frees the memory 
 ***********************************************/
void StructureFactor::deleteStrucFactor()
{
    delete [] strucFactor;
    delete [] k;
}

/********************************************//**
 *  Get the inverse distance k values corresponding to each
 structure factor value
 ***********************************************/
void StructureFactor::getK()
{
    // Loop through all bins
    for (int ibin=0; ibin<this->sbin; ibin++)
    {
        this->k[ibin] = (ibin+1)*this->kwidth;
    }
}


/********************************************//**
 *  Initialize the structure factor array to zero
 ***********************************************/
void StructureFactor::initToZero()
{
    // Loop over all bins
    for (int ibin=0; ibin < this->sbin; ibin++)
    {
        this->strucFactor[ibin] = 0.0;
    }
}

/********************************************//**
 *  Gets the number of bins for the structure factor,
 by calculating the width of the wave vector kwidth
 from the box dimensions using \f$\delta k = \frac{2 \pi}{L}\f$
 ***********************************************/
void StructureFactor::getBins(double box_length1, double box_lenth2)
{
    // Get the largest box length
    double length = this->largest(box_length1, box_lenth2);
    // kwidth is the should be such that the amplitude is not larger than the box length
    this->kwidth = 2*PI/length;
    this->sbin = 800; // temp
    this->sbin = this->sbin & ~1; // Round down to an even number for Simpson's rule
}

//-------------------------------------------------------------------------------------------------------
// CALCULATIONS
//-------------------------------------------------------------------------------------------------------

/********************************************//**
 *  Calculates the structure factor
 ***********************************************/
void StructureFactor::calcStrucFactor()
{
  //
}

/********************************************//**
 *  Calculates the value of the integral \f$\int\f$
 Takes the value of k as the argument 
 ***********************************************/
double StructureFactor::integrateSimpsons(double k)
{
    //
}

//-------------------------------------------------------------------------------------------------------
// OUTPUT FUNCTIONS
//-------------------------------------------------------------------------------------------------------

/********************************************//**
 *  Prints out the 2D RDF function to a file in the output folder
 ***********************************************/
void StructureFactor::printStrucFactor()
{
    // Prints the radial values and 3D RDF values to a file called rdf3D.txt
    this->printToFile(this->nbin, this->k, this->strucFactor, "strucFactor", "Inverse Distance Coordinate k", "S(k)");
} 

//-------------------------------------------------------------------------------------------------------
// HELPER FUNCTIONS
//-------------------------------------------------------------------------------------------------------

// Functions for returning the largest number
double StructureFactor::largest(double x, double y, double z)
{
  // return std::min({x,y,z}); // For C++11
  return std::min(std::min(x,y), z);
}

double StructureFactor::largest(double x, double y)
{
  return std::min(x,y);
}

