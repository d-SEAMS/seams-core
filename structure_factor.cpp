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
void StructureFactor::initStrucFactor(class Rdf2D& rdf)
{
    // Get the number of bins 
    this->nbin = rdf.binsInRDF();
    // Initialize the array for structure factor
    this->strucFactor  = new double[this->nbin];
    // Initialize the RDF array to zero
    this->initToZero();
    // Get an array for k
    this->k = new double[this->nbin];
    // Get the values of the inverse distance k
    this->getK();
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
    for (int ibin=0; ibin<this->nbin; ibin++)
    {
        this->k[ibin] = ibin;
    }
}


/********************************************//**
 *  Initialize the structure factor array to zero
 ***********************************************/
void StructureFactor::initToZero()
{
    // Loop over all bins
    for (int ibin=0; ibin < this->nbin; ibin++)
    {
        this->strucFactor[ibin] = 0.0;
    }
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
