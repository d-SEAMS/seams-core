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
    this->strucFactor = new double[this->nbin];
    // Initialize the RDF array to zero
    this->initToZero();
    // Get an array for k
    this->k = new double[this->nbin];
    // Get the values of the inverse distance k
    this->getK();
    // Get the value of S(k) at k=0
    this->transformAtZero(rdf);
    // Get S(k) for all other k
    this->fourierTransform(rdf);
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
        this->k[ibin] = ibin*this->kwidth;
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

/********************************************//**
 *  Gets the number of bins for the structure factor,
 by calculating the width of the wave vector kwidth
 from the box dimensions using \f$\delta k = \frac{2 \pi}{L}\f$
 ***********************************************/
void StructureFactor::getBins(double box_length1, double box_lenth2)
{
    // Get the largest box length
    double length = this->largest(box_length1, box_lenth2);
    // kwidth is the 
}

//-------------------------------------------------------------------------------------------------------
// CALCULATIONS
//-------------------------------------------------------------------------------------------------------

/********************************************//**
 *  Calculates the structure factor for \f$k=0\f$

 It uses the exponential form
 ***********************************************/
void StructureFactor::transformAtZero(class Rdf2D& rdf)
{
  double sum=0.0; // Variable used for calculating integral approximated as summation
  double rho=rdf.rho; // Number density
  // Loop through all the nbin points of g(r)
  // For k=0, the exponential term is 1
  for (int ibin=0; ibin < this->nbin; ibin++)
  {
    sum += (rdf.rdf2D[ibin] -1);
  }
  // Multiply by rho and add 1
  this->strucFactor[0] = 1 + rho*sum;
}

/********************************************//**
 *  Calculates the structure factor for \f$k \neq 0\f$

 It uses the sinc
 ***********************************************/
void StructureFactor::fourierTransform(class Rdf2D& rdf)
{
  double sum=0.0; // Variable used for calculating integral approximated as summation
  double rho=rdf.rho; // Number density
  double r; // Radial coordinate
  
  for (int ibin=0; ibin < this->nbin; ibin++)
  {
    // Loop through all the nbin points of g(r)
    for (int kbin=1; kbin < this->nbin; kbin++)
    {
      r = rdf.rVal[kbin]; // radial value 
      sum += ((rdf.rdf2D[kbin]-1)*r*sin(kbin*r)/(double(kbin)));
    }
    // Multiply by rho and add 1
    std::cout<<"sum is "<<sum<<"\n";
    this->strucFactor[ibin] = 1 + 4*PI*rho*sum;
    std::cout<<"s(k) "<<strucFactor[ibin]<<"\n";
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

