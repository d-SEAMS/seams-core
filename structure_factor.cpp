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
    this->calcStrucFactor(rdf);
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
        this->k[ibin] = this->k_min + (ibin+1)*this->kwidth;
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
void StructureFactor::getBins(double box_length1, double box_length2)
{
    // Get the largest possible length
    // double max_length = this->smallest(box_length1, box_length2);
    double max_length = sqrt(pow(box_length1,2) + pow(box_length2, 2));
    // kwidth is the should be such that the amplitude is not larger than the box length
    this->kwidth = 2*PI/max_length;
    // this->kwidth = 2*PI/(sqrt(3)*0.154);
    this->k_min = 2*PI/3.1589;
    // double k_max = 2*PI/0.075;
    double k_max = 20;
    this->sbin = int((k_max-k_min)/kwidth); // Using s_max
    std::cout<< "kwidth is " << this->kwidth << " and k_max is " << k_max << "\n";
}

//-------------------------------------------------------------------------------------------------------
// CALCULATIONS
//-------------------------------------------------------------------------------------------------------

/********************************************//**
 *  Calculates the structure factor
 ***********************************************/
void StructureFactor::calcStrucFactor(class Rdf2D& rdf)
{
    // Loop through all k 
    for (int kbin=0; kbin < this->sbin; kbin++)
    {
        this->strucFactor[kbin] = 1 + 4*PI*rdf.rho*this->integrateSimpsons(rdf, this->k[kbin]);
    }
}

/********************************************//**
 *  Calculates the value of the integral \f$\int\f$
 Takes the value of k as the argument 
 ***********************************************/
double StructureFactor::integrateSimpsons(class Rdf2D& rdf, double k)
{
    int nbin;
    double sum = 0.0;                   // Summation for the Simpson's rule
    double h = rdf.rVal[1]-rdf.rVal[0]; // Step size in r
    // Divide the range of g(r) into an even number of subintervals
    if (this->nbin%2 == 0){nbin = this->nbin-1;}
    else {nbin = this->nbin-1;}

    // Apply Simpson's 1/3 rule
    for (int ibin=1; ibin <= nbin/2; ibin++)
    {
        sum += (this->integrand(rdf, 2*ibin-2, k) + 4*this->integrand(rdf, 2*ibin-1, k) + this->integrand(rdf, 2*ibin, k));
    }

    return h*sum/3.0;
}

/********************************************//**
 *  Function to return the integrand for a particular \f$g(r)\f$ value
 ***********************************************/
double StructureFactor::integrand(class Rdf2D& rdf, int index, double k)
{
    double r = rdf.rVal[index]; // r at a particular index value
    double g_r = rdf.rdf2D[index]; // g(r)
    return ((g_r-1)*r*sin(k*r)/k);
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
    this->printToFile(this->sbin, this->k, this->strucFactor, "strucFactor", "Inverse Distance Coordinate k", "S(k)");
} 

//-------------------------------------------------------------------------------------------------------
// HELPER FUNCTIONS
//-------------------------------------------------------------------------------------------------------

// Functions for returning the largest number
double StructureFactor::largest(double x, double y, double z)
{
  // return std::min({x,y,z}); // For C++11
  return std::max(std::min(x,y), z);
}

double StructureFactor::largest(double x, double y)
{
  return std::max(x,y);
}

// Functions for returning the smallest number
double StructureFactor::smallest(double x, double y)
{
  return std::min(x,y);
}
