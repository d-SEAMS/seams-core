#ifndef _STRUCTURE_FACTOR_H
#define _STRUCTURE_FACTOR_H

#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "molecule.h"
#include "parameter.h"
#include "molecular_system.h"
#include "output.h"
#include "rdf2D.h"

/*! \brief Class for Structure Factor Calculation.
 *         This class creates an object for the structure factor, 
  and the output can be printed to a file
 *
 For example, if you want to find the in-plane RDF in the
 XY plane, you should use the  initRDFxy() function, followed by
 either the accumulateRDFxy() or singleRDFxy(), depending on whether
 you want to average over several frames or not. If you use the accumulateRDFxy
 function, you must also use the normalizeRDF2D() function. Don't call 
 call normalizeRDF2D() when using singleRDFxy(), since this function
 is called by singleRDFxy() itself. The printRDF2D() function prints
 the radial values and corresponding RDF values to a file called rdf2D.dat
 in the output folder. 

 The structure factor can be determined in two ways: directly from the coordinates
 or as a Fourier transform of the pair correlation function \f$g(r)\f$.

 The equation used for \f$k \neq 0\f$ is:
  \f[
  g^n(r) = \frac{1}{(\rho^n)^2 A \delta z} \Sigma_{i \neq j} \delta(r - r_{ij}) \left[ \Theta\left( \frac{\delta z}{2}-|z_i-z^2| \right) \times \Theta\left( \frac{\delta z}{2}-|z_j-z^n| \right) \right] 
  \f]

  For \f$k=0\f$ use the exponential term
  - \f$z^n\f$ is the z coordinate of the layer chosen. This can be determined
  by finding the peak in the density(or number) vs \f$z\f$ plot.
  - \f$\delta z\f$ is the thickness of the layer chosen. This can also be determined
  from the density(or number) vs \f$z\f$ plot by ascertaining the width of the peak corresponding to the desired
  layer. 
  - Here, \f$\rho\f$ corresponds to the bulk number density of the species chosen. This is determined in the code
  by calculating the number of atoms and the volume. By default, the program accepts the simulation box volume
  as the volume. The user should enter the volume for more complex systems. The accuracy of the RDF is sensitive to 
  the calculated value of density, so care should be taken when specifying/calculating the density. 
  - This calculation is computed upto a maximum radius, which can either be user-defined or is set by default to be half 
  the smallest box length of the \f$x\f$ or \f$y\f$ dimensions.
 */

class StructureFactor: public COutput {
	private:		
    // No. of bins of the RDF
		int nbin; 
    // No. of bins of the structure factor
    int sbin;
    // Wave vector increment, delta k
    double kwidth;

    // Calculate the number of bins, with box lengths as arguments
    void getBins(double, double);

    // Take the Fourier transform of g(r) when k is not zero
    void fourierTransform(class Rdf2D& rdf);

    // Calculate structure factor for k=0 
    void transformAtZero(class Rdf2D& rdf);
 
    // Initialize the structure factor array to zero 
    void initToZero();

    // Get values of the k inverse distance coordinate
    void getK();

    // Helper functions
    // Returns the largest value
    double largest(double, double, double);
    double largest(double, double); 

  public:
		  //the main object in which the structure factor is created and calculated
    	StructureFactor();
    	virtual ~StructureFactor();

    	// Dynamically allocated array for histogram values
    	// for RDF and radial values
		  double* strucFactor;
		  double* k;

    	// Initialize the histogram
    	void initStrucFactor(class Rdf2D& rdf, double box_length1, double box_lenth2);

    	// Print the structure factor to a file in the output folder
    	void printStrucFactor();

    	// Free the memory 
    	void deleteStrucFactor();

};

#endif