#ifndef _RDF2D_H
#define _RDF2D_H

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

/*! \brief Class for 2D RDF.
 *         This class creates an object for the 2D-RDF, 
  and the output of the RDF can be printed to a file
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

 The equation used for 2D-RDF for the \f$n^{th}\f$ layer is:
  \f[
  g^n(r) = \frac{1}{(\rho^n)^2 A \delta z} \Sigma_{i \neq j} \delta(r - r_{ij}) \left[ \Theta\left( \frac{\delta z}{2}-|z_i-z^2| \right) \times \Theta\left( \frac{\delta z}{2}-|z_j-z^n| \right) \right] 
  \f]
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

class Rdf2D: public COutput {
	private:
		// No. of snapshots for RDF
		int nframes;		
    	// No. of bins 
		int nbin; 
		// User-specified binwidth
    	double binwidth;
 		// Max distance upto which calculation is done
    	double max_radius;
    	// Total volume required for density calculation
    	double volume;
    	// Total no. of atoms
    	int nop;

    	// Calculate the histogram of the 2D RDF in the XY plane
    	void histogramRDFxy(class CMolecularSystem& molSys, double z_layer, double dz);

    	// Calculate the number of bins
    	void getBins();

        // Calculate the number of atoms in the box for the given frame and IDs
        int getNatoms(class CMolecularSystem& molSys, int, int);
        // Calculates the number of atoms in the XY plane
        int getNatomsXY(class CMolecularSystem& molSys, double, double);

        // Check whether the z-coordinate is within the layer
        bool atomInsideLayer(double z, double z_layer, double dz);

    	// Check to make sure that the user-defined max_radius is within limits
    	void checkParameterXY(class CMolecularSystem& molSys);
    	// Assigns volume
    	void checkVolume(class CMolecularSystem& molSys); 
    	// Initialize the 2D RDF array to zero before histogramming
    	void rdf2DInitToZero();
    	// Get absolute relative distance in the XY plane from wrapped coordinates
		double absDistanceXY(int, int, class CMolecularSystem& molSys);
		// Helper functions
		// Returns the smallest value
		double smallest(double, double, double); 
		double smallest(double, double);

        // Functions for YZ plane 
        // Calculate the histogram of the 2D RDF in the YZ plane
        void histogramRDFyz(class CMolecularSystem& molSys, double x_layer, double dx);
	    // Calculates the number of atoms in the YZ plane
        int getNatomsYZ(class CMolecularSystem& molSys, double, double);
    public:
		//the main object where all properties of all particles are saved
    	Rdf2D();
    	virtual ~Rdf2D();

    	// Dynamically allocated array for histogram values
    	// for RDF and radial values
		double* rdf2D;
		double* rVal;

        // Lammps trajectory IDs of the atoms to compute the RDF 
        // If not set, RDF for all atoms is calculated
        int typeA;
        int typeB; 

    	// Initialize the histogram
    	void initRDFxy(class CMolecularSystem& molSys, double binwidth, double volume=-1.0, double max_radius=-1.0);
    	// Calculates the RDF for a single snapshot
    	void singleRDFxy(class CMolecularSystem& molSys, double z_layer, double dz, int typeA=-1, int typeB=-1);
    	// Calculates the RDF over a number of snapshots
    	void accumulateRDFxy(class CMolecularSystem& molSys, double z_layer, double dz, int typeA=-1, int typeB=-1);
    	// Normalizes the RDF. You don't need to call this separately 
    	// for calculation of RDF for a single frame. You must call this 
    	// after using the accumulate RDF command for multiple snapshots
    	void normalizeRDF2D(double dr);
    	// Get the radial values corresponding to each radial bin
    	void getR();
        // Reintialize the histogram and number of frames to zero
        void clearRDF2D();

    	// Print the 3D RDF to a file in the output folder
    	void printRDF2D();

    	// Free the memory 
    	void deleteRDF2D();

};

#endif