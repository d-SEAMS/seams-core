#ifndef _RDF3D_H
#define _RDF3D_H

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

/*! \brief Class for 3D RDF.
 *         This class creates an object for the 3D-RDF, 
  and the output of the RDF can be printed to a file
 *
 Use initRDF3D() to create an object for the 3D RDF.
 Depending on the number of frames, use singleRDF3D() or accumulateRDF3D() .
 In case you have used the accumulate function, you will have to use the normalizeRDF3D()
 to normalize the RDF. 
 Finally, print the output to an output file, which is called rdf3D.txt by default using
 the print function
 */

class Rdf3D: public COutput {
	private:
		// No. of snapshots for RDF
		int nframes;		
    	// No. of bins 
		int nbin; 
		// User-specified binwidth
    	double binwidth;
 		// Max distance upto which calculation is done
    	double max_radius;
    	// Total volume in which 3D RDF is to be calculated
    	double volume;
        // Total number of central atoms (Itype)
        int n_iatoms;
        // Total number of distribution atoms (Jtype)
        int n_jatoms;

        // Dynamically allocated array for indices with particles of type I and J
        int* iIndex;
        int* jIndex;

    	// Calculate the histogram of the 3D RDF for I-I similar type calculations
    	void histogramRDF3Dii(class CMolecularSystem& molSys);
        // Calculate the histogram of the 3D RDF for I-J dissimilar type calculations
        void histogramRDF3Dij(class CMolecularSystem& molSys);

    	// Calculate the number of bins
    	void getBins();

        // Calculate the number of atoms in the box for the given frame and IDs
        void getNatoms(class CMolecularSystem& molSys, int, int);

    	// Check to make sure that the user-defined max_radius is within limits, assign volume
    	void checkParameter(class CMolecularSystem& molSys);
    	// Initialize the 3D RDF array to zero before histogramming
    	void rdf3DInitToZero();
    	// Get absolute relative distance from wrapped coordinates
		double getAbsDistance(int, int, class CMolecularSystem& molSys);
		// Helper functions
		// Returns the smallest value
		double smallest(double, double, double); 
		double smallest(double, double);
	public:
		//the main object where all properties of all particles are saved
    	Rdf3D();
    	virtual ~Rdf3D();

    	// Dynamically allocated array for histogram values
    	// for RDF and radial values
		double* rdf3D;
		double* rVal;

        // Lammps trajectory IDs of the atoms to compute the RDF 
        // If not set, RDF for all atoms is calculated
        int typeI;
        int typeJ; 

        // Volume bulk density of the system
        double rho;

    	// Initialize the histogram
    	void initRDF3D(class CMolecularSystem& molSys, double binwidth, double volume=-1.0, double max_radius=-1.0);
    	// Calculates the RDF for a single snapshot
    	void singleRDF3D(class CMolecularSystem& molSys, int typeA=-1, int typeB=-1);
    	// Calculates the RDF over a number of snapshots
    	void accumulateRDF3D(class CMolecularSystem& molSys, int typeA=-1, int typeB=-1);
    	// Normalizes the RDF. You don't need to call this separately 
    	// for calculation of RDF for a single frame. You must call this 
    	// after using the accumulate RDF command for multiple snapshots
    	void normalizeRDF3D();
    	// Get the radial values corresponding to each radial bin
    	void getR();
        // Reintialize the histogram and number of frames to zero
        void clearRDF3D();

    	// Print the 3D RDF to a file in the output folder
    	void printRDF3D();

    	// Free the memory 
    	void deleteRDF3D();

        // Returns the number of bins in the RDF array
        int binsInRDF();

};

#endif