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
 Use 
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

    	// Calculate the histogram of the 3D RDF
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
    	void normalizeRDF2D(double);
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