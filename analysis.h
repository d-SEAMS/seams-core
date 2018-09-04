#ifndef _ANALYSIS_H
#define _ANALYSIS_H

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

class CAnalysis: private CParameter, public COutput {
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
    	// Total no. of atoms
    	int nop;

    	// Calculate the histogram of the 3D RDF
    	void histogramRDF3D(class CMolecularSystem& molSys);

    	// Calculate the number of bins
    	void getBins();

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
    	CAnalysis();
    	virtual ~CAnalysis();
    	void readParameter(class CMolecularSystem& molSys);

    	// Dynamically allocated array for histogram values
    	// for RDF and radial values
		double* rdf3D;
		double* rVal; 

    	// Initialize the histogram
    	void initRDF3D(class CMolecularSystem& molSys);
    	// Calculates the RDF for a single snapshot
    	void singleRDF3D(class CMolecularSystem& molSys);
    	// Calculates the RDF over a number of snapshots
    	void accumulateRDF3D(class CMolecularSystem& molSys);
    	// Normalizes the RDF. You don't need to call this separately 
    	// for calculation of RDF for a single frame. You must call this 
    	// after using the accumulate RDF command for multiple snapshots
    	void normalizeRDF3D();
    	// Get the radial values corresponding to each radial bin
    	void getR();

    	// Print the 3D RDF to a file in the output folder
    	void printRDF3D();

    	// Free the memory 
    	void deleteRDF3D();

};

#endif