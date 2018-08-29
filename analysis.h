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

using namespace std;

class CAnalysis: private CParameter {
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

    	// Calculates the RDF over a number of snapshots
    	void accumulateRDF3D(class CMolecularSystem& molSys);
    	// Calculates the RDF for a single snapshot
    	void calcRDF3D(class CMolecularSystem& molSys);
    	// Normalizes the RDF 
    	void normalizeRDF3D();
	public:
		//the main object where all properties of all particles are saved
    	CAnalysis();
    	virtual ~CAnalysis();
    	void readParameter(class CMolecularSystem& molSys);

    	// Dynamically allocated array for histogram values
		int* rdf3D;

    	// Initialize the histogram
    	void initRDF3D(class CMolecularSystem& molSys);
    	// Check to make sure that the user-defined max_radius is within limits, assign volume
    	void checkParameter(class CMolecularSystem& molSys);
    	// Calculate the number of bins
    	void getBins();
    	void deleteRDF3D();

		// Get absolute relative distance from wrapped coordinates
		double getAbsDistance(int, int, class CMolecularSystem& molSys);

		// Helper functions
		// Returns the smallest value
		double smallest(double, double, double); 
		double smallest(double, double);
};

#endif