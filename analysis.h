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
		// Dynamically allocated array for histogram values
		int* rdf3D;
		// No. of snapshots for RDF
		int nframes;		
    	// No. of bins 
		int nbin; 
		// User-specified binwidth
    	double binwidth;
 		// Max distance upto which calculation is done
    	double max_radius;
	public:
		//the main object where all properties of all particles are saved
    	CAnalysis();
    	virtual ~CAnalysis();
    	void readParameter(class CMolecularSystem& molSys);

    	// Initialize the histogram
    	void initHistogram(class CMolecularSystem& molSys);
    	void checkRadius(class CMolecularSystem& molSys);

		// Get absolute relative distance from wrapped coordinates
		double getAbsDistance(int, int, class CMolecularSystem& molSys);

		// Helper functions
		// Returns the smallest value
		double smallest(double, double, double); 
		double smallest(double, double);
};

#endif