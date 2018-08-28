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

class CAnalysis {
	private:
		// Dynamically allocated array for histogram values
		int* rdf3D;
		// No. of snapshots for RDF
		int nframes;
	public:
		//the main object where all properties of all particles are saved
    	CAnalysis();
    	virtual ~CAnalysis();

    	// Initialize the histogram
    	void initHistogram();

    	// No. of bins 
		int nbin;

		// Get absolute relative distance from wrapped coordinates
		double getAbsDistance(int, int, class CMolecularSystem& molSys); 
};

#endif