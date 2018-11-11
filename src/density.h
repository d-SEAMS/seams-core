#ifndef _DENSITY_H
#define _DENSITY_H

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
#include "geometry.h"

class Density: public COutput, private CVolume {
	private:
		// No. of snapshots for density binning 
		int nframes;		
    	// No. of bins 
		int nbin;
		// User-specified binwidth
    	double binwidth;

		// Initializes the arrays etc
		void initNumberZ(class CMolecularSystem& molSys,double binwidth, int typeI, double xlo, double xhi, double ylo, double yhi, double zlo, double zhi);

		// Calculate the number of bins given the length in which binning is to be done
    	void getBins(class CMolecularSystem& molSys, double );

	public:
		Density();
		virtual ~Density();

		// Dynamically allocated array for number distribution and binned coordinate
		double* number;
		double* coord;

		// Lammps trajectory IDs of the atom
		int typeI;

		// Gets the number density of type I in the defined volume limits; binning done in the z dimension
		void NumberSingleFrameZ(class CMolecularSystem& molSys, double binwidth, int typeI=-1, double xlo = 0.0, double xhi = 0.0, double ylo=0.0, double yhi=0.0, double zlo=0.0, double zhi=0.0);
};

#endif