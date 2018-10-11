#ifndef _LAYER_H
#define _LAYER_H

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

class CLayer {
	public:
		CLayer();
		virtual ~CLayer();

		// Volume limits (may or may not be set)
		double xlo, xhi, ylo, yhi, zlo, zhi;
		// This is true if the volume has not been defined by the user (default)
		bool notset;
		// Flag to check if the atom is inside the given volume or not 
		int volFlag;

		// Check whether the atom is inside the user-defined volume
		bool atomInsideVol(class CMolecularSystem& molSys, int iatom, double xlo = 0.0, double xhi = 0.0, double ylo=0.0, double yhi=0.0, double zlo=0.0, double zhi=0.0);
		// Check whether the given coordinate (r) is within the layer (midpoint r_layer)
		bool atomInsideLayer(double r, double r_layer, double dr);

		// Check the volume limits
		void checkVolume();
		// Check if the atom is within a particular dimension range
		void atomCoordLimits(double r_atom, double r_min, double r_max);
};

#endif