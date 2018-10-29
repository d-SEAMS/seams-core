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
		// 

	public:
		Density();
		virtual ~Density();

		// Dynamically allocated array for number distribution and binned coordinate
		double* number;
		double* coord;
};

#endif