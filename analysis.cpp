#include "analysis.h"
#include "molecular_system.h"
#include "molecule.h"

CAnalysis::CAnalysis()
{
  //
}

CAnalysis::~CAnalysis()
{
  //
}

//-------------------------------------------------------------------------------------------------------
// DISTANCE CALCULATIONS
//-------------------------------------------------------------------------------------------------------

// Returns the absolute distance between two particles
// with particle indices iatom and jatom (x[iatom] - x[jatom])
double CAnalysis::getAbsDistance(int iatom, int jatom, class CMolecularSystem& molSys)
{
	double dr[3]; // Relative distance between wrapped coordinates
	double box[3] = {molSys.parameter->boxx, molSys.parameter->boxy, molSys.parameter->boxz};
	double r2 = 0.0; // Squared absolute distance

	// Get the relative distance in the x, y, z dim
	dr[0] = molSys.molecules[iatom].get_posx() - molSys.molecules[jatom].get_posx();
	dr[1] = molSys.molecules[iatom].get_posy() - molSys.molecules[jatom].get_posy();
	dr[2] = molSys.molecules[iatom].get_posz() - molSys.molecules[jatom].get_posz();

	// Get the squared absolute distance
	for (int k=0; k<3; k++)
	{
		// Correct for periodicity
		dr[k] -= round(dr[k]/box[k]);
		
		r2 += pow(dr[k],2.0);
	}
	

	return sqrt(r2);
}