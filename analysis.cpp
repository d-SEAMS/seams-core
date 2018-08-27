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

double CAnalysis::getAbsDistance(double x, double y, class CMolecularSystem& molSys)
{
	double box = molSys.parameter->boxx;
	return box;
}