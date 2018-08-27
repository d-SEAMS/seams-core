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

double CAnalysis::getAbsDistance(double x, double y)
{
	double box = CMolecularSystem::getBoxx();
	return box;
}