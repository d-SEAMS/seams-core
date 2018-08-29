#include "analysis.h"
#include "molecular_system.h"
#include "molecule.h"

const string PF_BINWIDTH = "binwidth";
const string PF_RADIUS = "max_radius";

// Constructor
CAnalysis::CAnalysis()
{
  this->binwidth = -1.0;
  this->max_radius = -1.0;
}

CAnalysis::~CAnalysis()
{
  delete [] rdf3D;
}

// Initialize the histogram array
void CAnalysis::initHistogram(class CMolecularSystem& molSys)
{
	//
}




//****************************************************************************************
//This procedure reads the parameter.txt file. The keywords are defined above with PF_...
//if a line starts with // it is handled as comment
//do not have spaces before or after =
//****************************************************************************************
void CAnalysis::readParameter(class CMolecularSystem& molSys)
{
	ifstream paraFile;
  	// Open the parameter file
  	paraFile.open("input/parameter.txt");
  	string line;
  	string::size_type pos;
  	int i = 0;
  	while (getline(paraFile,line))
  	{
  	  if(line.substr(0, 2).compare("//")!=0)
  	  {
  	    pos  = line.find('=');
  	    if (pos != string::npos)
  	    {
  	      this->rawParameter[i].name = line.substr(0, pos );
  	      this->rawParameter[i].value = line.substr(pos+1, string::npos );
        i += 1;
    	  } else {if (line.compare("")>0) {cerr << "malformed line in parameterfile :" << line << "\n";}}
    	}
  	}
  	for (int j = 0;j < i;j++)
  	{
    if (rawParameter[j].name.compare(PF_BINWIDTH) == 0) {this->binwidth = atof(rawParameter[j].value.c_str());}
    if (rawParameter[j].name.compare(PF_RADIUS) == 0) {this->max_radius = atof(rawParameter[j].value.c_str());}
  	}

  	// Check that the max_radius is within simulation box limits
  	this->checkRadius(molSys);
}

//-------------------------------------------------------------------------------------------------------
// CHECKS AND HELPER FUNCTIONS
//-------------------------------------------------------------------------------------------------------

// Checks that the max_radius entered is correct. If the max_radius is greater than half the simulation
// box length, by default it is set as half the smallest box length
// TODO: Modify for 2-D case
// TODO: Check binwidth
void CAnalysis::checkRadius(class CMolecularSystem& molSys)
{
  double boxx, boxy, boxz;
  double half_box; // Half the smallest box length
  double radius = this->max_radius;

  // Box lengths 
  boxx = molSys.parameter->boxx;
  boxy = molSys.parameter->boxy;
  boxz = molSys.parameter->boxz;

  half_box = 0.5*this->smallest(boxx, boxy, boxz);

  // Check if the max_radius is within bounds
  if (radius > half_box || radius <= 0.0)
  {
    cerr << "The maximum possible value of radius is " << half_box << " and in the parameter.txt file, it was set as " << radius << "\n";
    cerr << "I will now set the maximum radius to half the simulation box length " <<"\n";
    this->max_radius = half_box;
  }
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


//****************************************************************************************
// HELPER FUNCTIONS
//****************************************************************************************

double CAnalysis::smallest(double x, double y, double z)
{
  // return min({x,y,z}); // For C++11
  return min(min(x,y), z);
}

double CAnalysis::smallest(double x, double y)
{
  return min(x,y);
}
