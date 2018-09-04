#include "analysis.h"
#include "molecular_system.h"
#include "molecule.h"

const string PF_BINWIDTH = "binwidth";
const string PF_RADIUS = "max_radius";
const string PF_VOLUME = "volume";

// Constructor
CAnalysis::CAnalysis()
{
  this->binwidth = -1.0;
  this->max_radius = -1.0;
  this->volume = -1.0;
  this->nframes = 0;
}

CAnalysis::~CAnalysis()
{
  delete [] rdf3D;
  delete [] rVal;
}

// Initialize the histogram array
void CAnalysis::initRDF3D(class CMolecularSystem& molSys)
{
    // Read in the parameters from the parameter file
    this->readParameter(molSys);
    // Calculate the number of bins from user-defined parameters
    this->getBins();
    // Initialize the array for RDF
    this->rdf3D  = new double[this->nbin];
    // Initialize the RDF array to zero
    this->rdf3DInitToZero();
    // Get an array for R
    this->rVal = new double[this->nbin];
    // Get the values for radial distance (unchanged over nframes)
    this->getR();
}

//Frees the memory
void CAnalysis::deleteRDF3D()
{
    delete [] rdf3D;
    delete [] rVal;
}


// Calculate the number of bins from max_radius and binwidth
void CAnalysis::getBins()
{
    this->nbin = int(this->max_radius/this->binwidth);
}


//****************************************************************************************
//This procedure reads the parameter.txt file. The keywords are defined above with PF_...
//if a line starts with // it is handled as comment
//do not have spaces before or after =
//****************************************************************************************
void CAnalysis::readParameter(class CMolecularSystem& molSys)
{
    std::ifstream paraFile;
    // Open the parameter file
    paraFile.open("input/parameter.txt");
    std::string line;
    std::string::size_type pos;
    int i = 0;
    while (std::getline(paraFile,line))
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
    if (rawParameter[j].name.compare(PF_VOLUME) == 0) {this->volume = atof(rawParameter[j].value.c_str());}
    }

    // Check that the max_radius is within simulation box limits
    this->checkParameter(molSys);
    // Since the number of atoms is constant (TODO: Get nop)
    this->nop = molSys.parameter->nop;
}

//-------------------------------------------------------------------------------------------------------
// CALCULATIONS
//-------------------------------------------------------------------------------------------------------

// Calculates the 3D radial distribution function for a number of snapshots
// There is no need to use calcRDF3D if the RDF is to be calculated over a number of frames
// You will have to call the normalize function separately after accumulating to get the RDF
void CAnalysis::accumulateRDF3D(class CMolecularSystem& molSys)
{
    // Update the number of snapshots calculated
    this->nframes += 1;
    // Add to the RDF histogram
    this->histogramRDF3D(molSys);
    // Call the normalize function separately after accumulating 
    // the histogram over all the desired snapshots
}

// Calculates the 3D radial distribution function for a single snapshot
// Use this only if there is one frame only.
void CAnalysis::singleRDF3D(class CMolecularSystem& molSys)
{
    // There is only one snapshot
    this->nframes = 1;
    // Add to the RDF histogram
    this->histogramRDF3D(molSys);
    // Normalize the RDF 
    this->normalizeRDF3D();
}

// Updates 
void CAnalysis::histogramRDF3D(class CMolecularSystem& molSys)
{
    int natoms = this->nop; // Number of particles
    double dr;              // Relative distance between iatom and jatom (unwrapped)
    int ibin;               // Index of bin in which the particle falls wrt reference atom                              
    // Loop through every pair of particles
    for (int iatom = 0; iatom < natoms-1; iatom++)
    {
        for (int jatom = iatom+1; jatom < natoms; jatom++)
        {
            dr = this->getAbsDistance(iatom, jatom, molSys);
            // Only if dr is less than max_radius add to histogram
            if (dr < this->max_radius)
            {
                ibin = int(dr/this->binwidth); // Find which bin the particle falls in 
                this->rdf3D[ibin] += 2;        // Add to histogram for both iatom and jatom
            }
        }
    }
}

// Normalize the RDF 
// You will have to call this separately if you are averaging over
// several snapshots
void CAnalysis::normalizeRDF3D()
{
    double bin_volume;                      // Bin volume
    double nideal;                          // No. of ideal gas particles in each bin_volume
    double rho = this->nop/this->volume;    // Number density
    // Loop over all bins
    for (int ibin=0; ibin < this->nbin; ibin++)
    {
        // Volume between bin k+1 and k
        bin_volume = (pow(ibin+2, 3) - pow(ibin+1, 3)) * pow(this->binwidth, 3);
        // Assuming the total nop does not change with time 
        // Number of ideal gas particles in bin_volume
        nideal = (4.0/3.0)*PI*bin_volume*rho;
        // Normalization
        this->rdf3D[ibin] /= (this->nframes*this->nop*nideal);
    }
}

// Get the radial R values corresponding to each RDF value
// This will remain the same over all frames
void CAnalysis::getR()
{
    // Loop through all bins
    for (int ibin=0; ibin<this->nbin; ibin++)
    {
        this->rVal[ibin] = this->binwidth * (ibin + 1.5);
    }
}

// Initialize the RDF array to 0
void CAnalysis::rdf3DInitToZero()
{
    // Loop over all bins
    for (int ibin=0; ibin < this->nbin; ibin++)
    {
        this->rdf3D[ibin] = 0.0;
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
        dr[k] -= box[k]*round(dr[k]/box[k]);
        
        r2 += pow(dr[k],2.0);
    }
    

    return sqrt(r2);
}

//-------------------------------------------------------------------------------------------------------
// OUTPUT FUNCTION
//-------------------------------------------------------------------------------------------------------

// Prints out the 3D RDF function to a file in the output folder 
void CAnalysis::printRDF3D()
{
    // Prints the radial values and 3D RDF values to a file called rdf3D.txt
    this->printToFile(this->nbin, this->rVal, this->rdf3D, "rdf3D", "Radial Distance", "RDF");
} 

//-------------------------------------------------------------------------------------------------------
// CHECKS AND HELPER FUNCTIONS
//-------------------------------------------------------------------------------------------------------

// Checks that the max_radius entered is correct. If the max_radius is greater than half the simulation
// box length, by default it is set as half the smallest box length
// If the volume has not been set, set it as the simulation box volume
// TODO: Modify for 2-D case
// TODO: Check binwidth
void CAnalysis::checkParameter(class CMolecularSystem& molSys)
{
  double boxx, boxy, boxz;
  double half_box; // Half the smallest box length
  double radius = this->max_radius;
  double volume = this->volume;

  // Box lengths 
  boxx = molSys.parameter->boxx;
  boxy = molSys.parameter->boxy;
  boxz = molSys.parameter->boxz;

  half_box = 0.5*this->smallest(boxx, boxy, boxz);

  // Check if the max_radius is within bounds
  if (radius > half_box || radius <= 0.0)
  {
    std::cerr << "The maximum possible value of radius is " << half_box << " and in the parameter.txt file, it was set as " << radius << "\n";
    std::cerr << "I will now set the maximum radius to half the simulation box length " <<"\n";
    this->max_radius = half_box;
  }

  // Check if the volume has been entered. If not set it as the simulation box volume
  if (volume == -1.0)
  {
    std::cout << "Since the volume has not been set by the user, the simulation box volume will be used.\n";
    this->volume = boxx*boxy*boxz;
  }
  else if (volume <= 0.0 || volume > boxx*boxy*boxz)
  {
    std::cerr << "The volume entered cannot be used. I will use the simulation box volume instead. \n";
    this->volume = boxx*boxy*boxz;
  }
}

// Functions for returning the smallest number
double CAnalysis::smallest(double x, double y, double z)
{
  // return std::min({x,y,z}); // For C++11
  return std::min(std::min(x,y), z);
}

double CAnalysis::smallest(double x, double y)
{
  return std::min(x,y);
}
