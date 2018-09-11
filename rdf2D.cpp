#include "rdf2D.h"
#include "molecular_system.h"
#include "molecule.h"

// Constructor
Rdf2D::Rdf2D()
{
  this->binwidth = -1.0;
  this->max_radius = -1.0;
  this->volume = -1.0;
  this->nframes = 0;
  this->typeA = -1;
  this->typeB = -1;
  this->nop = -1;
}

Rdf2D::~Rdf2D()
{
  delete [] rdf2D;
  delete [] rVal;
}

/********************************************//**
 *  Initializes the histogram array for 2-D RDF 

 It takes the CMolecularSystem object reference and binwidth as arguments.
 Optional arguments include the maximum radius upto which 
 the 2-D RDF will be calculated and the desired volume. If not set, the default
 values are half the simulation box and the volume of the simulation box respectively
 ***********************************************/
void Rdf2D::initRDFxy(class CMolecularSystem& molSys, double binwidth, double volume, double max_radius)
{
    // Get the binwidth, max_radius and volume
    this->binwidth = binwidth; this->max_radius = max_radius; this->volume=volume;
    // Check the max_radius and parameters
    this->checkParameterXY(molSys);
    // Calculate the number of bins from user-defined parameters
    this->getBins();
    // Initialize the array for RDF
    this->rdf2D  = new double[this->nbin];
    // Initialize the RDF array to zero
    this->rdf2DInitToZero();
    // Get an array for R
    this->rVal = new double[this->nbin];
    // Get the values for radial distance (unchanged over nframes)
    this->getR();
}

/********************************************//**
 *  Frees the memory 
 ***********************************************/
void Rdf2D::deleteRDF2D()
{
    delete [] rdf2D;
    delete [] rVal;
}

/********************************************//**
 *  Calculates the number of bins from max_radius and binwidth
 ***********************************************/
void Rdf2D::getBins()
{
    this->nbin = int(this->max_radius/this->binwidth);
}

/********************************************//**
 *  Calculates the total number of atoms for the lammps IDs
 entered for the particular frame. The number of atoms is
 required for density calculations used for normalizing the RDF
 ***********************************************/
int Rdf2D::getNatoms(class CMolecularSystem& molSys, int typeA, int typeB)
{
  int nop=0; // No. of atoms 
  // If the lammps ID has not been set, then set nop as the total nop
  if (typeA==-1 || typeB==-1){return molSys.parameter->nop;}

  // Loop through all atoms
  for (int iatom = 0; iatom < molSys.parameter->nop; iatom++)
  {
    if (molSys.molecules[iatom].type==typeA || molSys.molecules[iatom].type==typeB){nop += 1;}
  }

  if (nop==0){
    std::cerr<<"You have entered incorrect type IDs\n"; 
    this->typeA = -1; 
    this->typeB = -1;
    return molSys.parameter->nop;
  }

  // Set the type IDs if they are correct
  this->typeA = typeA;
  this->typeB = typeB;
  return nop;
}

/********************************************//**
 *  Calculates the number of atoms for the lammps IDs
 entered for the particular frame in the XY plane defined
 by z_min and z_max. The number of atoms is for normalizing the RDF
 ***********************************************/
int Rdf2D::getNatomsXY(class CMolecularSystem& molSys, double z_min, double z_max)
{
  int nop=0; 	// No. of atoms
  double z; 	// z coordinate 
  // If the lammps ID has not been set, then set nop as the total nop
  if (typeA==-1 || typeB==-1){return molSys.parameter->nop;}

  // Loop through all atoms
  for (int iatom = 0; iatom < molSys.parameter->nop; iatom++)
  {
    if (molSys.molecules[iatom].type==typeA || molSys.molecules[iatom].type==typeB)
    	{
    		z = molSys.molecules[iatom].get_posz();
    		if (z>=z_min && z<=z_max){nop += 1;}
    	}
  }

  if (nop==0){
    std::cerr<<"There are no atoms of type"<< typeA<< "and" << typeB << "inside the z range given\n"; 
    return molSys.parameter->nop;
  }

  // Return the nop counted in the layer if greater than zero
  return nop;
}

//-------------------------------------------------------------------------------------------------------
// CALCULATIONS
//-------------------------------------------------------------------------------------------------------

/********************************************//**
 *  Calculates the 2D radial distribution function for a number of snapshots
 for a particular XY plane, specified by a range of z values

 It accepts the CMolecularSystem object reference and a pair
 of int type numbers corresponding to lammps type IDs in the trajectory
 file, and the z coordinate of the layer and width of the layer as arguments.
 The z coordinate of the layer can be chosen to be the midpoint of the peak
 in the density(or number) vs z plot. The width should be set as the width of the
 peak.
 If the integer type IDs are not set, then 
 the RDF is calculated for all the atoms in the frame, assuming they are all 
 of the same type.
 
 There is no need to use singleRDFxy() if the RDF is to be calculated over a number of frames.
 You will have to call the normalize function normalizeRDFxy() separately 
 after accumulating to get the RDF 
 ***********************************************/
void Rdf2D::accumulateRDFxy(class CMolecularSystem& molSys, double z_layer, double dz, int typeA, int typeB)
{
    // Check to make sure that the user has entered the correct type ID
    if (this->typeA!=-1 && this->typeA!=typeA && nframes>0){std::cerr<<"Type A cannot be changed after init\n";}
    if (this->typeB!=-1 && this->typeB!=typeB && nframes>0){std::cerr<<"Type B cannot be changed after init\n";}
    // Calculate the total number of particles in this particular frame
    this->nop = this->getNatoms(molSys, typeA, typeB);
    // Update the number of snapshots calculated
    this->nframes += 1;
    // Add to the RDF histogram
    this->histogramRDFxy(molSys, z_layer, dz);
    // Call the normalize function separately after accumulating 
    // the histogram over all the desired snapshots
}

/********************************************//**
 *  Calculates the 2D radial distribution function for a single snapshot for a
 particular XY plane, specified by a range of z values
 Use this only if there is one frame only.

 It accepts a pair of int type numbers corresponding to lammps type IDs in the trajectory
 file as arguments. If the integer type IDs are not set, then 
 the RDF is calculated for all the atoms in the frame, assuming they are all 
 of the same type.
 ***********************************************/
void Rdf2D::singleRDFxy(class CMolecularSystem& molSys, double z_layer, double dz, int typeA, int typeB)
{
    // There is only one snapshot
    this->nframes = 1;
    // Calculate the total number of particles in a particular frame
    this->nop = this->getNatoms(molSys, typeA, typeB);
    // Add to the RDF histogram
    this->histogramRDFxy(molSys, z_layer, dz);
    // Normalize the RDF 
    this->normalizeRDF2D(dz);
}

/********************************************//**
 *  Updates the 2D RDF histogram for a particular XY plane
 defined by z_layer and dz (layer thickness). The z_layer is the 
 z coordinate of the midpoint of the peak in the density vs z plot.
 dz is the width of the peak in the plot. 
 ***********************************************/
void Rdf2D::histogramRDFxy(class CMolecularSystem& molSys, double z_layer, double dz)
{
    int natoms = this->nop; // Number of particles
    double dr;              // Relative distance between iatom and jatom (unwrapped)
    int ibin;               // Index of bin in which the particle falls wrt reference atom                              
    int nop_layer; 			// Number of atoms in the XY plane
    double z_atom; 			// z coordinate of atom
    // Get the number of atoms in this layer
    // This is defined by z_min and z_max
    nop_layer = this->getNatomsXY(molSys, z_layer-0.5*dz, z_layer+0.5*dz);
    // Loop through every pair of particles
    for (int iatom = 0; iatom < natoms-1; iatom++)
    {
        // Only execute if the atom is of typeA
        if (molSys.molecules[iatom].type != typeA && typeA!= -1){continue;}
        z_atom = molSys.molecules[iatom].get_posz();
        if (this->atomInsideLayer(z_atom, z_layer, dz)==false){continue;}
        
        // Loop through the j^th atom
        for (int jatom = iatom+1; jatom < natoms; jatom++)
        {
            if (molSys.molecules[jatom].type != typeB && typeB!= -1){continue;}
            z_atom = molSys.molecules[iatom].get_posz();
            if (this->atomInsideLayer(z_atom, z_layer, dz)==false){continue;}
            // Test
            // double dz = molSys.molecules[jatom].get_posz() - molSys.molecules[iatom].get_posz();
            // if (abs(dz)>0.4){continue;}

            // dr = this->absDistanceXY(iatom, jatom, molSys);
            dr = this->getAbsDistance(iatom, jatom, molSys);
            // Only if dr is less than max_radius add to histogram
            if (dr < this->max_radius)
            {
                ibin = int(dr/this->binwidth); // Find which bin the particle falls in 
                // my intuition/Prerna
                this->rdf2D[ibin] += (2.0/nop_layer);        // Add to histogram for both iatom and jatom
            }
        }
    }
}

/********************************************//**
 *  Normalizes the RDF for the 2D RDF 

 You will have to call this after accumulateRDFxy() if you are averaging over
 several snapshots. This is automatically called inside singleRDFxy() and other single 
 frame RDF functions.
 This method also accepts the width of the layer as the argument
 The accuracy of the RDF calculation is sensitive to the layer width entered
 ***********************************************/

void Rdf2D::normalizeRDF2D(double dr)
{
    double bin_area;                      	// Bin area
    double nideal;                          // No. of ideal gas particles in each bin_volume
    double rho = this->nop/this->volume;    // Number density
    // Loop over all bins
    for (int ibin=0; ibin < this->nbin; ibin++)
    {
        // Area between bin k+1 and k
        bin_area = (pow(ibin+2, 2) - pow(ibin+1, 2)) * pow(this->binwidth, 2);
        // // Number of ideal gas particles in bin_volume
        nideal = PI*bin_area*dr*rho;
        // Normalization
        this->rdf2D[ibin] /= (this->nframes*nideal);
    }
}

/********************************************//**
 *  Get the radial R values corresponding to each RDF value
 This will remain the same over all frames
 ***********************************************/
void Rdf2D::getR()
{
    // Loop through all bins
    for (int ibin=0; ibin<this->nbin; ibin++)
    {
        this->rVal[ibin] = this->binwidth * (ibin + 1.5);
    }
}


/********************************************//**
 *  Initialize the RDF array to 0
 ***********************************************/
void Rdf2D::rdf2DInitToZero()
{
    // Loop over all bins
    for (int ibin=0; ibin < this->nbin; ibin++)
    {
        this->rdf2D[ibin] = 0.0;
    }
}

/********************************************//**
 *  Sets all the histogram values to zero and 
 sets the number of frames to zero so that the object can be 
 re-used. However, the binwidth and maximum radius 
 remain the same. This can be used to re-use the object 
 for a different frame etc. Use singleRDF3D() or accumulateRDF3D()
 after this function
 ***********************************************/
// TODO: Change comment
void Rdf2D::clearRDF2D()
{
    // Re-initialize the number of frames
    this->nframes = 0;

    // Loop over all bins
    for (int ibin=0; ibin < this->nbin; ibin++)
    {
        this->rdf2D[ibin] = 0.0;
    }
}

//-------------------------------------------------------------------------------------------------------
// DISTANCE CALCULATIONS
//-------------------------------------------------------------------------------------------------------

/********************************************//**
 *  Returns the absolute distance between two particles
 with particle indices iatom and jatom (x[iatom] - x[jatom])
 ***********************************************/
double Rdf2D::getAbsDistance(int iatom, int jatom, class CMolecularSystem& molSys)
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

/********************************************//**
 *  Returns the absolute distance between two particles
 with particle indices iatom and jatom (x[iatom] - x[jatom]) in 
 a specified XY plane
 ***********************************************/
double Rdf2D::absDistanceXY(int iatom, int jatom, class CMolecularSystem& molSys)
{
    double dr[2]; // Relative distance between wrapped coordinates
    double box[2] = {molSys.parameter->boxx, molSys.parameter->boxy};
    double r2 = 0.0; // Squared absolute distance

    // Get the relative distance in the x, y, z dim
    dr[0] = molSys.molecules[iatom].get_posx() - molSys.molecules[jatom].get_posx();
    dr[1] = molSys.molecules[iatom].get_posy() - molSys.molecules[jatom].get_posy();

    // Get the squared absolute distance
    for (int k=0; k<2; k++)
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

/********************************************//**
 *  Prints out the 2D RDF function to a file in the output folder
 ***********************************************/
void Rdf2D::printRDF2D()
{
    // Prints the radial values and 3D RDF values to a file called rdf3D.txt
    this->printToFile(this->nbin, this->rVal, this->rdf2D, "rdf2D", "Radial Distance", "RDF");
} 


//-------------------------------------------------------------------------------------------------------
// CHECKS AND HELPER FUNCTIONS
//-------------------------------------------------------------------------------------------------------

// Checks if the atom is inside the layer 
bool Rdf2D::atomInsideLayer(double r, double r_layer, double dr)
{
	if ((r-r_layer)<=0.5*dr){return true;}
	else {return false;}
}

/********************************************//**
 *  Checks that the max_radius entered is correct.
 If the max_radius is greater than half the simulation
 box length, by default it is set as half the smallest box length
 If the volume has not been set, set it as the simulation box volume 
 ***********************************************/
// TODO: Check binwidth
void Rdf2D::checkParameterXY(class CMolecularSystem& molSys)
{
  double boxx, boxy;
  double half_box; // Half the smallest box length
  double radius = this->max_radius;
  

  // Box lengths 
  boxx = molSys.parameter->boxx;
  boxy = molSys.parameter->boxy;

  half_box = 0.5*this->smallest(boxx, boxy);

  // Check if the max_radius is within bounds
  if (radius > half_box || radius <= 0.0)
  {
    std::cerr << "I will now set the maximum radius to half the simulation box length " <<"\n";
    this->max_radius = half_box;
  }

  // Check volume 
  this->checkVolume(molSys);

}

/********************************************//**
 *  Checks if the volume is correct.
 If the volume has not been set, set it as the simulation box volume 
 The volume will be used to calculate the density, which 
 will be used to calculate the normalization factor.
 ***********************************************/
void Rdf2D::checkVolume(class CMolecularSystem& molSys)
{
	double boxx, boxy, boxz;
	double volume = this->volume;

	// Box lengths 
  	boxx = molSys.parameter->boxx;
  	boxy = molSys.parameter->boxy;
  	boxz = molSys.parameter->boxz;
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
double Rdf2D::smallest(double x, double y, double z)
{
  // return std::min({x,y,z}); // For C++11
  return std::min(std::min(x,y), z);
}

double Rdf2D::smallest(double x, double y)
{
  return std::min(x,y);
}


//-------------------------------------------------------------------------------------------------------
// YZ PLANE SPECIFIC FUNCTIONS
//-------------------------------------------------------------------------------------------------------

/********************************************//**
 *  Initializes the histogram array for 2-D RDF in the YZ plane

 It takes the CMolecularSystem object reference and binwidth as arguments.
 Optional arguments include the maximum radius upto which 
 the 2-D RDF will be calculated and the desired volume. If not set, the default
 values are half the simulation box and the volume of the simulation box respectively
 ***********************************************/
void Rdf2D::initRDFyz(class CMolecularSystem& molSys, double binwidth, double volume, double max_radius)
{
    // Get the binwidth, max_radius and volume
    this->binwidth = binwidth; this->max_radius = max_radius; this->volume=volume;
    // Check the max_radius and parameters
    this->checkParameterYZ(molSys);
    // Calculate the number of bins from user-defined parameters
    this->getBins();
    // Initialize the array for RDF
    this->rdf2D  = new double[this->nbin];
    // Initialize the RDF array to zero
    this->rdf2DInitToZero();
    // Get an array for R
    this->rVal = new double[this->nbin];
    // Get the values for radial distance (unchanged over nframes)
    this->getR();
}

/********************************************//**
 *  Calculates the 2D radial distribution function for a single snapshot for a
 particular YZ plane, specified by a range of z values
 Use this only if there is one frame only.

 It accepts a pair of int type numbers corresponding to lammps type IDs in the trajectory
 file as arguments. If the integer type IDs are not set, then 
 the RDF is calculated for all the atoms in the frame, assuming they are all 
 of the same type.
 ***********************************************/
void Rdf2D::singleRDFyz(class CMolecularSystem& molSys, double x_layer, double dx, int typeA, int typeB)
{
    // There is only one snapshot
    this->nframes = 1;
    // Calculate the total number of particles in a particular frame
    this->nop = this->getNatoms(molSys, typeA, typeB);
    // Add to the RDF histogram
    this->histogramRDFyz(molSys, x_layer, dx);
    // Normalize the RDF 
    this->normalizeRDF2D(dx);
}

/********************************************//**
 *  Calculates the 2D radial distribution function for a number of snapshots
 for a particular YZ plane, specified by a range of x values

 It accepts the CMolecularSystem object reference and a pair
 of int type numbers corresponding to lammps type IDs in the trajectory
 file, and the x coordinate of the layer and width of the layer as arguments.
 The z coordinate of the layer can be chosen to be the midpoint of the peak
 in the density(or number) vs \f$x\f$ plot. The width should be set as the width of the
 peak.
 If the integer type IDs are not set, then 
 the RDF is calculated for all the atoms in the frame, assuming they are all 
 of the same type.
 
 There is no need to use singleRDFyz() if the RDF is to be calculated over a number of frames.
 You will have to call the normalize function normalizeRDF2D() separately 
 after accumulating to get the RDF 
 ***********************************************/
void Rdf2D::accumulateRDFyz(class CMolecularSystem& molSys, double x_layer, double dx, int typeA, int typeB)
{
    // Check to make sure that the user has entered the correct type ID
    if (this->typeA!=-1 && this->typeA!=typeA && nframes>0){std::cerr<<"Type A cannot be changed after init\n";}
    if (this->typeB!=-1 && this->typeB!=typeB && nframes>0){std::cerr<<"Type B cannot be changed after init\n";}
    // Calculate the total number of particles in this particular frame
    this->nop = this->getNatoms(molSys, typeA, typeB);
    // Update the number of snapshots calculated
    this->nframes += 1;
    // Add to the RDF histogram
    this->histogramRDFyz(molSys, x_layer, dx);
    // Call the normalize function separately after accumulating 
    // the histogram over all the desired snapshots
}

/********************************************//**
 *  Updates the 2D RDF histogram for a particular YZ plane
 defined by x_layer and dx (layer thickness). The x_layer is the 
 x coordinate of the midpoint of the peak in the density(or number) vs \f$x\f$ plot.
 dx is the width of the peak in the plot. 
 ***********************************************/
void Rdf2D::histogramRDFyz(class CMolecularSystem& molSys, double x_layer, double dx)
{
    int natoms = this->nop; // Number of particles
    double dr;              // Relative distance between iatom and jatom (unwrapped)
    int ibin;               // Index of bin in which the particle falls wrt reference atom                              
    int nop_layer;          // Number of atoms in the YZ plane
    double x_atom;          // x coordinate of atom
    // Get the number of atoms in this layer
    // This is defined by x_min and x_max
    nop_layer = this->getNatomsYZ(molSys, x_layer-0.5*dx, x_layer+0.5*dx);
    // Loop through every pair of particles
    for (int iatom = 0; iatom < natoms-1; iatom++)
    {
        // Only execute if the atom is of typeA
        if (molSys.molecules[iatom].type != typeA && typeA!= -1){continue;}
        x_atom = molSys.molecules[iatom].get_posx();
        if (this->atomInsideLayer(x_atom, x_layer, dx)==false){continue;}
        
        // Loop through the j^th atom
        for (int jatom = iatom+1; jatom < natoms; jatom++)
        {
            if (molSys.molecules[jatom].type != typeB && typeB!= -1){continue;}
            x_atom = molSys.molecules[iatom].get_posx();
            if (this->atomInsideLayer(x_atom, x_layer, dx)==false){continue;}

            dr = this->getAbsDistance(iatom, jatom, molSys);
            // Only if dr is less than max_radius add to histogram
            if (dr < this->max_radius)
            {
                ibin = int(dr/this->binwidth); // Find which bin the particle falls in 
                // my intuition/Prerna
                this->rdf2D[ibin] += (2.0/nop_layer);        // Add to histogram for both iatom and jatom
            }
        }
    }
}

/********************************************//**
 *  Calculates the number of atoms for the lammps IDs
 entered for the particular frame in the YZ plane defined
 by x_min and x_max. The number of atoms is for normalizing the RDF
 ***********************************************/
int Rdf2D::getNatomsYZ(class CMolecularSystem& molSys, double x_min, double x_max)
{
  int nop=0;  // No. of atoms
  double x;   // x coordinate 
  // If the lammps ID has not been set, then set nop as the total nop
  if (typeA==-1 || typeB==-1){return molSys.parameter->nop;}

  // Loop through all atoms
  for (int iatom = 0; iatom < molSys.parameter->nop; iatom++)
  {
    if (molSys.molecules[iatom].type==typeA || molSys.molecules[iatom].type==typeB)
      {
        x = molSys.molecules[iatom].get_posx();
        if (x>=x_min && x<=x_max){nop += 1;}
      }
  }

  if (nop==0){
    std::cerr<<"There are no atoms of type"<< typeA<< "and" << typeB << "inside the x range given\n"; 
    return molSys.parameter->nop;
  }

  // Return the nop counted in the layer if greater than zero
  return nop;
}

/********************************************//**
 *  Checks that the max_radius entered is correct.
 If the max_radius is greater than half the simulation
 box length, by default it is set as half the smallest box length
 If the volume has not been set, set it as the simulation box volume 
 ***********************************************/
// TODO: Check binwidth
void Rdf2D::checkParameterYZ(class CMolecularSystem& molSys)
{
  double boxz, boxy;
  double half_box; // Half the smallest box length
  double radius = this->max_radius;
  

  // Box lengths 
  boxy = molSys.parameter->boxy;
  boxz = molSys.parameter->boxz;

  half_box = 0.5*this->smallest(boxy, boxz);

  // Check if the max_radius is within bounds
  if (radius > half_box || radius <= 0.0)
  {
    std::cerr << "I will now set the maximum radius to half the simulation box length " <<"\n";
    this->max_radius = half_box;
  }

  // Check volume 
  this->checkVolume(molSys);

}

/********************************************//**
 *  Returns the absolute distance between two particles
 with particle indices iatom and jatom (r[iatom] - r[jatom]) in 
 a specified YZ plane
 ***********************************************/
double Rdf2D::absDistanceYZ(int iatom, int jatom, class CMolecularSystem& molSys)
{
    double dr[2]; // Relative distance between wrapped coordinates
    double box[2] = {molSys.parameter->boxy, molSys.parameter->boxz};
    double r2 = 0.0; // Squared absolute distance

    // Get the relative distance in the y, z dim
    dr[0] = molSys.molecules[iatom].get_posy() - molSys.molecules[jatom].get_posy();
    dr[1] = molSys.molecules[iatom].get_posz() - molSys.molecules[jatom].get_posz();

    // Get the squared absolute distance
    for (int k=0; k<2; k++)
    {
        // Correct for periodicity
        dr[k] -= box[k]*round(dr[k]/box[k]);
        
        r2 += pow(dr[k],2.0);
    }
    

    return sqrt(r2);
}

/********************************************//**
 *  Returns the number of bins in the RDF array
 ***********************************************/
int Rdf2D::binsInRDF()
{
    return this->nbin;
}