#include "molecular_system.h"
#include "molecule.h"


CMolecularSystem::CMolecularSystem()
{
  this->parameter = new CParameter;
  this->parameter->nop = -1;
}

CMolecularSystem::~CMolecularSystem()
{
  delete parameter;
  delete [] molecules;
}

//Initialize Simulation Box with number of Particles given
//the nop will be set to the given number
void CMolecularSystem::initializeMolecules(int numberOfParticles)
{
  this->molecules   = new CMolecule[numberOfParticles];
  this->parameter->nop = numberOfParticles;
}

//Assuming nop is allready set, the simulation box is initialized
void CMolecularSystem::initializeMolecules()
{
  try {
    if (this->parameter->nop <0) {throw NumberOfParticlesNotDefinedException();}
    else {
      this->molecules   = new CMolecule[this->parameter->nop];

    }
  } catch(NumberOfParticlesNotDefinedException& ) {}

}

//Frees the memory
void CMolecularSystem::deleteMolecules()
{
  delete [] molecules;
}


// // Initialize the System depending on the choice of xyzFile.
// void CMolecularSystem::InitializeSystem()
// {

// }

//-------------------------------------------------------------------------------------------------------
// INITIALIZATION FROM FILE, STRUCUTURE OR FROM ARGS
//-------------------------------------------------------------------------------------------------------

//Reads configuration from xyz file
// void CMolecularSystem::readParticleFile()
// {
//   double posx,posy,posz;
//   int nop;
//   string line;
//   ifstream confFile;
//   confFile.open(this->parameter->xyzFile.c_str(),ifstream::in);
//   if (confFile.is_open())
//   {
//     //the first line contains the number of particles
//     getline(confFile,line);
//     nop = atoi(line.c_str());
//     this->parameter->nop = nop;
//     this->initializeMolecules(nop);

//     //Followed by boxx,boxy,boxz
//     getline(confFile,line);
//     this->parameter->boxx = atof(line.c_str());
//     getline(confFile,line);
//     this->parameter->boxy = atof(line.c_str());
//     getline(confFile,line);
//     this->parameter->boxz = atof(line.c_str());

//     //so lets read the particles positions
//     for (int ti = 0;ti<nop;ti++)
//     {
//       getline(confFile,line);
//       string::size_type pos  = line.find(' ');
//       if (pos == string::npos) break;
//       posx = strtod(line.substr(0, pos ).c_str(),NULL);

//       string::size_type lpos = pos + 1;
//       pos  = line.find(' ', lpos);
//       if (pos == string::npos) break;
//       posy = strtod(line.substr(lpos, pos - lpos).c_str(),NULL);
//       lpos = pos+1;

//       posz = strtod(line.substr(lpos, line.length() - lpos).c_str(),NULL);

//       this->molecules[ti].posx = posx;
//       this->molecules[ti].posy = posy;
//       this->molecules[ti].posz = posz;
//     }
//   }
//   else
//   {
//     cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
//   }
// }


// Reads configuration from the lammps trajectory file 
// at a particular step number
void CMolecularSystem::readParticleFile(int step)
{
  double posx,posy,posz;
  string line;              // Current line being read in
  ifstream dumpFile;
  string::size_type pos;    // To find positions within the line
  string::size_type lpos;   // Next entry in the line             
  dumpFile.open(this->parameter->trajFile.c_str(),ifstream::in);
  // nop has already been saved 
  this->initializeMolecules(this->parameter->nop);

  // Error handling for an invalid step
  if (step > this->parameter->nsteps)
  {
    cerr << "The step number " << step << " is larger than the number of steps in the trajectory" << "\n";
    cerr << "Using the first snapshot in the lammps trajectory file by default"<< "\n";
    step = 1;
  }

  if (dumpFile.is_open())
  {
    // Loop through the number of snapshots in the traj file
    // until you reach the snapshot at step
    for (int istep = 1; istep <= this->parameter->nsteps; istep++)
    {
      // Lines before coordinates in every snapshot
      getline(dumpFile,line); // ITEM: TIMESTEP
      getline(dumpFile, line); // Timestep
      getline(dumpFile, line); // ITEM: NUMBER OF ATOMS
      getline(dumpFile, line); // No. of particles; already saved though
      getline(dumpFile, line); // ITEM: BOX BOUNDS pp pp pp
      // Skip the three lines with box dimenions
      for (int k=0; k<3; k++)
      {
        getline(dumpFile, line);
      } 
      getline(dumpFile, line); // ITEM: ATOMS id mol type x y z 
      // Now get the particle positions; only at the correct step
      for (int iatom=0; iatom < this->parameter->nop; iatom++)
      {
        getline(dumpFile, line);

        // Don't save the particle positions for the rest of the snapshots
        if (istep == step)
        { // Save the coordinates
          // Find the third space in the line
          pos = line.find(' '); // first space
          lpos = pos + 1;
          for (int k = 0; k < 2; k++) // second and third space
          {
            pos = line.find(' ', lpos);
            lpos = pos + 1; 
          }
          
          if (pos == string::npos) break;
          posx = strtod(line.substr(lpos, pos - lpos).c_str(),NULL);

          lpos = pos + 1;
          pos  = line.find(' ', lpos);
          if (pos == string::npos) break;
          posy = strtod(line.substr(lpos, pos - lpos).c_str(),NULL);
          // cout << posy << " ";
          lpos = pos+1;

          posz = strtod(line.substr(lpos, line.length() - lpos).c_str(),NULL);
          // cout << posz << "\n";

          this->molecules[iatom].posx = posx;
          this->molecules[iatom].posy = posy;
          this->molecules[iatom].posz = posz;

        }
        else // Skip lines for other snapshots
        {
          continue;
        }

      }   // End of loop through coordinate lines
    }     // End of loop through snapshots in the dump file
  
  }       // End of check for the file being open
  else
  {
    cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
  }
}

//-------------------------------------------------------------------------------------------------------
// HELPERS
//-------------------------------------------------------------------------------------------------------

//converts a int into a string
void CMolecularSystem::IntToString(int i, string& res)
{
    ostringstream tempel;
    tempel << i;
    res = tempel.str();
}

//-------------------------------------------------------------------------------------------------------
//DISTANCE CACULATIONS
//-------------------------------------------------------------------------------------------------------

void CMolecularSystem::makeperiodic(int ti)
{ 
  double boxx,boxy,boxz;
  boxx = this->parameter->boxx;
  boxy = this->parameter->boxy;
  boxz = this->parameter->boxz;
  
  if (this->molecules[ti].posx >=  boxx) {this->molecules[ti].posx -=boxx;}
  if (this->molecules[ti].posx <   0.0)  {this->molecules[ti].posx +=boxx;}
  if (this->molecules[ti].posy >=  boxy) {this->molecules[ti].posy -=boxy;}
  if (this->molecules[ti].posy <   0.0)  {this->molecules[ti].posy +=boxy;}
  if (this->molecules[ti].posz >=  boxz) {this->molecules[ti].posz -=boxz;}
  if (this->molecules[ti].posz <   0.0)  {this->molecules[ti].posz +=boxz;}

}


//Distance with nearest image convention
void CMolecularSystem::get_distance(int ti ,int tj ,double &diffx ,double &diffy,double &diffz)
{
  diffx = this->molecules[tj].posx - this->molecules[ti].posx;
  diffy = this->molecules[tj].posy - this->molecules[ti].posy;
  diffz = this->molecules[tj].posz - this->molecules[ti].posz;

  //nearest image
  if (diffx >  this->parameter->boxx/2.0) {diffx = diffx - this->parameter->boxx;};
  if (diffx < -this->parameter->boxx/2.0) {diffx = diffx + this->parameter->boxx;};
  if (diffy >  this->parameter->boxy/2.0) {diffy = diffy - this->parameter->boxy;};
  if (diffy < -this->parameter->boxy/2.0) {diffy = diffy + this->parameter->boxy;};
  if (diffz >  this->parameter->boxz/2.0) {diffz = diffz - this->parameter->boxz;};
  if (diffz < -this->parameter->boxz/2.0) {diffz = diffz + this->parameter->boxz;};
}


void CMolecularSystem::get_distancePosition(int ti ,double posx, double posy,double posz ,double &diffx,double &diffy,double &diffz)
{
  diffx = posx - this->molecules[ti].posx;
  diffy = posy - this->molecules[ti].posy;
  diffz = posz - this->molecules[ti].posz;

  //nearest image
  if (diffx >  this->parameter->boxx/2.0) {diffx = diffx - this->parameter->boxx;};
  if (diffx < -this->parameter->boxx/2.0) {diffx = diffx + this->parameter->boxx;};
  if (diffy >  this->parameter->boxy/2.0) {diffy = diffy - this->parameter->boxy;};
  if (diffy < -this->parameter->boxy/2.0) {diffy = diffy + this->parameter->boxy;};
  if (diffz >  this->parameter->boxz/2.0) {diffz = diffz - this->parameter->boxz;};
  if (diffz < -this->parameter->boxz/2.0) {diffz = diffz + this->parameter->boxz;};
}
    

double CMolecularSystem::get_absDistance(int ti ,int tj)
{
  double abs,diffx,diffy,diffz;
  diffx = this->molecules[tj].posx - this->molecules[ti].posx;
  diffy = this->molecules[tj].posy - this->molecules[ti].posy;
  diffz = this->molecules[tj].posz - this->molecules[ti].posz;
  //nearest image
  if (diffx >  this->parameter->boxx/2.0) {diffx = diffx - this->parameter->boxx;};
  if (diffx < -this->parameter->boxx/2.0) {diffx = diffx + this->parameter->boxx;};
  if (diffy >  this->parameter->boxy/2.0) {diffy = diffy - this->parameter->boxy;};
  if (diffy < -this->parameter->boxy/2.0) {diffy = diffy + this->parameter->boxy;};
  if (diffz >  this->parameter->boxz/2.0) {diffz = diffz - this->parameter->boxz;};
  if (diffz < -this->parameter->boxz/2.0) {diffz = diffz + this->parameter->boxz;};
  abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
  return abs;
}

//-------------------------------------------------------------------------------------------------------
//NEIGHBOURHOOD LIST
//-------------------------------------------------------------------------------------------------------

//Get all nearest neighbors within this->neighbordistance
//mainly used for the bond orders
void CMolecularSystem::get_AllNeighbors()
{
  double nd,d;
  nd = this->parameter->neighbordistance;
  int c,nop;
  nop = this->parameter->nop;

  for (int ti = 0;ti<nop;ti++)
  {
    for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++)
    {
      this->molecules[ti].neighbors[tn] = nilvalue;
      this->molecules[ti].neighbordist[tn] = -1.0;
    }
  }


  for (int ti = 0;ti<nop;ti++)
  {
    c = 0;
    for (int tj = 0;tj<nop;tj++)
    {

          if (tj != ti) { 
            d = get_absDistance(ti,tj); 
            if (d < nd) {
              this->molecules[ti].neighbors[c] = tj; 
              this->molecules[ti].neighbordist[c]=d; 
              c += 1;   
            }
         }
    }
    this->molecules[ti].n_neighbors = c;
    cout << ti << " " << c << "\n";
  }

}


//-------------------------------------------------------------------------------------------------------
//OUTPUTSECTION
//-------------------------------------------------------------------------------------------------------


















