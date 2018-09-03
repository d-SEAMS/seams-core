#include "molecular_system.h"
#include "molecule.h"

// For reading in lammps traj files
const string PF_ITEM = "ITEM:";
const string PF_ATOM = "ATOMS";

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


// // Initialize the System depending on the choice of file to read from.
void CMolecularSystem::InitializeSystem()
{
  if (this->parameter->xyzFile.compare("notset") == 0 && this->parameter->trajFile.compare("notset") == 0) {
    cerr << "The filename must be set in the input file parameter.txt\n";}
  else
  {
    if (this->parameter->xyzFile.compare("notset") == 0) // Read from the traj file
    {
      this->readWholeTrj();
    } else if (this->parameter->trajFile.compare("notset") == 0) // Read from the xyz file
    {
      this->readParticleFile();
    } 
  }
}

//-------------------------------------------------------------------------------------------------------
// INITIALIZATION FROM FILE, STRUCUTURE OR FROM ARGS
//-------------------------------------------------------------------------------------------------------

//****************************************************************************************
//This procedure reads in the entire trajectory file, finds the total number 
// of steps, natoms, checks the dump file 
//****************************************************************************************
void CMolecularSystem::readWholeTrj()
{
  int nop;
  int nsteps;
  double dr[3];             // Array for boxx, boxy, boxz             
  double rlo, rhi;          // To store xlo, xhi etc from the traj file
  string::size_type pos;    // To find positions within the line
  string::size_type lpos;   // Next entry in the line
  string line;              // Current line being read in
  string word;              // To check which line is being read in 
  ifstream dumpFile;
  bool atom_flag = true;    // To test if we are reading in a coordinate line or not
                            // set to true at first
  int natoms;               // To count no. of atoms in each timestep snapshot
  dumpFile.open(this->parameter->trajFile.c_str(),ifstream::in);


  nsteps = 1;
  //First line of traj file- ITEM: TIMESTEP
  getline(dumpFile,line);
  // Contains the timestep number
  getline(dumpFile,line);
  // ITEM: NUMBER OF ATOMS
  getline(dumpFile,line); // Skip!
    
  // The next line contains the number of particles 
  getline(dumpFile, line);
  nop = atoi(line.c_str());
  this->parameter->nop = nop;
  this->initializeMolecules(nop);
    
  // Ignore line- ITEM: BOX BOUNDS pp pp pp
  getline(dumpFile,line);
  // Followed by boxx,boxy,boxz
  for (int k = 0; k <= 2; k++) 
  {
    // line contains xlo xhi separated by a space
    getline(dumpFile,line);
    pos = line.find(' ');
    rlo = strtod(line.substr(0, pos ).c_str(),NULL); 
    lpos = pos + 1;
    rhi = strtod(line.substr(lpos, line.length() - lpos).c_str(),NULL);
    dr[k] = rhi-rlo;
  } 

  // for NVT these will remain constant for all frames
  this->parameter->boxx = dr[0]; 
  this->parameter->boxy = dr[1];
  this->parameter->boxz = dr[2];

  // Skip- ITEM: ATOMS id mol type x y z 
  getline(dumpFile,line);

  // Loop through the entire file
  while( getline(dumpFile,line) )
  {
    // Find the first space in the line
    pos = line.find(' ');
    if (pos != string::npos) // There's more than one word in the line!
    {                             
      word = line.substr(0, pos); // To find the first word      
    }
    // Now get the second word if the first word is ITEM:
    if (word.compare(PF_ITEM) == 0)
    {                                 // Test the second word if ITEM:
       lpos = pos + 1;
       pos  = line.find(' ', lpos);
       if (pos != string::npos)
       {                               // Get the second word if there is another space
         word = line.substr(lpos, pos - lpos);
         if (word.compare(PF_ATOM) == 0)
         {                             // If the second word is ATOMS, then turn on atom_flag
           atom_flag = true;
         }
       }
       else                            // If there are no more spaces it is the timestep line
       {
         nsteps += 1;
         atom_flag = false;
         natoms = 0;
       }                               // End check for TIMESTEP line 
    }                                 // End test of the second word if it is ITEM 

    // For particle coordinate lines
    if (atom_flag)
    {
      natoms += 1;
    }

  } // get line


  // If the dump file has an incomplete number of atom coordinates in a snapshot
  // throw an error
  if (natoms < nop)
  {
    cerr << "Snapshot number " << nsteps << " has an incomplete list of coordinates." << "\n";
  } 

  // Finally save no. of steps
  this->parameter->nsteps = nsteps;

}


//Reads configuration from xyz file
void CMolecularSystem::readParticleFile()
{
  double posx,posy,posz;
  int nop;
  string line;
  ifstream confFile;
  confFile.open(this->parameter->xyzFile.c_str(),ifstream::in);
  if (confFile.is_open())
  {
    //the first line contains the number of particles
    getline(confFile,line);
    nop = atoi(line.c_str());
    this->parameter->nop = nop;
    this->initializeMolecules(nop);

    //Followed by boxx,boxy,boxz
    getline(confFile,line);
    this->parameter->boxx = atof(line.c_str());
    getline(confFile,line);
    this->parameter->boxy = atof(line.c_str());
    getline(confFile,line);
    this->parameter->boxz = atof(line.c_str());

    //so lets read the particles positions
    for (int ti = 0;ti<nop;ti++)
    {
      getline(confFile,line);
      string::size_type pos  = line.find(' ');
      if (pos == string::npos) break;
      posx = strtod(line.substr(0, pos ).c_str(),NULL);

      string::size_type lpos = pos + 1;
      pos  = line.find(' ', lpos);
      if (pos == string::npos) break;
      posy = strtod(line.substr(lpos, pos - lpos).c_str(),NULL);
      lpos = pos+1;

      posz = strtod(line.substr(lpos, line.length() - lpos).c_str(),NULL);

      this->molecules[ti].set_position(posx, posy, posz);
    }
  }
  else
  {
    cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
  }
}


// Reads configuration from the lammps trajectory file 
// at a particular step number
void CMolecularSystem::readParticleFile(int step)
{
  double posx,posy,posz;
  string line;              // Current line being read in
  ifstream dumpFile;
  string::size_type pos;    // To find positions within the line
  string::size_type lpos;   // Next entry in the line 
  double number;            // Each number being read from the line
  vector<double> lineVal; // Vector containing all the elements in the line
          
  dumpFile.open(this->parameter->trajFile.c_str(),ifstream::in);

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
        { // Save the coordinates in the line 
          istringstream is( line );
          // Clear the contents of the vector
          lineVal.clear();
          while (is >> number)
          {
            lineVal.push_back(number);
          }

          posx = lineVal[3];
          posy = lineVal[4];
          posz = lineVal[5];

          this->molecules[iatom].set_position(posx, posy, posz);

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


