#include "molecular_system.h"
#include "molecule.h"

// For reading in lammps traj files
const std::string PF_ITEM = "ITEM:";
const std::string PF_ATOM = "ATOMS";

/********************************************//**
 *  Constructor
 ***********************************************/
CMolecularSystem::CMolecularSystem()
{
  this->parameter = new CParameter;
  this->parameter->nop = -1;
}

/********************************************//**
 *  Destructor
 ***********************************************/
CMolecularSystem::~CMolecularSystem()
{
  delete parameter;
  delete [] molecules;
}

/********************************************//**
 *  Initialize the simulation Box with number of Particles given.
 the nop will be set to the given number
 ***********************************************/
void CMolecularSystem::initializeMolecules(int numberOfParticles)
{
  this->molecules   = new CMolecule[numberOfParticles];
  this->parameter->nop = numberOfParticles;
}

/********************************************//**
 *  Assuming nop is allready set, the simulation box is initialized
 ***********************************************/
void CMolecularSystem::initializeMolecules()
{
  try {
    if (this->parameter->nop <0) {throw NumberOfParticlesNotDefinedException();}
    else {
      this->molecules   = new CMolecule[this->parameter->nop];

    }
  } catch(NumberOfParticlesNotDefinedException& ) {}

}

/********************************************//**
 *  Frees the memory
 ***********************************************/
void CMolecularSystem::deleteMolecules()
{
  delete [] molecules;
}

/********************************************//**
 *  Initialize the System depending on the choice of file to read from.
 ***********************************************/
void CMolecularSystem::InitializeSystem()
{
  if (this->parameter->xyzFile.compare("notset") == 0 && this->parameter->trajFile.compare("notset") == 0) {
    std::cerr << "The filename must be set in the input file parameter.txt\n";}
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

/********************************************//**
 *  This procedure reads in the entire trajectory file, finds the total number 
 of steps, natoms and checks the trajectory file
 ***********************************************/
void CMolecularSystem::readWholeTrj()
{
  int nop;
  int nsteps;
  double dr[3];             // Array for boxx, boxy, boxz             
  double rlo, rhi;          // To store xlo, xhi etc from the traj file
  std::string::size_type pos;    // To find positions within the line
  std::string::size_type lpos;   // Next entry in the line
  std::string line;              // Current line being read in
  std::string word;              // To check which line is being read in 
  std::ifstream dumpFile;
  bool atom_flag = true;    // To test if we are reading in a coordinate line or not
                            // set to true at first
  int natoms;               // To count no. of atoms in each timestep snapshot
  dumpFile.open(this->parameter->trajFile.c_str(), std::ifstream::in);


  nsteps = 1;
  //First line of traj file- ITEM: TIMESTEP
  std::getline(dumpFile,line);
  // Contains the timestep number
  std::getline(dumpFile,line);
  // ITEM: NUMBER OF ATOMS
  std::getline(dumpFile,line); // Skip!
    
  // The next line contains the number of particles 
  std::getline(dumpFile, line);
  nop = atoi(line.c_str());
  this->parameter->nop = nop;
  this->initializeMolecules(nop);
    
  // Ignore line- ITEM: BOX BOUNDS pp pp pp
  std::getline(dumpFile,line);
  // Followed by boxx,boxy,boxz
  for (int k = 0; k <= 2; k++) 
  {
    // line contains xlo xhi separated by a space
    std::getline(dumpFile,line);
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
  std::getline(dumpFile,line);

  // Loop through the entire file
  while( std::getline(dumpFile,line) )
  {
    // Find the first space in the line
    pos = line.find(' ');
    if (pos != std::string::npos) // There's more than one word in the line!
    {                             
      word = line.substr(0, pos); // To find the first word      
    }
    // Now get the second word if the first word is ITEM:
    if (word.compare(PF_ITEM) == 0)
    {                                 // Test the second word if ITEM:
       lpos = pos + 1;
       pos  = line.find(' ', lpos);
       if (pos != std::string::npos)
       {                               // Get the second word if there is another space
         word = line.substr(lpos, pos - lpos);
         if (word.compare(PF_ATOM) == 0)
         {                             // If the second word is ATOMS, then turn on atom_flag
           atom_flag = true;
         }
       }
       else                            // If there are no more spaces it is the timestep line
       {
         word = line.substr(lpos, std::string::npos); // Get the second word
         if (word.compare("TIMESTEP") == 0){nsteps += 1;}
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
    std::cerr << "Snapshot number " << nsteps << " has an incomplete list of coordinates." << "\n";
  } 

  // Finally save no. of steps
  this->parameter->nsteps = nsteps;

}

/********************************************//**
 *  Reads configuration from xyz file
 ***********************************************/
void CMolecularSystem::readParticleFile()
{
  double posx,posy,posz;
  int nop;
  std::string line;
  std::ifstream confFile;
  confFile.open(this->parameter->xyzFile.c_str(), std::ifstream::in);
  if (confFile.is_open())
  {
    //the first line contains the number of particles
    std::getline(confFile,line);
    nop = atoi(line.c_str());
    this->parameter->nop = nop;
    this->initializeMolecules(nop);

    //Followed by boxx,boxy,boxz
    std::getline(confFile,line);
    this->parameter->boxx = atof(line.c_str());
    std::getline(confFile,line);
    this->parameter->boxy = atof(line.c_str());
    std::getline(confFile,line);
    this->parameter->boxz = atof(line.c_str());

    //so lets read the particles positions
    for (int ti = 0;ti<nop;ti++)
    {
      std::getline(confFile,line);
      std::string::size_type pos  = line.find(' ');
      if (pos == std::string::npos) break;
      posx = strtod(line.substr(0, pos ).c_str(),NULL);

      std::string::size_type lpos = pos + 1;
      pos  = line.find(' ', lpos);
      if (pos == std::string::npos) break;
      posy = strtod(line.substr(lpos, pos - lpos).c_str(),NULL);
      lpos = pos+1;

      posz = strtod(line.substr(lpos, line.length() - lpos).c_str(),NULL);

      this->molecules[ti].set_position(posx, posy, posz);
    }
  }
  else
  {
    std::cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
  }
}

/********************************************//**
 *  Reads configuration from the lammps trajectory file 
 at a particular step number
 ***********************************************/
// TODO: Generalize for when mol ID has not been printed out
void CMolecularSystem::readParticleFile(int step)
{
  double posx,posy,posz;
  std::string line;              // Current line being read in
  std::ifstream dumpFile;
  std::string::size_type pos;    // To find positions within the line
  std::string::size_type lpos;   // Next entry in the line 
  double number;                 // Each number being read from the line
  std::vector<double> lineVal;   // Vector containing all the elements in the line
  int type;                      // Type ID of a particle
          
  dumpFile.open(this->parameter->trajFile.c_str(), std::ifstream::in);

  // Error handling for an invalid step
  if (step > this->parameter->nsteps)
  {
    std::cerr << "The step number " << step << " is larger than the number of steps in the trajectory" << "\n";
    std::cerr << "Using the first snapshot in the lammps trajectory file by default"<< "\n";
    step = 1;
  }

  if (dumpFile.is_open())
  {
    // Loop through the number of snapshots in the traj file
    // until you reach the snapshot at step
    for (int istep = 1; istep <= this->parameter->nsteps; istep++)
    {
      // Lines before coordinates in every snapshot
      std::getline(dumpFile,line); // ITEM: TIMESTEP
      std::getline(dumpFile, line); // Timestep
      std::getline(dumpFile, line); // ITEM: NUMBER OF ATOMS
      std::getline(dumpFile, line); // No. of particles; already saved though
      std::getline(dumpFile, line); // ITEM: BOX BOUNDS pp pp pp
      // Skip the three lines with box dimenions
      for (int k=0; k<3; k++)
      {
        std::getline(dumpFile, line);
      } 
      std::getline(dumpFile, line); // ITEM: ATOMS id mol type x y z 
      // Now get the particle positions; only at the correct step
      for (int iatom=0; iatom < this->parameter->nop; iatom++)
      {
        std::getline(dumpFile, line);

        // Don't save the particle positions for the rest of the snapshots
        if (istep == step)
        { // Save the coordinates in the line 
          std::istringstream is( line );
          // Clear the contents of the vector
          lineVal.clear();
          while (is >> number)
          {
            lineVal.push_back(number);
          }

          posx = lineVal[3];
          posy = lineVal[4];
          posz = lineVal[5];
          type = lineVal[2];

          this->molecules[iatom].set_position(posx, posy, posz);
          this->molecules[iatom].type = type;

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
    std::cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
  }
}


