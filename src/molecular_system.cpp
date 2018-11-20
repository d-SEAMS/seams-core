#include "molecular_system.h"
#include "molecule.h"
#include <bits/stdc++.h>

// For reading in lammps traj files
const std::string PF_ITEM = "ITEM:";
const std::string PF_ATOM = "ATOMS";

/********************************************/ /**
 *  Constructor
 ***********************************************/
CMolecularSystem::CMolecularSystem() {
  this->parameter = new CParameter;
  this->parameter->nop = -1;
  this->molecules = nullptr;
}

/********************************************/ /**
 *  Destructor
 ***********************************************/
CMolecularSystem::~CMolecularSystem() {
  delete parameter;
  delete[] molecules;
}

/********************************************/ /**
 *  Prepares frames for further processing.
 (Used by the transition module)
 ***********************************************/

void CMolecularSystem::initializeFrames(int nop, std::string fileName) {
  this->molecules = new CMolecule[nop];
  this->parameter->nop = nop;
  this->parameter->trajFile = fileName;
}

/********************************************/ /**
 *  Initialize the simulation Box with number of Particles given.
 the nop will be set to the given number
 ***********************************************/
void CMolecularSystem::initializeMolecules(int numberOfParticles) {
  this->molecules = new CMolecule[numberOfParticles];
  this->parameter->nop = numberOfParticles;
}

/********************************************/ /**
 *  Assuming nop is allready set, the simulation box is initialized
 ***********************************************/
void CMolecularSystem::initializeMolecules() {
  try {
    if (this->parameter->nop < 0) {
      throw NumberOfParticlesNotDefinedException();
    } else {
      this->molecules = new CMolecule[this->parameter->nop];
    }
  } catch (NumberOfParticlesNotDefinedException &) {
  }
}

/********************************************/ /**
 *  Frees the memory
 ***********************************************/
void CMolecularSystem::deleteMolecules() { delete[] molecules; }

/********************************************/ /**
 *  Initialize the System depending on the choice of file to read from.
 ***********************************************/
void CMolecularSystem::InitializeSystem() {
  if (this->parameter->xyzFile.compare("notset") == 0 &&
      this->parameter->trajFile.compare("notset") == 0) {
    std::cerr << "The filename must be set in the input file parameter.txt\n";
  } else {
    if (this->parameter->xyzFile.compare("notset") ==
        0) // Read from the traj file
    {
      this->readWholeTrj();
    } else if (this->parameter->trajFile.compare("notset") ==
               0) // Read from the xyz file
    {
      this->readParticleFile();
    }
  }
}

//-------------------------------------------------------------------------------------------------------
// INITIALIZATION FROM FILE, STRUCUTURE OR FROM ARGS
//-------------------------------------------------------------------------------------------------------

/********************************************/ /**
 *  This procedure reads in the entire trajectory file, finds the total number 
 of steps, natoms and checks the trajectory file
 ***********************************************/
void CMolecularSystem::readWholeTrj() {
  int nop;
  int nsteps;
  double dr[3];                // Array for boxx, boxy, boxz
  double rlo, rhi;             // To store xlo, xhi etc from the traj file
  std::string::size_type pos;  // To find positions within the line
  std::string::size_type lpos; // Next entry in the line
  std::string line;            // Current line being read in
  std::string word;            // To check which line is being read in
  std::ifstream dumpFile;
  bool atom_flag =
      true;   // To test if we are reading in a coordinate line or not
              // set to true at first
  int natoms; // To count no. of atoms in each timestep snapshot
  dumpFile.open(this->parameter->trajFile.c_str(), std::ifstream::in);

  nsteps = 1;
  //First line of traj file- ITEM: TIMESTEP
  std::getline(dumpFile, line);
  // Contains the timestep number
  std::getline(dumpFile, line);
  // ITEM: NUMBER OF ATOMS
  std::getline(dumpFile, line); // Skip!

  // The next line contains the number of particles
  std::getline(dumpFile, line);
  nop = atoi(line.c_str());
  this->parameter->nop = nop;
  this->initializeMolecules(nop);

  // Ignore line- ITEM: BOX BOUNDS pp pp pp
  std::getline(dumpFile, line);
  // Followed by boxx,boxy,boxz
  for (int k = 0; k <= 2; k++) {
    // line contains xlo xhi separated by a space
    std::getline(dumpFile, line);
    pos = line.find(' ');
    rlo = strtod(line.substr(0, pos).c_str(), nullptr);
    lpos = pos + 1;
    rhi = strtod(line.substr(lpos, line.length() - lpos).c_str(), nullptr);
    dr[k] = rhi - rlo;
  }

  // for NVT these will remain constant for all frames
  this->parameter->boxx = dr[0];
  this->parameter->boxy = dr[1];
  this->parameter->boxz = dr[2];

  // Skip- ITEM: ATOMS id mol type x y z
  std::getline(dumpFile, line);

  // Loop through the entire file
  while (std::getline(dumpFile, line)) {
    // Find the first space in the line
    pos = line.find(' ');
    if (pos != std::string::npos) // There's more than one word in the line!
    {
      word = line.substr(0, pos); // To find the first word
    }
    // Now get the second word if the first word is ITEM:
    if (word.compare(PF_ITEM) == 0) { // Test the second word if ITEM:
      lpos = pos + 1;
      pos = line.find(' ', lpos);
      if (pos !=
          std::string::npos) { // Get the second word if there is another space
        word = line.substr(lpos, pos - lpos);
        if (word.compare(PF_ATOM) ==
            0) { // If the second word is ATOMS, then turn on atom_flag
          atom_flag = true;
        }
      } else // If there are no more spaces it is the timestep line
      {
        word = line.substr(lpos, std::string::npos); // Get the second word
        if (word.compare("TIMESTEP") == 0) {
          nsteps += 1;
        }
        atom_flag = false;
        natoms = 0;
      } // End check for TIMESTEP line
    }   // End test of the second word if it is ITEM

    // For particle coordinate lines
    if (atom_flag) {
      natoms += 1;
    }

  } // get line

  // If the dump file has an incomplete number of atom coordinates in a snapshot
  // throw an error
  if (natoms < nop) {
    std::cerr << "Snapshot number " << nsteps
              << " has an incomplete list of coordinates."
              << "\n";
  }

  // Finally save no. of steps
  this->parameter->nsteps = nsteps;
}

/********************************************/ /**
 *  Reads configuration from xyz file
 ***********************************************/
void CMolecularSystem::readParticleFile() {
  double posx, posy, posz;
  int nop;
  std::string line;
  std::ifstream confFile;
  confFile.open(this->parameter->xyzFile.c_str(), std::ifstream::in);
  if (confFile.is_open()) {
    //the first line contains the number of particles
    std::getline(confFile, line);
    nop = atoi(line.c_str());
    this->parameter->nop = nop;
    this->initializeMolecules(nop);

    //Followed by boxx,boxy,boxz
    std::getline(confFile, line);
    this->parameter->boxx = atof(line.c_str());
    std::getline(confFile, line);
    this->parameter->boxy = atof(line.c_str());
    std::getline(confFile, line);
    this->parameter->boxz = atof(line.c_str());

    //so lets read the particles positions
    for (int ti = 0; ti < nop; ti++) {
      std::getline(confFile, line);
      std::string::size_type pos = line.find(' ');
      if (pos == std::string::npos)
        break;
      posx = strtod(line.substr(0, pos).c_str(), nullptr);

      std::string::size_type lpos = pos + 1;
      pos = line.find(' ', lpos);
      if (pos == std::string::npos)
        break;
      posy = strtod(line.substr(lpos, pos - lpos).c_str(), nullptr);
      lpos = pos + 1;

      posz = strtod(line.substr(lpos, line.length() - lpos).c_str(), nullptr);

      this->molecules[ti].set_position(posx, posy, posz);
    }
  } else {
    std::cerr << "Fatal Error : cannot open the file "
              << this->parameter->xyzFile << "\n";
  }
}

/********************************************/ /**
 *  Reads configuration from the lammps trajectory file 
 at a particular step number
 ***********************************************/
void CMolecularSystem::readParticleFile(int step) {
  double posx, posy, posz;
  std::string line; // Current line being read in
  std::ifstream dumpFile;
  double number;               // Each number being read from the line
  std::vector<double> lineVal; // Vector containing all the elements in the line
  int type;                    // Type ID of a particle
  int molID;                   // Molecule ID of each particle
  int atomID;                  // Stores lammps atom ID
  std::vector<std::string>
      lammpsLine; // Line that contains info about column numbers for type, ID etc.
  bool molFlag =
      false; // Flag that checks if the molecular ID has been entered or not
  std::string word;   // To store individual words
  int wordNumber = 0; // To store the word number
  int molNum;
  int typeNum;
  int xNum;

  dumpFile.open(this->parameter->trajFile.c_str(), std::ifstream::in);

  // Error handling for an invalid step
  // TODO: Do this better, maybe use <optional> wrt. https://stackoverflow.com/a/47677892/1895378
  if (this->parameter->nsteps > 1) {
    if (step > this->parameter->nsteps) {
      std::cerr << "The step number " << step
                << " is larger than the number of steps in the trajectory"
                << "\n";
      std::cerr
          << "Using the first snapshot in the lammps trajectory file by default"
          << "\n";
      step = 1;
    }
  }

  if (dumpFile.is_open()) {
    // Loop through the number of snapshots in the traj file
    // until you reach the snapshot at step
    for (int istep = 1; istep <= this->parameter->nsteps; istep++) {
      // Lines before coordinates in every snapshot
      std::getline(dumpFile, line); // ITEM: TIMESTEP
      std::getline(dumpFile, line); // Timestep
      std::getline(dumpFile, line); // ITEM: NUMBER OF ATOMS
      std::getline(dumpFile, line); // No. of particles; already saved though
      std::getline(dumpFile, line); // ITEM: BOX BOUNDS pp pp pp

      // Get the box lengths from the three lines with box dimenions
      for (int k = 0; k < 3; k++) {
        std::getline(dumpFile, line);
        getBoxLength(line);
      }

      // -----------------------
      std::getline(dumpFile, line); // ITEM: ATOMS id mol type x y z

      // Do this only for the first step!
      if (istep == 1) {
        // Check what the column number of the x,y,z coordinates are
        // Find out if the molecular ID has been entered!
        // breaking line into word using string stream
        std::stringstream ss(line); // Used for breaking words from line
        wordNumber = 0;
        while (ss >> word) {
          wordNumber++;
          if (word == "mol") {
            molFlag = true;
            molNum = wordNumber - 3;
          }
          if (word == "type") {
            typeNum = wordNumber - 3;
          }
          if (word == "x") {
            xNum = wordNumber - 3;
          }
        }
      }

      // ------------------------
      // Now get the particle positions; only at the correct step
      for (int iatom = 0; iatom < this->parameter->nop; iatom++) {
        std::getline(dumpFile, line);

        // Don't save the particle positions for the rest of the snapshots
        if (istep == step) { // Save the coordinates in the line
          std::istringstream is(line);
          // Clear the contents of the vector
          lineVal.clear();
          while (is >> number) {
            lineVal.push_back(number);
          }

          posx = lineVal[xNum];
          posy = lineVal[xNum + 1];
          posz = lineVal[xNum + 2];
          type = lineVal[typeNum];
          atomID = lineVal[0];
          if (molFlag == true) {
            molID = lineVal[molNum];
          } else {
            molID = iatom;
          }

          this->molecules[iatom].set_position(posx, posy, posz);
          this->molecules[iatom].type = type;
          this->molecules[iatom].molID = molID;
          this->molecules[iatom].atomID = atomID;

        } else // Skip lines for other snapshots
        {
          continue;
        }

      } // End of loop through coordinate lines
    }   // End of loop through snapshots in the dump file

  } // End of check for the file being open
  else {
    std::cerr << "Fatal Error : cannot open the file "
              << this->parameter->xyzFile << "\n";
  }
}

/********************************************/ /**
 *  Reads in and updates the box lengths at each step
 (in case the box volume is changing at each step; eg for NPT)
 ***********************************************/
// TODO: Write different logic if you want the exact box dimensions and
// not just the box lengths
double CMolecularSystem::getBoxLength(std::string line) {
  std::vector<double> lineVal; // Vector containing all the elements in the line
  double number;               // Each number being read from the line

  std::stringstream is(line); // Used for breaking words from line
  // Clear the contents of the vector
  lineVal.clear();
  while (is >> number) {
    lineVal.push_back(number);
  }

  // Get the box length
  return lineVal[1] - lineVal[0];
}
