#include "parameter.h"
#include<fstream>

const string PF_NUMBEROFPARTICLES = "NumberOfParticles";
const string PF_BOXX = "xBox";
const string PF_BOXY = "yBox";
const string PF_BOXZ = "zBox";
const string PF_XYZFILE = "XYZFile";
const string PF_TRAJFILE = "trajFile";
const string PF_NEIGHBORDISTANCE = "NeighborDistance";
// For reading in lammps traj files
const string PF_ITEM = "ITEM:";
const string PF_ATOM = "ATOMS";


//Constructor
CParameter::CParameter()
{
  this->nop = -1;
  this->boxx = -1.0;
  this->boxy = -1.0;
  this->boxz = -1.0;
  this->xyzFile = "notset";

}
// Destructor
CParameter::~CParameter()
{
}

//****************************************************************************************
//This procedure reads the parameter.txt file. The keywords are defined above with PF_...
//if a line starts with // it is handled as comment
//do not have spaces before or after =
//****************************************************************************************
void CParameter::readParameter(int myrank)
{
  char IntStr[80];
  sprintf( IntStr, "input/parameter.%d.txt", myrank);
  ifstream paraFile;
  if (myrank == 0)
  {
    paraFile.open("input/parameter.txt");
  } else {
    paraFile.open(IntStr);
  }
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
    if (rawParameter[j].name.compare(PF_NUMBEROFPARTICLES) == 0) {this->nop = atoi(rawParameter[j].value.c_str());}
    if (rawParameter[j].name.compare(PF_BOXX) == 0) {this->boxx = atof(rawParameter[j].value.c_str());}
    if (rawParameter[j].name.compare(PF_BOXY) == 0) {this->boxy = atof(rawParameter[j].value.c_str());}
    if (rawParameter[j].name.compare(PF_BOXZ) == 0) {this->boxz = atof(rawParameter[j].value.c_str());}
    if (rawParameter[j].name.compare(PF_XYZFILE) == 0) {this->xyzFile = rawParameter[j].value;}
    if (rawParameter[j].name.compare(PF_TRAJFILE) == 0) {this->trajFile = rawParameter[j].value; this->readWholeTrj();}
    if (rawParameter[j].name.compare(PF_NEIGHBORDISTANCE) == 0) {this->neighbordistance = atof(rawParameter[j].value.c_str());}
  }
}

//****************************************************************************************
//This procedure reads in the entire trajectory file, finds the total number 
// of steps, natoms, checks the dump file 
//****************************************************************************************
void CParameter::readWholeTrj()
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
  dumpFile.open(this->trajFile.c_str(),ifstream::in);


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
  this->nop = nop;
    
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

  this->boxx = dr[0];
  this->boxy = dr[1];
  this->boxz = dr[2];

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
  this->nsteps = nsteps;

}

