#include "parameter.h"
#include<fstream>

const std::string PF_NUMBEROFPARTICLES = "NumberOfParticles";
const std::string PF_BOXX = "xBox";
const std::string PF_BOXY = "yBox";
const std::string PF_BOXZ = "zBox";
const std::string PF_XYZFILE = "XYZFile";
const std::string PF_TRAJFILE = "trajFile";

/********************************************//**
 *  Connstructor
 Initialize parameters like the number of particles, 
 box dimensions, xyz and traj file char strings etc
 ***********************************************/
CParameter::CParameter()
{
  this->nop = -1;
  this->boxx = -1.0;
  this->boxy = -1.0;
  this->boxz = -1.0;
  this->xyzFile = "notset";
  this->trajFile = "notset";
  this->nsteps = -1.0;
}
// Destructor
CParameter::~CParameter()
{
}

/********************************************//**
 *  This procedure reads the parameter.txt file. The keywords are defined above with PF_...
 if a line starts with // it is handled as comment
 do not have spaces before or after =
 ***********************************************/
void CParameter::readParameter()
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
      if (pos != std::string::npos)
      {
        this->rawParameter[i].name = line.substr(0, pos );
        this->rawParameter[i].value = line.substr(pos+1, std::string::npos );
        i += 1;
      } else {if (line.compare("")>0) {std::cerr << "malformed line in parameterfile :" << line << "\n";}}
    }
  }
  for (int j = 0;j < i;j++)
  {
    if (rawParameter[j].name.compare(PF_NUMBEROFPARTICLES) == 0) {this->nop = atoi(rawParameter[j].value.c_str());}
    if (rawParameter[j].name.compare(PF_BOXX) == 0) {this->boxx = atof(rawParameter[j].value.c_str());}
    if (rawParameter[j].name.compare(PF_BOXY) == 0) {this->boxy = atof(rawParameter[j].value.c_str());}
    if (rawParameter[j].name.compare(PF_BOXZ) == 0) {this->boxz = atof(rawParameter[j].value.c_str());}
    if (rawParameter[j].name.compare(PF_XYZFILE) == 0) {this->xyzFile = rawParameter[j].value;}
    if (rawParameter[j].name.compare(PF_TRAJFILE) == 0) {this->trajFile = rawParameter[j].value;}
  }
}

