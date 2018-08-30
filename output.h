#ifndef _OUTPUT_H     
#define _OUTPUT_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <experimental/filesystem> // Remove experimental?

using namespace std;

class COutput {
  private:
    //The position of the particle
    void createOutputDir();
  public:
    COutput();
    virtual ~COutput();
 	// Print arrays to file
};

#endif