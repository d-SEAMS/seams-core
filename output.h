#ifndef _OUTPUT_H     
#define _OUTPUT_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

class COutput {
  private:
    //Checks whether the output directory exists, and creates it if it doesn't
    // The directory cannot be created for Windows
    void createOutputDir(const char *path);
  public:
    COutput();
    virtual ~COutput();
 	// Print arrays to file, with the size, and two x,y arrays 
 	// as arguments
 	void printToFile(int, double* x, double* y, const string& filename = "output");
};

#endif