#ifndef _OUTPUT_H
#define _OUTPUT_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

class COutput {
public:
  COutput();
  virtual ~COutput();
  // Print arrays to file, with the size, and two x,y arrays
  // as arguments. The names of the arrays can be printed out also, if provided
  // Otherwise default arguments are used
  void printToFile(int, double *x, double *y,
                   const std::string &filename = "output",
                   const std::string &xName = "Abscissa",
                   const std::string &yName = "Ordinate");
  //Checks whether the output directory exists, and creates it if it doesn't
  // The directory cannot be created for Windows
  void createOutputDir(const char *path);
};

#endif
