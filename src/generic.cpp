#include <generic.hpp>
#include <iostream>

/********************************************/ /**
 *  Function for printing out info in PairCorrel struct
 ***********************************************/
int gen::prettyPrintYoda(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::string outFile) {
  std::ofstream outputFile;
  // Create a new file in the output directory
  outputFile.open(outFile);

  if (outputFile.is_open()) {
    // First line
    outputFile << "# Frame\tAtomID\tx\ty\tz\tcij\ticeType\n";
    // Write out all the information out line by line
    for (int i = 0; i < yCloud->nop; i++) {
      outputFile << yCloud->currentFrame << "\t" << yCloud->pts[i].atomID
                 << "\t" << yCloud->pts[i].x << "\t" << yCloud->pts[i].y << "\t"
                 << yCloud->pts[i].z << "\t";
      // Print out cij
      // for(int c=0; c<yCloud->pts[i].c_ij.size(); c++){outputFile << yCloud->pts[i].c_ij[c]<<"\t";}
      // Print out the classifier
      outputFile << yCloud->pts[i].iceType << "\n";
    }
  }
  // Close the file
  outputFile.close();
  return 0;
}
