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

/********************************************/ /**
 *  Function for Converting to GSL and getting the angle
 ***********************************************/
double gen::gslVecAngle(std::vector<double> OO, std::vector<double> OH) {
  gsl_vector *gOO = gsl_vector_alloc(3);
  gsl_vector *gOH = gsl_vector_alloc(3);
  double norm_gOO, norm_gOH, xDummy, angle;
  for (int i = 0; i < 3; i++) {
    gsl_vector_set(gOO, i, OO[i]);
    gsl_vector_set(gOH, i, OH[i]);
  }
  norm_gOO = gsl_blas_dnrm2(gOO);
  norm_gOH = gsl_blas_dnrm2(gOH);
  gsl_blas_ddot(gOO, gOH, &xDummy);
  angle = acos(xDummy / (norm_gOO * norm_gOH));
  gsl_vector_free(gOO);
  gsl_vector_free(gOH);
  return angle;
}
