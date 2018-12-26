#include <generic.hpp>
#include <iostream>

// using namespace ;

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
 *  Function for printing out info in PairCorrel struct
 ***********************************************/
int gen::writeDump(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                   std::string outFile) {
  std::ofstream outputFile;
  // Create a new file in the output directory
  outputFile.open(outFile, std::ios_base::app);

  // Append stuff
  // -----------------------
  // Header

  // The format of the LAMMPS trajectory file is:
  // ITEM: TIMESTEP
  // 0
  // ITEM: NUMBER OF ATOMS
  // 4096
  // ITEM: BOX BOUNDS pp pp pp
  // -7.9599900000000001e-01 5.0164000000000001e+01
  // -7.9599900000000001e-01 5.0164000000000001e+01
  // -7.9599900000000001e-01 5.0164000000000001e+01
  // ITEM: ATOMS id type x y z
  // 1 1 0 0 0 etc
  outputFile << "ITEM: TIMESTEP\n";
  outputFile << yCloud->currentFrame << "\n";
  outputFile << "ITEM: NUMBER OF ATOMS\n";
  outputFile << yCloud->nop << "\n";
  outputFile << "ITEM: BOX BOUNDS pp pp pp\n";
  for (int k = 0; k < yCloud->boxLow.size(); k++) {
    outputFile << yCloud->boxLow[k] << " " << yCloud->boxLow[k] + yCloud->box[k]
               << "\n";
  } // end of printing box lengths
  outputFile << "ITEM: ATOMS id mol type x y z\n";
  // -----------------------
  // Atom lines
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    outputFile << yCloud->pts[iatom].atomID << " " << yCloud->pts[iatom].molID
               << " " << yCloud->pts[iatom].iceType << " "
               << yCloud->pts[iatom].x << " " << yCloud->pts[iatom].y << " "
               << yCloud->pts[iatom].z << "\n";
  } // end of loop through all atoms

  // Close the file
  outputFile.close();
  return 0;
}

/********************************************/ /**
 *  Function for printing out values of averaged Q6, averaged Q3 and
 Cij values
 ***********************************************/
int gen::writeHisto(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::vector<double> avgQ6) {
  std::ofstream cijFile;
  std::ofstream q3File;
  std::ofstream q6File;
  // Create a new file in the output directory
  int nNumNeighbours;
  double avgQ3;

  cijFile.open("cij.txt", std::ofstream::out | std::ofstream::app);
  q3File.open("q3.txt", std::ofstream::out | std::ofstream::app);
  q6File.open("q6.txt", std::ofstream::out | std::ofstream::app);

  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    if (yCloud->pts[iatom].type != 1) {
      continue;
    }
    // Check for slice later
    nNumNeighbours = yCloud->pts[iatom].neighList.size();
    avgQ3 = 0.0;
    for (int j = 0; j < nNumNeighbours; j++) {
      cijFile << yCloud->pts[iatom].c_ij[j].c_value << "\n";
      avgQ3 += yCloud->pts[iatom].c_ij[j].c_value;
    } // Loop through neighbours
    avgQ3 /= nNumNeighbours;
    q3File << avgQ3 << "\n";
    q6File << avgQ6[iatom] << "\n";
  } // loop through all atoms

  // Close the file
  cijFile.close();
  q3File.close();
  q6File.close();

  return 0;
}

/********************************************/ /**
 * Function to print out the largest ice cluster
 ***********************************************/
int gen::writeCluster(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                      std::string fileName, bool isSlice,
                      int largestIceCluster) {
  std::ofstream clusterFile;
  // Create a new file in the output directory
  clusterFile.open(fileName, std::ofstream::out | std::ofstream::app);
  clusterFile << yCloud->currentFrame << " " << largestIceCluster << "\n";
  // Close the file
  clusterFile.close();
  return 0;
}
