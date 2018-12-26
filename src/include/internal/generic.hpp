#ifndef __GENERIC_H_
#define __GENERIC_H_

#include <array>
#include <math.h>
#include <mol_sys.hpp>

namespace gen {

// Generic function for getting the unwrapped distance
inline double
periodicDist(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
             int iatom, int jatom) {
  std::array<double, 3> dr;
  double r2 = 0.0; // Squared absolute distance

  // Get x1-x2 etc
  dr[0] = fabs(yCloud->pts[iatom].x - yCloud->pts[jatom].x);
  dr[1] = fabs(yCloud->pts[iatom].y - yCloud->pts[jatom].y);
  dr[2] = fabs(yCloud->pts[iatom].z - yCloud->pts[jatom].z);

  // Get the squared absolute distance
  for (int k = 0; k < 3; k++) {
    // Correct for periodicity
    dr[k] -= yCloud->box[k] * round(dr[k] / yCloud->box[k]);
    r2 += pow(dr[k], 2.0);
  }

  return sqrt(r2);
}

// Generic function for getting the relative coordinates
inline std::array<double, 3>
relDist(molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatom,
        int jatom) {
  std::array<double, 3> dr;
  std::array<double, 3> box = {yCloud->box[0], yCloud->box[1], yCloud->box[2]};
  double r2 = 0.0; // Squared absolute distance

  // Get x1-x2 etc
  dr[0] = yCloud->pts[iatom].x - yCloud->pts[jatom].x;
  dr[1] = yCloud->pts[iatom].y - yCloud->pts[jatom].y;
  dr[2] = yCloud->pts[iatom].z - yCloud->pts[jatom].z;

  // Get the relative distance
  for (int k = 0; k < 3; k++) {
    //
    if (dr[k] < -box[k] * 0.5) {
      dr[k] = dr[k] + box[k];
    }
    if (dr[k] >= box[k] * 0.5) {
      dr[k] = dr[k] - box[k];
    }
  }

  return dr;
}

// Function for sorting according to atom ID
// Comparator for std::sort
inline bool compareByAtomID(const molSys::Point<double> &a,
                            const molSys::Point<double> &b) {
  return a.atomID < b.atomID;
}

// Generic function for printing all the struct information
int prettyPrintYoda(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::string outFile);

// Generic function for writing out to a dump file
int writeDump(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
              std::string outFile);

// Function for printing out Q6, Cij and averaged Q3 values as single columns to text files
// The file names are cij, q6, q3
int writeHisto(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
               std::vector<double> avgQ6);
// Function for printing the largest ice cluster
int writeCluster(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                 std::string fileName = "cluster.txt", bool isSlice = false,
                 int largestIceCluster = 0);
} // namespace gen

#endif // __NEIGHBOURS_H_
