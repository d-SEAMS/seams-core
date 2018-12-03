// For calculating the "correlation of orientational bond order parameter
// developed by ten Wolde et al. based on the Steinhardt parameters"
// Essentially the CHILL+ params by Molinero et al.
// This is a library.
#include <bop_correl.hpp>
#include <spherical_harmonics.h>

/********************************************/ /**
 *  Public Functions
 ***********************************************/

// void chill::test::sim() { std::cout << "sim"; };
int chill::bop::initBOP(int nop, int typeI, chill::initSlice<double> starter) {
  std::cout << "I am bop\n";
  // Set private variables
  this->nop = nop;
  this->typeI = typeI;
  this->filename = starter.filename;
  // Prepares the frame
  this->prepSnapshot(starter);
  // Use nsteps (dummy)
  this->snapshot->parameter->nsteps = starter.frameRange[1] + 1;
  // Read the file
  this->snapshot->readParticleFile(starter.frameRange[0]);
  // Now we populate the cloud (make int and do error handling)
  this->populateSnapshot(starter);
  return 1;
}

void chill::bop::prepSnapshot(chill::initSlice<double> starter) {
  this->snapshot->initializeFrames(nop, starter.filename);
}
void chill::bop::populateSnapshot(chill::initSlice<double> starter) {
  int filteredParticles(0);
  int nop(snapshot->parameter->nop);
  int dummy(0);
  // Do filtering
  for (int t; t < nop; t++) {
    // Accomodate one more point
    yCloud.pts.resize(yCloud.pts.size() + 1);
    // Match type
    if (snapshot->molecules[t].type == typeI) {
      // Check limits
      if (treeKNN::isThere(t, snapshot, starter.coordHigh, starter.coordLow)) {
        // Dump point
        yCloud.pts[t].x = snapshot->molecules[t].get_posx();
        yCloud.pts[t].y = snapshot->molecules[t].get_posy();
        yCloud.pts[t].z = snapshot->molecules[t].get_posz();
        yCloud.pts[t].inSlice = true;
      }
    }
  }
}

chill::yodaPoint<double> chill::bop::pointQ(int queryIndex) {
  chill::yodaPoint<double> resPoint;
  neigh::PointCloud<double> resultCloud;
  std::array<double, 3> delta;
  std::array<double, 2> angles;
  std::vector<std::complex<double>> ylm;
  int idx;
  treeKNN::initKNN(nop, filename, 1, typeI);
  resultCloud = treeKNN::byNumber(0, 4);
  for (int itr = 0; itr < 4; itr++) {
    idx = resultCloud.ret_index[itr];
    yCloud.pts[queryIndex].nearestID[itr] = idx;
    delta[0] = yCloud.pts[queryIndex].x - resultCloud.pts[idx].x;
    delta[1] = yCloud.pts[queryIndex].y - resultCloud.pts[idx].y;
    delta[2] = yCloud.pts[queryIndex].z - resultCloud.pts[idx].z;
    angles = trans::radialCoord(delta);
    if (itr == 0) {
      yCloud.pts[queryIndex].Q = trans::spheriHarmo(3, angles);
      continue;
    }
    ylm = trans::spheriHarmo(3, angles);
    for (int m = 0; m < 7; m++) {
      yCloud.pts[queryIndex].Q[m] += ylm[m];
    }
  }
  for (int i; i < 7; i++) {
    yCloud.pts[queryIndex].Q[i] = yCloud.pts[queryIndex].Q[i] / (double)4.0;
  }
  resPoint = yCloud.pts[queryIndex];
  return resPoint;
}
