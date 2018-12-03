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
  for (int t = 0; t < nop; t++) {
    // Accomodate one more point
    yCloud.pts.resize(yCloud.pts.size() + 1);
    // Match type
    if (snapshot->molecules[t].type == this->typeI) {
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
  // HAND Calculate at 0, to test
  chill::yodaPoint<double> resPointQ;
  neigh::PointCloud<double> resultCloud;
  std::array<double, 3> delta;
  std::array<double, 2> angles;
  std::vector<std::complex<double>> ylm;
  int idx;
  treeKNN::initKNN(nop, filename, 1, typeI);
  resultCloud = treeKNN::byNumber(queryIndex, 4);
  for (int itr = 0; itr < 4; itr++) {
    if (this->yCloud.pts[queryIndex].nearestID.size() < 4) {
      this->yCloud.pts[queryIndex].nearestID.resize(
          this->yCloud.pts[queryIndex].nearestID.size() + 1);
    }
    idx = resultCloud.ret_index[itr];
    this->yCloud.pts[queryIndex].nearestID[itr] = idx;
    delta[0] = this->yCloud.pts[queryIndex].x - resultCloud.pts[itr].x;
    delta[1] = this->yCloud.pts[queryIndex].y - resultCloud.pts[itr].y;
    delta[2] = this->yCloud.pts[queryIndex].z - resultCloud.pts[itr].z;
    angles = trans::radialCoord(delta);
    if (itr == 0) {
      this->yCloud.pts[queryIndex].Q = trans::spheriHarmo(3, angles);
      continue;
    }
    ylm = trans::spheriHarmo(3, angles);
    for (int m = 0; m < 7; m++) {
      this->yCloud.pts[queryIndex].Q[m] += ylm[m];
    }
  }
  for (int i = 0; i < 7; i++) {
    this->yCloud.pts[queryIndex].Q[i] =
        this->yCloud.pts[queryIndex].Q[i] * 0.25;
  }
  resPointQ = this->yCloud.pts[queryIndex];
  return resPointQ;
}

chill::yodaPoint<double> chill::bop::pointCij(int queryIndex) {
  chill::yodaPoint<double> resPointCij;
  int nearestID;
  std::complex<double> complexNumerator = {0, 0};
  std::complex<double> testI = {0, 0};
  std::complex<double> testJ = {0, 0};
  std::complex<double> complexDenominator = {0, 0};
  std::complex<double> complexDenominator1 = {0, 0};
  std::complex<double> complexDenominator2 = {0, 0};
  std::complex<double> complexDenominator3 = {0, 0};
  double realI, realJ, imI, imJ, realRes, imRes;
  realI = realJ = imI = imJ = realRes = imRes = 0;
  std::complex<double> cDum = {0, 0};
  for (int j = 0; j < 4; j++) {
    // Resize if not the right size
    if (this->yCloud.pts[queryIndex].nearestID.size() != 4) {
      this->yCloud.pts[queryIndex].nearestID.resize(4);
    }
    if (this->yCloud.pts[queryIndex].cij.size() != 4) {
      this->yCloud.pts[queryIndex].cij.resize(4);
    }
    // Use the ID of the nearest neighbor
    nearestID = this->yCloud.pts[queryIndex].nearestID[j];

    for (int m = 0; m < 7; m++) {
      testJ = std::conj(this->yCloud.pts[nearestID].Q[m]);
      testI = this->yCloud.pts[queryIndex].Q[m];
      complexNumerator = complexNumerator + (testI * testJ);
    }
    for (int m = 0; m < 7; m++) {
      testJ = this->yCloud.pts[nearestID].Q[m];
      testI = std::conj(this->yCloud.pts[nearestID].Q[m]);
      complexDenominator2 = complexDenominator2 + (testJ * testI);
    }
    for (int m = 0; m < 7; m++) {
      testJ = this->yCloud.pts[queryIndex].Q[m];
      testI = std::conj(this->yCloud.pts[queryIndex].Q[m]);
      complexDenominator3 = complexDenominator3 + (testJ * testI);
    }
    complexDenominator = sqrt(complexDenominator2 * complexDenominator3);

    cDum = complexNumerator / complexDenominator;
    std::cout << cDum << " C" << queryIndex << j << "\n";
    this->yCloud.pts[queryIndex].cij[j] = cDum.real();
  }
  resPointCij = this->yCloud.pts[queryIndex];
  return resPointCij;
}

// Assumes that you already have
chill::yodaPoint<double> chill::bop::atomVerdict(int queryIndex) {
  chill::yodaPoint<double> resPointFrame, temp;
  int nearestID = 0;
  this->yCloud.pts[queryIndex] = pointQ(queryIndex);
  for (int i = 0; i < 4; i++) {
    // Use the ID of the nearest neighbor
    nearestID = this->yCloud.pts[queryIndex].nearestID[i];
    this->yCloud.pts[nearestID] = pointQ(nearestID);
  }

  resPointFrame = pointCij(queryIndex);
  return resPointFrame;
}

void chill::bop::cleanUp() { delete snapshot; }

template class std::unique_ptr<chill::bop>;
