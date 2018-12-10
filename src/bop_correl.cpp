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
  this->frame = starter.frameRange[0];
  // Prepares the frame
  this->prepSnapshot(starter);
  // Use nsteps (dummy)
  this->snapshot->parameter->nsteps = starter.frameRange[1] + 1;
  // Read the file
  this->snapshot->readParticleFile(frame);
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

    // Roma version of Chill
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
  std::array<int, 4> bondIs = {-1, -1, -1, -1};
  int isEclipsed, isStaggered;
  int nearestID = 0;
  this->yCloud.pts[queryIndex] = pointQ(queryIndex);
  for (int i = 0; i < 4; i++) {
    // Use the ID of the nearest neighbor
    nearestID = this->yCloud.pts[queryIndex].nearestID[i];
    // Don't re-run this
    if (this->yCloud.pts[nearestID].Q.size() != 7) {
      this->yCloud.pts[nearestID] = pointQ(nearestID);
    }
  }

  resPointFrame = pointCij(queryIndex);

  // For chill+ (Molinero et. al.)
  isEclipsed = isStaggered = 0;
  for (int i = 0; i < 4; i++) {
    if (resPointFrame.cij[i] <= 0.25 && resPointFrame.cij[i] > -0.35) {
      bondIs[i] = 1; // Eclipsed
      isEclipsed++;
    } else if (resPointFrame.cij[i] <= -0.35 && resPointFrame.cij[i] >= -1) {
      bondIs[i] = 0; // Staggered
      isStaggered++;
    }
  }

  if (isEclipsed == 0 && isStaggered == 4) {
    resPointFrame.chillPlus.isCubic = true;
  } else if (isEclipsed == 1 && isStaggered == 3) {
    resPointFrame.chillPlus.isHexa = true;
  } else if (isStaggered == 2 && isEclipsed == 2) {
    resPointFrame.chillPlus.isInterfacial = true;
  } else if (isEclipsed == 3) {
    resPointFrame.chillPlus.isInterClathrate = true;
  } else if (isEclipsed == 4 && isStaggered == 0) {
    resPointFrame.chillPlus.isClathrate = true;
  } else {
    resPointFrame.chillPlus.isUndef = true;
  }

  // Update globally
  this->yCloud.pts[queryIndex] = resPointFrame;
  // Return locally
  return resPointFrame;
}

// Figure out relative percentages
chill::structurePercentage chill::bop::frameVerdict(int currentFrame) {
  chill::structurePercentage resultPercent;
  chill::yodaPoint<double> temp;
  int totalParticles = 0;
  for (int i = 0; i < 10; i++) {
    // Generate the raw values
    temp = atomVerdict(i);
    if (temp.chillPlus.isHexa) {
      this->yCloud.framePercent.hexa++;
    } else if (temp.chillPlus.isCubic) {
      this->yCloud.framePercent.cubic++;
    } else if (temp.chillPlus.isInterfacial) {
      this->yCloud.framePercent.interfacial++;
    } else if (temp.chillPlus.isClathrate) {
      this->yCloud.framePercent.clathrate++;
    } else if (temp.chillPlus.isWater) {
      this->yCloud.framePercent.water++;
    } else if (temp.chillPlus.isInterClathrate) {
      this->yCloud.framePercent.interClathrate++;
    } else { // isUndef is implied
      this->yCloud.framePercent.undef++;
    }
    // Track the number of iterations (for slicing)
    totalParticles++;
  }
  // Convert to percentage (local only)
  resultPercent.hexa = (this->yCloud.framePercent.hexa / totalParticles) * 100;
  resultPercent.cubic =
      (this->yCloud.framePercent.cubic / totalParticles) * 100;
  resultPercent.interfacial =
      (this->yCloud.framePercent.interfacial / totalParticles) * 100;
  resultPercent.clathrate =
      (this->yCloud.framePercent.clathrate / totalParticles) * 100;
  resultPercent.water =
      (this->yCloud.framePercent.water / totalParticles) * 100;
  resultPercent.interClathrate =
      (this->yCloud.framePercent.interClathrate / totalParticles) * 100;
  resultPercent.undef =
      (this->yCloud.framePercent.undef / totalParticles) * 100;
  return resultPercent;
}

void chill::bop::cleanUp() { delete snapshot; }

template class std::unique_ptr<chill::bop>;
