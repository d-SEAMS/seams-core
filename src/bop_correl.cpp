// For calculating the "correlation of orientational bond order parameter
// developed by ten Wolde et al. based on the Steinhardt parameters"
// Essentially the CHILL+ params by Molinero et al.
// This is a library.
#include <bop_correl.hpp>

/********************************************/ /**
 *  Public Functions
 ***********************************************/

// void chill::test::sim() { std::cout << "sim"; };
int chill::bop::initBOP(chill::initSlice<double> starter) {
  std::cout << "I am bop\n";
  return 1;
}

void chill::bop::prepSnapshot(int nop, chill::initSlice<double> starter) {
  this->snapshot->initializeFrames(nop, starter.filename);
}
void chill::bop::populateSnapshot(int typeI, chill::initSlice<double> starter) {
  int filteredParticles(0);
  int nop(snapshot->parameter->nop);
  int dummy(0);
  // Do filtering
  for (int t; t < nop; t++) {
    // Match type
    if (snapshot->molecules[t].type == typeI) {
      // Check limits
      if (treeKNN::isThere(t, snapshot, starter.coordHigh, starter.coordLow)) {
        // Accomodate one more point
        yCloud.pts.resize(yCloud.pts.size() + 1);
        // Dump point
        double coordX = snapshot->molecules[t].get_posx();
        double coordY = snapshot->molecules[t].get_posy();
        double coordZ = snapshot->molecules[t].get_posz();
        yCloud.pts[dummy].x = coordX;
        yCloud.pts[dummy].y = coordY;
        yCloud.pts[dummy].z = coordZ;
        // Update dummy
        dummy++;
      }
    }
  }
}
