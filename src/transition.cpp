#include "geometry.h"
#include "transition.h"

using namespace trans;

/********************************************/ /**
 *  Constructor
 ***********************************************/
TransitionSystem::TransitionSystem() {
  // Initialize things
  this->frameOne = new CMolecularSystem;
  this->frameTwo = new CMolecularSystem;
  // Volume Limits
  this->coordLow = this->coordHigh = {0, 0, 0};
  // Allocate size for the frame diffs
  this->currentDiff = nullptr;
}

/********************************************/ /**
 *  Destructor
 ***********************************************/
TransitionSystem::~TransitionSystem() {
  delete frameTwo;
  delete frameOne;
  delete[] currentDiff;
}

// Handles things pertaining to the phase transition
void TransitionSystem::mightTrans(int nop, int typeI, int frameNumOne,
                                  int frameNumTwo,
                                  std::array<double, 3> coordHigh,
                                  std::array<double, 3> coordLow,
                                  std::string fileName, int nstep) {
  int loopIndex = 1, frameNumber = frameNumOne + 1;
  double average = 0, currentAvg = 0;
  this->coordLow = coordLow;
  this->coordHigh = coordHigh;
  // Initialize the objects
  this->prepFrame(nop, fileName);
  this->frameOne->parameter->nsteps = nstep;
  this->frameTwo->parameter->nsteps = nstep;
  this->frameOne->readParticleFile(frameNumOne);
  while (1) {
    if (frameNumber > frameNumTwo) {
      break;
    }
    // Efficiently skip reading objects and work with indices
    if (this->isOdd(loopIndex)) {
      this->frameTwo->readParticleFile(frameNumber);
      // TODO: Sort wrt atomID using sort() [NOT qsort]
    } else {
      this->frameOne->readParticleFile(frameNumber);
    }
    // Diff with standard functions
    this->frameDiff(typeI, frameOne, frameTwo);
    currentAvg = this->timeAtomAvg(nop);
    average = (average + currentAvg) * 0.5;
    // Output stuff here
    std::cout << frameNumber << "  from mightTrans " << currentAvg << std::endl;
    frameNumber++;
    loopIndex++;
  }
  std::cout << "\nWe are done now\n"
            << "The final average is " << average << std::endl;
}

void TransitionSystem::prepFrame(int nop, std::string fileName) {
  this->frameOne->initializeFrames(nop, fileName);
  this->frameTwo->initializeFrames(nop, fileName);
  this->currentDiff = new double[nop];
  for (int i = 0; i < nop; i++) {
    this->currentDiff[i] = -1;
  }
}

inline int TransitionSystem::isOdd(int x) { return x & 1; }

void TransitionSystem::frameDiff(int typeI, CMolecularSystem *frameOne,
                                 CMolecularSystem *frameTwo) {

  int nop = frameOne->parameter->nop;

  for (int iatom = 0; iatom < nop; iatom++) {
    // By default we do not consider the atom
    this->currentDiff[iatom] = -1;
    // If the type matches
    if (frameOne->molecules[iatom].type == typeI) {
      // If the atom is within limits
      if (this->isThere(iatom, frameOne) && this->isThere(iatom, frameTwo)) {
        // This is in time. i.e this is the absolute difference in
        // distance for each valid atom over two frames
        this->currentDiff[iatom] =
            this->CGeneric::getAbsDistance(iatom, frameOne, frameTwo);
        // std::cout<<"FrameDiff "<<iatom<<" has absDist "<<this->currentDiff[iatom]<<std::endl;
      }
      // Do nothing
    }
  }
}

bool TransitionSystem::isThere(int iatom, CMolecularSystem *frame) {
  // TODO: Migrate to CGeneric
  double coordX = frame->molecules[iatom].get_posx();
  double coordY = frame->molecules[iatom].get_posy();
  double coordZ = frame->molecules[iatom].get_posz();
  std::array<double, 3> coord = {coordX, coordY, coordZ};
  // TODO: Handle non x-dimension things
  for (int i = 0; i < 3; i++) {
    if (coordHigh[i] == coordLow[i]) {
      return true;
    } else if (coord[i] >= coordLow[i] && coord[i] <= coordHigh[i]) {
      return true;
    } else {
      return false;
    }
  }
  return false;
}

double TransitionSystem::timeAtomAvg(int nop) {
  int norm = 0;
  double sum = 0;
  for (int i = 0; i <= nop; i++) {
    if (this->currentDiff[i] == -1) {
      continue;
    } else {
      norm++;
      sum += this->currentDiff[i];
      // for (auto l=0; l < nop; l++) {

      // std::cout<<"Current Diff for "<<l<<" "<<this->currentDiff[l]<<std::endl;
      std::cout << "Current Diff for " << i << " " << this->currentDiff[i]
                << std::endl;
      // }
    }
  }
  if (norm != 0) {
    return sum / norm;
  }
  // TODO: Error handling
  return -1000;
}
