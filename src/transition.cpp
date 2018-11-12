#include "transition.h"
#include "geometry.h"

/********************************************//**
 *  Constructor
 ***********************************************/
TransitionSystem::TransitionSystem()
{
    // Initialize things
  this->frameOne = new CMolecularSystem;
  this->frameTwo = new CMolecularSystem;
  // Volume Limits
  this->coordLow = this->coordHigh = {0,0,0};
  // Allocate size for the frame diffs
  this->frameAvg=NULL;
  this->currentDiff=NULL;
}

/********************************************//**
 *  Destructor
 ***********************************************/
TransitionSystem::~TransitionSystem()
{
  delete frameTwo;
  delete frameOne;
  delete [] frameAvg;
  delete [] currentDiff;
}

// Handles things pertaining to the phase transition
bool TransitionSystem::hasTrans(int nop, int typeI, int frameNumOne, int frameNumTwo, std::array<double,3> coordHigh, std::array<double,3> coordLow, std::string fileName){
    int loopIndex = 1, frameNumber=frameNumOne+1;
    this->coordLow = coordLow;
    this->coordHigh = coordHigh;
    // Initialize the objects
    this->prepFrame(nop, fileName);
    this->frameOne->readParticleFile(frameNumOne);
    while (1) {
        if (frameNumber > frameNumTwo) {
            break;
        }
// Efficiently skip reading objects and work with indices
        if (this->isOdd(loopIndex)) {
            this->frameTwo->readParticleFile(frameNumber);
            // TODO: Sort wrt atomID using sort() [NOT qsort]
        }
        else {
            this->frameOne->readParticleFile(frameNumber);
        }
        // Diff with standard functions
        frameNumber++;
        loopIndex++;
}
    return 1;
}

void TransitionSystem::prepFrame (int nop, std::string fileName) {
    this->frameOne->initializeFrames(nop, fileName);
    this->frameTwo->initializeFrames(nop, fileName);
    this->frameAvg = new double [nop];
    this->currentDiff= new double [nop];
    for (int i=0; i < nop; i++) {
        this->frameAvg[i]=-1;
        this->currentDiff[i]=-1;
    }
}

inline int TransitionSystem::isOdd(int x) { return x & 1; }

void TransitionSystem::frameDiff(int typeI, CMolecularSystem& frameOne, CMolecularSystem& frameTwo) {

    int nop = frameOne.parameter->nop;

    for (int iatom=0; iatom < nop; iatom++)
    {
        if (frameOne.molecules[iatom].type==typeI) {
            if (this->isThere(iatom, frameOne) && this->isThere(iatom, frameTwo)) {
                CGeneric::getAbsDistance(iatom, frameOne, frameTwo);
            }
            // Do nothing
        }
    }
}

bool TransitionSystem::isThere(int iatom, CMolecularSystem& frame) {
    // TODO: Migrate to CGeneric
                int iter=0;
                double coord = frame.molecules[iatom].get_posx();
                // TODO: Handle non x-dimension things
            for (auto const& value: this->coordHigh) {
                if (value==coordLow[iter]) {
                    return true;
                }
                else if (coord >= coordLow[iter] && coord <= value) {
                    return true;
                }
                else {
                    return false;
                }
            iter++;
            }
}
