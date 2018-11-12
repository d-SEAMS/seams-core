#ifndef _TRANSITION_H
#define _TRANSITION_H

#include <fstream>
#include <array>
#include <math.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "molecule.h"
#include "parameter.h"
#include "molecular_system.h"
#include "output.h"
#include "geometry.h"

/*! \brief Truncated minimal class for frame data.
 *          Blah
 *

 */

class TransitionSystem: public COutput, private CGeneric {
    double* frameAvg;
    double* currentDiff;
//  Super cool C++11
    std::array<double,3> coordHigh;
    std::array<double,3> coordLow;
    // Initializes Frames (CMolSys Obj)
    void prepFrame(int, std::string);
    // Calculates the relative change in position
    void frameDiff(int typeI, CMolecularSystem& frameOne, CMolecularSystem& frameTwo);
    // Somehow optimized (allegedly)
    static inline int isOdd(int);
    // Check if atom is within limits
    bool isThere(int iatom, CMolecularSystem& frame);

  public:
    //the main object where all properties of all particles are saved
    TransitionSystem();
    virtual ~TransitionSystem();
    // Frames we need
    CMolecularSystem* frameOne;
    CMolecularSystem* frameTwo;

    // Check for phase transitions
        bool hasTrans(int nop, int typeI, int frameNumOne, int frameNumTwo, std::array<double,3> coordHigh, std::array<double,3> coordLow, std::string fileName);
};

#endif
