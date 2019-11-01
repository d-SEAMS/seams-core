#ifndef __ORDER_PARAMETER_H_
#define __ORDER_PARAMETER_H_

#include <math.h>
#include <sys/stat.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <cage.hpp>
#include <mol_sys.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

namespace topoparam {

// Calculates the height%, an average measure of filled volume. The average
// height of a prism can be taken to be 2.75-2.85 Angstrom. (Koga et. al., 2001)
double normHeightPercent(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int nPrisms,
    double avgPrismHeight);

}  // namespace topoparam

#endif  // __ORDER_PARAMETER_H_
