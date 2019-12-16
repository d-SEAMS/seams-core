#ifndef __SHAPEMATCH_H_
#define __SHAPEMATCH_H_

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <math.h>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include <absOrientation.hpp>
#include <mol_sys.hpp>
#include <pntCorrespondence.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

namespace match {

// Shape-matching for a pair of polygon basal rings. Returns true if the pair of
// basal rings form a prism block.
bool matchPrismBlock(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                     std::vector<std::vector<int>> nList,
                     const Eigen::MatrixXd &refPoints, std::vector<int> *basal1,
                     std::vector<int> *basal2, std::vector<double> *rmsdPerAtom,
                     bool isPerfect = true);

} // namespace match

#endif // __SHAPEMATCH_H_
