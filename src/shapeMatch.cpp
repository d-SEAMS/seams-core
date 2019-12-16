#include <shapeMatch.hpp>

// Shape-matching for a pair of polygon basal rings. Returns true if the pair of
// basal rings form a prism block.
bool match::matchPrismBlock(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, const Eigen::MatrixXd &refPoints,
    std::vector<int> *basal1, std::vector<int> *basal2,
    std::vector<double> *rmsdPerAtom, bool isPerfect) {
  //
  int ringSize = (*basal1).size(); // Number of nodes in each basal ring
  std::vector<std::vector<int>> refPntToPnt; // Vector of vector of ints with
  // the connectivity information for the two basal rings
  std::vector<int> matchedBasal1,
      matchedBasal2; // Re-ordered basal rings 1 and 2

  // -----------------------
  // Getting the target Eigen vectors
  // Get the re-ordered matched basal rings, ordered with respect to each other
  pntToPnt::relOrderPrismBlock(yCloud, *basal1, *basal2, nList, &matchedBasal1,
                               &matchedBasal2);
  // -----------------------

  // Change this later
  return false;
} // end of function