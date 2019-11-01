#include <order_parameter.hpp>

/********************************************/ /**
*  The average height of prism blocks remains relatively constant. We have
observed a average prism heights of 2.7-2.85 Angstrom for prisms irrespective of
the number of nodes.
The equation is given by:

  \f[
  Height_{n}% = \frac{N_n}{N_{max}} \times 100
  \f]

Here, \f$N_{max} = H_{SWCT}/h_{avg}f$ and \f$N_{n}$ is the
number of prism blocks for n-sided prismatic phase.

This means that the normalization factor, is the same for
every node number \f$n\f$.
***********************************************/
double topoparam::normHeightPercent(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int nPrisms,
    double avgPrismHeight) {
  //
  double hPercent;        // Normalized height percent
  double nanoTubeHeight;  // Height of the SWCT
  double numberMax;  // Maximum number possible, given the average prism height

  // ---------------------------------------
  // Calculate the height of the SWCT
  // This is the longest dimension of the simulation box
  nanoTubeHeight = *max_element(yCloud->box.begin(), yCloud->box.end());
  // ---------------------------------------
  // Calculate the maximum possible height, given the average prism height
  // and the height of the nanotube
  numberMax = nanoTubeHeight / avgPrismHeight;
  // ---------------------------------------
  // Calculate the normalized height percentage
  hPercent = nPrisms / numberMax * 100.0;

  return hPercent;
}
