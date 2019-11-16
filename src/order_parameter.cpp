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
std::vector<double> topoparam::calcCoverageArea(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> rings, double sheetArea) {
  //
  double areaXY, areaXZ, areaYZ;    // Total coverage area
  std::vector<double> singleAreas;  // Area of single rings

  // ---------------------------------------
  // Initialization
  areaXY = 0.0;
  areaXZ = 0.0;
  areaYZ = 0.0;
  // ---------------------------------------
  // Loop through all the rings
  for (int iring = 0; iring < rings.size(); iring++) {
    // Get the coverage area for the current ring
    singleAreas = topoparam::projAreaSingleRing(yCloud, rings[iring]);
    // Add these to the total coverage area
    areaXY += singleAreas[0];
    areaXZ += singleAreas[1];
    areaYZ += singleAreas[2];
  }  // end of loop through all the rings
  // ---------------------------------------
  // Normalize the coverage area by the sheet area
  areaXY = areaXY / sheetArea * 100.0;
  areaXZ = areaXZ / sheetArea * 100.0;
  areaYZ = areaYZ / sheetArea * 100.0;

  return {areaXY, areaXZ, areaYZ};
}

/********************************************/ /**
*  Calculates the coverage area/ projected area of a single ring
given the ring and the PointCloud.
***********************************************/
std::vector<double> topoparam::projAreaSingleRing(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> ring) {
  //
  int iatomIndex, jatomIndex;  // Atom indices of the i^th and j^th atoms
  int ringSize = ring.size();  // Number of nodes in the ring
  double areaXY, areaXZ, areaYZ;
  // ----------------------------------------
  // Calculate projected area onto the XY, YZ and XZ planes for basal1

  // Init the projected area
  areaXY = 0.0;
  areaXZ = 0.0;
  areaYZ = 0.0;

  jatomIndex = ring[0];

  // All points except the first pair
  for (int k = 1; k < ringSize; k++) {
    iatomIndex = ring[k];  // Current vertex

    // Add to the polygon area
    // ------
    // XY plane
    areaXY += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
              (yCloud->pts[jatomIndex].y - yCloud->pts[iatomIndex].y);
    // ------
    // XZ plane
    areaXZ += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
              (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
    // ------
    // YZ plane
    areaYZ += (yCloud->pts[jatomIndex].y + yCloud->pts[iatomIndex].y) *
              (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
    // ------
    jatomIndex = iatomIndex;
  }

  // Closure point
  iatomIndex = ring[0];
  // ------
  // XY plane
  areaXY += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
            (yCloud->pts[jatomIndex].y - yCloud->pts[iatomIndex].y);
  // ------
  // XZ plane
  areaXZ += (yCloud->pts[jatomIndex].x + yCloud->pts[iatomIndex].x) *
            (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
  // ------
  // YZ plane
  areaYZ += (yCloud->pts[jatomIndex].y + yCloud->pts[iatomIndex].y) *
            (yCloud->pts[jatomIndex].z - yCloud->pts[iatomIndex].z);
  // ------
  // The actual projected area is half of this
  areaXY *= 0.5;
  areaXZ *= 0.5;
  areaYZ *= 0.5;

  // Absolute area
  areaXY = fabs(areaXY);
  areaXZ = fabs(areaXZ);
  areaYZ = fabs(areaYZ);

  return {areaXY, areaYZ, areaXZ};
}