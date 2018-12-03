#ifndef __BOP_CORREL_H_
#define __BOP_CORREL_H_

#include <complex>
#include <neighbors.hpp>

namespace chill {

// Struct for one point
template <typename T> struct yodaPoint {
  T x, y, z;
  std::vector<T> cij;
  std::vector<T> classifier;
  bool inSlice = false;
  std::vector<std::complex<double>> Q;
  std::vector<int> nearestID;
};

// Struct to hold coordinates, Cᵢⱼ, classifiers
template <typename T> struct yodaCloud { std::vector<T> pts; };

// Input object
template <typename T> struct initSlice {

  std::array<T, 3> coordHigh = {0, 0, 0};
  std::array<T, 3> coordLow = {0, 0, 0};
  std::string filename;
  std::array<int, 2> frameRange = {0, 0};
};

class bop : public neigh::treeKNN {
private:
  int nop;
  std::string filename;
  int typeI;
  // Volume slice
  chill::yodaCloud<yodaPoint<double>> yCloud;
  // Private snapshot
  CMolecularSystem *snapshot;

  void prepSnapshot(chill::initSlice<double> starter);
  void populateSnapshot(chill::initSlice<double> starter);

public:
  // Constructor
  bop() : treeKNN() { this->snapshot = new CMolecularSystem; };
  // Destructor (let the compiler figure it out DO NOT ALTER)
  virtual ~bop() { delete snapshot; }

  // Function to generate Q per point
  chill::yodaPoint<double> pointQ(int queryIndex);
  // Function to generate Cij per pair
  chill::yodaPoint<double> pointCij(int queryIndex);
  // Initializer to get stuff
  int initBOP(int nop, int typeI, chill::initSlice<double> starter);
  chill::yodaPoint<double> atomVerdict(int queryIndex);
};
} // namespace chill
#endif // __BOP_CORREL_H_
