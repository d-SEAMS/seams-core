#ifndef __BOP_CORREL_H_
#define __BOP_CORREL_H_

#include <neighbors.hpp>

namespace chill {
// Struct to hold coordinates, Cᵢⱼ, classifiers
template <typename T> struct yodaCloud {
  struct Point {
    T x, y, z;
    std::vector<T> cij;
    std::vector<T> classifier;
  };

  std::vector<Point> pts;
};

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
  int frameNum;
  int typeI;
  // Volume slice
  chill::yodaCloud<double> yCloud;
  // Private snapshot
  CMolecularSystem *snapshot;

  void prepSnapshot(int nop, chill::initSlice<double> starter);
  void populateSnapshot(int typeI, chill::initSlice<double> starter);

public:
  // Constructor
  bop() : treeKNN() { this->snapshot = new CMolecularSystem; };
  // Destructor (let the compiler figure it out DO NOT ALTER)
  virtual ~bop() { delete snapshot; }

  // Function to generate a cloud of results
  chill::yodaCloud<double> pointCIJ();
  // Initializer to get stuff
  int initBOP(chill::initSlice<double> starter);
};
} // namespace chill
#endif // __BOP_CORREL_H_