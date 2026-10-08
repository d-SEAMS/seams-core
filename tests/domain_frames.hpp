#ifndef SEAMS_TESTS_DOMAIN_FRAMES_HPP_
#define SEAMS_TESTS_DOMAIN_FRAMES_HPP_

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <random>
#include <vector>

// Uniform atoms in an L-edge box tilted by xy, xz, yz, stored as a LAMMPS
// dump stores it: bound spans in box[0..2] and the tilts in box[3..5]
inline molSys::PointCloud<molSys::Point<double>, double>
tiltedFrame(int n, double L, double xy, double xz, double yz, unsigned seed) {
  molSys::PointCloud<molSys::Point<double>, double> c;
  if (xy == 0.0 && xz == 0.0 && yz == 0.0) {
    c.box = {L, L, L};
    c.boxLow = {0.0, 0.0, 0.0};
  } else {
    const double xmin = std::min({0.0, xy, xz, xy + xz});
    const double xmax = std::max({0.0, xy, xz, xy + xz});
    c.box = {L + xmax - xmin, L + std::max(0.0, yz) - std::min(0.0, yz), L,
             xy, xz, yz};
    c.boxLow = {xmin, std::min(0.0, yz), 0.0};
  }
  c.nop = n;
  c.currentFrame = 1;
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int i = 0; i < n; i++) {
    const double sx = u(rng);
    const double sy = u(rng);
    const double sz = u(rng);
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    p.x = L * sx + xy * sy + xz * sz;
    p.y = L * sy + yz * sz;
    p.z = L * sz;
    c.pts.push_back(p);
    c.idIndexMap[i + 1] = i;
  }
  return c;
}

// The whole frame's rings on its cutoff graph with ascending rows
inline std::vector<std::vector<int>>
sortedRowRings(const molSys::PointCloud<molSys::Point<double>, double> &c,
               double cutoff, int depth) {
  auto nList = nneigh::getNewNeighbourListByIndex(c, cutoff);
  for (auto &row : nList) {
    if (row.size() > 2) {
      std::sort(row.begin() + 1, row.end());
    }
  }
  return primitive::ringNetwork(nList, depth);
}

#endif // SEAMS_TESTS_DOMAIN_FRAMES_HPP_
