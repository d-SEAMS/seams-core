// Tiled interface frames shared by the rank benchmark and the MPI tests.
// Ice is eight molecules per cubic cell. The sI cell is the 46-water
// genice frame. The upper half of z is shifted by a fixed generator so
// every rank builds the same cloud.

#ifndef SEAMS_TESTS_SLAB_FRAMES_HPP_
#define SEAMS_TESTS_SLAB_FRAMES_HPP_

#include <mol_sys.hpp>

#include <cmath>
#include <cstdint>

namespace seams::frames {

inline molSys::PointCloud<molSys::Point<double>, double> iceLattice(int mx, int my,
                                                                    int mz) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double bond = 2.75;
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  const double basis[8][3] = {{0, 0, 0},          {0.5, 0.5, 0},      {0.5, 0, 0.5},
                              {0, 0.5, 0.5},      {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  int id = 1;
  for (int i = 0; i < mx; i++) {
    for (int j = 0; j < my; j++) {
      for (int k = 0; k < mz; k++) {
        for (const auto &b : basis) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          p.x = (i + b[0]) * lattice;
          p.y = (j + b[1]) * lattice;
          p.z = (k + b[2]) * lattice;
          cloud.pts.push_back(p);
          cloud.idIndexMap[id] = id - 1;
          id++;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.currentFrame = 1;
  cloud.box = {mx * lattice, my * lattice, mz * lattice};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

// sI unit cell from input/traj/genice_sI.lammpstrj. 46 water molecules.
inline constexpr double kSI = 12.0034886;
inline constexpr int kSICount = 46;
inline constexpr double kSIWaters[kSICount][3] = {
    {10.64000000, 3.78000000, 0.00000000}, {3.75000000, 0.04000000, 1.43000000},
    {0.05000000, 3.00000000, 6.05000000},  {8.18000000, 3.79000000, 8.14000000},
    {2.18000000, 2.21000000, 2.26000000},  {0.05000000, 9.00000000, 5.96000000},
    {2.21000000, 9.82000000, 2.13000000},  {6.04000000, 9.70000000, 4.65000000},
    {7.43000000, 6.00000000, 9.66000000},  {9.79000000, 2.26000000, 9.82000000},
    {3.87000000, 3.82000000, 3.79000000},  {2.34000000, 4.58000000, 6.00000000},
    {3.82000000, 8.26000000, 3.79000000},  {3.79000000, 3.82000000, 8.26000000},
    {8.30000000, 11.96348860, 1.35000000}, {6.04000000, 2.25000000, 4.57000000},
    {9.82000000, 9.79000000, 2.26000000},  {5.96000000, 9.75000000, 7.44000000},
    {9.78000000, 4.64000000, 6.00000000},  {3.70000000, 11.96348860, 10.65000000},
    {5.96000000, 2.30000000, 7.35000000},  {3.79000000, 8.14000000, 8.18000000},
    {10.57000000, 8.26000000, 0.04000000}, {9.80000000, 9.82000000, 9.87000000},
    {11.96348860, 1.43000000, 3.75000000}, {9.80000000, 2.18000000, 2.13000000},
    {9.75000000, 7.44000000, 6.04000000},  {8.94000000, 6.00000000, 0.00000000},
    {11.96348860, 10.57000000, 8.26000000},{8.18000000, 8.26000000, 8.21000000},
    {6.05000000, 11.95348860, 3.00000000}, {11.96348860, 1.35000000, 8.30000000},
    {2.30000000, 7.35000000, 5.96000000},  {4.57000000, 6.04000000, 2.25000000},
    {4.64000000, 6.00000000, 9.78000000},  {2.18000000, 2.13000000, 9.80000000},
    {0.04000000, 10.65000000, 3.70000000}, {1.35000000, 8.30000000, 11.96348860},
    {1.43000000, 3.66000000, 0.00000000},  {8.26000000, 3.79000000, 3.82000000},
    {8.21000000, 8.14000000, 3.82000000},  {7.35000000, 5.96000000, 2.30000000},
    {2.18000000, 9.79000000, 9.75000000},  {5.96000000, 0.05000000, 9.00000000},
    {2.94000000, 6.00000000, 0.00000000},  {8.26000000, 0.04000000, 10.57000000},
};

inline molSys::PointCloud<molSys::Point<double>, double> tileSI(int nx, int ny, int nz) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {kSI * nx, kSI * ny, kSI * nz};
  cloud.boxLow = {0.0, 0.0, 0.0};
  int id = 1;
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        for (int w = 0; w < kSICount; w++) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          p.x = kSIWaters[w][0] + ix * kSI;
          p.y = kSIWaters[w][1] + iy * kSI;
          p.z = kSIWaters[w][2] + iz * kSI;
          cloud.pts.push_back(p);
          cloud.idIndexMap[id] = id - 1;
          id++;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.currentFrame = 1;
  return cloud;
}

// Deterministic shift of every molecule at or above the mid-plane.
// A coordinate that wraps into the lower half is folded back up.
inline void meltUpperHalf(molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  const double lx = cloud.box[0];
  const double ly = cloud.box[1];
  const double lz = cloud.box[2];
  const double mid = 0.5 * lz;
  auto wrap = [](double v, double length) {
    v -= length * std::floor(v / length);
    return v < 0.0 ? v + length : v;
  };
  for (auto &p : cloud.pts) {
    if (p.z < mid) {
      continue;
    }
    std::uint32_t state = static_cast<std::uint32_t>(p.atomID) * 1664525u + 1013904223u;
    p.x += (static_cast<int>(state % 1000) / 1000.0 - 0.5) * 3.5;
    state = state * 1664525u + 1013904223u;
    p.y += (static_cast<int>(state % 1000) / 1000.0 - 0.5) * 3.5;
    state = state * 1664525u + 1013904223u;
    p.z += (static_cast<int>(state % 1000) / 1000.0 - 0.5) * 2.0;
    p.x = wrap(p.x, lx);
    p.y = wrap(p.y, ly);
    p.z = wrap(p.z, lz);
    if (p.z < mid) {
      p.z += mid;
      if (p.z >= lz) {
        p.z = mid + 0.1;
      }
    }
  }
}

inline molSys::PointCloud<molSys::Point<double>, double> iceWaterSlab(int cells) {
  auto cloud = iceLattice(cells, cells, cells);
  meltUpperHalf(cloud);
  return cloud;
}

inline molSys::PointCloud<molSys::Point<double>, double> hydrateLiquidSlab(int nx, int ny,
                                                                           int nz) {
  auto cloud = tileSI(nx, ny, nz);
  meltUpperHalf(cloud);
  return cloud;
}

} // namespace seams::frames

#endif // SEAMS_TESTS_SLAB_FRAMES_HPP_
