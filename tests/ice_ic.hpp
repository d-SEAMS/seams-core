#pragma once

// Ideal ice Ic, the 8-atom cubic cell tiled on a cube. The oxygen bond
// is 2.75 A, the same length as the shell classifier. Coordinates stay
// in memory.

#include <mol_sys.hpp>

#include <cmath>
#include <vector>

inline molSys::PointCloud<molSys::Point<double>, double> iceIc(int molecules) {
  using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
  int side = 0;
  for (int s = 1; s <= 64; ++s) {
    if (8 * s * s * s == molecules) {
      side = s;
      break;
    }
  }
  Cloud cloud;
  if (side == 0) {
    return cloud;
  }
  const double bond = 2.75;
  const double a = 4.0 * bond / std::sqrt(3.0);
  const double length = a * static_cast<double>(side);
  const double basis[8][3] = {
      {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5},
      {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25}, {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75},
  };
  cloud.box = {length, length, length};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = molecules;
  cloud.pts.resize(static_cast<std::size_t>(molecules));
  int at = 0;
  for (int ix = 0; ix < side; ++ix) {
    for (int iy = 0; iy < side; ++iy) {
      for (int iz = 0; iz < side; ++iz) {
        for (int b = 0; b < 8; ++b) {
          auto &point = cloud.pts[static_cast<std::size_t>(at)];
          point.x = (static_cast<double>(ix) + basis[b][0]) / side * length;
          point.y = (static_cast<double>(iy) + basis[b][1]) / side * length;
          point.z = (static_cast<double>(iz) + basis[b][2]) / side * length;
          point.type = 1;
          point.atomID = at + 1;
          cloud.idIndexMap[point.atomID] = at;
          ++at;
        }
      }
    }
  }
  return cloud;
}

// Move every sixteenth oxygen by 1.2 A so its bond row changes, and
// wrap the coordinate back into the box. The stride leaves enough
// changed vertices that the carried update has a team of work.
inline molSys::PointCloud<molSys::Point<double>, double>
perturbIce(molSys::PointCloud<molSys::Point<double>, double> cloud) {
  const double length = cloud.box.empty() ? 1.0 : cloud.box[0];
  for (int i = 0; i < cloud.nop; i += 16) {
    auto &point = cloud.pts[static_cast<std::size_t>(i)];
    point.x += 1.2;
    if (point.x >= length) {
      point.x -= length;
    }
  }
  return cloud;
}
