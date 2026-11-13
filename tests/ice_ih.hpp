#pragma once

// Ideal ice Ih, eight oxygens in the orthorhombic cell of
// scripts/gen_ideal_ice.py, tiled by (nx, ny, nz). The oxygen bond is
// 2.75 A. Coordinates stay in memory.

#include <cmath>
#include <vector>

struct IceIhCell {
  std::vector<double> xyz;
  double box9[9]{};
  int n = 0;
};

inline IceIhCell iceIhOrtho(int nx, int ny, int nz) {
  const double bond = 2.75;
  const double z = 1.0 / 16.0;
  const double frac[4][3] = {{1.0 / 3.0, 2.0 / 3.0, z},
                             {2.0 / 3.0, 1.0 / 3.0, 0.5 + z},
                             {2.0 / 3.0, 1.0 / 3.0, 1.0 - z},
                             {1.0 / 3.0, 2.0 / 3.0, 0.5 - z}};
  const double c = bond / 0.375;
  const double a = c / std::sqrt(8.0 / 3.0);
  const double hv[3][3] = {{a, 0.0, 0.0},
                           {-0.5 * a, a * std::sqrt(3.0) / 2.0, 0.0},
                           {0.0, 0.0, c}};
  const double lengths[3] = {a, a * std::sqrt(3.0), c};
  double basis[8][3];
  for (int b = 0; b < 4; ++b) {
    for (int d = 0; d < 3; ++d) {
      basis[b][d] =
          frac[b][0] * hv[0][d] + frac[b][1] * hv[1][d] + frac[b][2] * hv[2][d];
      basis[b + 4][d] = basis[b][d] + hv[1][d];
    }
  }
  for (int b = 0; b < 8; ++b) {
    for (int d = 0; d < 3; ++d) {
      const double q = basis[b][d] / lengths[d] + 1e-8;
      basis[b][d] -= std::floor(q) * lengths[d];
    }
  }
  IceIhCell cell;
  cell.n = nx * ny * nz * 8;
  cell.xyz.resize(static_cast<std::size_t>(cell.n) * 3);
  cell.box9[3] = nx * lengths[0];
  cell.box9[4] = ny * lengths[1];
  cell.box9[5] = nz * lengths[2];
  int p = 0;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        for (int b = 0; b < 8; ++b, ++p) {
          cell.xyz[static_cast<std::size_t>(3 * p)] = i * lengths[0] + basis[b][0];
          cell.xyz[static_cast<std::size_t>(3 * p + 1)] = j * lengths[1] + basis[b][1];
          cell.xyz[static_cast<std::size_t>(3 * p + 2)] = k * lengths[2] + basis[b][2];
        }
      }
    }
  }
  return cell;
}
