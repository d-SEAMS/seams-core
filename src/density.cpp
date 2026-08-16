//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <density.hpp>

#include <neighbours.hpp>

#include <cmath>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

double axisCoord(const molSys::Point<double> &p, int axis) {
  if (axis == 0) {
    return p.x;
  }
  if (axis == 1) {
    return p.y;
  }
  return p.z;
}

// Face area conjugate to axis: |det H| / L_axis, i.e. ly*lz for x.
double faceArea(const Cloud &yCloud, int axis) {
  double lengths[3] = {0.0, 0.0, 0.0};
  nneigh::dumpCellLengths(yCloud.box, yCloud.boxLow, lengths);
  if (axis < 0 || axis > 2 || lengths[axis] == 0.0) {
    return 0.0;
  }
  return nneigh::dumpVolume(yCloud) / lengths[axis];
}

site::DensityZ fillHistogram(const Cloud &yCloud,
                             const std::vector<double> &coords, int nbin,
                             int axis, int type) {
  site::DensityZ out;
  out.type = type;
  if (nbin < 1) {
    nbin = 1;
  }
  if (axis < 0 || axis > 2) {
    axis = 2;
  }

  double lo = 0.0;
  double span = 0.0;
  if (static_cast<int>(yCloud.boxLow.size()) > axis) {
    lo = yCloud.boxLow[static_cast<std::size_t>(axis)];
  }
  if (static_cast<int>(yCloud.box.size()) > axis) {
    span = yCloud.box[static_cast<std::size_t>(axis)];
  }
  if (span < 0.0) {
    span = 0.0;
  }
  const double dz = span / static_cast<double>(nbin);
  const double area = faceArea(yCloud, axis);
  const double slab = area * dz;

  out.z.assign(static_cast<std::size_t>(nbin), 0.0);
  out.rho.assign(static_cast<std::size_t>(nbin), 0.0);
  for (int ibin = 0; ibin < nbin; ++ibin) {
    out.z[static_cast<std::size_t>(ibin)] =
        lo + dz * (static_cast<double>(ibin) + 0.5);
  }

  if (dz <= 0.0) {
    return out;
  }

  std::vector<int> count(static_cast<std::size_t>(nbin), 0);
  for (double c : coords) {
    int ibin = static_cast<int>(std::floor((c - lo) / dz));
    if (ibin == nbin && c <= lo + span) {
      ibin = nbin - 1;
    }
    if (ibin < 0 || ibin >= nbin) {
      continue;
    }
    count[static_cast<std::size_t>(ibin)] += 1;
  }

  if (slab <= 0.0) {
    return out;
  }
  for (int ibin = 0; ibin < nbin; ++ibin) {
    out.rho[static_cast<std::size_t>(ibin)] =
        static_cast<double>(count[static_cast<std::size_t>(ibin)]) / slab;
  }
  return out;
}

} // namespace

site::DensityZ site::densityZ(const Cloud &yCloud, int typeI, int nbin,
                              int axis) {
  if (axis < 0 || axis > 2) {
    axis = 2;
  }
  std::vector<double> coords;
  coords.reserve(yCloud.pts.size());
  for (const auto &pt : yCloud.pts) {
    if (typeI != 0 && pt.type != typeI) {
      continue;
    }
    coords.push_back(axisCoord(pt, axis));
  }
  return fillHistogram(yCloud, coords, nbin, axis, typeI);
}

site::DensityZ site::densityZ(const Cloud &yCloud, const Table &table,
                              Kind kind, int nbin, int axis) {
  if (axis < 0 || axis > 2) {
    axis = 2;
  }
  const auto idx = indicesOf(yCloud, table, kind);
  std::vector<double> coords;
  coords.reserve(idx.size());
  for (int i : idx) {
    if (i < 0 || i >= static_cast<int>(yCloud.pts.size())) {
      continue;
    }
    coords.push_back(axisCoord(yCloud.pts[static_cast<std::size_t>(i)], axis));
  }
  return fillHistogram(yCloud, coords, nbin, axis, 0);
}
