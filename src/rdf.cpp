//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <rdf.hpp>

#include <generic.hpp>
#include <neighbours.hpp>

/**
 * @details Partial g_IJ(r) for one frame. Distances are the dump
 *  minimum image. Each unordered I-J pair is counted once. The
 *  spherical-shell ideal-gas count uses nneigh::dumpVolume.
 */
rdf::PartialRdf rdf::partialRdf(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeI,
    int typeJ, double rmax, int nbins) {
  PartialRdf out;
  out.rmax = rmax;
  out.typeI = typeI;
  out.typeJ = typeJ;
  if (nbins < 1) {
    nbins = 1;
  }
  out.binwidth = (rmax > 0.0) ? rmax / static_cast<double>(nbins) : 0.0;
  out.volume = nneigh::dumpVolume(yCloud);
  out.r.assign(static_cast<std::size_t>(nbins), 0.0);
  out.g.assign(static_cast<std::size_t>(nbins), 0.0);
  out.count.assign(static_cast<std::size_t>(nbins), 0);
  for (int ibin = 0; ibin < nbins; ++ibin) {
    out.r[static_cast<std::size_t>(ibin)] =
        out.binwidth * (static_cast<double>(ibin) + 0.5);
  }

  for (const auto &pt : yCloud.pts) {
    if (pt.type == typeI) {
      ++out.nI;
    }
    if (pt.type == typeJ) {
      ++out.nJ;
    }
  }

  if (rmax <= 0.0 || yCloud.box.size() < 3 || yCloud.nop == 0) {
    return out;
  }

  const auto nList = nneigh::neighListPair(rmax, yCloud, typeI, typeJ);
  const bool like = typeI == typeJ;
  for (int iatom = 0; iatom < yCloud.nop; ++iatom) {
    if (yCloud.pts[static_cast<std::size_t>(iatom)].type != typeI) {
      continue;
    }
    if (static_cast<int>(nList.size()) <= iatom || nList[iatom].size() < 2) {
      continue;
    }
    for (std::size_t k = 1; k < nList[iatom].size(); ++k) {
      const int jid = nList[iatom][k];
      const auto found = yCloud.idIndexMap.find(jid);
      if (found == yCloud.idIndexMap.end()) {
        continue;
      }
      const int jatom = found->second;
      if (jatom < 0 || jatom >= yCloud.nop) {
        continue;
      }
      if (yCloud.pts[static_cast<std::size_t>(jatom)].type != typeJ) {
        continue;
      }
      if (like && jatom <= iatom) {
        continue;
      }
      const double rij = gen::periodicDist(yCloud, iatom, jatom);
      if (rij > rmax) {
        continue;
      }
      int ibin = (out.binwidth > 0.0)
                     ? static_cast<int>(rij / out.binwidth)
                     : 0;
      if (ibin < 0) {
        continue;
      }
      if (ibin >= nbins) {
        ibin = nbins - 1;
      }
      out.count[static_cast<std::size_t>(ibin)] += 1;
    }
  }

  for (int ibin = 0; ibin < nbins; ++ibin) {
    const double r0 = out.binwidth * static_cast<double>(ibin);
    const double r1 = out.binwidth * static_cast<double>(ibin + 1);
    const double shell =
        4.0 * gen::pi * (r1 * r1 * r1 - r0 * r0 * r0) / 3.0;
    if (shell <= 0.0 || out.volume <= 0.0) {
      continue;
    }
    double denom = 0.0;
    if (like) {
      if (out.nI < 2) {
        continue;
      }
      denom = 0.5 * static_cast<double>(out.nI) *
              static_cast<double>(out.nI - 1);
    } else {
      if (out.nI < 1 || out.nJ < 1) {
        continue;
      }
      denom = static_cast<double>(out.nI) * static_cast<double>(out.nJ);
    }
    const double ideal = denom * shell / out.volume;
    if (ideal <= 0.0) {
      continue;
    }
    out.g[static_cast<std::size_t>(ibin)] =
        static_cast<double>(out.count[static_cast<std::size_t>(ibin)]) / ideal;
  }
  return out;
}

/**
 * @details Discrete running integral 4 pi rho_J int_0^r s^2 g(s) ds.
 *  Each bin contributes a spherical shell of width binwidth. The
 *  value at bin i is the CN at the outer edge of that bin.
 */
std::vector<double> rdf::runningCN(const PartialRdf &h, double rhoJ) {
  std::vector<double> cn(h.g.size(), 0.0);
  if (h.binwidth <= 0.0 || rhoJ == 0.0) {
    return cn;
  }
  double acc = 0.0;
  for (std::size_t i = 0; i < h.g.size(); ++i) {
    const double r0 = h.binwidth * static_cast<double>(i);
    const double r1 = h.binwidth * static_cast<double>(i + 1);
    acc += 4.0 * gen::pi * rhoJ * h.g[i] * (r1 * r1 * r1 - r0 * r0 * r0) /
           3.0;
    cn[i] = acc;
  }
  return cn;
}

/**
 * @details First local maximum of g, then the first local minimum after
 *  that peak. A confirmed minimum needs a descent and a later rise.
 *  Monotonic or empty histograms return -1.
 */
int rdf::firstMinimumBin(const PartialRdf &h) {
  const int n = static_cast<int>(h.g.size());
  if (n < 3) {
    return -1;
  }
  int i = 0;
  while (i + 1 < n && h.g[static_cast<std::size_t>(i + 1)] >=
                           h.g[static_cast<std::size_t>(i)]) {
    ++i;
  }
  const int imax = i;
  if (imax >= n - 1) {
    return -1;
  }
  i = imax;
  while (i + 1 < n && h.g[static_cast<std::size_t>(i + 1)] <=
                           h.g[static_cast<std::size_t>(i)]) {
    ++i;
  }
  if (i == imax || i + 1 >= n) {
    return -1;
  }
  return i;
}

/**
 * @details Site-site CN up to rMax. Partial last bins are cut at rMax.
 */
double rdf::coordinationNumber(const PartialRdf &h, double rMax,
                               double rhoJ) {
  if (rMax <= 0.0 || h.binwidth <= 0.0 || rhoJ == 0.0 || h.g.empty()) {
    return 0.0;
  }
  double cn = 0.0;
  for (std::size_t i = 0; i < h.g.size(); ++i) {
    const double r0 = h.binwidth * static_cast<double>(i);
    if (r0 >= rMax) {
      break;
    }
    double r1 = h.binwidth * static_cast<double>(i + 1);
    if (r1 > rMax) {
      r1 = rMax;
    }
    cn += 4.0 * gen::pi * rhoJ * h.g[i] * (r1 * r1 * r1 - r0 * r0 * r0) / 3.0;
  }
  return cn;
}
