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

// Internal Libraries
#include <bond.hpp>
#include <generic.hpp>
#include <index_cast.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

namespace {

// Water populateHbonds* keep the two-H-per-molID set from hAtomMolList.
std::vector<int> donorHsFromMolList(const std::vector<std::vector<int>> &molIDlist) {
  std::vector<int> donorHs;
  for (const auto &row : molIDlist) {
    if (row.size() < 3) {
      continue;
    }
    donorHs.push_back(row[1]);
    donorHs.push_back(row[2]);
  }
  return donorHs;
}

} // namespace

bool bond::donatedHydrogenBond(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const molSys::PointCloud<molSys::Point<double>, double> &hCloud,
    int acceptorIndex, int donorIndex, const std::vector<int> &donorHs,
    double distCutoff, double angleCutoff) {
  if (acceptorIndex < 0 || donorIndex < 0 ||
      acceptorIndex >= yCloud.nop || donorIndex >= yCloud.nop) {
    return false;
  }
  const auto oo = gen::relDist(yCloud, acceptorIndex, donorIndex);
  const std::vector<double> ooVec{oo[0], oo[1], oo[2]};
  const double dist2 = distCutoff * distCutoff;
  for (int hAtomIndex : donorHs) {
    if (hAtomIndex < 0 || hAtomIndex >= hCloud.nop) {
      continue;
    }
    const auto oh = gen::relDistFromPoint(
        yCloud, acceptorIndex, hCloud.pts[static_cast<std::size_t>(hAtomIndex)].x,
        hCloud.pts[static_cast<std::size_t>(hAtomIndex)].y, hCloud.pts[static_cast<std::size_t>(hAtomIndex)].z);
    const double oh2 = oh[0] * oh[0] + oh[1] * oh[1] + oh[2] * oh[2];
    if (oh2 >= dist2) {
      continue;
    }
    const std::vector<double> ohVec{oh[0], oh[1], oh[2]};
    if (gen::radDeg(gen::eigenVecAngle(ooVec, ohVec)) <= angleCutoff) {
      return true;
    }
  }
  return false;
}

/**
 * @details Create a vector of vectors containing bond information (outputs
 * bonded atom IDs, not indices!) from the neighbour list vector of vectors
 * (which contains atom INDICES). Moreover, the first element corresponds to the
 * atom whose neighbours have been found.
 */
std::vector<std::vector<int>> bond::populateBonds(
    const std::vector<std::vector<int>> &nList,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  //
  std::vector<std::vector<int>> bonds; // Output vector of vectors
  std::vector<int> currentBond;        // Vector for the current bond
  int first, second;                   // Indices for the bonds
  int iatom, jatom;                    // Elements of the bond
  int iatomID, jatomID; // Atom IDs of the elements which are bonded

  // Error handling
  if (nList.empty()) {
    // There is some problem!
    std::cerr << "There are no bonds in the system!\n";
    return bonds;
  }

  // Form of the bonds vector of vectors:
  // 214    272
  // 1       2
  // Meaning that 272 and 214 are bonded; similarly 1 and 2 are bonded
  // To avoid repitition, discard bonds whose first value is greater than the
  // second value

  // Traverse the neighbour list

  // Loop through every atom in the neighbour list by index
  for (int i = 0; i < nList.size(); i++) {
    iatom = nList[static_cast<std::size_t>(i)][0]; // Index of the i^th atom
    // Get the neighbours of iatom
    for (int j = 1; j < nList[static_cast<std::size_t>(i)].size(); j++) {
      //
      jatom = nList[static_cast<std::size_t>(iatom)][static_cast<std::size_t>(j)]; // Index of the neighbour
      // To avoid duplicates, skip all bonds such
      // that iatom>jatom
      if (iatom > jatom) {
        continue;
      } // Skip to avoid duplicates

      // Clear the current bond vector
      currentBond.clear();
      // Fill the current bond vector with ATOM IDs
      iatomID = yCloud.pts[static_cast<std::size_t>(iatom)].atomID;
      jatomID = yCloud.pts[static_cast<std::size_t>(jatom)].atomID;
      currentBond.push_back(iatomID);
      currentBond.push_back(jatomID);
      // Add to the bonds vector of vectors
      bonds.push_back(currentBond);
    } // end of loop the neighbour list
    // The last pair is with the last element and the first element
    // Fill currentBond and update bonds

  } // end of loop through rings

  return bonds;
}

/**
 *  @details Create a vector of vectors containing bond information (outputs
 * bonded atom IDs, not indices!) from the neighbour list vector of vectors
 * (which contains atom INDICES). Bonds between dummy atoms, and between dummy
 * and ice atoms are not added. Moreover, the first element corresponds to the
 * atom whose neighbours have been found.
 *  @param[in] nList Row-ordered neighbour list by ID
 *  @param[in] yCloud The input molSys::PointCloud
 *  @param[in] atomTypes Contains an atom type for each particle in yCloud
 */
std::vector<std::vector<int>> bond::populateBonds(
    const std::vector<std::vector<int>> &nList,
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<cage::iceType> &atomTypes) {
  //
  std::vector<std::vector<int>> bonds; // Output vector of vectors
  std::vector<int> currentBond;        // Vector for the current bond
  int first, second;                   // Indices for the bonds
  int iatom, jatom;                    // Elements of the bond
  int iatomID, jatomID; // Atom IDs of the elements which are bonded

  // Error handling
  if (nList.empty()) {
    // There is some problem!
    std::cerr << "There are no bonds in the system!\n";
    return bonds;
  }

  // Form of the bonds vector of vectors:
  // 214    272
  // 1       2
  // Meaning that 272 and 214 are bonded; similarly 1 and 2 are bonded
  // To avoid repitition, discard bonds whose first value is greater than the
  // second value

  // Traverse the neighbour list

  // Loop through every atom in the neighbour list by index
  for (int i = 0; i < nList.size(); i++) {
    iatom = nList[static_cast<std::size_t>(i)][0]; // Index of the i^th atom
    // Skip for dummy atoms
    if (atomTypes[static_cast<std::size_t>(iatom)] == cage::iceType::dummy) {
      continue;
    } // Skip for dummy atoms
    // Get the neighbours of iatom
    for (int j = 1; j < nList[static_cast<std::size_t>(i)].size(); j++) {
      //
      jatom = nList[static_cast<std::size_t>(iatom)][static_cast<std::size_t>(j)]; // Index of the neighbour
      // Skip for dummy atoms
      if (atomTypes[static_cast<std::size_t>(jatom)] == cage::iceType::dummy) {
        continue;
      } // Skip for dummy atoms
      // To avoid duplicates, skip all bonds such
      // that iatom>jatom
      if (iatom > jatom) {
        continue;
      } // Skip to avoid duplicates

      // Clear the current bond vector
      currentBond.clear();
      // Fill the current bond vector with ATOM IDs
      iatomID = yCloud.pts[static_cast<std::size_t>(iatom)].atomID;
      jatomID = yCloud.pts[static_cast<std::size_t>(jatom)].atomID;
      currentBond.push_back(iatomID);
      currentBond.push_back(jatomID);
      // Add to the bonds vector of vectors
      bonds.push_back(currentBond);
    } // end of loop the neighbour list
    // The last pair is with the last element and the first element
    // Fill currentBond and update bonds

  } // end of loop through rings

  return bonds;
}

/********************************************/ /**
 *  Create a vector of vectors (similar to the neighbour list conventions). The
 output vector of vectors is row-ordered. The first element is the atom ID of
 the particle for which the neighbours are enumerated (the central atom),
 followed by the central atom's neighbour's IDs (not indices). Decides the
 existence of the hydrogen bond depending on the O--O and O--H vectors from the
 neighbour list (by ID) already constructed.
 *  @param[in] filename Filename of the trajectory, with the hydrogen and oxygen
 coordinates
 *  @param[in] yCloud The input molSys::PointCloud for the oxygen atoms only
 *  @param[in] nList Row-ordered neighbour list by atom ID
 *  @param[in] targetFrame The target or current frame number (starts from 1)
 and is not the timestep value
 *  @param[in] Htype The type ID of the hydrogen atoms
 ***********************************************/
std::vector<std::vector<int>>
bond::populateHbonds(std::string filename,
                     molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     const std::vector<std::vector<int>> &nList,
                     int targetFrame, int Htype, double distCutoff,
                     double angleCutoff) {
  // Read hydrogen atoms from the trajectory file
  molSys::PointCloud<molSys::Point<double>, double> hCloud;
  hCloud = sinp::readLammpsTrjreduced(filename, targetFrame, hCloud, Htype);
  // Delegate to the cloud-based overload
  return bond::populateHbondsWithInputClouds(yCloud, hCloud, nList, distCutoff,
                                             angleCutoff);
}

/********************************************/ /**
 *  Create a vector of vectors (similar to the neighbour list conventions). The
 output vector of vectors is row-ordered. The first element is the atom ID of
 the particle for which the neighbours are enumerated (the central atom),
 followed by the central atom's neighbour's IDs (not indices). Decides the
 existence of the hydrogen bond depending on the O--O and O--H vectors from the
 neighbour list (by ID) already constructed. The only difference between this
 function and the similar populateHbonds function is that this takes in the H
 atom PointCloud as input
 *  @param[in] yCloud The input molSys::PointCloud for the hydrogen atoms only,
 for the entire system (regardless of whether there is a slice or not)
 *  @param[in] hCloud The input molSys::PointCloud for the oxygen atoms only
 *  @param[in] nList Row-ordered neighbour list by atom ID
 ***********************************************/
std::vector<std::vector<int>> bond::populateHbondsWithInputClouds(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    molSys::PointCloud<molSys::Point<double>, double> &hCloud,
    const std::vector<std::vector<int>> &nList, double distCutoff,
    double angleCutoff) {
  auto molIDlist = molSys::hAtomMolList(hCloud, yCloud);
  return bond::populateHbondsFromDonors(yCloud, hCloud, nList,
                                        donorHsFromMolList(molIDlist),
                                        distCutoff, angleCutoff);
}

std::vector<std::vector<int>> bond::populateHbondsFromDonors(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    molSys::PointCloud<molSys::Point<double>, double> &hCloud,
    const std::vector<std::vector<int>> &nList,
    const std::vector<int> &donorHs, double distCutoff, double angleCutoff) {
  std::vector<std::vector<int>> hBondNet;
  std::unordered_map<int, std::vector<int>> donorHsByMol;
  for (int hAtomIndex : donorHs) {
    if (hAtomIndex < 0 || hAtomIndex >= hCloud.nop) {
      continue;
    }
    donorHsByMol[hCloud.pts[static_cast<std::size_t>(hAtomIndex)].molID].push_back(hAtomIndex);
  }

  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    hBondNet.push_back(std::vector<int>());
    hBondNet[static_cast<std::size_t>(iatom)].push_back(yCloud.pts[static_cast<std::size_t>(iatom)].atomID);
  }

  const std::vector<int> emptyHs;
  auto hsFor = [&](int donorIndex) -> const std::vector<int> & {
    auto it = donorHsByMol.find(yCloud.pts[static_cast<std::size_t>(donorIndex)].molID);
    if (it == donorHsByMol.end()) {
      return emptyHs;
    }
    return it->second;
  };

  for (size_t iatom = 0; iatom < nList.size(); iatom++) {
    if (nList[iatom].empty()) {
      continue;
    }
    const int iatomID = nList[iatom][0];
    const int nnumNeighbours = static_cast<int>(nList[iatom].size()) - 1;
    const int iatomIndex = static_cast<int>(iatom);
    for (int j = 1; j <= nnumNeighbours; j++) {
      const int jatomID = nList[iatom][static_cast<std::size_t>(j)];
      if (iatomID >= jatomID) {
        continue;
      }
      auto gotJ = yCloud.idIndexMap.find(jatomID);
      if (gotJ == yCloud.idIndexMap.end()) {
        continue;
      }
      const int jatomIndex = gotJ->second;
      if (bond::donatedHydrogenBond(yCloud, hCloud, iatomIndex, jatomIndex,
                                    hsFor(jatomIndex), distCutoff,
                                    angleCutoff) ||
          bond::donatedHydrogenBond(yCloud, hCloud, jatomIndex, iatomIndex,
                                    hsFor(iatomIndex), distCutoff,
                                    angleCutoff)) {
        hBondNet[static_cast<std::size_t>(iatomIndex)].push_back(jatomID);
        hBondNet[static_cast<std::size_t>(jatomIndex)].push_back(iatomID);
      }
    }
  }

  return hBondNet;
}

namespace {

double length3(const std::array<double, 3> &v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

} // namespace

std::vector<bond::DirectedHydrogenBond> bond::directedHydrogenBonds(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    double pairCutoff, double distCutoff, double angleCutoff,
    const std::vector<int> &donorTypes, const std::vector<int> &acceptorTypes,
    int hydrogenType, bool excludeSameMol) {
  std::vector<DirectedHydrogenBond> edges;
  if (donorTypes.empty() || acceptorTypes.empty() || hydrogenType <= 0) {
    return edges;
  }
  if (!(pairCutoff > 0.0) || !(distCutoff > 0.0) || !(angleCutoff > 0.0) ||
      !std::isfinite(pairCutoff) || !std::isfinite(distCutoff) ||
      !std::isfinite(angleCutoff)) {
    return edges;
  }
  const std::unordered_set<int> donors(donorTypes.begin(), donorTypes.end());
  const std::unordered_set<int> acceptors(acceptorTypes.begin(), acceptorTypes.end());
  if (donors.count(hydrogenType) != 0 || acceptors.count(hydrogenType) != 0) {
    return edges;
  }
  std::unordered_map<int, std::vector<int>> heaviesByMol;
  std::vector<int> hydrogens;
  hydrogens.reserve(static_cast<std::size_t>(std::max(cloud.nop, 0)));
  for (int i = 0; i < cloud.nop; ++i) {
    const auto &point = cloud.pts[static_cast<std::size_t>(i)];
    if (point.type == hydrogenType) {
      hydrogens.push_back(i);
      continue;
    }
    heaviesByMol[point.molID].push_back(i);
  }
  std::vector<std::vector<int>> hydrogensOnDonor(static_cast<std::size_t>(std::max(cloud.nop, 0)));
  for (const int hydrogen : hydrogens) {
    const auto &hpoint = cloud.pts[static_cast<std::size_t>(hydrogen)];
    const auto found = heaviesByMol.find(hpoint.molID);
    if (found == heaviesByMol.end()) {
      continue;
    }
    int best = -1;
    double best2 = 0.0;
    for (const int heavy : found->second) {
      const auto dr = gen::relDist(cloud, heavy, hydrogen);
      const double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
      if (best < 0 || d2 < best2) {
        best = heavy;
        best2 = d2;
      }
    }
    if (best >= 0 &&
        donors.count(cloud.pts[static_cast<std::size_t>(best)].type) != 0) {
      hydrogensOnDonor[static_cast<std::size_t>(best)].push_back(hydrogen);
    }
  }
  const auto neighbours = nneigh::getNewNeighbourListByIndex(cloud, pairCutoff);
  const double pair2 = pairCutoff * pairCutoff;
  const double dist2 = distCutoff * distCutoff;
  for (int donor = 0; donor < cloud.nop; ++donor) {
    const auto &owned = hydrogensOnDonor[static_cast<std::size_t>(donor)];
    if (owned.empty() || donor >= static_cast<int>(neighbours.size())) {
      continue;
    }
    const auto &dpoint = cloud.pts[static_cast<std::size_t>(donor)];
    const auto &row = neighbours[static_cast<std::size_t>(donor)];
    for (const int acceptor : row) {
      if (acceptor == donor || acceptor < 0 || acceptor >= cloud.nop) {
        continue;
      }
      const auto &apoint = cloud.pts[static_cast<std::size_t>(acceptor)];
      if (apoint.type == hydrogenType || acceptors.count(apoint.type) == 0) {
        continue;
      }
      if (excludeSameMol && apoint.molID == dpoint.molID) {
        continue;
      }
      const auto oo = gen::relDist(cloud, acceptor, donor);
      const double oo2 = oo[0] * oo[0] + oo[1] * oo[1] + oo[2] * oo[2];
      if (!(oo2 > 0.0) || oo2 >= pair2 || !std::isfinite(oo2)) {
        continue;
      }
      const std::vector<double> ooVec{oo[0], oo[1], oo[2]};
      for (const int hydrogen : owned) {
        const auto &hpoint = cloud.pts[static_cast<std::size_t>(hydrogen)];
        const auto oh = gen::relDistFromPoint(cloud, acceptor, hpoint.x, hpoint.y, hpoint.z);
        const double oh2 = oh[0] * oh[0] + oh[1] * oh[1] + oh[2] * oh[2];
        if (!(oh2 > 0.0) || oh2 >= dist2 || !std::isfinite(oh2)) {
          continue;
        }
        const std::vector<double> ohVec{oh[0], oh[1], oh[2]};
        const double angle = gen::radDeg(gen::eigenVecAngle(ooVec, ohVec));
        if (!std::isfinite(angle) || angle > angleCutoff) {
          continue;
        }
        DirectedHydrogenBond edge;
        edge.donor = donor;
        edge.hydrogen = hydrogen;
        edge.acceptor = acceptor;
        edge.d_AH = length3(oh);
        edge.d_DA = length3(oo);
        edge.angle = angle;
        edges.push_back(edge);
      }
    }
  }
  std::sort(edges.begin(), edges.end(),
            [](const DirectedHydrogenBond &a, const DirectedHydrogenBond &b) {
              if (a.donor != b.donor) {
                return a.donor < b.donor;
              }
              if (a.hydrogen != b.hydrogen) {
                return a.hydrogen < b.hydrogen;
              }
              return a.acceptor < b.acceptor;
            });
  return edges;
}

/**
*  Calculates the bond length between a Hydrogen and Oxygen
 atom of two different atoms, given their respective pointClouds and the indices
to each atom.
 *  @param[in] oCloud The molSys::PointCloud for the oxygen atoms only
 *  @param[in] hCloud The molSys::PointCloud for the hydrogen atoms only
 *  @param[in] oAtomIndex The index (in the oCloud) of the oxygen atom
 *  @param[in] hAtomIndex The index (in the hCloud) of the hydrogen atom
*/
double bond::getHbondDistanceOH(
    const molSys::PointCloud<molSys::Point<double>, double> &oCloud,
    const molSys::PointCloud<molSys::Point<double>, double> &hCloud,
    int oAtomIndex, int hAtomIndex) {
  const std::vector<double> hpt = {hCloud.pts[static_cast<std::size_t>(hAtomIndex)].x,
                                   hCloud.pts[static_cast<std::size_t>(hAtomIndex)].y,
                                   hCloud.pts[static_cast<std::size_t>(hAtomIndex)].z};
  return gen::unWrappedDistFromPoint(oCloud, oAtomIndex, hpt);
}

/**
 *  Create a vector of vectors containing bond information (bonded atom IDs, not
 vector or array indices!) from the ring vector of vectors and cageList
 *  @param[in] rings Row-ordered vector of vectors atom indices of ring
 information. Each row is a ring, containing the indices of the particles in
 that ring
 *  @param[in] cageList A vector of cage::Cage containing a list of HCs or DDCs
 *  @param[in] type The type of cage to get bonds for
 *  @param[in, out] nRings The total number of rings for all the cages, for the
 particular cage type
*/
std::vector<std::vector<int>>
bond::createBondsFromCages(const std::vector<std::vector<int>> &rings,
                           std::vector<cage::Cage> &cageList,
                           cage::cageType type, int &nRings) {
  std::vector<std::vector<int>> bonds; // Output vector of vectors
  std::vector<int> currentBond;        // Vector for the current bond

  // Error handling. rings[0] is not valid until this returns.
  if (rings.empty()) {
    // There is some problem!
    std::cerr << "There are no rings in the system!\n";
    return bonds;
  }
  // The closing bond of every ring uses the size of the first ring.
  const std::size_t ringSize = rings[0].size();

  // Form of the bonds vector of vectors:
  // 272    214
  // 1       2
  // Meaning that 272 and 214 are bonded; similarly 1 and 2 are bonded

  // Traverse the cageList vector of Cages

  nRings = 0; // init

  // Loop through all the cages
  for (const cage::Cage &cage : cageList) {
    // Skip if the cage is of a different type
    if (cage.type != type) {
      continue;
    }
    nRings += seams::checkedCast<int>(cage.rings.size()); // Add to the number of rings
    //
    // Now loop through a particular ring inside the i^th cage
    for (const int ringIndex : cage.rings) {
      const std::vector<int> &ring = rings[seams::toIndex(ringIndex)];
      // Get the first atom of each pair inside the ring
      for (std::size_t k = 0; k + 1 < ring.size(); k++) {
        // Clear the current bond vector
        currentBond.clear();
        // Fill the current bond vector
        currentBond.push_back(ring[k]);
        currentBond.push_back(ring[k + 1]);
        std::sort(currentBond.begin(), currentBond.end());
        // Add to the bonds vector of vectors
        bonds.push_back(currentBond);
      } // end of loop through ring elements, except the last one
      // The last pair is with the last element and the first element
      // Fill currentBond and update bonds
      currentBond.clear();
      currentBond.push_back(ring[ringSize - 1]);
      currentBond.push_back(ring[0]);
      std::sort(currentBond.begin(), currentBond.end());
      bonds.push_back(currentBond);
    } // end of loop through a particular ring
  } // end of loop through cages

  // This may have duplicates, so the duplicates should be removed
  std::sort(bonds.begin(), bonds.end());
  bonds.erase(std::unique(bonds.begin(), bonds.end()), bonds.end());

  return bonds;
}

/**
 *  The bond 1 2 and 2 1 are the same. To prevent multiple bonds between the
 same atoms, remove all bonds which are duplicates of the reversed vectors
 (denoting individual bonds) within the bonds vector of vectors
 *  @param[in, out] bonds Row-ordered vector of vectors of the bond matrix
 Each row is a ring, containing the indices of the particles in that ring
*/
std::vector<std::vector<int>>
bond::trimBonds(std::vector<std::vector<int>> bonds) {
  std::vector<int>
      reversedBond; // Vector for the current bond in reversed order
  std::vector<bool> isBondFlag;
  int temp = 0;

  std::sort(bonds.begin(), bonds.end());
  bonds.erase(std::unique(bonds.begin(), bonds.end()), bonds.end());

  // // Temp uncommented
  // // Resize the flag vector and init to true
  // isBondFlag.resize(bonds.size(), true);
  // // Form of the bonds vector of vectors:
  // // 272    214
  // // 1       2
  // // Meaning that 272 and 214 are bonded; similarly 1 and 2 are bonded

  // // Traverse the bonds vector of vectors

  // // Loop through all the bonds
  // for(int ibond=0; ibond<bonds.size(); ibond++){

  //   // Skip if the bond has been flagged as false
  //   if(isBondFlag[ibond] == false){continue;}

  //   reversedBond.clear();
  //   reversedBond.push_back( bonds[ibond][1] );
  //   reversedBond.push_back( bonds[ibond][0] );

  //   // Compare with other bonds by looping through
  //   // the rest
  //   for(int jbond=0; jbond<bonds.size(); jbond++){
  //     // Skip if jbond has been flagged
  //     if(isBondFlag[jbond] == false){continue;}
  //     // Compare reversed bond and jbond if jbond has not been flagged yet
  //     if(reversedBond==bonds[jbond]){
  //       isBondFlag[jbond] = false;
  //       temp++;
  //     } // end of check of comparison
  //   } // end of looping through all possible bonds
  // } // end of loop through all possible bonds
  // // temp

  return bonds;
}
