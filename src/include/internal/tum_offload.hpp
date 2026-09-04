//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#ifndef SEAMS_TUM_OFFLOAD_H_
#define SEAMS_TUM_OFFLOAD_H_

#include <vector>

/** @file tum_offload.hpp
 *  @brief TUM ice score: hop-bound six-rings and cage affiliation.
 *
 *  SEAMS_OFFLOAD=0 (or unset when no device) uses host Franzblau
 *  ringNetwork plus cageAffiliation. Any other SEAMS_OFFLOAD value
 *  selects the hop-bound six-ring enumerator and the same affiliation
 *  predicates, on the OpenMP target device when the build provides one
 *  and omp_get_num_devices() > 0, otherwise on the host. CHILL+ stays
 *  on the host. Steinhardt offload is unchanged.
 */

namespace tum {

struct CageCounts {
  int nSix = 0;
  int nHcRings = 0;
  int nDdcRings = 0;
  int nHcAtoms = 0;
  int nDdcAtoms = 0;
  int ringsDropped = 0;
  bool usedDevice = false;
  std::vector<int> atomHc;
  std::vector<int> atomDdc;
};

/** True when SEAMS_OFFLOAD is set and is not "0". */
[[nodiscard]] bool preferOffload();

/** Host Franzblau six-rings plus cageAffiliation on an index nList
 *  (leading self entry, as neighbourListByIndex produces). */
[[nodiscard]] CageCounts hostCageCounts(
    const std::vector<std::vector<int>> &nList);

/** Hop-bound six-ring enumerator plus affiliation predicates. Uses the
 *  OpenMP target device when offload is compiled and a device exists. */
[[nodiscard]] CageCounts specializedCageCounts(
    const std::vector<std::vector<int>> &nList);

/** preferOffload() selects specializedCageCounts, otherwise hostCageCounts. */
[[nodiscard]] CageCounts cageCounts(const std::vector<std::vector<int>> &nList);

} // namespace tum

#endif
