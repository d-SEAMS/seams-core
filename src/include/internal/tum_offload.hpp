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
 *  and omp_get_num_devices() > 0, otherwise on the host. The seeded
 *  CLI path uses that enumerator on the union four-nearest graph.
 *  CHILL+ stays on the host. Steinhardt offload is unchanged.
 *
 *  The enumerator writes into fixed buffers: a row of ring ids per atom,
 *  a ring table and a basal-pair table. Every counter is incremented
 *  before the bound check, so a pass reports the size it needed. A pass
 *  that overflowed any buffer is discarded and rerun with that buffer
 *  grown to the reported size. A row longer than tum::device::kRowCap
 *  does not fit the kernels' local arrays, and that graph is evaluated
 *  by hostCageCounts instead.
 */

namespace tum {

struct CageCounts {
  int nSix = 0;
  int nHcRings = 0;
  int nDdcRings = 0;
  int nHcAtoms = 0;
  int nDdcAtoms = 0;
  int retries = 0;           // passes discarded after a buffer overflowed
  bool hostFallback = false; // a row passed kRowCap; hostCageCounts ran
  bool usedDevice = false;
  std::vector<int> atomSix; // six-rings through each atom
  std::vector<int> atomHc;
  std::vector<int> atomDdc;
};

/** Initial buffer sizes of the hop-bound path. Each grows on overflow. */
struct SpecializedCaps {
  int rowCap = 32; // ring ids per atom row, at most tum::device::kRowCap
  int ringCap = 0; // rings; 0 is 16 per atom
  int pairCap = 0; // basal pairs; 0 is 8 per ring of ringCap
  // Run the OpenMP target branch when no device is present; it then
  // executes on the initial device. Has no effect without offload.
  bool forceTarget = false;
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

/** As above, starting from the given buffer sizes. */
[[nodiscard]] CageCounts
specializedCageCounts(const std::vector<std::vector<int>> &nList,
                      const SpecializedCaps &caps);

/** preferOffload() selects specializedCageCounts, otherwise hostCageCounts. */
[[nodiscard]] CageCounts cageCounts(const std::vector<std::vector<int>> &nList);

} // namespace tum

#endif
