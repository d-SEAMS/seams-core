<!-- https://docs.dseams.info/api/namespace_traj.md -->

# namespace traj {#namespace-traj}

## Overview {#overview}

Multi-frame PBC unwrapping, molecular orientation/dipole vectors, and streaming time-correlation accumulators. [More…](#details-namespacetraj)

## Detailed Documentation {#detailed-documentation}

Multi-frame PBC unwrapping, molecular orientation/dipole vectors, and streaming time-correlation accumulators.

### Global Variables {#global-variables}

Largest particle or vector count n whose 3 * n doubles one std::vector can hold.

Larger counts are refused, so 3 * n never wraps.

Largest box count n whose 9 * n doubles fit in a single array.

Larger counts are refused, so 9 * n and pointer offsets never overflow.

Largest max_lag whose max_lag + 1 lag bins one std::vector can hold.

Larger lags are refused, so max_lag + 1 never wraps to zero.

### Global Functions {#global-functions}

out = a * b; false, leaving out alone, when that overflows size_t.

Unwrap a multi-frame coordinate array \[nFrames * 3 * nAtoms\].

**Parameters:**

| boxes | Either 1 cell (9 doubles) or nFrames cells (9 * nFrames doubles). |
| --- | --- |

Compute intact PBC-assembled per-molecule centres for one frame.

Molecules appear in first-seen order of molIds. If masses is nullptr, equal atom weights are used.

Cloud overload using stored atom coordinates and molIDs.

Compute and unwrap molecule centres across nFrames frames.

Extract per-water O->midpoint(H1,H2) and H1->H2 vectors under MIC.

Pairs each oxygen (oxygenType) with its two hydrogens (hydrogenType) by shared molID when multi-atom molecule IDs exist, falling back to the two nearest hydrogens within covalentCutoff Angstrom when molIDs do not group O and H atoms.

Compute per-molecule dipole or geometric orientation vectors.

When charges is non-null, computes sum_a q_a (r_a^intact - R_geom) in e*A (origin-independent for neutral molecules and translation-invariant for ions). When charges is null, computes the mean MIC bond vector from the first atom of each molecule to its remaining atoms (reducing to O->H-midpoint for O-H-H).
