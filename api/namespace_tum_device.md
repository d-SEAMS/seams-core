<!-- https://docs.dseams.info/api/namespace_tum_device.md -->

# namespace tum::device {#namespace-tum-device}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Variables {#global-variables}

Longest atom row the kernels accept.

Candidate buffers are sized by it and ddcFrom keeps one bit per entry in an unsigned. A liquid atom of a 3.5 A cutoff graph can lie in more than 16 six-rings; a graph with a longer row goes to the host (`tum_offload.hpp`).

### Global Functions {#global-functions}

Graph distance from a to b when it is at most cap, else -1.

cap is 1 or 2: the Franzblau reject bound for a six-ring.

Hop-bounded Franzblau SP test on a six-cycle of the bond graph.

beta(b1, b2) for rings sharing no atom: some alternating triplet {l_p, l_p+2, l_p+4} of b1 bonds one-to-one onto an alternating triplet {m_q, m_q+2, m_q+4} of b2, and no other atom of b2 is bonded to b1.

All 2 x 2 x 6 alignments are tried, so rotation, direction and ring order do not matter.

Enumerate hop-bound primitive six-rings whose lowest vertex is i.

nRings counts every ring found; only the first maxRings are stored.

throughCount counts every stored ring through the atom; only the first maxPer ids are kept.

True when ring j is in the through row of a bond neighbour of bi\[0\] or bi\[1\] visited before (slot, n), so each candidate is tested once.

Ordered basal pairs (i, j).

Every pair passing beta has a bond from bi\[0\] or bi\[1\] into bj, so the neighbours of those two atoms reach every candidate.

epsilon(i) for a ring that is not HC, as [ring::equatorialConditions](namespace_ring.md#doxid-namespacering-1a01243f4581fb55b417184d98f9ac70aa) : every atom in 3 other rings, and pairwise distinct P_k in T_k with P_0, P_2, P_4 sharing an atom, P_1, P_3, P_5 sharing an atom and, for some starting atom and direction of i, the alternate pairs sharing 3 atoms.

On a pass marks i and every P_k of every valid assignment.
