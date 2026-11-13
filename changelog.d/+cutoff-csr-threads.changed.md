The orthorhombic cutoff neighbour list is built by a parallel cell sort.
Each thread counts its atoms, a prefix sum places them, and the pairs are
written straight into one CSR. The same buffers are reused on the next frame.
The pair set does not depend on the thread count.
q_l reads that CSR. The neighbour average is skipped when only q_l is requested.
