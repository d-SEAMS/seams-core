# PLUMED 2.10

Pinned ABAB, 20 samples. Serial is pinned to one core. Eight threads are pinned to cores 0-7. The machine has 8 cores. OpenMP uses close binding and a passive wait, and the process affinity is restored after the libraries load. The raw log is benchmarks/h2h-plumed.txt.

Ice Ih at 4096, 32768 and 110592 molecules. PLUMED 2.10.0 COORDINATIONNUMBER and Q6 use `SWITCH={CUSTOM FUNC=step(3.5-x) R_0=1 D_MAX=4.2}`. The step is the same 3.5 A neighbour cutoff as d-SEAMS. `D_MAX` is only the neighbour-list skin: on this crystal the first oxygen shell is at 2.76 A and the next is at 4.51 A, so 4.2 A contains those four neighbours and no others. Without `D_MAX`, PLUMED builds a dense adjacency matrix, which matched the 4096-molecule cell and then faulted inside `AdjacencyMatrixBase::updateNeighbourList` at 32768. Q6 is scaled by `sqrt(4*pi/13)` so it is the local Steinhardt parameter.

Coordination matches the cutoff degree on every atom (maximum absolute difference 0). Local q6 matches to at most 5.6e-16.

Every timed comparison wins on the median and on the minimum. At one core and 110592 molecules, coordination is 9.29 ms (IQR 1.16, minimum 8.68) against PLUMED at 272.09 ms (IQR 2.73, minimum 270.35), and q6 is 34.86 ms (IQR 1.09, minimum 29.79) against 2814.60 ms (IQR 6.83, minimum 2799.57). At eight threads the same q6 is 16.92 ms (IQR 2.08, minimum 14.31) against 2799.48 ms (IQR 16.21, minimum 2783.21).

## In-tree colvar

`engines/plumed/SeamsColvar.cpp` builds as a PLUMED action named SEAMS. The shared library is produced next to the release build and is not part of the tree. PLUMED passes nanometres; the action scales positions and the cell to angstroms and calls the C API. Derivatives stay off unless a restraint allocates them. The raw log is benchmarks/h2h-plumed-plugin.txt.

On the same three ice cells the mean q6, the CHILL+ cubic fraction, the primitive six-ring fraction, the hexagonal count and the six-ring count match the C API. The largest mean-q6 gap is 5.4e-13 at 110592 molecules. Ideal ice is entirely hexagonal, so the cubic fraction is 0 and every molecule sits on a six-ring.

The timed plugin call is one C++ evaluation. The API arm is the same work as three ctypes calls (q6, CHILL+, rings). The plugin is faster on that pairing: at eight threads and 110592 molecules the API median is 95.49 ms (IQR 2.76, minimum 91.82) and the plugin median is 63.33 ms (IQR 1.94, minimum 60.00). That gap is the ctypes boundary, not a second implementation.
