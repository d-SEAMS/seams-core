# Triclinic cells

Pinned ABAB, 20 samples. Serial is pinned to one core. Eight threads are pinned to cores 0-7. The machine has 8 cores. OpenMP uses close binding and a passive wait, and the process affinity is restored after libraries import.

Tilted ice Ih is the orthorhombic generator with LAMMPS tilt factors 0.15 ly, 0.05 lz and 0.08 lz, at 4096, 32768 and 110592 molecules. sH is the genice triclinic cell. 136 does not divide those counts, so the supercell is the nearest multiple whose edges stay compact: 2x5x3 = 4080, 4x10x6 = 32640, 6x15x9 = 110160. The raw log is benchmarks/h2h-triclinic.txt.

Cutoff pairs match freud, OVITO and mdapy exactly on every frame. k nearest neighbours (k=4) match atom for atom. q6 and q6bar match OVITO and mdapy to about 1e-14. Freud's q6, from a float32 cell, differs by at most 3.3e-5 on the largest sH cell. CHILL+ in ovito mode matches the OVITO modifier and mdapy on every atom. Paper mode, which is the published four-neighbour average, gives the same label as ovito mode on these crystals because every oxygen's cutoff shell is those four neighbours. Primitive rings match rings-statistics-matter as sets, not only as counts.

Every timed comparison wins on the median and on the minimum. At one core and 110592 tilted ice molecules the cutoff is 7.87 ms (IQR 0.08, minimum 7.74) against OVITO CreateBonds at 41.07 ms (IQR 0.38, minimum 40.70). The same cutoff at eight threads is faster than freud, OVITO and mdapy on every cell, including the 4080-molecule sH cell.
