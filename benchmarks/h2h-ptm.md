# Polyhedral template matching

Pinned ABAB, 20 samples. Serial is pinned to one core. Eight threads are pinned to cores 0-7. The machine has 8 cores. OpenMP uses close binding and a passive wait, and the process affinity is restored after the libraries load. The raw log is benchmarks/h2h-ptm.txt.

No remote branch holds the PTM/IRA tree from cloud agent bc-035d2d78. The standalone PTM library is not installed. The reference is OVITO 3.17.2 polyhedral template matching with the cubic-diamond and hexagonal-diamond structures enabled and `rmsd_cutoff` 0.1. d-SEAMS is `seams_system_symmetry`. OVITO codes 6 and 7 are cubic and hexagonal diamond; every other code is other. d-SEAMS codes are hexagonal, cubic, stacking, and other.

Ideal ice Ih at 4096 and 32768 molecules is hexagonal on every atom in both codes. Ideal ice Ic at the same sizes is cubic on every atom. The 512-molecule stacking frame agrees on all 512 atoms: 128 cubic, 256 hexagonal, 128 other.

The 512-molecule interface agrees on 492 atoms (128 hexagonal on both, 364 other on both, 9 cubic that PTM leaves as other, 11 hexagonal that PTM leaves as other). The melted frame agrees on 477 atoms (315 hexagonal on both, 162 other on both, 32 hexagonal that PTM leaves as other, 2 other that PTM calls hexagonal, and 1 stacking atom that PTM leaves as other). PTM matches a local diamond template. d-SEAMS labels the core environment. Those leftovers are that definition, not a threshold.

At one core, ice Ih at 32768 molecules is 19.99 ms (IQR 1.40, minimum 18.18) against OVITO PTM at 852.61 ms (IQR 16.75, minimum 840.11). Ice Ic at the same size is 21.48 ms (IQR 1.58, minimum 19.70) against 953.78 ms (IQR 7.55, minimum 932.84). At eight threads those two cells are 13.63 ms (IQR 0.44, minimum 13.12) against 113.51 ms (IQR 3.67, minimum 109.99), and 14.11 ms (IQR 0.38, minimum 13.67) against 122.26 ms (IQR 5.60, minimum 119.35). Every comparison wins on the median and on the minimum.
