# d-SEAMS, freud 3.6.1, and OVITO 3.17.2

One process, eight samples, eight cores. Orthorhombic ice Ih at 4096, 32768, and 110592 oxygens, plus the mixed interface, stacking, and melted frames (512 each). Medians are milliseconds. The raw lines are `benchmarks/h2h-freud-ovito.txt`. CHILL+ is the ovito definition (every neighbour inside 3.5 A). Freud has no CHILL+ and no diamond identifier. OVITO's polyhedral template matching leaves the diamond templates off by default; they were enabled here, with an RMSD cutoff of 0.1, because that is the diamond counterpart of `seams_system_symmetry`.

## Agreement

Cutoff pairs match exactly on every frame, against freud and against OVITO CreateBonds. Counts are 8192, 65536, and 221184 on the ice cells, and 1208, 960, and 1085 on interface, stacking, and melted.

On the ice cells the four nearest neighbours match freud and OVITO exactly. On the mixed frames d-SEAMS stores the symmetrized one-way graph: an edge is kept when either oxygen lists the other among its four nearest inside 3.5 + 1.5 A, so a row can hold more than four neighbours. Freud and OVITO return exactly four directed neighbours. Taking the four shortest bonds of each d-SEAMS row, the remaining swaps sit 1e-15 A apart. That is a tie between the fourth and fifth neighbour, not a different shell. The swap counts are 7 (interface), 55 (stacking), and 10 (melted).

q3, q6, and both Lechner-Dellago averages match OVITO's Steinhardt modifier to 2e-15 on every frame. The modifier's columns, with `max_l = 6`, are q3 at index 1 and q6 at index 4. Freud stores `ql` as float32. On the ice cells the largest absolute gap is 3.4e-6 for q6 and 3.0e-6 for the q3 average, a few tens of float32 units in the last place of a value near 0.6. On stacking, 32 oxygens have an l = 3 sum that cancels: d-SEAMS and OVITO return about 1e-16, and freud's float32 residual is 1.04e-3. The same pattern, at a few times 1e-4, is the whole of the mixed-frame q gap.

CHILL+ in ovito mode matches the OVITO Chill+ modifier on every oxygen of every frame, including the three mixed frames.

The ice cells are entirely hexagonal under d-SEAMS symmetry (label 0), OVITO PTM (`HEX_DIAMOND`, 7), and identify diamond (`HEX_DIAMOND`, 4). On the mixed frames the crystal cores agree and the shells do not:

| frame | PTM | identify diamond |
|---|---|---|
| interface | hex/hex 128, hex/other 2, cubic/other 4, other/other 378 | the 128 hex cores stay hex; 4 cubic cores and 2 hex cores are tagged first neighbours; 122 further oxygens are first or second neighbours of hex |
| stacking | cubic/cubic 128, hex/hex 256, other/other 128 | the same 128 cubic and 256 hex cores; the 128 others are split evenly into cubic and hex first and second neighbours |
| melted | hex/hex 315, hex/other 17, other/hex 2, other/other 178 | 308 of the hex cores stay hex, 24 are tagged first neighbours, and 106 further oxygens are first or second neighbours |

d-SEAMS labels the crystal core (ice Ih, ice Ic, or stacking disorder). PTM at RMSD 0.1 drops a few distorted hex cores that the shell residual still accepts, and names two melted oxygens hex that d-SEAMS leaves unlabelled. Identify diamond then paints the first and second neighbours of those cores. Those tags are not a second crystal class in d-SEAMS.

## Timings

| Molecules | Operation | d-SEAMS 1 | d-SEAMS 8 | freud 1 | freud 8 | OVITO 1 | OVITO 8 |
|---:|---|---:|---:|---:|---:|---:|---:|
| 4096 | cutoff | 0.202 | 0.159 | 2.959 | 1.015 | 1.768 | 1.552 |
| 4096 | knn | 0.556 | 0.434 | 3.552 | 1.156 | 3.581 | 3.956 |
| 4096 | q3 | 0.905 | 0.469 | 10.22 | 2.188 | — | — |
| 4096 | q6 | 0.829 | 0.461 | 14.97 | 2.839 | — | — |
| 4096 | q3 and q6 | 1.714 | 0.980 | — | — | 3.901 | 3.797 |
| 4096 | CHILL+ ovito | 0.637 | 0.323 | — | — | 9.012 | 10.22 |
| 4096 | symmetry / PTM | 0.700 | 0.295 | — | — | 31.34 | 31.79 |
| 4096 | symmetry / diamond | 0.700 | 0.295 | — | — | 2.284 | 2.678 |
| 32768 | cutoff | 2.081 | 1.400 | 26.33 | 7.296 | 3.709 | 3.144 |
| 32768 | knn | 3.234 | 3.453 | 29.69 | 7.656 | 6.169 | 7.230 |
| 32768 | q3 | 7.767 | 3.311 | 81.34 | 18.39 | — | — |
| 32768 | q6 | 7.196 | 3.306 | 119.1 | 20.28 | — | — |
| 32768 | q3 and q6 | 14.91 | 6.581 | — | — | 9.759 | 9.163 |
| 32768 | CHILL+ ovito | 5.679 | 2.034 | — | — | 21.27 | 33.89 |
| 32768 | symmetry / PTM | 5.969 | 1.857 | — | — | 118.0 | 150.8 |
| 32768 | symmetry / diamond | 5.969 | 1.857 | — | — | 10.37 | 9.190 |
| 110592 | cutoff | 7.343 | 4.430 | 106.0 | 28.52 | 8.358 | 9.124 |
| 110592 | knn | 10.88 | 10.67 | 127.9 | 31.60 | 24.90 | 24.51 |
| 110592 | q3 | 28.55 | 11.61 | 296.9 | 57.41 | — | — |
| 110592 | q6 | 26.00 | 11.81 | 443.9 | 74.37 | — | — |
| 110592 | q3 and q6 | 54.68 | 25.52 | — | — | 40.38 | 42.46 |
| 110592 | CHILL+ ovito | 19.88 | 6.277 | — | — | 83.58 | 94.22 |
| 110592 | symmetry / PTM | 22.68 | 6.404 | — | — | 422.7 | 483.5 |
| 110592 | symmetry / diamond | 22.68 | 6.404 | — | — | 35.01 | 35.11 |

OVITO's Steinhardt arm is one modifier with `max_l = 6`, so it computes q2 through q6 and the averages. d-SEAMS q3 and q6 are two calls. Freud's q arm is one `Steinhardt` with `average=True`.

The serial cutoff at 110592 is 7.343 ms against OVITO CreateBonds at 8.358 ms, and 4.430 ms against 9.124 ms at eight threads. The half stencil already tested each pair once and kept both ends in a fixed stash; a second geometry walk was not where the time went. The stash walk now tests four candidates at a time (AVX2), the degree-4 row is sorted with a five-compare network while it is copied into the CSR, and an in-box coordinate is binned with a multiply instead of a division and a floor. `SEAMS_NEIGH_PROFILE` on the steady serial call is about 1.05 ms in the bin, 4.5 ms in that stencil, 0.09 ms in the prefix, and 0.61 ms in the copy. `perf stat` on this guest reports the hardware counters as not supported, so that phase split is the profile. The neighbour pairs are unchanged: the cutoff tests still match the minimum image at one, two, four, and eight threads. On 512-atom frames an eight-thread kNN call can lose a few hundredths of a millisecond to freud; the team setup is the whole of that gap, and at 110592 the same call is 10.67 ms against freud's 31.60 ms.
