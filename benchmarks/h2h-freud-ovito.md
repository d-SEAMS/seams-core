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
| 4096 | cutoff | 0.274 | 0.240 | 2.921 | 0.992 | 1.429 | 1.311 |
| 4096 | knn | 0.502 | 0.507 | 3.483 | 1.048 | 3.366 | 3.390 |
| 4096 | q3 | 0.950 | 0.590 | 9.803 | 2.338 | — | — |
| 4096 | q6 | 0.889 | 0.602 | 14.49 | 2.906 | — | — |
| 4096 | q3 and q6 | 1.827 | 1.181 | — | — | 3.313 | 3.127 |
| 4096 | CHILL+ ovito | 0.713 | 0.424 | — | — | 7.304 | 7.295 |
| 4096 | symmetry / PTM | 0.687 | 0.292 | — | — | 28.70 | 29.17 |
| 4096 | symmetry / diamond | 0.687 | 0.292 | — | — | 2.475 | 2.204 |
| 32768 | cutoff | 2.611 | 1.519 | 24.09 | 8.511 | 3.221 | 3.071 |
| 32768 | knn | 3.436 | 3.432 | 28.75 | 8.558 | 4.989 | 5.144 |
| 32768 | q3 | 8.219 | 4.425 | 79.76 | 15.26 | — | — |
| 32768 | q6 | 7.583 | 4.420 | 114.4 | 20.15 | — | — |
| 32768 | q3 and q6 | 15.75 | 8.720 | — | — | 9.942 | 6.373 |
| 32768 | CHILL+ ovito | 6.198 | 2.810 | — | — | 20.38 | 20.38 |
| 32768 | symmetry / PTM | 5.846 | 2.087 | — | — | 111.4 | 113.1 |
| 32768 | symmetry / diamond | 5.846 | 2.087 | — | — | 8.454 | 8.483 |
| 110592 | cutoff | 8.347 | 5.602 | 91.34 | 28.05 | 7.148 | 8.028 |
| 110592 | knn | 10.52 | 11.13 | 107.0 | 29.14 | 20.36 | 20.78 |
| 110592 | q3 | 27.07 | 15.18 | 270.7 | 51.34 | — | — |
| 110592 | q6 | 25.33 | 13.09 | 399.6 | 67.70 | — | — |
| 110592 | q3 and q6 | 52.14 | 29.42 | — | — | 31.71 | 32.21 |
| 110592 | CHILL+ ovito | 20.49 | 8.704 | — | — | 71.82 | 72.10 |
| 110592 | symmetry / PTM | 21.06 | 7.405 | — | — | 391.3 | 394.3 |
| 110592 | symmetry / diamond | 21.06 | 7.405 | — | — | 31.45 | 30.81 |

OVITO's Steinhardt arm is one modifier with `max_l = 6`, so it computes q2 through q6 and the averages. d-SEAMS q3 and q6 are two calls. Freud's q arm is one `Steinhardt` with `average=True`.

The one place a large cell is slower in d-SEAMS is the serial cutoff at 110592: 8.347 ms against OVITO CreateBonds at 7.148 ms. `SEAMS_NEIGH_PROFILE` puts 5.50 ms of the steady serial call in the pair-count walk (`passFast`), 1.62 ms in the bin, 0.09 ms in the prefix, and 1.09 ms in the fill. That walk is the specialized cell-list kernel, and the neighbour pairs it writes already match both codes exactly, so the count loop was not retuned. At eight threads the same cutoff is 5.602 ms against OVITO's 8.028 ms and freud's 28.05 ms. On 512-atom frames an eight-thread kNN call can lose a few hundredths of a millisecond to freud; the team setup is the whole of that gap, and at 110592 the same call is 11.13 ms against freud's 29.14 ms.
