# Paper-review head-to-head

One session on 8 cores. Eight alternating samples. Serial calls are pinned to core 0. Eight-thread calls use cores 0-7, with `OMP_PROC_BIND=close`, `OMP_PLACES=cores`, and `OMP_WAIT_POLICY=PASSIVE`. A win is a strictly smaller median and a strictly smaller minimum. OVITO is 3.17.2. The raw log is `benchmarks/h2h-paper-review.txt`.

The recorded losses below are the numbers from the review. They are not from this session. This session did not reproduce them.

## Before and after, milliseconds

Median, then minimum. The ring and symmetry rows are this host against OVITO. The cage rows are the library signature against the `seams` command.

| Comparison | Recorded | This session, 1 thread | This session, 8 threads |
|---|---|---|---|
| Rings, 136 molecules | 865 vs 511 | 0.312 / 0.301 vs 0.796 / 0.750 | 0.307 / 0.290 vs 0.902 / 0.816 |
| Rings, 4096 molecules | 865 vs 511 | 4.511 / 4.417 vs 16.700 / 15.817 | 2.466 / 2.264 vs 17.310 / 16.401 |
| Rings, 110592 molecules | 865 vs 511 | 129.317 / 127.326 vs 612.537 / 593.639 | 47.643 / 43.895 vs 638.421 / 604.509 |
| Symmetry, 4096, identify-diamond | 40.6 vs 36.1 | 0.770 / 0.740 vs 5.016 / 4.790 | 0.986 / 0.691 vs 5.096 / 4.884 |
| Symmetry, 110592, identify-diamond | 40.6 vs 36.1 | 24.704 / 23.146 vs 144.620 / 143.422 | 11.472 / 9.717 vs 152.911 / 147.840 |
| Cages, sII 69632, signature 51264 | library 282 vs command 213 | 267.195 / 260.982 vs 339.867 / 330.757 | 90.277 / 86.860 vs 170.621 / 168.468 |

Every row on this host is a win. Ideal ice is hexagonal on every atom for both symmetry codes (4096 of 4096, 110592 of 110592). The cage census is 4096 cages at both thread counts.

## Rings and symmetry

The ring path and the symmetry classifier are unchanged. On this host the existing Franzblau census and `seams_system_symmetry` (shell RMSD ceiling 0.15, the caller passes 0) are already faster than OVITO FindRings and IdentifyDiamond at 136, 4096, and 110592 molecules, at one thread and at eight. The recorded 865 ms against 511 ms, and 40.6 ms against 36.1 ms, did not recur.

## Library cages

`seams_cages_signature` copied the whole frame, including the atom-id map, and rebuilt the cutoff graph through atom ids. A typed call now keeps the caller's cloud, and an orthorhombic single-type frame builds the graph with the cell-list CSR. One thread searches the packed rows. A team searches the per-atom vectors, which is the layout the command uses. A tilted cell still uses the neighbour-list builder, because the tilted cell list is not the same pair set.

The sII 8x8x8 frame (69632 molecules, orthorhombic) returns 4096 cages of signature `51264` from the library and from `seams cages`. The sI unit cell returns 2 cages of `512`. The tilted sH cell, with the LAMMPS edges and the xy tilt, returns 12 cages of `512` and 4 of `51268`, the same counts as the command.

## CHILL+

Thresholds are unchanged. The published rule is Nguyen and Molinero, JPCB 119, 9369 (2015): the four nearest neighbours, staggered when `c_ij <= -0.35`, eclipsed when `-0.8 <= c_ij <= 0.25`. The OVITO rule uses every neighbour inside 3.5 A and the open edges (`c_ij < -0.35`, `-0.8 < c_ij < 0.25`). That option is `SEAMS_CHILL_OVITO` / `rule="ovito"`, already on this branch. `tests/test_ovito_chill_atoms.py` compares it with OVITO 3.17.2 `ChillPlusModifier` on every atom.

The mW recall pair 0.902 against 0.862, and the TIP4P pairs 0.975/0.970 against 0.966/0.944, are not recomputed here. Those labeled frames are not in the tree. On the frames that are, the published rule and the OVITO rule agree on ideal Ih, ideal Ic, and the stacking slab. They differ on the interface and the melted slab, which is the four-nearest rule against the 3.5 A neighbour rule. The OVITO rule matches the modifier on every atom of all five frames (0 mismatches).

Rows are `paper label / ovito-rule label = count`, then `ovito-rule / modifier`.

| Frame | Paper against OVITO rule | OVITO rule against the modifier |
|---|---|---|
| Ih 512 | hexagonal/hexagonal=512 | hexagonal/hexagonal=512 |
| Ic 512 | cubic/cubic=512 | cubic/cubic=512 |
| stacking 512 | cubic/cubic=128 hexagonal/hexagonal=320 water/water=64 | cubic/cubic=128 hexagonal/hexagonal=320 water/water=64 |
| interface 512 | clathrate/clathrate=1 clathrate/water=1 cubic/cubic=1 hexagonal/hexagonal=138 hexagonal/interfacial=3 interclathrate/clathrate=4 interclathrate/interclathrate=6 interclathrate/water=4 interfacial/hexagonal=1 interfacial/interfacial=44 water/clathrate=1 water/interclathrate=6 water/water=302 | clathrate/clathrate=6 cubic/cubic=1 hexagonal/hexagonal=139 interclathrate/interclathrate=12 interfacial/interfacial=47 water/water=307 |
| melted 512 | clathrate/interclathrate=1 hexagonal/hexagonal=364 hexagonal/interfacial=5 interclathrate/water=4 interfacial/hexagonal=1 interfacial/interfacial=19 interfacial/water=4 water/interclathrate=2 water/water=112 | hexagonal/hexagonal=365 interclathrate/interclathrate=3 interfacial/interfacial=24 water/water=120 |
