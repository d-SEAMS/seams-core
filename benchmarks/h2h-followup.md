# Follow-up to the paper review

One session, 8 cores, eight alternating samples. Serial work is pinned to core 0. Eight-thread work uses cores 0-7. A win is a strictly smaller median and a strictly smaller minimum. OVITO is 3.17.2. The raw log is `benchmarks/h2h-followup.txt`. The ring samples are `benchmarks/h2h-recorded-rings.json`.

## Where the recorded ring loss came from

`paper/tex/figures/data/ring-head2head.json` and `paper/tex/figures/data/plot-serial-rings.dat`, git `c5da2d9f`, Intel Xeon, 8 cores, no host name. One thread. Each sample is a fresh `bench_stages` process. The neighbour list is built before the timer. The median is the type-7 median of 12 repeats. OVITO FindRings also runs on one thread.

| Molecules | Frame | d-SEAMS (ms) | OVITO (ms) |
|---:|---|---:|---:|
| 136 | `genice_sH.lammpstrj` | 1.078 | 0.697 |
| 4096 | `mW_cubic.lammpstrj` | 31.617 | 15.831 |
| 110592 | that cubic cell, replicated 3 by 3 by 3 | 864.701 | 511.024 |

The 865 vs 511 pair is the 110592 row, rounded. The same record has OVITO faster at 136 and 4096. The earlier paper-review session timed a different call: `seams_rings_find`, which includes the neighbour build, on ideal Ih and an sII cell, in one process.

This session repeats the recorded boundary (graph outside the timer, those three frames, census checked). d-SEAMS is faster at every size, at one thread and at eight.

| Molecules | Threads | d-SEAMS median / min | OVITO median / min |
|---:|---:|---:|---:|
| 136 | 1 | 0.260 / 0.226 | 0.551 / 0.508 |
| 4096 | 1 | 4.272 / 3.820 | 14.150 / 13.722 |
| 110592 | 1 | 116.044 / 111.428 | 627.160 / 613.367 |
| 136 | 8 | 0.385 / 0.344 | 0.611 / 0.554 |
| 4096 | 8 | 2.620 / 2.093 | 14.633 / 14.279 |
| 110592 | 8 | 33.728 / 30.295 | 334.846 / 315.568 |

The ring census matches on every frame (sH 12/120/68 of sizes 4/5/6, cubic ice 8192 and 221184 hexagons).

## Symmetry, 40.6 ms vs 36.1 ms

That pair is not a labelled median in `benchmarks/` or the paper data. The nearest identify-diamond record is `benchmarks/logs/symmetry-v20-session.txt`, ideal ice, 110592 molecules, one thread: d-SEAMS kernel 21.743 ms, identify-diamond 36.564 ms. On this host the same classifier against identify-diamond on the cubic frames above is a win: 4096 molecules, one thread, 1.750 / 1.664 vs 5.560 / 5.274; 110592 molecules, one thread, 62.650 / 59.515 vs 171.670 / 167.184. Eight threads win as well. No symmetry threshold was changed.

## Carried rings

`nucleation-churn.txt` counts 147 transitions. Every one exceeds the old full-pass rules (changed atoms above one fifth of the frame, or changed bonds above 0.15 per molecule). Those rules used to replace the frame with a search from every vertex, so the carried update cost the same as a full search.

A ring of at most 7 members changes only by meeting the vertices in T, the changed rows and, for a deeper search, the vertices within reach of them. For a maximum ring of 6 the reach is 0, so T is the changed rows. The atom-fraction and bond-rate switches no longer replace that search. A whole-frame enumeration remains only when T is every vertex.

On a 4096-molecule cell with 3670 changed bonds (the old rule fires) the update searches 2881 vertices and matches `ringNetwork`. One thread: carry 20.073 / 19.763 vs full 20.602 / 20.097, a win. Eight threads: 6.590 / 5.670 vs 5.225 / 5.103. The eight-thread partial search is slower than the tuned whole-frame search on this high-churn frame.

A local move of 1200 molecules also trips the old rule (850 vertices recomputed, not 4096) and is faster at both widths: one thread 9.350 / 9.028 vs 18.688 / 17.862; eight threads 3.758 / 3.209 vs 5.420 / 4.621.

## CHILL+ recall, without the circular truth

`benchmarks/paper-cells/symmetry-holdout.csv` is the source of 0.902 vs 0.862 and of 0.975/0.970 vs 0.966/0.944. The higher number is CHILL+ and the lower number is the symmetry gate fitted on the 252 K TIP4P slab. The truth for those `both_recall` columns is PTM diamond and CHILL+ ice together, so CHILL+ is scored against a label that includes itself, and the gate was chosen on the training slab. Those figures are not recomputed that way.

The frames here are built in the session from `ice_positions`. mW uses the 2.76 A hexagonal lattice (512 molecules). TIP4P uses the same lattice scaled to 2.75 A. A second frame of each adds 256 liquid points beyond the lattice, through the periodic boundary, so one face of the crystal meets the liquid. An atom is ice because it was placed on the lattice. Published CHILL+ thresholds only: staggered at `c_ij <= -0.35`, eclipsed on `-0.8 <= c_ij <= 0.25`, four nearest neighbours. The OVITO rule uses the open edges and every neighbour inside 3.5 A. Ice in the recall is cubic or hexagonal.

| Frame | Paper recall | OVITO-rule recall | Modifier recall | False ice in the liquid |
|---|---:|---:|---:|---:|
| mW lattice | 1 | 1 | 1 | 0 |
| mW lattice plus liquid | 0.775 (397/512) | 0.785 (402/512) | 0.785 | 0 |
| TIP4P lattice | 1 | 1 | 1 | 0 |
| TIP4P lattice plus liquid | 0.773 (396/512) | 0.783 (401/512) | 0.783 | 0 |

On the pure lattices every atom is hexagonal under the paper rule, the OVITO rule, and OVITO 3.17.2. On the seeded frames the paper rule and the OVITO rule differ by the five or six interfacial atoms the open neighbour rule calls hexagonal. The liquid contributes no cubic or hexagonal label. Thresholds were not changed.

The equilibrated holdout trajectories that produced the old table are built by `scripts/paper_symmetry_holdout.py` when `lmp` is on `PATH`. They are not committed. `input/mw215/fetch.py` already fetches the public set `d-SEAMS/mw-liquid-215K` when a run needs those frames.

## SQQ and the triclinic sH cell

SQQ 0.5.8 accepts `genice_sI.lammpstrj` and rejects `genice_sH.lammpstrj` with `Triclinic/tilted LAMMPS dump boxes are not supported`. Setting `pbc.box_mode` to `triclinic` raises `pbc.box_mode must be orthorhombic`. The citations are in `benchmarks/h2h-sqq-box.md`. `tests/test_sqq_sh.py` runs the same check.
