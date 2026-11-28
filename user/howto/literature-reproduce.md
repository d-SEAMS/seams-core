<!-- https://docs.dseams.info/user/howto/literature-reproduce.md -->

# Reproduce the literature head-to-head {#reproduce-the-literature-head-to-head}

`scripts/h2h_rc3_sitting.py` writes the six-call sitting. One warm-up,
then six timed calls. The median of that even sample is the mean of
the two central values, and the interquartile range is Hyndman-Fan
type 7. The medians are the long table on
[Literature head-to-head](literature-h2h.md). This page is how
each guest code receives a thread width, and what that width did to
the median.

The driver starts every row under `taskset`. One thread is
`taskset -c 0`. Eight threads are `taskset -c 0-7`. It sets
`OMP_NUM_THREADS` to that width and `OMP_DYNAMIC=FALSE`. d-SEAMS
follows that OpenMP width. Each guest then receives the same width
through its own control, except PLUMED and ClavisCage.

| code | what is timed | thread control |
| --- | --- | --- |
| freud | local q_4 and q_6 | freud.parallel.set_num_threads(n) |
| OVITO | CHILL+, primitive rings, polyhedral template matching | OVITO_THREAD_COUNT=n before import ovito |
| mdapy | CHILL+, q_4, q_6, rings | mdapy.set_num_threads(n) |
| vesin | full cell-list neighbour pairs, copy included | NeighborList(..., n_threads=n, algorithm”cell_list”)= |
| PLUMED | symfunc q_4 and q_6 | OMP_NUM_THREADS is set. PLUMED_NUM_THREADS is not |
| rings-statistics-matter | primitive rings | ring_statistics(..., num_threads=n, backend”rust”)= |
| fastatomstruct | King’s shortest paths through size 6 | fastatomstruct.set_num_threads(n) |
| ClavisCage | TRACE cages | the analyze_waters call takes no thread count |

Polyhedral template matching is the OVITO modifier. Structure ids 6
and 7 are enabled and the other templates are left off. The d-SEAMS
cell on that row is `seams_frame_symmetry` with tolerance 0.12. Both
are labeled PTM. rings-statistics-matter is imported as
`rings_statistics_matter`.

The worker opens `bbdir-release/bindings/libseams_lang.so`. `PLUMED_LIB`
names `libplumed.so`; otherwise the script takes that library from the
`plumed` executable on `PATH`.

```bash
python3 scripts/h2h_rc3_sitting.py
```

The vesin probe runs first.

## What the medians do {#what-the-medians-do}

On ideal ice Ic at 110592 and 884736 molecules, freud, OVITO, mdapy,
vesin, and rings-statistics-matter are faster at eight threads than
at one. OVITO’s template matching and its ring finder fall with the
same `OVITO_THREAD_COUNT` as its CHILL+. `t1` and `t8` are the
medians in the threads=1 and threads=8 columns. `ratio` is `t1 / t8`.

PLUMED is the exception in this file. The driver sets
`OMP_NUM_THREADS` and leaves `PLUMED_NUM_THREADS` unset. In PLUMED
2.10 `PLMD::OpenMP::getNumThreads` reads `PLUMED_NUM_THREADS`, and an
unset variable is one thread. The `threads=1` column and the
`threads=8` column are that same measurement, so the ratio stays
near 1. The sitting that sets `PLUMED_NUM_THREADS` to the row width
is on [the performance overview](../../explanation/performance-overview.md).
ClavisCage also has no thread argument. On the sI unit cell the
medians are 51.122 ms and 51.206 ms. On the 368-oxygen replicate they
are 301.514 ms and 302.471 ms. On sII they are 69.332 ms and 69.657 ms.
The first number of each pair is one thread and the second is eight.

| code | operation | molecules | t1 ms | t8 ms | ratio | note |
| --- | --- | --- | --- | --- | --- | --- |
| freud | q4 | 110592 | 276.580 | 51.302 | 5.39 |  |
| freud | q4 | 884736 | 2358.091 | 431.376 | 5.47 |  |
| freud | q6 | 110592 | 362.326 | 68.760 | 5.27 |  |
| freud | q6 | 884736 | 2973.372 | 550.161 | 5.40 |  |
| OVITO | CHILL+ | 110592 | 622.424 | 88.685 | 7.02 |  |
| OVITO | CHILL+ | 884736 | 4890.789 | 666.288 | 7.34 |  |
| OVITO | PTM | 110592 | 2790.624 | 370.341 | 7.54 |  |
| OVITO | PTM | 884736 | 21304.364 | 2849.382 | 7.48 |  |
| OVITO | rings | 110592 | 583.334 | 354.605 | 1.65 |  |
| OVITO | rings | 884736 | 4883.899 | 2582.748 | 1.89 |  |
| mdapy | CHILL+ | 110592 | 61.938 | 12.379 | 5.00 |  |
| mdapy | CHILL+ | 884736 | 744.487 | 125.578 | 5.93 |  |
| mdapy | q4 | 110592 | 64.687 | 12.400 | 5.22 |  |
| mdapy | q4 | 884736 | 807.101 | 96.238 | 8.39 |  |
| mdapy | q6 | 110592 | 79.232 | 14.572 | 5.44 |  |
| mdapy | q6 | 884736 | 921.574 | 107.347 | 8.59 |  |
| mdapy | rings | 110592 | 2276.747 | 375.967 | 6.06 |  |
| mdapy | rings | 884736 | 19221.718 | 3361.913 | 5.72 |  |
| vesin | neighbours | 110592 | 51.194 | 22.269 | 2.30 |  |
| vesin | neighbours | 884736 | 1112.702 | 196.866 | 5.65 |  |
| rings-statistics-matter | rings | 110592 | 1151.862 | 264.953 | 4.35 |  |
| rings-statistics-matter | rings | 884736 | 11504.451 | 2665.750 | 4.32 |  |
| PLUMED | q4 | 110592 | 1310.025 | 1316.423 | 1.00 | one thread; PLUMED_NUM_THREADS unset |
| PLUMED | q4 | 884736 | 10053.680 | 10122.158 | 0.99 | one thread; PLUMED_NUM_THREADS unset |
| PLUMED | q6 | 110592 | 1685.454 | 1680.281 | 1.00 | one thread; PLUMED_NUM_THREADS unset |
| PLUMED | q6 | 884736 | 12947.767 | 12818.837 | 1.01 | one thread; PLUMED_NUM_THREADS unset |
