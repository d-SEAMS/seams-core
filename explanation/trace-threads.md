<!-- https://docs.dseams.info/explanation/trace-threads.md -->

# Signature-search thread controls {#signature-search-thread-controls}

The signature search reads one OpenMP width, and the neighbour list may
narrow it. The tables record that width on hydrate cells and on ice Ic.
`threads` is the width that was requested. On the signature and kernel
tables, `cores` is the number of CPUs the process was allowed. `n` is
the number of timed calls that were kept.

## How the rows were timed {#how-the-rows-were-timed}

The signature and kernel rows at width $n$ were run with this environment:

| variable | value |
| --- | --- |
| OMP_NUM_THREADS | n |
| OMP_DYNAMIC | FALSE |
| OMP_PROC_BIND | close |
| OMP_PLACES | cores |
| RAYON_NUM_THREADS | n |
| OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, NUMEXPR_NUM_THREADS | 1 |

The process was pinned to cpu 0 when $n = 1$, and to cpus 0 through
$n - 1$ otherwise. One warm-up was dropped. Eight timed calls remain.
The median of an even count is the mean of the two central samples. The
interquartile range is the type-7 upper quartile minus the lower one.

## The width the search reads {#the-width-the-search-reads}

`seams_set_num_threads` sets the OpenMP width later calls read, and the
linkcell width with it. Zero restores the width saved at first use.
`seams_get_num_threads` reports the width the next call will request.

`seams::teamFor` returns that OpenMP width. The atom count is not read.
`seams::TeamWidth` applies the width for the rest of the call and
restores the previous width afterwards. Bond order, published CHILL+,
the k-nearest list, the shell label and local $q_l$ use `TeamWidth`.
Those calls do not start a Rayon pool.

`seams::adoptThreads` calls `prepareCompute`, then `lc_set_threads` with
the OpenMP width. `prepareCompute` calls `omp_set_dynamic(0)` and leaves
`OMP_WAIT_POLICY` as the process started it. The environment above does
not set that variable. The linkcell workers sleep between searches.
Forcing a passive OpenMP wait on top of that sleep made the next region
late by tens of milliseconds, so the OpenMP team is left running.

`lc_set_threads` stores the width of later linkcell searches. A count
below one is stored as one thread, and one thread is a serial walk. A
width keeps its pool, and a later search at the same width reuses it.
A different width uses its own pool.

The C calls `seams_signature_cages` and `seams_signature_cages_edge`
enter through `seams::FairWidth`. When `SEAMS_FAIR_WIDTH` is unset, two
forks in a row that each take more than a millisecond halve the team,
down to a quarter of the requested width. A narrowed team tries twice
its width again at most once a second. `SEAMS_FAIR_WIDTH=0` keeps the
requested width. `seams_frame_signature` and `seams cages --signature`
do not enter through `FairWidth`.

[C](../stable/c.md) documents `seams_set_num_threads` and `seams_get_num_threads`.

## The neighbour list {#the-neighbour-list}

The timed hydrate sample is one cutoff list from `neighListLinkcell`,
one depth-6 primitive-ring network, and one search of $5^{12}$,
$5^{12}6^{2}$ and $5^{12}6^{4}$ together. Below 4096 molecules that
list caps the team at four and sets the OpenMP width to the cap. The
cap is left in place, so the ring search that follows runs at four
threads as well. A request for eight threads on sI $2\times2\times2$
(368 molecules), sI $4\times4\times4$ (2944) and sII $2\times2\times2$
(1088) therefore has the same width as the four-thread row, and the
medians sit together. At 4096 molecules and above the requested width
is kept.

Several threads in that list build an orthorhombic or triclinic cutoff
with the cell-list CSR. One thread uses linkcell’s buffered walk. If the
cell list cannot be built, the call uses that walk as well. Both builders
write the same sorted cutoff set, so the ring search does not depend on
which builder ran.

`seams cages --signature` on a cutoff graph calls `neighListO`, which
does not apply the 4096 cap. The C API builds an orthorhombic cutoff
with the cell list when every atom has the requested type, and uses
`neighListO` for a triclinic frame or a frame that mixes types.

`seams_frame_signature` on a single census uses the largest face as the
depth, so `512` alone is a depth-5 network. The timed sample above
searches all three censuses on one depth-6 network.

When the pairs come from a linkcell pair list, the rows stay on the
calling thread until the frame has at least 4096 atoms and the pair
list has at least $2^{20}$ pairs. Waking a team that slept through
the search costs more than writing those rows on the caller until the
list is that large. The several-thread cell-list path does not use
that placement.

| Signature cages on the sI and sII cells. threads is the width that was requested. n is the number of timed calls. The median of an even count is the mean of the two central samples. iqr_ms is the type-7 interquartile range. |
| --- |

### signature cages, sI {#signature-cages-si}

| cell | molecules | threads | cores | n | median_ms | iqr_ms |
| --- | --- | --- | --- | --- | --- | --- |
| sI 2x2x2 | 368 | 1 | 8 | 8 | 1.671145 | 0.042545 |
| sI 2x2x2 | 368 | 2 | 8 | 8 | 1.194630 | 0.144788 |
| sI 2x2x2 | 368 | 4 | 8 | 8 | 0.845288 | 0.093749 |
| sI 2x2x2 | 368 | 8 | 8 | 8 | 0.889254 | 0.067923 |
| sI 4x4x4 | 2944 | 1 | 8 | 8 | 11.155500 | 0.052250 |
| sI 4x4x4 | 2944 | 2 | 8 | 8 | 7.572065 | 0.876455 |
| sI 4x4x4 | 2944 | 4 | 8 | 8 | 5.198810 | 0.470322 |
| sI 4x4x4 | 2944 | 8 | 8 | 8 | 5.295710 | 0.664862 |
| sI 6x6x6 | 9936 | 1 | 8 | 8 | 34.513400 | 0.145950 |
| sI 6x6x6 | 9936 | 2 | 8 | 8 | 29.893900 | 0.793650 |
| sI 6x6x6 | 9936 | 4 | 8 | 8 | 15.499400 | 0.382550 |
| sI 6x6x6 | 9936 | 8 | 8 | 8 | 11.406800 | 1.692000 |
| sI 8x8x8 | 23552 | 1 | 8 | 8 | 79.230800 | 0.313300 |
| sI 8x8x8 | 23552 | 2 | 8 | 8 | 50.058800 | 0.456100 |
| sI 8x8x8 | 23552 | 4 | 8 | 8 | 33.261250 | 4.790825 |
| sI 8x8x8 | 23552 | 8 | 8 | 8 | 26.431450 | 0.435175 |

### signature cages, sII {#signature-cages-sii}

| cell | molecules | threads | cores | n | median_ms | iqr_ms |
| --- | --- | --- | --- | --- | --- | --- |
| sII 2x2x2 | 1088 | 1 | 8 | 8 | 4.815100 | 0.054395 |
| sII 2x2x2 | 1088 | 2 | 8 | 8 | 3.748070 | 0.253115 |
| sII 2x2x2 | 1088 | 4 | 8 | 8 | 2.464180 | 0.202713 |
| sII 2x2x2 | 1088 | 8 | 8 | 8 | 2.543545 | 0.405943 |
| sII 4x4x4 | 8704 | 1 | 8 | 8 | 35.907600 | 0.137650 |
| sII 4x4x4 | 8704 | 2 | 8 | 8 | 22.485750 | 0.262925 |
| sII 4x4x4 | 8704 | 4 | 8 | 8 | 14.209000 | 0.490550 |
| sII 4x4x4 | 8704 | 8 | 8 | 8 | 11.173800 | 3.512300 |
| sII 6x6x6 | 29376 | 1 | 8 | 8 | 121.799000 | 0.416000 |
| sII 6x6x6 | 29376 | 2 | 8 | 8 | 74.774150 | 1.357250 |
| sII 6x6x6 | 29376 | 4 | 8 | 8 | 45.948400 | 0.754300 |
| sII 6x6x6 | 29376 | 8 | 8 | 8 | 31.094250 | 0.210175 |
| sII 8x8x8 | 69632 | 1 | 8 | 8 | 292.810500 | 2.154000 |
| sII 8x8x8 | 69632 | 2 | 8 | 8 | 180.097000 | 1.709250 |
| sII 8x8x8 | 69632 | 4 | 8 | 8 | 111.554000 | 3.843000 |
| sII 8x8x8 | 69632 | 8 | 8 | 8 | 73.108700 | 0.648400 |

## The ring update {#the-ring-update}

`sourceTeam` sizes the enumeration of changed ring sources, and the
carried update’s full pass. The calling thread does the work until
there are about 80 sources for every worker. `teamFor` is a different
rule: it is the requested width at every frame size.

On one thread the carried update copies the new rings and sorts them
directly, then merges the rings that stay with a direct walk. Several
threads sort indices, scan for rings that meet the changed set, and
merge the rings that are kept, using `teamFor`. Dropped ring buffers
are kept and reused by the next frame.

The cage index is merged from per-thread buckets in ring order.
Membership is the sorted atom ids of each face. The census at every
width is the one-thread census.

## Eight kernels {#eight-kernels}

Eight kernels are timed on ideal ice Ic at 4096, 32768, 110592 and
884736 molecules, each at 1, 2, 4 and 8 threads. They are primitive
rings, local and cutoff $q_4$ and $q_6$, CHILL+, the shell label,
and one carried hop update. The neighbour graph is built before the
timer starts. The rings kernel is `primitive::ringNetwork` on that
graph. Local $q_4$ and $q_6$ read the same graph. Cutoff $q_4$,
cutoff $q_6$ and published CHILL+ build their own cutoff list inside
the timed call. The hop kernel times the second update, after the first
update has filled the carried state. A repeat whose fingerprint differs
from the first call aborts the run. On 4096 molecules the kernels match
the one-thread result.

`seams_frame_rings` builds the cutoff list and then searches, so its
time is not this rings kernel.

| Eight kernels on ideal ice Ic, each at 1, 2, 4 and 8 threads. n is the number of timed calls. The median of an even count is the mean of the two central samples. iqr_ms is the type-7 interquartile range. |
| --- |

### ice Ic, 4096 molecules {#ice-ic-4096-molecules}

| kernel | threads | cores | n | median_ms | iqr_ms |
| --- | --- | --- | --- | --- | --- |
| rings | 1 | 8 | 8 | 2.466905 | 0.075503 |
| rings | 2 | 8 | 8 | 1.495700 | 0.050533 |
| rings | 4 | 8 | 8 | 0.940597 | 0.151879 |
| rings | 8 | 8 | 8 | 0.633803 | 0.038622 |
| q4-local | 1 | 8 | 8 | 0.667118 | 0.018324 |
| q4-local | 2 | 8 | 8 | 0.389290 | 0.007018 |
| q4-local | 4 | 8 | 8 | 0.267276 | 0.001623 |
| q4-local | 8 | 8 | 8 | 0.262831 | 0.002778 |
| q6-local | 1 | 8 | 8 | 0.771595 | 0.003968 |
| q6-local | 2 | 8 | 8 | 0.451925 | 0.007984 |
| q6-local | 4 | 8 | 8 | 0.380479 | 0.008200 |
| q6-local | 8 | 8 | 8 | 0.391622 | 0.020511 |
| q4-cutoff | 1 | 8 | 8 | 0.612442 | 0.006985 |
| q4-cutoff | 2 | 8 | 8 | 0.378497 | 0.008923 |
| q4-cutoff | 4 | 8 | 8 | 0.242802 | 0.039543 |
| q4-cutoff | 8 | 8 | 8 | 0.243986 | 0.037453 |
| q6-cutoff | 1 | 8 | 8 | 0.637176 | 0.018393 |
| q6-cutoff | 2 | 8 | 8 | 0.412152 | 0.013879 |
| q6-cutoff | 4 | 8 | 8 | 0.230091 | 0.008668 |
| q6-cutoff | 8 | 8 | 8 | 0.242103 | 0.006373 |
| chill | 1 | 8 | 8 | 0.907933 | 0.045465 |
| chill | 2 | 8 | 8 | 0.552909 | 0.012054 |
| chill | 4 | 8 | 8 | 0.348419 | 0.008053 |
| chill | 8 | 8 | 8 | 0.349171 | 0.004840 |
| symmetry | 1 | 8 | 8 | 0.913122 | 0.011245 |
| symmetry | 2 | 8 | 8 | 0.691176 | 0.007351 |
| symmetry | 4 | 8 | 8 | 0.414814 | 0.012065 |
| symmetry | 8 | 8 | 8 | 0.414618 | 0.021285 |
| hop | 1 | 8 | 8 | 0.978623 | 0.022944 |
| hop | 2 | 8 | 8 | 0.597030 | 0.034203 |
| hop | 4 | 8 | 8 | 0.410767 | 0.017664 |
| hop | 8 | 8 | 8 | 0.381499 | 0.009265 |

### ice Ic, 32768 molecules {#ice-ic-32768-molecules}

| kernel | threads | cores | n | median_ms | iqr_ms |
| --- | --- | --- | --- | --- | --- |
| rings | 1 | 8 | 8 | 19.360550 | 0.379825 |
| rings | 2 | 8 | 8 | 10.805550 | 0.215250 |
| rings | 4 | 8 | 8 | 6.449795 | 0.082868 |
| rings | 8 | 8 | 8 | 4.263060 | 0.236065 |
| q4-local | 1 | 8 | 8 | 5.905020 | 0.029900 |
| q4-local | 2 | 8 | 8 | 3.253615 | 0.053870 |
| q4-local | 4 | 8 | 8 | 1.909015 | 0.111800 |
| q4-local | 8 | 8 | 8 | 1.299955 | 0.023633 |
| q6-local | 1 | 8 | 8 | 6.897740 | 0.058900 |
| q6-local | 2 | 8 | 8 | 3.728225 | 0.037745 |
| q6-local | 4 | 8 | 8 | 3.084135 | 0.023317 |
| q6-local | 8 | 8 | 8 | 1.956415 | 0.062893 |
| q4-cutoff | 1 | 8 | 8 | 5.263250 | 0.058237 |
| q4-cutoff | 2 | 8 | 8 | 2.753380 | 0.027987 |
| q4-cutoff | 4 | 8 | 8 | 1.770700 | 0.079300 |
| q4-cutoff | 8 | 8 | 8 | 1.213175 | 0.046030 |
| q6-cutoff | 1 | 8 | 8 | 5.426560 | 0.037195 |
| q6-cutoff | 2 | 8 | 8 | 2.854520 | 0.034102 |
| q6-cutoff | 4 | 8 | 8 | 1.822180 | 0.086745 |
| q6-cutoff | 8 | 8 | 8 | 1.316010 | 0.059367 |
| chill | 1 | 8 | 8 | 7.429445 | 0.132305 |
| chill | 2 | 8 | 8 | 3.889890 | 0.125662 |
| chill | 4 | 8 | 8 | 2.437750 | 0.107562 |
| chill | 8 | 8 | 8 | 1.680220 | 0.097135 |
| symmetry | 1 | 8 | 8 | 7.109115 | 0.078567 |
| symmetry | 2 | 8 | 8 | 4.891060 | 0.148090 |
| symmetry | 4 | 8 | 8 | 2.897955 | 0.084083 |
| symmetry | 8 | 8 | 8 | 1.865990 | 0.062373 |
| hop | 1 | 8 | 8 | 7.778165 | 0.065072 |
| hop | 2 | 8 | 8 | 4.798770 | 0.679093 |
| hop | 4 | 8 | 8 | 2.998450 | 0.069390 |
| hop | 8 | 8 | 8 | 2.032610 | 0.068672 |

### ice Ic, 110592 molecules {#ice-ic-110592-molecules}

| kernel | threads | cores | n | median_ms | iqr_ms |
| --- | --- | --- | --- | --- | --- |
| rings | 1 | 8 | 8 | 62.012050 | 0.158025 |
| rings | 2 | 8 | 8 | 35.865750 | 0.459375 |
| rings | 4 | 8 | 8 | 21.042200 | 0.582800 |
| rings | 8 | 8 | 8 | 13.537600 | 0.173700 |
| q4-local | 1 | 8 | 8 | 20.345300 | 0.172975 |
| q4-local | 2 | 8 | 8 | 11.096000 | 0.052775 |
| q4-local | 4 | 8 | 8 | 6.376390 | 0.068050 |
| q4-local | 8 | 8 | 8 | 4.051045 | 0.039383 |
| q6-local | 1 | 8 | 8 | 23.794950 | 0.063875 |
| q6-local | 2 | 8 | 8 | 12.896700 | 0.082325 |
| q6-local | 4 | 8 | 8 | 7.310575 | 0.056100 |
| q6-local | 8 | 8 | 8 | 4.613615 | 0.274352 |
| q4-cutoff | 1 | 8 | 8 | 16.288450 | 0.182575 |
| q4-cutoff | 2 | 8 | 8 | 8.978580 | 0.105928 |
| q4-cutoff | 4 | 8 | 8 | 5.604980 | 0.102873 |
| q4-cutoff | 8 | 8 | 8 | 3.591890 | 0.234900 |
| q6-cutoff | 1 | 8 | 8 | 16.811850 | 0.136250 |
| q6-cutoff | 2 | 8 | 8 | 9.404800 | 0.225812 |
| q6-cutoff | 4 | 8 | 8 | 5.770415 | 0.166370 |
| q6-cutoff | 8 | 8 | 8 | 3.852335 | 0.112190 |
| chill | 1 | 8 | 8 | 24.159000 | 0.204625 |
| chill | 2 | 8 | 8 | 13.322850 | 0.096675 |
| chill | 4 | 8 | 8 | 8.032445 | 0.149675 |
| chill | 8 | 8 | 8 | 5.275310 | 0.284212 |
| symmetry | 1 | 8 | 8 | 24.393350 | 0.072975 |
| symmetry | 2 | 8 | 8 | 17.057900 | 0.214850 |
| symmetry | 4 | 8 | 8 | 13.666650 | 0.453375 |
| symmetry | 8 | 8 | 8 | 6.160520 | 0.194820 |
| hop | 1 | 8 | 8 | 26.381800 | 0.228900 |
| hop | 2 | 8 | 8 | 16.760500 | 2.426175 |
| hop | 4 | 8 | 8 | 9.924360 | 0.561755 |
| hop | 8 | 8 | 8 | 7.194500 | 0.529503 |

### ice Ic, 884736 molecules {#ice-ic-884736-molecules}

| kernel | threads | cores | n | median_ms | iqr_ms |
| --- | --- | --- | --- | --- | --- |
| rings | 1 | 8 | 8 | 495.117000 | 2.365750 |
| rings | 2 | 8 | 8 | 299.999500 | 2.423750 |
| rings | 4 | 8 | 8 | 184.109500 | 3.346250 |
| rings | 8 | 8 | 8 | 112.990000 | 2.134500 |
| q4-local | 1 | 8 | 8 | 201.293000 | 6.906000 |
| q4-local | 2 | 8 | 8 | 108.775500 | 2.981250 |
| q4-local | 4 | 8 | 8 | 61.171350 | 2.056050 |
| q4-local | 8 | 8 | 8 | 40.111950 | 2.153775 |
| q6-local | 1 | 8 | 8 | 239.785500 | 4.282750 |
| q6-local | 2 | 8 | 8 | 129.450500 | 4.643500 |
| q6-local | 4 | 8 | 8 | 80.323450 | 6.257450 |
| q6-local | 8 | 8 | 8 | 51.277600 | 0.623550 |
| q4-cutoff | 1 | 8 | 8 | 154.053000 | 7.017500 |
| q4-cutoff | 2 | 8 | 8 | 80.075400 | 1.599300 |
| q4-cutoff | 4 | 8 | 8 | 43.635150 | 1.151825 |
| q4-cutoff | 8 | 8 | 8 | 25.552200 | 0.728500 |
| q6-cutoff | 1 | 8 | 8 | 158.820000 | 2.469750 |
| q6-cutoff | 2 | 8 | 8 | 85.462600 | 0.926100 |
| q6-cutoff | 4 | 8 | 8 | 45.467300 | 0.915550 |
| q6-cutoff | 8 | 8 | 8 | 26.620450 | 0.664350 |
| chill | 1 | 8 | 8 | 232.527500 | 2.841500 |
| chill | 2 | 8 | 8 | 122.900000 | 0.568750 |
| chill | 4 | 8 | 8 | 67.026200 | 0.327775 |
| chill | 8 | 8 | 8 | 39.195350 | 0.346225 |
| symmetry | 1 | 8 | 8 | 214.114000 | 11.427750 |
| symmetry | 2 | 8 | 8 | 138.153500 | 7.535750 |
| symmetry | 4 | 8 | 8 | 79.281800 | 2.953025 |
| symmetry | 8 | 8 | 8 | 48.350600 | 1.473875 |
| hop | 1 | 8 | 8 | 256.551500 | 6.978250 |
| hop | 2 | 8 | 8 | 171.991500 | 8.024750 |
| hop | 4 | 8 | 8 | 103.765500 | 4.727250 |
| hop | 8 | 8 | 8 | 62.761750 | 2.417975 |

## The full cutoff ring call {#the-full-cutoff-ring-call}

The full call is the neighbour list, the index conversion and the
depth-6 search, timed together on ice Ic at 884736 molecules.
`cell-list` is the orthorhombic list inside that call. `linkcell` is
the buffered walk. `link`, `convert` and `search` are the three parts.
`total` is the whole call. This is not the rings kernel above.

Several threads take the half stencil one z layer at a time, in one
team. A layer stores a pair only on its own atoms. The partner row is
filled after the team has stopped, because a forward step in x reaches
back one z layer. An interior cell has no periodic image. Its thirteen
forward neighbours are plain index steps and the shift is zero, so
that cell skips the periodic-shift arithmetic. A boundary cell still
applies the shift. One thread uses the same half stencil without a
team. Below 2048 atoms the orthorhombic list runs on one thread and
then restores the caller’s width. Published CHILL+ and cutoff $q_6$
read this list. The one-thread and several-thread pair lists match.

The 4096-molecule cap and the $2^{20}$-pair row placement belong to
the signature search above. `teamFor` does not replace them.

| Ice Ic at 884736 molecules. cell-list is the orthorhombic cutoff list. linkcell is the buffered walk. total is the neighbour list, the index conversion and the depth-6 search. link, convert and search are those three parts. n is the number of timed calls. The median of an even count is the mean of the two central samples. iqr_ms is the type-7 interquartile range, filled on the total rows. |
| --- |

### Full cutoff ring call {#full-cutoff-ring-call}

| code | operation | threads | n | median_ms | iqr_ms |
| --- | --- | --- | --- | --- | --- |
| cell-list | total | 1 | 8 | 654.896 | 8.296 |
| cell-list | total | 2 | 8 | 418.741 | 3.510 |
| cell-list | total | 4 | 8 | 252.877 | 1.991 |
| cell-list | total | 8 | 8 | 173.215 | 2.554 |
| linkcell | total | 1 | 8 | 653.639 | 9.726 |
| linkcell | total | 2 | 8 | 521.980 | 4.294 |
| linkcell | total | 4 | 8 | 323.101 | 1.161 |
| linkcell | total | 8 | 8 | 228.632 | 5.324 |
| cell-list | link | 1 | 8 | 132.279 |  |
| cell-list | link | 2 | 8 | 86.208 |  |
| cell-list | link | 4 | 8 | 52.179 |  |
| cell-list | link | 8 | 8 | 35.622 |  |
| linkcell | link | 1 | 8 | 131.875 |  |
| linkcell | link | 2 | 8 | 206.762 |  |
| linkcell | link | 4 | 8 | 137.551 |  |
| linkcell | link | 8 | 8 | 107.620 |  |
| cell-list | convert | 1 | 8 | 28.719 |  |
| cell-list | convert | 2 | 8 | 15.368 |  |
| cell-list | convert | 4 | 8 | 7.783 |  |
| cell-list | convert | 8 | 8 | 4.628 |  |
| linkcell | convert | 1 | 8 | 28.935 |  |
| linkcell | convert | 2 | 8 | 15.060 |  |
| linkcell | convert | 4 | 8 | 7.941 |  |
| linkcell | convert | 8 | 8 | 4.544 |  |
| cell-list | search | 1 | 8 | 492.962 |  |
| cell-list | search | 2 | 8 | 316.537 |  |
| cell-list | search | 4 | 8 | 192.863 |  |
| cell-list | search | 8 | 8 | 132.576 |  |
| linkcell | search | 1 | 8 | 493.166 |  |
| linkcell | search | 2 | 8 | 298.425 |  |
| linkcell | search | 4 | 8 | 177.345 |  |
| linkcell | search | 8 | 8 | 115.945 |  |
