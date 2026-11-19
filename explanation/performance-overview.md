<!-- https://docs.dseams.info/explanation/performance-overview.md -->

# Performance overview {#performance-overview}

Each median is printed on one page.

| page | comparison |
| --- | --- |
| Release-candidate sitting | Ideal ice Ic at 110592 and 884736 molecules, one thread and eight. n is 6. Median and type-7 interquartile range for neighbours, q_4, q_6, CHILL+, rings, cages, the row stored as PTM, and King’s count. |
| Head-to-head | Ice Ic and a TIP4P/ice liquid at 884736 molecules, at 1, 2, 4 and 8 threads. n is 6 or 8. Median and interquartile range recomputed from the samples. |
| 4096 to 110592 | Neighbour and ring rows for 4096 to 110592 molecules, kept when n is at least 6, with the recorded interval. |
| Signature thread controls | The hydrate signature search, the eight kernels, and the cutoff ring call split into link, convert and search. |

The operation string is the string the table stored. It is a different call on different pages.

On the release-candidate comparison, $q_4$ and $q_6$ are one Steinhardt call, and the rings row is `seams_frame_rings` on the cutoff graph: the neighbour list and the search together. The 884736-molecule rings row is that same call in the other comparison. The rings kernel on the thread-controls page is `ringNetwork` after the neighbour graph already exists. Those ring times are separate measurements.

The release-candidate PLUMED worker sets `OMP_NUM_THREADS` and does not set `PLUMED_NUM_THREADS`. The head-to-head page holds the other PLUMED $q_6$ rows.
