host bench-3f64f2db4e0d; core 1; Intel(R) Xeon(R) Processor; family 6 model 207 stepping 2; 8 cpus; 1 thread per core; L2 16 MiB (8 instances); L3 320 MiB (1 instance); hypervisor KVM; kernel 6.12.113+
build release -O3; code of each run:

- fas23-ice-ic: code_revision a726636b8fbf0baa368117c3a904c0a373abd4b0; base_revision c0208da6fba7265ac44e279a17c4dffa24e4aef4; engine_tree 9ecd000c053a07c1498f9b41092f16996a58b8e2; core 1

arms:
- fastatomstruct: fastatomstruct 2.3.0 (build release, commit 4a833a871dd26063381b2183a574b5011ba56eed, tag 2.3.0)
- fastatomstruct_2_2: fastatomstruct 2.2.0 (source pypi)

## rings: wall ms, median (Q1-Q3), n per cell

| system | molecules | arm | n | wall ms | x seams | A/A band | peak MiB | delta MiB | cpu ms | threads | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ice_ic_1000 | 1000 | fastatomstruct | 8 | 24.829 (24.723-25.077) | 15.24 | 0.601-1.091 | 61.4 | 3.1 | 24.817 | 3 | seams faster |
| ice_ic_1000 | 1000 | fastatomstruct_2_2 | 8 | 24.847 (24.808-24.949) | 15.25 | 0.601-1.091 | 62.2 | 3.0 | 24.844 | 3 | seams faster |
| ice_ic_1000 | 1000 | seams | 8 | 1.630 (1.624-1.667) | 1.00 | 0.601-1.091 | 43.9 | 0.9 | 1.628 | 1 |  |
| ice_ic_1000 | 1000 | seams_aa | 8 | 1.666 (1.650-1.692) | 1.02 | 0.601-1.091 | 43.9 | 0.9 | 1.664 | 1 |  |
| ice_ic_10648 | 10648 | fastatomstruct | 8 | 1764.147 (1757.507-1803.799) | 123.82 | 0.969-1.020 | 77.9 | 18.2 | 1763.995 | 3 | seams faster |
| ice_ic_10648 | 10648 | fastatomstruct_2_2 | 8 | 1827.344 (1805.975-1854.690) | 128.26 | 0.969-1.020 | 78.8 | 18.2 | 1827.244 | 3 | seams faster |
| ice_ic_10648 | 10648 | seams | 8 | 14.247 (14.091-14.455) | 1.00 | 0.969-1.020 | 50.7 | 6.2 | 14.245 | 1 |  |
| ice_ic_10648 | 10648 | seams_aa | 8 | 14.429 (14.246-14.587) | 1.01 | 0.969-1.020 | 50.7 | 6.2 | 14.426 | 1 |  |
| ice_ic_110592 | 110592 | fastatomstruct | 6 | 565560.528 (561836.539-571198.541) | 3820.38 | 0.551-1.116 | 249.2 | 175.7 | 565483.798 | 3 | seams faster |
| ice_ic_110592 | 110592 | fastatomstruct_2_2 | 6 | 586837.255 (574290.749-593809.651) | 3964.11 | 0.551-1.116 | 250.0 | 175.7 | 586746.077 | 3 | seams faster |
| ice_ic_110592 | 110592 | seams | 8 | 148.038 (146.383-151.271) | 1.00 | 0.551-1.116 | 118.2 | 58.9 | 148.031 | 1 |  |
| ice_ic_110592 | 110592 | seams_aa | 8 | 152.626 (147.753-156.756) | 1.03 | 0.551-1.116 | 118.2 | 58.9 | 152.624 | 1 |  |

## q6: wall ms, median (Q1-Q3), n per cell

| system | molecules | arm | n | wall ms | x seams | A/A band | peak MiB | delta MiB | cpu ms | threads | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ice_ic_1000 | 1000 | fastatomstruct | 8 | 2.749 (2.422-4.472) | 4.97 | 0.124-1.114 | 59.2 | 0.9 | 2.747 | 3 | seams faster |
| ice_ic_1000 | 1000 | fastatomstruct_2_2 | 8 | 2.741 (2.595-3.370) | 4.95 | 0.124-1.114 | 60.1 | 1.0 | 2.739 | 3 | seams faster |
| ice_ic_1000 | 1000 | seams | 8 | 0.553 (0.347-1.353) | 1.00 | 0.124-1.114 | 43.8 | 0.8 | 0.551 | 1 |  |
| ice_ic_1000 | 1000 | seams_aa | 8 | 0.364 (0.350-0.371) | 0.66 | 0.124-1.114 | 43.8 | 0.8 | 0.362 | 1 |  |
| ice_ic_10648 | 10648 | fastatomstruct | 8 | 66.052 (65.498-66.264) | 22.15 | 0.926-1.054 | 61.5 | 1.8 | 66.049 | 3 | seams faster |
| ice_ic_10648 | 10648 | fastatomstruct_2_2 | 8 | 90.612 (89.058-91.815) | 30.38 | 0.926-1.054 | 62.4 | 1.8 | 90.589 | 3 | seams faster |
| ice_ic_10648 | 10648 | seams | 8 | 2.983 (2.865-3.090) | 1.00 | 0.926-1.054 | 49.5 | 5.0 | 2.980 | 1 |  |
| ice_ic_10648 | 10648 | seams_aa | 8 | 2.989 (2.906-3.090) | 1.00 | 0.926-1.054 | 49.5 | 5.0 | 2.986 | 1 |  |
| ice_ic_110592 | 110592 | fastatomstruct | 8 | 5322.907 (5301.341-5364.740) | 177.17 | 0.819-1.224 | 83.0 | 9.6 | 5322.214 | 3 | seams faster |
| ice_ic_110592 | 110592 | fastatomstruct_2_2 | 8 | 7282.628 (7196.485-7382.057) | 242.40 | 0.819-1.224 | 84.0 | 9.6 | 7282.347 | 3 | seams faster |
| ice_ic_110592 | 110592 | seams | 8 | 30.044 (29.588-33.331) | 1.00 | 0.819-1.224 | 106.1 | 46.8 | 30.037 | 1 |  |
| ice_ic_110592 | 110592 | seams_aa | 8 | 32.846 (30.587-33.348) | 1.09 | 0.819-1.224 | 106.1 | 46.8 | 32.843 | 1 |  |

## agreement with seams (same cell, round 0 output)

| workload | system | arm | equal | detail |
|---|---|---|---|---|
| q6 | ice_ic_1000 | fastatomstruct | yes | n=1000; mean=0.62853936; max_abs=3.378e-09; atol=1e-06; vs 2.2.0 yes (max_abs=0.000e+00) |
| q6 | ice_ic_1000 | fastatomstruct_2_2 | yes | n=1000; mean=0.62853936; max_abs=3.378e-09; atol=1e-06 |
| q6 | ice_ic_10648 | fastatomstruct | yes | n=10648; mean=0.62853936; max_abs=3.378e-09; atol=1e-06; vs 2.2.0 yes (max_abs=0.000e+00) |
| q6 | ice_ic_10648 | fastatomstruct_2_2 | yes | n=10648; mean=0.62853936; max_abs=3.378e-09; atol=1e-06 |
| q6 | ice_ic_110592 | fastatomstruct | yes | n=110592; mean=0.62853936; max_abs=3.378e-09; atol=1e-06; vs 2.2.0 yes (max_abs=0.000e+00) |
| q6 | ice_ic_110592 | fastatomstruct_2_2 | yes | n=110592; mean=0.62853936; max_abs=3.378e-09; atol=1e-06 |
| rings | ice_ic_1000 | fastatomstruct | yes | sizes [3, 4, 5, 6]; vs 2.2.0 yes |
| rings | ice_ic_1000 | fastatomstruct_2_2 | yes | sizes [3, 4, 5, 6] |
| rings | ice_ic_10648 | fastatomstruct | yes | sizes [3, 4, 5, 6]; vs 2.2.0 yes |
| rings | ice_ic_10648 | fastatomstruct_2_2 | yes | sizes [3, 4, 5, 6] |
| rings | ice_ic_110592 | fastatomstruct | yes | sizes [3, 4, 5, 6]; vs 2.2.0 yes |
| rings | ice_ic_110592 | fastatomstruct_2_2 | yes | sizes [3, 4, 5, 6] |

## samples that did not complete

| workload | system | arm | samples | status | detail |
|---|---|---|---|---|---|
| rings | ice_ic_110592 | fastatomstruct | 2 | skipped | capped at 6 rounds |
| rings | ice_ic_110592 | fastatomstruct_2_2 | 2 | skipped | capped at 6 rounds |

## what the timed call covers

- thread controls: every sample process is pinned to one CPU core via taskset with single-thread env vars (OMP, MKL, OPENBLAS, NUMBA, TBB, RAYON, POLARS, TOKIO, NUMEXPR, OVITO set to 1) and the tool's 1-thread API knob.
- threads: observed OS threads in /proc/self/status after the timed call: 96/188 good samples had 1 OS thread; 92/188 good samples had >1 OS thread.
- multi-thread tools (os_threads > 1 despite 1-thread controls and taskset): fastatomstruct (3).
- seams/q6: sharp 3.5 A cutoff local Steinhardt q_6 via seams_steinhardt (qlbar skipped)
- fastatomstruct/q6: Steinhardt q_l with near-hard radial weight p1=48 inside 3.5 A
- fastatomstruct_2_2/q6: fastatomstruct 2.2.0 Steinhardt q_l with near-hard radial weight p1=48 inside 3.5 A
- fastatomstruct_2_2/rings: fastatomstruct 2.2.0 King's shortest-path rings for continuity
