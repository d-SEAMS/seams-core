host bench-3f64f2db4e0d; core 0; Intel(R) Xeon(R) Processor; family 6 model 207 stepping 2; 8 cpus; 1 thread per core; L2 16 MiB (8 instances); L3 320 MiB (1 instance); hypervisor KVM; kernel 6.12.113+
build release -O3; code of each run:

- fas23-ice-ic-r4: code_revision 8fb28dc54dc4e20a67691cd83bae4f8854aaa707; base_revision c0208da6fba7265ac44e279a17c4dffa24e4aef4; engine_tree 753e0a3319ce6847568a878a3085f928f5a80d15; core 0

arms:
- fastatomstruct: fastatomstruct 2.3.0 (build release, commit 4a833a871dd26063381b2183a574b5011ba56eed, tag 2.3.0)
- fastatomstruct_2_2: fastatomstruct 2.2.0 (source pypi)

## rings: wall ms, median (Q1-Q3), n per cell

| system | molecules | arm | n | wall ms | x seams | A/A band | peak MiB | delta MiB | cpu ms | threads | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ice_ic_1000 | 1000 | fastatomstruct | 8 | 24.986 (24.752-25.239) | 15.18 (not single-thread) | 0.879-1.034 | 62.1 | 3.1 | 24.984 | 3 | seams faster |
| ice_ic_1000 | 1000 | fastatomstruct_2_2 | 8 | 24.995 (24.909-25.633) | 15.19 (not single-thread) | 0.879-1.034 | 62.9 | 3.0 | 24.991 | 3 | seams faster |
| ice_ic_1000 | 1000 | seams | 8 | 1.646 (1.611-1.690) | 1.00 (not single-thread) | 0.879-1.034 | 44.7 | 0.9 | 1.644 | 1 |  |
| ice_ic_1000 | 1000 | seams_aa | 8 | 1.641 (1.627-1.668) | 1.00 (not single-thread) | 0.879-1.034 | 44.6 | 0.9 | 1.639 | 1 |  |
| ice_ic_10648 | 10648 | fastatomstruct | 8 | 1771.772 (1765.180-1782.253) | 124.00 (not single-thread) | 0.909-1.029 | 78.6 | 18.2 | 1771.680 | 3 | seams faster |
| ice_ic_10648 | 10648 | fastatomstruct_2_2 | 8 | 1831.757 (1815.645-1851.246) | 128.19 (not single-thread) | 0.909-1.029 | 79.5 | 18.1 | 1831.657 | 3 | seams faster |
| ice_ic_10648 | 10648 | seams | 8 | 14.289 (14.224-14.579) | 1.00 (not single-thread) | 0.909-1.029 | 51.4 | 6.2 | 14.286 | 1 |  |
| ice_ic_10648 | 10648 | seams_aa | 8 | 14.494 (14.268-14.626) | 1.01 (not single-thread) | 0.909-1.029 | 51.4 | 6.2 | 14.490 | 1 |  |
| ice_ic_110592 | 110592 | fastatomstruct | 6 | 603174.076 (586538.102-609036.252) | 4029.61 (not single-thread) | 0.695-1.103 | 249.9 | 175.8 | 603101.132 | 3 | seams faster |
| ice_ic_110592 | 110592 | fastatomstruct_2_2 | 6 | 590552.738 (584425.311-596202.938) | 3945.29 (not single-thread) | 0.695-1.103 | 250.8 | 175.7 | 590484.225 | 3 | seams faster |
| ice_ic_110592 | 110592 | seams | 8 | 149.685 (146.426-150.260) | 1.00 (not single-thread) | 0.695-1.103 | 118.9 | 58.9 | 149.684 | 1 |  |
| ice_ic_110592 | 110592 | seams_aa | 8 | 147.394 (145.057-150.438) | 0.98 (not single-thread) | 0.695-1.103 | 118.9 | 58.9 | 147.392 | 1 |  |

## q6: wall ms, median (Q1-Q3), n per cell

| system | molecules | arm | n | wall ms | x seams | A/A band | peak MiB | delta MiB | cpu ms | threads | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ice_ic_1000 | 1000 | fastatomstruct | 8 | 2.413 (2.399-2.586) | 6.83 (not single-thread) | 0.278-1.050 | 59.9 | 0.9 | 2.411 | 3 | seams faster |
| ice_ic_1000 | 1000 | fastatomstruct_2_2 | 8 | 2.565 (2.550-2.639) | 7.26 (not single-thread) | 0.278-1.050 | 60.9 | 1.0 | 2.563 | 3 | seams faster |
| ice_ic_1000 | 1000 | seams | 8 | 0.353 (0.343-0.391) | 1.00 (not single-thread) | 0.278-1.050 | 44.5 | 0.8 | 0.351 | 1 |  |
| ice_ic_1000 | 1000 | seams_aa | 8 | 0.354 (0.346-0.358) | 1.00 (not single-thread) | 0.278-1.050 | 44.5 | 0.8 | 0.352 | 1 |  |
| ice_ic_10648 | 10648 | fastatomstruct | 8 | 65.730 (65.588-65.923) | 23.17 (not single-thread) | 0.974-1.042 | 62.2 | 1.8 | 65.727 | 3 | seams faster |
| ice_ic_10648 | 10648 | fastatomstruct_2_2 | 8 | 90.862 (89.211-91.326) | 32.02 (not single-thread) | 0.974-1.042 | 63.1 | 1.8 | 90.860 | 3 | seams faster |
| ice_ic_10648 | 10648 | seams | 8 | 2.837 (2.798-2.930) | 1.00 (not single-thread) | 0.974-1.042 | 50.2 | 5.0 | 2.835 | 1 |  |
| ice_ic_10648 | 10648 | seams_aa | 8 | 2.812 (2.788-3.028) | 0.99 (not single-thread) | 0.974-1.042 | 50.1 | 5.0 | 2.809 | 1 |  |
| ice_ic_110592 | 110592 | fastatomstruct | 8 | 5331.071 (5304.075-5443.660) | 172.25 (not single-thread) | 0.873-1.150 | 83.8 | 9.6 | 5330.774 | 3 | seams faster |
| ice_ic_110592 | 110592 | fastatomstruct_2_2 | 8 | 7391.091 (7341.476-7526.439) | 238.81 (not single-thread) | 0.873-1.150 | 84.7 | 9.6 | 7390.775 | 3 | seams faster |
| ice_ic_110592 | 110592 | seams | 8 | 30.950 (28.994-32.945) | 1.00 (not single-thread) | 0.873-1.150 | 106.8 | 46.8 | 30.941 | 1 |  |
| ice_ic_110592 | 110592 | seams_aa | 8 | 32.059 (29.792-32.968) | 1.04 (not single-thread) | 0.873-1.150 | 106.8 | 46.8 | 32.056 | 1 |  |

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
- x seams is one compute thread where compute_threads is 1. The count is taken over the sample from scheduled time in /proc/<pid>/task. An OS thread that stays parked does not count as compute.
- x seams is marked not single-thread on every figure when a competitor did not stay on one compute thread.
- seams/q6: sharp 3.5 A cutoff local Steinhardt q_6 via seams_steinhardt (qlbar skipped)
- fastatomstruct/q6: Steinhardt q_l with near-hard radial weight p1=48 inside 3.5 A
- fastatomstruct_2_2/q6: fastatomstruct 2.2.0 Steinhardt q_l with near-hard radial weight p1=48 inside 3.5 A
- fastatomstruct_2_2/rings: fastatomstruct 2.2.0 King's shortest-path rings for continuity
