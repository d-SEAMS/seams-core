# Cutoff neighbour list

Ice Ic, lattice constant 6.35 Å, cutoff 3.5 Å, orthorhombic periodic cell.
Each conventional cell holds 8 molecules. The list is the full both-direction
shell: four neighbours per molecule, `4N` pairs. A pair is kept when the
squared minimum image is strictly below the cutoff squared. Each row is
sorted by cloud index. The same checksum is produced at 1, 2, 4 and 8 threads.

`d-SEAMS` is `orthoFullCsr`: a parallel counting sort into cells, a prefix
sum, and a direct CSR write. Scratch buffers stay allocated and are refilled
on the next call. `rows` in the raw log is `neighListO`, which builds that
CSR and then the atom-ID rows. vesin is 0.6.2 (`NeighborList`, `full_list`,
`n_threads`, quantities `ij`, `copy=False`). freud is 3.6.1 (`AABBQuery`,
ball query, `r_max=3.5`, `exclude_ii`, tree build included). matscipy is
1.3.1 (`neighbour_list('ij')`) with ASE 3.29.0. Every code returned `4N` pairs.

## Hardware

- Kernel `6.12.113+`, KVM guest, one socket
- Intel Xeon, family 6, model 207, stepping 2, 8 cores, 1 thread per core
- L1d 384 KiB, L1i 256 KiB, L2 16 MiB, L3 320 MiB, 47 GiB RAM
- GCC 15.2.0, release build (`meson --buildtype=release`, `-O3`), OpenMP

The guest does not expose hardware performance counters. `perf stat` on the
8-thread, 110592-molecule binary (6 repeats, cores 0–7) reported
`task-clock` 435.64 ms, 3.227 CPUs utilized, 0.135 s elapsed, and
`<not supported>` for cycles, instructions, cache-references and
cache-misses. The phase lines below are the breakdown.

## Median milliseconds

Eight timed calls after two warmups, pinned with `taskset` to cores `0..T-1`.
The median is the 5th of the 8 sorted samples. Lower is faster.

| molecules | threads | d-SEAMS CSR | d-SEAMS rows | vesin | freud | matscipy |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10648 | 1 | 3.531 | 3.859 | 5.674 | 9.508 | 8.458 |
| 10648 | 2 | 1.859 | 2.178 | 3.789 | 9.553 | 8.468 |
| 10648 | 4 | 1.068 | 1.242 | 3.896 | 9.638 | 8.599 |
| 10648 | 8 | 0.737 | 0.889 | 2.375 | 9.575 | 9.215 |
| 110592 | 1 | 25.150 | 29.865 | 52.801 | 104.725 | 76.528 |
| 110592 | 2 | 13.311 | 16.449 | 28.762 | 101.544 | 75.869 |
| 110592 | 4 | 12.746 | 15.758 | 27.467 | 105.529 | 76.589 |
| 110592 | 8 | 4.427 | 5.928 | 17.474 | 101.768 | 75.541 |
| 884736 | 1 | 236.403 | 289.309 | 1073.681 | 1067.929 | 585.680 |
| 884736 | 2 | 121.671 | 163.047 | 576.814 | 1055.421 | 583.512 |
| 884736 | 4 | 61.670 | 91.716 | 323.580 | 1055.389 | 585.517 |
| 884736 | 8 | 33.579 | 48.668 | 176.231 | 1068.058 | 590.335 |

At 8 threads on 110592 molecules the CSR phases were about 1.5 ms to bin,
2.2 ms to count, 0.2–0.7 ms for the prefix, and 2.7 ms to fill and sort.
The atom-ID rows added about 2.7 ms. Binning, the count, and the fill run
on the whole team. The prefix is one pass over the cell totals.

## Rerun

```
meson setup build-release --buildtype=release -Dwith_tests=true \
  -Dwith_cli=true -Dwith_mpi=disabled -Dwith_openmp_offload=disabled
ninja -C build-release tests/bench_cutoff_list
python3 -m venv /tmp/neigh-venv
/tmp/neigh-venv/bin/pip install 'vesin==0.6.2' freud-analysis matscipy ase numpy
NEIGH_PYTHON=/tmp/neigh-venv/bin/python \
  NEIGH_BENCH=$PWD/build-release/tests/bench_cutoff_list \
  /tmp/neigh-venv/bin/python scripts/bench_neighbours_threads.py --reps 8
```

Phase lines, same binary, cores 0–7:

```
OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores \
  SEAMS_NEIGH_PROFILE=1 \
  taskset -c 0-7 build-release/tests/bench_cutoff_list 110592 8 2
```

`bench_cutoff_list` arguments are molecules, threads, repeats.
`SEAMS_NEIGH_PROFILE=1` prints `bin_ms`, `count_ms`, `prefix_ms`,
`fill_ms` and `rows_ms`.

## Raw log

```
molecules=10648 threads=1 pairs=42592 checksum=13782790659416074331 csr_median_ms=3.5307 rows_median_ms=3.8589 csr_samples=3.5773,3.5307,3.4458,3.5823,3.5760,3.4455,3.4481,3.5184 rows_samples=4.0212,4.0064,3.7980,3.7772,3.7964,3.9416,3.8516,3.8589
kind=vesin molecules=10648 threads=1 pairs=42592 median_ms=5.6735 samples_ms=5.5034,5.5333,5.5587,5.6038,5.6735,5.7046,5.7110,5.9183
kind=freud molecules=10648 threads=1 pairs=42592 median_ms=9.5077 samples_ms=9.3886,9.3924,9.4444,9.4478,9.5077,9.5211,9.5658,12.0163
kind=matscipy molecules=10648 threads=1 pairs=42592 median_ms=8.4582 samples_ms=8.3164,8.3390,8.4425,8.4490,8.4582,8.4634,8.4681,8.4955
molecules=10648 threads=2 pairs=42592 checksum=13782790659416074331 csr_median_ms=1.8591 rows_median_ms=2.1776 csr_samples=1.9952,1.8591,1.9889,1.8534,1.8550,1.8552,1.8537,1.9879 rows_samples=2.3440,2.1776,2.2380,2.1423,2.2820,2.1572,2.1478,2.1592
kind=vesin molecules=10648 threads=2 pairs=42592 median_ms=3.7890 samples_ms=3.5138,3.5245,3.6509,3.7695,3.7890,5.6180,5.7539,5.8868
kind=freud molecules=10648 threads=2 pairs=42592 median_ms=9.5526 samples_ms=9.4330,9.4365,9.5329,9.5511,9.5526,9.6211,9.6587,9.6917
kind=matscipy molecules=10648 threads=2 pairs=42592 median_ms=8.4682 samples_ms=8.3435,8.3759,8.4426,8.4594,8.4682,8.4853,8.4856,8.5105
molecules=10648 threads=4 pairs=42592 checksum=13782790659416074331 csr_median_ms=1.0681 rows_median_ms=1.2423 csr_samples=1.1991,1.0153,1.0396,1.0344,1.2065,1.0359,1.0681,1.0838 rows_samples=1.6745,1.2535,1.2423,1.2360,1.2422,1.2823,1.2352,1.2232
kind=vesin molecules=10648 threads=4 pairs=42592 median_ms=3.8958 samples_ms=2.8777,2.8801,3.0309,3.7719,3.8958,4.0966,4.2211,4.3461
kind=freud molecules=10648 threads=4 pairs=42592 median_ms=9.6375 samples_ms=9.5529,9.6001,9.6284,9.6288,9.6375,9.6922,9.7302,9.9327
kind=matscipy molecules=10648 threads=4 pairs=42592 median_ms=8.5988 samples_ms=8.4818,8.4831,8.5278,8.5683,8.5988,8.6281,8.6367,8.6830
molecules=10648 threads=8 pairs=42592 checksum=13782790659416074331 csr_median_ms=0.7371 rows_median_ms=0.8891 csr_samples=0.7371,0.6868,0.6757,0.7702,0.7787,0.7132,0.7398,0.7220 rows_samples=1.9136,0.8768,0.8777,0.8858,0.8852,0.8891,0.9157,0.9022
kind=vesin molecules=10648 threads=8 pairs=42592 median_ms=2.3754 samples_ms=2.0775,2.2333,2.3533,2.3753,2.3754,2.8719,5.2864,5.5782
kind=freud molecules=10648 threads=8 pairs=42592 median_ms=9.5748 samples_ms=9.4788,9.4873,9.5231,9.5580,9.5748,9.5771,9.5871,9.7791
kind=matscipy molecules=10648 threads=8 pairs=42592 median_ms=9.2152 samples_ms=8.7930,8.8023,8.9622,9.0832,9.2152,9.2929,9.3009,9.4009
molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 csr_median_ms=25.1502 rows_median_ms=29.8652 csr_samples=26.5322,28.2799,25.2378,25.1094,24.9398,25.1502,25.1400,24.8695 rows_samples=32.7555,32.2879,29.8652,29.5509,29.4537,29.8547,29.8791,29.4933
kind=vesin molecules=110592 threads=1 pairs=442368 median_ms=52.8007 samples_ms=49.6497,52.1663,52.4496,52.6253,52.8007,54.5748,55.2475,58.6171
kind=freud molecules=110592 threads=1 pairs=442368 median_ms=104.7253 samples_ms=101.7118,102.4326,102.6009,103.0954,104.7253,108.5395,110.1960,110.9510
kind=matscipy molecules=110592 threads=1 pairs=442368 median_ms=76.5279 samples_ms=75.7622,75.8363,75.9660,76.3268,76.5279,76.7310,79.4174,79.9316
molecules=110592 threads=2 pairs=442368 checksum=8222705064348099955 csr_median_ms=13.3107 rows_median_ms=16.4493 csr_samples=13.8094,13.2682,12.9802,13.0746,13.3107,13.3408,13.3729,13.0102 rows_samples=19.0968,16.1053,16.4493,16.3297,16.0897,16.6797,16.3769,16.4805
kind=vesin molecules=110592 threads=2 pairs=442368 median_ms=28.7622 samples_ms=28.3802,28.4056,28.6377,28.6575,28.7622,28.8184,29.0074,29.3780
kind=freud molecules=110592 threads=2 pairs=442368 median_ms=101.5435 samples_ms=100.6077,100.9050,100.9881,101.1397,101.5435,102.3015,102.3181,103.4454
kind=matscipy molecules=110592 threads=2 pairs=442368 median_ms=75.8692 samples_ms=75.7009,75.7485,75.8378,75.8527,75.8692,76.2013,76.6119,78.1870
molecules=110592 threads=4 pairs=442368 checksum=8222705064348099955 csr_median_ms=12.7459 rows_median_ms=15.7581 csr_samples=12.7459,13.1396,12.5321,12.6463,13.2115,13.3651,12.6128,7.6835 rows_samples=18.5690,15.5339,15.5194,15.6537,15.8130,16.3161,15.7581,9.6465
kind=vesin molecules=110592 threads=4 pairs=442368 median_ms=27.4665 samples_ms=26.6767,27.3047,27.3272,27.3398,27.4665,27.5308,27.5656,27.6411
kind=freud molecules=110592 threads=4 pairs=442368 median_ms=105.5292 samples_ms=104.3341,104.6533,105.3081,105.4954,105.5292,105.7560,108.0639,108.1220
kind=matscipy molecules=110592 threads=4 pairs=442368 median_ms=76.5885 samples_ms=76.0297,76.0418,76.2561,76.4843,76.5885,76.6478,76.7380,76.8230
molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 csr_median_ms=4.4267 rows_median_ms=5.9277 csr_samples=4.5809,5.7288,4.4129,4.3229,4.4267,4.3235,4.5058,4.2508 rows_samples=10.0468,7.5728,5.9277,5.8058,5.4980,7.3309,5.5193,5.4890
kind=vesin molecules=110592 threads=8 pairs=442368 median_ms=17.4740 samples_ms=16.1485,16.5963,16.7936,17.0000,17.4740,17.5138,17.5542,19.1365
kind=freud molecules=110592 threads=8 pairs=442368 median_ms=101.7675 samples_ms=100.7947,101.1494,101.1664,101.5857,101.7675,101.9745,102.1099,102.7809
kind=matscipy molecules=110592 threads=8 pairs=442368 median_ms=75.5413 samples_ms=75.3469,75.4681,75.5288,75.5298,75.5413,75.6579,75.9823,76.9681
molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 csr_median_ms=236.4034 rows_median_ms=289.3091 csr_samples=235.8964,236.4034,237.4392,236.8145,235.2647,235.0685,245.5089,234.2122 rows_samples=298.6153,289.3091,286.2176,285.3170,285.5766,293.7344,295.2892,287.5968
kind=vesin molecules=884736 threads=1 pairs=3538944 median_ms=1073.6805 samples_ms=1071.8321,1071.8360,1073.2135,1073.2625,1073.6805,1074.1097,1076.2066,1078.5157
kind=freud molecules=884736 threads=1 pairs=3538944 median_ms=1067.9293 samples_ms=1056.2978,1058.9365,1063.2766,1064.7146,1067.9293,1069.4364,1072.3539,1074.2250
kind=matscipy molecules=884736 threads=1 pairs=3538944 median_ms=585.6797 samples_ms=580.6717,583.4579,585.1193,585.4606,585.6797,585.7660,585.8462,585.9333
molecules=884736 threads=2 pairs=3538944 checksum=9316061812838145587 csr_median_ms=121.6709 rows_median_ms=163.0474 csr_samples=125.9350,147.7060,119.3425,123.3162,121.0451,119.3861,119.4294,121.6709 rows_samples=171.5983,158.4196,163.5711,157.9367,164.6845,157.7758,163.0474,160.9659
kind=vesin molecules=884736 threads=2 pairs=3538944 median_ms=576.8143 samples_ms=570.5083,574.0647,576.3226,576.5264,576.8143,577.0379,578.3187,579.5403
kind=freud molecules=884736 threads=2 pairs=3538944 median_ms=1055.4212 samples_ms=1049.6374,1050.5901,1051.4172,1055.4062,1055.4212,1059.5470,1063.1849,1067.4778
kind=matscipy molecules=884736 threads=2 pairs=3538944 median_ms=583.5115 samples_ms=580.5139,580.6043,581.1850,582.5127,583.5115,584.7398,585.5645,588.8628
molecules=884736 threads=4 pairs=3538944 checksum=9316061812838145587 csr_median_ms=61.6701 rows_median_ms=91.7156 csr_samples=60.5885,60.7146,61.5494,61.3836,63.3247,62.4305,61.6701,61.7357 rows_samples=102.7189,90.1946,88.5884,94.1006,117.2673,91.6014,89.1074,91.7156
kind=vesin molecules=884736 threads=4 pairs=3538944 median_ms=323.5803 samples_ms=322.6266,322.8458,322.9298,323.5248,323.5803,323.8637,325.1847,326.1428
kind=freud molecules=884736 threads=4 pairs=3538944 median_ms=1055.3889 samples_ms=1044.2061,1049.1768,1049.4911,1051.0277,1055.3889,1056.0477,1056.7326,1057.5254
kind=matscipy molecules=884736 threads=4 pairs=3538944 median_ms=585.5168 samples_ms=580.6758,582.6148,583.9991,584.3550,585.5168,585.7258,590.7514,590.9497
molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 csr_median_ms=33.5789 rows_median_ms=48.6675 csr_samples=41.9785,34.1191,33.5789,33.0317,32.7467,33.2250,32.3763,34.3708 rows_samples=73.0009,57.1652,47.0358,48.4149,48.6675,46.4784,52.0003,47.7116
kind=vesin molecules=884736 threads=8 pairs=3538944 median_ms=176.2310 samples_ms=171.9261,172.4380,173.2194,173.8490,176.2310,178.1876,180.0601,199.2978
kind=freud molecules=884736 threads=8 pairs=3538944 median_ms=1068.0580 samples_ms=1062.5830,1065.3066,1065.3198,1065.8386,1068.0580,1069.4406,1073.7921,1082.6855
kind=matscipy molecules=884736 threads=8 pairs=3538944 median_ms=590.3347 samples_ms=585.8217,586.1255,586.7858,587.2295,590.3347,591.0257,593.1093,614.8044
```

Phase sample, 110592 molecules, 8 threads:

```
neigh_profile n=110592 threads=8 pairs=442368 bin_ms=0.9450 count_ms=2.1503 prefix_ms=0.1514 fill_ms=3.1931 rows_ms=2.7216
```
