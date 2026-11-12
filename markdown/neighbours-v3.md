# Triclinic list, CHILL+, diamond CNA, and the hop ring search

Ice Ic, lattice constant 6.35 Å. The triclinic cell keeps that cube and
adds xy = 0.25 L, xz = 0.05 L, yz = 0.12 L. A pair is kept when the
squared minimum image is strictly below the cutoff squared. Every list
here uses 3.5 Å. Eight timed calls after one warmup. The median is the
5th of the 8 sorted samples. Threads are pinned with `taskset` to cores
`0..T-1`, with `OMP_PROC_BIND=close`, `OMP_PLACES=cores`,
`OMP_DYNAMIC=FALSE`, and `OMP_WAIT_POLICY=PASSIVE`.

Same machine as `markdown/neighbours.md`: kernel 6.12.113+, KVM guest,
Intel Xeon family 6 model 207, 8 cores, 2.4 GHz, L2 16 MiB, L3 320 MiB,
GCC 15.2.0, release `-O3`. vesin 0.6.2, freud 3.6.1, mdapy 2.0.0, OVITO
3.17.2. OVITO is run with `LD_LIBRARY_PATH` pointed at the system libGL
directories. That variable is unset for d-SEAMS. Each OVITO sample
builds a new pipeline, because a pipeline caches its output.

## Why the triclinic list was slow

At 110592 molecules the perpendicular widths give 42 x 43 x 43 cells
against 43 cubed for the orthorhombic box. The reach is 1 on every
axis, so both lists use the 27-cell stencil. The walks test 4,297,540
and 4,204,334 candidate pairs. Occupancy matches: most cells hold one
or two molecules, and none holds more than five. The cell count and the
skew of the bins are not the gap.

The previous 43 ms was the cell matrix, multiplied for every candidate.
With `SEAMS_NEIGH_PROFILE=1` that split was about 18 ms in the count
pass and 25 ms in the fill pass, against about 10 ms and 14 ms for the
orthorhombic loop. The 27 cartesian images of a fractional shift
`(sx, sy, sz)` are now filled once. Both boxes call that same loop.
On this binary the orthorhombic list at 110592 molecules and 1 thread
is 12.09 ms and the triclinic list is 32.99 ms. At 8 threads the pair
is 4.00 ms and 6.58 ms.

At 884736 molecules the triclinic list is 462 ms against 116 ms
orthorhombic. The same stencil is in use. The working set of the three
coordinate streams no longer fits in L2, and the sheared coordinates
miss more often.

## Triclinic neighbour list

Every code returned `4N` pairs. The index checksum matches across thread
counts and matches the orthorhombic list: this shear does not change
which molecules are neighbours.

| molecules | threads | d-SEAMS | vesin | freud | orthorhombic |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10648 | 1 | 3.614 | 6.909 | 13.563 | 1.591 |
| 10648 | 8 | 0.905 | 6.497 | 3.878 | 0.500 |
| 110592 | 1 | 32.986 | 58.369 | 176.267 | 12.089 |
| 110592 | 8 | 6.580 | 61.736 | 49.393 | 4.002 |
| 884736 | 1 | 461.745 | 1284.579 | 1425.870 | 115.950 |
| 884736 | 8 | 75.689 | 1335.178 | 459.751 | 31.276 |

The orthorhombic column is the same binary, run after the triclinic
session. At 110592 molecules and 8 threads the orthorhombic cutoff
list is 4.00 ms and vesin on the same orthorhombic cell is 50.20 ms.

## CHILL+

Each d-SEAMS sample builds the CSR and then runs the analysis on those
rows. CHILL+ keeps the four nearest neighbours of each row and
accumulates the real spherical harmonics at l = 3, the same contraction
as the fused kernel. Every molecule of this ice is cubic. freud 3.6.1
has no CHILL+.

| molecules | threads | d-SEAMS | mdapy | OVITO |
| ---: | ---: | ---: | ---: | ---: |
| 110592 | 1 | 28.635 | 67.061 | 658.972 |
| 110592 | 8 | 7.512 | 88.905 | 665.308 |
| 884736 | 1 | 325.990 | 809.027 | 5371.775 |
| 884736 | 8 | 100.510 | 895.189 | 5359.389 |

## Diamond common-neighbour test

The symmetry column is the diamond common-neighbour test on the 3.5 Å
CSR. A cubic shell is the twelve outward bonds of the four nearest
neighbours. It is not the template overlay. PTM has no cubic-ice
template: on this crystal every structure code is 0, in mdapy and in
OVITO, which is why those sinks are 0. Identify diamond and mdapy
`DiamondStructure` label every molecule, and so does this kernel.

| molecules | threads | d-SEAMS | mdapy diamond | OVITO identify diamond | mdapy PTM | OVITO PTM |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 110592 | 1 | 39.342 | 137.849 | 143.585 | 1268.823 | 1500.891 |
| 110592 | 8 | 11.491 | 186.123 | 151.964 | 1321.445 | 1489.324 |
| 884736 | 1 | 378.629 | 1134.471 | 1155.594 | 9906.975 | 11728.495 |
| 884736 | 8 | 95.424 | 1264.562 | 1184.624 | 10100.536 | 11854.570 |

## Rings

The census is the short-hop search: for rings of size at most 7 the
adjacency answers the primitivity test, and the ball table is not
built. The CSR is handed to that search. The counts match mdapy:
221184 rings at 110592 molecules and 1769472 at 884736. The OVITO ring
sink is the particle count. `FindRings` does not write a structure
code. freud has no ring census.

| molecules | threads | d-SEAMS | mdapy | OVITO |
| ---: | ---: | ---: | ---: | ---: |
| 110592 | 1 | 92.590 | 2764.487 | 562.649 |
| 110592 | 8 | 37.075 | 2996.590 | 594.795 |
| 884736 | 1 | 868.278 | 22539.341 | 4408.334 |
| 884736 | 8 | 288.024 | 21541.297 | 4552.499 |

## How to rerun

From the repository root, with the release build and the comparison
environment on `PATH` via `NEIGH_PYTHON`:

```
NEIGH_PYTHON="$TMPDIR/neigh-venv/bin/python" \
NEIGH_BENCH=bbdir-release/tests/bench_v2 \
python scripts/bench_neighbours_v2.py
```

The orthorhombic rows are `bench_v2 ortho <molecules> <threads> 8`
under the same pin.

## Raw log

```
mode=tilt molecules=10648 threads=1 pairs=42592 checksum=13782790659416074331 sink=0 median_ms=3.6136 samples_ms=3.5290,3.5116,3.6462,3.5746,3.6136,3.7374,3.6805,3.5818
kind=tilt-vesin molecules=10648 threads=1 sink=42592 median_ms=6.9089 samples_ms=6.1481,6.3096,6.3132,6.5718,6.9089,6.9213,7.0270,7.2705
kind=tilt-freud molecules=10648 threads=1 sink=42592 median_ms=13.5628 samples_ms=13.4089,13.4204,13.5068,13.5177,13.5628,13.5660,13.9384,15.0522
mode=tilt molecules=10648 threads=8 pairs=42592 checksum=13782790659416074331 sink=0 median_ms=0.9050 samples_ms=1.0675,0.8831,0.9349,0.9165,0.8946,0.8785,0.8835,0.9050
kind=tilt-vesin molecules=10648 threads=8 sink=42592 median_ms=6.4968 samples_ms=6.3332,6.4523,6.4796,6.4959,6.4968,6.5081,6.5400,6.5757
kind=tilt-freud molecules=10648 threads=8 sink=42592 median_ms=3.8783 samples_ms=3.7471,3.7902,3.7990,3.8378,3.8783,3.9088,4.1633,4.1646
mode=tilt molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=0 median_ms=32.9859 samples_ms=33.6284,33.6802,32.9859,32.7138,32.1201,31.9763,32.8132,33.7346
kind=tilt-vesin molecules=110592 threads=1 sink=442368 median_ms=58.3694 samples_ms=56.6844,57.0274,57.0646,57.5344,58.3694,59.3172,59.3271,60.0721
kind=tilt-freud molecules=110592 threads=1 sink=442368 median_ms=176.2666 samples_ms=166.8276,169.7186,170.0775,171.6822,176.2666,176.9997,182.3942,182.9749
mode=tilt molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=0 median_ms=6.5796 samples_ms=7.2042,6.6091,6.5796,6.2714,6.2174,7.4520,6.1244,6.3126
kind=tilt-vesin molecules=110592 threads=8 sink=442368 median_ms=61.7357 samples_ms=58.1809,58.8971,59.0076,61.1365,61.7357,63.5179,68.6842,75.1982
kind=tilt-freud molecules=110592 threads=8 sink=442368 median_ms=49.3934 samples_ms=48.4543,48.6769,48.8791,49.2567,49.3934,49.4176,50.1489,50.2504
mode=tilt molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=0 median_ms=461.7447 samples_ms=490.6386,467.6978,471.7271,461.7447,453.8932,447.8837,456.4112,458.8929
kind=tilt-vesin molecules=884736 threads=1 sink=3538944 median_ms=1284.5788 samples_ms=1239.3421,1239.5664,1247.1746,1272.9547,1284.5788,1290.4438,1297.5934,1312.3276
kind=tilt-freud molecules=884736 threads=1 sink=3538944 median_ms=1425.8696 samples_ms=1380.6571,1382.9046,1398.9587,1421.4590,1425.8696,1426.7651,1429.9934,1455.9621
mode=tilt molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=0 median_ms=75.6887 samples_ms=82.9710,75.6887,73.3067,78.5547,71.9985,74.4027,71.1172,77.1862
kind=tilt-vesin molecules=884736 threads=8 sink=3538944 median_ms=1335.1776 samples_ms=1283.1001,1304.3266,1308.8527,1324.4148,1335.1776,1347.1455,1357.0829,1453.4672
kind=tilt-freud molecules=884736 threads=8 sink=3538944 median_ms=459.7513 samples_ms=452.9957,453.2066,454.1927,454.5660,459.7513,459.9613,473.1465,487.2203
mode=chill molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=110592 median_ms=28.6346 samples_ms=31.6623,28.6346,28.3681,27.9607,27.8725,27.7006,29.2106,28.6726
kind=mdapy-chill molecules=110592 threads=1 sink=110592 median_ms=67.0605 samples_ms=65.6662,66.3215,66.3976,66.7692,67.0605,67.3278,67.4323,67.4897
kind=ovito-chill molecules=110592 threads=1 sink=110592 median_ms=658.9716 samples_ms=645.9547,647.1954,652.1019,654.3331,658.9716,663.6213,681.9221,693.8544
kind=freud-none molecules=110592 threads=1 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=chill molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=110592 median_ms=7.5116 samples_ms=14.6228,9.8933,7.2321,7.3881,7.5116,8.2223,6.7308,6.4281
kind=mdapy-chill molecules=110592 threads=8 sink=110592 median_ms=88.9052 samples_ms=87.1070,88.0382,88.1475,88.2820,88.9052,89.1076,89.9463,92.0046
kind=ovito-chill molecules=110592 threads=8 sink=110592 median_ms=665.3082 samples_ms=656.8273,657.1579,662.1448,662.8558,665.3082,673.7192,675.3397,699.7180
kind=freud-none molecules=110592 threads=8 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=chill molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=884736 median_ms=325.9905 samples_ms=325.9905,326.5068,325.5837,328.0363,323.2055,318.3777,318.1229,332.4789
kind=mdapy-chill molecules=884736 threads=1 sink=884736 median_ms=809.0274 samples_ms=793.1692,793.7185,799.6084,807.9949,809.0274,813.6752,825.7008,833.6796
kind=ovito-chill molecules=884736 threads=1 sink=884736 median_ms=5371.7746 samples_ms=5143.1507,5149.6017,5201.3209,5326.0036,5371.7746,5492.5485,5592.9248,5683.6245
kind=freud-none molecules=884736 threads=1 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=chill molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=884736 median_ms=100.5104 samples_ms=108.0509,100.6428,99.3612,96.6862,96.4581,98.3457,100.5104,101.9906
kind=mdapy-chill molecules=884736 threads=8 sink=884736 median_ms=895.1890 samples_ms=862.0503,879.1258,880.0615,882.4211,895.1890,898.4240,904.8394,950.4902
kind=ovito-chill molecules=884736 threads=8 sink=884736 median_ms=5359.3892 samples_ms=5275.1901,5313.7701,5337.1574,5355.1598,5359.3892,5557.7074,5582.1964,5824.3119
kind=freud-none molecules=884736 threads=8 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=rings molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=221184 median_ms=92.5899 samples_ms=92.5899,93.3141,90.9893,89.9044,91.7381,91.8897,93.3897,92.8667
kind=mdapy-rings molecules=110592 threads=1 sink=221184 median_ms=2764.4871 samples_ms=2707.5034,2733.1687,2742.4816,2759.1961,2764.4871,2806.3040,2812.8644,2813.4336
kind=ovito-rings molecules=110592 threads=1 sink=110592 median_ms=562.6486 samples_ms=551.8617,554.5727,556.7448,558.4628,562.6486,572.7330,579.3828,611.7632
mode=rings molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=221184 median_ms=37.0745 samples_ms=39.0109,38.1296,37.2926,36.4884,35.9699,37.0745,36.1229,36.2109
kind=mdapy-rings molecules=110592 threads=8 sink=221184 median_ms=2996.5902 samples_ms=2898.5053,2920.9022,2934.1534,2962.4218,2996.5902,3068.3538,3110.0272,3141.9184
kind=ovito-rings molecules=110592 threads=8 sink=110592 median_ms=594.7947 samples_ms=574.2233,578.0247,580.0843,585.8687,594.7947,599.6228,601.2701,605.5221
mode=rings molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=1769472 median_ms=868.2783 samples_ms=870.1243,868.2783,847.8250,888.6081,876.7560,828.6685,857.1269,850.7157
kind=mdapy-rings molecules=884736 threads=1 sink=1769472 median_ms=22539.3410 samples_ms=21728.2251,22223.0674,22253.7336,22386.8878,22539.3410,22657.5039,23195.5240,23229.4331
kind=ovito-rings molecules=884736 threads=1 sink=884736 median_ms=4408.3337 samples_ms=4291.5874,4324.2477,4329.1406,4366.5833,4408.3337,4425.7142,4586.1416,4652.5332
mode=rings molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=1769472 median_ms=288.0241 samples_ms=281.1083,287.9801,285.7423,288.0241,289.6980,287.3453,294.1276,288.8472
kind=mdapy-rings molecules=884736 threads=8 sink=1769472 median_ms=21541.2973 samples_ms=21404.8060,21405.5801,21455.0717,21537.3186,21541.2973,21563.3928,21577.5517,21735.2793
kind=ovito-rings molecules=884736 threads=8 sink=884736 median_ms=4552.4987 samples_ms=4397.8127,4481.7712,4488.4843,4495.8304,4552.4987,4574.8770,4650.0623,4685.9715
mode=symmetry molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=110592 median_ms=39.3419 samples_ms=39.3419,39.2426,39.2951,38.7307,38.9384,39.3683,39.8771,39.5112
kind=mdapy-symmetry molecules=110592 threads=1 sink=0 median_ms=1268.8226 samples_ms=1260.8496,1265.0466,1265.9495,1267.3792,1268.8226,1273.2153,1276.3912,1322.8777
kind=ovito-symmetry molecules=110592 threads=1 sink=0 median_ms=1500.8914 samples_ms=1471.2850,1473.5622,1474.2113,1479.8385,1500.8914,1501.0793,1501.7458,1502.4293
kind=mdapy-diamond molecules=110592 threads=1 sink=110592 median_ms=137.8485 samples_ms=137.4763,137.6174,137.6323,137.7066,137.8485,137.8677,137.9829,138.8264
kind=ovito-diamond molecules=110592 threads=1 sink=110592 median_ms=143.5853 samples_ms=142.3339,142.9965,143.0891,143.1635,143.5853,143.6241,143.6911,143.7238
mode=symmetry molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=110592 median_ms=11.4913 samples_ms=11.6319,11.2097,11.4913,11.6499,11.4546,11.5239,11.4674,11.1924
kind=mdapy-symmetry molecules=110592 threads=8 sink=0 median_ms=1321.4445 samples_ms=1316.3119,1318.7064,1319.2738,1319.3577,1321.4445,1321.9608,1322.8685,1324.7801
kind=ovito-symmetry molecules=110592 threads=8 sink=0 median_ms=1489.3236 samples_ms=1481.4217,1486.6203,1488.0240,1489.2451,1489.3236,1490.2243,1496.5122,1503.4503
kind=mdapy-diamond molecules=110592 threads=8 sink=110592 median_ms=186.1225 samples_ms=185.1389,185.3903,185.7019,186.0050,186.1225,186.2390,186.8598,187.5554
kind=ovito-diamond molecules=110592 threads=8 sink=110592 median_ms=151.9645 samples_ms=143.9273,144.4579,145.2344,151.9187,151.9645,152.0143,152.0567,152.7948
mode=symmetry molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=884736 median_ms=378.6286 samples_ms=372.1380,374.8579,596.4697,581.9226,378.6286,378.3131,379.9081,374.4870
kind=mdapy-symmetry molecules=884736 threads=1 sink=0 median_ms=9906.9748 samples_ms=9890.7389,9893.5175,9898.6585,9901.6417,9906.9748,9907.7464,9916.7445,9917.9189
kind=ovito-symmetry molecules=884736 threads=1 sink=0 median_ms=11728.4951 samples_ms=11701.2404,11707.4471,11707.4914,11716.6456,11728.4951,11731.2614,11733.3916,11748.2163
kind=mdapy-diamond molecules=884736 threads=1 sink=884736 median_ms=1134.4712 samples_ms=1131.4855,1131.5384,1132.8975,1132.8991,1134.4712,1135.0274,1135.6548,1136.4672
kind=ovito-diamond molecules=884736 threads=1 sink=884736 median_ms=1155.5936 samples_ms=1145.1700,1149.0310,1151.0698,1152.1847,1155.5936,1156.4332,1156.6671,1159.4389
mode=symmetry molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=884736 median_ms=95.4243 samples_ms=96.4681,95.8060,91.9194,93.5268,94.1117,95.4243,92.7922,97.0989
kind=mdapy-symmetry molecules=884736 threads=8 sink=0 median_ms=10100.5363 samples_ms=10062.3637,10065.3372,10071.5438,10078.8063,10100.5363,10100.7818,10132.3576,10166.5051
kind=ovito-symmetry molecules=884736 threads=8 sink=0 median_ms=11854.5700 samples_ms=11781.5125,11803.5680,11832.3351,11845.5900,11854.5700,11858.3931,11961.7474,12075.3104
kind=mdapy-diamond molecules=884736 threads=8 sink=884736 median_ms=1264.5616 samples_ms=1249.2075,1264.2742,1264.4386,1264.4511,1264.5616,1268.6039,1268.6147,1269.5840
kind=ovito-diamond molecules=884736 threads=8 sink=884736 median_ms=1184.6236 samples_ms=1164.9931,1171.5005,1176.8059,1184.3809,1184.6236,1186.1084,1187.7768,1204.5497
mode=ortho molecules=10648 threads=1 pairs=42592 checksum=13782790659416074331 sink=0 median_ms=1.5914 samples_ms=1.6117,1.5914,1.5154,1.5930,1.5746,1.5902,1.5726,1.6246
mode=ortho molecules=10648 threads=8 pairs=42592 checksum=13782790659416074331 sink=0 median_ms=0.5004 samples_ms=0.6153,0.4965,0.5037,0.4854,0.4928,0.5013,0.4849,0.5004
mode=ortho molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=0 median_ms=12.0891 samples_ms=13.1306,11.7482,11.9661,12.0891,12.0511,12.1623,11.7855,12.3495
mode=ortho molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=0 median_ms=4.0022 samples_ms=4.6685,5.0282,3.9284,4.0329,4.0022,3.7875,3.8247,3.9604
mode=ortho molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=0 median_ms=115.9495 samples_ms=123.1739,117.0919,115.7952,116.1937,115.9495,114.4051,114.5626,115.6036
mode=ortho molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=0 median_ms=31.2759 samples_ms=35.4312,31.1635,30.8659,31.2759,32.1812,30.4012,31.5164,30.5219
```

Orthorhombic list, same binary, after that session:

```
mode=ortho molecules=10648 threads=1 pairs=42592 checksum=13782790659416074331 sink=0 median_ms=3.4731 samples_ms=3.4731,3.5137,3.3797,3.4386,3.5110,3.4071,3.5545,3.4138
mode=ortho molecules=10648 threads=8 pairs=42592 checksum=13782790659416074331 sink=0 median_ms=0.8825 samples_ms=0.9532,0.8837,0.8543,0.8825,0.8500,0.8588,0.8689,0.9265
mode=ortho molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=0 median_ms=25.4262 samples_ms=25.5751,25.2739,25.3977,25.2691,25.4661,25.1220,25.4753,25.4262
mode=ortho molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=0 median_ms=4.8818 samples_ms=5.7702,4.8028,5.0368,4.6458,4.8818,4.7980,4.8668,4.9628
mode=ortho molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=0 median_ms=299.0097 samples_ms=302.4505,299.2290,297.6589,298.1825,299.0097,298.1216,298.8550,302.7927
mode=ortho molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=0 median_ms=45.3569 samples_ms=49.2256,45.3569,46.3919,45.3136,45.0866,45.2903,45.4778,44.8664
```
