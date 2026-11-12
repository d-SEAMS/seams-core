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
is 25.43 ms and the triclinic list is 30.57 ms. At 8 threads the pair
is 4.88 ms and 5.95 ms.

At 884736 molecules the triclinic list is 443 ms against 299 ms
orthorhombic. The same stencil is in use. The working set of the three
coordinate streams no longer fits in L2, and the sheared coordinates
miss more often.

## Triclinic neighbour list

Every code returned `4N` pairs. The index checksum matches across thread
counts and matches the orthorhombic list: this shear does not change
which molecules are neighbours.

| molecules | threads | d-SEAMS | vesin | freud | orthorhombic |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10648 | 1 | 3.485 | 5.964 | 10.373 | 3.473 |
| 10648 | 8 | 1.005 | 4.601 | 10.434 | 0.883 |
| 110592 | 1 | 30.572 | 51.994 | 122.628 | 25.426 |
| 110592 | 8 | 6.093 | 24.911 | 127.350 | 4.882 |
| 884736 | 1 | 443.024 | 1141.327 | 1191.634 | 299.010 |
| 884736 | 8 | 65.550 | 203.597 | 1153.981 | 45.357 |

The orthorhombic column is the same binary, run after the session below.
The 110592 rows were measured back to back with a repeat of the
triclinic list (30.57 ms and 5.95 ms).

## CHILL+

Each d-SEAMS sample builds the CSR and then runs the analysis on those
rows. CHILL+ keeps the four nearest neighbours of each row and
accumulates the real spherical harmonics at l = 3, the same contraction
as the fused kernel. Every molecule of this ice is cubic. freud 3.6.1
has no CHILL+.

| molecules | threads | d-SEAMS | mdapy | OVITO |
| ---: | ---: | ---: | ---: | ---: |
| 110592 | 1 | 40.862 | 65.618 | 628.011 |
| 110592 | 8 | 9.070 | 88.058 | 92.206 |
| 884736 | 1 | 508.358 | 758.383 | 4942.828 |
| 884736 | 8 | 92.831 | 794.323 | 659.641 |

## Diamond common-neighbour test

The symmetry column is the diamond common-neighbour test on the 3.5 Å
CSR. A cubic shell is the twelve outward bonds of the four nearest
neighbours. It is not the template overlay. PTM has no cubic-ice
template: on this crystal every structure code is 0, in mdapy and in
OVITO, which is why those sinks are 0. Identify diamond and mdapy
`DiamondStructure` label every molecule, and so does this kernel.

| molecules | threads | d-SEAMS | mdapy diamond | OVITO identify diamond | mdapy PTM | OVITO PTM |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 110592 | 1 | 79.765 | 147.656 | 148.962 | 1341.099 | 1548.133 |
| 110592 | 8 | 19.395 | 188.508 | 37.889 | 1344.673 | 208.699 |
| 884736 | 1 | 720.704 | 1190.894 | 1186.497 | 9953.016 | 11783.641 |
| 884736 | 8 | 100.959 | 1295.796 | 290.571 | 10196.825 | 1620.720 |

## Rings

The census is the short-hop search: for rings of size at most 7 the
adjacency answers the primitivity test, and the ball table is not
built. The CSR is handed to that search. The counts match mdapy:
221184 rings at 110592 molecules and 1769472 at 884736. The OVITO ring
sink is the particle count. `FindRings` does not write a structure
code. freud has no ring census.

| molecules | threads | d-SEAMS | mdapy | OVITO |
| ---: | ---: | ---: | ---: | ---: |
| 110592 | 1 | 86.970 | 2330.859 | 560.793 |
| 110592 | 8 | 32.497 | 2342.910 | 282.745 |
| 884736 | 1 | 831.572 | 17972.932 | 4855.287 |
| 884736 | 8 | 273.613 | 18066.436 | 2601.361 |

## How to rerun

From the repository root, with the release build and the comparison
environment on `PATH` via `NEIGH_PYTHON`:

```
NEIGH_PYTHON=/tmp/neigh-venv/bin/python \
NEIGH_BENCH=build-release/tests/bench_v2 \
python scripts/bench_neighbours_v2.py
```

The orthorhombic rows are `bench_v2 ortho <molecules> <threads> 8`
under the same pin.

## Raw log

```
mode=tilt molecules=10648 threads=1 pairs=42592 checksum=13782790659416074331 sink=0 median_ms=3.4845 samples_ms=3.4845,3.5386,3.4384,3.4343,3.5647,3.3983,3.4383,3.5188
kind=tilt-vesin molecules=10648 threads=1 sink=42592 median_ms=5.9642 samples_ms=5.8079,5.8388,5.8820,5.9372,5.9642,5.9650,6.0336,6.0468
kind=tilt-freud molecules=10648 threads=1 sink=42592 median_ms=10.3732 samples_ms=10.3084,10.3655,10.3669,10.3693,10.3732,10.4011,10.4332,10.4731
mode=tilt molecules=10648 threads=8 pairs=42592 checksum=13782790659416074331 sink=0 median_ms=1.0046 samples_ms=1.0062,0.9460,0.9495,1.0046,0.9722,1.0144,1.0139,0.9479
kind=tilt-vesin molecules=10648 threads=8 sink=42592 median_ms=4.6014 samples_ms=3.0579,3.1368,3.5329,3.6727,4.6014,4.6093,4.6589,4.7537
kind=tilt-freud molecules=10648 threads=8 sink=42592 median_ms=10.4342 samples_ms=10.3497,10.4138,10.4159,10.4184,10.4342,10.4482,10.5116,10.5148
mode=tilt molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=0 median_ms=30.5720 samples_ms=30.8342,30.3009,30.5775,30.5255,30.5955,30.5720,30.3354,30.4279
kind=tilt-vesin molecules=110592 threads=1 sink=442368 median_ms=51.9942 samples_ms=51.6919,51.6956,51.7981,51.8662,51.9942,52.0650,52.0916,52.4943
kind=tilt-freud molecules=110592 threads=1 sink=442368 median_ms=122.6279 samples_ms=120.5045,121.3189,121.3590,121.8853,122.6279,122.7369,122.7974,127.0188
mode=tilt molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=0 median_ms=6.0934 samples_ms=6.9065,6.2982,6.1149,5.9821,6.0174,6.0934,5.9758,6.0216
kind=tilt-vesin molecules=110592 threads=8 sink=442368 median_ms=24.9108 samples_ms=24.6145,24.7194,24.7664,24.8324,24.9108,24.9483,25.1186,26.0167
kind=tilt-freud molecules=110592 threads=8 sink=442368 median_ms=127.3502 samples_ms=126.4897,127.0935,127.1074,127.1738,127.3502,127.9802,128.4043,137.9679
mode=tilt molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=0 median_ms=443.0236 samples_ms=446.1631,449.7723,439.4564,443.0236,438.4797,453.6831,437.5949,438.9384
kind=tilt-vesin molecules=884736 threads=1 sink=3538944 median_ms=1141.3268 samples_ms=1118.4384,1124.0209,1129.3905,1134.7615,1141.3268,1141.7548,1147.7133,1191.9243
kind=tilt-freud molecules=884736 threads=1 sink=3538944 median_ms=1191.6338 samples_ms=1156.5232,1159.8109,1183.4468,1190.8607,1191.6338,1193.9519,1202.9605,1222.5416
mode=tilt molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=0 median_ms=65.5497 samples_ms=70.4154,65.5063,66.1948,65.5497,64.3381,64.5287,65.4321,85.5898
kind=tilt-vesin molecules=884736 threads=8 sink=3538944 median_ms=203.5968 samples_ms=181.5102,185.9666,201.5857,203.2913,203.5968,219.4238,221.3960,221.9107
kind=tilt-freud molecules=884736 threads=8 sink=3538944 median_ms=1153.9814 samples_ms=1137.9867,1142.0892,1142.1750,1148.5169,1153.9814,1157.1931,1157.7741,1165.3098
mode=chill molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=110592 median_ms=40.8617 samples_ms=47.7965,44.4743,41.1771,40.7971,40.8617,40.1922,39.9350,39.9228
kind=mdapy-chill molecules=110592 threads=1 sink=110592 median_ms=65.6183 samples_ms=62.2688,62.3768,62.8975,63.2062,65.6183,66.0663,66.7361,68.0841
kind=ovito-chill molecules=110592 threads=1 sink=110592 median_ms=628.0114 samples_ms=625.8838,626.1325,626.5083,627.8998,628.0114,632.9294,653.2898,694.5793
kind=freud-none molecules=110592 threads=1 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=chill molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=110592 median_ms=9.0697 samples_ms=12.2989,10.0419,9.6273,8.8737,9.0697,8.9163,8.8115,8.8248
kind=mdapy-chill molecules=110592 threads=8 sink=110592 median_ms=88.0583 samples_ms=85.3887,85.6530,86.9503,87.2410,88.0583,88.2143,96.7992,100.5640
kind=ovito-chill molecules=110592 threads=8 sink=110592 median_ms=92.2060 samples_ms=89.5883,89.8483,91.9378,91.9847,92.2060,93.0686,96.6500,98.1130
kind=freud-none molecules=110592 threads=8 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=chill molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=884736 median_ms=508.3578 samples_ms=484.4607,516.6233,525.1612,497.2228,490.2562,508.3578,528.8858,488.2780
kind=mdapy-chill molecules=884736 threads=1 sink=884736 median_ms=758.3834 samples_ms=722.9076,731.4055,735.8126,748.7469,758.3834,764.2421,776.5018,816.5757
kind=ovito-chill molecules=884736 threads=1 sink=884736 median_ms=4942.8283 samples_ms=4932.6113,4935.4970,4937.0440,4938.0039,4942.8283,4992.5839,5046.1470,5447.1433
kind=freud-none molecules=884736 threads=1 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=chill molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=884736 median_ms=92.8312 samples_ms=99.8517,92.8312,91.6297,95.5418,91.6615,94.1322,91.5696,92.1825
kind=mdapy-chill molecules=884736 threads=8 sink=884736 median_ms=794.3232 samples_ms=737.9534,748.0695,763.8569,773.6553,794.3232,795.0764,823.1204,1002.9260
kind=ovito-chill molecules=884736 threads=8 sink=884736 median_ms=659.6414 samples_ms=649.9426,651.5619,652.0473,653.7416,659.6414,661.1631,661.5111,679.6672
kind=freud-none molecules=884736 threads=8 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=rings molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=221184 median_ms=86.9703 samples_ms=87.5309,86.9088,86.5027,86.9703,86.4979,86.8492,87.3591,87.1577
kind=mdapy-rings molecules=110592 threads=1 sink=221184 median_ms=2330.8591 samples_ms=2318.7842,2319.6780,2323.2265,2324.9518,2330.8591,2338.2140,2351.8763,2358.0680
kind=ovito-rings molecules=110592 threads=1 sink=110592 median_ms=560.7932 samples_ms=536.0512,537.1700,540.6751,553.8979,560.7932,564.6724,575.9788,614.6378
mode=rings molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=221184 median_ms=32.4965 samples_ms=33.3868,32.4965,32.5108,31.3206,35.8687,29.8305,29.7614,28.7368
kind=mdapy-rings molecules=110592 threads=8 sink=221184 median_ms=2342.9103 samples_ms=2331.9907,2333.2140,2335.7342,2337.7846,2342.9103,2343.5751,2345.8716,2349.1580
kind=ovito-rings molecules=110592 threads=8 sink=110592 median_ms=282.7452 samples_ms=257.5110,264.3925,273.1819,280.8424,282.7452,329.3495,329.6736,330.2904
mode=rings molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=1769472 median_ms=831.5721 samples_ms=827.3854,824.1449,849.1644,832.0921,841.6491,830.8284,831.5721,825.6962
kind=mdapy-rings molecules=884736 threads=1 sink=1769472 median_ms=17972.9315 samples_ms=17783.4366,17840.4640,17842.1896,17884.9131,17972.9315,18078.4878,18138.7402,18359.2966
kind=ovito-rings molecules=884736 threads=1 sink=884736 median_ms=4855.2865 samples_ms=4740.1998,4804.8986,4810.2427,4834.6359,4855.2865,4884.2288,4917.7673,4943.5882
mode=rings molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=1769472 median_ms=273.6129 samples_ms=274.8406,276.6449,269.1452,273.6886,267.8904,263.6135,273.6129,267.8700
kind=mdapy-rings molecules=884736 threads=8 sink=1769472 median_ms=18066.4357 samples_ms=17981.4100,18032.5574,18035.3208,18038.2780,18066.4357,18341.2014,18440.5632,18662.7000
kind=ovito-rings molecules=884736 threads=8 sink=884736 median_ms=2601.3606 samples_ms=2518.6888,2549.2956,2564.4558,2591.2438,2601.3606,2635.2182,2773.4944,2785.3184
mode=symmetry molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=110592 median_ms=79.7653 samples_ms=77.5893,77.2043,77.0657,81.1700,80.1033,79.7653,80.3813,78.8864
kind=mdapy-symmetry molecules=110592 threads=1 sink=0 median_ms=1341.0993 samples_ms=1323.0447,1324.2901,1329.9375,1330.5075,1341.0993,1366.7108,1373.8571,1389.8981
kind=ovito-symmetry molecules=110592 threads=1 sink=0 median_ms=1548.1332 samples_ms=1521.5844,1531.7906,1533.7089,1547.3274,1548.1332,1553.4330,1553.8319,1575.5963
kind=mdapy-diamond molecules=110592 threads=1 sink=110592 median_ms=147.6562 samples_ms=145.4016,145.4931,146.4816,147.4334,147.6562,147.8950,149.5743,150.0196
kind=ovito-diamond molecules=110592 threads=1 sink=110592 median_ms=148.9624 samples_ms=146.7945,146.9238,147.9757,148.2574,148.9624,149.4918,152.8464,154.1139
mode=symmetry molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=110592 median_ms=19.3948 samples_ms=19.4088,19.3604,19.3490,19.3948,18.8877,18.9297,19.4076,19.6276
kind=mdapy-symmetry molecules=110592 threads=8 sink=0 median_ms=1344.6730 samples_ms=1336.2287,1338.0815,1339.0507,1342.4747,1344.6730,1348.4066,1350.0285,1354.9662
kind=ovito-symmetry molecules=110592 threads=8 sink=0 median_ms=208.6991 samples_ms=206.5527,206.9600,207.4455,207.6940,208.6991,209.2907,210.8938,211.8429
kind=mdapy-diamond molecules=110592 threads=8 sink=110592 median_ms=188.5075 samples_ms=187.1940,187.5423,187.9940,188.4358,188.5075,188.6299,188.9083,191.8159
kind=ovito-diamond molecules=110592 threads=8 sink=110592 median_ms=37.8889 samples_ms=30.4845,30.9932,31.9385,32.0141,37.8889,43.2822,44.7059,45.1765
mode=symmetry molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=884736 median_ms=720.7038 samples_ms=720.7038,724.1510,720.8560,715.9333,734.8162,716.4234,720.4348,720.3075
kind=mdapy-symmetry molecules=884736 threads=1 sink=0 median_ms=9953.0155 samples_ms=9919.1566,9926.2854,9928.6786,9928.8757,9953.0155,9954.5970,9961.8055,9987.2578
kind=ovito-symmetry molecules=884736 threads=1 sink=0 median_ms=11783.6407 samples_ms=11707.3837,11711.7610,11715.0135,11737.1384,11783.6407,11803.6552,11909.0205,11917.8488
kind=mdapy-diamond molecules=884736 threads=1 sink=884736 median_ms=1190.8939 samples_ms=1182.7263,1187.0285,1188.4473,1189.6744,1190.8939,1192.0739,1192.8388,1194.7090
kind=ovito-diamond molecules=884736 threads=1 sink=884736 median_ms=1186.4967 samples_ms=1178.6761,1180.9197,1183.4787,1185.7268,1186.4967,1189.3460,1199.3876,1203.6319
mode=symmetry molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=884736 median_ms=100.9590 samples_ms=101.6876,97.3104,96.8658,101.2449,98.5655,103.9143,97.7471,100.9590
kind=mdapy-symmetry molecules=884736 threads=8 sink=0 median_ms=10196.8250 samples_ms=10123.0603,10135.7217,10179.8746,10193.8160,10196.8250,10204.9248,10234.5713,10556.0501
kind=ovito-symmetry molecules=884736 threads=8 sink=0 median_ms=1620.7203 samples_ms=1586.5863,1599.6344,1604.5385,1619.9022,1620.7203,1625.6959,1632.5223,1634.3473
kind=mdapy-diamond molecules=884736 threads=8 sink=884736 median_ms=1295.7963 samples_ms=1264.7065,1289.5569,1292.5730,1295.6444,1295.7963,1298.8485,1304.4819,1315.7347
kind=ovito-diamond molecules=884736 threads=8 sink=884736 median_ms=290.5705 samples_ms=270.2464,278.3434,282.0129,285.5709,290.5705,330.5988,342.5576,343.6381
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
