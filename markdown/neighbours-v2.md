# Cutoff list, triclinic cells, and the analyses that use it

Ice Ic, lattice constant 6.35 Å. The triclinic cell keeps that cube and
adds xy = 0.25 L, xz = 0.05 L, yz = 0.12 L. A pair is kept when the
squared minimum image is strictly below the cutoff squared. CHILL+ and
the rings use 3.5 Å. The template overlay uses 6.0 Å so a 12-neighbour
shell exists. Eight timed calls after one warmup, pinned with `taskset`
to cores `0..T-1`. The median is the 5th of the 8 sorted samples.

The same machine as `markdown/neighbours.md`: kernel 6.12.113+, KVM
guest, Intel Xeon family 6 model 207, 8 cores, GCC 15.2.0, release
`-O3`. vesin 0.6.2, freud 3.6.1, mdapy 2.0.0, OVITO 3.17.2.

## Triclinic neighbour list

Every code returned `4N` pairs. The index checksum is the same at 1 and
8 threads, and the same as the orthorhombic list: this shear does not
change which molecules are neighbours.

| molecules | threads | d-SEAMS | vesin | freud |
| ---: | ---: | ---: | ---: | ---: |
| 10648 | 1 | 5.469 | 6.140 | 10.570 |
| 10648 | 8 | 1.409 | 3.146 | 10.544 |
| 110592 | 1 | 43.250 | 52.412 | 138.003 |
| 110592 | 8 | 10.873 | 22.582 | 137.605 |
| 884736 | 1 | 643.776 | 1161.736 | 1291.377 |
| 884736 | 8 | 97.572 | 194.677 | 1228.686 |

## End to end

Each d-SEAMS sample builds the CSR and then runs the analysis on those
rows. CHILL+ labels every molecule of this ice cubic (110592 and
884736). The ring census is 221184 rings at 110592 molecules and
1769472 at 884736, the same counts mdapy reports. freud 3.6.1 has no
CHILL+, no ring census, and no polyhedral template match.

OVITO caches a pipeline. The first session's OVITO samples were that
cache (a few microseconds) and are not the analysis. The OVITO column
is a second pinned run, same binary and same inputs, with a new
pipeline on every sample.

| analysis | molecules | threads | d-SEAMS | mdapy | OVITO |
| --- | ---: | ---: | ---: | ---: | ---: |
| CHILL+ | 110592 | 1 | 66.261 | 69.247 | 623.611 |
| CHILL+ | 110592 | 8 | 16.797 | 87.292 | 89.249 |
| CHILL+ | 884736 | 1 | 645.498 | 836.341 | 4896.235 |
| CHILL+ | 884736 | 8 | 171.240 | 976.795 | 671.322 |
| rings | 110592 | 1 | 127.569 | 2444.762 | 530.515 |
| rings | 110592 | 8 | 45.258 | 2408.506 | 263.796 |
| rings | 884736 | 1 | 1143.636 | 18713.615 | 5447.092 |
| rings | 884736 | 8 | 387.861 | 18813.651 | 2974.254 |
| symmetry | 110592 | 1 | 5271.913 | 1320.978 | 1492.315 |
| symmetry | 110592 | 8 | 752.255 | 1383.292 | 215.615 |
| symmetry | 884736 | 1 | 42315.438 | 10203.100 | 11799.265 |
| symmetry | 884736 | 8 | 5489.950 | 10255.042 | 1663.210 |

The template overlay is the d-SEAMS symmetry classifier. mdapy and
OVITO time polyhedral template matching. Those are different kernels;
the table is the end-to-end cost of each code's own classifier.

## Rerun

```
meson setup build-release --buildtype=release -Dwith_tests=true \
  -Dwith_cli=true -Dwith_mpi=disabled -Dwith_openmp_offload=disabled
ninja -C build-release tests/bench_v2
python3 -m venv /tmp/neigh-venv
/tmp/neigh-venv/bin/pip install 'vesin==0.6.2' freud-analysis mdapy ovito
# OVITO needs libGL on the loader path inside that process only.
NEIGH_PYTHON=/tmp/neigh-venv/bin/python \
  NEIGH_BENCH=$PWD/build-release/tests/bench_v2 \
  /tmp/neigh-venv/bin/python scripts/bench_neighbours_v2.py
```

`bench_v2` arguments are `tilt|chill|rings|symmetry`, molecules, threads,
repeats.

## Raw log

First session. OVITO lines here are pipeline-cache hits, not the modifier.

```
mode=tilt molecules=10648 threads=1 pairs=42592 checksum=13782790659416074331 sink=0 median_ms=5.4694 samples_ms=5.8035,5.8004,5.4693,5.4672,5.3196,5.4694,5.4569,5.5500
kind=tilt-vesin molecules=10648 threads=1 sink=42592 median_ms=6.1398 samples_ms=5.8734,5.9420,6.0128,6.1336,6.1398,6.1417,6.1440,8.2617
kind=tilt-freud molecules=10648 threads=1 sink=42592 median_ms=10.5696 samples_ms=10.4876,10.5216,10.5335,10.5593,10.5696,10.6212,10.6959,10.8410
mode=tilt molecules=10648 threads=8 pairs=42592 checksum=13782790659416074331 sink=0 median_ms=1.4086 samples_ms=1.5940,1.2814,1.3005,1.4171,1.3959,1.5946,1.4086,1.3566
kind=tilt-vesin molecules=10648 threads=8 sink=42592 median_ms=3.1458 samples_ms=3.0275,3.0735,3.0789,3.0918,3.1458,3.1564,3.1837,3.3098
kind=tilt-freud molecules=10648 threads=8 sink=42592 median_ms=10.5442 samples_ms=10.4068,10.4725,10.4864,10.4943,10.5442,10.5640,10.5913,10.7294
mode=tilt molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=0 median_ms=43.2496 samples_ms=44.1532,42.9910,43.2191,43.2496,42.9686,43.1869,43.5140,43.7562
kind=tilt-vesin molecules=110592 threads=1 sink=442368 median_ms=52.4123 samples_ms=52.1721,52.2641,52.3041,52.3959,52.4123,52.4749,52.5368,53.5748
kind=tilt-freud molecules=110592 threads=1 sink=442368 median_ms=138.0031 samples_ms=131.0809,134.6335,135.2706,137.6814,138.0031,140.1396,140.5574,143.0089
mode=tilt molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=0 median_ms=10.8734 samples_ms=11.3421,10.8734,10.5368,10.7526,11.5577,10.8917,10.7687,10.7668
kind=tilt-vesin molecules=110592 threads=8 sink=442368 median_ms=22.5815 samples_ms=22.2321,22.3464,22.3494,22.5538,22.5815,22.7213,22.7905,27.6723
kind=tilt-freud molecules=110592 threads=8 sink=442368 median_ms=137.6046 samples_ms=129.6424,134.0783,134.1220,135.8662,137.6046,138.2686,142.4777,147.7156
mode=tilt molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=0 median_ms=643.7760 samples_ms=678.4454,656.2402,643.7760,642.0629,634.8649,633.8881,704.8107,643.1859
kind=tilt-vesin molecules=884736 threads=1 sink=3538944 median_ms=1161.7356 samples_ms=1128.4045,1132.3492,1142.3155,1143.6869,1161.7356,1200.0479,1211.3307,1286.8518
kind=tilt-freud molecules=884736 threads=1 sink=3538944 median_ms=1291.3770 samples_ms=1196.9410,1242.5377,1259.7852,1274.0986,1291.3770,1298.8091,1306.9458,1670.1540
mode=tilt molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=0 median_ms=97.5721 samples_ms=111.6611,96.1457,93.9973,102.0276,108.7273,94.5705,97.5721,93.0815
kind=tilt-vesin molecules=884736 threads=8 sink=3538944 median_ms=194.6769 samples_ms=191.4052,191.5190,191.7221,191.8472,194.6769,194.8288,196.0686,211.0110
kind=tilt-freud molecules=884736 threads=8 sink=3538944 median_ms=1228.6862 samples_ms=1163.8773,1171.1296,1174.2214,1177.6055,1228.6862,1228.8607,1234.3629,1305.5577
mode=chill molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=995328 median_ms=66.2609 samples_ms=66.2609,68.2940,67.1449,65.7767,65.5065,64.7987,69.1611,64.6082
kind=mdapy-chill molecules=110592 threads=1 sink=110592 median_ms=69.2474 samples_ms=66.3048,67.2739,67.5456,68.5590,69.2474,69.6208,70.7197,84.0183
kind=ovito-chill molecules=110592 threads=1 sink=110592 median_ms=0.0050 samples_ms=0.0043,0.0044,0.0044,0.0049,0.0050,0.0062,0.0067,0.0112
kind=freud-none molecules=110592 threads=1 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=chill molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=995328 median_ms=16.7972 samples_ms=21.9749,19.4653,17.0059,16.3259,16.7458,16.6689,16.7972,15.6978
kind=mdapy-chill molecules=110592 threads=8 sink=110592 median_ms=87.2924 samples_ms=86.3223,86.4517,86.6365,87.1728,87.2924,87.8812,87.8878,88.5378
kind=ovito-chill molecules=110592 threads=8 sink=110592 median_ms=0.0047 samples_ms=0.0038,0.0039,0.0039,0.0040,0.0047,0.0054,0.0071,0.0089
kind=freud-none molecules=110592 threads=8 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=chill molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=7962624 median_ms=645.4979 samples_ms=637.8857,651.5962,660.0987,630.3249,623.7931,621.9199,652.2614,645.4979
kind=mdapy-chill molecules=884736 threads=1 sink=884736 median_ms=836.3412 samples_ms=800.9123,810.1987,812.8663,834.5773,836.3412,838.6643,847.9760,899.2196
kind=ovito-chill molecules=884736 threads=1 sink=884736 median_ms=0.0074 samples_ms=0.0066,0.0068,0.0070,0.0073,0.0074,0.0086,0.0099,0.0182
kind=freud-none molecules=884736 threads=1 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=chill molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=7962624 median_ms=171.2404 samples_ms=169.9806,209.9289,205.4207,181.1393,168.2621,171.2404,168.4129,165.9268
kind=mdapy-chill molecules=884736 threads=8 sink=884736 median_ms=976.7947 samples_ms=881.5986,941.2456,966.9429,973.4144,976.7947,985.0446,1013.1856,1057.7554
kind=ovito-chill molecules=884736 threads=8 sink=884736 median_ms=0.0070 samples_ms=0.0047,0.0051,0.0053,0.0061,0.0070,0.0071,0.0073,0.0099
kind=freud-none molecules=884736 threads=8 sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm
mode=rings molecules=110592 threads=1 pairs=442368 checksum=8222705064348099955 sink=1990656 median_ms=127.5689 samples_ms=130.9310,130.2538,128.3202,126.4972,125.1778,124.3259,127.3946,127.5689
kind=mdapy-rings molecules=110592 threads=1 sink=221184 median_ms=2444.7624 samples_ms=2350.0490,2393.7827,2433.4134,2437.0296,2444.7624,2454.6790,2464.5144,2501.5082
kind=ovito-rings molecules=110592 threads=1 sink=110592 median_ms=0.0064 samples_ms=0.0054,0.0054,0.0054,0.0058,0.0064,0.0065,0.0099,0.0116
mode=rings molecules=110592 threads=8 pairs=442368 checksum=8222705064348099955 sink=1990656 median_ms=45.2579 samples_ms=47.3595,47.3958,44.9227,45.6678,45.2579,44.0494,42.9542,42.4028
kind=mdapy-rings molecules=110592 threads=8 sink=221184 median_ms=2408.5059 samples_ms=2363.5271,2378.8398,2402.6982,2404.2948,2408.5059,2423.6773,2445.5161,2550.8374
kind=ovito-rings molecules=110592 threads=8 sink=110592 median_ms=0.0062 samples_ms=0.0042,0.0046,0.0051,0.0057,0.0062,0.0068,0.0080,0.0106
mode=rings molecules=884736 threads=1 pairs=3538944 checksum=9316061812838145587 sink=15925248 median_ms=1143.6355 samples_ms=1149.7267,1143.6355,1135.4000,1132.2211,1151.7669,1138.1171,1140.4480,1181.7630
kind=mdapy-rings molecules=884736 threads=1 sink=1769472 median_ms=18713.6149 samples_ms=18232.3542,18319.5140,18589.6939,18632.7612,18713.6149,18916.7269,19242.5190,20350.7528
kind=ovito-rings molecules=884736 threads=1 sink=884736 median_ms=0.0108 samples_ms=0.0060,0.0060,0.0064,0.0066,0.0108,0.0173,0.0184,0.0190
mode=rings molecules=884736 threads=8 pairs=3538944 checksum=9316061812838145587 sink=15925248 median_ms=387.8613 samples_ms=387.8613,385.5491,404.1719,410.7128,366.0673,406.4983,372.9953,382.5144
kind=mdapy-rings molecules=884736 threads=8 sink=1769472 median_ms=18813.6509 samples_ms=18337.4482,18368.5645,18576.4948,18756.5549,18813.6509,18826.3282,18866.1382,19134.2245
kind=ovito-rings molecules=884736 threads=8 sink=884736 median_ms=0.0063 samples_ms=0.0039,0.0044,0.0049,0.0057,0.0063,0.0074,0.0097,0.0139
mode=symmetry molecules=110592 threads=1 pairs=3096576 checksum=7146905601224204987 sink=211941 median_ms=5271.9131 samples_ms=5270.1061,5294.6938,5271.9131,5336.1892,5283.8265,5247.8390,5212.0431,5210.0729
kind=mdapy-symmetry molecules=110592 threads=1 sink=110592 median_ms=1320.9778 samples_ms=1285.2821,1303.6750,1304.2746,1304.3065,1320.9778,1346.8384,1352.2243,1372.0633
kind=ovito-symmetry molecules=110592 threads=1 sink=110592 median_ms=0.0072 samples_ms=0.0045,0.0059,0.0059,0.0067,0.0072,0.0074,0.0086,0.0138
mode=symmetry molecules=110592 threads=8 pairs=3096576 checksum=7146905601224204987 sink=211941 median_ms=752.2546 samples_ms=791.4657,730.3517,772.7716,752.2546,744.7039,741.3751,737.8743,764.3317
kind=mdapy-symmetry molecules=110592 threads=8 sink=110592 median_ms=1383.2924 samples_ms=1353.5930,1354.8185,1360.7211,1364.5653,1383.2924,1434.9549,1446.4086,1495.6604
kind=ovito-symmetry molecules=110592 threads=8 sink=110592 median_ms=0.0066 samples_ms=0.0050,0.0052,0.0055,0.0061,0.0066,0.0072,0.0081,0.0103
mode=symmetry molecules=884736 threads=1 pairs=24772608 checksum=12092942904217230779 sink=1565235 median_ms=42315.4380 samples_ms=43369.9231,42582.9106,42115.1657,42315.4380,44814.5362,42095.6510,41736.4238,42075.9051
kind=mdapy-symmetry molecules=884736 threads=1 sink=884736 median_ms=10203.0997 samples_ms=10053.7081,10085.6693,10130.3299,10178.4926,10203.0997,10223.5510,10234.0972,10824.0909
kind=ovito-symmetry molecules=884736 threads=1 sink=884736 median_ms=0.0095 samples_ms=0.0051,0.0059,0.0062,0.0083,0.0095,0.0123,0.0129,0.0243
mode=symmetry molecules=884736 threads=8 pairs=24772608 checksum=12092942904217230779 sink=1565235 median_ms=5489.9496 samples_ms=5794.4638,5372.6124,5396.9478,5489.9496,5580.1773,5567.7455,5473.1703,5441.1910
kind=mdapy-symmetry molecules=884736 threads=8 sink=884736 median_ms=10255.0416 samples_ms=10221.3171,10231.6093,10246.5506,10247.3593,10255.0416,10307.7585,10347.1446,10357.6808
kind=ovito-symmetry molecules=884736 threads=8 sink=884736 median_ms=0.0082 samples_ms=0.0055,0.0060,0.0062,0.0070,0.0082,0.0083,0.0133,0.0149
```

OVITO, new pipeline on every sample.

```
kind=ovito-chill molecules=110592 threads=1 sink=110592 median_ms=623.6111 samples_ms=617.0180,618.1420,618.2556,618.6737,623.6111,623.6662,627.1549,632.3069
kind=ovito-chill molecules=110592 threads=8 sink=110592 median_ms=89.2492 samples_ms=85.1719,85.7006,86.6478,88.6278,89.2492,91.4434,91.6826,93.2933
kind=ovito-chill molecules=884736 threads=1 sink=884736 median_ms=4896.2353 samples_ms=4884.1540,4885.9815,4888.6953,4890.4226,4896.2353,4899.0701,4899.3618,4900.8383
kind=ovito-chill molecules=884736 threads=8 sink=884736 median_ms=671.3220 samples_ms=642.4324,649.6557,651.7378,670.2815,671.3220,671.9047,678.8567,688.9371
kind=ovito-rings molecules=110592 threads=1 sink=110592 median_ms=530.5154 samples_ms=512.6235,513.2112,514.3848,525.6075,530.5154,531.1163,538.2856,544.0843
kind=ovito-rings molecules=110592 threads=8 sink=110592 median_ms=263.7959 samples_ms=253.5340,256.0265,257.9678,262.2769,263.7959,290.6575,313.1720,332.0823
kind=ovito-rings molecules=884736 threads=1 sink=884736 median_ms=5447.0917 samples_ms=4873.5515,5036.1850,5149.0061,5289.5472,5447.0917,5535.2471,5895.1242,6011.9412
kind=ovito-rings molecules=884736 threads=8 sink=884736 median_ms=2974.2537 samples_ms=2689.5709,2883.6877,2944.4229,2961.4580,2974.2537,3110.4597,3189.1869,3385.7535
kind=ovito-symmetry molecules=110592 threads=1 sink=110592 median_ms=1492.3153 samples_ms=1482.4881,1485.5177,1486.5317,1487.3412,1492.3153,1513.0945,1516.5200,1537.4541
kind=ovito-symmetry molecules=110592 threads=8 sink=110592 median_ms=215.6151 samples_ms=209.4045,210.6129,211.7443,212.4900,215.6151,216.4951,219.6904,227.5904
kind=ovito-symmetry molecules=884736 threads=1 sink=884736 median_ms=11799.2654 samples_ms=11728.4484,11730.6762,11753.6352,11758.3147,11799.2654,11893.3023,11894.5733,12167.3642
kind=ovito-symmetry molecules=884736 threads=8 sink=884736 median_ms=1663.2098 samples_ms=1616.6216,1635.8872,1650.1107,1653.7339,1663.2098,1686.1690,1728.3562,1800.0034
```
