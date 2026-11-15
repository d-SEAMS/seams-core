# Pinned interleaved timings

One process per thread count. Serial is pinned to core 0. Eight threads are pinned to cores 0-7. OpenMP uses close binding and a passive wait. The Qt pool is idle before each sample, so a library does not time the other library's leftover workers.

Order is ABAB: each of 20 d-SEAMS samples is followed at once by the same sample of freud 3.6.1, OVITO 3.17.2, or mdapy 2.0.0. A win is a strictly smaller median and a strictly smaller minimum. Every one of the 66 comparisons won.

q6 is the local Steinhardt parameter and the Lechner-Dellago average. OVITO computes q2 through q6 in one modifier. mdapy computes the local value and the average as two calls. CHILL+ is the ovito definition. Freud has no CHILL+.

## 110592 molecules

| Threads | Operation | Other | d-SEAMS median | IQR | min | other median | IQR | min |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | cutoff | freud | 10.41 | 2.12 | 8.41 | 98.24 | 4.36 | 94.53 |
| 1 | cutoff | ovito | 7.05 | 0.67 | 6.65 | 23.33 | 0.40 | 22.98 |
| 1 | cutoff | mdapy | 7.21 | 0.69 | 6.79 | 41.64 | 1.08 | 40.70 |
| 1 | knn | freud | 11.37 | 0.87 | 10.38 | 116.98 | 7.05 | 111.72 |
| 1 | knn | ovito | 10.93 | 1.01 | 9.95 | 95.28 | 1.54 | 94.40 |
| 1 | knn | mdapy | 10.50 | 0.55 | 9.51 | 112.25 | 0.61 | 111.64 |
| 1 | q6 | freud | 36.91 | 1.06 | 35.91 | 1169.45 | 4.97 | 1161.58 |
| 1 | q6 | ovito | 30.51 | 3.29 | 27.47 | 79.82 | 4.49 | 77.92 |
| 1 | q6 | mdapy | 32.44 | 2.98 | 29.35 | 152.55 | 1.99 | 150.35 |
| 1 | chill_plus | ovito | 26.94 | 1.24 | 26.02 | 524.09 | 8.85 | 520.55 |
| 1 | chill_plus | mdapy | 20.51 | 1.13 | 19.59 | 58.25 | 1.14 | 57.51 |
| 8 | cutoff | freud | 4.81 | 0.70 | 4.18 | 32.92 | 0.73 | 31.60 |
| 8 | cutoff | ovito | 4.34 | 0.59 | 3.64 | 6.20 | 0.70 | 5.67 |
| 8 | cutoff | mdapy | 3.70 | 0.32 | 3.42 | 6.71 | 0.22 | 6.58 |
| 8 | knn | freud | 9.76 | 0.50 | 9.32 | 33.95 | 1.47 | 32.81 |
| 8 | knn | ovito | 9.65 | 0.55 | 9.44 | 19.90 | 0.29 | 19.51 |
| 8 | knn | mdapy | 10.16 | 0.59 | 9.55 | 16.47 | 0.64 | 15.99 |
| 8 | q6 | freud | 12.01 | 0.50 | 11.28 | 65.38 | 1.42 | 64.33 |
| 8 | q6 | ovito | 12.75 | 1.40 | 11.66 | 30.56 | 14.28 | 25.20 |
| 8 | q6 | mdapy | 12.11 | 0.88 | 11.20 | 23.29 | 1.35 | 22.35 |
| 8 | chill_plus | ovito | 8.44 | 0.61 | 6.95 | 71.81 | 0.69 | 71.42 |
| 8 | chill_plus | mdapy | 7.02 | 0.71 | 6.55 | 10.13 | 0.23 | 9.87 |

The same pattern holds at 4096 and 32768. Serial cutoff at 110592 is 7.05 ms (IQR 0.67, min 6.65) against OVITO CreateBonds at 23.33 ms (IQR 0.40, min 22.98) on the same core. The raw lines are benchmarks/h2h-robust-1.txt and benchmarks/h2h-robust-8.txt.
