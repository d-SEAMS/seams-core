# Scale-up, 442368 molecules

Tilted ice Ih, four times the 110592-molecule cell (24 by 24 by 96). Cutoff neighbours, both CHILL+ modes, and primitive rings. Eight samples. Serial is pinned to one core. Eight threads are pinned to cores 0-7. The machine has 8 cores. No host name is recorded. The raw log is benchmarks/h2h-scale.txt.

| Threads | Operation | Other | d-SEAMS median | IQR | min | other median | IQR | min |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | cutoff | freud | 55.96 | 1.90 | 55.57 | 509.34 | 26.52 | 496.38 |
| 1 | cutoff | ovito | 42.45 | 2.83 | 41.62 | 130.94 | 2.88 | 127.21 |
| 1 | cutoff | mdapy | 49.50 | 3.69 | 44.76 | 223.53 | 20.62 | 208.11 |
| 1 | chill_plus | ovito | 130.66 | 2.10 | 104.55 | 2375.70 | 41.01 | 2340.01 |
| 1 | chill_plus | mdapy | 133.78 | 2.83 | 131.78 | 291.33 | 8.83 | 282.16 |
| 1 | rings | ovito | 502.90 | 15.03 | 484.38 | 2260.12 | 37.71 | 2233.38 |
| 1 | rings | rings-statistics-matter | 495.81 | 18.16 | 484.52 | 6729.69 | 38.39 | 6684.91 |
| 8 | cutoff | freud | 26.71 | 1.76 | 25.44 | 179.96 | 9.59 | 172.67 |
| 8 | cutoff | ovito | 27.52 | 3.51 | 25.06 | 142.30 | 7.93 | 136.30 |
| 8 | cutoff | mdapy | 26.08 | 4.52 | 21.02 | 31.03 | 3.39 | 29.11 |
| 8 | chill_plus | ovito | 44.10 | 3.66 | 39.89 | 2461.74 | 55.49 | 2381.22 |
| 8 | chill_plus | mdapy | 34.82 | 12.41 | 32.30 | 38.58 | 8.26 | 37.71 |
| 8 | rings | ovito | 201.99 | 18.45 | 174.48 | 2378.19 | 81.14 | 2330.96 |
| 8 | rings | rings-statistics-matter | 187.71 | 6.97 | 176.30 | 1556.91 | 30.50 | 1521.39 |

Paper-mode CHILL+ is timed on its own: 107.59 ms serial (IQR 0.68, minimum 106.40) and 44.52 ms at eight threads (IQR 2.39, minimum 41.44). Every comparison wins on the median and on the minimum.
