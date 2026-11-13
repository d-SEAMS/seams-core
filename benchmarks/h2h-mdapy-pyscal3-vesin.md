# d-SEAMS, mdapy, pyscal3, and vesin

One process, eight samples, eight cores. Orthorhombic ice Ih, oxygen molecules only. Medians are milliseconds. The raw lines are `benchmarks/h2h-mdapy-pyscal3-vesin.txt`.

Unordered pairs inside 3.5 Å match exactly (8192, 65536, and 221184). Each oxygen's four nearest neighbours match exactly, and on this crystal they are that cutoff shell. q6 and q6bar agree within 1e-8 (the largest absolute difference is 2.4e-15). Published CHILL, published CHILL+, and mdapy CHILL+ label every oxygen hexagonal. pyscal3 has no CHILL; its extended-CNA diamond count is hexagonal diamond for every oxygen. vesin builds a cutoff list and has no k-nearest or ice identifier. Thresholds were not changed.

| Molecules | Operation | d-SEAMS 1 | d-SEAMS 8 | mdapy 1 | mdapy 8 | pyscal3 1 | pyscal3 8 | vesin 1 | vesin 8 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 4096 | cutoff | 0.282 | 0.233 | 1.939 | 0.553 | 1.939 | 1.615 | 1.672 | 0.768 |
| 4096 | knn | 0.570 | 0.509 | 3.788 | 0.821 | 15.49 | 5.525 | — | — |
| 4096 | q6 | 0.904 | 0.556 | 4.630 | 0.990 | 3.473 | 2.156 | — | — |
| 4096 | CHILL | 1.207 | 1.171 | — | — | — | — | — | — |
| 4096 | CHILL+ | 0.760 | 0.431 | 2.625 | 0.720 | — | — | — | — |
| 4096 | diamond | — | — | — | — | 4.442 | 2.368 | — | — |
| 32768 | cutoff | 2.703 | 2.036 | 12.42 | 2.795 | 17.66 | 5.779 | 13.50 | 5.798 |
| 32768 | knn | 3.417 | 3.369 | 30.89 | 5.281 | 109.4 | 25.82 | — | — |
| 32768 | q6 | 7.630 | 4.136 | 32.90 | 5.722 | 29.61 | 6.688 | — | — |
| 32768 | CHILL | 10.12 | 5.711 | — | — | — | — | — | — |
| 32768 | CHILL+ | 6.520 | 2.613 | 17.54 | 3.972 | — | — | — | — |
| 32768 | diamond | — | — | — | — | 28.66 | 7.551 | — | — |
| 110592 | cutoff | 8.418 | 4.619 | 39.91 | 6.402 | 60.26 | 12.68 | 45.78 | 15.50 |
| 110592 | knn | 10.35 | 10.11 | 110.0 | 15.54 | 464.9 | 86.59 | — | — |
| 110592 | q6 | 25.63 | 11.64 | 108.2 | 15.10 | 99.90 | 19.39 | — | — |
| 110592 | CHILL | 34.13 | 16.76 | — | — | — | — | — | — |
| 110592 | CHILL+ | 22.02 | 6.969 | 56.64 | 9.069 | — | — | — | — |
| 110592 | diamond | — | — | — | — | 96.37 | 22.88 | — | — |

q6 is the local Steinhardt parameter and the Lechner-Dellago average, timed together. pyscal3's ice row is extended CNA (`diamond_structure`), not CHILL. Original CHILL still builds the bond list; at eight threads that list is slower than mdapy's CHILL+ (16.8 ms against 9.1 ms at 110592 molecules). Published CHILL+ on the same shell is faster than mdapy at every size.
