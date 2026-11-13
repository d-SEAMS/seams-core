# d-SEAMS, mdapy, pyscal3, and vesin

One process, eight samples, eight cores. Orthorhombic ice Ih, oxygen molecules only. Medians are milliseconds. The raw lines are `benchmarks/h2h-mdapy-pyscal3-vesin.txt`.

Unordered pairs inside 3.5 Å match exactly (8192, 65536, and 221184). Each oxygen's four nearest neighbours match exactly, and on this crystal they are that cutoff shell. q6 and q6bar agree within 1e-8 (the largest absolute difference is 2.4e-15). Published CHILL, published CHILL+, and mdapy CHILL+ label every oxygen hexagonal. pyscal3 has no CHILL; its extended-CNA diamond count is hexagonal diamond for every oxygen. vesin builds a cutoff list and has no k-nearest or ice identifier. Thresholds were not changed.

| Molecules | Operation | d-SEAMS 1 | d-SEAMS 8 | mdapy 1 | mdapy 8 | pyscal3 1 | pyscal3 8 | vesin 1 | vesin 8 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 4096 | cutoff | 0.320 | 0.228 | 2.154 | 0.571 | 2.051 | 2.137 | 1.847 | 1.157 |
| 4096 | knn | 0.685 | 0.665 | 4.610 | 0.890 | 17.01 | 6.971 | — | — |
| 4096 | q6 | 0.987 | 0.810 | 5.090 | 1.075 | 3.761 | 2.945 | — | — |
| 4096 | CHILL | 0.914 | 0.488 | — | — | — | — | — | — |
| 4096 | CHILL+ | 0.866 | 0.481 | 2.779 | 0.776 | — | — | — | — |
| 4096 | diamond | — | — | — | — | 4.581 | 2.598 | — | — |
| 32768 | cutoff | 2.813 | 1.728 | 12.82 | 3.808 | 18.38 | 7.196 | 14.62 | 6.160 |
| 32768 | knn | 4.028 | 3.511 | 32.52 | 5.682 | 137.1 | 30.48 | — | — |
| 32768 | q6 | 8.106 | 4.087 | 37.45 | 5.841 | 30.43 | 9.608 | — | — |
| 32768 | CHILL | 7.365 | 2.557 | — | — | — | — | — | — |
| 32768 | CHILL+ | 6.899 | 2.950 | 18.48 | 4.106 | — | — | — | — |
| 32768 | diamond | — | — | — | — | 29.58 | 9.399 | — | — |
| 110592 | cutoff | 8.709 | 5.805 | 41.35 | 8.171 | 70.29 | 21.24 | 46.72 | 18.66 |
| 110592 | knn | 12.43 | 11.29 | 112.7 | 18.19 | 590.1 | 107.2 | — | — |
| 110592 | q6 | 27.43 | 14.83 | 123.7 | 20.45 | 126.3 | 23.89 | — | — |
| 110592 | CHILL | 24.08 | 10.41 | — | — | — | — | — | — |
| 110592 | CHILL+ | 22.57 | 9.307 | 74.18 | 13.18 | — | — | — | — |
| 110592 | diamond | — | — | — | — | 101.6 | 22.51 | — | — |

q6 is the local Steinhardt parameter and the Lechner-Dellago average, timed together. pyscal3's ice row is extended CNA (`diamond_structure`), not CHILL. Published CHILL and CHILL+ both read the 3.5 Å cell-list shell and do not store a bond vector per atom. At eight threads and 110592 molecules, CHILL is 10.41 ms and CHILL+ is 9.31 ms; mdapy CHILL+ is 13.18 ms.
