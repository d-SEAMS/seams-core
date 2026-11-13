# d-SEAMS, mdapy, pyscal3, vesin, and OVITO

One process, eight samples, eight cores. Orthorhombic ice Ih, oxygen molecules only. Medians are milliseconds. The raw lines are `benchmarks/h2h-mdapy-pyscal3-vesin.txt`.

Unordered pairs inside 3.5 A match exactly (8192, 65536, and 221184). Each oxygen's four nearest neighbours match exactly, and on this crystal they are that cutoff shell. q6 and q6bar agree within 1e-8 (the largest absolute difference is 2.442e-15). Paper CHILL+ averages q_lm over the four closest neighbours. Ovito CHILL+ averages q_lm over every neighbour inside 3.5 A. On this ice the two shells are the same, and both modes, mdapy CHILL+, and the OVITO Chill+ modifier label every oxygen hexagonal. The same ovito mode matches mdapy and the modifier atom for atom on the mixed interface, stacking, and melted frames. pyscal3 has no CHILL; its extended-CNA diamond count is hexagonal diamond for every oxygen. vesin builds a cutoff list and has no k-nearest or ice identifier. Thresholds were not changed.

| Molecules | Operation | d-SEAMS 1 | d-SEAMS 8 | mdapy 1 | mdapy 8 | pyscal3 1 | pyscal3 8 | vesin 1 | vesin 8 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 4096 | cutoff | 0.2462 | 0.2063 | 2.211 | 0.4723 | 2.138 | 1.620 | 1.829 | 0.8214 |
| 4096 | knn | 0.5211 | 0.6070 | 4.135 | 0.8684 | 17.16 | 6.143 | — | — |
| 4096 | q6 | 0.9719 | 0.5143 | 5.253 | 0.9478 | 3.816 | 2.247 | — | — |
| 4096 | CHILL | 0.8562 | 0.3830 | — | — | — | — | — | — |
| 4096 | CHILL+ paper | 1.123 | 0.3779 | — | — | — | — | — | — |
| 4096 | CHILL+ ovito | 1.102 | 0.3350 | 2.874 | 0.6360 | — | — | — | — |
| 4096 | diamond | — | — | — | — | 4.901 | 2.246 | — | — |
| 32768 | cutoff | 2.397 | 1.715 | 13.81 | 2.657 | 20.61 | 5.924 | 15.25 | 6.446 |
| 32768 | knn | 3.637 | 3.959 | 35.68 | 5.432 | 150.9 | 28.77 | — | — |
| 32768 | q6 | 8.619 | 4.112 | 43.36 | 6.047 | 32.94 | 8.550 | — | — |
| 32768 | CHILL | 7.590 | 2.510 | — | — | — | — | — | — |
| 32768 | CHILL+ paper | 9.621 | 2.791 | — | — | — | — | — | — |
| 32768 | CHILL+ ovito | 9.103 | 2.471 | 19.83 | 3.710 | — | — | — | — |
| 32768 | diamond | — | — | — | — | 31.21 | 8.081 | — | — |
| 110592 | cutoff | 7.065 | 5.784 | 47.20 | 8.864 | 75.12 | 18.40 | 49.84 | 20.21 |
| 110592 | knn | 11.77 | 11.95 | 131.8 | 21.75 | 626.5 | 119.5 | — | — |
| 110592 | q6 | 29.17 | 15.31 | 127.5 | 23.55 | 140.8 | 24.11 | — | — |
| 110592 | CHILL | 24.70 | 6.825 | — | — | — | — | — | — |
| 110592 | CHILL+ paper | 22.24 | 6.386 | — | — | — | — | — | — |
| 110592 | CHILL+ ovito | 25.94 | 6.342 | 71.51 | 14.59 | — | — | — | — |
| 110592 | diamond | — | — | — | — | 131.9 | 23.96 | — | — |

q6 is the local Steinhardt parameter and the Lechner-Dellago average, timed together. pyscal3's ice row is extended CNA (`diamond_structure`), not CHILL. The mdapy column on the CHILL+ ovito row is mdapy CHILL+. Paper mode is the default. The serial cutoff at 110592 is 7.065 ms. At eight threads and 110592 molecules, paper CHILL+ is 6.386 ms and ovito CHILL+ is 6.342 ms; mdapy CHILL+ is 14.59 ms.
