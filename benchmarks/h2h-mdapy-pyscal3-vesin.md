# d-SEAMS, mdapy, pyscal3, vesin, and OVITO

One process, eight samples, eight cores. Orthorhombic ice Ih, oxygen molecules only. Medians are milliseconds. The raw lines are `benchmarks/h2h-mdapy-pyscal3-vesin.txt`.

Unordered pairs inside 3.5 A match exactly (8192, 65536, and 221184). Each oxygen's four nearest neighbours match exactly, and on this crystal they are that cutoff shell. q6 and q6bar agree within 1e-8 (the largest absolute difference is 2.442e-15). Paper CHILL+ averages q_lm over the four closest neighbours. Ovito CHILL+ averages q_lm over every neighbour inside 3.5 A. On this ice the two shells are the same, and both modes, mdapy CHILL+, and the OVITO Chill+ modifier label every oxygen hexagonal. The same ovito mode matches mdapy and the modifier atom for atom on the mixed interface, stacking, and melted frames. pyscal3 has no CHILL; its extended-CNA diamond count is hexagonal diamond for every oxygen. vesin builds a cutoff list and has no k-nearest or ice identifier. Thresholds were not changed.

| Molecules | Operation | d-SEAMS 1 | d-SEAMS 8 | mdapy 1 | mdapy 8 | pyscal3 1 | pyscal3 8 | vesin 1 | vesin 8 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 4096 | cutoff | 0.2785 | 0.2549 | 1.978 | 0.6336 | 1.925 | 1.605 | 1.681 | 0.8267 |
| 4096 | knn | 0.5871 | 0.6319 | 3.809 | 0.8289 | 15.31 | 6.065 | — | — |
| 4096 | q6 | 0.9753 | 0.6191 | 4.655 | 1.069 | 3.588 | 2.246 | — | — |
| 4096 | CHILL | 0.8209 | 0.4655 | — | — | — | — | — | — |
| 4096 | CHILL+ paper | 0.7715 | 0.4569 | — | — | — | — | — | — |
| 4096 | CHILL+ ovito | 0.7173 | 0.4224 | 2.639 | 0.8103 | — | — | — | — |
| 4096 | diamond | — | — | — | — | 4.433 | 2.754 | — | — |
| 32768 | cutoff | 2.659 | 1.613 | 12.67 | 3.041 | 17.79 | 5.443 | 13.54 | 7.842 |
| 32768 | knn | 3.764 | 3.847 | 30.90 | 4.941 | 111.1 | 28.72 | — | — |
| 32768 | q6 | 8.214 | 4.171 | 32.91 | 5.877 | 29.39 | 7.937 | — | — |
| 32768 | CHILL | 6.977 | 2.755 | — | — | — | — | — | — |
| 32768 | CHILL+ paper | 6.540 | 2.657 | — | — | — | — | — | — |
| 32768 | CHILL+ ovito | 6.195 | 2.711 | 17.59 | 4.171 | — | — | — | — |
| 32768 | diamond | — | — | — | — | 28.24 | 9.078 | — | — |
| 110592 | cutoff | 8.543 | 4.833 | 40.22 | 8.599 | 60.45 | 15.56 | 45.70 | 24.60 |
| 110592 | knn | 11.23 | 11.67 | 110.3 | 17.03 | 492.8 | 101.0 | — | — |
| 110592 | q6 | 27.12 | 12.45 | 108.9 | 17.45 | 99.43 | 21.87 | — | — |
| 110592 | CHILL | 23.32 | 7.929 | — | — | — | — | — | — |
| 110592 | CHILL+ paper | 21.99 | 7.723 | — | — | — | — | — | — |
| 110592 | CHILL+ ovito | 20.72 | 7.492 | 57.21 | 11.39 | — | — | — | — |
| 110592 | diamond | — | — | — | — | 95.09 | 24.78 | — | — |

q6 is the local Steinhardt parameter and the Lechner-Dellago average, timed together. pyscal3's ice row is extended CNA (`diamond_structure`), not CHILL. The mdapy column on the CHILL+ ovito row is mdapy CHILL+. Paper mode is the default. At eight threads and 110592 molecules, paper CHILL+ is 7.723 ms and ovito CHILL+ is 7.492 ms; mdapy CHILL+ is 11.39 ms.
