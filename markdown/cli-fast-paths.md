# CLI times

Ice Ic, orthorhombic, one oxygen per molecule, one frame of 110592 molecules. The machine reported 8 cores. Each row is one warm-up and eight samples, median and interquartile range, with the process pinned to that many cores and the thread count set to the same number. CLI rows include process start. Library and package rows time the call after the frame is loaded.

mdapy 2.0.0, OVITO 3.17.2, freud 3.6.1, rings-statistics-matter 0.3.0. Measured in this session on this frame.

CHILL+ labels every molecule cubic in `seams chill-plus`, `seams_frame_chill_plus`, mdapy and OVITO. `seams symmetry` labels every molecule ice Ic, and OVITO's identify-diamond modifier labels every molecule cubic diamond. Polyhedral template matching with cubic and hexagonal diamond (PTM types 6 and 7) labels every molecule cubic diamond. The default modifier leaves both templates off and then reports Other. Local q4 and q6 from `seams_frame_q4_q6` and freud differ by at most 1.565e-07. Primitive six-rings: the library and rings-statistics-matter both report 221184. `seams cages` and `seams_frame_cages` label every molecule cubic.

| command | molecules | cores | median (ms) | IQR (ms) |
| --- | ---: | ---: | ---: | ---: |
| `seams read` | 110592 | 1 | 24.83 | 0.73 |
| `seams chill-plus` | 110592 | 1 | 62.48 | 0.74 |
| `seams symmetry` | 110592 | 1 | 62.83 | 1.48 |
| `seams steinhardt --l 4,6` | 110592 | 1 | 90.83 | 0.42 |
| `seams cages` | 110592 | 1 | 556.87 | 3.00 |
| `library chill-plus` | 110592 | 1 | 22.06 | 0.10 |
| `library q4 q6` | 110592 | 1 | 32.79 | 0.14 |
| `library cages` | 110592 | 1 | 871.14 | 2.42 |
| `library rings` | 110592 | 1 | 92.64 | 0.90 |
| `mdapy ChillPlus` | 110592 | 1 | 61.74 | 0.23 |
| `OVITO ChillPlus` | 110592 | 1 | 615.33 | 4.53 |
| `OVITO PTM` | 110592 | 1 | 1524.19 | 3.70 |
| `OVITO identify diamond` | 110592 | 1 | 153.19 | 2.84 |
| `freud q4` | 110592 | 1 | 259.96 | 4.79 |
| `freud q6` | 110592 | 1 | 344.94 | 3.81 |
| `rings-statistics-matter` | 110592 | 1 | 1101.20 | 16.00 |
| `seams read` | 110592 | 8 | 25.30 | 0.57 |
| `seams chill-plus` | 110592 | 8 | 49.70 | 2.08 |
| `seams symmetry` | 110592 | 8 | 40.57 | 1.12 |
| `seams steinhardt --l 4,6` | 110592 | 8 | 59.17 | 0.99 |
| `seams cages` | 110592 | 8 | 212.85 | 5.09 |
| `library chill-plus` | 110592 | 8 | 5.64 | 0.14 |
| `library q4 q6` | 110592 | 8 | 8.60 | 0.75 |
| `library cages` | 110592 | 8 | 282.25 | 10.54 |
| `library rings` | 110592 | 8 | 31.31 | 4.13 |
| `mdapy ChillPlus` | 110592 | 8 | 83.73 | 0.84 |
| `OVITO ChillPlus` | 110592 | 8 | 86.92 | 2.56 |
| `OVITO PTM` | 110592 | 8 | 209.28 | 1.42 |
| `OVITO identify diamond` | 110592 | 8 | 36.08 | 1.17 |
| `freud q4` | 110592 | 8 | 259.15 | 3.93 |
| `freud q6` | 110592 | 8 | 342.18 | 3.50 |
| `rings-statistics-matter` | 110592 | 8 | 1122.83 | 9.67 |
