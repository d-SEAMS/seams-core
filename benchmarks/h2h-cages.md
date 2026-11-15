# Cages

Pinned samples, 20 for d-SEAMS. Serial is pinned to one core. Eight threads are pinned to cores 0-7. The machine has 8 cores. The command is `seams cages` with cutoff 3.5 A and signatures `5:12`, `5:12,6:2`, `5:12,6:4`, `5:12,6:8` and `4:3,5:6,6:3`. The raw log is benchmarks/h2h-cages.txt.

SQQ 0.5.8 is the python cage engine. ClavisCage is absent: `clavis_engine` 0.3.18 is not on PyPI, and the package named clavis is a different client. SQQ stores an orthorhombic box, so the sH cell, which is tilted, is counted by d-SEAMS only.

Unit cells agree. sI (46 molecules): d-SEAMS finds 2 of `5:12` and 6 of `5:12,6:2`; SQQ finds 2 of `512` and 6 of `51262`. sII (136): 16 of `5:12` and 8 of `5:12,6:4` against SQQ 16 of `512` and 8 of `51264`. sH (136): 12 of `5:12`, 4 of `5:12,6:8` and 8 of `4:3,5:6,6:3`.

The large cells are 14x14x14 sI (126224), 10x10x9 sII (122400) and 7x14x9 sH (119952). Each cage count is the unit-cell count times the number of replicas.

At one core the medians are 625.32 ms on sI (IQR 26.20, minimum 607.73), 740.89 ms on sII (IQR 8.72, minimum 726.51) and 725.57 ms on sH (IQR 8.82, minimum 718.16). At eight threads they are 306.92 ms (IQR 13.51, minimum 299.27), 328.92 ms (IQR 10.52, minimum 316.90) and 322.36 ms (IQR 7.20, minimum 311.60).

SQQ does not reach these cells in a sample. On the 46-molecule sI cell, four samples at one core have median 474910 ms (IQR 3416, minimum 472564). One sample with eight workers took 471112 ms. The python walk stayed on one core, so the extra workers did not change it. A separate one-shot of SQQ's C++ engine on that same cell took 80.4 s and then 80.9 s, with the same two counts, so the repeated timing stayed on the python engine.
