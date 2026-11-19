<!-- https://docs.dseams.info/user/howto/cutoff-padding.md -->

# Choosing a cutoff and padding {#choosing-a-cutoff-and-padding}

The neighbour cutoff is the `--cutoff` of `seams chill-plus`.
The default is 3.5 Å. `SEAMS_CUTOFF` sets that default when the flag
is omitted.
Open-face padding is the `openFace` argument of `symm::classifyIce`,
`symm::classifyIceCloud`, and `symm::classifyIceArrays`.
It completes a short second shell on a bonded neighbour that has fewer
than four bonds. A full coordination shell is never padded.
The argument defaults to false.
`seams symmetry` calls `classifyIce` with that default, so the command
leaves the completion off.
CHILL+ stays the published rule. It does not take `openFace`.

Pass `openFace` true when the bulk-liquid rows record no added false
positive. On an interface, read the layers whose label changed before
treating the extra shell labels as ice. The classify timings are the
same bulk liquid with the completion off and on.

A false positive is a non-negative shell label on an atom CHILL+ does
not call cubic or hexagonal. `dseams_off` and `fp_off` are the
classifier with the completion off. `dseams_on` and `fp_on` are the
same frame with it on.

## The bulk-liquid gate {#the-bulk-liquid-gate}

padding adds no false positive on the bulk liquid; the open-face completion stays on

False positives off are `0` and false positives on are `0`.

## Summary rows {#summary-rows}

Each row is one `summary` record. `frame` is the number of
frames that row adds.

| lane | temperature_K | frame | molecules | chill_ice | dseams_off | dseams_on | fp_off | fp_on | kind_changed | agree_off | agree_on | source |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mw-bulk | 230 | 10 | 40960 | 108 | 0 | 0 | 0 | 0 | 0 | 40852 | 40852 | sum of the ten hold frames |
| mw-bulk | 300 | 10 | 40960 | 19 | 0 | 0 | 0 | 0 | 0 | 40941 | 40941 | sum of the ten hold frames |
| tip4p-basal | 230 | 10 | 14640 | 9390 | 8491 | 8669 | 13 | 16 | 178 | 13715 | 13887 | sum of the ten basal frames |
| T230_m0_r0_t0 | 230 | 1 | 1536 | 776 | 743 | 745 | 0 | 1 | 2 | 1503 | 1503 | one pinned brine frame |
| T230_m0_r0_t1 | 230 | 1 | 1536 | 844 | 511 | 511 | 2 | 2 | 0 | 1199 | 1199 | one pinned brine frame |
| T230_m0_r0_t2 | 230 | 1 | 1536 | 832 | 539 | 539 | 0 | 0 | 0 | 1243 | 1243 | one pinned brine frame |
| T230_m0_r0_tend | 230 | 1 | 1536 | 856 | 577 | 577 | 1 | 1 | 0 | 1255 | 1255 | one pinned brine frame |
| T230_m20_r0_t0 | 230 | 1 | 1496 | 734 | 693 | 706 | 2 | 4 | 13 | 1451 | 1460 | one pinned brine frame |
| T230_m20_r0_t1 | 230 | 1 | 1496 | 698 | 420 | 426 | 0 | 0 | 6 | 1218 | 1224 | one pinned brine frame |
| T230_m20_r0_t2 | 230 | 1 | 1496 | 698 | 447 | 447 | 2 | 2 | 0 | 1241 | 1241 | one pinned brine frame |
| T230_m20_r0_tend | 230 | 1 | 1496 | 741 | 479 | 481 | 2 | 2 | 2 | 1230 | 1232 | one pinned brine frame |
| T230_m40_r0_t0 | 230 | 1 | 1456 | 707 | 671 | 687 | 3 | 5 | 16 | 1414 | 1426 | one pinned brine frame |
| T230_m40_r0_t1 | 230 | 1 | 1456 | 664 | 408 | 413 | 2 | 2 | 5 | 1196 | 1201 | one pinned brine frame |
| T230_m40_r0_t2 | 230 | 1 | 1456 | 682 | 369 | 375 | 0 | 0 | 6 | 1143 | 1149 | one pinned brine frame |
| T230_m40_r0_tend | 230 | 1 | 1456 | 667 | 414 | 415 | 0 | 0 | 1 | 1203 | 1204 | one pinned brine frame |

## Layers whose label changed {#layers-whose-label-changed}

These are the `layer` records whose `kind_changed` cell is
not `0`. `z_lo_A` is the lower edge of the slab.

| lane | frame | layer | z_lo_A | molecules | dseams_off | dseams_on | fp_off | fp_on | kind_changed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tip4p-basal | 8 | 1 | 2.000000 | 128 | 20 | 33 | 2 | 2 | 13 |
| tip4p-basal | 8 | 14 | 28.000000 | 127 | 14 | 43 | 1 | 2 | 29 |
| tip4p-basal | 16 | 1 | 2.000000 | 128 | 41 | 43 | 1 | 1 | 2 |
| tip4p-basal | 16 | 14 | 28.000000 | 128 | 21 | 42 | 1 | 1 | 21 |
| tip4p-basal | 24 | 1 | 2.000000 | 128 | 51 | 52 | 0 | 0 | 1 |
| tip4p-basal | 24 | 14 | 28.000000 | 128 | 29 | 38 | 0 | 0 | 9 |
| tip4p-basal | 32 | 14 | 28.000000 | 127 | 26 | 41 | 1 | 1 | 15 |
| tip4p-basal | 40 | 14 | 28.000000 | 128 | 26 | 39 | 0 | 0 | 13 |
| tip4p-basal | 48 | 13 | 26.000000 | 1 | 0 | 1 | 0 | 0 | 1 |
| tip4p-basal | 48 | 14 | 28.000000 | 127 | 25 | 42 | 1 | 2 | 17 |
| tip4p-basal | 56 | 14 | 28.000000 | 127 | 24 | 38 | 1 | 2 | 14 |
| tip4p-basal | 64 | 1 | 2.000000 | 128 | 66 | 68 | 0 | 0 | 2 |
| tip4p-basal | 64 | 14 | 28.000000 | 126 | 22 | 32 | 1 | 1 | 10 |
| tip4p-basal | 72 | 14 | 28.000000 | 128 | 27 | 42 | 2 | 2 | 15 |
| tip4p-basal | 80 | 14 | 28.000000 | 128 | 23 | 39 | 0 | 0 | 16 |
| T230_m0_r0_t0 | 1 | 0 | 0.000000 | 64 | 46 | 47 | 0 | 1 | 1 |
| T230_m0_r0_t0 | 1 | 12 | 24.000000 | 64 | 52 | 53 | 0 | 0 | 1 |
| T230_m20_r0_t0 | 1 | 0 | 0.000000 | 64 | 27 | 28 | 2 | 2 | 1 |
| T230_m20_r0_t0 | 1 | 1 | 2.000000 | 64 | 61 | 62 | 0 | 0 | 1 |
| T230_m20_r0_t0 | 1 | 11 | 22.000000 | 64 | 61 | 64 | 0 | 0 | 3 |
| T230_m20_r0_t0 | 1 | 12 | 24.000000 | 64 | 30 | 38 | 0 | 2 | 8 |
| T230_m20_r0_t1 | 1 | 0 | 0.000000 | 54 | 12 | 13 | 0 | 0 | 1 |
| T230_m20_r0_t1 | 1 | 1 | 2.000000 | 47 | 18 | 19 | 0 | 0 | 1 |
| T230_m20_r0_t1 | 1 | 2 | 4.000000 | 37 | 15 | 18 | 0 | 0 | 3 |
| T230_m20_r0_t1 | 1 | 12 | 24.000000 | 54 | 16 | 17 | 0 | 0 | 1 |
| T230_m20_r0_tend | 1 | 0 | 0.000000 | 56 | 22 | 23 | 0 | 0 | 1 |
| T230_m20_r0_tend | 1 | 25 | 50.000000 | 69 | 4 | 5 | 2 | 2 | 1 |
| T230_m40_r0_t0 | 1 | 0 | 0.000000 | 64 | 21 | 28 | 0 | 2 | 7 |
| T230_m40_r0_t0 | 1 | 1 | 2.000000 | 64 | 61 | 64 | 0 | 0 | 3 |
| T230_m40_r0_t0 | 1 | 11 | 22.000000 | 64 | 59 | 63 | 0 | 0 | 4 |
| T230_m40_r0_t0 | 1 | 12 | 24.000000 | 64 | 18 | 20 | 3 | 3 | 2 |
| T230_m40_r0_t1 | 1 | 2 | 4.000000 | 53 | 24 | 25 | 1 | 1 | 1 |
| T230_m40_r0_t1 | 1 | 3 | 6.000000 | 61 | 48 | 49 | 1 | 1 | 1 |
| T230_m40_r0_t1 | 1 | 12 | 24.000000 | 65 | 13 | 15 | 0 | 0 | 2 |
| T230_m40_r0_t1 | 1 | 13 | 26.000000 | 65 | 0 | 1 | 0 | 0 | 1 |
| T230_m40_r0_t2 | 1 | 0 | 0.000000 | 43 | 7 | 8 | 0 | 0 | 1 |
| T230_m40_r0_t2 | 1 | 1 | 2.000000 | 45 | 10 | 13 | 0 | 0 | 3 |
| T230_m40_r0_t2 | 1 | 13 | 26.000000 | 63 | 1 | 3 | 0 | 0 | 2 |
| T230_m40_r0_tend | 1 | 12 | 24.000000 | 64 | 18 | 19 | 0 | 0 | 1 |

## Classify timings {#classify-timings}

`padding` is `0` when the open-face completion is off and
`1` when it is on. `seconds` is one replicate.

| molecules | padding | rep | seconds | cores_online | cores_affinity | cores_used | source |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 4096 | 0 | 1 | 0.000398 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 1 | 1 | 0.000389 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 0 | 2 | 0.000375 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 1 | 2 | 0.000370 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 0 | 3 | 0.000365 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 1 | 3 | 0.000367 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 0 | 4 | 0.000379 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 1 | 4 | 0.000384 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 0 | 5 | 0.000366 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 1 | 5 | 0.000375 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 0 | 6 | 0.000362 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 1 | 6 | 0.000367 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 0 | 7 | 0.000391 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 1 | 7 | 0.000396 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 0 | 8 | 0.000398 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 4096 | 1 | 8 | 0.000371 | 8 | 8 | 8 | 300 K bulk-liquid frame |
| 110592 | 0 | 1 | 0.039944 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 1 | 1 | 0.036897 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 0 | 2 | 0.041904 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 1 | 2 | 0.039008 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 0 | 3 | 0.040180 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 1 | 3 | 0.035358 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 0 | 4 | 0.034795 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 1 | 4 | 0.018831 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 0 | 5 | 0.009178 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 1 | 5 | 0.028915 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 0 | 6 | 0.028953 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 1 | 6 | 0.038089 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 0 | 7 | 0.020041 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 1 | 7 | 0.025629 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 0 | 8 | 0.031073 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
| 110592 | 1 | 8 | 0.015290 | 8 | 8 | 8 | 3x3x3 tile of the 300 K bulk-liquid frame |
