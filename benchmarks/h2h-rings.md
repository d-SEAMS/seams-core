# Primitive rings

Franzblau shortest-path rings of size 3 to 6, cutoff 3.5 A. One session, eight samples, one thread and eight threads. Codes are d-SEAMS, rings-statistics-matter 0.3.0 (Rust backend), cycless 0.7, and OVITO 3.17.2 FindRings. The raw log is benchmarks/h2h-rings.txt.

`seams_rings_find` returns the ring list and the size label. Cage and prism codes are filled the first time `seams_rings_labels` reads symmetry, from the cloud kept on the handle. The timed call is the find plus the histogram.

## Counts

| Frame | Molecules | Size | d-SEAMS | rings-statistics-matter | OVITO | cycless |
|---|---:|---:|---:|---:|---:|---:|
| ice | 4096 | 6 | 8192 | 8192 | 8192 | 8192 |
| ice | 32768 | 6 | 65536 | 65536 | 65536 | not repeated |
| ice | 110592 | 6 | 221184 | 221184 | 221184 | not repeated |
| interface | 512 | 3 | 637 | 637 | -- | 637 |
| interface | 512 | 4 | 52 | 52 | 52 | 52 |
| interface | 512 | 5 | 55 | 55 | 55 | 48 |
| interface | 512 | 6 | 407 | 407 | 407 | 403 |
| stacking | 512 | 6 | 800 | 800 | 800 | 800 |
| melted | 512 | 3 | 191 | 191 | -- | 191 |
| melted | 512 | 4 | 21 | 21 | 21 | 21 |
| melted | 512 | 5 | 25 | 25 | 25 | 20 |
| melted | 512 | 6 | 807 | 807 | 807 | 807 |
| sII | 136 | 5 | 144 | 144 | 144 | 144 |
| sII | 136 | 6 | 16 | 16 | 16 | 16 |

Sizes that are zero on the ideal ice and on sII are omitted. d-SEAMS, rings-statistics-matter, and OVITO agree on every size OVITO reports.

OVITO's ring-size histogram starts at 4, so the 637 triangles on the interface and the 191 triangles on the melted frame are absent from its table. d-SEAMS, rings-statistics-matter, and cycless agree on those triangles.

cycless is short by 7 pentagons and 4 hexagons on the interface, and by 5 pentagons on the melted frame. `cycles_iter` grows a ring only as a shortest path between one pair of neighbours of a root, and that walk does not emit every primitive ring. The rings it misses are cycles in the same 3.5 A graph, and no pair of members has a graph distance shorter than the distance along the ring, so they are shortest-path rings. On the interface the witness pentagon, in cycle order, is atoms 123, 298, 260, 337, 383. `tests/test_cycless_pentagon.py` rebuilds that frame and checks the d-SEAMS census contains it. On ice, stacking, and sII, where the cycless walk does not drop a ring, the four codes match.

## Time, median of eight samples, milliseconds

| Frame | Threads | d-SEAMS | rings-statistics-matter | OVITO | cycless |
|---|---:|---:|---:|---:|---:|
| ice 4096 | 1 | 4.63 | 59.49 | 16.71 | -- |
| ice 4096 | 8 | 2.49 | 15.57 | 22.72 | -- |
| ice 32768 | 1 | 33.71 | 440.91 | 147.88 | -- |
| ice 32768 | 8 | 12.71 | 98.68 | 148.06 | -- |
| ice 110592 | 1 | 113.82 | 1398.53 | 669.81 | -- |
| ice 110592 | 8 | 48.02 | 361.45 | 651.43 | -- |
| interface | 1 | 1.81 | 13.52 | 3.30 | 7600.31 |
| interface | 8 | 0.68 | 3.71 | 3.15 | 7677.03 |
| stacking | 1 | 0.81 | 6.89 | 2.55 | 1310.43 |
| stacking | 8 | 0.54 | 2.20 | 2.41 | 1313.86 |
| melted | 1 | 1.28 | 9.54 | 3.03 | 2384.49 |
| melted | 8 | 0.51 | 2.74 | 2.84 | 2421.18 |
| sII | 1 | 0.27 | 2.20 | 1.03 | 131.08 |
| sII | 8 | 0.25 | 1.46 | 0.77 | 130.66 |

cycless is timed through 512 molecules. A 512-molecule interface already spends about 7.6 s in `cycles_iter`, and the same search is not repeated on the larger ice. d-SEAMS is the fastest census on every frame and at both thread counts.
