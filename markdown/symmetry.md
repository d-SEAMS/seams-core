# Local symmetry of ice, cages and rings

Copyright (c) d-SEAMS developers

d-SEAMS can match a local environment to an ideal ice shell, a clathrate
cage or a ring, and can name the point group of that cluster and the space
group of a crystal. The match is optional at build time.

## Which code the recollection points to

The Japanese code is spglib, by Atsushi Togo.

The other codes in this note were checked against that recollection.

- spglib is the crystal-symmetry library written by Atsushi Togo (National
  Institute for Materials Science). Its methods paper is Togo, Shinohara and
  Tanaka, Science and Technology of Advanced Materials: Methods 4, 2384822
  (2024), doi:10.1080/27660400.2024.2384822. pymatgen's `SpacegroupAnalyzer`
  is documented as an interface to that library, not as a separate symmetry
  engine.
- moyo, and its Python binding moyopy, are also by a Japanese author, Kohei
  Shinohara, who coauthors the 2024 spglib paper. The moyo repository was
  created in October 2023 inside the spglib organisation, and the citation
  the project asks for is a January 2026 figshare deposit
  (doi:10.6084/m9.figshare.31081162). A rewrite that recent is named moyo.
  The code one remembers as "a Japanese code" for space groups, from the
  years this engine has called external symmetry tools, is spglib.
- GenIce and cycless are Masakazu Matsumoto's ice-network codes (Okayama).
  This repository already cites them as ring and lattice generators
  (doi:10.1002/jcc.25077, doi:10.1063/5.0198056). They build hydrogen-bond
  networks and count cycles. They do not find point groups or space groups,
  so they are not the symmetry code in that recollection.
- SHAPE and the continuous shape measure are not Japanese. The measure is
  Pinsky and Avnir, Inorganic Chemistry 37, 5575 (1998),
  doi:10.1021/ic9804925. The Alvarez program is the Barcelona implementation.

## What this engine calls

| Piece | Role here | Licence we found | Vendored? |
| --- | --- | --- | --- |
| IRA and SOFI | Rotation, permutation and RMSD; point groups of cages and rings | GPL-3.0-or-later OR Apache-2.0 (MAMMASMIAS, `LICENSE.txt`) | No. Linked when `with_ira` finds libira. The Apache-2.0 choice is compatible with MIT. Taking the GPL-3.0 choice instead would conflict with shipping the combined work under MIT only, so that choice is not used. |
| spglib | Space groups of crystalline regions | BSD-3-Clause | No. Linked when `with_spglib` finds libsymspg. BSD-3-Clause is compatible with MIT. |
| moyo / moyopy | Same crystal-symmetry job as spglib, in Rust | Apache-2.0 for moyo; MIT OR Apache-2.0 for moyopy | No. Not linked. It duplicates spglib. |
| pymatgen `SpacegroupAnalyzer` and `PointGroupAnalyzer` | Python interface over spglib, plus a molecular point-group search | MIT | No. Not linked. The space-group call is the spglib call. |
| continuous-symmetry-measure/shape | Continuous shape measure reference program | GPL-2.0 | Optional `with_csm` (default auto). Built when OpenBabel is present, as its own executable. Not linked into libyodaLib. The in-tree measure is the published formula. |
| SHAPE (Llunell, Casanova, Cirera, Alemany, Alvarez) | Continuous shape measures of polyhedra | No OSI licence in the binary distribution we could inspect | No. |
| OVITO | Polyhedral template matching in the GUI and Python module | GPL-3.0 | No. Not linked. A benchmark may call an installed module. |
| polyhedral-template-matching (P. M. Larsen) | The PTM reference implementation | MIT | No. Not linked. |
| GenIce2 | Ice and clathrate lattices | MIT | No. |
| cycless | Ring cycles used by GenIce | No `LICENSE` file at the repository root we fetched | No. |

Ideal templates are generated in `local_symmetry.cpp`. No trajectory frames are stored in the tree.

## Templates

Oxygen–oxygen nearest-neighbour distance 2.75 Å.

- ice Ic: the central oxygen, its four tetrahedral neighbours, and the twelve second neighbours of the cubic diamond net.
- ice Ih: the same seventeen-atom shell with one bond eclipsed, the geometry of hexagonal ice.
- ice Isd: the same shell with two eclipsed bonds, the local geometry of a stacking fault. CHILL+ has no label for that shell. It can call the oxygen cubic, hexagonal, interfacial or water, and it does not name the fault.
- 5^12: regular dodecahedron, the small cage of both sI and sII.
- 5^12 6^2: the D6d tetrakaidecahedron, the large cage of sI.
- 5^12 6^4: the Td hexakaidecahedron, the large cage of sII. Three orbit seeds are expanded by the tetrahedral group.
- ring-5 and ring-6: regular pentagon and hexagon.

A bulk oxygen is compared with the ice shells by taking itself and its sixteen nearest periodic neighbours, and only when the sixteenth neighbour lies inside 5 Å. sI and sII are told apart by which large-cage template wins. CHILL+ labels both cage vertex sets from the bond list alone and does not carry a cage name.

## Build

```
meson setup build -Dwith_ira=auto -Dwith_spglib=auto
```

`auto` links the library when pkg-config or the linker can see it. `enabled` fails the configure when the library is missing. This build links IRA under the Apache-2.0 option of its dual licence. spglib stays optional. The per-atom ice label does not call either library: a missing libira or libsymspg does not change those labels. Point groups still need SOFI, and space groups still need spglib.

The command is `seams symmetry`. The C entry points are `seams_symmetry_match`, `seams_symmetry_point_group`, `seams_symmetry_space_group`, `seams_symmetry_template`, `seams_symmetry_crystal` and `seams_symmetry_classify` in `seams_c_api.h`. `seams_symmetry_classify_tol` takes the eclipsed-cosine cutoff. Python is `bindings/python/symmetry.py`.

## Kernel

The production label of a bulk oxygen is the ice template with the lowest scale-free residual. A periodic cell list gathers neighbours inside 5.6 Å, and unwrapped coordinates are folded into the primary cell before the stencil runs. The four nearest atoms inside 3.85 Å are the tetrahedron. The twelve second neighbours are the atoms bonded to that tetrahedron (within 3.8 Å of a first neighbour), not merely the twelve nearest to the centre: a thermal shell otherwise drops a real bond and keeps an atom that belongs to no arm. Each proper orientation of the tetrahedron is a quaternion Kabsch (Coutsias, Seok and Dill 2004; Theobald 2005). A mirror flips the scalar triple product and is skipped. Around each bond the three second neighbours are ordered by azimuth and by the three cyclic shifts of that order, which is the graph automorphism once the tetrahedron is fixed (Larsen, Schmidt and Schiøtz 2016). One 17-point quaternion then supplies the residual after the optimal rotation and the optimal scale, divided by the template radius. The ceiling is 0.12. It was chosen on the 13824-molecule seeds 91001 and 92002, against the majority of this kernel, CHILL+ and PTM on the steepest-descent quench, and seed 93003 was not used to pick it. A shell that is still far from the first template after three proper orientations is not ice, so a clathrate or a liquid does not pay for the other templates. The argument of `classifyIce` is that scale-free ceiling. Zero is a perfect shape.

A named shell, cage or ring still returns a rotation, a permutation and an RMSD. That alignment is a 3×3 SVD of the covariance (Kabsch), tried proper and improper, on a permutation taken from the bond graph rather than from N!. Rings use the cyclic order. Ice uses the 24 assignments of the first shell and then the bonded second neighbours. Cages grow a degree-3 matching from each vertex. The environments are independent, so `matchAll` runs them in parallel. The continuous shape measure is the scaled residual of that overlay, 100 times the mean square distance to the template divided by the mean square radius of the template (Pinsky and Avnir 1998; Alvarez et al. 2005). It is computed in-tree. The GPL-2 shape program is a separate executable, built when OpenBabel is present, and is not on this path.

## Speed

Release `-O3`, one warm-up, then six calls. The median of an even sample is the mean of the two central values. Quartiles are linear (type 7). Every oxygen of ideal ice Ic is labelled ice Ic. OVITO 3.17.2 PTM has cubic and hexagonal diamond enabled. It does not speed up when `OMP_NUM_THREADS` goes from 1 to 8. The kernel does.

| Molecules | Threads | Kernel median (ms) | Kernel per atom (µs) | PTM median (ms) | PTM per atom (µs) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 216 | 1 | 0.250 | 1.16 | 8.447 | 39.11 |
| 216 | 8 | 0.060 | 0.279 | 8.572 | 39.69 |
| 10648 | 1 | 15.239 | 1.43 | 72.986 | 6.85 |
| 10648 | 8 | 2.268 | 0.213 | 73.003 | 6.86 |
| 110592 | 1 | 88.015 | 0.796 | 545.459 | 4.93 |
| 110592 | 8 | 17.167 | 0.155 | 539.262 | 4.88 |

Raw kernel samples, 1 thread, milliseconds, 110592: 85.6238 86.7279 87.5105 88.5188 89.0185 92.2193. Eight threads, 110592: 16.5615 16.9311 16.9432 17.3903 17.5342 20.3109. The PTM column is the previous measurement of these ideal cells. The shell is the four nearest neighbours on the 3.5 Å bond list and their bonded neighbours. Before that walk, a tetrahedron whose bond-orientational order q is below 0.88 is left unlabelled. q is 1 − (3/8) Σ (cos θ + 1/3)² over the six angles of the four bonds. When the twelve second neighbours sit on one radius, the staggered and eclipsed count labels the shell directly: four staggered bonds are cubic diamond (ice Ic) and one eclipsed bond is hexagonal diamond (ice Ih). That is the second-shell test used by OVITO's diamond identifier and by pyscal3. A shell that is not that clean goes to the permutation fit. One thread at 110592 is 0.796 µs/atom. Eight threads are 5.1 times faster.

IRA on the same 216-atom cell, one shell per oxygen through libira, is serial. Median 296.31 ms (raw 284.602 288.204 291.98 300.632 325.416 335.735), about 1.37 ms per atom, three orders of magnitude above the kernel. spglib on the 8-atom ice Ic cell returns 227 `Fd-3m` in a median 4.72 ms. That is one cell, not a per-atom cost.

The naming path is the SVD and the bond-graph permutation, timed on rotated copies of the templates rather than on a neighbour search. IRA is the same shells through libira, serial. PTM does not name a cage or a ring. The PTM number below is structure identification of a 216-atom ice Ic crystal, from the same protocol.

| Path | Shells | Threads | Median (ms) | Per shell (µs) |
| --- | ---: | ---: | ---: | ---: |
| Kernel, ice Ic | 216 | 1 | 35.215 | 163 |
| Kernel, ice Ic | 216 | 8 | 5.290 | 24.5 |
| Kernel, ice Ic | 2048 | 1 | 334.85 | 164 |
| Kernel, ice Ic | 2048 | 8 | 44.712 | 21.8 |
| IRA, ice Ic | 216 | 1 | 159.26 | 737 |
| Kernel, 5^12 | 216 | 1 | 1.041 | 4.82 |
| Kernel, 5^12 | 216 | 8 | 19.502 | 90.3 |
| IRA, 5^12 | 64 | 1 | 22.427 | 350 |
| Kernel, 5^12 6^2 | 216 | 1 | 66.773 | 309 |
| Kernel, 5^12 6^2 | 216 | 8 | 9.177 | 42.5 |
| Kernel, ring-6 | 216 | 1 | 0.822 | 3.80 |
| Kernel, ring-6 | 216 | 8 | 0.588 | 2.72 |
| Kernel, ring-6 | 2048 | 1 | 43.738 | 21.4 |
| Kernel, ring-6 | 2048 | 8 | 5.576 | 2.72 |
| IRA, ring-6 | 216 | 1 | 3.130 | 14.5 |
| In-tree CShM, ice Ic | 216 | 1 | 35.197 | 163 |
| In-tree CShM, ice Ic | 216 | 8 | 4.479 | 20.7 |
| Shape program, one ice shell, identity permutation | 1 | 1 | 9.884 | 9884 |
| Shape program, tetrahedron, all 24 permutations | 1 | 1 | 10.360 | 10360 |
| PTM structure id, ice Ic crystal | 216 atoms | 1 | 8.573 | 39.7 per atom |

On one thread the hexagon is 3.80 µs a shell and the dodecahedron is 4.82 µs, against 14.5 µs and 357 µs for IRA. The quaternion replaces the 3×3 SVD, and a cage returns on the first graph match inside 0.08 Å. The shape-program row is one process, OpenBabel startup included. The in-tree measure on the same distorted ice shell agrees with that program to 0.05. The program printed `SHAPE = 0.043800` for a shell with one oxygen moved by a few tenths of an angstrom. Raw kernel samples, ice Ic, 1 thread, milliseconds: 216 shells 35.1267 35.1646 35.1953 35.2352 35.4733 37.1024; 2048 shells 333.275 333.691 334.234 335.466 335.909 336.319. Eight threads: 216 shells 5.0707 5.16434 5.1825 5.39789 5.43898 5.46815; 2048 shells 42.8295 43.3195 43.9354 45.4883 47.708 47.8798. IRA ice, 216 shells, 1 thread: 158.408 159.203 159.251 159.268 160.189 163.331. In-tree CShM, 216 shells, 1 thread: 35.0172 35.1007 35.1804 35.214 35.275 35.9081. Eight threads: 4.42846 4.44697 4.45494 4.50364 4.54158 4.61616. Shape program, ice identity: 9.768 9.851 9.882 9.886 9.970 10.429. Tetrahedron: 9.886 10.043 10.251 10.469 10.508 10.929. PTM, 216 atoms: 8.515148 8.631621 8.632679 8.681955 8.433761 8.488970.

## What PTM does not name

On the stacking-fault cell, PTM returns OTHER for all 8 oxygens. The IRA shell match returns ice Ic for 4 and ice Isd for 4. CHILL+ on those 8 oxygens taken out of the periodic cell, so the crystal neighbours are missing, returns water for all 8. It has no stacking-fault label.

On the sI oxygen lattice (46 waters, cell 12.239 Å) PTM returns OTHER for every water and CHILL+ returns clathrate for every water. The ring census finds 48 five-rings and 54 six-rings, and the closed cages are 2 of 5^12 and 6 of 5^12 6^2. On the sII lattice (136 waters, nearest neighbour scaled to 2.80 Å) PTM is again OTHER throughout and CHILL+ is again clathrate throughout. The cages are 16 of 5^12 and 8 of 5^12 6^4. Neither PTM nor CHILL+ separates those cage types.

A 64-atom ice Ic crystallite in a box with 128 uncorrelated points is the small stand-in for a nucleus. The kernel and PTM both mark 8 core atoms (ice Ic, cubic diamond) and leave the 128 points unlabelled. CHILL+ marks 20 atoms cubic, and also 7 interfacial, 7 inter-clathrate and 1 clathrate, because four neighbours are enough for its bond test. The kernel and PTM ask for the second shell, so they agree with each other and disagree with those extra 12 cubic labels.

## Thermal frames

The frames are mW water (Molinero and Moore, J. Phys. Chem. B 113, 4008, 2009) in LAMMPS 7 Feb 2024, `pair_style sw`, real units, timestep 5 fs. They are not stored in the tree. `scripts/mw_frames.py` rebuilds them.

On the bulk crystals the quench is labelled with this kernel at ceiling 0.06, and that label is the whole cell. A thermal label is a true positive only when it names that same template. CHILL+ maps cubic and reCubic to ice Ic, and hexagonal and reHex to ice Ih. PTM maps cubic diamond to ice Ic and hexagonal diamond to ice Ih. At 260 K and 270 K the kernel's precision is 1 and its recall is at or above PTM. A neighbour just past 3.5 Å is kept when the next atom sits in a gap of at least 0.12 Å, and only for an atom that already has three crystal neighbours. An all-eclipsed shell in that same neighbourhood is searched as a real template instead of being dropped.

| Frame | Kernel precision | Kernel recall | PTM precision | PTM recall |
| --- | ---: | ---: | ---: | ---: |
| Ice Ih, 260 K, 256 | 1 | 1 | 1 | 1 |
| Ice Ih, 270 K, 256 | 1 | 1 | 1 | 0.996 |
| Ice Ic, 260 K, 512 | 1 | 0.998 | 1 | 0.996 |
| Ice Ic, 270 K, 512 | 1 | 0.982 | 1 | 0.961 |

One thread on the previous kernel, median of six calls after a warm-up: the 4096-molecule cell is 16.908 ms (4.13 µs/atom) against PTM 42.885 ms (10.47 µs/atom); the 5750-water sI cell is 34.309 ms (5.97 µs/atom) against PTM 39.783 ms (6.92 µs/atom); the 220 K stack is 4.845 ms (18.9 µs/atom) against PTM 12.023 ms (47.0 µs/atom). Eight threads: 2.246 ms, 4.636 ms and 0.907 ms, against PTM 43.590 ms, 40.536 ms and 12.264 ms. PTM does not speed up from one thread to eight. Raw kernel samples, milliseconds, 4096 then 5750, one thread: 16.8127 16.9228 16.9441 16.8937 17.0251 16.8902 and 34.4555 34.1342 35.4087 34.2692 34.3481 34.2526. Eight threads: 2.21406 2.20469 2.27712 2.30962 2.31932 2.20973 and 4.58636 4.64197 4.54961 4.62909 4.99974 6.17401.

The nucleus is the 40 ps melt at 500 K with the 71-oxygen ice Ic core put back and held for 2 ps at 250 K. A 2 ps time average with the core fixed smears the liquid through the second shell, and the average is not ice. The ground truth used for the table below freezes the liquid at its instantaneous coordinates and quenches only the core with steepest descent, then labels that quench with the earlier kernel. That inherent structure has 12 ice oxygens, all in the core. The instantaneous frame is scored against it. Interface atoms are the non-core oxygens within 3.5 Å of a core oxygen. The table is that earlier kernel, not the production ceiling.

| Region | n | Ground-truth ice | Kernel P | Kernel R | PTM P | PTM R | CHILL+ P | CHILL+ R |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Core | 71 | 12 | 0.700 | 0.583 | 0.467 | 0.583 | 0.667 | 0.667 |
| Interface | 60 | 0 | 1 | 1 | 1 | 1 | 1 | 1 |
| Liquid | 316 | 0 | 1 | 1 | 1 | 1 | 1 | 1 |

The kernel and PTM have the same core recall (7 of 12). The kernel adds 3 false core labels and PTM adds 8. Neither labels the interface or the liquid.

No d-SEAMS ice-nucleation or clathrate dataset is published on Hugging Face from this machine. Three 1728-molecule ice Ic cells were melted, a 71-oxygen core put back, and each cell held at 220 K for 120 ps with that core stationary. The liquid is about 1670 waters. None of the seeds melts. The plot is `markdown/nucleation-120ps.svg`. Ice counts for seeds 11001, 22002 and 33003:

| Time | Kernel | CHILL+ crystal | PTM |
| --- | --- | --- | --- |
| 0 | 6, 4, 11 | 9, 10, 17 | 9, 4, 8 |
| 20 ps | 58, 58, 58 | 42, 45, 57 | 26, 28, 28 |
| 40 ps | 63, 69, 76 | 42, 60, 73 | 32, 34, 30 |
| 60 ps | 66, 74, 90 | 62, 55, 72 | 30, 30, 28 |
| 80 ps | 72, 79, 93 | 64, 60, 77 | 33, 40, 33 |
| 100 ps | 96, 99, 84 | 75, 73, 69 | 32, 40, 31 |
| 120 ps | 112, 92, 89 | 81, 75, 61 | 44, 44, 30 |

Those 1728-molecule counts are the earlier kernel. By 120 ps the nucleus is still a few percent of the cell. The counts below are the production kernel.

The same protocol on a 13824-molecule cell with a 300-oxygen seed, held at 220 K for 200 ps, three seeds (91001, 92002, 93003). About 13660 waters remain after the core is put back. The seed grows in all three and does not melt. Counts are the three seeds in that order.

| Time | Kernel | CHILL+ crystal | PTM |
| --- | --- | --- | --- |
| 0 | 92, 90, 82 | 156, 146, 140 | 101, 99, 102 |
| 40 ps | 211, 214, 213 | 374, 370, 369 | 198, 198, 202 |
| 80 ps | 243, 235, 243 | 410, 401, 406 | 231, 218, 226 |
| 120 ps | 291, 283, 287 | 493, 446, 465 | 269, 263, 265 |
| 160 ps | 324, 296, 325 | 518, 497, 467 | 294, 275, 296 |
| 200 ps | 335, 316, 353 | 544, 488, 492 | 306, 295, 329 |

On the 200 ps frames, one thread, the kernel takes 85.827, 87.216 and 75.212 ms (6.28, 6.39 and 5.51 µs/atom). Raw samples, milliseconds: 84.541 85.338 85.448 86.205 87.365 88.169, then 85.529 86.750 86.925 87.506 87.605 87.662, then 74.151 74.531 75.181 75.243 75.373 75.901. PTM on the same three frames, a fresh pipeline each sample and only `compute` timed, is 89.434, 91.753 and 90.721 ms (6.54, 6.72 and 6.65 µs/atom). Raw: 87.451 87.921 89.254 89.613 90.041 92.565, then 88.598 90.707 91.028 92.478 93.142 93.679, then 88.666 89.514 89.940 91.501 92.141 104.149. One thread is under PTM on all three seeds.

The rebuilt 200 ps cells are 13553 molecules: 12×12×12 ice Ic melted 40 ps at 500 K, the 300 central molecules put back, liquid inside 2.4 Å of that core removed, the core held fixed, and the rest run at 220 K. The first test is the Appendix A fingerprint of Maras, Trushin, Stukowski, Albe and Nordlund, Comput. Phys. Commun. 205, 13 (2016). The neighbour tables are allocated once and reused. A first shell whose bonds lie outside 2.30–3.35 Å does not open the second shell. A shell with q at least 0.95, bonds within 15%, and a second shell within a factor 1.6 is passed to the permutation fit on the twelve neighbours already in hand. Callgrind on one disordered frame puts 69% of the instructions in the permutation fit, 14% in the neighbour build, 13% in the second-shell assembly and the fingerprint, and 4% in the hexagonal recovery. The neighbour buffers are reused, so allocation does not appear. A tetrahedron with q below 0.90 does not open the second shell. The fit runs only for a shell that already touches a labelled crystal. Ideal ice, equal bonds and one second-shell radius, is labelled by the raw cosine of each second-shell bond with the other first-shell bonds, 1/3 when staggered and about 0.78 when eclipsed, and does not build the fingerprint. The neighbour search uses 3.35 Å cells and extends a short list out to 3.5 Å. One thread on these frames takes 10.198, 9.604 and 10.331 ms. Raw: 10.174 10.189 10.190 10.205 10.210 10.428, then 9.549 9.592 9.602 9.606 9.687 9.746, then 10.155 10.295 10.329 10.333 10.337 10.373. OVITO identify diamond in that session is 5.847, 5.243 and 5.903 ms. pyscal3 is 101.142, 102.530 and 101.895 ms. PTM is 34.830, 34.826 and 34.299 ms. Eight threads: kernel 4.634, 4.234 and 5.556 ms; identify diamond 6.506, 6.418 and 5.487 ms; pyscal3 23.023, 22.872 and 21.448 ms; PTM 38.706, 37.960 and 36.327 ms. On ideal ice Ic, 110592 molecules, one thread is 25.241 ms against identify diamond at 36.973 ms, and every oxygen is ice Ic. Raw kernel: 25.060 25.115 25.216 25.265 25.366 25.995. Raw identify diamond: 30.789 35.367 36.672 37.274 37.819 38.157. The slowest kernel sample is under the fastest identify-diamond sample. Eight threads on that cell: kernel 7.478 ms, identify diamond 36.662 ms, pyscal3 26.653 ms, PTM 206.688 ms. Raw lines are `benchmarks/logs/symmetry-v17-session.txt`. Raw lines are `benchmarks/logs/symmetry-v17-session.txt`.

The earlier table labelled the quenched frames with the kernel alone at ceiling 0.06. That truth favoured the kernel: its recall sat near 0.99 and its precision near 0.42, because the kernel's own loose labels were the answer. The quenches used below are labelled by each method with its own production rule. The kernel's rule is ceiling 0.12, the q gate, the eclipsed reject, and one crystal neighbour. CHILL+ and PTM are not retuned.

Two definitions use the 200 ps quench. Majority is the template named by at least two of the three methods. Unanimous is the template all three name. The third definition uses the quenches of the 190 ps frame and the 200 ps frame: an atom counts only when the majority template is the same at both times. A thermal label is a true positive only when it names that template.

| Seed | Truth | n | Kernel P | Kernel R | Kernel F1 | CHILL+ P | CHILL+ R | PTM P | PTM R | PTM F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 91001 | majority | 382 | 0.988 | 0.866 | 0.923 | 0.667 | 0.950 | 0.997 | 0.798 | 0.887 |
| 91001 | unanimous | 345 | 0.967 | 0.939 | 0.953 | 0.616 | 0.971 | 0.987 | 0.875 | 0.928 |
| 91001 | 190 and 200 ps | 353 | 0.955 | 0.907 | 0.930 | 0.625 | 0.963 | 0.977 | 0.847 | 0.907 |
| 92002 | majority | 350 | 1 | 0.903 | 0.949 | 0.686 | 0.957 | 0.993 | 0.837 | 0.909 |
| 92002 | unanimous | 320 | 0.953 | 0.941 | 0.947 | 0.633 | 0.966 | 0.976 | 0.900 | 0.937 |
| 92002 | 190 and 200 ps | 330 | 0.965 | 0.924 | 0.944 | 0.656 | 0.970 | 0.980 | 0.876 | 0.925 |
| 93003 | majority | 391 | 0.983 | 0.887 | 0.933 | 0.750 | 0.944 | 0.964 | 0.811 | 0.881 |
| 93003 | unanimous | 351 | 0.932 | 0.937 | 0.935 | 0.687 | 0.963 | 0.933 | 0.875 | 0.903 |
| 93003 | 190 and 200 ps | 366 | 0.955 | 0.921 | 0.937 | 0.709 | 0.954 | 0.945 | 0.850 | 0.895 |

Seeds 91001 and 92002 are the training pair. Pooled against the majority, PTM is precision 0.995 and recall 0.817. The ceiling sweep, with one crystal neighbour kept on, is:

| Ceiling | Precision | Recall | F1 |
| ---: | ---: | ---: | ---: |
| 0.08 | 1 | 0.675 | 0.806 |
| 0.10 | 0.998 | 0.796 | 0.886 |
| 0.12 | 0.994 | 0.883 | 0.935 |
| 0.14 | 0.960 | 0.947 | 0.953 |
| 0.16 | 0.880 | 0.960 | 0.918 |
| 0.18 | 0.750 | 0.967 | 0.845 |

0.12 is the tightest ceiling whose recall stays above PTM. Its precision is 0.994, just under PTM's 0.995. Requiring two crystal neighbours instead of one puts the training majority at precision 0.992 and recall 0.858, which is above PTM on that definition, and puts the held-out recall under PTM, so the production rule stays at one neighbour.

Seed 93003 was held out. The q level for the fit was chosen on seeds 91001 and 92002. A partial cubic fingerprint with fewer than eight (4,2,1) signatures is not labelled directly; that shell is passed to the fit, and the fit rejects the liquid false positives of that class. The held-out majority is precision 0.982 and recall 0.839 (328 true positives, 6 false positives), against PTM at 0.964 and 0.811. The labels that persist from 190 ps to 200 ps are precision 0.955 and recall 0.872, against PTM at 0.945 and 0.850. Precision and recall on both are above PTM.

pyscal3 4.1.0 diamond structure and OVITO 3.17.2 identify diamond were scored on the same three truths. Only cubic diamond and hexagonal diamond count; the first-neighbour and second-neighbour tags are not crystal. The two programs return the same labels. On the held-out majority their precision is 0.967 and their recall is 0.742, F1 0.839, against the kernel at 0.933. They miss more of the nucleus than PTM does.

| Seed | Truth | pyscal/OVITO diamond P | R | F1 | Kernel F1 | PTM F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 91001 | majority | 0.986 | 0.746 | 0.849 | 0.923 | 0.887 |
| 91001 | unanimous | 0.976 | 0.817 | 0.890 | 0.953 | 0.928 |
| 91001 | 190 and 200 ps | 0.969 | 0.793 | 0.872 | 0.930 | 0.907 |
| 92002 | majority | 0.992 | 0.751 | 0.855 | 0.949 | 0.909 |
| 92002 | unanimous | 0.977 | 0.809 | 0.885 | 0.947 | 0.937 |
| 92002 | 190 and 200 ps | 0.966 | 0.776 | 0.861 | 0.944 | 0.925 |
| 93003 | majority | 0.967 | 0.742 | 0.839 | 0.933 | 0.881 |
| 93003 | unanimous | 0.950 | 0.812 | 0.876 | 0.935 | 0.903 |
| 93003 | 190 and 200 ps | 0.953 | 0.781 | 0.859 | 0.937 | 0.895 |

The kernel's false positives against the majority are interfacial, not a second shell of the wrong atoms. On 91001 there are 4, all within 3.5 Å of a majority-crystal oxygen, and CHILL+ calls them cubic or hexagonal. On 92002 there are none. On 93003 there are 6, all within 3.5 Å of a majority-crystal oxygen. Their scale-free residual sits higher than the true positives: of the 4 on 91001, none pass a ceiling of 0.08 and only 1 passes 0.10, while 298 of 330 true positives pass 0.10. The true-positive q has 10th percentile 0.97 and median 1.00. The liquid that the kernel leaves unlabelled has 10th percentile 0.65 and median 0.87, so the q gate removes disordered liquid and does not separate the interfacial false positives from the crystal. The twelve topological second neighbours of a true positive are the twelve geometric neighbours past 3.2 Å (Jaccard median 1, radius median 4.49 Å, the ideal second-neighbour distance at a 2.75 Å bond). The false positives match that shell too. Liquid's Jaccard median is about 0.65, and only about half of those atoms even have twelve topological second neighbours. The false positive is a crystal-like shell at the surface of the nucleus, admitted because its residual is between 0.10 and 0.12.

Guest-free sI (2×2×2, 368 waters) and sII (2×2×2, 1088 waters) were run for 10 ps at 220 K and at 250 K and quenched the same way. Cages come from the primitive ring network at a 3.5 Å bond cutoff and `cage::findBySignature`, kept only when the shape residual of the cage against its template is under 1.5 Å. At both temperatures the thermal sI cell has 16 of 5^12 and 48 of 5^12 6^2, and every cage matches an inherent-structure cage of the same kind within 2 Å. The thermal sII cell has 128 of 5^12 and 64 of 5^12 6^4, again with every cage matched. Widening the bond cutoff to 3.7 Å drops 4 small cages and 3 large cages from the 250 K sII cell, so 3.5 Å is the cutoff that agrees with the quench. The ice classifier names none of the 5750-water sI cell as ice. PTM names that cell other throughout. CHILL+ names it clathrate.

## Paragraph

Local ice is labelled by the scale-free residual of the seventeen-atom shell to ice Ic, ice Ih or stacking-disordered ice. A quaternion overlay of the tetrahedron, with the second shell ordered by azimuth, supplies that residual after the optimal rotation and scale. The atoms are split across threads. On ideal ice Ic this is faster, per atom, than polyhedral template matching in OVITO 3.17.2 from 216 molecules through 110592, on one thread and on eight. A general IRA overlay of the same shells is kept as a reference and is about a thousand times slower per atom. The overlay that is returned for a named shell, cage or ring uses a 3×3 SVD and a permutation fixed by the bond graph, and the continuous shape measure is the scaled residual of that overlay. The same machinery names the cages of the sI and sII lattices, two and six in sI and sixteen and eight in sII, which template matching reports only as other and which CHILL+ reports only as clathrate. On a stacking-fault cell, template matching returns other for every oxygen, while the shell match splits the cell evenly between ice Ic and ice Isd.

The label kept for production is the scale-free residual of the seventeen-atom shell, ceiling 0.12, and a crystal atom must have one crystal neighbour. The neighbour list is the orthorhombic cutoff CSR: cells no wider than the 3.5 Å bond, and one thread walks each pair once. A tetrahedron with q below 0.94 does not open the second shell. A clean second shell is labelled by the common-neighbour test before the fit. The fit is a 4-neighbour Horn quaternion, and the bond-angle fingerprint names the one template it tries; a near miss on the hexagonal recovery still sees the other templates. On the held-out seed 93003 the majority recall is 0.839, above PTM at 0.811, and the precision is 0.982, above PTM at 0.964. The persistent recall is 0.872, above PTM at 0.850, and the precision is 0.955, above PTM at 0.945. The Horn fit uses two Jacobi sweeps. A deferred shell that is still far after four correspondences is not crystal, and a shell with one crystal neighbour stops after six. Hexagonal recovery is one pass. On the rebuilt 13553-molecule frames one thread takes 5.672, 5.287 and 5.472 ms, against identify diamond at 5.703, 5.427 and 5.273 ms in the same session. Eight threads take 2.277, 2.104 and 2.133 ms, against 6.246, 5.532 and 5.033 ms. Ideal ice at 110592 molecules is 22.953 ms against 36.453 ms on one thread and 6.604 ms against 34.603 ms on eight, and every oxygen is cubic. The raw samples are in `benchmarks/logs/symmetry-v19-session.txt`. `scripts/symmetry_frame_labels.py` writes `index,id,kind` for `input/traj/nucleation_f1.lammpstrj`, `input/traj/interface_f0.lammpstrj`, and the three 200 ps seeds. The two in-tree frames are in `benchmarks/labels/`. pyscal3 and the OVITO diamond identifier agree with each other and recall less of the nucleus than the kernel or PTM. A 3.5 Å ring network plus a shape ceiling of 1.5 Å names every thermal sI and sII cage and matches each one to the inherent structure. One thread names a hexagon in 3.80 µs and a dodecahedron in 4.82 µs, below IRA.

## Rerun

The machine is an Intel Xeon, 8 CPUs, affinity 0–7, GCC 15.2.0. The build is Meson `release` (`-O3`, and `-mavx2 -mfma` on the shell translation unit). One warm-up, then six calls. The median of an even sample is the mean of the two central values.

```
meson setup build --buildtype=release
meson compile -C build tests/bench_shell_scale tests/test_local_symmetry
./build/tests/bench_shell_scale 1
./build/tests/bench_shell_scale 8
```

The head-to-head builds the same 24×24×24 ice Ic cell (110592 molecules) and times `symmetry.classify`, `pyscal3.diamond_structure` and `IdentifyDiamondModifier` in one process at `OMP_NUM_THREADS=1`, then again at 8. One thread, milliseconds: kernel 92.652 93.212 93.530 96.134 96.333 101.512; pyscal3 98.769 98.941 100.183 100.297 100.678 101.754; OVITO identify diamond 29.440 30.122 35.399 36.767 37.263 37.803. Eight threads: kernel 20.567 20.905 21.426 22.350 23.332 23.641; pyscal3 23.719 24.013 24.026 24.102 24.168 26.548; OVITO identify diamond 30.363 32.463 37.536 37.901 40.083 43.576. On one thread the kernel is ahead of pyscal3 and behind OVITO's diamond identifier. On eight threads it is ahead of both. Every oxygen is cubic diamond for all three.

## Citations

- Gunde, Salles, Hémeryck, Martin-Samos. IRA: A shape matching approach for recognition and comparison of generic atomic patterns. J. Chem. Inf. Model. 61, 5446 (2021). doi:10.1021/acs.jcim.1c00567
- Gunde, Salles, Grisanti, Martin-Samos, Hémeryck. SOFI: Finding point group symmetries in atomic clusters as finding the set of degenerate solutions in a shape-matching problem. J. Chem. Phys. 161, 062503 (2024). doi:10.1063/5.0215689
- Togo, Shinohara, Tanaka. Spglib: a software library for crystal symmetry search. Sci. Technol. Adv. Mater. Methods 4, 2384822 (2024). doi:10.1080/27660400.2024.2384822
- Shinohara. moyo: A fast and robust crystal symmetry finder, written in Rust. figshare (2026). doi:10.6084/m9.figshare.31081162
- Ong, Richards, Jain, Hautier, Kocher, Cholia, Gunter, Chevrier, Persson, Ceder. Python Materials Genomics (pymatgen). Comput. Mater. Sci. 68, 314 (2013). doi:10.1016/j.commatsci.2012.10.028
- Pinsky, Avnir. Continuous symmetry measures. 5. The classical polyhedra. Inorg. Chem. 37, 5575 (1998). doi:10.1021/ic9804925
- Alvarez, Alemany, Casanova, Cirera, Llunell, Avnir. Shape maps and polyhedral interconversion paths in transition metal chemistry. Coord. Chem. Rev. 249, 1693 (2005). doi:10.1016/j.ccr.2005.03.031
- Coutsias, Seok, Dill. Using quaternions to calculate RMSD. J. Comput. Chem. 25, 1849 (2004). doi:10.1002/jcc.20110
- Larsen, Schmidt, Schiøtz. Robust structural identification via polyhedral template matching. Modell. Simul. Mater. Sci. Eng. 24, 055007 (2016). doi:10.1088/0965-0393/24/5/055007
- Nguyen, Molinero. Identification of clathrate hydrates, hexagonal ice, cubic ice, and liquid water in simulations: the CHILL+ algorithm. J. Phys. Chem. B 119, 9369 (2015). doi:10.1021/jp510289t
- Matsumoto, Yagasaki, Tanaka. GenIce: Hydrogen-disordered ice generator. J. Comput. Chem. 39, 61 (2018). doi:10.1002/jcc.25077
- Matsumoto, Yagasaki, Tanaka. GenIce-core. J. Chem. Phys. 160, 094101 (2024). doi:10.1063/5.0198056
- Molinero, Moore. Water modeled as an intermediate element between carbon and silicon. J. Phys. Chem. B 113, 4008 (2009). doi:10.1021/jp805227c
