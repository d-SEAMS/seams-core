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

The production label of a bulk oxygen is not a general permutation search. A periodic cell list gathers the sixteen nearest neighbours inside 5 Å. The four nearest of those are the tetrahedral shell. Each bond is eclipsed or staggered according to the azimuthal alignment of the outer neighbours with the other three bonds of the centre. The ideal eclipsed cosine is 1 and the ideal staggered cosine is 0.5. The default cutoff, 0.75, sits between them. It is an argument of `classifyIce` and of `seams_symmetry_classify_tol`. Zero eclipsed bonds is ice Ic, one is ice Ih, and two is ice Isd. Anything else is left unlabelled by this pass. On the ideal cells the label matches a full IRA overlay of the shell. Raising the cutoff above 1 drops the hexagonal bond, so ideal ice Ih is labelled ice Ic. Lowering it below 0.5 counts the staggered bonds, so ideal ice Ic is left unlabelled.

A named shell, cage or ring still returns a rotation, a permutation and an RMSD. That alignment is a 3×3 SVD of the covariance (Kabsch), tried proper and improper, on a permutation taken from the bond graph rather than from N!. Rings use the cyclic order. Ice uses the 24 assignments of the first shell and then the bonded second neighbours. Cages grow a degree-3 matching from each vertex. The environments are independent, so `matchAll` runs them in parallel. The continuous shape measure is the scaled residual of that overlay, 100 times the mean square distance to the template divided by the mean square radius of the template (Pinsky and Avnir 1998; Alvarez et al. 2005). It is computed in-tree. The GPL-2 shape program is a separate executable, built when OpenBabel is present, and is not on this path.

## Speed

Release `-O3`, one warm-up, then six calls. The median of an even sample is the mean of the two central values. Quartiles are linear (type 7). Every oxygen of ideal ice Ic is labelled ice Ic. OVITO 3.17.2 PTM has cubic and hexagonal diamond enabled. It does not speed up when `OMP_NUM_THREADS` goes from 1 to 8. The kernel does.

| Molecules | Threads | Kernel median (ms) | Kernel per atom (µs) | PTM median (ms) | PTM per atom (µs) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 216 | 1 | 0.256 | 1.18 | 8.447 | 39.11 |
| 216 | 8 | 0.050 | 0.230 | 8.572 | 39.69 |
| 10648 | 1 | 14.176 | 1.33 | 72.986 | 6.85 |
| 10648 | 8 | 2.009 | 0.189 | 73.003 | 6.86 |
| 110592 | 1 | 130.380 | 1.18 | 545.459 | 4.93 |
| 110592 | 8 | 18.728 | 0.169 | 539.262 | 4.88 |

Raw kernel samples, 1 thread, milliseconds: 216 is 0.246768 0.254731 0.255038 0.256716 0.260282 0.272244; 10648 is 14.0435 14.1211 14.1751 14.1775 14.2388 14.2674; 110592 is 129.171 130.027 130.162 130.599 130.698 140.118. Eight threads: 216 is 0.041396 0.047861 0.049413 0.049966 0.050533 0.054676; 10648 is 1.94715 1.96203 2.00624 2.0112 2.07868 2.10591; 110592 is 18.5042 18.6201 18.7258 18.7293 18.7547 19.6979.

IRA on the same 216-atom cell, one shell per oxygen through libira, is serial. Median 283.84 ms (raw 282.74 282.963 283.56 284.119 285.766 286.336), about 1.31 ms per atom, three orders of magnitude above the kernel. spglib on the 8-atom ice Ic cell returns 227 `Fd-3m` in a median 6.13 ms. That is one cell, not a per-atom cost.

The naming path is the SVD and the bond-graph permutation, timed on rotated copies of the templates rather than on a neighbour search. IRA is the same shells through libira, serial. PTM does not name a cage or a ring. The PTM number below is structure identification of a 216-atom ice Ic crystal, from the same protocol.

| Path | Shells | Threads | Median (ms) | Per shell (µs) |
| --- | ---: | ---: | ---: | ---: |
| Kernel, ice Ic | 216 | 1 | 35.215 | 163 |
| Kernel, ice Ic | 216 | 8 | 5.290 | 24.5 |
| Kernel, ice Ic | 2048 | 1 | 334.85 | 164 |
| Kernel, ice Ic | 2048 | 8 | 44.712 | 21.8 |
| IRA, ice Ic | 216 | 1 | 159.26 | 737 |
| Kernel, 5^12 | 216 | 1 | 150.65 | 697 |
| Kernel, 5^12 | 216 | 8 | 19.502 | 90.3 |
| IRA, 5^12 | 64 | 1 | 22.427 | 350 |
| Kernel, 5^12 6^2 | 216 | 1 | 66.773 | 309 |
| Kernel, 5^12 6^2 | 216 | 8 | 9.177 | 42.5 |
| Kernel, ring-6 | 216 | 1 | 4.549 | 21.1 |
| Kernel, ring-6 | 216 | 8 | 0.588 | 2.72 |
| Kernel, ring-6 | 2048 | 1 | 43.738 | 21.4 |
| Kernel, ring-6 | 2048 | 8 | 5.576 | 2.72 |
| IRA, ring-6 | 216 | 1 | 3.130 | 14.5 |
| In-tree CShM, ice Ic | 216 | 1 | 35.197 | 163 |
| In-tree CShM, ice Ic | 216 | 8 | 4.479 | 20.7 |
| Shape program, one ice shell, identity permutation | 1 | 1 | 9.884 | 9884 |
| Shape program, tetrahedron, all 24 permutations | 1 | 1 | 10.360 | 10360 |
| PTM structure id, ice Ic crystal | 216 atoms | 1 | 8.573 | 39.7 per atom |

On one thread a ring is slightly slower than IRA, and a dodecahedron is slower than IRA. Eight threads reverse both. The shape-program row is one process, OpenBabel startup included. The in-tree measure on the same distorted ice shell agrees with that program to 0.05. The program printed `SHAPE = 0.043800` for a shell with one oxygen moved by a few tenths of an angstrom. Raw kernel samples, ice Ic, 1 thread, milliseconds: 216 shells 35.1267 35.1646 35.1953 35.2352 35.4733 37.1024; 2048 shells 333.275 333.691 334.234 335.466 335.909 336.319. Eight threads: 216 shells 5.0707 5.16434 5.1825 5.39789 5.43898 5.46815; 2048 shells 42.8295 43.3195 43.9354 45.4883 47.708 47.8798. IRA ice, 216 shells, 1 thread: 158.408 159.203 159.251 159.268 160.189 163.331. In-tree CShM, 216 shells, 1 thread: 35.0172 35.1007 35.1804 35.214 35.275 35.9081. Eight threads: 4.42846 4.44697 4.45494 4.50364 4.54158 4.61616. Shape program, ice identity: 9.768 9.851 9.882 9.886 9.970 10.429. Tetrahedron: 9.886 10.043 10.251 10.469 10.508 10.929. PTM, 216 atoms: 8.515148 8.631621 8.632679 8.681955 8.433761 8.488970.

## What PTM does not name

On the stacking-fault cell, PTM returns OTHER for all 8 oxygens. The IRA shell match returns ice Ic for 4 and ice Isd for 4. CHILL+ on those 8 oxygens taken out of the periodic cell, so the crystal neighbours are missing, returns water for all 8. It has no stacking-fault label.

On the sI oxygen lattice (46 waters, cell 12.239 Å) PTM returns OTHER for every water and CHILL+ returns clathrate for every water. The ring census finds 48 five-rings and 54 six-rings, and the closed cages are 2 of 5^12 and 6 of 5^12 6^2. On the sII lattice (136 waters, nearest neighbour scaled to 2.80 Å) PTM is again OTHER throughout and CHILL+ is again clathrate throughout. The cages are 16 of 5^12 and 8 of 5^12 6^4. Neither PTM nor CHILL+ separates those cage types.

A 64-atom ice Ic crystallite in a box with 128 uncorrelated points is the small stand-in for a nucleus. The kernel and PTM both mark 8 core atoms (ice Ic, cubic diamond) and leave the 128 points unlabelled. CHILL+ marks 20 atoms cubic, and also 7 interfacial, 7 inter-clathrate and 1 clathrate, because four neighbours are enough for its bond test. The kernel and PTM ask for the second shell, so they agree with each other and disagree with those extra 12 cubic labels.

## Thermal frames

The frames are mW water (Molinero and Moore, J. Phys. Chem. B 113, 4008, 2009) in LAMMPS 7 Feb 2024, `pair_style sw`, real units, timestep 5 fs, after a conjugate-gradient minimisation. They are not stored in the tree. `scripts/mw_frames.py` rebuilds them.

The stacking cell is 4×4 copies of the hexagonal basal cell and eight bilayers (256 oxygens). The site sequence is ABABCACB, with A = (1/3, 2/3), B = (2/3, 1/3) and C = (0, 0). A step between two different sites keeps the oxygen–oxygen bond at 2.75 Å. The ideal cell, before any dynamics, is 128 ice Ic and 128 ice Ih, and no ice Isd. The layer string, one token per bilayer of 32 oxygens, is `hhhcchcc`. CHILL+ returns 128 cubic and 128 hexagonal. PTM, with cubic and hexagonal diamond enabled, returns the same 128 and 128 at every RMSD cutoff from 0.05 to 0.30. All three agree on every oxygen. The minimised energy is −12.29 kcal/mol per molecule.

Each thermal run is 20 ps of NVT. At 220 K the kernel reads 55 ice Ic, 97 ice Ih, 30 ice Isd and 74 unlabelled. CHILL+ is still 128 cubic and 128 hexagonal. PTM at cutoff 0.10 is the same 128 and 128, and it stays there through 0.30. At cutoff 0.05 it keeps one hexagonal oxygen and drops the rest. At 270 K the kernel reads 47 ice Ic, 82 ice Ih, 29 ice Isd and 98 unlabelled. CHILL+ is 128 cubic, 126 hexagonal and 2 interfacial. PTM at 0.10 is 128 cubic and 125 hexagonal.

The nucleus is a 4×4×4 ice Ic cell (512 oxygens). Oxygens more than 8 Å from the centre are held at 360 K for 10 ps while the core is held at 230 K, and the whole cell is then held at 250 K for 10 ps. The kernel labels 64 of the 71 core oxygens as ice and 177 of the 441 outer oxygens as ice (154 ice Ic, 61 ice Ih, 26 ice Isd, 271 unlabelled in the whole cell). CHILL+ labels all 512 cubic. PTM at cutoff 0.10 also labels all 512 cubic, and at 0.05 it keeps one.

A contested oxygen is one where the kernel, CHILL+ and PTM do not all name the same family. The second shell decides it. If the bond-graph residual to one of the three templates is below 0.45 Å, that template is the shell. Otherwise the four bond cosines decide, and only when every cosine lies outside [0.65, 0.85], so the 0.75 cut is not what the label hinges on. On the ideal cell all 256 shells lock and all three methods match them. At 220 K, 112 shells lock: all three methods match 61 of them, and CHILL+ and PTM match the shell on a further 51 where the eclipsed count has already moved. The other 144 oxygens have a cosine inside that window, so the cutoff itself is the dispute and no hard label is uniquely right. At 270 K, 65 shells lock (36 agreed by all three, 27 by CHILL+ and PTM alone) and 191 are inside the window. In the nucleus, 239 shells lock: 85 agreed by all three, 153 by CHILL+ and PTM where the kernel has moved, and 273 inside the window.

Moving the kernel cutoff off 0.75 changes the 220 K labels by 143 oxygens at 0.65 and by 56 at 0.85. The same frame's PTM labels do not change between 0.10 and 0.30, and they collapse at 0.05. The default cutoff stays at 0.75, between the ideal values 0.5 and 1. The argument is there because a 20 ps mW frame puts many cosines into the gap.

## Paragraph

Local ice is labelled from the eclipsed-bond count of the tetrahedral shell, which is the same geometric distinction as the seventeen-atom templates of ice Ic, ice Ih and stacking-disordered ice. The count is evaluated from a periodic cell list, with no heap traffic in the atom loop and with the atoms split across threads. On ideal ice Ic this is faster, per atom, than polyhedral template matching in OVITO 3.17.2 from 216 molecules through 110592, on one thread and on eight. A general IRA overlay of the same shells is kept as a reference and is about a thousand times slower per atom. The overlay that is returned for a named shell, cage or ring uses a 3×3 SVD and a permutation fixed by the bond graph, and the continuous shape measure is the scaled residual of that overlay. The same machinery names the cages of the sI and sII lattices, two and six in sI and sixteen and eight in sII, which template matching reports only as other and which CHILL+ reports only as clathrate. On a stacking-fault cell, template matching returns other for every oxygen, while the shell match splits the cell evenly between ice Ic and ice Isd.

On a 256-oxygen stacking-disordered mW cell the same count, CHILL+ and PTM agree on every oxygen before any motion: 128 ice Ic and 128 ice Ih, in the layer order `hhhcchcc`, and no ice Isd. After 20 ps at 220 K and at 270 K, PTM at its default RMSD cutoff of 0.1 and CHILL+ still return that ideal split, while the eclipsed count moves, because many bond cosines now lie between the ideal staggered value 0.5 and the ideal eclipsed value 1. Where the second-shell overlay still locks onto one template, CHILL+ and PTM are the ones that match it; where a cosine sits in [0.65, 0.85], the cutoff is the whole of the dispute. The cutoff is an argument and stays at 0.75. Naming the same shells, cages and rings with the SVD and the bond-graph permutation is parallel. At eight threads it is faster per shell than a serial IRA overlay for ice, for both cages and for the hexagon. The in-tree shape measure on a distorted ice shell agrees with the GPL shape program to 0.05, and the program is a separate executable.

## Citations

- Gunde, Salles, Hémeryck, Martin-Samos. IRA: A shape matching approach for recognition and comparison of generic atomic patterns. J. Chem. Inf. Model. 61, 5446 (2021). doi:10.1021/acs.jcim.1c00567
- Gunde, Salles, Grisanti, Martin-Samos, Hémeryck. SOFI: Finding point group symmetries in atomic clusters as finding the set of degenerate solutions in a shape-matching problem. J. Chem. Phys. 161, 062503 (2024). doi:10.1063/5.0215689
- Togo, Shinohara, Tanaka. Spglib: a software library for crystal symmetry search. Sci. Technol. Adv. Mater. Methods 4, 2384822 (2024). doi:10.1080/27660400.2024.2384822
- Shinohara. moyo: A fast and robust crystal symmetry finder, written in Rust. figshare (2026). doi:10.6084/m9.figshare.31081162
- Ong, Richards, Jain, Hautier, Kocher, Cholia, Gunter, Chevrier, Persson, Ceder. Python Materials Genomics (pymatgen). Comput. Mater. Sci. 68, 314 (2013). doi:10.1016/j.commatsci.2012.10.028
- Pinsky, Avnir. Continuous symmetry measures. 5. The classical polyhedra. Inorg. Chem. 37, 5575 (1998). doi:10.1021/ic9804925
- Alvarez, Alemany, Casanova, Cirera, Llunell, Avnir. Shape maps and polyhedral interconversion paths in transition metal chemistry. Coord. Chem. Rev. 249, 1693 (2005). doi:10.1016/j.ccr.2005.03.031
- Larsen, Schmidt, Schiøtz. Robust structural identification via polyhedral template matching. Modell. Simul. Mater. Sci. Eng. 24, 055007 (2016). doi:10.1088/0965-0393/24/5/055007
- Nguyen, Molinero. Identification of clathrate hydrates, hexagonal ice, cubic ice, and liquid water in simulations: the CHILL+ algorithm. J. Phys. Chem. B 119, 9369 (2015). doi:10.1021/jp510289t
- Matsumoto, Yagasaki, Tanaka. GenIce: Hydrogen-disordered ice generator. J. Comput. Chem. 39, 61 (2018). doi:10.1002/jcc.25077
- Matsumoto, Yagasaki, Tanaka. GenIce-core. J. Chem. Phys. 160, 094101 (2024). doi:10.1063/5.0198056
- Molinero, Moore. Water modeled as an intermediate element between carbon and silicon. J. Phys. Chem. B 113, 4008 (2009). doi:10.1021/jp805227c
