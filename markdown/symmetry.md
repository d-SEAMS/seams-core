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

The production label of a bulk oxygen is the ice template with the lowest residual on the sixteen neighbours inside 5 Å. A periodic cell list gathers those neighbours. On a skewed cell the stencil reaches far enough to include every neighbour the shell search would, which a one-cell stencil misses. The four nearest bonds inside 2.35–3.25 Å are the tetrahedron. Each of the 24 assignments of those bonds is laid on the template by a quaternion Kabsch (Coutsias, Seok and Dill 2004; Theobald 2005), and the second neighbours are scored only against the three atoms bonded to the matching first neighbour. If that score is already under 0.55 Å the atom is labelled and the search stops. Otherwise the same shell is passed to the full bond-graph overlay used as the reference. The argument of `classifyIce` is that 0.55 Å ceiling. Zero residual to ice Ic, ice Ih or ice Isd is the label. A shell past the ceiling is left unlabelled.

A named shell, cage or ring still returns a rotation, a permutation and an RMSD. That alignment is a 3×3 SVD of the covariance (Kabsch), tried proper and improper, on a permutation taken from the bond graph rather than from N!. Rings use the cyclic order. Ice uses the 24 assignments of the first shell and then the bonded second neighbours. Cages grow a degree-3 matching from each vertex. The environments are independent, so `matchAll` runs them in parallel. The continuous shape measure is the scaled residual of that overlay, 100 times the mean square distance to the template divided by the mean square radius of the template (Pinsky and Avnir 1998; Alvarez et al. 2005). It is computed in-tree. The GPL-2 shape program is a separate executable, built when OpenBabel is present, and is not on this path.

## Speed

Release `-O3`, one warm-up, then six calls. The median of an even sample is the mean of the two central values. Quartiles are linear (type 7). Every oxygen of ideal ice Ic is labelled ice Ic. OVITO 3.17.2 PTM has cubic and hexagonal diamond enabled. It does not speed up when `OMP_NUM_THREADS` goes from 1 to 8. The kernel does.

| Molecules | Threads | Kernel median (ms) | Kernel per atom (µs) | PTM median (ms) | PTM per atom (µs) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 216 | 1 | 0.371 | 1.72 | 8.447 | 39.11 |
| 216 | 8 | 0.065 | 0.301 | 8.572 | 39.69 |
| 10648 | 1 | 25.586 | 2.40 | 72.986 | 6.85 |
| 10648 | 8 | 2.670 | 0.251 | 73.003 | 6.86 |
| 110592 | 1 | 174.431 | 1.58 | 545.459 | 4.93 |
| 110592 | 8 | 26.954 | 0.244 | 539.262 | 4.88 |

Raw kernel samples, 1 thread, milliseconds: 216 is 0.348139 0.353123 0.369240 0.373251 0.373592 0.397438; 10648 is 17.9323 23.7384 25.5109 25.6612 25.7260 25.9182; 110592 is 173.079 173.774 173.782 175.080 177.467 188.816. Eight threads: 216 is 0.049742 0.054405 0.057775 0.059024 0.062307 0.075262; 10648 is 2.53141 2.64807 2.65663 2.68348 2.85126 2.94446; 110592 is 24.1951 24.7195 27.4514 32.5204 32.605 32.7736. On a skewed thermal cell the first overlay often stays above 0.55 Å, so those atoms take the full bond-graph search. The 256-oxygen 220 K frame then takes 55.5 ms. The ideal cells above stay on the fast path.

IRA on the same 216-atom cell, one shell per oxygen through libira, is serial. Median 283.84 ms (raw 282.74 282.963 283.56 284.119 285.766 286.336), about 1.31 ms per atom, three orders of magnitude above the kernel. spglib on the 8-atom ice Ic cell returns 227 `Fd-3m` in a median 6.13 ms. That is one cell, not a per-atom cost.

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

The stacking cell is 4×4 copies of the hexagonal basal cell and eight bilayers (256 oxygens), site sequence ABABCACB. Before dynamics it is 128 ice Ic and 128 ice Ih, layer order `hhhcchcc`. CHILL+ and PTM agree on every oxygen. A shell is locked when the full bond-graph residual to one template is below 0.45 Å. That residual is the ground truth for the rows below. A false label is an ice name on a shell that does not lock.

| Decision | 220 K locked (of 94) | 270 K locked (of 57) | False labels, 220 K | Ideal ice Ic, 1 thread (µs/atom) |
| --- | ---: | ---: | ---: | ---: |
| Eclipsed count, cutoff 0.75 | 68 | 43 |  | 1.18 |
| Soft distance to the ideal cosine pattern | 68 | 43 |  | 1.2 |
| Quaternion fit of the tetrahedron, no full search | 47 | 26 | 1 | 6.6 |
| Fit, then the full bond-graph search past 0.55 Å | 94 | 56 | 0 and 1 | 1.58 |
| PTM, RMSD cutoff 0.10 | 94 | 56 |  | 4.93 |

The kept decision is the last kernel row. On the locked shells it matches PTM, 94 of 94 at 220 K and 56 of 57 at 270 K, and it does not swap ice Ic with ice Ih. At 220 K the frame as a whole is 54 ice Ic, 58 ice Ih, 3 ice Isd and 141 unlabelled. PTM still returns 128 cubic and 128 hexagonal, because its 0.10 cutoff accepts shells the 0.45 Å lock does not. At 270 K the kernel is 42 ice Ic, 39 ice Ih, 4 ice Isd and 171 unlabelled. PTM is 128 cubic, 125 hexagonal and 3 other. CHILL+ stays on the ideal split at 220 K (128 cubic, 128 hexagonal) and is 128 cubic, 126 hexagonal and 2 interfacial at 270 K.

The old nucleus held the outside at 360 K for 10 ps. CHILL+ and PTM still called all 512 oxygens cubic, so the crystal had not melted. The frame used here is a 40 ps melt at 500 K, the 71 oxygens that started inside 8 Å put back on the ice Ic lattice, liquid oxygens that overlapped them removed, and 2 ps at 250 K with the core held fixed. 447 oxygens remain. CHILL+ calls 386 of them water. PTM calls 432 other, 14 cubic and 1 hexagonal. The kernel calls 7 ice Ic, all inside the core, and leaves the other 440 unlabelled.

## Paragraph

Local ice is labelled by the residual of the seventeen-atom shell to ice Ic, ice Ih or stacking-disordered ice. A quaternion overlay of the tetrahedron handles the shells that sit close to a template, and the full bond-graph search handles the rest. The atoms are split across threads. On ideal ice Ic this is faster, per atom, than polyhedral template matching in OVITO 3.17.2 from 216 molecules through 110592, on one thread and on eight. A general IRA overlay of the same shells is kept as a reference and is about a thousand times slower per atom. The overlay that is returned for a named shell, cage or ring uses a 3×3 SVD and a permutation fixed by the bond graph, and the continuous shape measure is the scaled residual of that overlay. The same machinery names the cages of the sI and sII lattices, two and six in sI and sixteen and eight in sII, which template matching reports only as other and which CHILL+ reports only as clathrate. On a stacking-fault cell, template matching returns other for every oxygen, while the shell match splits the cell evenly between ice Ic and ice Isd.

The label kept for production is the template residual, not the eclipsed count. A count with the cutoff at 0.75 agrees with the locked second shell on 68 of 94 oxygens at 220 K and 43 of 57 at 270 K. A quaternion fit of the tetrahedron, followed by the full bond-graph search only when that fit stays above 0.55 Å, agrees on 94 of 94 and 56 of 57, which is PTM at its default RMSD cutoff of 0.10. On ideal ice Ic the fast path is enough, and the same label is faster per atom than PTM from 216 molecules through 110592, on one thread and on eight. A 40 ps melt at 500 K, with the core put back, leaves a liquid shell: CHILL+ calls 386 of 447 oxygens water, where the earlier 360 K hold had left the whole cell cubic. One thread names a hexagon in 3.80 µs and a dodecahedron in 4.82 µs, below IRA.

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
