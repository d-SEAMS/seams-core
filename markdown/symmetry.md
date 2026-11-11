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
| continuous-symmetry-measure/shape | Continuous shape measure reference program | GPL-2.0 | Optional `with_csm` (default auto). Off unless OpenBabel is present. Not linked into libyodaLib. This build leaves it off. The in-tree measure is the published formula. |
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

The command is `seams symmetry`. The C entry points are `seams_symmetry_match`, `seams_symmetry_point_group`, `seams_symmetry_space_group`, `seams_symmetry_template` and `seams_symmetry_crystal` in `seams_c_api.h`. Python is `bindings/python/symmetry.py`.

## Kernel

The production label of a bulk oxygen is not a general permutation search. A periodic cell list gathers the sixteen nearest neighbours inside 5 Å. The four nearest of those are the tetrahedral shell. Each bond is eclipsed or staggered according to the azimuthal alignment of the outer neighbours with the other three bonds of the centre: eclipsed when the best cosine exceeds 0.75, which separates the ideal values 1 and 0.5. Zero eclipsed bonds is ice Ic, one is ice Ih, and two is ice Isd. Anything else is left unlabelled by this pass. The same count is what the three seventeen-atom templates were built from, so on the ideal cells the label matches a full IRA overlay of the shell.

A named shell, cage or ring still returns a rotation, a permutation and an RMSD. That alignment is a 3×3 SVD of the covariance (Kabsch), tried proper and improper, on a permutation taken from the bond graph rather than from N!. Rings use the cyclic order. Ice uses the 24 assignments of the first shell and then the bonded second neighbours. Cages grow a degree-3 matching from each vertex. The continuous shape measure is the scaled residual of that overlay, 100 times the mean square distance to the template divided by the mean square radius of the template (Pinsky and Avnir 1998; Alvarez et al. 2005). It is computed in-tree. The GPL-2 shape program is a separate optional build and is not on this path.

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

IRA on the same 216-atom cell, one shell per oxygen through libira, is serial. Median 283.84 ms (raw 282.74 282.963 283.56 284.119 285.766 286.336), about 1.31 ms per atom, three orders of magnitude above the kernel. spglib on the 8-atom ice Ic cell returns 227 `Fd-3m` in a median 6.13 ms. That is one cell, not a per-atom cost. The in-tree shape measure on one ideal shell, scaled to 216 atoms, is a median 14.3 ms. The GPL shape program was not timed: OpenBabel is absent, so `with_csm` stays off.

## What PTM does not name

On the stacking-fault cell, PTM returns OTHER for all 8 oxygens. The IRA shell match returns ice Ic for 4 and ice Isd for 4. CHILL+ on those 8 oxygens taken out of the periodic cell, so the crystal neighbours are missing, returns water for all 8. It has no stacking-fault label.

On the sI oxygen lattice (46 waters, cell 12.239 Å) PTM returns OTHER for every water and CHILL+ returns clathrate for every water. The ring census finds 48 five-rings and 54 six-rings, and the closed cages are 2 of 5^12 and 6 of 5^12 6^2. On the sII lattice (136 waters, nearest neighbour scaled to 2.80 Å) PTM is again OTHER throughout and CHILL+ is again clathrate throughout. The cages are 16 of 5^12 and 8 of 5^12 6^4. Neither PTM nor CHILL+ separates those cage types.

A 64-atom ice Ic crystallite in a box with 128 uncorrelated points is the small stand-in for a nucleus. The kernel and PTM both mark 8 core atoms (ice Ic, cubic diamond) and leave the 128 points unlabelled. CHILL+ marks 20 atoms cubic, and also 7 interfacial, 7 inter-clathrate and 1 clathrate, because four neighbours are enough for its bond test. The kernel and PTM ask for the second shell, so they agree with each other and disagree with those extra 12 cubic labels.

## Paragraph

Local ice is labelled from the eclipsed-bond count of the tetrahedral shell, which is the same geometric distinction as the seventeen-atom templates of ice Ic, ice Ih and stacking-disordered ice. The count is evaluated from a periodic cell list, with no heap traffic in the atom loop and with the atoms split across threads. On ideal ice Ic this is faster, per atom, than polyhedral template matching in OVITO 3.17.2 from 216 molecules through 110592, on one thread and on eight. A general IRA overlay of the same shells is kept as a reference and is about a thousand times slower per atom. The overlay that is returned for a named shell, cage or ring uses a 3×3 SVD and a permutation fixed by the bond graph, and the continuous shape measure is the scaled residual of that overlay. The same machinery names the cages of the sI and sII lattices, two and six in sI and sixteen and eight in sII, which template matching reports only as other and which CHILL+ reports only as clathrate. On a stacking-fault cell, template matching returns other for every oxygen, while the shell match splits the cell evenly between ice Ic and ice Isd.

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
