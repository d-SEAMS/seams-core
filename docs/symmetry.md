# Local symmetry

The engine names a local environment as ice Ih, ice Ic, stacking-disordered ice, a hydrate cage, or a five- or six-ring, and it labels every oxygen in a cell. The per-atom ice label is the production path. Cage and ring names, the point group, and the space group are the same templates asked one cluster at a time.

The Sphinx page is `docs/orgmode/explanation/symmetry.org`. This file is the same account.

## What is labelled

A crystal oxygen is ice Ih (0), ice Ic (1), or stacking-disordered ice (2). Anything else is −1. Only those two diamond structures count as crystal when the same frame is scored with OVITO. The first-neighbour and second-neighbour tags from identify diamond are not crystal, and polyhedral template matching counts only cubic diamond and hexagonal diamond.

The seventeen-atom shell is the central oxygen, its four nearest neighbours, and the twelve neighbours of those four. That is the shell polyhedral template matching uses for diamond (Larsen, Schmidt and Schiøtz, *Modelling and Simulation in Materials Science and Engineering* 24, 055007, 2016, [10.1088/0965-0393/24/5/055007](https://doi.org/10.1088/0965-0393/24/5/055007)).

## Neighbours

The cell is orthorhombic. Atoms are wrapped into the box, histogrammed into cells of the bond cutoff, and stored as three coordinate columns. The bond cutoff is 3.5 Å and a neighbour closer than 2.15 Å is dropped. One thread walks each pair once and writes both ends of the bond. Interior cells use a fixed half of the 27-cell stencil, so that walk does not build a neighbour list for the cell. The four nearest bonds are kept by distance and then sorted, because the common-neighbour test and the fit both read them in that order. Each accepted pair is offered to those four seats, and a closer pair replaces the farthest. A frame with more than sixteen accepted pairs per atom still keeps the nearer neighbours. Eight threads scan the full stencil and write one end. The coordinate columns and the bond rows are allocated once and reused on the next frame of the same size.

A centre whose four bonds do not all lie between 2.30 Å and 3.35 Å does not open a second shell.

## Tetrahedral order

The tetrahedral order parameter *q* is the Chau–Hardwick / Errington–Debenedetti invariant of the four bonds (Errington and Debenedetti, *Nature* 409, 318, 2001, [10.1038/35053024](https://doi.org/10.1038/35053024)). A shell with *q* below 0.94 does not open the second shell. That gate was chosen on the training seeds below. A shell that does open and then measures *q* below 0.88 is liquid.

## Common-neighbour fingerprint

The second shell is the adaptive-cutoff common-neighbour signature of Maras, Trushin, Stukowski, Albe and Nordlund, *Computer Physics Communications* 205, 13 (2016), Appendix A ([10.1016/j.cpc.2016.04.001](https://doi.org/10.1016/j.cpc.2016.04.001)). The cutoff for a bond of length *r* is *r*(1+√2)/2. The twelve neighbour distances are compared with that cutoff in one vector pass.

An exact cubic signature, twelve (4,2,1) patterns, is ice Ic. An exact hexagonal signature, six (4,2,1) and six (4,2,2), is ice Ih. A partial cubic signature with at least eight (4,2,1) patterns and *q* at least 0.97 is ice Ic. A partial hexagonal signature with *q* at least 0.97 is ice Ih. Six (4,2,1) patterns are not labelled here: on the training seeds that count is the liquid false positive, and the shell is passed to the fit instead. An ideal shell, equal bonds and one second-shell radius, is named from the raw cosine of each second-shell bond with another first-shell bond and does not build this signature. The bond vectors are not projected. On the ideal cell that cosine is exactly 1/3 for a staggered arm and exactly 7/9 for an eclipsed arm. The fast path compares the squares, 1/9 and 49/81: an arm is eclipsed when the best square is at least 0.36, and staggered when that square is at most 0.25. Four staggered arms are ice Ic. One eclipsed arm and three staggered arms are ice Ih.

## The fit

A shell is fitted only when it already has one crystal neighbour, *q* is at least 0.95, the four bonds differ by at most 15 percent, and the second-shell radius varies by at most a factor 1.6. The residual is scale-free. Both clouds are centred, the optimal rotation is a quaternion (Horn, *Journal of the Optical Society of America A* 4, 629, 1987, [10.1364/JOSAA.4.000629](https://doi.org/10.1364/JOSAA.4.000629); Coutsias, Seok and Dill, *Journal of Computational Chemistry* 25, 1849, 2004, [10.1002/jcc.20110](https://doi.org/10.1002/jcc.20110)), and the residual is divided by the template radius, which is the normalisation Larsen, Schmidt and Schiøtz use. Zero is a perfect shape. The ceiling is 0.12.

The bond-angle fingerprint names the one template the fit tries. This cosine is the azimuthal one: both bonds are projected onto the plane perpendicular to the central bond. Ideal staggered ice lies at 60 degrees, cosine 1/2, and ideal eclipsed ice lies at 0 degrees, cosine 1. An arm above 0.90 is eclipsed. The ideal-cell fast path uses the raw cosine above. Four eclipsed arms are liquid, except on the recovery pass, where that pattern is hexagonal. The eclipsed arm is anchored to the eclipsed arm of the template. The four-point orientation uses one Jacobi sweep and the seventeen-atom residual uses two. A shell that is still far after a handful of correspondences is not crystal, and a shell with only one crystal neighbour stops earlier than a shell with two. Hexagonal recovery is one pass over the atoms the first fit left unlabelled, and it tries the other templates.

After the labels are written, a crystal atom that has no crystal neighbour among its four bonds is returned to liquid. The tests read the labels from the fit, before that filter.

## Cages, rings, and the other codes

The same in-tree matcher names a cage or a ring. It is a 3×3 SVD with the permutation fixed by the bond graph. The continuous shape measure of that overlay is the scaled residual of Pinsky and Avnir (*Inorganic Chemistry* 37, 5575, 1998, [10.1021/ic9804925](https://doi.org/10.1021/ic9804925)) and of the Alvarez shape program (*Coordination Chemistry Reviews* 249, 1693, 2005, [10.1016/j.ccr.2005.03.031](https://doi.org/10.1016/j.ccr.2005.03.031)). The shape program is not in the tree.

IRA (Gunde, Salles, Hémeryck and Martin-Samos, *Journal of Chemical Information and Modeling* 61, 5446, 2021, [10.1021/acs.jcim.1c00567](https://doi.org/10.1021/acs.jcim.1c00567)) is an optional second matcher for the same templates. SOFI, from the same group (*Journal of Chemical Physics* 161, 062503, 2024, [10.1063/5.0215689](https://doi.org/10.1063/5.0215689)), names the point group. spglib (Togo, Shinohara and Tanaka, *Science and Technology of Advanced Materials: Methods* 4, 2384822, 2024, [10.1080/27660400.2024.2384822](https://doi.org/10.1080/27660400.2024.2384822)) names the space group. CHILL+ (Nguyen and Molinero, *Journal of Physical Chemistry B* 119, 9369, 2015, [10.1021/jp510289t](https://doi.org/10.1021/jp510289t)) is a different label and is not this shell.

## Calling it

C++, in `symm`:

- `classifyIce(cartesian, vectors, kind)` writes one label per row. `kind` is ice Ih, ice Ic, ice Isd, or −1. The default ceiling is `kShellRmsd` (0.12).
- `classifyIceArrays(x, y, z, n, lattice, kind)` is the same call on three coordinate columns.
- `best`, `matchKernel` and `matchAll` name one environment. Template ids are ice Ih, ice Ic, ice Isd, 5^12, 5^12 6^2, 5^12 6^4, ring-5 and ring-6.
- `pointGroup` needs SOFI. `spaceGroup` needs spglib. Either returns status 2 when that library was not linked.
- `seams symmetry` on the command line prints the space group, one shell overlay, and the label tally.

C, in `seams_c_api.h`:

- `seams_symmetry_classify` and `seams_symmetry_classify_tol` take row-major coordinates and three lattice rows, and write the same labels.
- `seams_symmetry_match`, `seams_symmetry_template`, `seams_symmetry_crystal`, `seams_symmetry_point_group`, `seams_symmetry_space_group` and `seams_symmetry_backend` cover one cluster.

Python, `bindings/python/symmetry.py`, loads `libyodaLib` through ctypes. Set `SEAMS_CORE_LIB` to the built library. `classify(cartesian, lattice, shell_rmsd=0.12)` is the per-atom label.

## Build options

`with_ira` is a feature option, `auto` by default. When libira links, the matcher can use IRA and SOFI. IRA is the Apache-2.0 option of its dual licence, and it is not required to build.

`with_csm` is a feature option, `disabled` by default. The continuous shape measure program is GPL-2.0. It is an external executable (`csm_shape` or `shape`), it is not vendored, and it is not linked into `libyodaLib`. The Catch2 case that compares the in-tree measure with that program runs only when `CSM_SHAPE` is set. `with_spglib` stays `auto` and links libsymspg when it is present.

On GCC and Clang the shell translation unit is a separate static library built with `-mavx2 -mfma`. Other compilers compile it into `libyodaLib` with the ordinary flags. A release build is `-O3`.

## Held-out protocol

The frames are the steepest-descent quenches of three mW nucleation seeds, 91001, 92002 and 93003, each a 12×12×12 ice Ic cell that was melted, reseeded, and run at 220 K. Seeds 91001 and 92002 are the training pair. They chose the *q* gate and the ceiling 0.12. Seed 93003 was not used to choose either, and it is the held-out frame.

Three methods label each quench with their own production rule: this kernel, CHILL+, and polyhedral template matching at an RMSD cutoff of 0.1 with cubic and hexagonal diamond enabled. The majority truth is the template named by at least two of the three. The persistent truth keeps an atom only when the majority template is the same on the 190 ps quench and the 200 ps quench. A prediction is a true positive when the truth is crystal and the prediction names that same template, a false positive when the prediction is crystal and is not the truth, and a false negative when the truth is crystal and the prediction differs.

On the held-out quench the kernel at ceiling 0.12 is:

| truth | true positive | false positive | false negative | precision | recall |
| --- | ---: | ---: | ---: | ---: | ---: |
| majority | 328 | 6 | 63 | 0.9820 | 0.8389 |
| persistent | 319 | 15 | 47 | 0.9551 | 0.8716 |

Polyhedral template matching on the same two truths is majority precision 0.9635 and recall 0.8107 (317 true positives, 12 false positives, 74 false negatives), and persistent precision 0.9453 and recall 0.8497. The kernel is ahead of template matching on both precision and recall. Raising the ceiling to 0.125 on the training pair adds false positives as well as true positives, and an embedded-only relaxation of the ceiling did the same, so the production ceiling stays 0.12. The held-out recall is the number above. It was not pushed up by looking at seed 93003.

## Real frames

The classifier rule is the one above. These frames were not used to choose it. An mW site is one molecule. Crystal for the kernel is ice Ih, ice Ic, or stacking-disordered ice. Crystal for polyhedral template matching and for identify diamond is cubic or hexagonal diamond only. The first-neighbour and second-neighbour tags from identify diamond are counted apart from crystal. CHILL+ is the 2015 rule in this tree, cutoff 3.5 Å: four nearest neighbours, the staggered and eclipsed windows of Nguyen and Molinero, and the neighbour test in footnote a of their Table 1. It is not the OVITO Chill+ modifier.

The interface frame is `interface_f0.lammpstrj` from `science/interface-order` at `f5ca8ffc`, 864 mW molecules in a slab. The nucleation file is figshare [10.6084/m9.figshare.11448702.v1](https://doi.org/10.6084/m9.figshare.11448702.v1), md5 `c28fc209a2854d2ba49e8f9b55accbe8`, 148 frames of 4096 molecules, the trajectory read by `science/nucleation-pathways` at `170a6c18`. Its first frame is the shipped `input/traj/nucleation_f1.lammpstrj`.


Interface slab, kernel against polyhedral template matching:

| kernel | Ih | Ic | other |
| --- | ---: | ---: | ---: |
| Ih | 233 | 0 | 12 |
| Ic | 0 | 0 | 0 |
| Isd | 0 | 0 | 0 |
| other | 66 | 0 | 553 |

Interface slab, kernel against identify diamond:

| kernel | Ih | Ic | Ih neighbour | Ic neighbour | Ih second | Ic second | other |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Ih | 232 | 0 | 13 | 0 | 0 | 0 | 0 |
| other | 81 | 0 | 67 | 0 | 72 | 0 | 399 |

Interface slab, kernel against CHILL+:

| kernel | hexagonal | cubic | interfacial | clathrate | interfacial clathrate | water |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Ih | 245 | 0 | 0 | 0 | 0 | 0 |
| other | 114 | 1 | 67 | 13 | 69 | 355 |

Nucleation frame 1, kernel against polyhedral template matching:

| kernel | Ih | Ic | other |
| --- | ---: | ---: | ---: |
| Ih | 7 | 0 | 1 |
| Ic | 0 | 9 | 3 |
| Isd | 0 | 0 | 1 |
| other | 4 | 4 | 4067 |

Nucleation frame 1, kernel against identify diamond:

| kernel | Ih | Ic | Ih neighbour | Ic neighbour | Ih second | Ic second | other |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Ih | 7 | 0 | 0 | 1 | 0 | 0 | 0 |
| Ic | 0 | 8 | 2 | 1 | 0 | 1 | 0 |
| Isd | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| other | 7 | 5 | 28 | 18 | 68 | 36 | 3913 |

Nucleation frame 1, kernel against CHILL+:

| kernel | hexagonal | cubic | interfacial | clathrate | interfacial clathrate | water |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Ih | 8 | 0 | 0 | 0 | 0 | 0 |
| Ic | 0 | 12 | 0 | 0 | 0 | 0 |
| Isd | 0 | 0 | 1 | 0 | 0 | 0 |
| other | 59 | 36 | 272 | 49 | 305 | 3354 |

The one stacking-disordered molecule on that frame is CHILL+ interfacial. Across all 148 frames the 1037 stacking-disordered molecule-frames are all template-matching liquid and all CHILL+ interfacial.

Over the 148 frames, 606208 molecule-frames, the kernel calls 7342 hexagonal, 10020 cubic and 1037 stacking-disordered. Template matching calls 7838 hexagonal and 10771 cubic. Identify diamond calls 8423 hexagonal and 11336 cubic, and a further 19417 first-neighbour plus 27904 second-neighbour tags. CHILL+ calls 19760 hexagonal and 17869 cubic, plus 46992 interfacial, 6468 clathrate and 50164 interfacial clathrate. Kernel hexagonal is CHILL+ hexagonal on 7264 of 7342 molecule-frames, and kernel cubic is CHILL+ cubic on 9980 of 10020. The extra CHILL+ crystal, 12495 hexagonal and 7889 cubic, is kernel liquid.

The disagreements follow from what each rule is allowed to name. CHILL+ cubic and hexagonal need four bond correlations in the staggered or eclipsed windows. A molecule can pass that table and still miss the kernel: the kernel also asks for a second-shell common-neighbour signature, a scale-free residual at most 0.12, and one crystal neighbour. On the slab, 114 of the 359 CHILL+ hexagonal molecules are kernel liquid. On the trajectory, the same gap is the 12495 hexagonal and 7889 cubic molecule-frames above. CHILL+ interfacial ice is two or three staggered bonds. The kernel has no interfacial label, so those molecules are liquid, except the stacking-disordered shells, which land in the interfacial row of CHILL+ because they are not four bonds of one kind. Template matching has cubic and hexagonal diamond templates and no stacking-disordered template, so every stacking-disordered label is template-matching liquid. The remaining template-matching mismatches are shells on opposite sides of the two cutoffs: an isolated template-matching diamond can fail the kernel's crystal-neighbour test, and a shell the kernel fits can sit outside the 0.1 RMSD cutoff. Identify diamond's neighbour tags are the molecules touching a diamond molecule. The kernel does not promote them to crystal, and that tag is most of the off-diagonal with identify diamond. Where both methods call the molecule itself crystal, hexagonal stays with hexagonal and cubic stays with cubic.

IRA and SOFI are optional. A build links them only when `with_ira` finds libira; the default build does not. The numbers in this paragraph are from one build that did link them. In that build the matcher backend is IRA. SOFI is asked for the point group of the same 17-atom shell the ice templates use: the centre, its four neighbours between 2.15 Å and 3.5 Å, and three further neighbours of each of those four, in the minimum image. A shell that cannot be filled is left unmatched (68 molecules on the slab, 377 on nucleation frame 1). IRA returns a nearest ice template for every complete shell, because the only 17-atom templates are the three ice shells. On both frames every kernel crystal molecule is given that same template. The IRA RMSD, in ångströms, separates them from the liquid: on the slab the hexagonal median is 0.314 and the liquid median is 1.194; on frame 1 the hexagonal, cubic and stacking-disordered medians are 0.342, 0.377 and 0.458, and the liquid median is 1.112. The liquid still has a low tail, so the nearest-template name alone is not a crystal label. SOFI at its default threshold 0.05 returns C1 for every complete shell on both frames, crystal and liquid alike, so the point group does not separate ice Ih from ice Ic here.

## Timings

One sitting, 11 November 2026, an Intel Xeon with 8 cores. `nproc` was read with `OMP_NUM_THREADS` unset and reported 8. The compiler is g++ 15.2.0. The library is a Meson release build, `-O3`, with `-mavx2 -mfma` on the shell translation unit. Eight samples, one warm-up discarded, and the median of an even count is the mean of the two central samples. OVITO is 3.17.2. Identify diamond is `IdentifyDiamondModifier`. Template matching is `PolyhedralTemplateMatchingModifier` with the RMSD cutoff 0.1 and only cubic diamond and hexagonal diamond enabled.

The bond stage was the largest piece of a one-thread call, about 3 ms of 5.3 ms at 13553 molecules and almost the whole 24 ms call on the ideal 110592-molecule cell. The pair walk now keeps the minimum-image vector it has already computed, and the later pass that rebuilt those vectors from the coordinates is gone. Inside the template fit, the second-shell arms are ordered by the unnormalised cross product, and an eclipsed arm is a squared-cosine test. Common-neighbour rows are built one at a time, and a shell that can no longer reach the ambiguous band does not build the rest.

The timed frames are the 13553-molecule cells, not the quenches. One thread, eight samples, milliseconds, kernel then identify diamond in one process. The range is the slowest and fastest sample.

| seed | kernel | identify diamond |
| --- | --- | --- |
| 91001 | 5.436 (5.325–5.485) | 6.618 (5.179–7.773) |
| 92002 | 5.099 (4.987–5.281) | 5.939 (5.043–6.472) |
| 93003 | 5.251 (5.163–5.342) | 6.176 (5.225–7.059) |

The eight kernel samples are 5.429 5.485 5.392 5.484 5.364 5.442 5.325 5.463, then 5.092 5.281 5.016 5.106 4.987 5.172 5.071 5.144, then 5.298 5.166 5.341 5.219 5.342 5.163 5.283 5.200. Identify diamond on the same frames is 6.726 6.509 7.119 7.773 6.057 5.179 6.934 5.550, then 5.815 6.457 5.178 6.472 5.043 6.062 5.605 6.314, then 5.739 6.272 7.059 5.918 6.314 5.225 6.247 6.104. On every seed the kernel median is below the identify-diamond median, and the kernel's own spread is under 0.3 ms. The fastest identify-diamond sample on each seed is still below the slowest kernel sample, by 0.31 ms on 91001, 0.24 ms on 92002 and 0.12 ms on 93003.

Eight threads, same cells, measured in the preceding minute:

| seed | kernel | identify diamond |
| --- | --- | --- |
| 91001 | 2.222 (2.153–2.296) | 6.087 (5.590–13.502) |
| 92002 | 1.940 (1.912–2.076) | 5.787 (4.925–6.548) |
| 93003 | 2.009 (1.958–2.056) | 5.714 (5.075–6.425) |

The kernel's slowest eight-thread sample on each seed is under identify diamond's fastest sample.

Template matching was measured in that same preceding run, kernel and template matching in one process. It does not speed up at eight threads. Per-atom times are the median divided by the atom count.

| frame | threads | kernel µs/atom | template matching µs/atom |
| --- | --- | --- | --- |
| 91001, 13553 | 1 | 0.406 | 4.953 |
| 92002, 13553 | 1 | 0.380 | 4.655 |
| 93003, 13553 | 1 | 0.399 | 4.698 |
| ideal ice Ic, 110592 | 1 | 0.218 | 3.644 |
| 91001, 13553 | 8 | 0.177 | 4.644 |
| 92002, 13553 | 8 | 0.146 | 4.658 |
| 93003, 13553 | 8 | 0.151 | 4.713 |
| ideal ice Ic, 110592 | 8 | 0.071 | 3.536 |

The classify workspace, the coordinate columns and the canonical shell are one copy per calling thread, so two frames can be labelled at once. A later sitting on the same eight cores timed the C call again, two warm-ups and then eight samples. One thread, microseconds per atom: 91001 0.401, 92002 0.362, 93003 0.376, ideal ice Ic 0.215. Eight threads: 0.169, 0.137, 0.137, and the ideal cell remeasured at 0.071. The one-thread samples, milliseconds, are 5.416 5.168 5.291 5.294 5.446 6.678 5.777 5.561, then 5.054 4.869 4.758 4.831 5.066 4.729 4.950 5.120, then 5.145 5.134 4.954 5.331 5.265 5.053 4.868 5.034, and on the ideal cell 25.693 23.961 23.874 23.780 23.734 24.298 23.393 23.321. Eight threads: 2.276 2.293 2.287 3.296 2.398 4.829 2.170 2.266, then 2.006 2.123 1.913 1.765 1.747 1.815 1.891 1.829, then 1.857 1.868 1.872 1.766 1.801 1.816 1.846 1.897. The seeds sit under the per-atom figures above. The ideal cell is 0.215 against 0.218 at one thread. At eight threads that cell was measured again on 12 November 2026. `nproc` was read with `OMP_NUM_THREADS` unset and reported 8. One process ran `tests/bench_shell_scale` on the release build, one warm-up and then eight calls. The median of an even count is the mean of the two central samples, and the interquartile range is the upper quartile minus the lower. The median is 7.804 ms, the interquartile range is 0.148 ms, and that is 0.071 µs/atom. Template matching was timed in that same sitting, eight threads, the same warm-up and eight calls, RMSD cutoff 0.1, only cubic diamond and hexagonal diamond. Its median is 391.012 ms, its interquartile range is 4.280 ms, and that is 3.536 µs/atom. The kernel is faster per atom.

The one-thread template-matching samples, milliseconds, are 63.493 64.348 65.309 71.131 81.903 72.515 68.933 63.321, then 64.195 62.642 65.019 63.439 64.351 61.660 62.748 61.902, then 64.377 62.895 63.863 63.724 63.622 64.526 63.316 63.037, and on the ideal cell 411.364 406.241 396.372 399.726 407.844 396.695 399.634 415.364. The ideal cell is 24×24×24 ice Ic, 110592 molecules, and every atom is labelled cubic. The one-thread kernel on that cell is 24.262 24.381 23.859 24.137 24.262 24.060 24.121 24.174. At eight threads the ideal kernel is 7.661 7.767 7.790 7.792 7.816 7.906 8.012 8.862 and template matching is 394.550 390.626 390.898 400.186 389.956 395.303 391.127 377.933.
