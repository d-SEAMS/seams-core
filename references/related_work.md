# Related work for d-SEAMS 2.0

Survey of ice and clathrate structure analysis from 2022 onward, plus the codes named in the brief (freud, OVITO, MDAnalysis, pyscal, GenIce, and IcePic). Bibliographic data is in `related_work.bib`. DOIs were read back from Crossref. The head-to-head below is new; it does not replace the campaign against OVITO, freud, pyscal, matscipy and the CHILL+ library.

## Codes

| Source | What it does | Against d-SEAMS |
| --- | --- | --- |
| [mdapy 2.0](https://github.com/mushroomfire/mdapy) ([Wu and Shao, CPC 2023](https://doi.org/10.1016/j.cpc.2023.108764)) | Cutoff CHILL+ in `mdapy.identify.ChillPlus` | Competes on CHILL+. Timed below. |
| [fastatomstruct 2.2](https://zivgitlab.uni-muenster.de/ag-salinga/fastatomstruct) ([Holle et al., JOSS 2025](https://doi.org/10.21105/joss.08106)) | King and Guttman rings; Steinhardt `q_l` | Competes on primitive-ring counts and on `q_6`. Timed below. |
| [TRACE](https://github.com/junweihsu/TRACE) ([Hsu and Lin, JCTC 2025](https://doi.org/10.1021/acs.jctc.5c01459)) | DFS rings, cups, complete and incomplete hydrate cages, including one additive. One C++ program. Reads GROMACS `.gro` | Competes on hydrate cage topology. Complements d-SEAMS on additive-coordinated cages. No library API, so it is not in the harness. A full search each frame. |
| [RingDetect-HPC 1.4](https://github.com/ouielba90/RingDetect-HPC) | Python `ringdetect.find_rings` and a C/Fortran CLI. DFS toward a minimum cycle basis. Bonds are covalent radii. The wrapper asks for every OpenMP thread | Does not see a 3.5 Å oxygen hydrogen-bond graph, so it is not a matched ice-ring code. Left out of the harness. |
| [OVITO](https://doi.org/10.1088/0965-0393/18/1/015012), [Find rings](https://www.ovito.org/manual/reference/pipelines/modifiers/find_rings.html) | CHILL+, rings, PTM, diamond, adaptive CNA | Already in the paper's matched comparison. |
| [freud](https://doi.org/10.1016/j.cpc.2020.107275), [repository](https://github.com/glotzerlab/freud) | Neighbour lists and Steinhardt `q_l` | Already matched on `q_6`. Complements cages. |
| [pyscal3](https://doi.org/10.21105/joss.01824), [repository](https://github.com/pyscal/pyscal3) | `q_l`, CNA, diamond structure | Already matched on `q_6` and timed on diamond structure. |
| [MDAnalysis](https://doi.org/10.1002/jcc.21787), [SciPy 2016](https://doi.org/10.25080/Majora-629e541a-00e), [repository](https://github.com/MDAnalysis/mdanalysis) | Trajectories, distances, hydrogen bonds | Complements. No ice or cage classifier in the library. |
| [GenIce](https://doi.org/10.1002/jcc.25077), [GenIce-core](https://doi.org/10.1063/5.0198056), [cycless](https://github.com/vitroid/cycless) | Builds hydrogen-disordered ice. cycless counts cycles | Generation complements d-SEAMS. cycless rings are already in the campaign and disagree on a confined frame. |
| [matscipy](https://doi.org/10.21105/joss.05668) | Primitive rings | Already matched. |
| [atomes](https://doi.org/10.1016/j.commatsci.2025.113805) | Ring statistics, editing and visualization | Competes on ring counts. A desktop program, not a library called from this harness. |
| [RingsStatisticsMatter](https://github.com/MorrowChem/RingsStatisticsMatter) ([0.2.0](https://pypi.org/project/rings-statistics-matter/0.2.0/), [Yuan and Cormack](https://doi.org/10.1016/s0927-0256(01)00256-7)) | Primitive rings. Rust kernel and a Python wrapper; Julia is optional | Competes on the census. Timed in the afternoon section. A full search each frame. |
| [Palmer CHILL+ library](https://github.com/palmergroupUH/Order_Parameters_Library) | Fortran CHILL+ | Already matched. |
| [mdevaluate `chill.py`](https://gitea.pkm.physik.tu-darmstadt.de/IPKM/mdevaluate) | Staggered and eclipsed bond counts | A non-periodic KDTree. Not the periodic cutoff list, so it was not timed. |
| [ClavisCage 0.3.18](https://doi.org/10.5281/zenodo.24112918) | CHILL+, TRACE, HTR, GRADE, F3/F4, averaged Steinhardt, MCG-3, tetrahedral order. Browser and a Python wheel | Post-processing. CHILL+ is timed in the 12 November section. |
| [GRADE](https://doi.org/10.1016/j.cpc.2019.06.004) | Hydrate structure from cage composition | Earlier cage code. TRACE is the 2025 topological successor. |
| [IceCoder](https://doi.org/10.1021/acs.jctc.4c01298), [repository](https://github.com/TeamSuman/IceCoder) | SOAP plus a variational autoencoder over many ice phases | Competes as a phase tag. Needs a trained model. Does not return cages. |
| [IcePic](https://doi.org/10.1073/pnas.2205347119), [repository](https://github.com/mbdavies13/IcePic) | CNN for the heterogeneous nucleation temperature of a room-temperature contact layer | Complements. It does not label molecules in a trajectory. |

The medicinal-chemistry program IcePick (1999) is a different tool. The ice-nucleation code is IcePic.

## Papers, 2022 onward

Learned classifiers, none of which return hexagonal or double-diamond cages:

- Ishiai, Endo and Yasuoka, [JCP 2023](https://doi.org/10.1063/5.0156203). TeaNet graph network. Nine water molecules suffice for ice versus liquid.
- Ishiai, Endo, Brumby, Sum and Yasuoka, [JCP 2024](https://doi.org/10.1063/5.0177815). The same network on sI, sII and sH.
- Ishiai, Yasuda, Endo and Yasuoka, [JCTC 2024](https://doi.org/10.1021/acs.jctc.3c00995). Unsupervised temporal contrast. Ice, liquid, and ice–liquid and ice–vapour interfaces.
- Moradzadeh, Oliaei and Aluru, [JPCC 2023](https://doi.org/10.1021/acs.jpcc.2c07423). Edge-conditioned graph network for bulk, interface and confined water, including a nanotube.
- Takahashi and Hiratsuka, [Cryst. Growth Des. 2023](https://doi.org/10.1021/acs.cgd.2c01519). MALIO search over local order parameters for ice Ih and cyclopropane sI/sII.
- Sun, Hamel, Hsu, Sadigh, Lordi and Zhou, [JCIM 2024](https://doi.org/10.1021/acs.jcim.4c00822). Score-based denoiser plus SOAP. Ice polymorphs and sI.
- Faure Beaulieu, Deringer and Martelli, [JCP 2024](https://doi.org/10.1063/5.0193340). Neural-network classifiers for amorphous ices.
- Maity and Chakrabarty, IceCoder, above.
- Davies, Fitzner and Michaelides, IcePic, above.

PointNet classification of local ice and hydrate environments is DeFever, Targonski, Hall, Smith and Sarupria, [Chem. Sci. 2019](https://doi.org/10.1039/c9sc02097g). No 2022–2026 PointNet library turned up with a drop-in CHILL+ or cage API. The later papers above use graph networks, a denoiser, or a variational autoencoder.

## Claims a finding can move

1. **First code to return a full ring search from a local update** (abstract and the introduction). fastatomstruct, TRACE, RingDetect-HPC, atomes, cycless, OVITO Find rings and RingsStatisticsMatter.jl recompute the frame. None of those sources describe carrying rings whose bonds did not change. The sentence still stands. A revision should name fastatomstruct and TRACE next to RINGS, matscipy, cycless and OVITO.

2. **Faster than OVITO, freud, pyscal, matscipy and the CHILL+ library in every matched comparison** (abstract and the comparison section). That list is unchanged: mdapy and fastatomstruct were not in it. On this protocol mdapy's cutoff CHILL+ is faster than d-SEAMS at 1 000, 10 648 and 110 592 molecules, with the same cubic label on every molecule. A sentence widened to every CHILL+ library would be false. fastatomstruct `q_6` is faster at 1 000 molecules and slower at 10 648 and at 110 592. King's rings are slower at all three sizes and return the same census (two six-rings per molecule). At 110 592 molecules the King median is about 410 times the d-SEAMS median, and the six samples spread from 447 s to 544 s.

3. **Ring-statistics codes stop at counts, while d-SEAMS builds cages** (introduction). True of RINGS, cycless, matscipy, fastatomstruct and OVITO Find rings. False if read as "no other code builds cages": GRADE and TRACE assemble hydrate cages. TRACE also keeps incomplete cages and one additive, which the d-SEAMS signature cages do not.

4. **OVITO, freud and pyscal have no nanotube-ice class** (introduction). Still true of those three. Moradzadeh et al. label confined water liquid versus solid, including inside a nanotube. That is a phase tag, not the prism blocks.

5. **A LAMMPS compute assigns CHILL+ during the run** (program summary and the implementation section). The manuscript does not claim a GROMACS or OpenMM fix. Clavis and TRACE read trajectories afterwards. mdapy is a post-processing library. Nothing surveyed here assigns CHILL+ or cages from inside an MD step except the d-SEAMS LAMMPS compute.

6. **CHILL+ follows Nguyen and Molinero, with the OVITO modifier as a second option.** mdapy uses the same staggered (`c < -0.8`) and eclipsed (`-0.35 < c < 0.25`) cuts on neighbours inside 3.5 Å. On these ideal ice Ic frames the second shell is beyond 3.5 Å, so that cutoff list is the four nearest neighbours. Agreement is 1 on every molecule.

7. **Cage labels come from the topology of the network, and local classifiers answer a different question** (comparison and conclusions). The graph-network, denoiser, autoencoder and MALIO papers are further local or learned tags. They do not replace the hexagonal-cage and double-diamond-cage counts.

## Head-to-head

Ideal ice Ic, oxygen bond 2.75 Å, cutoff 3.5 Å, rings of at most six. One core (`taskset -c 0`), one thread, one warm-up, six timed calls. Median and type-7 interquartile range. Frames were built in a temporary directory and not stored. d-SEAMS is `libseams_lang` 2.13.0. mdapy is 2.0.0. fastatomstruct is 2.2.0.

CHILL+ is the cutoff graph. `q_6` for d-SEAMS is the mutual four-neighbour list; fastatomstruct uses `q_l` with the near-hard weight `p_1 = 48`. King's rings are compared with the cutoff-graph primitive rings. Agreement is the fraction of molecules with the same merged CHILL+ class, the fraction of `q_6` values within `1e-6`, or 1 when the ring census matches.

| Molecules | Operation | d-SEAMS median (IQR) ms | Other median (IQR) ms | Agreement |
| --- | --- | --- | --- | --- |
| 1000 | CHILL+ | 4.376 (0.732) | mdapy 0.903 (0.026) | 1 |
| 1000 | `q_6` | 2.718 (0.111) | fastatomstruct 2.369 (0.127) | 1, max abs 3.4e-9 |
| 1000 | rings | 10.297 (0.158) | fastatomstruct 24.609 (0.260) | 1; both 2000 six-rings |
| 10648 | CHILL+ | 46.807 (1.144) | mdapy 7.128 (0.130) | 1 |
| 10648 | `q_6` | 27.034 (0.535) | fastatomstruct 94.520 (5.535) | 1, max abs 3.4e-9 |
| 10648 | rings | 111.135 (1.260) | fastatomstruct 1886.682 (2.615) | 1; both 21296 six-rings |
| 110592 | CHILL+ | 469.734 (16.309) | mdapy 63.300 (3.223) | 1 |
| 110592 | `q_6` | 279.387 (3.831) | fastatomstruct 7703.540 (32.572) | 1, max abs 3.4e-9 |
| 110592 | rings | 1132.576 (6.633) | fastatomstruct 464469.420 (57784.377) | 1; both 221184 six-rings |

The samples are `benchmarks/h2h-watch.csv`. `snakemake related_watch` regenerates that file. `benchmarks/h2h_watch.py` is the one-call Airspeed Velocity arm of the same commands.

## Matched cases added in v2

pyscal3 4.1.0 exports `steinhardt_parameter` and `diamond_structure`. It has no CHILL+ and no ring census, so those two operations are not a matched case. The Palmer library at commit `ee8a3cb` does not link: `CHILL_mod` imports `gen_cluster_nb_list`, `gen_cluster_ID_list`, `get_largest_cluster` and `single_linkage_cluster_by_Jeremy`, and none of those symbols is in `single_linkage_cluster.f90`. OVITO 3.17.2 from PyPI ran `ChillPlusModifier` and `FindRingsModifier` with no licence prompt once `libGL` was on the loader path.

Same protocol, ideal ice Ic, cutoff 3.5 Å, rings of at most six. OVITO's ring timer includes `CreateBondsModifier`. matscipy uses `ring_statistics(..., maxlength=7)`, which counts rings of size at most six. cycless is `pairlist.pairs_iter` plus `cycles_iter`. Agreement is 1 wherever both codes returned a result: every CHILL+ label is cubic, and both ring censuses are two six-rings per molecule.

| Molecules | Operation | d-SEAMS median (IQR) ms | Other median (IQR) ms |
| --- | --- | --- | --- |
| 1000 | CHILL+ | 4.140 (0.066) | OVITO 6.137 (0.165) |
| 1000 | rings | 10.199 (0.090) | OVITO 4.760 (0.193) |
| 1000 | rings | 10.199 (0.090) | matscipy 15.359 (0.138) |
| 1000 | rings | 10.199 (0.090) | cycless 3366.593 (153.502) |
| 10648 | CHILL+ | 46.009 (0.136) | OVITO 64.663 (0.763) |
| 10648 | rings | 111.588 (1.962) | OVITO 48.489 (9.627) |
| 10648 | rings | 111.588 (1.962) | matscipy 2642.873 (581.054) |
| 10648 | rings | 111.588 (1.962) | cycless 36801.907 (634.919) |
| 110592 | CHILL+ | 481.726 (3.800) | OVITO 629.005 (3.107) |
| 110592 | rings | 1132.111 (11.657) | OVITO 564.356 (3.353) |
| 110592 | rings | 1132.111 (11.657) | matscipy signal 11 |
| 110592 | rings | 1132.111 (11.657) | cycless 398433.762 (3764.967) |

OVITO Find rings is about twice as fast as d-SEAMS at every size, with the same census. OVITO CHILL+ is slower, and the labels match. matscipy agrees at 1 000 and 10 648 molecules and is slower; at 110 592 molecules the worker dies with signal 11, the same abort the paper already records for large frames. cycless agrees on this crystal and is slower by about 330 times at 1 000 molecules and about 350 times at 110 592. Samples are `benchmarks/h2h-watch-v2.csv`. `snakemake related_watch_v2` regenerates them.

The abstract's speed sentence is scoped to the main campaign. It must not be read as every later OVITO ring timing: on this watch, Find rings is the faster code.

## 2025–2026

No new code in this window assigns CHILL+ or cages from inside a PLUMED, LAMMPS or OpenMM step. LAMMPS `fix plumed` and [openmm-plumed](https://github.com/openmm/openmm-plumed) still pass collective variables. The manuscript's in-situ claim remains the LAMMPS compute.

Order parameters and phase tags, none of which return hexagonal or double-diamond cages:

- Sinaeian and Haji-Akbari, [JCP 2025](https://doi.org/10.1063/5.0263587). For Steinhardt coordinates of homogeneous and heterogeneous ice nucleation, the distance cutoff and the \(q_6\) cutoff dominate. Chain exclusion and a hydration shell are not uniformly helpful. This does not touch the speed sentence. It does warn against treating one \(q_6\) threshold as the reaction coordinate.
- Kadre, Ghoshdastidar, Dhankar and Punnathanam, [JCP 2025](https://doi.org/10.1063/5.0273484). XGBoost labels of gas-hydrate polymorphs from a learned set of generic order parameters. A phase tag, not a cage finder.
- Harless, Sindelar and Gezelter, [JCP 2025](https://doi.org/10.1063/5.0282695). Ice-active solutes at the ice–water interface. Complements salt and solute work. No CHILL+ replacement.
- Wu and others, [JFM 2025](https://doi.org/10.1017/jfm.2025.10633). Salt exclusion during directional freezing of brine. Not a per-molecule order parameter.
- Wang, Du, Xie, Calzavarini and Sun, [Phys. Rev. Fluids 2026](https://doi.org/10.1103/8ghn-yrlr). Freezing and aging of saline water under convection. Same limitation.
- Zhu and others, [JPCB 2025](https://doi.org/10.1021/acs.jpcb.5c01346). Phase change of a confined bilayer of ice. Not a nanotube prism label.
- Zhu, Zhang, Wang and Su, [Mod. Phys. Lett. B 2025](https://doi.org/10.1142/s0217984925502690). Terahertz driving of water inside a carbon nanotube. Not a classifier.
- Torres-Arenas, Fernández-Fernández, Pérez-Rodríguez and Piñeiro, [arXiv:2605.14533](https://arxiv.org/abs/2605.14533) (2026). Methane hydrate between silica slits. F3 and F4 show in-plane order where a three-dimensional hydrate does not form. Preprint; Crossref has no DOI for it. Complements confined-hydrate analysis. No cage code.

## Proposed one-line fixes

Introduction, the ring-code sentence. Name fastatomstruct and TRACE:

> Ring statistics codes for network solids count these rings, among them RINGS, ISAACS, atomes, matscipy and fastatomstruct, the cycless package that GenIce uses, the ring finder of OVITO, and TRACE.

Speed sentence. Keep the five named codes. The mdapy caveat below was the slow tree (`f3e9e394`). On the fast tree, cutoff CHILL+ is faster than mdapy 2.0 as well:

> In all matched comparisons of the main campaign, d-SEAMS 2.0 is faster than OVITO, freud, pyscal, matscipy and a CHILL library. On the fast tree it is also faster than mdapy 2.0 cutoff CHILL+ on ideal ice Ic, at one core and at eight.

Hydrate cages. The following sentence currently says d-SEAMS is the code that assembles rings into cages:

> Unlike the ring-count codes, d-SEAMS assembles the hexagonal and double-diamond cages of ice, while TRACE assembles hydrate cages, including incomplete cages and cages coordinated to one additive.

A further scope line. The slow-tree release build lost the 110 592 ring row to OVITO Find rings. The fast tree does not:

> That speed sentence is the main campaign. On the fast tree, cutoff primitive rings of size at most six are faster than OVITO 3.17.2 Find rings at 1 000, 10 648 and 110 592 molecules, on one core and on eight, with the same census of two six-rings per molecule.

## Release rebuild

The d-SEAMS columns in the earlier tables are a debug build (`buildtype=debug`). The watch binary is now `bbdir-release/tests/bench_watch`, compiled with meson `-Dbuildtype=release` (`-O3`, `SEAMS_WATCH_RELEASE`). It refuses to build otherwise. CHILL+ is `seams_chill_plus` (four mutual neighbours; on this crystal that is the 3.5 Å shell). Rings are cutoff primitive rings of at most six. One core, one warm-up, six calls. Raw lines are `benchmarks/h2h-watch-v3.log`. Agreement is 1 wherever both codes returned a result.

This base is `f3e9e394`, not the later main whose release CHILL+ is about 45 ms at 110 592 molecules. The measured release median here is 147 ms, about three times the debug 482 ms on the previous tree, not ten.

| Molecules | Operation | d-SEAMS median (IQR) ms | Other median (IQR) ms |
| --- | --- | --- | --- |
| 1000 | CHILL+ | 1.344 (0.070) | OVITO 6.145 (0.055) |
| 1000 | CHILL+ | 1.344 (0.070) | mdapy 0.914 (0.036) |
| 1000 | rings | 3.044 (0.059) | OVITO 4.586 (0.223) |
| 1000 | rings | 3.044 (0.059) | matscipy 16.260 (0.669) |
| 1000 | rings | 3.044 (0.059) | cycless 3295.845 (162.218) |
| 10648 | CHILL+ | 13.561 (0.112) | OVITO 70.920 (13.741) |
| 10648 | CHILL+ | 13.561 (0.112) | mdapy 7.613 (0.686) |
| 10648 | rings | 45.673 (2.337) | OVITO 44.577 (0.025) |
| 10648 | rings | 45.673 (2.337) | matscipy 1309.019 (20.284) |
| 10648 | rings | 45.673 (2.337) | cycless 35583.400 (233.363) |
| 110592 | CHILL+ | 146.964 (4.363) | OVITO 635.465 (12.525) |
| 110592 | CHILL+ | 146.964 (4.363) | mdapy 63.639 (1.843) |
| 110592 | rings | 2206.976 (5.026) | OVITO 553.569 (1.689) |
| 110592 | rings | 2206.976 (5.026) | matscipy signal 11 |
| 110592 | rings | 2206.976 (5.026) | cycless 385129.916 (3580.078) |

mdapy's cutoff CHILL+ stays faster at every size, with the same cubic label, so the speed-sentence caveat stands on the release build. d-SEAMS CHILL+ is faster than OVITO. Primitive rings are faster than matscipy and cycless where those finish, and the census matches. OVITO Find rings is faster at 110 592 molecules. matscipy still dies with signal 11 at that size.

## Symmetry and brine, 2025–2026

IRA ([Software Impacts 2022](https://doi.org/10.1016/j.simpa.2022.100264)) and SOFI ([JCP 2024](https://doi.org/10.1063/5.0215689)) match a finite cluster to a reference. This tree can link them (`with_ira`) for point-group labels of a cage. They do not classify a thermal ice trajectory. spglib ([Togo, Shinohara and Tanaka, 2024](https://doi.org/10.1080/27660400.2024.2384822)) and [moyo](https://github.com/spglib/moyo) find the space group of a periodic crystal. That complements a bulk ice phase call and does not replace CHILL+ or a ring update. No 2025–2026 journal DOI for moyo turned up.

Brine and electrolyte models, not per-molecule classifiers:

- Blazquez, Sedano and Vega, [JCP 2024](https://doi.org/10.1063/5.0241233). Madrid-2019 ions with TIP4P/Ice. Freezing temperatures of salt solutions, not an order parameter.
- Blazquez, de Lucas, Vega and Gámez, [JCP 2025](https://doi.org/10.1063/5.0267223). A rigid H3O+ added to Madrid-2019.
- Cruz-Sánchez, Blazquez, Vega and Trejos, [Condensed Matter Physics 2026](https://doi.org/10.5488/cmp.29.13601). Perchloric acid on the extended Madrid-2019 field.
- Sedano, Espinosa, Tejedor, Vega and Noya, [JCP 2026](https://doi.org/10.1063/5.0314384). Ice III–water interfacial free energy with TIP4P/Ice. Espinosa and Vega. Not a brine classifier and not an in-run CV.
- Romero-Guzmán, Algaba, Vega and Blas, [JCP 2026](https://doi.org/10.1063/5.0326799). Melting points of TIP4P/2005 and TIP4P/Ice when dispersion uses particle-mesh Ewald.
- Sedano, Vega, Noya and Sanz, [JCP 2025](https://doi.org/10.1063/5.0247832). A two-state water model, TIP4P/2005-Ice.

None of these assign CHILL+, rings or cages during a PLUMED, LAMMPS or OpenMM step. The in-situ claim is unchanged.

## Fast tree, one core and eight

The head-to-head above is `f3e9e394`, where `seams_chill_plus` is the mutual four-neighbour path and the ring driver is `bench_watch`. This section is `perf/chillplus-fast` (`cf8bcc3d`; the requested `a0cc5e6f` is not an object on that branch). CHILL+ is `seams_frame_chill_plus` with the cutoff graph. Rings are `seams_frame_rings` with the same cutoff and a maximum depth of six. The release library is meson `-Dbuildtype=release -Doptimization=3`. One warm-up, six calls, median and type-7 IQR. One core is `taskset -c 0` with one OpenMP thread and one Rayon thread. Eight cores are `taskset -c 0-7`. d-SEAMS reads `OMP_NUM_THREADS` and `RAYON_NUM_THREADS`. mdapy uses `set_num_threads`. OVITO reads `OVITO_THREAD_COUNT` before import. matscipy and cycless have no parallel ring API; they run at the same OpenMP width. Agreement is 1 wherever both codes returned a result: every CHILL+ label is cubic, and the ring census is two six-rings per molecule (2 000, 21 296, 221 184).

| Threads | Molecules | Operation | d-SEAMS median (IQR) ms | Other median (IQR) ms |
| --- | --- | --- | --- | --- |
| 1 | 1000 | CHILL+ | 0.447 (0.037) | mdapy 0.985 (0.028) |
| 1 | 1000 | CHILL+ | 0.447 (0.037) | OVITO 6.122 (0.070) |
| 1 | 1000 | rings | 1.181 (0.075) | mdapy 27.622 (0.473) |
| 1 | 1000 | rings | 1.181 (0.075) | OVITO 4.616 (0.079) |
| 1 | 1000 | rings | 1.181 (0.075) | matscipy 15.734 (0.929) |
| 1 | 1000 | rings | 1.181 (0.075) | cycless 3226.987 (31.319) |
| 1 | 10648 | CHILL+ | 4.918 (0.132) | mdapy 7.652 (0.014) |
| 1 | 10648 | CHILL+ | 4.918 (0.132) | OVITO 64.190 (0.928) |
| 1 | 10648 | rings | 10.563 (0.061) | mdapy 265.497 (4.453) |
| 1 | 10648 | rings | 10.563 (0.061) | OVITO 44.160 (5.085) |
| 1 | 10648 | rings | 10.563 (0.061) | matscipy 1280.353 (16.751) |
| 1 | 10648 | rings | 10.563 (0.061) | cycless 35693.574 (239.725) |
| 1 | 110592 | CHILL+ | 50.614 (1.184) | mdapy 67.816 (0.833) |
| 1 | 110592 | CHILL+ | 50.614 (1.184) | OVITO 626.181 (3.018) |
| 1 | 110592 | rings | 105.412 (0.584) | mdapy 2243.860 (6.125) |
| 1 | 110592 | rings | 105.412 (0.584) | OVITO 596.963 (6.928) |
| 1 | 110592 | rings | 105.412 (0.584) | matscipy signal 11 |
| 1 | 110592 | rings | 105.412 (0.584) | cycless 382279.618 (2827.667) |
| 8 | 1000 | CHILL+ | 0.430 (0.049) | mdapy 23.682 (0.051) |
| 8 | 1000 | CHILL+ | 0.430 (0.049) | OVITO 6.577 (0.067) |
| 8 | 1000 | rings | 0.934 (0.068) | mdapy 46.017 (0.587) |
| 8 | 1000 | rings | 0.934 (0.068) | OVITO 4.937 (0.357) |
| 8 | 1000 | rings | 0.934 (0.068) | matscipy 18.160 (0.800) |
| 8 | 1000 | rings | 0.934 (0.068) | cycless 3224.250 (31.427) |
| 8 | 10648 | CHILL+ | 4.278 (0.088) | mdapy 30.087 (0.174) |
| 8 | 10648 | CHILL+ | 4.278 (0.088) | OVITO 14.656 (0.258) |
| 8 | 10648 | rings | 5.673 (0.407) | mdapy 292.194 (0.419) |
| 8 | 10648 | rings | 5.673 (0.407) | OVITO 40.472 (1.161) |
| 8 | 10648 | rings | 5.673 (0.407) | matscipy 1259.609 (20.215) |
| 8 | 10648 | rings | 5.673 (0.407) | cycless 35561.341 (178.566) |
| 8 | 110592 | CHILL+ | 35.188 (1.966) | mdapy 89.400 (0.604) |
| 8 | 110592 | CHILL+ | 35.188 (1.966) | OVITO 93.037 (2.953) |
| 8 | 110592 | rings | 52.740 (2.533) | mdapy 2624.489 (7.136) |
| 8 | 110592 | rings | 52.740 (2.533) | OVITO 353.332 (20.949) |
| 8 | 110592 | rings | 52.740 (2.533) | matscipy signal 11 |
| 8 | 110592 | rings | 52.740 (2.533) | cycless 383652.478 (7729.991) |

d-SEAMS is the fastest code in every row that both sides finished, including cycless. The 110 592 ring census is the fast-tree one (105.412 ms on one core), not the slow-tree `bench_watch` census (about 2 207 ms). OVITO Find rings at that size is 596.963 ms. mdapy `RingStatistics` (Franzblau shortest-path rings, maximum size 6, cutoff 3.5 Å) agrees and is slower at every size. matscipy still dies with signal 11 at 110 592 molecules. cycless agrees and is about 3 600 times slower at that size. Raw lines are `benchmarks/h2h-watch-v4.csv`. `snakemake related_watch_v4` regenerates them.

## Rings, cages and in-situ, 2026

- Maharjan, Williamson and Fennell, [JPCB 2026](https://doi.org/10.1021/acs.jpcb.5c08791). A ring summation factor on common-angle-pruned hydrogen-bond rings, for liquid, supercooled and crystalline water. No public code turned up, so it is not in the table. It is a ring metric, not a CHILL+ label and not a hexagonal or double-diamond cage.
- [CageIdentification](https://github.com/UnconvRS/CageIdentification) (Adibifard, 2024) counts 5^12, 5^12 6^2, 5^12 6^3 and 5^12 6^4 cages from a GROMACS trajectory. That is a hydrate-cage frequency, not the ice Ic primitive-ring census, so it is not timed.
- [vmd-hydrate-mcp](https://github.com/wjgoarxiv/vmd-hydrate-mcp) (2026) counts clathrate cages and F3/F4 through a local VMD session. VMD is not installed here, so it is not timed.

No 2025–2026 PLUMED, LAMMPS or OpenMM plugin turned up that assigns CHILL+ or cages during the MD step. The in-situ claim remains the LAMMPS compute.

## Symmetry head-to-head, 110592

`59d85155` is not an object on `feat/symmetry`. The kernel timed here is the branch tip `071489584c70e3bfdf224012c53476547175b5df`. The call is `symm::classifyIce` (scale-free 17-atom shell, ceiling 0.18) from a release `-O3` build of that tree. Ideal ice Ic, 110592 molecules. One warm-up, six calls, median and type-7 IQR. One core is `taskset -c 0`. Eight cores are `taskset -c 0-7`. The kernel reads `OMP_NUM_THREADS`. OVITO reads `OVITO_THREAD_COUNT` before import. pyscal3 uses `set_num_threads`.

Every oxygen is ice Ic under the kernel. OVITO PTM with cubic and hexagonal diamond enabled (the default build leaves both off) labels every oxygen cubic diamond, type 6. Turning every PTM template on does the same and costs more. OVITO Identify diamond structure, the extended common-neighbour analysis, labels every oxygen cubic diamond, type 1. pyscal3 `diamond_structure` does the same under the name "cubic diamond". Adaptive common-neighbour analysis, in OVITO and in pyscal3, labels every oxygen other: the first-shell CNA signatures are fcc, hcp, bcc and icosahedral, and a diamond lattice is none of those. Ackland–Jones does the same. pyscal3 `ace` is the Drautz atomic-cluster-expansion descriptor, not a structure label, so it is not a row. freud 3.6.1 has environment motif matching and bond-order parameters and no polyhedral template matching class.

| Threads | Code | Operation | Median (IQR) ms | Cubic ice |
| --- | --- | --- | --- | --- |
| 1 | d-SEAMS | classifyIce | 189.056 (1.380) | 1.000000 |
| 1 | OVITO 3.17.2 | PTM, diamond templates | 2801.560 (61.830) | 1.000000 |
| 1 | OVITO 3.17.2 | PTM, every template | 4012.979 (33.409) | 1.000000 |
| 1 | OVITO 3.17.2 | identify diamond | 136.690 (0.901) | 1.000000 |
| 1 | OVITO 3.17.2 | adaptive CNA | 275.435 (0.918) | 0.000000 |
| 1 | pyscal3 4.1.0 | diamond structure | 99.743 (1.685) | 1.000000 |
| 1 | pyscal3 4.1.0 | adaptive CNA | 382.516 (2.195) | 0.000000 |
| 1 | pyscal3 4.1.0 | Ackland–Jones | 561.481 (3.372) | 0.000000 |
| 8 | d-SEAMS | classifyIce | 35.520 (1.460) | 1.000000 |
| 8 | OVITO 3.17.2 | PTM, diamond templates | 388.972 (5.139) | 1.000000 |
| 8 | OVITO 3.17.2 | PTM, every template | 533.050 (5.525) | 1.000000 |
| 8 | OVITO 3.17.2 | identify diamond | 35.660 (3.581) | 1.000000 |
| 8 | OVITO 3.17.2 | adaptive CNA | 51.267 (1.626) | 0.000000 |
| 8 | pyscal3 4.1.0 | diamond structure | 22.760 (0.132) | 1.000000 |
| 8 | pyscal3 4.1.0 | adaptive CNA | 73.227 (1.323) | 0.000000 |
| 8 | pyscal3 4.1.0 | Ackland–Jones | 165.939 (4.557) | 0.000000 |

On the matched diamond labels, pyscal3 is faster than the kernel at both widths, and OVITO's diamond identifier is faster at one thread. The kernel is faster than both PTM settings at both widths. Raw lines are `benchmarks/h2h-watch-v5.csv`. The first `ptm` lines in `benchmarks/h2h-watch-v5.log` are every template; the `ptm-diamond` lines are cubic and hexagonal diamond only.

## matscipy rings at 110592

matscipy 1.3.1 `ring_statistics` dies with signal 11 at 110592 molecules. The fault is in `py_distances_on_graph` (`c/ring_statistics.cpp`). That wrapper builds an n by n int32 distance map and `std::fill` writes every entry. At 110592 molecules the map is 48.92 GB. The store is a `vmovdqu64` in libc, called from `py_distances_on_graph`, and the process receives SIGSEGV. A NumPy allocation of the same shape returns without touching the pages, so the failure is the fill, not a Python `MemoryError`.

The same fill completes at 44064 molecules (18 by 18 by 17 cells, 7.77 GB) and the full call with `maxlength=7` returns 88128 six-rings, two per molecule, in 34.8 s. The next cell tried, 46512 molecules (19 by 18 by 17, 8.65 GB), dies in the same fill. 39304 molecules (17 cubed cells) also completes, with 78608 six-rings.

## Graph-partitioned cages, 2026

- Pang and Sun, [JPCA 2026](https://doi.org/10.1021/acs.jpca.6c05769). SQQ implements the graph-guided ring overlap workflow: an oxygen or hydrogen-bond graph, chordless rings as faces, and cages grown along exposed edges. It partitions hydrate domains from shared faces. The code is [pimooni/sqq](https://github.com/pimooni/sqq). It names cages and phases, not a PTM or CHILL+ label, so it is not in the symmetry table.
- Chen, Huang, Xu, Li, Song, Wang, Tang, Liu and Deng, [arXiv:2609.37036](https://arxiv.org/abs/2609.37036) (2026). HyTopo classifies water as liquid-like, precursor-like, hydrate-like or ice-like, then searches rings and cages and builds a cage-transition network. Preprint; Crossref has no DOI. No public source turned up, so it is not timed.

## Kernel 59d85155

`59d85155` is on `feat/symmetry`. The call is still `symm::classifyIce`, now with the production ceiling 0.12, the tetrahedral q gate, and the all-eclipsed reject. Release `-O3`. Ideal ice Ic, 110592 molecules, one warm-up, six calls. Every oxygen is ice Ic. pyscal3 diamond structure remains the fastest matched label. The kernel is faster than both PTM settings and, at eight threads, level with OVITO identify diamond.

| Threads | Code | Operation | Median (IQR) ms | Cubic ice |
| --- | --- | --- | --- | --- |
| 1 | d-SEAMS 59d85155 | classifyIce | 166.575 (0.911) | 1.000000 |
| 1 | OVITO 3.17.2 | PTM, diamond templates | 2711.014 (35.722) | 1.000000 |
| 1 | OVITO 3.17.2 | PTM, every template | 3953.667 (6.929) | 1.000000 |
| 1 | OVITO 3.17.2 | identify diamond | 134.629 (1.275) | 1.000000 |
| 1 | OVITO 3.17.2 | adaptive CNA | 272.391 (0.352) | 0.000000 |
| 1 | pyscal3 4.1.0 | diamond structure | 96.051 (0.362) | 1.000000 |
| 1 | pyscal3 4.1.0 | adaptive CNA | 358.423 (4.680) | 0.000000 |
| 1 | pyscal3 4.1.0 | Ackland–Jones | 529.562 (5.703) | 0.000000 |
| 8 | d-SEAMS 59d85155 | classifyIce | 32.868 (0.382) | 1.000000 |
| 8 | OVITO 3.17.2 | PTM, diamond templates | 379.501 (46.228) | 1.000000 |
| 8 | OVITO 3.17.2 | PTM, every template | 593.633 (14.970) | 1.000000 |
| 8 | OVITO 3.17.2 | identify diamond | 36.833 (3.302) | 1.000000 |
| 8 | OVITO 3.17.2 | adaptive CNA | 51.414 (2.428) | 0.000000 |
| 8 | pyscal3 4.1.0 | diamond structure | 25.946 (0.724) | 1.000000 |
| 8 | pyscal3 4.1.0 | adaptive CNA | 71.561 (0.668) | 0.000000 |
| 8 | pyscal3 4.1.0 | Ackland–Jones | 153.488 (0.775) | 0.000000 |

Raw lines are `benchmarks/h2h-watch-v6-ideal.csv`.

The 200 ps frames are not in the tree. `scripts/nucleation_200ps.py` rebuilds them: 12×12×12 ice Ic, 40 ps at 500 K, the 300 central molecules put back, liquid inside 2.4 Å of that core dropped (13553 atoms remain), core fixed, the rest at 220 K for 200 ps. Velocity seeds 91001, 92002, 93003. LAMMPS 22 Jul 2025, `pair_style sw`, timestep 5 fs. Counts below are Ih / Ic / Isd for the kernel, cubic diamond / hexagonal diamond for PTM, and the same two bulk labels for identify diamond and pyscal3. Neighbour shells of the diamond identifiers are in the csv.

| Seed | Threads | Kernel ms (Ih, Ic, Isd) | PTM ms (cubic, hex) | Identify diamond ms (cubic, hex) | pyscal3 diamond ms (cubic, hex) |
| --- | --- | --- | --- | --- | --- |
| 91001 | 1 | 84.099 (57, 275, 1) | 378.401 (288, 39) | 18.721 (247, 37) | 96.381 (247, 37) |
| 91001 | 8 | 12.742 (57, 275, 1) | 60.352 (288, 39) | 6.255 (247, 37) | 20.789 (247, 37) |
| 92002 | 1 | 83.266 (33, 261, 0) | 378.093 (267, 19) | 18.762 (225, 25) | 98.157 (225, 25) |
| 92002 | 8 | 14.584 (33, 261, 0) | 63.615 (267, 19) | 6.531 (225, 25) | 22.274 (225, 25) |
| 93003 | 1 | 84.842 (28, 237, 0) | 379.238 (255, 19) | 18.659 (216, 20) | 96.333 (216, 20) |
| 93003 | 8 | 12.660 (28, 237, 0) | 60.298 (255, 19) | 5.938 (216, 20) | 18.681 (216, 20) |

On these frames the kernel is faster than PTM and slower than OVITO identify diamond. pyscal3's bulk diamond count matches OVITO's. Raw lines are `benchmarks/h2h-watch-v6-frames.csv`.

## SQQ on the 110592 ring row

SQQ 0.5.8 is public ([pimooni/sqq](https://github.com/pimooni/sqq), PyPI). The row is the Python engine, O–O cutoff 0.35 nm, `ring.definition=shortest_path`, sizes 4–6. Cages and order parameters are off. SQQ does not search three-rings. The census is 0, 0, 221184, the same two six-rings per molecule. One frame does not split across `parallel.workers`, so eight cores do not change the time.

| Threads | d-SEAMS rings ms | SQQ median (IQR) ms |
| --- | --- | --- |
| 1 | 105.412 | 18234.926 (590.433) |
| 8 | 52.740 | 17907.744 (430.095) |

The d-SEAMS times are the fast-tree ring row. Raw SQQ lines are `benchmarks/h2h-watch-v6-sqq.csv`.

All 42 journal DOIs in `references/related_work.bib` match Crossref. HyTopo and the metatensor preprint remain without a Crossref DOI.

## Kernel b0aeaf41

`b0aeaf41` adds a clean-diamond path: a second shell on one radius with four staggered bonds is ice Ic, and one eclipsed bond is ice Ih, before the permutation fit. Ideal ice Ic, 110592 molecules, is still entirely ice Ic.

| Threads | Code | Operation | Median (IQR) ms | Cubic ice |
| --- | --- | --- | --- | --- |
| 1 | d-SEAMS b0aeaf41 | classifyIce | 102.752 (7.446) | 1.000000 |
| 1 | OVITO 3.17.2 | PTM, diamond templates | 2721.777 (31.526) | 1.000000 |
| 1 | OVITO 3.17.2 | identify diamond | 135.320 (1.283) | 1.000000 |
| 1 | pyscal3 4.1.0 | diamond structure | 96.876 (1.059) | 1.000000 |
| 8 | d-SEAMS b0aeaf41 | classifyIce | 21.687 (0.386) | 1.000000 |
| 8 | OVITO 3.17.2 | PTM, diamond templates | 371.530 (4.535) | 1.000000 |
| 8 | OVITO 3.17.2 | identify diamond | 40.453 (1.648) | 1.000000 |
| 8 | pyscal3 4.1.0 | diamond structure | 24.472 (1.750) | 1.000000 |

At one thread pyscal3 is ahead of the kernel. At eight threads the kernel is ahead of pyscal3 and of OVITO identify diamond. Adaptive CNA and Ackland–Jones still label every oxygen other. Raw lines are `benchmarks/h2h-watch-v7-ideal.csv`.

The same three 200 ps frames, 13553 atoms. Counts are Ih, Ic, Isd for the kernel and cubic, hexagonal diamond for the others.

| Seed | Threads | Kernel ms (Ih, Ic, Isd) | PTM ms (cubic, hex) | Identify diamond ms (cubic, hex) | pyscal3 ms (cubic, hex) |
| --- | --- | --- | --- | --- | --- |
| 91001 | 1 | 79.158 (57, 278, 1) | 375.379 (288, 39) | 18.717 (247, 37) | 96.228 (247, 37) |
| 91001 | 8 | 12.768 (57, 278, 1) | 57.719 (288, 39) | 5.963 (247, 37) | 21.558 (247, 37) |
| 92002 | 1 | 79.782 (33, 262, 0) | 382.575 (267, 19) | 18.757 (225, 25) | 96.465 (225, 25) |
| 92002 | 8 | 12.409 (33, 262, 0) | 57.607 (267, 19) | 5.851 (225, 25) | 19.323 (225, 25) |
| 93003 | 1 | 81.514 (29, 237, 0) | 375.969 (255, 19) | 18.647 (216, 20) | 96.059 (216, 20) |
| 93003 | 8 | 13.147 (29, 237, 0) | 57.492 (255, 19) | 5.817 (216, 20) | 19.232 (216, 20) |

Raw lines are `benchmarks/h2h-watch-v7-frames.csv`.

## q4, q6 and neighbour lists

Ideal ice Ic. d-SEAMS `seams_frame_steinhardt` uses the four nearest neighbours inside a 5 Å window. On this crystal that shell is the 3.5 Å cutoff shell. Freud, pyscal3 and OVITO use the cutoff directly. OVITO's modifier fills every order from 2 through `max_l` in the timed call. q4 mean is 0.509175 and q6 mean is 0.628539.

| Threads | Order | d-SEAMS ms | freud 3.6.1 ms | pyscal3 4.1.0 ms | OVITO 3.17.2 ms | max abs vs d-SEAMS |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | q4 | 66.990 (0.298) | 263.283 (1.454) | 96.236 (2.077) | 49.989 (1.078) | 7.464e-08, 3.331e-16, 3.331e-16 |
| 1 | q6 | 72.020 (0.621) | 349.392 (0.491) | 108.280 (3.691) | 58.710 (1.554) | 1.416e-07, 7.772e-16, 8.882e-16 |
| 8 | q4 | 27.273 (4.288) | 52.032 (1.747) | 18.128 (1.515) | 11.354 (0.336) | 7.464e-08, 3.331e-16, 3.331e-16 |
| 8 | q6 | 29.327 (7.914) | 61.005 (1.846) | 22.435 (1.489) | 13.113 (0.388) | 1.416e-07, 7.772e-16, 8.882e-16 |

Neighbour lists are the 3.5 Å periodic shell, self excluded. Every code returns 4N directed pairs (42592 and 442368).

| Threads | Molecules | d-SEAMS ms | freud AABBQuery ms | matscipy ms | ASE ms | MDAnalysis ms |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 10648 | 2.469 (0.033) | 9.618 (0.117) | 8.465 (0.028) | 112.474 (3.253) | 60.567 (0.333) |
| 1 | 110592 | 25.778 (0.390) | 111.369 (2.834) | 75.898 (0.129) | 1604.503 (191.021) | 616.682 (1.798) |
| 8 | 10648 | 2.804 (0.050) | 3.604 (0.162) | 8.446 (0.090) | 111.092 (4.157) | 60.761 (0.533) |
| 8 | 110592 | 23.709 (0.190) | 35.486 (1.534) | 76.010 (0.570) | 1555.691 (218.647) | 640.507 (14.143) |

d-SEAMS is the fastest neighbour list in every row. OVITO's Steinhardt modifier is the fastest q4 and q6. Raw lines are `benchmarks/h2h-watch-v7.csv`.

Bigi and others, [arXiv:2508.15704](https://arxiv.org/abs/2508.15704) (2025), present vesin, a cutoff neighbour list used by metatensor. Crossref has no record for that preprint. It is not in the timed neighbour table.

## 12 November 2026

Window 6–12 November. Journal DOIs below were read back from Crossref. The ClavisCage DOI resolves at DataCite; Crossref has no work record for it. An arXiv query for ice, clathrate or hydrate titles, and for CHILL, Steinhardt, ring statistics and bond-orientational order, in `physics.chem-ph`, `cond-mat.soft` and `physics.comp-ph`, submitted 8–12 November, returned no hits. Computer Physics Communications published no article in 8–13 November. The four JOSS papers in that window are Minion, VERNIER, LightCurveLynx and exoatlas.

ClavisCage 0.3.18 ([Go, Mok and Seo](https://doi.org/10.5281/zenodo.24112918), 10 November) is a Rust engine, compiled to WebAssembly and to the Python module `clavis_engine`. One install runs TRACE, HTR, HTR+, GRADE, CHILL+, F3/F4, locally averaged Steinhardt q4/q6, MCG-3 and the tetrahedral order. HTR is [Liu and others, Nanotechnology Reviews 2022](https://doi.org/10.1515/ntrev-2022-0044). HTR+ is [Shi and others, Journal of Physics: Condensed Matter 2024](https://doi.org/10.1088/1361-648X/ad52df). The engine recomputes each frame.

CHILL+ uses the Nguyen and Molinero cuts, staggered below -0.8 and eclipsed on (-0.35, 0.25), on oxygens inside 0.35 nm, and the OVITO class codes. Ideal ice Ic, oxygen bond 2.75 Å. Hydrogens sit on the oxygen, the placement in the engine's own tests; the classifier reads the oxygens. d-SEAMS is `seams_frame_chill_plus` on the cutoff graph, release `-O3`, version 2.13.0. ClavisCage is `analyze_waters`, which parses the JSON water list and then classifies. One warm-up, six calls, median and type-7 IQR. One core is `taskset -c 0`. Eight cores are `taskset -c 0-7`. Raw lines are `benchmarks/h2h-watch-v11.csv`. `snakemake related_watch_v11` regenerates them. Agreement is the fraction of molecules with the same merged class.

| Threads | Molecules | d-SEAMS median (IQR) ms | ClavisCage median (IQR) ms | Cubic ice |
| --- | --- | --- | --- | --- |
| 1 | 1000 | 0.457 (0.052) | 1.471 (0.016) | 1 |
| 1 | 10648 | 4.933 (0.141) | 16.787 (0.081) | 1 |
| 1 | 110592 | 50.430 (0.603) | 176.836 (2.432) | 1 |
| 8 | 1000 | 0.429 (0.042) | 1.476 (0.041) | 1 |
| 8 | 10648 | 4.125 (0.552) | 16.829 (0.227) | 1 |
| 8 | 110592 | 38.508 (3.858) | 175.856 (1.470) | 1 |

Every oxygen is cubic in both codes. Eight cores leave the ClavisCage time where the one-core time is: `chill.rs` has no parallel loop. d-SEAMS is faster in every row.

Other code movement, with no new timing:

- freud's default branch, after 3.6.1, computes `Hexatic` and `StaticStructureFactorDirect` with Euler's formula where the real part of the exponent is zero. The change is not on the Steinhardt `q_l` path, and it is not in the 22 September release. The q4 and q6 rows stand as recorded.
- vesin 0.6.2 is still the release. Commits on 9–10 November add CUDA streams in `vesin_neighbors` and a LAMMPS column in that project's benchmarks. The CPU neighbour table was left as recorded.
- mdapy 2.0.0 was already timed. Commits on 12 November add phonon group velocities, thermal displacements and conductivity.
- SQQ 0.5.8, OVITO 3.17.2, pyscal3 4.1.0, matscipy 1.3.1, fastatomstruct 2.2.0, cycless 0.7, GenIce 2.2.13.3 and PLUMED 2.10.2 have no newer release. PLUMED's commits since 8 November are documentation. No new ice collective variable turned up.

Papers in the window, none of which return hexagonal or double-diamond cages:

- Yerby, Chen, Bertho, Zhu, Francisco and Yang, [Science Advances 2026](https://doi.org/10.1126/sciadv.aej1123). Molecular simulations of trace gases at the ice interface. The gases adsorb, slow the ice front, and suppress ice-like ordering in the quasi-liquid layer.
- Villacampa and others, [JACS 2026](https://doi.org/10.1021/jacs.6c14020). Antifreeze glycoproteins bind specific crystallographic planes. Synthesis, NMR, activity and molecular dynamics. A plane-recognition model.
- Najm Alsayed, Pabst, Cassone, Hassanali and Novelli, [JCP 2026](https://doi.org/10.1063/5.0362915). The population of waters with unequal donated and accepted hydrogen bonds tracks a terahertz Drude–Smith weight.
- Kocherbitov, [Scientific Reports 2026](https://doi.org/10.1038/s41598-026-79646-0). TIP4P/2005 cooling scans. The glass transition sits on the same temperature range as the high-density to low-density structural evolution.
- Blossey, [JCP 2026](https://doi.org/10.1063/5.0358448). A short account of Onsager–Dupuis structural polarization in ice, and a fourth-order equation for the electrostatic potential.

## Claims after 12 November

1. **First code to return a full ring search from a local update.** Stands. ClavisCage, the freud hexatic change, and the papers above recompute a frame or do not search rings. None describes keeping rings whose bonds did not change.

2. **Faster than OVITO, freud, pyscal, matscipy and the CHILL library in every matched comparison.** Stands for that list. ClavisCage is outside the list and is slower than d-SEAMS 2.13.0 cutoff CHILL+ at 1 000, 10 648 and 110 592 molecules, on one core and on eight. The mdapy caveat already recorded for cutoff CHILL+ is unchanged by this engine.

3. **Ring-statistics codes stop at counts, while d-SEAMS builds cages.** Still false if read as "no other code builds cages." ClavisCage 0.3.18 now ships TRACE, HTR, HTR+ and GRADE together, plus CHILL+. A revision should name it next to TRACE. Cage counts were not timed: this crystal has no clathrate cages.

4. **OVITO, freud and pyscal have no nanotube-ice class.** Stands. Nothing this week adds a prism-block label.

5. **A LAMMPS compute assigns CHILL+ during the run.** Stands. ClavisCage runs in the browser and as a Python wheel after the trajectory exists. PLUMED's recent commits are documentation. No surveyed source assigns CHILL+ or cages from inside a PLUMED, LAMMPS or OpenMM step.

6. **CHILL+ follows Nguyen and Molinero, with the OVITO modifier as a second option.** ClavisCage uses those cuts and the OVITO codes. On these frames the merged label agrees for every molecule, and the label is cubic. d-SEAMS is the faster of the two calls. The ClavisCage median includes the JSON parse inside `analyze_waters`.

7. **Cage labels come from the topology of the network, and local classifiers answer a different question.** Stands. The five papers are interface adsorption, plane recognition, a hydrogen-bond population, a two-state glass scan, and a polarization theory. None returns the hexagonal-cage or double-diamond-cage count.

## After 06:15 UTC, 12 November 2026

The morning sweep closed at 00:15 CST. This section is what appeared after 2026-11-12T06:15:00Z.

arXiv `submittedDate` from 2026-11-11 00:00 through 2026-11-12 23:59 returned no works. The newest records in `physics.chem-ph`, `cond-mat.soft`, `physics.comp-ph` and `cond-mat.mtrl-sci` are 10 November. Crossref titles created on 12 November contain no clathrate, Steinhardt or CHILL paper. The ice titles are geophysical, food or classroom notes. The hydrate titles after 06:15 UTC are a CO2 sequestration experiment, creep of hydrate-bearing sediment, calcium silicate hydrate, and a medicinal hydrate. None identifies ice or cages in a trajectory. No new citation.

[rings-statistics-matter 0.2.0](https://pypi.org/project/rings-statistics-matter/0.2.0/) was uploaded at 2026-11-12T11:48:09Z. The repository is [MorrowChem/RingsStatisticsMatter](https://github.com/MorrowChem/RingsStatisticsMatter) (created 2023; the morning note said it had no Python library). The default backend is the Rust kernel. `ring_statistics` takes an ASE atoms object, a cutoff, and `maxlvl`. Neighbours come from vesin 0.6.2 inside the call. The kernel tries each node closer than `maxlvl / 2` as the midpoint of a ring, so a six-ring, whose midpoint sits at distance 3, needs `maxlvl` 8. `maxlvl` 7 returns no rings on this ice. The package text says rings of up to `maxlvl - 1`, which would have chosen 7. The algorithm is Yuan and Cormack, [Comput. Mater. Sci. 24, 343 (2002)](https://doi.org/10.1016/s0927-0256(01)00256-7). It searches the whole frame. It returns the count and the atom indices of each ring. It does not build a hexagonal or double-diamond cage.

Ideal ice Ic, oxygen bond 2.75 Å, cutoff 3.5 Å, `maxlvl` 8, `supercell` off, deduplicated rings. d-SEAMS is `seams_frame_rings` on the cutoff graph, depth 6, release `-O3`, version 2.13.0. One warm-up, six calls, median and type-7 IQR. One core is `taskset -c 0`. Eight cores are `taskset -c 0-7`. The Rust call is given `num_threads` of that width. Raw lines are `benchmarks/h2h-watch-v12.csv`. `snakemake related_watch_v12` regenerates them. Agreement is 1 when the census of rings of size 3 through 6 matches. It matches in every row: 0, 0, 0, and two six-rings per molecule (2 000, 21 296, 221 184).

| Threads | Molecules | d-SEAMS median (IQR) ms | rings-statistics-matter median (IQR) ms |
| --- | --- | --- | --- |
| 1 | 1000 | 1.132 (0.008) | 12.034 (0.138) |
| 1 | 10648 | 10.269 (0.099) | 123.192 (0.233) |
| 1 | 110592 | 107.461 (0.855) | 1073.553 (16.738) |
| 8 | 1000 | 4.017 (1.053) | 3.555 (0.431) |
| 8 | 10648 | 5.498 (0.171) | 30.724 (1.153) |
| 8 | 110592 | 56.196 (4.606) | 258.706 (16.674) |

At one core d-SEAMS is faster at every size, by a factor of ten to twelve. At eight cores it is faster at 10 648 and 110 592 molecules. At eight cores and 1 000 molecules the rival's median is lower. That d-SEAMS batch had not settled: the six samples run from 0.771 ms to 4.242 ms. The same d-SEAMS call, measured again in the next rows, has median 0.753 ms (IQR 0.197).

The default `supercell=True` repeats a cell narrower than `2 * maxlvl * cutoff` = 56 Å. The 1 000-molecule side is 31.754 Å, and the repeat is 2×2×2. Counts are mapped back onto the original cell and still match. The 10 648-molecule and 110 592-molecule cells are wider than 56 Å, so the default does not repeat them.

| Threads | d-SEAMS median (IQR) ms | default supercell median (IQR) ms |
| --- | --- | --- |
| 1 | 1.121 (0.016) | 88.356 (2.047) |
| 8 | 0.753 (0.197) | 24.090 (1.604) |

Other code, with no new timing:

- mdapy 2.0.0 is still the release. After 06:15 UTC, master gained dislocation-analysis output, not a change to cutoff CHILL+: [pull 75](https://github.com/mushroomfire/mdapy/pull/75) writes an OVITO crystal-analysis file, [pull 76](https://github.com/mushroomfire/mdapy/pull/76) reads one back, [pull 77](https://github.com/mushroomfire/mdapy/pull/77) accepts hcp with a low c/a, and [pull 78](https://github.com/mushroomfire/mdapy/pull/78) adds line densities and a selected-atom system. The CHILL+ rows stand.
- PLUMED 2.10.2 is still the release. [Pull 1599](https://github.com/plumed/plumed2/pull/1599), merged 09:56 UTC, gives `LOAD` an `INCLUDE_PATH`. It is not an ice collective variable.
- fastatomstruct 2.2.0 is still the release. Two merges at 06:22 UTC add a test for `parse_composition` ([!277](https://zivgitlab.uni-muenster.de/ag-salinga/fastatomstruct/-/merge_requests/277)) and the missing autodoc entry for the existing `MACECalculatorASE` ([!276](https://zivgitlab.uni-muenster.de/ag-salinga/fastatomstruct/-/merge_requests/276)). The `q_6` and King-ring rows stand.
- freud 3.6.1, OVITO 3.17.2, MDAnalysis 2.11.0, pyscal3 4.1.0, GenIce 2.2.13.3, vesin 0.6.2, SQQ 0.5.8, ClavisCage 0.3.18, matscipy 1.3.1 and cycless 0.7 have no commit after 06:15 UTC and no newer release. The `stuko/ovito` GitLab project path returns 404; the PyPI upload of OVITO 3.17.2 remains 2 November.

## Claims after the afternoon sweep

1. **First code to return a full ring search from a local update.** Stands. rings-statistics-matter 0.2.0 searches the whole frame. Nothing in this window keeps rings whose bonds did not change.

2. **Faster than OVITO, freud, pyscal, matscipy and the CHILL library in every matched comparison.** Stands for that list. rings-statistics-matter is outside it. On the matched rows here, d-SEAMS 2.13.0 is faster except at eight cores and 1 000 molecules, where the first six-call median is slower than the rival and the repeated measurement of the same call is faster. The mdapy CHILL+ caveat is unchanged: the new mdapy commits are dislocation output.

3. **Ring-statistics codes stop at counts, while d-SEAMS builds cages.** The new package returns the atom indices of each ring as well as the counts. It still does not assemble hexagonal or double-diamond cages. The morning sentence that it had no Python library is no longer true.

4. **OVITO, freud and pyscal have no nanotube-ice class.** Stands.

5. **A LAMMPS compute assigns CHILL+ during the run.** Stands. The PLUMED change is a `LOAD` search path. mdapy's new dislocation output is post-processing on master. fastatomstruct's merges are a test and a documentation entry. No surveyed source assigns CHILL+ or cages from inside a PLUMED, LAMMPS or OpenMM step.

6. **CHILL+ follows Nguyen and Molinero, with the OVITO modifier as a second option.** Unchanged. No new CHILL+ implementation appeared after 06:15 UTC.

7. **Cage labels come from the topology of the network, and local classifiers answer a different question.** Stands. No new classifier from this window returns the hexagonal-cage or double-diamond-cage count.
