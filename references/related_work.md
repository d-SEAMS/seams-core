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
| [RingsStatisticsMatter.jl](https://github.com/MorrowChem/RingsStatisticsMatter.jl) | Julia primitive-ring counts | Competes on the census. No Python or C++ API. |
| [Palmer CHILL+ library](https://github.com/palmergroupUH/Order_Parameters_Library) | Fortran CHILL+ | Already matched. |
| [mdevaluate `chill.py`](https://gitea.pkm.physik.tu-darmstadt.de/IPKM/mdevaluate) | Staggered and eclipsed bond counts | A non-periodic KDTree. Not the periodic cutoff list, so it was not timed. |
| [Clavis CHILL+](https://docs.claviscage.app/analysis/chill-plus) | CHILL+ on loaded GROMACS and LAMMPS frames | Post-processing. Complements the LAMMPS compute. |
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

Speed sentence. Keep the five named codes, and do not let mdapy be swept into them:

> In all matched comparisons of the main campaign, d-SEAMS 2.0 is faster than OVITO, freud, pyscal, matscipy and a CHILL library; mdapy 2.0 is outside that list, and its cutoff CHILL+ is faster on ideal ice Ic.

Hydrate cages. The following sentence currently says d-SEAMS is the code that assembles rings into cages:

> Unlike the ring-count codes, d-SEAMS assembles the hexagonal and double-diamond cages of ice, while TRACE assembles hydrate cages, including incomplete cages and cages coordinated to one additive.

A further scope line, from the release ring timings below:

> That speed sentence is the main campaign. On this tree's release build, OVITO 3.17.2 Find rings is faster than d-SEAMS at 110 592 molecules and tied at 10 648, with the same six-ring census; d-SEAMS is faster at 1 000.

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
