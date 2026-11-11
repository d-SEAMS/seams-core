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
