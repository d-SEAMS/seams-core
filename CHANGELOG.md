# Changelog

All notable changes to this project are documented in this file.

The format is [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Unreleased notes live in [`changelog.d/`](changelog.d/) and are assembled
by [towncrier](https://towncrier.readthedocs.io/).

<!-- towncrier release notes start -->

## [2.14.0] - 2026-11-13

### Added

- A nucleation-pathway census reports primitive rings by size, both CHILL+ rules, Lechner-Dellago q6-bar, and diamond symmetry labels on every frame.
  The CHILL+ staggered and eclipsed edges are unchanged.
  Largest-cluster size and the Ic/Ih fraction are recorded for each definition.
  `nucleation-pathway` writes the table. `input/nucleation/regenerate.py` rebuilds the short mW quenches outside the work tree.
- A running molecular dynamics calculation can call the in-situ engine on
  each frame it chooses.
  The engine reports the primitive ring census.
  It also reports the published common neighbour ice labelling (CHILL)
  labels.
  The cage flags cover hexagonal and double-diamond cages.
  LAMMPS, GROMACS, OpenMM, i-PI, and the Atomic Simulation Environment
  attach through optional adapters.
  -Dplumed=enabled now errors. Live analysis does not use PLUMED.
- Add Rodger F4 and host-only Steinhardt l=12 so the Zeron q3/q12 hydrate pair is a per-atom field.
- Add TUM/rings stacking planes (HC-basal and DDC-equatorial) next to the CHILL+ I_sd molecule-bin reference.
- Add a C ABI seams_chill_plus for native and emscripten builds.
- Add a pinned head-to-head of cutoff neighbours, q4, q6, diamond, rings and local symmetry against mdapy, pyscal3, vesin, IRA and SQQ.
- Add extras/lammps compute dseams that writes seams_chill_plus per atom.
- Add optional IRA, SOFI and spglib symmetry reports for ice shells, hydrate cages and rings.
- Add pinned head-to-heads of 5^12 cages against ClavisCage TRACE and of q4, q6 and coordination against PLUMED.
- Assign guests by ray-parity inside fan-triangulated cage faces, and count ions per ice cluster on the oxygen graph.
- Bin CHILL+ cubic and hexagonal molecules into basal layers and emit cubicity Phi_c plus the H/C stacking string.
- Build CHILL+, cages, seeded ions and k-NN on an explicit water-type mask (`--water-types`) so substrate and ions never enter the four-neighbour list. The flag loads the mixed dump; ice and cage counts are the water types only.
- C API: `seams_last_error_message`, `seams_system_natoms`, `seams_neighbours_nrows`, `seams_rings_nrings`, and the length-checked `seams_chill_plus_n`, `seams_q4q6_n`, `seams_rings_labels_n`, and `seams_rings_histogram_n`. Every public declaration in `seams.h` now has a Doxygen block with parameters, ownership, thread safety, and an example.
- C API: `seams_system_set_positions` overwrites the coordinates, and optionally the cell, of a stored frame, and `seams_system_box`, `seams_system_positions`, `seams_system_types`, and `seams_system_ids` copy the stored frame back out through a const handle. `seams_system_set_atoms` reuses the storage of the frame it replaces, and keeps the id index when the ids are 1 to n, so loading 110592 atoms costs about the same as building the frame in C++.
- Direct-coexistence mW growth is scored on basal and secondary-prism slabs of 8064 molecules at 1 bar.
  The 2015 CHILL+ thresholds are unchanged.
- Each OpenMP thread owns a block of the linkcell weighted cell cut and writes the neighbour rows of those cells. The halo is painted in that owned walk. On 884736 molecules, eight threads build the neighbour rows in 5.72 times the one-thread time.
- Each language front reports rings, cages, signature cages and symmetry, and the binding CI runs those calls.
- Enumerate cups and incomplete cages on the ring graph; closed signature counts stay the closed path.
- Fortran callers use a `seams` module of derived types with final procedures, allocatable results, and an optional `stat`/`errmsg` pair, and Lua callers get a Fennel module, a rockspec, and LDoc comments. The numbers match the existing tests.
- Hexagonal channels are stacked six-ring prisms, not a raw six-ring count. Ice XXI (Lee et al. 2026, Z=152 BCT) also requires a tetrahedral 3.5 A graph and a primitive six-ring. Hydrogen MSD uses the minimum image. Glass labels use local-density windows for ice/LDA/MDA/HDA. `compute dseams` ships with `water.data` and `pair_style zero`.
- Ice XXI library matches the Lee et al. 2026 I-42d oxygen cell (SI CIF, Z=152).
- Julia, R, and C++ call the analyses `seams.h` now exports. Julia uses `ccall` on `libseams`. R keeps the `libseams_lang` frame API and adds a C shim. `seams.hpp` is a header-only RAII wrapper whose exceptions carry `seams_handle_error`.
- Live analysis for LAMMPS, GROMACS, and OpenMM calls the library.
  The LAMMPS plugin and the in-tree package share one fix.
  GROMACS uses a module inside mdrun.
  OpenMM uses a reporter.
  Chemfiles reads the saved trajectory.
- Local q4 and q6 on the cutoff shell use the solid-harmonic kernel from `seams steinhardt`, `seams_q4_q6`, and the Python, R, and Julia `q4_q6` fronts.
- Name 51264, 51268 and sH cages, and report a per-cage occupancy histogram.
- One sitting pairs ice Ic at 884736 molecules with freud, mdapy, vesin, OVITO PTM and CHILL+, and rings-statistics-matter, at 1, 2, 4 and 8 threads. The detail column records the agreement check on that frame. The rings row is the cutoff list and the ring search together. The pool-audit rings kernel is ringNetwork alone.
- One sitting retimes neighbours, q4, q6, CHILL+, rings, symmetry and hydrate cages at 110592 and 884736 molecules on one thread and on eight. The raw lines are benchmarks/rc3-sitting.csv.
- OpenMP target offload of the TUM ice score: hop-bound primitive six-rings and HC/DDC cage affiliation. `SEAMS_OFFLOAD=1` cage counts match the host on mW cubic, including `seams cages --graph seeded` (union 4-NN graph). Device CHILL+ is not this path.
- Orthorhombic frames are cut by linkcell's weighted recursive coordinate bisection of the cell grid. The halo is the stencil painted while walking the cells a rank owns, and the same cut is reused on the next frame of that size. Neighbour rows, CHILL+, q6 and primitive rings keep their scratch across that repeat, and the ring search answers a length-2 shortcut that passes below the source.
- Python, Julia, Rust and Fortran call the C API. Python reads NumPy arrays in place and accepts ASE and MDAnalysis frames. Julia uses ccall. Rust splits the raw FFI from the safe crate. Fortran binds the same symbols with iso_c_binding.
- R (`bindings/r/dseams`) and Julia (`bindings/julia/Dseams`) call
  `libseams_lang`. Python stays `pydseamslib`. The executable tutorials
  under `docs/tutorials` import `pydseams`.
- Saved trajectories other than CON are read with chemfiles, including xyz, a LAMMPS dump, gro, pdb, dcd, xtc and Amber NetCDF. A selection expression names water oxygen and hydrogen. The dependency is optional and on by default.
- TIP4P/2005 ring updates at 0.01, 0.1 and 1 ps stay inside the hop ball and match a full search. The strict seeded graph is 8.3, 3.3 and 2.6 times a rebuild.
- The 884736-molecule sitting adds pyscal3 neighbours and q6, PLUMED q6 with PLUMED_NUM_THREADS set to the row width, and OVITO coordination. Those rows agree with d-SEAMS on the same ice Ic frame. scripts/plot_sitting_large_time.py draws time against threads and does not store the PDF.
- The PLUMED SEAMS collective variable reports q6, the CHILL+ ice fraction and the six-ring fraction, and a restraint can bias q6.
- The Rust, Fortran, and Python ctypes bindings call the analyses `seams.h` now exports, including frame reading, CHILL, Steinhardt parameters, cages, neighbour lists, the RDF, hydrogen bonds, density, the fingerprint, prisms, and the monolayer.
- The binding-overhead session times each language wrapper of read, k-nearest neighbours, q6 with q6bar, CHILL, and primitive rings on ice Ih, and records the median and interquartile range.
- The ice label is the residual of the seventeen-atom shell after the optimal rotation and scale. The shell is the four bonded neighbours and each of their bonded neighbours.
- The ice label takes an eclipsed-cosine cutoff, and cage, ring and shell naming runs in parallel.
- The installed C library reads frames from arrays and files. Callers can build neighbour lists and primitive rings, evaluate q3, q4, q6, and q6-bar, and run CHILL, CHILL+, cages, the shell classifier, and the brine front. `seams.pc` and `seams-config.cmake` ship with the library, and `examples/c` builds against that install.
- The mW front speed was measured again under the other readings of L and tau, under Lechner-Dellago q6, on the three faces, and at a 5 fs step.
  The one-front speed is a quarter of Espinosa et al. Table I at both 240 K and 260 K, and it matches the primary-prism markers in Shi et al., J. Chem. Phys. 162, 054714 (2025), Figure 3.
  Equation 2 is unchanged.
- The mW one-front speed was measured on the basal, primary-prism and secondary-prism faces from 220 K to 265 K.
  A Wilson-Frenkel curve is fit to the count-slope estimator on each face.
- The secondary-prism speed is reported under equation 2 and under the solid-count slope, on the 4000-molecule cell.
  The 2015 CHILL+ thresholds are unchanged.
- Under mpiexec, an MPI build of `seams` shares a LAMMPS dump's `--frame`/`--last` frames among the ranks in serpentine rounds (`sinp::forEachLammpsFrame` takes a part and a part count), and rank 0 prints every frame in order. Other readers stay on rank 0, `--per-atom` needs a single rank, ranks that differ in arguments, directory or `SEAMS_` settings are refused, `OMP_NUM_THREADS` and `RAYON_NUM_THREADS` default to 1 on a launch of several ranks, and a run without a launcher, or with `SEAMS_MPI=0`, never starts MPI. `chill::setSteinhardtAtomSplit(false)` keeps `steinhardtQl` from splitting atoms over ranks that hold different frames, and a dump frame cut short keeps the atoms it read, with its warning on stderr. With Open MPI's libfabric probe skipped (`--mca btl ^ofi`), eight ranks match eight `--jobs` workers on one node for `cages` over 32 frames of 32768-atom ice (4.2 against 4.3 s), and keep load balance 0.99 where frame cost steps or alternates along the run (2.6 s against 2.8 and 3.1 s for `--jobs`; contiguous blocks reach 0.59 on the step). On two four-core core sets standing in for nodes, with MPI over TCP only, eight ranks take 4.4 s against 8.7 s for threads on one, at POP parallel efficiency 0.98 (release builds).
- Update topology keys only on the hop-ball of atoms whose neighbourhood changed.
- `-Dcapi_sanitize=address,undefined` builds the C API suite under ASan and UBSan at once, and every UBSan build also checks float-to-int overflow, which GCC leaves out of `-fsanitize=undefined`.
- `bench_rings_h2h` times primitive rings against OVITO FindRings with the cutoff graph outside the timer. A small-graph test checks that census against a brute-force cycle list.
- `capi_abi` compares libseams type by type with the 2.14.0 ABI 1 corpus through libabigail's `abidiff` when it is installed, and passes only added functions. The `abi` pixi environment provides `abidiff`, and the `capi-abi` CI job uses it. Without `abidiff`, or with one too old to read the corpus, that step is skipped.
- `capi_reload_fuzzer` drives `seams_system_set_atoms`, `seams_system_set_positions` and `seams_system_stage_positions` with random atom counts, NULL pointers, wrong lengths, NaN and infinite values in both reload modes, and reads the frame back after every call. Under clang the libFuzzer harnesses are built with `-fsanitize=fuzzer` and the `capi_sanitize` sanitizer.
- `interface_order wc` reports a Willard-Chandler width on the crystal centre of mass.
  `input/interface/campaign.py` grows the shipped mW slab at several temperatures and seeds.
  The 2015 CHILL+ thresholds are unchanged.
- `interface_order` compares the Nguyen and Molinero 2015 CHILL+ rule with the OVITO 3.17.2 rule on one frame or a short trajectory.
  The thresholds are unchanged.
  A square, pentagonal or hexagonal nanotube is labelled by its prism blocks, and a honeycomb sheet by its polygon or bilayer rings.
- `require("dseams")` is a Lua 5.4 and LuaJIT module over `seams.h`, with userdata handles and a Fennel wrapper on the same module. Reading, CHILL, Steinhardt parameters, cages, neighbour lists, the RDF, hydrogen bonds, density, the fingerprint, prisms, and the monolayer call that header.
- `seams::domain` splits a periodic frame over ranks: cells ordered along a Hilbert curve are cut into runs of equal atom count, and each rank also holds every atom within a halo of its own. `domain::rings` finds the primitive rings whose lowest-indexed member a rank owns from its share alone, `primitive::ringNetwork` takes a mask of the sources to enumerate from, and `domain::gatherRings` collects every rank's rings on every MPI rank in the order `ringNetwork` lists them. On one eight-core node threads still win (a 65536-atom frame takes 58 ms on one rank of eight threads and 67 ms on eight single-thread ranks), and every rank rebuilding the whole list costs 72 ms of 251 at 262144 atoms on eight ranks; `bench_strong` reports it as `domain/ms` in MPI builds, and CI runs the gather under `mpiexec`.
- science::layersOf bins published CHILL+, the shell label and q6 along the interface normal. science::criticalSize is the mean size of the band around the smallest largest cluster. That number is an upper bound when the minimum falls in the first two frames. It is not a committor.
- seams rings prints the primitive-ring census by size. seams prisms counts perfect and imperfect prisms of size 3 to 6. seams cluster prints the largest domain, the subset, and whether it percolates.
- seams.h is ABI 1. A frame carries an orthorhombic or triclinic cell, and the library returns neighbours, both CHILL+ rules, q4 and q6, primitive rings, and ring symmetry labels.
- seams_frame_q4_q6_chill and the Python wrapper write q4, q6, q4-bar and q6-bar together with the published CHILL+ labels.

### Changed

- A threaded cutoff neighbour list on an orthorhombic cell uses the scaling cell list instead of linkcell's buffered walk. The rows stay the same sorted cutoff set, so the ring search matches the one-thread list.
- Below half the narrowest face separation the cutoff list skips minimage's pair reduction, which only sorted there, and sorts its rows on every thread; `getNewNeighbourListByIndex` and `neighListO` fill their rows in parallel, and `neighListO` takes the cutoff list first at every thread count. On `bench_strong`'s 65536-atom frame the index list goes from 25.8 to 13.9 ms on one thread and from 22.7 to 7.8 ms on eight (release builds).
- C API: `seams_system_set_atoms` and `seams_system_set_positions` check the coordinates with a vector scan that runs from the end of the input, so the write pass that follows finds the start of it still in cache. A default reload in a trajectory loop costs 2 to 5% more than writing the coordinates into the C++ frame directly at 4096 and 110592 molecules, against 13 to 15% before, and about 12% at 884736.
- C API: `seams_system_stage_positions` makes `seams_system_set_positions` check and write the coordinates in one pass, into a second set of per-atom records that it swaps in only when every value is finite, at about 72 bytes per atom. By default a reload checks every coordinate first and writes in place, with no extra memory; either way a failed reload keeps the old frame.
- CHILL and CHILL+ default to the mutual four-nearest graph. vesin still owns cutoff pairs; linkcell still owns k-nearest.
- CHILL and CHILL+ use one default neighbour list on every front: the cutoff list (pairs closer than `cutoff`, 3.5 Å by default). The mutual four-nearest list is `graph = "mutual"` under that name everywhere: `seams chill` and `seams chill-plus --graph mutual`, R `chill_plus(graph = "mutual")`, Julia `chill_plus(; graph = "mutual")`, pydseamslib `Frame.chill_plus(graph="mutual")`, the C ABI `seams_chill_plus_graph`, and the LAMMPS `compute dseams ... graph mutual`. `seams chill`, R, Julia, `seams_chill_plus`, and `compute dseams` used the mutual list by default before. `tests/chill_plus_fronts.py` checks that the CLI, Python, R, and Julia agree on `exampleTraj.lammpstrj` with both lists.
- CON frames are read with readcon-core's frame iterator. Lengths and angles become the minimage cell, and that cell is the dump box linkcell already uses. A readcon-db corpus, including a sharded campaign root, is selected and decoded with the same reader.
- Cutoff neighbour lists (`neighListPair`, `halfNeighList`, `SkinNeighborList`, and `neighListO` / `getNewNeighbourListByIndex` below 2048 atoms) come from linkcell's `lc_pairs_within` with cells of about 20 atoms, and vesin is the fallback. `halfNeighList` on 8192 atoms at water density with a 3.5 Å cutoff goes from 8.9 to 4.7 ms on one thread and 4.0 ms on eight. The first linkcell cutoff list keeps freed pages on the heap for the rest of the process (`M_TRIM_THRESHOLD`, `M_MMAP_THRESHOLD`).
- Cutoff neighbour lists use linkcell pairs_within with cells of about 20 atoms rather than one cutoff, ahead of the threaded cell rows. The neighbour list still keeps one index per atom. k-nearest stays knearest.
- Document the two LAMMPS type-filtered readers: ~readLammpsTrjO~
  flags ~inSlice~ and keeps every atom of the type;
  ~readLammpsTrjreduced~ drops atoms outside the slice. An axis with
  ~lo == hi~ is unconstrained.
- Inside the one-image ball the in-plane RDF walks Rapaport cell pairs: cells span the cutoff across each face separation, a half stencil visits each neighbouring pair of cells once, and one lattice shift serves the pair. It replaces the packed grid and the vesin list; 8192 atoms at a 6 Å cutoff take 2.2 ms on one thread and 0.4 ms on eight.
- Interior cells of the one-thread orthorhombic cutoff skip the periodic-shift branches. A parallel ring search writes each neighbour row, already sorted, into one array, so the walk after the cell list matches the walk after linkcell.
- MPI ring analyses can reduce ring counts and per-atom ring types with `seams::domain::reduceRings`, which keeps each rank on its own domain and sends only owned atoms to the root instead of gathering every ring on every rank.
- Neighbour, order-parameter, ring and shell calls keep one worker pool and leave the OpenMP team running between a Rayon search and the next region. A per-atom walk under 8192 atoms stays at four threads. Ring sources still take one thread per eighty sources. The one-thread and several-thread results match.
- Neighbour-list q3 uses the one-bond kernel. The eight-bond batch stays on degrees 4 and 6.
- Neighbours, CHILL, CHILL+, rings and cages can keep water oxygens by atom type or by molecule, and leave ions out of that graph.
  `seams --mol water` uses the atom that appears once in a multi-atom molecule.
  `seams --mol 1,5-12` keeps those molecule ids.
  A type or molecule that is not in the frame is an error.
  On a LAMMPS dump with multi-atom molecules and no `--type`, the same choice replaces the old preference for type 2, which labelled hydrogens or sodium.
  `compute dseams` takes `type I` and `mol water`. A mixed frame that sets neither is an error.
  CHILL+ then reports `ion-ice`, `ion-liquid` and `ion-front` when ions remain in the frame.
- Pair distances use the fractional wrap below half the narrowest face separation, and call the Euclidean image only past that. The in-plane RDF bins the cell-list distance inside that ball, and the direct loop runs in parallel.
- Release builds define NDEBUG, so the libstdc++ assertions that meson adds when b_ndebug is false no longer run in them. The closed-cage search on the sI 2x2x2 crystal takes 13.2 ms instead of 16.0 ms, with the same output.
- The C API version macros come from the project version at build time. This release is 2.14.0.
- The SOAP power spectrum evaluates each neighbour's spherical harmonics once and reuses them for every radial function, since its basis is a product of the two; the values are unchanged. On 8192 atoms at water density with nMax 8, lMax 6 and a 6 Å cutoff, `soapSpectrumAll` goes from 831 to 199 ms on one thread (release builds).
- The carried ring update sorts ring indices and merges the kept rings in parallel, and it keeps dropped ring buffers instead of freeing them on the update. The orthorhombic cutoff list counts each pair once, a z layer at a time, so several threads beat the single-thread walk. One-thread and several-thread results still match.
- The clathrate face-edge window is `SEAMS_CAGE_EDGE` and `--cage-edge`, default `2.65,3.15` in the coordinate units of the file. `none` turns the window off, and the cage bound is then the face census alone. `seams_frame_signature_edge` takes the same string. Python is `frame.cages_by_signature("512", edge="none")` and Lua is `dseams.cages_by_signature(frame, "512", {edge = "none"})`.
- The closed-cage search skips a ring that shares two consecutive edges with a chosen face or passes through a vertex whose faces already close, and keeps the vertex and edge counts inside V - E + F = 2. No cage is lost. The search runs on dense edge ids, adds the forced faces and picks the open edge with the fewest rings in one pass, and lists the cages of a signature by their sorted faces.
- The cutoff list now covers restricted triclinic boxes.
  A half list keeps each pair on the lower index.
  A skin reuses the wide list while every atom stays inside the Verlet shell.
  CHILL+, the ring census, and the template overlay read that CSR.
- The cutoff lists bin linkcell's pairs by row on every thread: the pairs are dealt into blocks of rows, four per thread, and one thread counts, places and sorts each block, so no pass waits on shared counters and the rows do not depend on the thread count. linkcell fills a buffer through `lc_pairs_within_rows`, where `linkcell::pairs_within` zeroes its vector on one thread first, and `neighListO` allocates each row once, on the thread that fills it. On a jittered 65536-atom lattice at water density the cutoff list goes from 17.6 to 14.5 ms on one thread, and at 262144 atoms from 67.7 to 57.6 ms on one thread and from 35.1 to 15.2 ms on eight (medians of three interleaved rounds, release builds, both on the same linkcell build).
- The hexagonal-cage sweep in `seams cages` runs across the threads. One thread and several threads mark the same rings.
- The hydrate signature search scales with threads. The cutoff list is the threaded linkcell list, built on one reused worker pool, and crystals under a few thousand molecules stay at four threads so a wider team does not slow them down. Neighbour rows are placed on the calling thread until the pair list is large enough that a parallel sort is the smaller cost. The cage index is merged from per-thread buckets in ring order. Census and cage membership match the one-thread result.
- The ice cage search counts distinct atoms in a small buffer and builds only the least ring reading.
- The in-plane RDF builds an explicit grid whose cells are the cutoff and bins packed coordinates from the neighbouring cells. Each thread keeps its own histogram until the frame finishes.
- The linkcell wrap pin is 3700ea3 (0.3.5).
- The linkcell wrap tracks v0.3.9 and the minimage wrap tracks v0.1.4. The linkcell patch builds that crate against the minimage wrap checked out beside it.
- The minimage wrap pin is v0.1.2. Distances use mi_displacement_euclidean, and a failed call keeps the in-tree formula.
- The one-thread carried ring update copies and sorts its rings directly again. Several threads still sort indices and merge in parallel. A team narrower than the process is bound once to that many cores of the allowed set. Interior cells of the orthorhombic cutoff skip the periodic-shift arithmetic. One-thread and several-thread results still match.
- The one-thread ring search packs sorted neighbour rows in one pass, and interior cells of a triclinic cutoff skip the image-wrap branches. Both stencils keep the same directed pairs as the branchy walk.
- The orthorhombic cutoff neighbour list is built by a parallel cell sort.
  Each thread counts its atoms, a prefix sum places them, and the pairs are
  written straight into one CSR. The same buffers are reused on the next frame.
  The pair set does not depend on the thread count.
  q_l reads that CSR. The neighbour average is skipped when only q_l is requested.
- The ring search builds each source's shortest paths once, in breadth-first order into one buffer, keeps the graph in two packed arrays, and moves the per-source results out on every thread. With the fix above, `bench_strong`'s ring stage at 65536 atoms goes from 1408 to 265 ms on one thread and from 190 to 40 ms on eight (release builds), and lists the same rings in the same order.
- The seeded cage binding searches one graph when the mutual and union edge sets match. A dirty cage frame starts the backoff before it returns, so the following frames are warm sweeps instead of another diff.
- The seeded cage call searches once when the mutual and union edge sets match. fellBack() runs before a dirty frame returns, so the following frames are warm sweeps.
- The triclinic cutoff list uses the same stencil as the orthorhombic list.
  CHILL+ evaluates the real spherical harmonics on that CSR.
  The diamond common-neighbour test reads the same rows.
  The ring census uses the short-hop search on that CSR.
- Users who run chill-plus on a cutoff list now see different counts.
  This is a bug fix in the shell count.
  The old code counted bonds after the list had been cut to four.
  A crystal row needs four neighbours inside the cutoff.
  On the mixed TIP4P fixture, interClathrate moves from 12 to 9.
  Water on that fixture moves from 238 to 241.
  Default q3 uses up to four nearest neighbours in the cutoff.
  Interfacial rows keep the neighbour test from the published table.
  Users who want the OVITO labels use the ovito-chill+ command.
  That command averages q3 on the full shell and skips the test.
  Plain CHILL keeps four closest neighbours and Moore's test.
  A reader error now returns exit code 2 instead of aborting the process.
  We print reCubic and reHex as plain text.
  A path that ends in xyz and cannot be opened returns an error.
  An exception outside std::exception returns exit code 2.
- `cageAffiliation` runs its hexagonal- and double-diamond-cage sweeps on every thread, and `basalConditions` and `commonElementsInThreeRings` stop allocating per call. Seeded cage affiliation of a 32768-atom jittered cubic-ice frame goes from 1046 to 849 ms on one thread and to 132 ms on eight, and `cages` over 32 such frames from 36.2 to 6.3 s on eight threads (release builds).
- `nearestUnlike` and the threaded neighbour rows score candidates with minimage's `mi_dist2_many` over positions packed once per call, and take the Euclidean image on that call's cell; the in-plane RDF's direct loop past the one-image ball takes each row's Euclidean images in one `mi_dist2_euclidean_many` call. On eight threads `nearestUnlike` over 8192 atoms goes from 24.6 to 3.7 ms, the orthorhombic direct loop from 4.2 to 1.6 ms, and tilted threaded rows from 3.7 to 2.4 ms. Under tilt the row batch takes the direct loop on 2048 atoms at a 30 Å cutoff from 60.4 to 19.5 ms on one thread, against `mi_dist2_many` plus one Euclidean call per pair past the Smith ball, with the same distances bit for bit on default compiler flags.
- `seams chill-plus` on a published cutoff classifies from the fused CSR kernel. `seams steinhardt` computes the degrees named by `--l` (default `4,6`). `seams symmetry` builds the reported shell from the cells inside 5 Å. Primitive rings and the orthorhombic cutoff list stay on one thread below 2048 atoms. The vendored vesin meson project reports 0.5.7.
- `steinhardtQl` builds its bond table with the box computed once rather than once per bond, and counts and fills its rows on every thread; the values are unchanged. On `bench_strong`'s 65536-atom frame the stage goes from 11.7 to 4.0 ms on eight threads and from 23.1 to 19.7 ms on one (release builds).
- findDDC on the mixed-site frame returns 0, not 605243. A double-diamond cage is kept only when its seven rings span 14 molecules and each peripheral shares three equatorial atoms, and every hexagon of that frame fails the test. Clathrate faces are rings whose edges all lie between 2.65 and 3.15 angstroms, so the closure search stays inside one cage instead of the whole ring complex.
- q4 and q6 share one cutoff list. A multi-thread cutoff walk records each pair once. Seeded cages with identical neighbour rows search rings and cages once, and the k-nearest nomination uses the same thread count as the rest of the command.

### Fixed

- A LAMMPS dump that grows while it is read is kept. The language fronts follow the build Python and skip a library mismatch.
- A tilted cell's cage edge window uses the lattice minimum image, so a GenIce sH dump keeps every 5^12, 4^3 5^6 6^3 and 5^12 6^8 cage. The head-to-head against ClavisCage and cycless is benchmarks/cage-h2h-v8.txt.
- Accept a cup that fills the signature face count with dangling edges.
- An orthorhombic displacement wraps every period, not one. Unfolded coordinates 16 apart in a 10 box are 4 apart, not 6.
- Build the Sphinx book with warnings as errors from the committed sources. Doxygen XML and the C++ pages are optional, and each page is listed in one toctree.
- C API: `seams_ql`, `seams_q4q6`, `seams_q4q6_n`, `seams_q6bar`, and `seams_steinhardt` kept their q_lm scratch in function-local statics, so two of them on different handles in different threads corrupted each other's results or the heap. The scratch now belongs to the calling thread, and the one-time affinity setup they reach is guarded, so q calls on different handles may run at the same time.
- C API: calls on a const handle no longer write the handle on success, and take a lock to record a failure, so several threads may query one handle at once. The ring label cache that `seams_rings_labels` fills on the first symmetry request is now built once under a lock and marked done only when complete; before, a second thread arriving during that request, or any request after a failed one, could read unfinished labels.
- C API: libseams exports only the `seams_*` functions in `seams.h` (hidden visibility plus a linker export list). No call aborts or exits: the MIS self-test throws instead of calling `exit`, every allocation in the wrappers is caught, a missing trajectory returns `SEAMS_ERR_ARG`, and infinite cutoffs are rejected by `seams_hbonds` and `seams_monolayer`.
- Cage object counts no longer follow the wrapped coordinates. `findDDC` reports a maximum set of non-conflicting double-diamond cages, and both searches order rings by sorted atom indices plus the unwrapped shift of the cycle. A periodic translation or a permutation of the atoms leaves the count unchanged.
- Cage search drops rings whose minimum-image lattice shifts wind through the periodic cell, and refuses a cutoff that can see two images of one pair. The R, Julia and C fronts keep the seeded label on the host rule when offload is requested. `seams cages` and the language fronts report disjoint hexagonal, cubic, both and water counts. `scripts/rerun_main_text.py` reruns the engine counts behind the main-text cage numbers and writes the archive hash.
- Clang builds record the directory of clang's `libomp.so` as a runpath when it lies outside `/usr`, so the libraries, tests and libFuzzer harnesses start without `LD_LIBRARY_PATH`.
- Coordination and homopolar counts skip the leading self entry on an index neighbour row.
- Evaluate Steinhardt l=12 qlBar on the host so the 25-component average does not overflow the device barRe[17] buffer.
- Host OpenMP cell-list loops no longer fork under nvc++ `-mp=gpu`, so `[bulkTUM][offload]` reaches `usedDevice`.
- In OpenMP builds the threaded neighbour rows gave an atom missing from `idIndexMap` its neighbours' IDs as a row header, and gave its partners -1. That atom now keeps an empty row.
- Keep spherical-harmonic scratch local to each call so concurrent frames and changing thread counts do not share buffers.
- Keep the 2.13 signatures of `getCorrelPlus`, `getIceTypePlus`, and `getIceTypePlusNoPrint` as overloads. The longer forms take the open-edge flag and the OVITO neighbour test.
- Load oxygen and hydrogen together for `seams f4` so Rodger F4 is finite when mol IDs exist.
- On a crowded shell the diamond label uses the four nearest neighbours, the list mdapy uses, and the labels stay those of the cutoff list.
- Place sI hydrogens with TIP3P HOH and a 15 degree libration so Rodger F4 sits near 0.7, score Zeron q12bar against a disordered liquid, and evaluate host l=12 Ylm when sphericart is off.
- Rename the site-site pair flag from `--types` to `--pair-types`, so `--type` is not a prefix of another flag. Give the gpulite and vesin wraps their own option files so `meson configure` does not copy the engine options onto them.
- Require readcon-db 0.1.7 so static consumers receive only the archive and shared consumers receive the runtime search path.
- SOAP now fills every l through lMax. classifyBonds reuses a Ylm buffer. Fingerprint builds one hop graph per atom.
- The 2D RDF histogram, `nearestUnlike`, and `shellSeparation` trust the fractional wrap only below half the narrowest face separation of a tilted cell. Half an edge or half a bound span let a pair inside the cutoff wrap to a longer image and drop out, and let the RDF grid put a pair two cells apart.
- The C API builds q4 and q6 from one cutoff shell and does not fill an atom-index map on that path.
- The PLUMED q6 derivative is expressed per nanometre. A 0.001 angstrom step matched an independent difference to 3.99e-11 on a shifted eight-molecule cell.
- The cutoff neighbour list sizes its row buffer from the pairs found. A shell denser than the ideal-gas estimate no longer writes past that buffer.
- The cutoff neighbour search clamps its cell count and stencil reach before converting them to integers. A cell edge many orders of magnitude shorter or longer than the cutoff overflowed both, which is undefined behaviour and in practice left the reach at one cell. Found by `capi_reload_fuzzer` under UBSan.
- The library builds with clang again: the prism code in `topo_one_dim.cpp` casts its signed ring and atom indices explicitly, as the rest of the tree does, so the `-Werror=sign-conversion` and `-Werror=shorten-64-to-32` that clang builds use no longer stop it.
- The multi-thread cutoff list records neighbour ids during the count walk, and the linkcell pool is sized to the OpenMP width.
- The nucleation notebook rule declares ``figshare-incremental.json``
  as an output. A later ``test -s`` stamp cannot see a side-effect
  file because Snakemake deletes declared outputs before the job.
- The reproducibility DAG builds yodaStruct as ``require("dseams")`` and
  runs the five figshare deposits through that library. There is no
  in-tree ``yodaStruct`` binary.
- The ring search no longer clears the level field over every lower index for each source, a pass quadratic in the frame: on one thread `bench_strong`'s ring stage grew about threefold with each doubling of the frame, to 16.9 s at 262144 atoms, and now doubles.
- Use linkcell 0.3.10 and minimage 0.1.4 without a downstream source patch. The with_linkcell and with_minimage options can require these dependencies.
- Use the one-cell neighbour stencil from two cells per edge, so a short box is not left on the general walk.
- `cages --signature` counts a face set as a cage only when its faces tile a sphere: every vertex lies on at least three faces that close around it, and V - E + F = 2. Rings whose minimum-image bonds add up to a lattice vector wind through the periodic cell, bound no face, and are dropped first (`ring::contractibleRings`). On frame 1 of the figshare mW nucleation run (4096 molecules, 3.5 Å cutoff) the search found 3 5¹²6², 20 5¹²6⁴, 270 5¹²6⁸ and 143 4³5⁶6³ cages in the liquid; it now finds none. On the GenIce sI unit cell it found 90 5¹²6², 4248 5¹²6⁴ and 218305 5¹²6⁸; it now finds 2 5¹² and 6 5¹²6².
- `seams_rings_find` classifies prisms from the neighbour index instead of
  an all-pairs scan, and the ice Ih stacking-plane pass runs across
  threads. Ring lists and labels on the 4096- and 32768-molecule cells
  are unchanged.
- `tests/check_invariance` exits 2 when the trajectory cannot be read, the frame holds no atom of the requested type, or the graph has no six-ring. It used to compare empty results, print PASS and exit 0. The meson tests `invariance_missing_input` and `invariance_empty_frame` expect that failure.
- q4 and q6 build the multi-thread cutoff list in one pass again, instead of copying each linkcell block.

## [2.13.0] - 2026-11-06

### Changed

- The bond test used a dense bit for every atom pair.
  Each row is now the sorted neighbour ids.
- The cage search hashes the sorted undirected edges and compares the lists.
  One search is reused only when the edge sets match.
  Each graph is searched on its own when the lists differ.
  The default cages command uses that check.
- The paper records the release comparison on one core.
  The build rule writes the ring times and the neighbour times.
  It writes the two commands, the resident set, and the abort.
  The Open Visualization Tool call includes startup and draws no cages.
  The diamond finder gives one label to each point.
- The ring search reads one offset and column table.
  A single thread writes one output list.
  The scan of those rows is the parallel prefix in par_team.hpp.
- The rings of one atom are one block of a single list.
  The search stops when a candidate shares an atom with the ring.
  A clash is stored in a row that grows as further clashes are found.
- v2.12.0 was tagged on a tree whose CI was red.

### Fixed

- The Ozaki tests link libmpfr and libgmp.
  Linux pixi environments did not install either library.
  The pixi environment installs both, so the test build file finds them.
- The bindings job fetches Julia 1.11.6 from the published x64 archive.
- The conda libc++ for macOS 11 deletes floating std::from_chars.
  Cage edge windows and runtime double options call that function.
  A C locale conversion accepts the same decimal text on every toolchain.
- The coverage build updates arc counters, and 8 cores bounce them.
  The bulk and cage tests then exceed the existing time limit.
  One thread and a non-atomic update keep those counters and those tests.
  The cage fuzz pool follows that same one-thread request.
- The ignore pattern run* matched runtests.jl in every directory.
  The pattern is anchored at the repository root, and the bindings job runs that file.
- The runner ignores sudo -E, so a root install cannot see the engine paths.
  cpp11, testthat, and dseams install into the user library, and the configure script sees those paths.
  testthat's fs package builds against libuv headers.
- The version pin howto names the version written in the other version sources.
- pydseamslib chill_plus takes no graph argument.
  That call labels the cutoff neighbour list.
  The agreement test builds the mutual list itself and checks the cutoff counts.
  No graph argument is added to the Python package.


## [2.12.0] - 2026-11-01

### Added

- R (`bindings/r/dseams`) and Julia (`bindings/julia/Dseams`) call `libseams_lang`.
  Python stays `pydseamslib`.
  The executable tutorials under `docs/tutorials` import `pydseams`.

### Changed

- On the all-site niu2f1 frame, one thread, one run of `findHC`, `findDDC` and `findBySignatures` together takes 176.4 s.
  TRACE on that frame takes 1762.064 s.

Eleven release runs of `bench_cage_h2h` on this tree.
The median is the 6th value.
The bracket is the 3rd and the 9th.

| input | threads | d-SEAMS ms |
|---|---|---|
| mw | 1 | 10.2957 [10.2441, 10.3796] |
| mw | 4 | 5.1094 [5.0790, 5.1514] |
| mw | 8 | 4.4331 [4.1860, 4.9120] |
| nuc148 | 1 | 7.9022 [7.8521, 7.9838] |
| nuc148 | 4 | 5.9419 [5.8551, 6.0224] |
| nuc148 | 8 | 5.4825 [5.2546, 8.8501] |
| sI | 1 | 0.2567 [0.2517, 0.2630] |
| sI | 4 | 0.2555 [0.2534, 0.2673] |
| sI | 8 | 0.2586 [0.2521, 0.2645] |
| sII | 1 | 0.5342 [0.5293, 0.5411] |
| sII | 4 | 0.4253 [0.4163, 0.4297] |
| sII | 8 | 0.3941 [0.3901, 0.4181] |
| niu3f148 | 1 | 5.1213 [5.0898, 5.1600] |
| niu3f148 | 4 | 4.6676 [4.6178, 4.7560] |
| niu3f148 | 8 | 4.4719 [4.4143, 4.5055] |
| sI222 | 1 | 1.3870 [1.3747, 1.3927] |
| sI222 | 4 | 0.9477 [0.9215, 0.9606] |
| sI222 | 8 | 0.8243 [0.8139, 0.8585] |
- findDDC on the mixed-site frame returns 0, not 605243.
  A double-diamond cage is kept only when its seven rings span 14 molecules and each peripheral shares three equatorial atoms, and every hexagon of that frame fails the test.
  Clathrate faces are rings whose edges all lie between 2.65 and 3.15 angstroms, so the closure search stays inside one cage instead of the whole ring complex.
- The clathrate face-edge window is `SEAMS_CAGE_EDGE` and `--cage-edge`, default `2.65,3.15` in the coordinate units of the file.
  `none` turns the window off, and the cage bound is then the face census alone.
  `seams_frame_signature_edge` takes the same string.
  Python is `frame.cages_by_signature("512", edge="none")` and Lua is `dseams.cages_by_signature(frame, "512", {edge = "none"})`.
- A peripheral set is kept only when the seven rings span 14 molecules and each peripheral shares exactly three atoms with the equator.
  On a 5508-frame census the seeded partition changes on 498 frames, all from this rule.
  The deposited main-text seeded counts do not move.
- CHILL and CHILL+ use one default neighbour list on every front: the cutoff list (pairs closer than `cutoff`, 3.5 Å by default).
  The mutual four-nearest list is `graph = "mutual"` under that name everywhere: `seams chill` and `seams chill-plus --graph mutual`, R `chill_plus(graph = "mutual")`, Julia `chill_plus(; graph = "mutual")`, pydseamslib `Frame.chill_plus(graph="mutual")`, the C ABI `seams_chill_plus_graph`, and the LAMMPS `compute dseams ... graph mutual`.
  `seams chill`, R, Julia, `seams_chill_plus`, and `compute dseams` used the mutual list by default before.
  `tests/chill_plus_fronts.py` checks that the CLI, Python, R, and Julia agree on `exampleTraj.lammpstrj` with both lists.
- Release builds define NDEBUG, so the libstdc++ assertions that meson adds when b_ndebug is false no longer run in them.
  The closed-cage search keeps the same output.
- The closed-cage search skips a ring that shares two consecutive edges with a chosen face or passes through a vertex whose faces already close, and keeps the vertex and edge counts inside V - E + F = 2.
  No cage is lost.
  The search runs on dense edge ids, adds the forced faces and picks the open edge with the fewest rings in one pass, and lists the cages of a signature by their sorted faces.
- The orthorhombic cutoff list folds coordinates into the cell, bins them one cutoff wide, and shares the 27-cell stencil across the atoms in a cell.
  Interior cells skip the periodic wrap.
  A squared distance within 1e-8 of the cutoff squared falls through to the existing list.
  The walk is used when the box is orthogonal, the subset has at least 2048 atoms, and every edge is at least three cutoffs.
  CHILL and CHILL+ on the default left-to-right fp64 path evaluate the l=3 harmonics as real and imaginary parts under OpenMP.
  On all 5508 frames the neighbour rows match the pre-patch rows and the CHILL+ ice labels match.
- When the product of rings through the even or odd triplets exceeds 256, `equatorialConditions` builds the chains from edges that share at least three atoms and keeps the lexicographically least assignment that spans 14 molecules and gives each peripheral three equatorial atoms.
  Smaller products stay on the triple loop.
  Cubic mW stays at 2048 double-diamond cages, certificate method 1.
  Nucleation frame 148 stays at 102 hexagonal cages, 121 double-diamond cages, and 223 cages in total, with seeded counts 255, 363, and 377.
- `findDDC` keeps a maximum non-conflicting set and certifies it.
  The mutual-conflict graph is built once.
  A one-cage component is kept.
  A component whose neighbourhoods are the two cliques of a bipartite line graph is a maximum matching, certified by a König cover.
  Any other component is a branch-and-bound maximum independent set, certified by a clique partition.
  On cubic mW the certificate is method 1 and the double-diamond count is 2048; the unpatched library reported 2040.
  On nucleation frame 148 the double-diamond count is 121 and the hexagonal count stays 102; the unpatched library reported 122 double-diamond cages.
  Seeded labels on that frame stay 255 hexagonal, 363 cubic, and 377 both, and the deposited mW seeded row stays 0 / 4096 / 0 / 0.

### Fixed

- The equatorial search keeps a degree-one candidate and drops only that one conflict.
  K4 minus an edge keeps two vertices.
- Cage object counts no longer follow the wrapped coordinates.
  `findDDC` reports a maximum set of non-conflicting double-diamond cages, and both searches order rings by sorted atom indices plus the unwrapped shift of the cycle.
  A periodic translation or a permutation of the atoms leaves the count unchanged.
- Cage search drops rings whose minimum-image lattice shifts wind through the periodic cell, and refuses a cutoff that can see two images of one pair.
  The R, Julia and C fronts keep the seeded label on the host rule when offload is requested.
  `seams cages` and the language fronts report disjoint hexagonal, cubic, both and water counts.
  `scripts/rerun_main_text.py` reruns the engine counts behind the main-text cage numbers and writes the archive hash.
- `cages --signature` counts a face set as a cage only when its faces tile a sphere: every vertex lies on at least three faces that close around it, and V - E + F = 2.
  Rings whose minimum-image bonds add up to a lattice vector wind through the periodic cell, bound no face, and are dropped first (`ring::contractibleRings`).
  On frame 1 of the figshare mW nucleation run (4096 molecules, 3.5 Å cutoff) the search found 3 5¹²6², 20 5¹²6⁴, 270 5¹²6⁸ and 143 4³5⁶6³ cages in the liquid; it now finds none.
  On the GenIce sI unit cell it found 90 5¹²6², 4248 5¹²6⁴ and 218305 5¹²6⁸; it now finds 2 5¹² and 6 5¹²6².
- `tests/check_invariance` exits 2 when the trajectory cannot be read, the frame holds no atom of the requested type, or the graph has no six-ring.
  It used to compare empty results, print PASS and exit 0.
  The meson tests `invariance_missing_input` and `invariance_empty_frame` expect that failure.
- The Sphinx book builds with no warnings.
  The CHILL tables are verbatim, the book title underline matches the title, the nanotube page is in the tutorial contents, and each generated API link has a target.


## [2.11.0] - 2026-10-09

### Added

- Ice XXI library matches the Lee et al. 2026 I-42d oxygen cell (SI CIF, Z=152).
- Under mpiexec, an MPI build of `seams` shares a LAMMPS dump's `--frame`/`--last` frames among the ranks in serpentine rounds (`sinp::forEachLammpsFrame` takes a part and a part count), and rank 0 prints every frame in order. Other readers stay on rank 0, `--per-atom` needs a single rank, ranks that differ in arguments, directory or `SEAMS_` settings are refused, `OMP_NUM_THREADS` and `RAYON_NUM_THREADS` default to 1 on a launch of several ranks, and a run without a launcher, or with `SEAMS_MPI=0`, never starts MPI. `chill::setSteinhardtAtomSplit(false)` keeps `steinhardtQl` from splitting atoms over ranks that hold different frames, and a dump frame cut short keeps the atoms it read, with its warning on stderr. With Open MPI's libfabric probe skipped (`--mca btl ^ofi`), eight ranks match eight `--jobs` workers on one node for `cages` over 32 frames of 32768-atom ice (4.2 against 4.3 s), and keep load balance 0.99 where frame cost steps or alternates along the run (2.6 s against 2.8 and 3.1 s for `--jobs`; contiguous blocks reach 0.59 on the step). On two four-core core sets standing in for nodes, with MPI over TCP only, eight ranks take 4.4 s against 8.7 s for threads on one, at POP parallel efficiency 0.98 (release builds).
- `seams::domain` splits a periodic frame over ranks: cells ordered along a Hilbert curve are cut into runs of equal atom count, and each rank also holds every atom within a halo of its own. `domain::rings` finds the primitive rings whose lowest-indexed member a rank owns from its share alone, `primitive::ringNetwork` takes a mask of the sources to enumerate from, and `domain::gatherRings` collects every rank's rings on every MPI rank in the order `ringNetwork` lists them. On one eight-core node threads still win (a 65536-atom frame takes 58 ms on one rank of eight threads and 67 ms on eight single-thread ranks), and every rank rebuilding the whole list costs 72 ms of 251 at 262144 atoms on eight ranks; `bench_strong` reports it as `domain/ms` in MPI builds, and CI runs the gather under `mpiexec`.

### Changed

- Below half the narrowest face separation the cutoff list skips minimage's pair reduction, which only sorted there, and sorts its rows on every thread; `getNewNeighbourListByIndex` and `neighListO` fill their rows in parallel, and `neighListO` takes the cutoff list first at every thread count. On `bench_strong`'s 65536-atom frame the index list goes from 25.8 to 13.9 ms on one thread and from 22.7 to 7.8 ms on eight (release builds).
- CON frames are read with readcon-core's frame iterator. Lengths and angles become the minimage cell, and that cell is the dump box linkcell already uses. A readcon-db corpus, including a sharded campaign root, is selected and decoded with the same reader.
- Cutoff neighbour lists use linkcell pairs_within with cells of about 20 atoms rather than one cutoff, ahead of the threaded cell rows. The neighbour list still keeps one index per atom. k-nearest stays knearest.
- Inside the one-image ball the in-plane RDF walks Rapaport cell pairs: cells span the cutoff across each face separation, a half stencil visits each neighbouring pair of cells once, and one lattice shift serves the pair. It replaces the packed grid and the vesin list; 8192 atoms at a 6 Å cutoff take 2.2 ms on one thread and 0.4 ms on eight.
- Pair distances use the fractional wrap below half the narrowest face separation, and call the Euclidean image only past that. The in-plane RDF bins the cell-list distance inside that ball, and the direct loop runs in parallel.
- The SOAP power spectrum evaluates each neighbour's spherical harmonics once and reuses them for every radial function, since its basis is a product of the two; the values are unchanged. On 8192 atoms at water density with nMax 8, lMax 6 and a 6 Å cutoff, `soapSpectrumAll` goes from 831 to 199 ms on one thread (release builds).
- The cutoff lists bin linkcell's pairs by row on every thread: the pairs are dealt into blocks of rows, four per thread, and one thread counts, places and sorts each block, so no pass waits on shared counters and the rows do not depend on the thread count. linkcell fills a buffer through `lc_pairs_within_rows`, where `linkcell::pairs_within` zeroes its vector on one thread first, and `neighListO` allocates each row once, on the thread that fills it. On a jittered 65536-atom lattice at water density the cutoff list goes from 17.6 to 14.5 ms on one thread, and at 262144 atoms from 67.7 to 57.6 ms on one thread and from 35.1 to 15.2 ms on eight (medians of three interleaved rounds, release builds, both on the same linkcell build).
- The linkcell wrap tracks v0.3.9 and the minimage wrap tracks v0.1.4. The linkcell patch builds that crate against the minimage wrap checked out beside it.
- The ring search builds each source's shortest paths once, in breadth-first order into one buffer, keeps the graph in two packed arrays, and moves the per-source results out on every thread. With the fix above, `bench_strong`'s ring stage at 65536 atoms goes from 1408 to 265 ms on one thread and from 190 to 40 ms on eight (release builds), and lists the same rings in the same order.
- `cageAffiliation` runs its hexagonal- and double-diamond-cage sweeps on every thread, and `basalConditions` and `commonElementsInThreeRings` stop allocating per call. Seeded cage affiliation of a 32768-atom jittered cubic-ice frame goes from 1046 to 849 ms on one thread and to 132 ms on eight, and `cages` over 32 such frames from 36.2 to 6.3 s on eight threads (release builds).
- `nearestUnlike` and the threaded neighbour rows score candidates with minimage's `mi_dist2_many` over positions packed once per call, and take the Euclidean image on that call's cell; the in-plane RDF's direct loop past the one-image ball takes each row's Euclidean images in one `mi_dist2_euclidean_many` call. On eight threads `nearestUnlike` over 8192 atoms goes from 24.6 to 3.7 ms, the orthorhombic direct loop from 4.2 to 1.6 ms, and tilted threaded rows from 3.7 to 2.4 ms. Under tilt the row batch takes the direct loop on 2048 atoms at a 30 Å cutoff from 60.4 to 19.5 ms on one thread, against `mi_dist2_many` plus one Euclidean call per pair past the Smith ball, with the same distances bit for bit on default compiler flags.
- `steinhardtQl` builds its bond table with the box computed once rather than once per bond, and counts and fills its rows on every thread; the values are unchanged. On `bench_strong`'s 65536-atom frame the stage goes from 11.7 to 4.0 ms on eight threads and from 23.1 to 19.7 ms on one (release builds).

### Fixed

- An orthorhombic displacement wraps every period, not one. Unfolded coordinates 16 apart in a 10 box are 4 apart, not 6.
- Coordination and homopolar counts skip the leading self entry on an index neighbour row.
- Host OpenMP cell-list loops no longer fork under nvc++ `-mp=gpu`, so `[bulkTUM][offload]` reaches `usedDevice`.
- In OpenMP builds the threaded neighbour rows gave an atom missing from `idIndexMap` its neighbours' IDs as a row header, and gave its partners -1. That atom now keeps an empty row.
- SOAP now fills every l through lMax. classifyBonds reuses a Ylm buffer. Fingerprint builds one hop graph per atom.
- The 2D RDF histogram, `nearestUnlike`, and `shellSeparation` trust the fractional wrap only below half the narrowest face separation of a tilted cell. Half an edge or half a bound span let a pair inside the cutoff wrap to a longer image and drop out, and let the RDF grid put a pair two cells apart.
- The cutoff neighbour list sizes its row buffer from the pairs found. A shell denser than the ideal-gas estimate no longer writes past that buffer.
- The ring search no longer clears the level field over every lower index for each source, a pass quadratic in the frame: on one thread `bench_strong`'s ring stage grew about threefold with each doubling of the frame, to 16.9 s at 262144 atoms, and now doubles.
- Use linkcell 0.3.10 and minimage 0.1.4 without a downstream source patch. The with_linkcell and with_minimage options can require these dependencies.


## [2.10.0] - 2026-09-06

### Added

- Add Rodger F4 and host-only Steinhardt l=12 so the Zeron q3/q12 hydrate pair is a per-atom field.
- Add TUM/rings stacking planes (HC-basal and DDC-equatorial) next to the CHILL+ I_sd molecule-bin reference.
- Add a C ABI seams_chill_plus for native and emscripten builds.
- Add extras/lammps compute dseams that writes seams_chill_plus per atom.
- Assign guests by ray-parity inside fan-triangulated cage faces, and count ions per ice cluster on the oxygen graph.
- Bin CHILL+ cubic and hexagonal molecules into basal layers and emit cubicity Phi_c plus the H/C stacking string.
- Build CHILL+, cages, seeded ions and k-NN on an explicit water-type mask (`--water-types`) so substrate and ions never enter the four-neighbour list. The flag loads the mixed dump; ice and cage counts are the water types only.
- Enumerate cups and incomplete cages on the ring graph; closed signature counts stay the closed path.
- Hexagonal channels are stacked six-ring prisms, not a raw six-ring count. Ice XXI (Lee et al. 2026, Z=152 BCT) also requires a tetrahedral 3.5 A graph and a primitive six-ring. Hydrogen MSD uses the minimum image. Glass labels use local-density windows for ice/LDA/MDA/HDA. `compute dseams` ships with `water.data` and `pair_style zero`.
- Name 51264, 51268 and sH cages, and report a per-cage occupancy histogram.
- OpenMP target offload of the TUM ice score: hop-bound primitive six-rings and HC/DDC cage affiliation. `SEAMS_OFFLOAD=1` cage counts match the host on mW cubic, including `seams cages --graph seeded` (union 4-NN graph). Device CHILL+ is not this path.
- Update topology keys only on the hop-ball of atoms whose neighbourhood changed.

### Changed

- CHILL and CHILL+ default to the mutual four-nearest graph. vesin still owns cutoff pairs; linkcell still owns k-nearest.
- Dump MIC and vesin pair reduction go through the minimage wrap. The simd Catch2 binary links it when the wrap is present.

### Fixed

- Accept a cup that fills the signature face count with dangling edges.
- Evaluate Steinhardt l=12 qlBar on the host so the 25-component average does not overflow the device barRe[17] buffer.
- Load oxygen and hydrogen together for `seams f4` so Rodger F4 is finite when mol IDs exist.
- Place sI hydrogens with TIP3P HOH and a 15 degree libration so Rodger F4 sits near 0.7, score Zeron q12bar against a disordered liquid, and evaluate host l=12 Ylm when sphericart is off.


## [2.9.2] - 2026-09-02

### Added

- The cutoff neighbour list and the index-ordered list are built by a threaded cell list (`nneigh::cellListRowsThreaded`) from 2048 atoms when OpenMP is compiled: atoms are binned in fractional coordinates of the recovered triclinic cell, every axis needs three cells of perpendicular width one cutoff, and each row is one thread's task; the rows equal the minimum-image reference. `neighbourListByIndex` runs one row per thread. `bench_strong` times every host stage of the ring pipeline, prints their sum and takes the best of five runs, so the strong-scaling figure covers the whole pipeline rather than the ring stage alone.

### Changed

- Document the two LAMMPS type-filtered readers: ~readLammpsTrjO~
  flags ~inSlice~ and keeps every atom of the type;
  ~readLammpsTrjreduced~ drops atoms outside the slice. An axis with
  ~lo == hi~ is unconstrained.

### Fixed

- The nucleation notebook rule declares ``figshare-incremental.json``
  as an output. A later ``test -s`` stamp cannot see a side-effect
  file because Snakemake deletes declared outputs before the job.
- The reproducibility DAG builds yodaStruct as ``require("dseams")`` and
  runs the five figshare deposits through that library. There is no
  in-tree ``yodaStruct`` binary.


## [2.9.1] - 2026-09-02

### Added

- OpenMP target offload for the Steinhardt kernel compiles and links with nvc++ 23.7 (`-mp=gpu`) against CUDA 12.2 on an NVIDIA A100-PCIE-40GB. The Catch2 suite is green, including the identity check that the device path matches the serial and threaded host paths bit for bit on the FCC lattice and on `input/traj/mW_cubic.lammpstrj`. `nsys stats --force-export=true` writes a non-empty kernel and memcpy summary. `SEAMS_OFFLOAD=0` still forces the host path.

## [2.9.0] - 2026-09-02

### Added

- `cage::Signature` names a polyhedron by its ring-size census (`sodalite` is `{4:6, 6:8}`). `cage::findBySignature` grows face-sharing primitive rings that close (every edge shared by exactly two faces). `seams cages --signature 4:6,6:8` or a named table entry (`sodalite|alpha|512|51262|hc|ddc`) prints the cage and atom counts. Named `hc` and `ddc` call `findHC` / `findDDC` so the vertex sets match those finders; a raw list `4:6,6:2` is the hexagonal prism.
- `topo::matchLibraries` names atoms against key libraries at several hop counts, the deepest library that holds an atom's key winning, and reports the depth that named each atom, so a molecule whose wide neighbourhood is disturbed still gets its name from the inner shells. Libraries record whether their keys carry vertex colours (`colours` in the header; an absent field reads as uncoloured) and a coloured library is refused against plain keys. `seams fingerprint --library` accepts a comma separated list and prints the per-depth counts. `site::guestOccupancy` places guests (methane, THF, ions) in enumerated cages by the periodic centroid of each cage's vertices and counts occupied, multiply occupied and free, the occupancy of a clathrate hydrate; `site::periodicCentroid` is the helper. `site::IonEnvironment` lists each ion's shell molecules and `site::shellRingCensus` counts the rings of the water network that pass through a shell, by size: how far the network survives around an ion. `seams cages --signature SPEC --guest-types T,U` places guests in the found cages and reports occupied, multiply occupied and free. `--per-atom FILE` on `cages`, `fingerprint` and `ions` appends a LAMMPS dump frame with one extra column (cage membership, topology class or library label, water and ion state) so OVITO or VMD colour the trajectory by the engine's decision.

## [2.8.0] - 2026-09-02

### Added

- `topo::localKey` names the isomorphism class of an atom's rooted bonded neighbourhood within a number of hops (the nauty certificate with the centre in its own colour cell when nauty is linked, a Weisfeiler-Lehman refinement hash otherwise); `topo::fingerprint` keys a frame by the sorted refinement hashes and the primitive ring census, so relabelled configurations share a key. `site::ionEnvironment` classes ions by their first water shell against a per-atom ice flag. `seams fingerprint` and `seams ions` expose both; `seams cages --graph seeded --complete` turns on the ring completion; `--hops`, `--ion-types` and `--ion-cutoff` are new options. Vertex colours (atom types) partition the keys by species (`--colour-types`). A `topo::KeyLibrary` collects the keys of reference structures under labels and names any atom whose key it holds; a key shared by several references carries all their labels. `seams fingerprint --emit-library LABEL` writes a frame's keys as library lines and `--library FILE` names the atoms of a frame by one. The readcon-core fallback follows the engine's library kind, so a static engine embeds the CON reader; CON and chemfiles input without their reader is refused instead of returning an empty frame.

## [2.7.0] - 2026-09-02

### Added

- `ring::seededCageAffiliation` takes an optional completion flag: `ring::ringAdjacentCompletion` fills the last vertex of any six-ring whose other vertices carry a cage label, iterated to a fixed point, separately for the HC and DDC labels. The rule is stated and proven in `lean/` (Mathlib): the completion is the least fixed point above the seed, empty on an empty seed, sound, and independent of visiting order; the edge-sharing rule it replaces is shown unsound on a five-vertex instance. `tests/walk_compare` walks a trajectory and prints per-frame CHILL+ and cage labels with their largest clusters. The reproducibility workflow moved to the `dseams2_repro` package.

## [2.6.0] - 2026-08-17

### Changed

- Meson links FlexiBLAS when pkg-config or libflexiblas is present, and falls back to the conda-forge libblas/liblapack ABI packages. The reproducibility campaign no longer configures in-tree Python; bindings come from pydseamslib on PyPI. The nucleation notebook reads CageScore fields. Figshare Lua demos resolve example_lua from yodaStruct.

## [2.5.0] - 2026-08-16

### Changed

- `--family` (default `waterIce`) is an input. Ice scores refuse non-water families and name the family. `seams cn --ions` is cage degree on `site::ionCloud`; `seams pairs` is the mutual-nearest contact-pair count, not ionicity. `seams domains` is the largest polar or apolar Stoddard component. `seams density-z` is type-resolved `rho(z)` with slab volume from dump H. Running CN and the first minimum of `g_IJ` come from `rdf::runningCN` / `firstMinimumBin`. Dump writers emit `xy, xz, yz` when the cloud carries tilt.

## [2.4.0] - 2026-08-16

### Changed

- Partial 3D site-site `g_IJ(r)` and coordination numbers under the dump MIC, with `seams rdf` and `seams cn`. Site chemistry is a `site::Table` of LAMMPS types and optional atom-ID overrides; `site::ionCloud` collapses each ion `molID` to one unwrapped COM vertex. Hydrogen bonds accept an explicit donor-H index set (`populateHbondsFromDonors`). `seams hbonds --donors` uses every hydrogen as a donor candidate; water cutoffs stay on that command and are not an ionic-liquid criterion. In-plane RDF normalization uses the dump-cell volume `|det H|` and the restricted-triclinic in-plane area `lx*ly`, not the particle AABB or bound-span product.

## [2.3.4] - 2026-08-16

### Changed

- Steinhardt, hydrogen-bond, and cluster wraps use the dump minimum image, not an independent-axis wrap on bound spans.

## [2.3.3] - 2026-08-16

### Changed

- Cutoff fallbacks and the skin Verlet refresh use the dump 3x3 minimum image, not a diagonal of bound spans.

## [2.3.2] - 2026-08-16

### Changed

- Cutoff lists take the LAMMPS dump 3x3, not a diagonal of bound spans. Device ice-score batches pass `lammpsBoxToLinkcell` when the dump tilt is present. `linkcell` wrap pin is v0.3.1.

## [2.3.1] - 2026-08-16

### Changed

- `neighList`, `halfNeighList`, and the in-plane RDF sampler use the vesin cell list (same helper as `neighListO`). Brute force remains the fallback. The GPU ice-score workspace is unchanged.

## [2.3.0] - 2026-08-16

### Changed

- The TUM ice score (mutual four-nearest graph, primitive hexagons, HC/DDC) runs as a device-resident frame batch when gpulite and `linkcell` 0.3.0 are present. CHILL and `q_{lm}` stay on the host. Runtime knobs are twelve-factor: `SEAMS_CONFIG` or `./seams.env`, then the environment, then CLI flags. `seams --print-config` dumps the table. Examples and the README use `seams` / `pydseams` / `require("dseams")`; the 2020 `yodaStruct -c` / `conf.yaml` surface is gone.
- The C++ API in the book is Doxygen via doxyrest (`api/index`).

## [2.2.5] - 2026-08-16

### Changed

- `subprojects/linkcell.wrap` is v0.2.4. A static `libyodaLib.a` no longer passes `liblinkcell.so` to `ar`. Wrap consumers link the static archive, so `pydseams.yoda` does not need `liblinkcell` at import time. Dump doubles parse with `strtod` on Apple libc++.

## [2.2.4] - 2026-08-16

### Changed

- Periodic *k*-nearest graphs via `d-SEAMS/linkcell` v0.2.2. `kNearestNeighbourList` / `kNearestNeighbourPair` write packed `n*k` nominations. Mutual and union come from one walk. The LAMMPS dump box is bound spans plus `xy, xz, yz`; `periodicDistSq` recovers H. Empty frames no longer abort the walker.
- Docs orgmode covers the `seams` CLI, the Nix flake, and the three-repo split. Tutorials use the live APIs.

## [2.2.3] - 2026-08-15

### Changed

- The `seams` CLI uses Argum, the same parser as eonclient. Help, errors, ice-type counts, and `--features` are colorized. `NO_COLOR` turns the colors off.

## [2.2.2] - 2026-08-15

### Changed

- The compiled Python surface is `pydseams.yoda`. Docs and the repro scripts use that name. `_core` remains an alias.

## [2.2.1] - 2026-08-15

### Added

- Flake-based Nix package for `libyodaLib` and the `seams` CLI
- (meson, not the CMake-era `yodaStruct` derivation).
- `nLammpsFrames` / `dropLammpsDumpIndex`: live dump session with a
- lazy `ITEM: TIMESTEP` offset table (LAMMPS `ReaderNative` cursor, chemfiles `read_step`, readcon frame offsets). Sequential `load_frame` walks no longer rescan prior snapshots.
- `forEachLammpsFrame`: OpenMP walk over a frame range. Each worker
- opens its own handle and seeks. `seams --frame N --last M --jobs J`.
- `SkinNeighborList`: vesin candidates at cutoff+skin, rebuilt on
- the Verlet trigger. `BondGraph` is chosen at runtime (`cutoff`, `knn`, `knn-union`). `seams cages --graph` adds `seeded`.

### Fixed

- LAMMPS readers bind `xu`/`yu`/`zu` and `xs`/`ys`/`zs`
- when `x y z` are absent (Niu/Parrinello TIP4P/Ice dumps).

## [2.2.0] - 2026-08-15

### Added

- `seams` CLI: `read`, `chill`, `chill-plus`, `cages`.
- This is the engine command line. Lua is the `dseams` library (yodaStruct repo). Python is `pydseams`.

## [2.1.2] - 2026-08-15

### Removed

- CMake-era Nix entrypoints (`default.nix`, `shell.nix`,
- `nix/yodaStruct.nix`). They still named the product `yodaStruct`. Build with pixi.

## [2.1.1] - 2026-08-15

### Fixed

- CHILL `isInterfacial` walks the four nearest neighbours, the same
- star as `c_ij`.
- CHILL `getIceType` / `getIceTypeNoPrint` send atoms without four
- recorded bonds to water, matching CHILL+.
- Seeded affiliation floods HC and DDC separately so an HC seed does
- not keep a DDC-only atom in the same H-bond component.

## [2.1.0] - 2026-08-15

### Changed

- The C++ engine is this repository. Front ends are not.

### Changed

- The `yodaStruct` Lua and Fennel CLI moved to
- https://github.com/d-SEAMS/yodaStruct . `-Dwith_lua=enabled` is now an error that names that repository.
- Python bindings remain in
- https://github.com/d-SEAMS/PydSEAMSlib (moved in 2.0.1).
- C++ `getCorrel`, `getCorrelPlus`, `getIceType`,
- `getIceTypeNoPrint`, `getIceTypePlus`, `getIceTypePlusNoPrint` and `reclassifyWater` return `void`. They already took the cloud by reference; the extra copy was the return. PydSEAMSlib still returns the object from `_core`.

### Added

- Incremental Franzblau rings, order-free HC/DDC affiliation, seeded
- hysteresis, Voronoi `c/2` certificate, scalar Steinhardt parameters, and the neighbour-list reverse-index that makes those paths linear.

## [2.0.1] - 2026-04-07

### Changed

- Build and CI fixes following v2.0.0 release.

### Fixed

- Python bindings moved to
- https://github.com/d-SEAMS/PydSEAMSlib . This repository is the C++ engine. The Lua CLI later moved to yodaStruct (2.1.0).
- Fix Eigen install path in wheel builds
- Fix macOS deployment target (bump to 14), remove --enable-new-dtags on macOS
- Regenerate pixi.lock for chemfiles compatibility
- Increase test timeout to 120s for bulkTUM tests on macOS
- Drop macOS x86_64 wheel builds (no macos-13 runners available)

### Changed

- Added release workflow for automatic Zenodo DOI minting on tag push

## [2.0.0] - 2026-03-22

### Changed

- This is a major release that modernizes the entire codebase.

### Changed

- Replaced Lua scripting interface with Python bindings via nanobind
- All C++ APIs converted from pointer to reference parameters
- Removed Nix build system in favor of Meson + pixi (conda-forge)
- License changed from GPL to MIT

### Added

- Python package `pydseamslib` with high-level `Trajectory` class
- O(n) cell-list neighbor search via vesin (optional, with brute-force fallback)
- SIMD-vectorized distance computation via Highway
- chemfiles integration for reading PDB, GRO, DCD, and 40+ trajectory formats
- readcon-core integration for reading .con format (eOn trajectories)
- Formal verification pipeline: SymPy symbolic proofs, Coq machine-checked
- proofs, Hypothesis property-based tests (32 properties)
- Comprehensive test suite: 21 Catch2 test binaries, 100% function coverage
- Binary wheels for Linux and macOS via cibuildwheel
- Cross-platform CI (GitHub Actions) for Linux, macOS, and coverage

### Fixed

- Fixed 5/9 incorrect elements in quaternion rotation matrix (quat2RotMatrix)
- Fixed division by zero in getAverageWithoutOutliers
- Fixed off-by-one in getAverageWithoutOutliers fallback loop
- Fixed XZ/YZ return order swap in projAreaSingleRing
- Fixed acos domain vulnerability in eigenVecAngle and angDistDegQuaternions
- Fixed iceCloud++ instead of iatom++ in cluster.cpp loop
- Fixed out-of-bounds in populateHbonds when molecule ID not in H-atom list
- Fixed out-of-bounds in populateHbonds when molecule has <2 hydrogen atoms
- Fixed writeLAMMPSdata accessing rings[0] before empty guard
- Fixed moleculesInSingleSlice overwriting inSlice from molecule-mate
- Fixed matchPrism/matchUntetheredPrism crash on failed basal ordering
- Fixed relOrderPrismBlock crash on symmetric (equidistant) prisms
- Fixed H-bond network duplicate entries from symmetric pair processing
- Fixed topoUnitMatchingBulk hardcoded template path (now configurable)

### Changed

- Required: Eigen 3.4+, BLAS, LAPACK, Catch2 3.6+, Meson 1.3+, nanobind 2.0+
- Optional: Highway 1.3+ (SIMD), vesin (cell-list), chemfiles 0.10+ (formats),
- readcon-core 0.5+ (.con format)
- Version 1.0.0 (2020)
- Initial release. See Goswami et al., J. Chem. Inf. Model. 2020, 60, 2169-2177.
