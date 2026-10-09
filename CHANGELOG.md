# Changelog

All notable changes to this project are documented in this file.

The format is [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Unreleased notes live in [`changelog.d/`](changelog.d/) and are assembled
by [towncrier](https://towncrier.readthedocs.io/).

<!-- towncrier release notes start -->

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
