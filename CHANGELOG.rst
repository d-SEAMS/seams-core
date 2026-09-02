=========
Changelog
=========

Unreleased
==========

``topo::matchLibraries`` names atoms against key libraries at several
hop counts, the deepest library that holds an atom's key winning, and
reports the depth that named each atom, so a molecule whose wide
neighbourhood is disturbed still gets its name from the inner shells.
Libraries record whether their keys carry vertex colours (``colours``
in the header; an absent field reads as uncoloured) and a coloured
library is refused against plain keys. ``seams fingerprint --library``
accepts a comma separated list and prints the per-depth counts.
``site::guestOccupancy`` places guests (methane, THF, ions) in
enumerated cages by the periodic centroid of each cage's vertices and
counts occupied, multiply occupied and free, the occupancy of a
clathrate hydrate; ``site::periodicCentroid`` is the helper.

Version 2.8.0 (2026-09-02)
===========================

``topo::localKey`` names the isomorphism class of an atom's rooted bonded
neighbourhood within a number of hops (the nauty certificate with the
centre in its own colour cell when nauty is linked, a Weisfeiler-Lehman
refinement hash otherwise); ``topo::fingerprint`` keys a frame by the
sorted refinement hashes and the primitive ring census, so relabelled
configurations share a key. ``site::ionEnvironment`` classes ions by
their first water shell against a per-atom ice flag. ``seams
fingerprint`` and ``seams ions`` expose both; ``seams cages --graph
seeded --complete`` turns on the ring completion; ``--hops``,
``--ion-types`` and ``--ion-cutoff`` are new options. Vertex colours
(atom types) partition the keys by species (``--colour-types``). A
``topo::KeyLibrary`` collects the keys of reference structures under
labels and names any atom whose key it holds; a key shared by several
references carries all their labels. ``seams fingerprint
--emit-library LABEL`` writes a frame's keys as library lines and
``--library FILE`` names the atoms of a frame by one. The readcon-core
fallback follows the engine's library kind, so a static engine embeds
the CON reader; CON and chemfiles input without their reader is refused
instead of returning an empty frame.

Version 2.7.0 (2026-09-02)
===========================

``ring::seededCageAffiliation`` takes an optional completion flag:
``ring::ringAdjacentCompletion`` fills the last vertex of any six-ring
whose other vertices carry a cage label, iterated to a fixed point,
separately for the HC and DDC labels. The rule is stated and proven in
``lean/`` (Mathlib): the completion is the least fixed point above the
seed, empty on an empty seed, sound, and independent of visiting
order; the edge-sharing rule it replaces is shown unsound on a
five-vertex instance. ``tests/walk_compare`` walks a trajectory and
prints per-frame CHILL+ and cage labels with their largest clusters.
The reproducibility workflow moved to the ``dseams2_repro`` package.

Version 2.6.0 (2026-08-17)
===========================

Meson links FlexiBLAS when pkg-config or libflexiblas is present,
and falls back to the conda-forge libblas/liblapack ABI packages.
The reproducibility campaign no longer configures in-tree Python;
bindings come from pydseamslib on PyPI. The nucleation notebook
reads CageScore fields. Figshare Lua demos resolve example_lua
from yodaStruct.

Version 2.5.0 (2026-08-16)
===========================

``--family`` (default ``waterIce``) is an input. Ice scores refuse
non-water families and name the family. ``seams cn --ions`` is cage
degree on ``site::ionCloud``; ``seams pairs`` is the mutual-nearest
contact-pair count, not ionicity. ``seams domains`` is the largest
polar or apolar Stoddard component. ``seams density-z`` is type-resolved
``rho(z)`` with slab volume from dump H. Running CN and the first
minimum of ``g_IJ`` come from ``rdf::runningCN`` /
``firstMinimumBin``. Dump writers emit ``xy, xz, yz`` when the cloud
carries tilt.

Version 2.4.0 (2026-08-16)
===========================

Partial 3D site-site ``g_IJ(r)`` and coordination numbers under
the dump MIC, with ``seams rdf`` and ``seams cn``. Site chemistry
is a ``site::Table`` of LAMMPS types and optional atom-ID
overrides; ``site::ionCloud`` collapses each ion ``molID`` to one
unwrapped COM vertex. Hydrogen bonds accept an explicit donor-H
index set (``populateHbondsFromDonors``). ``seams hbonds --donors``
uses every hydrogen as a donor candidate; water cutoffs stay on
that command and are not an ionic-liquid criterion. In-plane RDF
normalization uses the dump-cell volume ``|det H|`` and the
restricted-triclinic in-plane area ``lx*ly``, not the particle AABB
or bound-span product.

Version 2.3.4 (2026-08-16)
===========================

Steinhardt, hydrogen-bond, and cluster wraps use the dump
minimum image, not an independent-axis wrap on bound spans.

Version 2.3.3 (2026-08-16)
===========================

Cutoff fallbacks and the skin Verlet refresh use the dump 3x3
minimum image, not a diagonal of bound spans.

Version 2.3.2 (2026-08-16)
===========================

Cutoff lists take the LAMMPS dump 3x3, not a diagonal of bound
spans. Device ice-score batches pass ``lammpsBoxToLinkcell`` when
the dump tilt is present. ``linkcell`` wrap pin is v0.3.1.

Version 2.3.1 (2026-08-16)
===========================

``neighList``, ``halfNeighList``, and the in-plane RDF sampler use
the vesin cell list (same helper as ``neighListO``). Brute force
remains the fallback. The GPU ice-score workspace is unchanged.

Version 2.3.0 (2026-08-16)
===========================

The TUM ice score (mutual four-nearest graph, primitive hexagons,
HC/DDC) runs as a device-resident frame batch when gpulite and
``linkcell`` 0.3.0 are present. CHILL and :math:`q_{lm}` stay on the
host. Runtime knobs are twelve-factor: ``SEAMS_CONFIG`` or
``./seams.env``, then the environment, then CLI flags.
``seams --print-config`` dumps the table. Examples and the README
use ``seams`` / ``pydseams`` / ``require("dseams")``; the 2020
``yodaStruct -c`` / ``conf.yaml`` surface is gone.

The C++ API in the book is Doxygen via doxyrest (``api/index``).

Version 2.2.5 (2026-08-16)
===========================

``subprojects/linkcell.wrap`` is v0.2.4. A static ``libyodaLib.a``
no longer passes ``liblinkcell.so`` to ``ar``. Wrap consumers link
the static archive, so ``pydseams.yoda`` does not need
``liblinkcell`` at import time. Dump doubles parse with ``strtod``
on Apple libc++.

Version 2.2.4 (2026-08-16)
===========================

Periodic *k*-nearest graphs via ``d-SEAMS/linkcell`` v0.2.2.
``kNearestNeighbourList`` / ``kNearestNeighbourPair`` write packed
``n*k`` nominations. Mutual and union come from one walk. The LAMMPS
dump box is bound spans plus ``xy, xz, yz``; ``periodicDistSq``
recovers H. Empty frames no longer abort the walker.

Docs orgmode covers the ``seams`` CLI, the Nix flake, and the
three-repo split. Tutorials use the live APIs.

Version 2.2.3 (2026-08-15)
===========================

The ``seams`` CLI uses Argum, the same parser as eonclient. Help,
errors, ice-type counts, and ``--features`` are colorized. ``NO_COLOR``
turns the colors off.

Version 2.2.2 (2026-08-15)
===========================

The compiled Python surface is ``pydseams.yoda``. Docs and the
repro scripts use that name. ``_core`` remains an alias.

Version 2.2.1 (2026-08-15)
===========================

Added
-----
- Flake-based Nix package for ``libyodaLib`` and the ``seams`` CLI
  (meson, not the CMake-era ``yodaStruct`` derivation).
- ``nLammpsFrames`` / ``dropLammpsDumpIndex``: live dump session with a
  lazy ``ITEM: TIMESTEP`` offset table (LAMMPS ``ReaderNative`` cursor,
  chemfiles ``read_step``, readcon frame offsets). Sequential
  ``load_frame`` walks no longer rescan prior snapshots.
- ``forEachLammpsFrame``: OpenMP walk over a frame range. Each worker
  opens its own handle and seeks. ``seams --frame N --last M --jobs J``.
- ``SkinNeighborList``: vesin candidates at cutoff+skin, rebuilt on
  the Verlet trigger. ``BondGraph`` is chosen at runtime
  (``cutoff``, ``knn``, ``knn-union``). ``seams cages --graph``
  adds ``seeded``.

Fixed
-----
- LAMMPS readers bind ``xu``/``yu``/``zu`` and ``xs``/``ys``/``zs``
  when ``x y z`` are absent (Niu/Parrinello TIP4P/Ice dumps).

Version 2.2.0 (2026-08-15)
===========================

Added
-----
- ``seams`` CLI: ``read``, ``chill``, ``chill-plus``, ``cages``.
  This is the engine command line. Lua is the ``dseams`` library
  (yodaStruct repo). Python is ``pydseams``.

Version 2.1.2 (2026-08-15)
===========================

Removed
-------
- CMake-era Nix entrypoints (``default.nix``, ``shell.nix``,
  ``nix/yodaStruct.nix``). They still named the product ``yodaStruct``.
  Build with pixi.

Version 2.1.1 (2026-08-15)
===========================

Fixed
-----
- CHILL ``isInterfacial`` walks the four nearest neighbours, the same
  star as ``c_ij``.
- CHILL ``getIceType`` / ``getIceTypeNoPrint`` send atoms without four
  recorded bonds to water, matching CHILL+.
- Seeded affiliation floods HC and DDC separately so an HC seed does
  not keep a DDC-only atom in the same H-bond component.

Version 2.1.0 (2026-08-15)
===========================

The C++ engine is this repository. Front ends are not.

Breaking
--------
- The ``yodaStruct`` Lua and Fennel CLI moved to
  https://github.com/d-SEAMS/yodaStruct . ``-Dwith_lua=enabled`` is now
  an error that names that repository.
- Python bindings remain in
  https://github.com/d-SEAMS/PydSEAMSlib (moved in 2.0.1).
- C++ ``getCorrel``, ``getCorrelPlus``, ``getIceType``,
  ``getIceTypeNoPrint``, ``getIceTypePlus``, ``getIceTypePlusNoPrint``
  and ``reclassifyWater`` return ``void``. They already took the cloud
  by reference; the extra copy was the return. PydSEAMSlib still
  returns the object from ``_core``.

Added
-----
- Incremental Franzblau rings, order-free HC/DDC affiliation, seeded
  hysteresis, Voronoi ``c/2`` certificate, scalar Steinhardt parameters,
  and the neighbour-list reverse-index that makes those paths linear.

Version 2.0.1 (2026-04-07)
===========================

Build and CI fixes following v2.0.0 release.

Fixed
------
- Python bindings moved to
  https://github.com/d-SEAMS/PydSEAMSlib . This repository is the C++
  engine. The Lua CLI later moved to yodaStruct (2.1.0).
- Fix Eigen install path in wheel builds
- Fix macOS deployment target (bump to 14), remove --enable-new-dtags on macOS
- Regenerate pixi.lock for chemfiles compatibility
- Increase test timeout to 120s for bulkTUM tests on macOS
- Drop macOS x86_64 wheel builds (no macos-13 runners available)

Developer
----------
- Added release workflow for automatic Zenodo DOI minting on tag push

Version 2.0.0 (2026-03-22)
===========================

This is a major release that modernizes the entire codebase.

Breaking Changes
-----------------
- Replaced Lua scripting interface with Python bindings via nanobind
- All C++ APIs converted from pointer to reference parameters
- Removed Nix build system in favor of Meson + pixi (conda-forge)
- License changed from GPL to MIT

New Features
-------------
- Python package ``pydseamslib`` with high-level ``Trajectory`` class
- O(n) cell-list neighbor search via vesin (optional, with brute-force fallback)
- SIMD-vectorized distance computation via Highway
- chemfiles integration for reading PDB, GRO, DCD, and 40+ trajectory formats
- readcon-core integration for reading .con format (eOn trajectories)
- Formal verification pipeline: SymPy symbolic proofs, Coq machine-checked
  proofs, Hypothesis property-based tests (32 properties)
- Comprehensive test suite: 21 Catch2 test binaries, 100% function coverage
- Binary wheels for Linux and macOS via cibuildwheel
- Cross-platform CI (GitHub Actions) for Linux, macOS, and coverage

Bug Fixes
----------
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

Dependencies
-------------
- Required: Eigen 3.4+, BLAS, LAPACK, Catch2 3.6+, Meson 1.3+, nanobind 2.0+
- Optional: Highway 1.3+ (SIMD), vesin (cell-list), chemfiles 0.10+ (formats),
  readcon-core 0.5+ (.con format)

Version 1.0.0 (2020)
=====================

Initial release. See Goswami et al., J. Chem. Inf. Model. 2020, 60, 2169-2177.
