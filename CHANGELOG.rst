=========
Changelog
=========

Version 2.0.1 (2026-04-07)
===========================

Build and CI fixes following v2.0.0 release.

Fixed
------
- Python extension builds against the CPython 3.12 stable ABI
  (=limited_api: '3.12'=). Wheels are one abi3 artifact per platform
  (GIL). Free-threaded CPython has no limited ABI and is not a wheel
  target.
- =Frame= / =read= / =from_ase= / =from_arrays= are the public Python
  API. CHILL and CHILL+ classify without writing scratch files.
  =Trajectory= remains an alias of =Frame=.
- Pin Python to 3.12 for the abi3 build host
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
