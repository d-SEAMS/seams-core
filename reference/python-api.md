<!-- https://docs.dseams.info/reference/python-api.md -->

# Python {#python}

`pydseamslib` provides `import dseams` and `import pydseams`. Signatures, parameters, returns, and examples are read from the docstrings in [d-SEAMS/PydSEAMSlib](https://github.com/d-SEAMS/PydSEAMSlib). The C column is taken from the nanobind definitions in `src/dseams/_core.cpp` and from `seams_*` calls in the function.

## Native API {#native-api}

Generated from the sources named below. A name on this page is a public call in that language.

### dseams.frames.cages {#dseams-frames-cages}

```python
dseams.frames.cages(positions, box, cutoff=3.5, graph='seeded', k=4)
```

Per-atom cage code and the four population counts.

### dseams.frames.rings {#dseams-frames-rings}

```python
dseams.frames.rings(positions, box, cutoff=3.5, max_depth=6, graph='cutoff')
```

Primitive-ring histogram. Index ``i`` counts rings of size ``i``.

### dseams.frames.signature_cages {#dseams-frames-signature-cages}

```python
dseams.frames.signature_cages(positions, box, spec='512', cutoff=3.5, graph='cutoff')
```

Closed cages whose ring census matches ``spec``. Returns (n_cages, n_atoms).

### dseams.frames.symmetry {#dseams-frames-symmetry}

```python
dseams.frames.symmetry(positions, box, shell_rmsd=0.12)
```

Per-atom shell label. 0 ice Ih, 1 ice Ic, 2 stacking-disordered, -1 otherwise.

### dseams.userframe.CageScore {#dseams-userframe-cagescore}

```python
dseams.userframe.CageScore
```

Per-atom hexagonal-cage and double-diamond flags. ``n_ih`` counts atoms in a hexagonal cage, ``n_ic`` atoms in a double-diamond cage, and ``n_water`` atoms in neither. An atom in both contributes to ``n_ih`` and ``n_ic``.

### dseams.userframe.CageScore.n_both {#dseams-userframe-cagescore-n-both}

```python
dseams.userframe.CageScore.n_both()
```

### dseams.userframe.CageScore.n_ic {#dseams-userframe-cagescore-n-ic}

```python
dseams.userframe.CageScore.n_ic()
```

### dseams.userframe.CageScore.n_ih {#dseams-userframe-cagescore-n-ih}

```python
dseams.userframe.CageScore.n_ih()
```

### dseams.userframe.CageScore.n_water {#dseams-userframe-cagescore-n-water}

```python
dseams.userframe.CageScore.n_water()
```

### dseams.userframe.Frame {#dseams-userframe-frame}

```python
dseams.userframe.Frame
```

One configuration. Classification calls libseams_lang.

### dseams.userframe.Frame.box {#dseams-userframe-frame-box}

```python
dseams.userframe.Frame.box()
```

### dseams.userframe.Frame.cages {#dseams-userframe-frame-cages}

```python
dseams.userframe.Frame.cages(seeded=True, k=4, ring_adjacent=False)
```

Hexagonal and double-diamond membership from seams_frame_cages.

### dseams.userframe.Frame.cages_by_signature {#dseams-userframe-frame-cages-by-signature}

```python
dseams.userframe.Frame.cages_by_signature(signature, max_ring_size=None, edge=None)
```

Closed cages of one ring census. Each dict has ``signature`` and ``vertices`` (atom indices of that cage). ``n_atoms`` on the list is the number of distinct vertices across the cages. ``edge`` is the face-edge window. ``None`` keeps ``2.65,3.15`` in the coordinate units of the file. ``”none”`` leaves the cage bound to the face census, the same string as ``–cage-edge``.

### dseams.userframe.Frame.chill_plus {#dseams-userframe-frame-chill-plus}

```python
dseams.userframe.Frame.chill_plus(graph='cutoff')
```

Published CHILL+ histogram. Index names match the C label order.

### dseams.userframe.Frame.close {#dseams-userframe-frame-close}

```python
dseams.userframe.Frame.close()
```

### dseams.userframe.Frame.find_prisms {#dseams-userframe-frame-find-prisms}

```python
dseams.userframe.Frame.find_prisms(output_dir='output/', max_depth=6, shape_matching=False)
```

Prism census from seams_frame_prisms. The return is the C rows. A frame with no atoms returns one zero row for each ring size from 3 through ``max_depth``. The engine reads the first point while removing axial drift, so this method does not call it.

### dseams.userframe.Frame.monolayer_rings {#dseams-userframe-frame-monolayer-rings}

```python
dseams.userframe.Frame.monolayer_rings(output_dir, sheet_area, max_depth=4)
```

Quasi-2D ring counts from seams_frame_monolayer. Coverage is the number the C writer stores in ``topoMonolayer/coverageAreaXY.dat``.

### dseams.userframe.Frame.n_atoms {#dseams-userframe-frame-n-atoms}

```python
dseams.userframe.Frame.n_atoms()
```

### dseams.userframe.Frame.positions {#dseams-userframe-frame-positions}

```python
dseams.userframe.Frame.positions()
```

Coordinates, shape ``(n, 3)``. A frame with no atoms is shape ``(0, 3)``.

### dseams.userframe.Frame.rdf_2d {#dseams-userframe-frame-rdf-2d}

```python
dseams.userframe.Frame.rdf_2d(output_dir, cutoff=12.0, binwidth=0.05)
```

In-plane g(r) from seams_frame_rdf2d. A frame with no atoms returns the bin centres and a zero ``g(r)`` of that length. The engine reads the first point of an empty cloud, so this method does not call it.

### dseams.userframe.Frame.ring_census {#dseams-userframe-frame-ring-census}

```python
dseams.userframe.Frame.ring_census(max_size=8, graph=None)
```

Primitive rings by size. ``counts\[i\]`` from seams_frame_rings. ``graph`` is ``”hbond”`` or ``”cutoff”``. The default is the frame’s bonded list: hydrogen bonds for a type-2 oxygen dump, and the cutoff list otherwise.

### dseams.userframe.Frame.steinhardt {#dseams-userframe-frame-steinhardt}

```python
dseams.userframe.Frame.steinhardt(order_l=6)
```

Local q_l and the Lechner-Dellago average from seams_frame_steinhardt.

### dseams.userframe.Frame.to_ase {#dseams-userframe-frame-to-ase}

```python
dseams.userframe.Frame.to_ase()
```

ASE Atoms for this frame, with ice labels and cage flags when set.

### dseams.userframe.Frame.to_solvis {#dseams-userframe-frame-to-solvis}

```python
dseams.userframe.Frame.to_solvis(expand_box=True)
```

solvis System around :meth:`to_ase`. solvis-tools must be installed.

### dseams.userframe.IceCounts {#dseams-userframe-icecounts}

```python
dseams.userframe.IceCounts
```

CHILL+ histogram. Missing names read as 0.

### dseams.userframe.from_ase {#dseams-userframe-from-ase}

```python
dseams.userframe.from_ase(atoms, select='O', cutoff=None, bonded='auto')
```

Frame for one ASE species. The coordinates go through seams_frame_from_arrays.

### dseams.userframe.read {#dseams-userframe-read}

```python
dseams.userframe.read(filename, frame=1, atom_type=None, all_atoms=False, **kwargs)
```

Read one frame with seams_frame_read. ``atom_type`` keeps that LAMMPS species. ``all_atoms`` keeps every atom. Otherwise the library guesses oxygen, which is what the tutorials do when they omit the type.

### dseams.userframe.to_solvis {#dseams-userframe-to-solvis}

```python
dseams.userframe.to_solvis(frame, expand_box=True)
```

### pydseams.aseio.frame_from_ase {#pydseams-aseio-frame-from-ase}

```python
pydseams.aseio.frame_from_ase(atoms, select='O', cutoff=3.5, bonded='auto')
```

Construct ``cls`` (a :class:`~pydseams.frame.Frame`) from ASE ``Atoms``. Parameters cls : type :class:`~pydseams.frame.Frame` or a subclass. atoms : ase.Atoms Configuration with a nonsingular cell periodic in all three directions. select : str, int or sequence, optional Chemical symbol or atomic number kept for analysis. Default ``”O”``. ``None`` keeps every atom. A sequence such as ``(“O”, “Na”, “Cl”)`` keeps every listed species in the cloud and analyses the first; the others stay as their atomic numbers in ``c_type``, which is what :func:`pydseams.features.ion_environment` reads. cutoff : float, optional Neighbour cutoff in Angstroms. Default ``3.5``. bonded : {“auto”, “hbond”, “cutoff”}, optional Graph for rings. ``”auto”`` becomes ``”hbond”`` when the ``Atoms`` contain hydrogen and the selected analysis cloud excludes hydrogen; otherwise it becomes ``”cutoff”``. Explicit ``”hbond”`` also requires a heavy-atom selection. Returns Frame Selected-species cloud plus an optional hydrogen cloud for hydrogen-bond analysis. Raises ImportError If ASE is not installed. TypeError If ``atoms`` has no ``get_positions``. ValueError If the cell is singular, is not fully periodic, or ``select`` matches no atom, or if ``bonded=”hbond”`` selects hydrogen. Notes For a heavy-atom selection, hydrogens stay in a side cloud so the analysed species remain the CHILL / ring particles. General cells are rotated into LAMMPS restricted form for analysis and restored by :func:`frame_to_ase`. An ASE ``mol-id`` array supplies molecule ownership. Without it, each hydrogen belongs to the nearest selected atom under the minimum-image convention.

### pydseams.aseio.frame_to_ase {#pydseams-aseio-frame-to-ase}

```python
pydseams.aseio.frame_to_ase(frame)
```

Build an ASE ``Atoms`` from a :class:`~pydseams.frame.Frame`. Parameters frame : Frame Configuration to export. Returns ase.Atoms Periodic cell in its imported orientation. Symbols come from the ASE import when present; LAMMPS-only frames fall back to ``O``. Notes After :meth:`~pydseams.frame.Frame.chill_plus` (or :meth:`~pydseams.frame.Frame.chill`), ``atoms.arrays\[‘ice_type’\]`` holds the per-atom labels. After :meth:`~pydseams.frame.Frame.cages`, ``atoms.arrays\[‘hc’\]`` and ``atoms.arrays\[‘ddc’\]`` hold the cage flags. ``atoms.info\[‘dseams_n_atoms’\]`` is the analysed particle count.

### pydseams.band.band_rows {#pydseams-band-band-rows}

```python
pydseams.band.band_rows(frames, max_ring_size=8)
```

One row per frame: bonds that entered or left, shortest same-type distances, rings. Parameters frames : sequence of Frame Images in path order with the same atom IDs and types in the same order. Each frame uses its own cutoff graph and periodicity. max_ring_size : int, optional Largest primitive ring. Default 8. Returns list of dict ``index``, ``n_bonds``, ``entered``, ``left``, ``shortest_within_type``, ``rings``. The first row has empty ``entered`` and ``left``.

### pydseams.band.bond_pairs {#pydseams-band-bond-pairs}

```python
pydseams.band.bond_pairs(frame)
```

Undirected pairs in the cutoff graph, as index tuples ``(i, j)`` with ``i < j``.

### pydseams.band.ring_census {#pydseams-band-ring-census}

```python
pydseams.band.ring_census(frame, max_size=8)
```

Count primitive rings by size, up to ``max_size``.

### pydseams.band.shortest_within_type {#pydseams-band-shortest-within-type}

```python
pydseams.band.shortest_within_type(frame, type_code)
```

Shortest distance between two atoms of ``type_code``. Periodic frames use the nearest image in the simulation cell; open frames use the Cartesian distance. Returns ``None`` when fewer than two atoms carry that type.

### pydseams.bond_network.components {#pydseams-bond-network-components}

```python
pydseams.bond_network.components(n_sites, pairs)
```

Return per-site component labels and component sizes. Each bond joins two zero-based site indices regardless of direction. Labels follow the first site encountered in each component. Isolated sites have size one.

### pydseams.bond_network.residence_counts {#pydseams-bond-network-residence-counts}

```python
pydseams.bond_network.residence_counts(membership, lags)
```

Return lag, at-risk origin, continuous and intermittent counts. The boolean matrix has rows at equally spaced saved frames and columns for stable site identities. An origin counts only when its endpoint exists. Continuous residence requires membership in every saved frame in the interval. Intermittent residence requires both endpoints.

### pydseams.config.apply {#pydseams-config-apply}

```python
pydseams.config.apply()
```

Load the backing file once. Already-set env vars stay put.

### pydseams.config.apply_file {#pydseams-config-apply-file}

```python
pydseams.config.apply_file(path, required=True)
```

### pydseams.config.cutoff {#pydseams-config-cutoff}

```python
pydseams.config.cutoff(default=3.5)
```

### pydseams.config.frame {#pydseams-config-frame}

```python
pydseams.config.frame(default=1)
```

### pydseams.config.get_float {#pydseams-config-get-float}

```python
pydseams.config.get_float(name, default)
```

### pydseams.config.get_int {#pydseams-config-get-int}

```python
pydseams.config.get_int(name, default)
```

### pydseams.config.graph {#pydseams-config-graph}

```python
pydseams.config.graph(default='seeded')
```

### pydseams.config.k_neighbors {#pydseams-config-k-neighbors}

```python
pydseams.config.k_neighbors(default=4)
```

### pydseams.config.reset {#pydseams-config-reset}

```python
pydseams.config.reset()
```

Test helper: allow apply() to run again.

### pydseams.features.IceFeaturizer {#pydseams-features-icefeaturizer}

```python
pydseams.features.IceFeaturizer
```

Turn a trajectory into per-frame features and per-molecule states. Parameters trajectory : Trajectory An open :class:`~pydseams.frame.Trajectory`; frames are visited with :meth:`~pydseams.frame.Frame.load_frame`, so the incremental updaters stay warm. k : int, optional Neighbours in the four-nearest graphs. Default ``4``. ring_adjacent : bool, optional Ring-adjacent completion of the seeded assignment. Default ``True``. chill : bool, optional Also compute CHILL+ on the cutoff graph. Default ``True``. ion_types : iterable of int, optional ``c_type`` codes of ions present in the cloud. When given, the vector gains :data:`ION_FEATURE_NAMES` from :func:`ion_features`.

### pydseams.features.IceFeaturizer.feature_names {#pydseams-features-icefeaturizer-feature-names}

```python
pydseams.features.IceFeaturizer.feature_names()
```

### pydseams.features.IceFeaturizer.frame_features {#pydseams-features-icefeaturizer-frame-features}

```python
pydseams.features.IceFeaturizer.frame_features()
```

Feature vector and per-molecule states of the loaded frame. Returns (features, states) ``features`` is a float array of length ``len(self.feature_names)``; ``states`` an ``int8`` array with one of ``STATE_*`` per molecule.

### pydseams.features.IceFeaturizer.transform {#pydseams-features-icefeaturizer-transform}

```python
pydseams.features.IceFeaturizer.transform(frames=None)
```

Features and states for a range of frames. Parameters frames : iterable of int or None 1-based frame numbers; ``None`` walks the stored frames from the current one to the last. Returns (X, S) ``X`` has shape ``(n_frames, len(self.feature_names))``; ``S`` has shape ``(n_frames, n_molecules)`` with ``int8`` states.

### pydseams.features.count_frames {#pydseams-features-count-frames}

```python
pydseams.features.count_frames(filename)
```

Stored frames in a LAMMPS dump (``ITEM: TIMESTEP`` headers).

### pydseams.features.discretize_nmax {#pydseams-features-discretize-nmax}

```python
pydseams.features.discretize_nmax(n_max, edges)
```

Map a largest-cluster series onto integer states by bin edges. Parameters n_max : array_like Largest cage cluster per frame (``X\[:, 1\]``). edges : sequence of float Increasing bin edges; state ``i`` is ``edges\[i\] <= n_max < edges\[i+1\]``, below the first edge is state 0, at or above the last is the last. Returns numpy.ndarray of int A discrete trajectory for :class:`deeptime.markov.msm.MaximumLikelihoodMSM` or ``pyemma.msm.estimate_markov_model``.

### pydseams.features.ion_environment {#pydseams-features-ion-environment}

```python
pydseams.features.ion_environment(frame, states, ion_types, cutoff=None)
```

Ice content of the first water shell of every ion. Ions are not part of the hydrogen-bond network; the cage assignment runs on the water alone and the ions are read against it. A frame built with ``all_atoms=True`` or :meth:`Frame.from_arrays` carries every species in its cloud; ``frame.atom_type`` names the water. The engine’s ``ionEnvironment`` does the work when the compiled module exposes it; the NumPy path below gives the same answer on older engines. Parameters frame : Frame Frame whose cloud holds water and ions. states : array of int8 Per-atom states from :meth:`IceFeaturizer.frame_features`; non-water entries are ignored. ion_types : iterable of int ``c_type`` codes of the ions. cutoff : float, optional First-shell radius in Angstrom. Default ``frame.cutoff``. Returns (indices, shell, fraction, ion_states) Cloud indices of the ions, the water count in each first shell, the fraction of that shell carrying an ice label, and one of ``ION_ICE`` (whole shell ice), ``ION_LIQUID`` (no ice) or ``ION_FRONT`` per ion. An ion with an empty shell is liquid.

### pydseams.features.ion_features {#pydseams-features-ion-features}

```python
pydseams.features.ion_features(frame, states, ion_types, cutoff=None)
```

Per-frame ion summary in the order of :data:`ION_FEATURE_NAMES`.

### pydseams.features.to_deeptime {#pydseams-features-to-deeptime}

```python
pydseams.features.to_deeptime(X, lagtime, dim=None)
```

TICA of a per-frame feature array with deeptime. Returns the fitted :class:`deeptime.decomposition.TICA` model; call ``model.transform(X)`` for the projection.

### pydseams.features.to_pyemma_featurizer {#pydseams-features-to-pyemma-featurizer}

```python
pydseams.features.to_pyemma_featurizer(featurizer, topology)
```

Register the per-frame vector as a PyEMMA custom feature. PyEMMA’s featurizer wants a function of an ``mdtraj`` trajectory; this ignores the coordinates it is given and reads the frames through the engine instead, so the feature order is the one of :attr:`FEATURE_NAMES`.

### pydseams.frame.CageScore {#pydseams-frame-cagescore}

```python
pydseams.frame.CageScore
```

Per-atom hexagonal-cage (HC) and double-diamond-cage (DDC) flags. A molecule in an HC is ice Ih; a molecule in a DDC is ice Ic. Membership is a boolean per analysed atom. Parameters hc : sequence of bool True when the atom belongs to at least one hexagonal cage. ddc : sequence of bool True when the atom belongs to at least one double-diamond cage. Attributes hc : list of bool Per-atom HC membership. ddc : list of bool Per-atom DDC membership.

### pydseams.frame.CageScore.n_ic {#pydseams-frame-cagescore-n-ic}

```python
pydseams.frame.CageScore.n_ic()
```

Number of atoms flagged DDC (ice Ic).

### pydseams.frame.CageScore.n_ih {#pydseams-frame-cagescore-n-ih}

```python
pydseams.frame.CageScore.n_ih()
```

Number of atoms flagged HC (ice Ih).

### pydseams.frame.CageScore.n_water {#pydseams-frame-cagescore-n-water}

```python
pydseams.frame.CageScore.n_water()
```

Number of atoms in neither cage.

### pydseams.frame.ContactPairs {#pydseams-frame-contactpairs}

```python
pydseams.frame.ContactPairs
```

Mutual nearest unlike ion pairs returned by :meth:`Frame.pairs`.

### pydseams.frame.DensityProfile {#pydseams-frame-densityprofile}

```python
pydseams.frame.DensityProfile
```

Cartesian number-density profile returned by :meth:`Frame.density`.

### pydseams.frame.DomainStats {#pydseams-frame-domainstats}

```python
pydseams.frame.DomainStats
```

Largest connected site domain returned by :meth:`Frame.domain`.

### pydseams.frame.DualOccupancy {#pydseams-frame-dualoccupancy}

```python
pydseams.frame.DualOccupancy
```

Radius and inside occupancy of the same guests in the same cages.

### pydseams.frame.Frame {#pydseams-frame-frame}

```python
pydseams.frame.Frame
```

One configuration: neighbours, rings, CHILL(+), and cage membership. Load a LAMMPS dump, an ASE ``Atoms``, or raw arrays, then call :meth:`chill_plus` or :meth:`cages`. Classification does not write files. Prefer :func:`pydseams.io.read`, :meth:`from_ase`, or :meth:`from_arrays` over constructing this class by filename. Parameters filename : path-like, optional LAMMPS dump. Types 1 and 2 are treated as hydrogen and oxygen unless ``atom_type`` is set. frame : int, optional 1-indexed frame. Default ``1``. atom_type : int or None, optional Species to analyse. ``None`` picks oxygen (type 2) if that type is present, otherwise type 1 (mW-style single-site dumps). cutoff : float, optional Neighbour cutoff in Angstroms. Default ``3.5``. bonded : {“auto”, “hbond”, “cutoff”}, optional Graph for rings. ``”auto”`` uses hydrogen bonds when hydrogens are available, otherwise the cutoff neighbour list. region : ((xlo, ylo, zlo), (xhi, yhi, zhi)) or None, optional Rectangular slice passed to :func:`pydseams.yoda.readLammpsTrjreduced`. An axis with ``lo == hi`` is unconstrained. The cloud’s ``nop`` is the kept count. all_atoms : bool, optional Keep every LAMMPS atom type instead of filtering to ``atom_type``. ``atom_type`` still selects the species used by neighbour and ice analyses. Cannot be combined with ``region``. Raises ValueError If ``bonded`` is not one of ``auto``, ``hbond``, ``cutoff``, or if a LAMMPS dump has no atoms of type 1 or 2. TypeError If neither a filename nor a pre-built cloud is supplied. Notes Keyword-only ``cloud``, ``h_cloud``, and ``symbols`` are the constructor plumbing used by :meth:`from_ase` and :meth:`from_arrays`.

### pydseams.frame.Frame.bonds_by_index {#pydseams-frame-frame-bonds-by-index}

```python
pydseams.frame.Frame.bonds_by_index()
```

Index-based bonded graph used for rings. Hydrogen-bond list when :attr:`bonded` is ``”hbond”``, otherwise the cutoff neighbour list. Converted with :func:`pydseams.yoda.neighbourListByIndex`.

### pydseams.frame.Frame.box {#pydseams-frame-frame-box}

```python
pydseams.frame.Frame.box()
```

Engine box as lengths or LAMMPS restricted-triclinic values. Orthorhombic boxes have three lengths. Triclinic boxes have dump bound spans followed by ``xy``, ``xz``, and ``yz`` tilts.

### pydseams.frame.Frame.cage_affiliation {#pydseams-frame-frame-cage-affiliation}

```python
pydseams.frame.Frame.cage_affiliation()
```

Order-free per-ring HC/DDC flags on this frame’s six-rings. Uses :class:`pydseams.yoda.AffiliationUpdater` on the cutoff-or-hbond six-rings. Returns dict ``six_rings`` (list of 6-cycles), ``hc`` and ``ddc`` (per-ring bools), and ``reclassified`` (updater delta).

### pydseams.frame.Frame.cages {#pydseams-frame-frame-cages}

```python
pydseams.frame.Frame.cages(seeded=True, k=4, candidate_cutoff=None, ring_adjacent=False)
```

Ice score: HC = Ih, DDC = Ic, neither = water. Parameters seeded : bool, optional ``True`` (default) is the hysteresis construction: mutual four-nearest seeds, union-graph completion (:meth:`seeded_affiliation`). ``False`` is cutoff-graph affiliation on this frame’s six-rings (:meth:`cage_affiliation`). k : int, optional Neighbours kept in the seeded k-nearest graphs. Default ``4``. ring_adjacent : bool, optional Ring completion of the seeded assignment (fill the last vertex of a six-ring whose other vertices carry a label). Default ``False``. candidate_cutoff : float or None, optional Candidate-list cutoff for the k-nearest graphs. ``None`` uses :attr:`cutoff` ``+ 1.5``. Returns CageScore Per-atom HC/DDC flags. Cached as ``_cages`` for :meth:`to_ase`.

### pydseams.frame.Frame.cages_by_signature {#pydseams-frame-frame-cages-by-signature}

```python
pydseams.frame.Frame.cages_by_signature(signature, max_ring_size=None)
```

Closed polyhedra matching a ring-size census. Parameters signature : str Comma list (``4:6,6:8``) or named table entry (``sodalite``, ``alpha``, ``512``, ``51262``, ``hc``, ``ddc``). Named ``hc`` and ``ddc`` use the TUM finders. max_ring_size : int or None, optional Franzblau depth. ``None`` uses 8, or the largest size in a comma list. Returns list of dict Each dict has ``signature``, ``faces``, ``vertices``, and ``certificate``.

### pydseams.frame.Frame.chill {#pydseams-frame-frame-chill}

```python
pydseams.frame.Frame.chill()
```

CHILL labels for every analysed atom. Does not write a file. Calls :func:`pydseams.yoda.getCorrel` then :func:`pydseams.yoda.getIceTypeNoPrint` on :attr:`neighbor_list`. Mutates ``cloud.pts\[\].iceType``. Returns IceCounts Histogram of CHILL labels on this frame.

### pydseams.frame.Frame.chill_plus {#pydseams-frame-frame-chill-plus}

```python
pydseams.frame.Frame.chill_plus()
```

CHILL+ labels for every analysed atom. Does not write a file. Calls :func:`pydseams.yoda.getCorrelPlus` then :func:`pydseams.yoda.getIceTypePlusNoPrint` on :attr:`neighbor_list`. Mutates ``cloud.pts\[\].iceType``. Returns IceCounts Histogram of CHILL+ labels on this frame.

### pydseams.frame.Frame.classify_chill {#pydseams-frame-frame-classify-chill}

```python
pydseams.frame.Frame.classify_chill()
```

Alias of :meth:`chill`.

### pydseams.frame.Frame.classify_chill_plus {#pydseams-frame-frame-classify-chill-plus}

```python
pydseams.frame.Frame.classify_chill_plus()
```

Alias of :meth:`chill_plus`.

### pydseams.frame.Frame.classify_templates {#pydseams-frame-frame-classify-templates}

```python
pydseams.frame.Frame.classify_templates(k_neigh=12)
```

IRA/Horn overlay onto FCC, HCP, BCC, and SC neighbour shells. Parameters k_neigh : int, optional Neighbours in each template shell. Default ``12``. Returns list of dict Each hit has ``name``, ``rmsd``, and ``kind``.

### pydseams.frame.Frame.classify_topology {#pydseams-frame-frame-classify-topology}

```python
pydseams.frame.Frame.classify_topology(library, hops=2, max_ring_size=7, colour_types=False)
```

Name every analysed atom by one key library, or by several at different hop counts. Parameters library : pydseams.yoda.KeyLibrary, str, or a sequence of them A library or its text form from ``yoda.writeLibrary``. A sequence of libraries built at different ``hops`` names each atom by the deepest library that holds its key, so a molecule whose wide neighbourhood is disturbed still gets a name from its inner shells; ``hops`` is then ignored. hops, max_ring_size, colour_types As in :meth:`fingerprint`; ``colour_types`` must match the colouring the libraries were built with. Returns pydseams.yoda.LibraryMatch ``labels`` per atom (``””`` when no reference matches), ``counts`` per label, ``depth`` per atom (the hops of the library that named it, ``0`` when none) and ``matched``.

### pydseams.frame.Frame.cn {#pydseams-frame-frame-cn}

```python
pydseams.frame.Frame.cn(type_i, type_j, cutoff, binwidth=0.05)
```

Site-site coordination number integrated to ``cutoff``. Uses :func:`pydseams.yoda.partialRdfHist` and :func:`pydseams.yoda.coordinationNumber` with ``rho_J = nJ / dumpVolume``.

### pydseams.frame.Frame.density {#pydseams-frame-frame-density}

```python
pydseams.frame.Frame.density(bins=None, axis='z', atom_type=0)
```

Cartesian number density by particle type or mapped site kind. Parameters bins : int or None, optional Number of equal-width slabs. ``None`` uses about 0.1 Angstrom spacing along the selected restricted-cell axis. axis : {“x”, “y”, “z”, 0, 1, 2}, optional Profile axis. Default ``”z”``. atom_type : int, optional Particle type to count; ``0`` counts every particle. Ignored when ``table`` and ``kind`` are supplied. table, kind : optional A :class:`pydseams.yoda.SiteTable` and matching :class:`pydseams.yoda.SiteKind` for a chemistry-resolved profile. Returns DensityProfile

### pydseams.frame.Frame.domain {#pydseams-frame-frame-domain}

```python
pydseams.frame.Frame.domain(table, kind, cutoff=None)
```

Largest cutoff-connected component of a mapped site subset.

### pydseams.frame.Frame.find_prisms {#pydseams-frame-frame-find-prisms}

```python
pydseams.frame.Frame.find_prisms(output_dir='output/', max_depth=6, shape_matching=False)
```

Identify prism blocks and write engine output under ``output_dir``. Parameters output_dir : path-like, optional Directory passed to :func:`pydseams.yoda.prismAnalysis`. Default ``”output/”``. max_depth : int, optional Largest ring size considered. Default ``6``. shape_matching : bool, optional Enable TUM shape matching inside the prism search. Notes This path writes files. :meth:`chill_plus` and :meth:`cages` do not.

### pydseams.frame.Frame.fingerprint {#pydseams-frame-frame-fingerprint}

```python
pydseams.frame.Frame.fingerprint(hops=2, max_ring_size=7, colour_types=False)
```

Label-independent topology keys of the bonded graph. Parameters hops : int, optional Bonds from the centre in each local key. Default ``2``. max_ring_size : int, optional Largest primitive ring counted in the census. Default ``7``. colour_types : bool, optional Colour vertices by ``c_type`` (atom type or atomic number), so species never match across types. Default ``False``. Returns pydseams.yoda.FrameFingerprint ``key`` names the frame (same for any relabelling of the same bonded graph), ``atomKeys`` one class per atom, ``classes`` the histogram, ``ringCensus\[s\]`` the primitive rings of size ``s``. ``method`` is ``”nauty”`` when the engine links nauty, else ``”wl”``.

### pydseams.frame.Frame.fit_classifier {#pydseams-frame-frame-fit-classifier}

```python
pydseams.frame.Frame.fit_classifier(X, y, labels=None)
```

Fit a :class:`pydseams.yoda.LinearClassifier` on feature rows. Parameters X : sequence of sequence of float Feature matrix. y : sequence of int Integer class labels. labels : sequence of str, optional Human-readable class names stored on the classifier. Returns pydseams.yoda.LinearClassifier Fitted classifier, also stored on the frame.

### pydseams.frame.Frame.from_arrays {#pydseams-frame-frame-from-arrays}

```python
pydseams.frame.Frame.from_arrays(positions, cell, numbers=None, cutoff=None, bonded='cutoff', box_low=None)
```

Build a frame from ``(N, 3)`` positions and a periodic box. Parameters positions : sequence of (x, y, z) Cartesian coordinates. cell : sequence of float Orthorhombic lengths ``\[lx, ly, lz\]`` or LAMMPS restricted triclinic bound spans and tilts ``\[xspan, yspan, zspan, xy, xz, yz\]``. numbers : sequence of int, optional Per-atom type codes stored as ``c_type``. Default ``1``. cutoff : float, optional Neighbour cutoff in Angstroms. Default ``3.5``. bonded : {“auto”, “hbond”, “cutoff”}, optional Graph for rings. Default ``”cutoff”``. box_low : sequence of float, optional Box origin. Default ``(0, 0, 0)``. mol_ids : sequence of int, optional Molecule ID for each analysed atom. Defaults to 1 through N. hydrogen_positions : sequence of (x, y, z), optional Hydrogen coordinates in the same Cartesian frame as positions. hydrogen_mol_ids : sequence of int, optional Parent molecule IDs, one per hydrogen. Required with hydrogen_positions. Returns Frame Raises ValueError If ``positions`` is empty or ``cell`` is not three or six values.

### pydseams.frame.Frame.from_ase {#pydseams-frame-frame-from-ase}

```python
pydseams.frame.Frame.from_ase(atoms, select='O', cutoff=None, bonded='auto')
```

Load an ASE ``Atoms``. ``select`` is a symbol, atomic number, sequence, or None. Parameters atoms : ase.Atoms Configuration with a nonsingular cell periodic in all three directions. select : str, int or sequence, optional Chemical symbol or atomic number of the species to analyse. Default ``”O”``. ``None`` keeps every atom. A sequence such as ``(“O”, “Na”, “Cl”)`` keeps the listed species and analyses the first, so ions ride along for :func:`pydseams.features.ion_environment`. cutoff : float, optional Neighbour cutoff in Angstroms. Default ``3.5``. bonded : {“auto”, “hbond”, “cutoff”}, optional Graph for rings. ``”auto”`` uses hydrogen bonds when the ``Atoms`` contain H and the selected analysis cloud excludes H; otherwise it uses the cutoff neighbour list. Explicit ``”hbond”`` also requires a heavy-atom selection. Returns Frame

### pydseams.frame.Frame.from_chemfiles {#pydseams-frame-frame-from-chemfiles}

```python
pydseams.frame.Frame.from_chemfiles(filename, frame=None, type_filter=-1, cutoff=None, bonded='cutoff', atom_type=None)
```

Load PDB/GRO/DCD (or any chemfiles format) when chemfiles is linked. Parameters filename : path-like Trajectory chemfiles can read. frame : int, optional 1-indexed frame. Default ``1``. type_filter : int, optional Chemfiles type filter. ``-1`` keeps every type. cutoff : float, optional Neighbour cutoff in Angstroms. Default ``3.5``. bonded : {“auto”, “hbond”, “cutoff”}, optional Graph for rings. Default ``”cutoff”``. atom_type : int or None, optional Species to analyse. ``None`` uses the first particle’s ``c_type``. Returns Frame Raises RuntimeError If chemfiles is not linked in this build of seams-core.

### pydseams.frame.Frame.from_con {#pydseams-frame-frame-from-con}

```python
pydseams.frame.Frame.from_con(filename, frame=None, cutoff=None, bonded='cutoff', atom_type=None)
```

Load an eOn ``.con`` file when readcon-core is linked. Parameters filename : path-like eOn ``.con`` trajectory. frame : int, optional 1-indexed frame. Default ``1``. cutoff : float, optional Neighbour cutoff in Angstroms. Default ``3.5``. bonded : {“auto”, “hbond”, “cutoff”}, optional Graph for rings. Default ``”cutoff”``. atom_type : int or None, optional Species to analyse. ``None`` uses the first particle’s ``c_type``. Returns Frame Raises RuntimeError If readcon-core is not linked in this build of seams-core.

### pydseams.frame.Frame.from_file {#pydseams-frame-frame-from-file}

```python
pydseams.frame.Frame.from_file(filename, frame=None, atom_type=None, cutoff=None, bonded='auto', region=None, all_atoms=False)
```

Load a LAMMPS dump. Parameters filename : path-like LAMMPS dump (``.lammpstrj``, ``.dump``, ``.lammps``). frame : int, optional 1-indexed frame. Default ``1``. atom_type : int or None, optional Species to keep. ``None`` tries type 2 (oxygen) then type 1. cutoff : float, optional Neighbour cutoff in Angstroms. Default ``3.5``. bonded : {“auto”, “hbond”, “cutoff”}, optional Graph for rings. region : ((xlo, ylo, zlo), (xhi, yhi, zhi)) or None, optional Optional rectangular slice. all_atoms : bool, optional Keep every LAMMPS atom type. ``atom_type`` still selects the species used by neighbour and ice analyses. Cannot be combined with ``region``. Returns Frame

### pydseams.frame.Frame.from_xyz {#pydseams-frame-frame-from-xyz}

```python
pydseams.frame.Frame.from_xyz(filename, cutoff=None, bonded='cutoff', atom_type=None, periodic=True)
```

Load an XYZ file through :func:`pydseams.yoda.readXYZ`. Parameters filename : path-like XYZ structure. cutoff : float, optional Neighbour cutoff in Angstroms. Default ``3.5``. bonded : {“auto”, “hbond”, “cutoff”}, optional Graph for rings. Default ``”cutoff”``. atom_type : int or None, optional Species recorded for single-type routines. ``None`` uses the first particle’s ``c_type``. Several elements still enter the neighbour list. periodic : bool, optional When ``False``, the cutoff graph uses free-space distances. The stored box stays the bounding span. Default ``True`` keeps the minimum image inside that span. Returns Frame Raises RuntimeError If this build has no ``readXYZ``.

### pydseams.frame.Frame.guest_occupancy {#pydseams-frame-frame-guest-occupancy}

```python
pydseams.frame.Frame.guest_occupancy(cages, guest_types, radius=4.0)
```

Place guests (methane, THF, ions) in enumerated cages. Parameters cages : sequence of sequences of int Vertex indices into ``self.cloud.pts`` for each cage. guest_types : iterable of int ``c_type`` codes of the guests (LAMMPS types, or atomic numbers for frames built through ASE). radius : float, optional A guest belongs to the nearest cage centroid within this distance in Angstrom. Default ``4.0``, about the radius of a 5^12 cage. Returns pydseams.yoda.GuestOccupancy ``guestsPerCage``, ``cageOfGuest`` (``-1`` when free), ``centreDistance`` and the counts ``occupied``, ``multiply`` and ``free``.

### pydseams.frame.Frame.guest_occupancy_both {#pydseams-frame-frame-guest-occupancy-both}

```python
pydseams.frame.Frame.guest_occupancy_both(signature, guest_types, radius=4.0, max_ring_size=None)
```

Radius and ray-parity occupancy of guests in signature cages. The cages are :meth:`cages_by_signature` on the bonded graph. The radius assignment is :func:`pydseams.yoda.guestOccupancy` (nearest periodic cage centroid within ``radius``); the inside assignment is :func:`pydseams.yoda.guestOccupancyInside` (the guest lies inside the cage’s fan-triangulated faces). Parameters signature : str Cage signature, for example ``”512”``. guest_types : iterable of int ``c_type`` codes of the guests. radius : float, optional Radius of the centroid assignment in Angstrom. Default ``4.0``. max_ring_size : int, optional Largest primitive ring searched. Default ``8``. Returns DualOccupancy ``radius`` and ``inside``, each a :class:`pydseams.yoda.GuestOccupancy` on the same cages, in the order :func:`pydseams.yoda.findBySignature` returns them.

### pydseams.frame.Frame.hbonds {#pydseams-frame-frame-hbonds}

```python
pydseams.frame.Frame.hbonds()
```

Hydrogen-bond neighbour list. Uses :meth:`hbonds_from_donors` when an ASE hydrogen cloud is attached, otherwise :func:`pydseams.yoda.populateHbonds` on the LAMMPS dump. Raises ValueError If no hydrogens are available (arrays-only frame, or ``bonded=’hbond’`` without H).

### pydseams.frame.Frame.hbonds_from_donors {#pydseams-frame-frame-hbonds-from-donors}

```python
pydseams.frame.Frame.hbonds_from_donors(donor_hs, h_cloud=None, dist=2.42, angle=30.0)
```

Hydrogen-bond network from an explicit donor-H index list. Parameters donor_hs : sequence of int ``hCloud`` indices. Each H is paired with the heavy atom that shares its ``molID``. h_cloud : PointCloudDouble, optional Hydrogen cloud. Defaults to the ASE hydrogen cloud. dist, angle : float, optional Acceptor-H distance and acceptor-centered O-O-H angle.

### pydseams.frame.Frame.hydration_shell_rings {#pydseams-frame-frame-hydration-shell-rings}

```python
pydseams.frame.Frame.hydration_shell_rings(ion_types, k=4, ring_adjacent=True, cutoff=None)
```

Rings of the water network through each ion’s first shell. Parameters ion_types : iterable of int ``c_type`` codes of the ions, as in :meth:`ion_environment`. k, ring_adjacent, cutoff As in :meth:`ion_environment`. Returns tuple ``(env, census)``: the :class:`pydseams.yoda.IonEnvironment` and a list with one row per ion, ``census\[i\]\[s\]`` counting the primitive rings of size ``s`` (up to six, the rings :attr:`rings` keeps) with a vertex in ion ``i``’s shell. An ion is not a vertex of the network, so the rings it would have closed are gone; the shell census measures how far the network survives around it.

### pydseams.frame.Frame.ion_cloud {#pydseams-frame-frame-ion-cloud}

```python
pydseams.frame.Frame.ion_cloud(table)
```

One COM vertex per ion molecule. Parameters table : pydseams.yoda.SiteTable Type-to-kind map. Cation molecules restamp to type 1, anions to type 2.

### pydseams.frame.Frame.ion_environment {#pydseams-frame-frame-ion-environment}

```python
pydseams.frame.Frame.ion_environment(ion_types, k=4, ring_adjacent=True, cutoff=None)
```

Class every ion by its first water shell against the seeded cages. Parameters ion_types : iterable of int ``c_type`` codes of the ions (LAMMPS types, or atomic numbers for frames built through ASE). k : int, optional Neighbours kept in the seeded k-nearest graphs. Default ``4``. ring_adjacent : bool, optional Passed to :meth:`seeded_affiliation`. Default ``True``. cutoff : float, optional First-shell radius in Angstrom. Default ``self.cutoff``. Returns pydseams.yoda.IonEnvironment ``ion`` (cloud indices), ``shell``, ``iceFraction``, ``state`` (:class:`pydseams.yoda.IonState`) and the counts ``nIce``, ``nFront``, ``nLiquid``.

### pydseams.frame.Frame.load_frame {#pydseams-frame-frame-load-frame}

```python
pydseams.frame.Frame.load_frame(frame)
```

Reload a later (or earlier) frame from the same trajectory file. Parameters frame : int 1-indexed frame to read. Raises ValueError If this ``Frame`` was built from arrays or ASE and has no trajectory path. Notes Clears cached neighbours, hydrogen bonds, rings, and cages.

### pydseams.frame.Frame.local_structure_index {#pydseams-frame-frame-local-structure-index}

```python
pydseams.frame.Frame.local_structure_index(cutoff=3.7)
```

Per-site variance of successive neighbour-distance gaps, in Å². The last gap reaches the closest site outside the shell. Sites with no interior or no exterior neighbour return NaN. Compute on water oxygens for the Shiratani–Sasai LSI. Requires the local-order extra.

### pydseams.frame.Frame.monolayer_rings {#pydseams-frame-frame-monolayer-rings}

```python
pydseams.frame.Frame.monolayer_rings(output_dir, sheet_area, max_depth=4)
```

Classify quasi-2D polygon rings and read coverage back. Parameters output_dir : path-like Directory for :func:`pydseams.yoda.polygonRingAnalysis`. sheet_area : float Sheet area passed to the engine. max_depth : int, optional Largest ring size. Default ``4``. Returns dict ``{ring_size: {“count”: int, “coverage_xy”: float}}`` parsed from ``topoMonolayer/coverageAreaXY.dat``.

### pydseams.frame.Frame.n_atoms {#pydseams-frame-frame-n-atoms}

```python
pydseams.frame.Frame.n_atoms()
```

Number of particles in the analysed cloud (``cloud.nop``).

### pydseams.frame.Frame.neighbor_list {#pydseams-frame-frame-neighbor-list}

```python
pydseams.frame.Frame.neighbor_list()
```

Cutoff neighbour list. One atom type uses :func:`pydseams.yoda.neighListO`. Several types use the union of like and unlike lists from :func:`pydseams.yoda.neighList`. Rows are atom IDs within :attr:`cutoff`. :attr:`periodic` ``False`` widens the box for that call so a molecule is not wrapped across its bounding span, then restores the span.

### pydseams.frame.Frame.pairs {#pydseams-frame-frame-pairs}

```python
pydseams.frame.Frame.pairs(table)
```

Mutual nearest cation-anion pairs in an ion COM cloud. Pair indices refer to the ion cloud returned by :meth:`ion_cloud`.

### pydseams.frame.Frame.positions {#pydseams-frame-frame-positions}

```python
pydseams.frame.Frame.positions()
```

List of ``(x, y, z)`` coordinates for the analysed particles.

### pydseams.frame.Frame.predict_class {#pydseams-frame-frame-predict-class}

```python
pydseams.frame.Frame.predict_class(x)
```

Predict a class for one feature row. Parameters x : sequence of float Feature vector matching the last :meth:`fit_classifier` fit. Returns int Predicted class index. Raises RuntimeError If :meth:`fit_classifier` has not been called.

### pydseams.frame.Frame.prism_counts {#pydseams-frame-frame-prism-counts}

```python
pydseams.frame.Frame.prism_counts(max_depth=6, shape_matching=False)
```

Return (perfect, imperfect) prism counts by basal ring size. Uses the frame’s bonded graph and enumerates primitive rings through max_depth. Shape matching enables recognition of deformed prisms. Results stay in memory.

### pydseams.frame.Frame.rdf {#pydseams-frame-frame-rdf}

```python
pydseams.frame.Frame.rdf(type_i, type_j, cutoff=12.0, binwidth=0.05)
```

Partial 3D radial distribution function. Parameters type_i : int First species type code (``c_type``). type_j : int Second species type code (``c_type``). cutoff : float, optional RDF cutoff in Angstroms. Default ``12.0``. binwidth : float, optional Histogram width. Default ``0.05``. Returns r, g : list of float Bin centres and ``g_IJ(r)`` from :func:`pydseams.yoda.partialRdf`.

### pydseams.frame.Frame.rdf_2d {#pydseams-frame-frame-rdf-2d}

```python
pydseams.frame.Frame.rdf_2d(output_dir, cutoff=12.0, binwidth=0.05)
```

2D radial distribution function for identical atom types. Parameters output_dir : path-like Directory for :func:`pydseams.yoda.rdf2Danalysis_AA`. cutoff : float, optional RDF cutoff in Angstroms. Default ``12.0``. binwidth : float, optional Histogram width. Default ``0.05``. Returns r, g : list of float Bin centres and ``g(r)`` parsed from ``topoMonolayer/rdf.dat``.

### pydseams.frame.Frame.rings {#pydseams-frame-frame-rings}

```python
pydseams.frame.Frame.rings()
```

Primitive rings up to size 6 on :attr:`bonds_by_index`. Computed by :class:`pydseams.yoda.RingUpdater` and cached.

### pydseams.frame.Frame.rings_recomputed_sources {#pydseams-frame-frame-rings-recomputed-sources}

```python
pydseams.frame.Frame.rings_recomputed_sources()
```

Sources recomputed by the last :class:`~pydseams.yoda.RingUpdater` pass.

### pydseams.frame.Frame.running_cn {#pydseams-frame-frame-running-cn}

```python
pydseams.frame.Frame.running_cn(type_i, type_j, cutoff=12.0, binwidth=0.05)
```

Running site-site coordination number. Parameters type_i : int First species type code (``c_type``). type_j : int Second species type code (``c_type``). cutoff : float, optional Histogram cutoff in Angstroms. Default ``12.0``. binwidth : float, optional Histogram width. Default ``0.05``. Returns list of float Running CN at each bin outer edge from :func:`pydseams.yoda.partialRdfHist` and :func:`pydseams.yoda.runningCN` with ``rho_J = nJ / dumpVolume``.

### pydseams.frame.Frame.seeded_affiliation {#pydseams-frame-frame-seeded-affiliation}

```python
pydseams.frame.Frame.seeded_affiliation(k=4, candidate_cutoff=None, ring_adjacent=False)
```

Seeded (hysteresis) per-atom cage flags. Strict-graph seeds from a mutual k-nearest list, permissive completion on the union k-nearest list. Delegates to :func:`pydseams.yoda.seededCageAffiliation`. Parameters k : int, optional Neighbours kept in each k-nearest graph. Default ``4``. candidate_cutoff : float or None, optional Candidate-list cutoff. ``None`` uses :attr:`cutoff` ``+ 1.5``. ring_adjacent : bool, optional Fill the last vertex of any union-graph six-ring whose other vertices carry a label, repeated to a fixed point. A frame with no accepted ring stays empty. Default ``False``. Returns CageScore

### pydseams.frame.Frame.soap {#pydseams-frame-frame-soap}

```python
pydseams.frame.Frame.soap(iatom=None, n_max=3, l_max=6, rcut=None)
```

SOAP power spectrum of one particle, or of every particle. Parameters iatom : int or None, optional Particle index. ``None`` (default) computes every particle via :func:`pydseams.yoda.soapSpectrumAll`. n_max : int, optional Radial basis size. Default ``3``. l_max : int, optional Angular momentum cutoff. Default ``6``. rcut : float or None, optional SOAP cutoff. ``None`` uses :attr:`cutoff`. Returns list of float or list of list of float One spectrum, or one spectrum per particle.

### pydseams.frame.Frame.steinhardt {#pydseams-frame-frame-steinhardt}

```python
pydseams.frame.Frame.steinhardt(order_l=6)
```

Local and neighbour-averaged Steinhardt parameters. Parameters order_l : int, optional Degree ``l`` (3, 4, or 6). Default ``6``. Returns dict ``ql`` and ``ql_bar`` lists from :func:`pydseams.yoda.steinhardtQl`.

### pydseams.frame.Frame.steinhardt_voronoi {#pydseams-frame-frame-steinhardt-voronoi}

```python
pydseams.frame.Frame.steinhardt_voronoi(order_l=6, cutoff=None)
```

Voronoi facet-area weighted Steinhardt parameters. Parameters order_l : int, optional Degree ``l``. Default ``6``. cutoff : float or None, optional Candidate cutoff for the Voronoi pass. ``None`` uses :attr:`cutoff`. Returns dict ``ql`` and ``ql_bar`` lists from :func:`pydseams.yoda.steinhardtQlVoronoi`.

### pydseams.frame.Frame.tetrahedral_q {#pydseams-frame-frame-tetrahedral-q}

```python
pydseams.frame.Frame.tetrahedral_q(cutoff=3.5)
```

Per-site angular order of the four closest selected particles. Returns NaN for sites with fewer than four neighbours inside the cutoff (angstroms) or a coincident neighbour. Pass None to use the four closest sites at any distance. Requires the local-order extra.

### pydseams.frame.Frame.to_ase {#pydseams-frame-frame-to-ase}

```python
pydseams.frame.Frame.to_ase()
```

ASE ``Atoms`` for this frame. Returns ase.Atoms Orthorhombic cell, ``pbc=True``. After :meth:`chill_plus` (or :meth:`chill`), ``arrays\[‘ice_type’\]`` holds the labels. After :meth:`cages`, ``arrays\[‘hc’\]`` and ``arrays\[‘ddc’\]`` hold the cage flags.

### pydseams.frame.Frame.to_solvis {#pydseams-frame-frame-to-solvis}

```python
pydseams.frame.Frame.to_solvis(expand_box=True)
```

``solvis.System`` for this frame. Parameters expand_box : bool, optional Passed to ``solvis.system.System``. Default ``True``. Returns solvis.system.System Notes Optional extra: ``pip install ‘pydseamslib\[solvis\]’``.

### pydseams.frame.Frame.topology_library {#pydseams-frame-frame-topology-library}

```python
pydseams.frame.Frame.topology_library(label, hops=2, max_ring_size=7, colour_types=False, library=None)
```

Add this frame’s local keys to a key library under ``label``. Returns the :class:`pydseams.yoda.KeyLibrary` (a new one unless ``library`` is given). ``yoda.writeLibrary`` turns it into text.

### pydseams.frame.Frame.voronoi_features {#pydseams-frame-frame-voronoi-features}

```python
pydseams.frame.Frame.voronoi_features(cutoff=None)
```

Per-atom ``\[q4, q6, q8\]`` from one Voronoi-weighted pass. Parameters cutoff : float or None, optional Candidate cutoff. ``None`` uses :attr:`cutoff`. Returns list of list of float One ``\[q4, q6, q8\]`` row per particle from :func:`pydseams.yoda.voronoiFeatures`.

### pydseams.frame.IceCounts {#pydseams-frame-icecounts}

```python
pydseams.frame.IceCounts
```

CHILL / CHILL+ histogram of ice labels on one frame. Keys are the ``AtomStateType`` names written by the classifier (``cubic``, ``hexagonal``, ``water``, ``interfacial``, ``clathrate``, ``interClathrate``, ``unclassified``, ``reCubic``, ``reHex``). Missing keys read as ``0`` via attribute access, so ``counts.cubic`` and ``counts\[‘cubic’\]`` are equivalent. Notes ``repr`` omits zero-count labels.

### pydseams.frame.read {#pydseams-frame-read}

```python
pydseams.frame.read(filename, **kwargs)
```

Load a trajectory. Format follows the suffix. Thin wrapper around :func:`pydseams.io.read`. Parameters filename : path-like Trajectory or structure file. **kwargs Forwarded to :func:`pydseams.io.read` (``frame``, ``cutoff``, ``bonded``, ``atom_type``, ``region``, …). Returns Frame

### pydseams.io.available_readers {#pydseams-io-available-readers}

```python
pydseams.io.available_readers()
```

Report which compiled I/O readers this build linked. Returns dict of str to bool ``lammps`` is always ``True``. ``xyz``, ``chemfiles``, and ``readcon`` are ``True`` when the matching :mod:`pydseams.yoda` symbol exists (``readXYZ``, ``readChemfiles``, ``readCon``).

### pydseams.io.has_chemfiles {#pydseams-io-has-chemfiles}

```python
pydseams.io.has_chemfiles()
```

True when this build can read PDB, GRO, and DCD through chemfiles.

### pydseams.io.has_readcon {#pydseams-io-has-readcon}

```python
pydseams.io.has_readcon()
```

True when this build can read an eOn ``.con`` file through readcon-core.

### pydseams.io.iter_trajectory {#pydseams-io-iter-trajectory}

```python
pydseams.io.iter_trajectory(path)
```

Yield frames selected by atom and residue names through readcon. Coordinate-only formats such as XTC and TRR take a labelled topology with the same atom count and order. Indices follow Python slicing. Each selected residue contributes one centre and its named hydrogens. Coordinates are in angstroms. Frames retain the full periodic cell. Install the readcon extra to use this reader.

### pydseams.io.read {#pydseams-io-read}

```python
pydseams.io.read(path, frame=1, **kwargs)
```

Load one configuration. Format follows the file suffix. Parameters path : path-like Trajectory or structure file. frame : int, optional 1-indexed frame for multi-frame formats (LAMMPS, chemfiles, ``.con``). Ignored for XYZ. Default ``1``. **kwargs Forwarded to the matching :class:`~pydseams.frame.Frame` constructor (``cutoff``, ``bonded``, ``atom_type``, ``region``, ``all_atoms``, …). ``all_atoms=True`` retains every LAMMPS type for mixed-site analyses such as ionic pairs and domains. Returns Frame Analysable configuration. Notes Suffix dispatch: * ``.xyz`` – :meth:`Frame.from_xyz` * ``.con`` – :meth:`Frame.from_con` * ``.pdb``, ``.gro``, ``.dcd`` – :meth:`Frame.from_chemfiles` * otherwise – :meth:`Frame.from_file` (LAMMPS dump) ``.lammpstrj``, ``.dump``, and ``.lammps`` take the LAMMPS path. Builds without chemfiles or readcon raise ``RuntimeError`` when those suffixes are used.

### pydseams.liquid_order.local_structure_index {#pydseams-liquid-order-local-structure-index}

```python
pydseams.liquid_order.local_structure_index(frame, cutoff)
```

### pydseams.liquid_order.tetrahedral_q {#pydseams-liquid-order-tetrahedral-q}

```python
pydseams.liquid_order.tetrahedral_q(frame, cutoff)
```

### pydseams.solvis.to_solvis {#pydseams-solvis-to-solvis}

```python
pydseams.solvis.to_solvis(frame, expand_box=True)
```

Wrap a :class:`~pydseams.frame.Frame` as a ``solvis.System``. Parameters frame : Frame Configuration to view. Converted through :meth:`~pydseams.frame.Frame.to_ase`. expand_box : bool, optional Passed to ``solvis.system.System``. Default ``True``. Returns solvis.system.System Raises ImportError If solvis is not installed (``pip install ‘pydseamslib\[solvis\]’``).

### seams.errors.SeamsError {#seams-errors-seamserror}

```python
seams.errors.SeamsError
```

A failing ``libseams`` call. Parameters status : int ``seams_status`` value. Zero is never stored here. message : str Text from ``seams_handle_error``. The library owns that string. This object keeps a Python copy. Examples >>> err = SeamsError(5, “no frame”) >>> err.status 5 >>> “no frame” in str(err) True

### seams.results.BrineFronts {#seams-results-brinefronts}

```python
seams.results.BrineFronts
```

The two planar crossings of an ice slab. Notes z_enter, z_exit, thickness : float Angstrom, in the frame cell. found : bool True when both crossings exist.

### seams.results.BrinePockets {#seams-results-brinepockets}

```python
seams.results.BrinePockets
```

Brine pockets trapped in the ice. Notes n_pockets, reservoir_waters : int Pocket count, and waters in the liquid reservoir.

### seams.results.BrineProfile {#seams-results-brineprofile}

```python
seams.results.BrineProfile
```

Concentration profile along the front axis. Notes z_mid : numpy.ndarray Bin centers, angstrom. molality : numpy.ndarray mol/kg. ice_fraction : numpy.ndarray Dimensionless.

### seams.results.BrineRejection {#seams-results-brinerejection}

```python
seams.results.BrineRejection
```

Molal rejection coefficient ``k = c_ice / c_liquid``. Notes c_ice, c_liquid, k : float Molality is formula units per kilogram of water. defined : bool True when the ratio is defined.

### seams.results.CageCounts {#seams-results-cagecounts}

```python
seams.results.CageCounts
```

Disjoint cage populations. The four counts sum to the kept atoms. Notes n_hex, n_cubic, n_water, n_both : int Hexagonal only, double-diamond only, neither, both.

### seams.results.ChillResult {#seams-results-chillresult}

```python
seams.results.ChillResult
```

Per-atom labels and the nine population counts. Notes labels : numpy.ndarray One ``seams_label`` code per atom. counts : numpy.ndarray Length 9, in label order. ``counts\[i\]`` is the population of label ``i``.

### seams.results.Contacts {#seams-results-contacts}

```python
seams.results.Contacts
```

Mutual nearest unlike contacts. Notes n_pairs, n_cation, n_anion : int Pair count, and molecule counts in the ion cloud.

### seams.results.Density {#seams-results-density}

```python
seams.results.Density
```

Number density along one cell axis. Notes z, rho : numpy.ndarray ``z`` is angstrom. ``rho`` is number density.

### seams.results.Domains {#seams-results-domains}

```python
seams.results.Domains
```

Largest polar or apolar domain. Notes n_subset, largest : int Atoms in the subset, and the largest domain. percolation : float Dimensionless percolation score.

### seams.results.NeighbourCSR {#seams-results-neighbourcsr}

```python
seams.results.NeighbourCSR
```

Caller-owned CSR. Column indices omit the self entry. Notes offsets : numpy.ndarray Length ``n + 1``. columns : numpy.ndarray Neighbour indices, length ``nnz``.

### seams.results.OrderResult {#seams-results-orderresult}

```python
seams.results.OrderResult
```

Local Steinhardt parameters. Notes ql : numpy.ndarray One value per atom. Dimensionless. qlbar : numpy.ndarray or None Lechner-Dellago neighbour average, or ``None`` when it was not requested.

### seams.results.Prisms {#seams-results-prisms}

```python
seams.results.Prisms
```

Perfect and imperfect prism counts. Index is the ring size. Notes perfect, imperfect : numpy.ndarray Length 7. Sizes 3 through 6 are filled.

### seams.results.Q4Q6 {#seams-results-q4q6}

```python
seams.results.Q4Q6
```

Local q4 and q6 on the cutoff shell. Notes q4, q6 : numpy.ndarray One dimensionless value per atom.

### seams.results.Q6Bar {#seams-results-q6bar}

```python
seams.results.Q6Bar
```

Lechner-Dellago q6-bar and the summary above a cut. Notes values : numpy.ndarray or None Per-atom q6-bar when it was requested. mean, median : float Over finite values. An empty finite set is 0. n_ice : int Count of finite values strictly above the cut.

### seams.results.Radial {#seams-results-radial}

```python
seams.results.Radial
```

Partial radial distribution. Notes r : numpy.ndarray Bin centers, angstrom. g : numpy.ndarray Dimensionless g(r). count : numpy.ndarray Pair counts. volume : float Cell volume, cubic angstrom.

### seams.results.RingLabels {#seams-results-ringlabels}

```python
seams.results.RingLabels
```

Size and symmetry of each primitive ring. Notes sizes, symmetry : numpy.ndarray ``symmetry`` holds ``seams_ring_symmetry`` codes.

### seams.results.RunningCN {#seams-results-runningcn}

```python
seams.results.RunningCN
```

Running coordination number. Notes r, cn : numpy.ndarray ``r`` is angstrom. ``cn`` is dimensionless.

### seams.results.Signature {#seams-results-signature}

```python
seams.results.Signature
```

Closed cages of one ring census. Notes n_cages, n_atoms : int Closed cages, and the distinct vertices in those cages.

### seams.system.Brine {#seams-system-brine}

```python
seams.system.Brine
```

One planar ice/brine front. Lengths are angstrom in the frame cell. Molality is formula units per kilogram of water.

### seams.system.Brine.fronts {#seams-system-brine-fronts}

```python
seams.system.Brine.fronts()
```

The two planar crossings of the ice slab. Returns BrineFronts

### seams.system.Brine.pockets {#seams-system-brine-pockets}

```python
seams.system.Brine.pockets()
```

Brine pockets, excluding the liquid reservoir. Returns BrinePockets

### seams.system.Brine.profile {#seams-system-brine-profile}

```python
seams.system.Brine.profile()
```

Concentration profile along the front axis. Returns BrineProfile

### seams.system.Brine.rejection {#seams-system-brine-rejection}

```python
seams.system.Brine.rejection()
```

Molal rejection coefficient. Returns BrineRejection

### seams.system.CageSet {#seams-system-cageset}

```python
seams.system.CageSet
```

Hexagonal and double-diamond membership. The four population counts are disjoint.

### seams.system.CageSet.counts {#seams-system-cageset-counts}

```python
seams.system.CageSet.counts()
```

Disjoint population counts. Returns CageCounts

### seams.system.CageSet.labels {#seams-system-cageset-labels}

```python
seams.system.CageSet.labels(n)
```

One cage code per atom. Parameters n : int Frame count stored when the set was built. Atoms the species filter dropped are ``CAGE_WATER``. Returns numpy.ndarray ``seams_cage_label`` codes.

### seams.system.Fingerprint {#seams-system-fingerprint}

```python
seams.system.Fingerprint
```

Topological fingerprint of one bonded frame.

### seams.system.Fingerprint.atom {#seams-system-fingerprint-atom}

```python
seams.system.Fingerprint.atom(index)
```

Local key of one atom. Parameters index : int Atom index, from 0. Returns str

### seams.system.Fingerprint.key {#seams-system-fingerprint-key}

```python
seams.system.Fingerprint.key()
```

Frame key text. Returns str

### seams.system.Fingerprint.rings {#seams-system-fingerprint-rings}

```python
seams.system.Fingerprint.rings(length=13)
```

Primitive-ring census stored on the fingerprint. Parameters length : int ``counts\[k\]`` is the number of rings of size ``k``. Returns numpy.ndarray

### seams.system.Neighbours {#seams-system-neighbours}

```python
seams.system.Neighbours
```

A neighbour list copied out of a frame. The list does not borrow the :class:`System`. Closing the system leaves this object valid. :meth:`copy` returns NumPy arrays the caller owns. The deprecated ``seams_neighbours_csr`` borrow is not exposed. Examples Use it as a context manager so the handle is freed on the way out.

### seams.system.Neighbours.copy {#seams-system-neighbours-copy}

```python
seams.system.Neighbours.copy()
```

Copy the CSR into new arrays. Returns NeighbourCSR Offsets of length ``n + 1`` and columns that omit self.

### seams.system.Neighbours.count {#seams-system-neighbours-count}

```python
seams.system.Neighbours.count()
```

Number of neighbour rows. Returns int

### seams.system.Rings {#seams-system-rings}

```python
seams.system.Rings
```

Primitive rings (Franzblau) for one frame. The caller owns the handle. Copied arrays do not borrow it.

### seams.system.Rings.count {#seams-system-rings-count}

```python
seams.system.Rings.count()
```

Number of primitive rings. Returns int

### seams.system.Rings.fraction {#seams-system-rings-fraction}

```python
seams.system.Rings.fraction(size=0)
```

Fraction of atoms that sit on a ring. Parameters size : int Ring size. Zero counts atoms that sit on any ring. Returns float A value in ``\[0, 1\]``.

### seams.system.Rings.histogram {#seams-system-rings-histogram}

```python
seams.system.Rings.histogram(length=13)
```

Count rings by size. Parameters length : int ``counts\[k\]`` is the number of rings of size ``k``. Returns numpy.ndarray

### seams.system.Rings.labels {#seams-system-rings-labels}

```python
seams.system.Rings.labels()
```

Size and symmetry of every ring. Returns RingLabels

### seams.system.Rings.vertices {#seams-system-rings-vertices}

```python
seams.system.Rings.vertices(index=0)
```

Atom indices of one ring. Parameters index : int Ring index, from 0. Returns numpy.ndarray Positions in the frame.

### seams.system.System {#seams-system-system}

```python
seams.system.System
```

One frame. The constructor calls ``seams_system_create``. :meth:`close` calls ``seams_system_destroy``. Neighbour, ring, cage, fingerprint, and brine handles keep their own copies, so they stay valid after this frame is closed. A failing call raises :class:`seams.errors.SeamsError`. ``status`` is the ``seams_status`` code and ``message`` is the string from ``seams_handle_error`` on this handle. Parameters path : str, optional When set, the frame is read before the constructor returns. frame : int 1-based frame index. Used only with ``path``. type_filter : int Negative keeps every atom. Zero guesses oxygen. A positive value keeps that species. Used only with ``path``. Examples >>> err = SeamsError(5, “no frame”) >>> err.status == 5 True

### seams.system.System.brine {#seams-system-system-brine}

```python
seams.system.System.brine(ice)
```

Planar ice/brine front from an ice assignment. Parameters ice : array_like Per-atom ice flags, length equal to the frame count. Non-zero means ice. water_type : int Species kept as water, or 0 to keep every atom that is not ``ion_type``. ion_type : int Species kept as ions, or 0 when the frame has no ions. n_bins : int Bin count along ``axis``. Zero selects 20. axis : int 0, 1, or 2. Returns Brine

### seams.system.System.cages {#seams-system-system-cages}

```python
seams.system.System.cages(cutoff=3.5)
```

Hexagonal and double-diamond cage membership. Parameters cutoff : float Positive cutoff, angstrom. The seeded graph uses ``cutoff + 1.5`` as the candidate shell. type : int Species filter. graph : int ``CUTOFF``, ``MUTUAL``, or ``SEEDED``. k : int Neighbour count for the seeded graph. Zero selects 4. complete : bool Run ring-adjacent completion on the seeded labels. Returns CageSet

### seams.system.System.chill {#seams-system-system-chill}

```python
seams.system.System.chill(cutoff=3.5)
```

Original CHILL labels. Parameters cutoff : float Positive cutoff, angstrom. type : int Species filter. Zero and negative keep every atom. graph : int ``CUTOFF`` or ``MUTUAL``. rule : int ``PUBLISHED`` (inclusive edges) or ``OVITO`` (open edges). Returns ChillResult Labels and nine population counts.

### seams.system.System.chill_plus {#seams-system-system-chill-plus}

```python
seams.system.System.chill_plus(cutoff=3.5)
```

CHILL+ labels. Parameters cutoff : float Positive cutoff, angstrom. type : int Species filter. graph : int ``CUTOFF`` or ``MUTUAL``. rule : int ``PUBLISHED`` is the Nguyen and Molinero table. ``OVITO`` is the Chill+ modifier. Returns ChillResult

### seams.system.System.contact_pairs {#seams-system-system-contact-pairs}

```python
seams.system.System.contact_pairs(spec)
```

Mutual nearest unlike contacts of a site table. Parameters spec : str Text such as ``”1=cationHead,2=anion”``. Returns Contacts

### seams.system.System.coordination {#seams-system-system-coordination}

```python
seams.system.System.coordination(rmax=3.5, bins=8)
```

Coordination number integrated to ``rmax``. Parameters rmax : float Cutoff, angstrom. bins : int Bin count used by the integral. type_i, type_j : int Species pair. Returns float

### seams.system.System.count {#seams-system-system-count}

```python
seams.system.System.count()
```

Number of points in the frame. Returns int Zero when the handle holds no frame.

### seams.system.System.density {#seams-system-system-density}

```python
seams.system.System.density()
```

Number density along one cell axis. Parameters type : int Species filter. Zero keeps every atom. bins : int Bin count. Zero selects about one bin per 0.1 angstrom. axis : int 0, 1, or 2. Returns Density

### seams.system.System.domains {#seams-system-system-domains}

```python
seams.system.System.domains(spec, cutoff=3.5)
```

Largest domain of a polar or apolar site subset. Parameters spec : str Site-table text. cutoff : float Positive cutoff, angstrom. subset : int 0 is apolar, 1 is polar. Returns Domains

### seams.system.System.fingerprint {#seams-system-system-fingerprint}

```python
seams.system.System.fingerprint(cutoff=3.5)
```

Topological fingerprint of the bonded graph. Parameters cutoff : float Positive cutoff, angstrom. type : int Species filter. k : int Neighbour count. graph : int ``CUTOFF`` or ``MUTUAL``. hops : int Neighbourhood depth. Zero selects 2. color_types : bool Color vertices by LAMMPS type. Returns Fingerprint

### seams.system.System.hbonds {#seams-system-system-hbonds}

```python
seams.system.System.hbonds(cutoff=3.5)
```

Hydrogen-bond count. Parameters cutoff : float Heavy-atom neighbour cutoff, angstrom. oxygen_type : int Heavy-atom species. Zero keeps every atom that is not ``hydrogen_type``. hydrogen_type : int Positive hydrogen species. dist : float Acceptor-H distance, angstrom. angle : float O-O-H angle, degrees. Returns int

### seams.system.System.knn {#seams-system-system-knn}

```python
seams.system.System.knn(k=4, cutoff=3.5)
```

Build a k-nearest neighbour list. Parameters k : int Positive neighbour count. cutoff : float Nominal cutoff, angstrom. Candidates lie inside ``cutoff + 1.5``. type : int Species filter. Zero keeps every atom. mutual : bool True keeps the mutual graph. Returns Neighbours

### seams.system.System.monolayer {#seams-system-system-monolayer}

```python
seams.system.System.monolayer(cutoff=3.5)
```

Primitive-ring histogram of a monolayer. Parameters cutoff : float Positive cutoff, angstrom. type : int Species filter. max_depth : int Largest ring size. sheet_area : float Sheet area, square angstrom. Must be positive. Returns numpy.ndarray ``counts\[k\]`` is the number of rings of size ``k``.

### seams.system.System.neighbours {#seams-system-system-neighbours}

```python
seams.system.System.neighbours(cutoff=3.5)
```

Build a cutoff or mutual neighbour list. Parameters cutoff : float Positive cutoff, angstrom. The mutual graph keeps the four nearest neighbours inside ``cutoff + 1.5``. type : int Species to keep. Zero and negative keep every atom. graph : int ``CUTOFF`` or ``MUTUAL``. ``SEEDED`` is rejected. Returns Neighbours

### seams.system.System.prisms {#seams-system-system-prisms}

```python
seams.system.System.prisms(cutoff=3.5)
```

Perfect and imperfect prism counts by ring size. Parameters cutoff : float Positive cutoff, angstrom. type : int Species filter. Returns Prisms

### seams.system.System.q4q6 {#seams-system-system-q4q6}

```python
seams.system.System.q4q6(cutoff=3.5)
```

Local q4 and q6 on the cutoff shell. Parameters cutoff : float Positive cutoff, angstrom. type : int Species filter. Returns Q4Q6

### seams.system.System.q6bar {#seams-system-system-q6bar}

```python
seams.system.System.q6bar(cutoff=3.5)
```

Lechner-Dellago q6-bar and the ice count above a cut. Parameters cutoff : float Positive cutoff, angstrom. type : int Species filter. ice_cut : float Finite threshold. The nucleation census uses 0.36. values : bool When true, keep the per-atom array. Returns Q6Bar

### seams.system.System.ql {#seams-system-system-ql}

```python
seams.system.System.ql(order, cutoff=3.5)
```

Local ``q_l``, and the Lechner-Dellago average on request. Parameters order : int Spherical-harmonic degree. 3, 4, 6, 8, or 12. cutoff : float Positive cutoff, angstrom. type : int Species filter. graph : int ``CUTOFF`` or ``MUTUAL``. average : bool When true, also return the neighbour mean of ``q_lm``. Returns OrderResult Notes The scratch belongs to the calling thread, so separate frames may run this in separate threads.

### seams.system.System.rdf {#seams-system-system-rdf}

```python
seams.system.System.rdf(rmax=6.0, bins=0)
```

Partial radial distribution. Parameters rmax : float Maximum radius, angstrom. bins : int Bin count. Zero selects about ``rmax / 0.1``. type_i, type_j : int Species pair. Zero keeps every atom. Returns Radial

### seams.system.System.read {#seams-system-system-read}

```python
seams.system.System.read(path)
```

Read one frame from a trajectory. Parameters path : path-like CON, LAMMPS dump, or another format the linked readers accept. frame : int 1-based index, in file order. type_filter : int Negative keeps every atom. Zero guesses oxygen. A positive value keeps that species. Returns int Atom count after the read. Raises seams.SeamsError On a rejected path, index, or reader failure. The previous frame is kept.

### seams.system.System.rings {#seams-system-system-rings}

```python
seams.system.System.rings(cutoff=3.5)
```

Find primitive rings up to a maximum size. Parameters cutoff : float Positive cutoff, angstrom. type : int Species filter. graph : int ``CUTOFF`` or ``MUTUAL``. max_depth : int Largest ring size, from 3 through 12. Returns Rings

### seams.system.System.running_cn {#seams-system-system-running-cn}

```python
seams.system.System.running_cn(rmax=6.0, bins=0)
```

Running coordination number along the partial RDF. Parameters rmax : float Maximum radius, angstrom. bins : int Bin count. Zero selects the library default. type_i, type_j : int Species pair. Returns RunningCN

### seams.system.System.set_frame {#seams-system-system-set-frame}

```python
seams.system.System.set_frame(positions, box)
```

Store points and a cell. Parameters positions : array_like Shape ``(n, 3)``, angstrom. The library copies them. box : array_like Three lengths or nine LAMMPS numbers. types, ids : array_like, optional Length ``n``. Omitted types are 1. Omitted ids are ``1..n``. Ids must be unique. Raises ValueError If the arrays have the wrong shape. seams.SeamsError If the cell has no volume or an id is repeated.

### seams.system.System.set_ortho {#seams-system-system-set-ortho}

```python
seams.system.System.set_ortho(positions, lengths)
```

Store points in an orthorhombic cell. Parameters positions : array_like Shape ``(n, 3)``, angstrom. lengths : tuple of float Edge lengths ``lx, ly, lz``, angstrom. origin : tuple of float Origin ``ox, oy, oz``, angstrom. The tilt is zero.

### seams.system.System.signature {#seams-system-system-signature}

```python
seams.system.System.signature(spec='512', cutoff=3.5)
```

Count closed cages whose ring census matches ``spec``. Parameters spec : str A census such as ``”512”``, ``”51262”``, or the names ``”hc”`` and ``”ddc”``. cutoff : float Positive cutoff, angstrom. type : int Species filter. graph : int ``CUTOFF`` or ``MUTUAL``. Returns Signature

### seams.system.System.site_table {#seams-system-system-site-table}

```python
seams.system.System.site_table(spec)
```

Per-atom site kind. Parameters spec : str Site-table text such as ``”1=waterO”``. Returns numpy.ndarray One ``seams_site_kind`` per atom.

### seams.system.System.steinhardt {#seams-system-system-steinhardt}

```python
seams.system.System.steinhardt(order, cutoff=3.5)
```

Local ``q_l`` on the cutoff shell. Parameters order : int 3, 4, 6, 8, or 12. cutoff : float Positive cutoff, angstrom. type : int Species filter. average : bool When true, include the Lechner-Dellago average. Returns OrderResult

### seams.system.System.symmetry {#seams-system-system-symmetry}

```python
seams.system.System.symmetry(shell_rmsd=0.0)
```

Per-atom ice shell label. Parameters shell_rmsd : float Positive scale-free residual ceiling. A non-positive value selects the library ceiling. Returns numpy.ndarray 0 ice Ih, 1 ice Ic, 2 stacking disorder, or -1.

### seams.system.abi {#seams-system-abi}

```python
seams.system.abi()
```

ABI number compiled into this library. Returns int ``SEAMS_ABI``, which is 1.

### seams.system.as_box9 {#seams-system-as-box9}

```python
seams.system.as_box9(box)
```

Nine LAMMPS cell numbers. Parameters box : array_like Three edge lengths, or ``ox, oy, oz, lx, ly, lz, xy, xz, yz``. Lengths are angstrom. Three lengths place the origin at zero and set the tilt to zero. Returns numpy.ndarray Shape ``(9,)``, dtype ``float64``. Raises ValueError If ``box`` does not have 3 or 9 numbers. Examples >>> as_box9((4.0, 5.0, 6.0)).tolist() \[0.0, 0.0, 0.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0\]

### seams.system.from_pydseams {#seams-system-from-pydseams}

```python
seams.system.from_pydseams(frame)
```

Copy a pydseamslib frame into a :class:`System`. ``pydseamslib`` (``import pydseams``) is the trajectory package. A frame exposes ``positions`` as ``(x, y, z)`` rows and ``box`` as three lengths or the LAMMPS restricted-triclinic numbers. This package is the ctypes binding of ``libseams``. The in-tree nanobind module ``dseams`` covers ``chill_plus``, ``q4q6``, neighbours, and rings. Use this function when a pydseams frame should run an analysis that module does not wrap. Parameters frame : object Any object with ``positions`` and ``box``. Returns System A new handle. The caller closes it. Raises TypeError If ``positions`` or ``box`` is missing.

### seams.system.label_name {#seams-system-label-name}

```python
seams.system.label_name(code)
```

Name of a CHILL label. Parameters code : int A ``seams_label`` value. Returns str ``”unknown”`` when the code is not a shipped label.

### seams.system.nframes {#seams-system-nframes}

```python
seams.system.nframes(path)
```

Number of frames in a trajectory path. Parameters path : path-like Returns int

### seams.system.num_threads {#seams-system-num-threads}

```python
seams.system.num_threads()
```

Team width the next analysis call will request. Returns int

### seams.system.set_num_threads {#seams-system-set-num-threads}

```python
seams.system.set_num_threads(n)
```

Set the OpenMP team used by later analysis calls. Parameters n : int Positive width, or 0 to restore the width captured at first use. Raises seams.SeamsError If ``n`` is negative. Notes The width is process-global. Do not call this concurrently with analysis.

### seams.system.status_string {#seams-system-status-string}

```python
seams.system.status_string(code)
```

Short name of a status code. Parameters code : int Returns str

### seams.system.version {#seams-system-version}

```python
seams.system.version()
```

Library version text. Returns str For example ``”2.14.3”``. The caller does not free the C string. This function returns a Python copy.

### seams.system.version_tuple {#seams-system-version-tuple}

```python
seams.system.version_tuple()
```

Major, minor, and patch. Returns tuple of int

## C API {#c-api}

Each native call and the `seams.h` or `seams_lang.h` function it uses. The C column is the stable ABI.

| Native | C |
| --- | --- |
| dseams.frames.cages | seams_frame_cages, seams_frame_free |
| dseams.frames.rings | seams_frame_free, seams_frame_rings, seams_rings_count, seams_rings_destroy, seams_rings_find, seams_rings_fraction, seams_rings_histogram, seams_rings_labels |
| dseams.frames.signature_cages | seams_frame_free, seams_frame_signature |
| dseams.frames.symmetry | seams_frame_free, seams_frame_symmetry |
| dseams.userframe.CageScore | no direct C call |
| dseams.userframe.CageScore.n_both | no direct C call |
| dseams.userframe.CageScore.n_ic | no direct C call |
| dseams.userframe.CageScore.n_ih | no direct C call |
| dseams.userframe.CageScore.n_water | no direct C call |
| dseams.userframe.Frame | seams_frame_box, seams_frame_cages, seams_frame_chill_plus, seams_frame_free, seams_frame_monolayer, seams_frame_nop, seams_frame_positions, seams_frame_prisms, seams_frame_rdf2d, seams_frame_rdf2d_nbin, seams_frame_rings, seams_frame_signature, seams_frame_signature_edge, seams_frame_signature_vertices, seams_frame_signature_vertices_edge, seams_frame_steinhardt, seams_last_error |
| dseams.userframe.Frame.box | seams_frame_box |
| dseams.userframe.Frame.cages | seams_frame_cages |
| dseams.userframe.Frame.cages_by_signature | seams_frame_signature, seams_frame_signature_edge, seams_frame_signature_vertices, seams_frame_signature_vertices_edge |
| dseams.userframe.Frame.chill_plus | seams_chill_plus, seams_frame_chill_plus, seams_system_count |
| dseams.userframe.Frame.close | seams_frame_free |
| dseams.userframe.Frame.find_prisms | seams_frame_prisms |
| dseams.userframe.Frame.monolayer_rings | seams_frame_monolayer |
| dseams.userframe.Frame.n_atoms | seams_frame_nop |
| dseams.userframe.Frame.positions | seams_frame_positions |
| dseams.userframe.Frame.rdf_2d | seams_frame_rdf2d, seams_frame_rdf2d_nbin |
| dseams.userframe.Frame.ring_census | seams_frame_rings |
| dseams.userframe.Frame.steinhardt | seams_frame_steinhardt |
| dseams.userframe.Frame.to_ase | no direct C call |
| dseams.userframe.Frame.to_solvis | no direct C call |
| dseams.userframe.IceCounts | no direct C call |
| dseams.userframe.from_ase | seams_frame_from_arrays |
| dseams.userframe.read | seams_frame_read |
| dseams.userframe.to_solvis | no direct C call |
| pydseams.aseio.frame_from_ase | no direct C call |
| pydseams.aseio.frame_to_ase | no direct C call |
| pydseams.band.band_rows | no direct C call |
| pydseams.band.bond_pairs | no direct C call |
| pydseams.band.ring_census | no direct C call |
| pydseams.band.shortest_within_type | no direct C call |
| pydseams.bond_network.components | no direct C call |
| pydseams.bond_network.residence_counts | no direct C call |
| pydseams.config.apply | no direct C call |
| pydseams.config.apply_file | no direct C call |
| pydseams.config.cutoff | no direct C call |
| pydseams.config.frame | no direct C call |
| pydseams.config.get_float | no direct C call |
| pydseams.config.get_int | no direct C call |
| pydseams.config.graph | no direct C call |
| pydseams.config.k_neighbors | no direct C call |
| pydseams.config.reset | no direct C call |
| pydseams.features.IceFeaturizer | seams_chill_plus, seams_system_count |
| pydseams.features.IceFeaturizer.feature_names | no direct C call |
| pydseams.features.IceFeaturizer.frame_features | seams_chill_plus, seams_system_count |
| pydseams.features.IceFeaturizer.transform | no direct C call |
| pydseams.features.count_frames | no direct C call |
| pydseams.features.discretize_nmax | no direct C call |
| pydseams.features.ion_environment | no direct C call |
| pydseams.features.ion_features | no direct C call |
| pydseams.features.to_deeptime | no direct C call |
| pydseams.features.to_pyemma_featurizer | no direct C call |
| pydseams.frame.CageScore | no direct C call |
| pydseams.frame.CageScore.n_ic | no direct C call |
| pydseams.frame.CageScore.n_ih | no direct C call |
| pydseams.frame.CageScore.n_water | no direct C call |
| pydseams.frame.ContactPairs | no direct C call |
| pydseams.frame.DensityProfile | no direct C call |
| pydseams.frame.DomainStats | no direct C call |
| pydseams.frame.DualOccupancy | no direct C call |
| pydseams.frame.Frame | seams_chill_plus, seams_system_count |
| pydseams.frame.Frame.bonds_by_index | no direct C call |
| pydseams.frame.Frame.box | no direct C call |
| pydseams.frame.Frame.cage_affiliation | no direct C call |
| pydseams.frame.Frame.cages | no direct C call |
| pydseams.frame.Frame.cages_by_signature | no direct C call |
| pydseams.frame.Frame.chill | no direct C call |
| pydseams.frame.Frame.chill_plus | seams_chill_plus, seams_system_count |
| pydseams.frame.Frame.classify_chill | no direct C call |
| pydseams.frame.Frame.classify_chill_plus | seams_chill_plus, seams_system_count |
| pydseams.frame.Frame.classify_templates | no direct C call |
| pydseams.frame.Frame.classify_topology | no direct C call |
| pydseams.frame.Frame.cn | no direct C call |
| pydseams.frame.Frame.density | no direct C call |
| pydseams.frame.Frame.domain | no direct C call |
| pydseams.frame.Frame.find_prisms | no direct C call |
| pydseams.frame.Frame.fingerprint | no direct C call |
| pydseams.frame.Frame.fit_classifier | no direct C call |
| pydseams.frame.Frame.from_arrays | no direct C call |
| pydseams.frame.Frame.from_ase | no direct C call |
| pydseams.frame.Frame.from_chemfiles | no direct C call |
| pydseams.frame.Frame.from_con | no direct C call |
| pydseams.frame.Frame.from_file | no direct C call |
| pydseams.frame.Frame.from_xyz | no direct C call |
| pydseams.frame.Frame.guest_occupancy | no direct C call |
| pydseams.frame.Frame.guest_occupancy_both | no direct C call |
| pydseams.frame.Frame.hbonds | no direct C call |
| pydseams.frame.Frame.hbonds_from_donors | no direct C call |
| pydseams.frame.Frame.hydration_shell_rings | no direct C call |
| pydseams.frame.Frame.ion_cloud | no direct C call |
| pydseams.frame.Frame.ion_environment | no direct C call |
| pydseams.frame.Frame.load_frame | no direct C call |
| pydseams.frame.Frame.local_structure_index | no direct C call |
| pydseams.frame.Frame.monolayer_rings | no direct C call |
| pydseams.frame.Frame.n_atoms | no direct C call |
| pydseams.frame.Frame.neighbor_list | no direct C call |
| pydseams.frame.Frame.pairs | no direct C call |
| pydseams.frame.Frame.positions | no direct C call |
| pydseams.frame.Frame.predict_class | no direct C call |
| pydseams.frame.Frame.prism_counts | no direct C call |
| pydseams.frame.Frame.rdf | no direct C call |
| pydseams.frame.Frame.rdf_2d | no direct C call |
| pydseams.frame.Frame.rings | seams_rings_count, seams_rings_destroy, seams_rings_find, seams_rings_fraction, seams_rings_histogram, seams_rings_labels |
| pydseams.frame.Frame.rings_recomputed_sources | no direct C call |
| pydseams.frame.Frame.running_cn | no direct C call |
| pydseams.frame.Frame.seeded_affiliation | no direct C call |
| pydseams.frame.Frame.soap | no direct C call |
| pydseams.frame.Frame.steinhardt | no direct C call |
| pydseams.frame.Frame.steinhardt_voronoi | no direct C call |
| pydseams.frame.Frame.tetrahedral_q | no direct C call |
| pydseams.frame.Frame.to_ase | no direct C call |
| pydseams.frame.Frame.to_solvis | no direct C call |
| pydseams.frame.Frame.topology_library | no direct C call |
| pydseams.frame.Frame.voronoi_features | no direct C call |
| pydseams.frame.IceCounts | no direct C call |
| pydseams.frame.read | no direct C call |
| pydseams.io.available_readers | no direct C call |
| pydseams.io.has_chemfiles | no direct C call |
| pydseams.io.has_readcon | no direct C call |
| pydseams.io.iter_trajectory | no direct C call |
| pydseams.io.read | no direct C call |
| pydseams.liquid_order.local_structure_index | no direct C call |
| pydseams.liquid_order.tetrahedral_q | no direct C call |
| pydseams.solvis.to_solvis | no direct C call |
| seams.errors.SeamsError | no direct C call |
| seams.results.BrineFronts | no direct C call |
| seams.results.BrinePockets | no direct C call |
| seams.results.BrineProfile | no direct C call |
| seams.results.BrineRejection | no direct C call |
| seams.results.CageCounts | no direct C call |
| seams.results.ChillResult | no direct C call |
| seams.results.Contacts | no direct C call |
| seams.results.Density | no direct C call |
| seams.results.Domains | no direct C call |
| seams.results.NeighbourCSR | no direct C call |
| seams.results.OrderResult | no direct C call |
| seams.results.Prisms | no direct C call |
| seams.results.Q4Q6 | no direct C call |
| seams.results.Q6Bar | no direct C call |
| seams.results.Radial | no direct C call |
| seams.results.RingLabels | no direct C call |
| seams.results.RunningCN | no direct C call |
| seams.results.Signature | no direct C call |
| seams.system.Brine | no direct C call |
| seams.system.Brine.fronts | no direct C call |
| seams.system.Brine.pockets | no direct C call |
| seams.system.Brine.profile | no direct C call |
| seams.system.Brine.rejection | no direct C call |
| seams.system.CageSet | no direct C call |
| seams.system.CageSet.counts | no direct C call |
| seams.system.CageSet.labels | no direct C call |
| seams.system.Fingerprint | no direct C call |
| seams.system.Fingerprint.atom | no direct C call |
| seams.system.Fingerprint.key | no direct C call |
| seams.system.Fingerprint.rings | seams_rings_count, seams_rings_destroy, seams_rings_find, seams_rings_fraction, seams_rings_histogram, seams_rings_labels |
| seams.system.Neighbours | no direct C call |
| seams.system.Neighbours.copy | no direct C call |
| seams.system.Neighbours.count | seams_system_count |
| seams.system.Rings | seams_system_count |
| seams.system.Rings.count | seams_system_count |
| seams.system.Rings.fraction | no direct C call |
| seams.system.Rings.histogram | no direct C call |
| seams.system.Rings.labels | seams_system_count |
| seams.system.Rings.vertices | no direct C call |
| seams.system.System | seams_chill_plus, seams_q4q6, seams_system_count |
| seams.system.System.brine | seams_system_count |
| seams.system.System.cages | no direct C call |
| seams.system.System.chill | seams_system_count |
| seams.system.System.chill_plus | seams_chill_plus, seams_system_count |
| seams.system.System.contact_pairs | no direct C call |
| seams.system.System.coordination | no direct C call |
| seams.system.System.count | seams_system_count |
| seams.system.System.density | no direct C call |
| seams.system.System.domains | no direct C call |
| seams.system.System.fingerprint | no direct C call |
| seams.system.System.hbonds | no direct C call |
| seams.system.System.knn | no direct C call |
| seams.system.System.monolayer | no direct C call |
| seams.system.System.neighbours | seams_neighbours_build, seams_neighbours_copy, seams_neighbours_destroy |
| seams.system.System.prisms | no direct C call |
| seams.system.System.q4q6 | seams_q4q6, seams_system_count |
| seams.system.System.q6bar | seams_system_count |
| seams.system.System.ql | seams_system_count |
| seams.system.System.rdf | no direct C call |
| seams.system.System.read | seams_system_count |
| seams.system.System.rings | seams_rings_count, seams_rings_destroy, seams_rings_find, seams_rings_fraction, seams_rings_histogram, seams_rings_labels |
| seams.system.System.running_cn | no direct C call |
| seams.system.System.set_frame | seams_system_set_frame |
| seams.system.System.set_ortho | no direct C call |
| seams.system.System.signature | no direct C call |
| seams.system.System.site_table | seams_system_count |
| seams.system.System.steinhardt | seams_system_count |
| seams.system.System.symmetry | seams_system_count |
| seams.system.abi | no direct C call |
| seams.system.as_box9 | no direct C call |
| seams.system.from_pydseams | seams_system_set_frame |
| seams.system.label_name | no direct C call |
| seams.system.nframes | no direct C call |
| seams.system.num_threads | seams_get_num_threads |
| seams.system.set_num_threads | seams_get_num_threads, seams_set_num_threads |
| seams.system.status_string | no direct C call |
| seams.system.version | seams_get_num_threads, seams_set_num_threads, seams_version_string |
| seams.system.version_tuple | no direct C call |
