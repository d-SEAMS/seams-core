# Python

The package `dseams` is a nanobind module over `seams.h`. NumPy arrays are
C-contiguous and are read in place.

```python
import dseams
labels, counts = dseams.chill_plus(positions, box)
q4, q6 = dseams.q4q6(positions, box)
offsets, cols = dseams.neighbours(positions, box)
sizes, symmetry, hist, fraction = dseams.rings(positions, box)
per_atom, (n_hex, n_cubic, n_water, n_both) = dseams.cages(positions, box)
n_cages, n_atoms = dseams.signature_cages(positions, box, spec="512")
shell = dseams.symmetry(positions, box)
```

`dseams.cages` is the six-ring ice score from `frames.cages`: a per-atom code
and the hexagonal, cubic, water, and both counts. `dseams.signature_cages`
counts closed polyhedra whose ring census matches `spec` and returns
`(n_cages, n_atoms)`. `dseams.symmetry` is the shell label: 0 ice Ih, 1 ice
Ic, 2 stacking-disordered, and -1 otherwise.

`box` is three lengths or nine LAMMPS numbers. `dseams.PUBLISHED` and
`dseams.OVITO` select the CHILL+ rule. `dseams.CUTOFF` and `dseams.MUTUAL`
select the graph. `dseams.set_num_threads` and `dseams.num_threads` match the
C calls. `dseams.version` is the library version. `Frame.set_frame` stores one
cell. `Frame.count` is the number of points.

`docs/examples/minimal_py.py` reads
`input/traj/genice_sI.lammpstrj` and calls the C API for neighbours, q6,
CHILL+, and rings. `meson test docs_examples` runs it.

```{literalinclude} ../../../docs/examples/minimal_py.py
:language: python
```

`dseams.from_ase(atoms)` rotates an ASE cell into the LAMMPS frame and returns
positions and the nine-number box. `dseams.from_mdanalysis(universe)` reads
`universe.dimensions` (lengths and angles) and `universe.atoms.positions`.

`python -m dseams` labels one cubic cell.

`cages`, `signature_cages`, and `symmetry` load `libseams_lang` through
`frames.py`. Set `SEAMS_LANG_LIB` to that shared library. The nanobind
calls (`chill_plus`, `q4q6`, `neighbours`, `rings`, `Frame`) use `libseams`.

`meson test docs_api_python` runs the eight-atom program below. It calls
`dseams.chill_plus` and `dseams.q4q6` when the nanobind module imports.
The ctypes fallback calls those two and also `seams_system_read`,
`seams_chill`, `seams_steinhardt`, `seams_cages`, `seams_signature_cages`,
`seams_neighbours_knn`, `seams_rdf`, `seams_running_cn`, `seams_coordination`,
`seams_hbonds`, `seams_density_axis`, `seams_fingerprint`,
`seams_contact_pairs`, `seams_domains`, `seams_prisms`, `seams_monolayer`,
`seams_rings_vertices`, and `seams_site_table`.
`meson test python_dseams` runs `bindings/python/tests/test_dseams.py`.
`tests/test_frames_abi.py` is the `frames_abi` check inside
`tests/run_language_fronts.sh` (`meson test language_fronts`).

```{literalinclude} ../../../docs/examples/python_api.py
:language: python
```

```{literalinclude} ../../../bindings/python/tests/test_dseams.py
:language: python
:lines: 41-74,80-93
```

```{literalinclude} ../../../tests/test_frames_abi.py
:language: python
:lines: 50-57
```

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### CUTOFF

Purpose: Graph id of the cutoff neighbour list.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_py.py`

### Frame

Purpose: One stored cell and the classifiers that read it.

Parameters: none

Units: object

Example: `docs/examples/minimal_py.py`

### MUTUAL

Purpose: Graph id of the mutual neighbour list.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_py.py`

### OVITO

Purpose: CHILL+ rule id matching OVITO.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_py.py`

### PUBLISHED

Purpose: CHILL+ rule id of the 2015 paper.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_py.py`

### cages

Purpose: Score hexagonal and double-diamond cage membership.

Parameters: positions; box; or a Frame

Units: coordinates and cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_py.py`

### chill_plus

Purpose: Label each atom with CHILL+.

Parameters: positions; box; rule; graph; cutoff

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### count

Purpose: Number of points stored on a Frame.

Parameters: frame

Units: dimensionless count

Example: `docs/examples/minimal_py.py`

### from_ase

Purpose: Turn an ASE atoms object into positions and a LAMMPS cell.

Parameters: atoms

Units: positions and box lengths in angstrom

Example: `docs/examples/minimal_py.py`

### from_mdanalysis

Purpose: Turn an MDAnalysis universe into positions and a LAMMPS cell.

Parameters: universe; select

Units: positions in angstrom; angles inside the universe are degrees

Example: `docs/examples/minimal_py.py`

### neighbours

Purpose: Cutoff or mutual neighbour list as CSR.

Parameters: positions; box; cutoff; graph

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### num_threads

Purpose: OpenMP width the next call will use.

Parameters: none

Units: dimensionless thread count

Example: `docs/examples/minimal_py.py`

### q4q6

Purpose: Local q4 and q6 on the cutoff shell.

Parameters: positions; box; cutoff

Units: coordinates and cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_py.py`

### rings

Purpose: Primitive-ring sizes, symmetry, histogram, and fraction.

Parameters: positions; box; cutoff; max_depth

Units: coordinates and cutoff in angstrom; max_depth is a dimensionless ring size

Example: `docs/examples/minimal_py.py`

### set_frame

Purpose: Store one cell on a Frame.

Parameters: frame; positions; box

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### set_num_threads

Purpose: Set the OpenMP width for the next call.

Parameters: n

Units: dimensionless thread count; 0 restores the library default

Example: `docs/examples/minimal_py.py`

### signature_cages

Purpose: Count closed cages matching a ring-census spec.

Parameters: positions; box; spec

Units: coordinates and cutoff in angstrom; spec is text such as 512

Example: `docs/examples/minimal_py.py`

### symmetry

Purpose: Shell label for Ih, Ic, stacking disorder, or other.

Parameters: positions; box

Units: coordinates and cutoff in angstrom; labels are dimensionless codes

Example: `docs/examples/minimal_py.py`

### version

Purpose: Library version text.

Parameters: none

Units: text

Example: `docs/examples/minimal_py.py`
<!-- API_SYMBOLS_END -->
