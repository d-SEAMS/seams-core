# Python

Three packages call the same engine.

`seams` is the ctypes package in `bindings/python`. It wraps `include/seams.h`.
NumPy arrays go in and come back out. Results are frozen dataclasses.
A failing call raises `seams.SeamsError`. `status` is the `seams_status`
code and `message` is the Python copy of `seams_handle_error`.

`dseams` is the in-tree nanobind module over `seams.h`. It covers
`chill_plus`, `q4q6`, neighbours, rings, and `Frame`. NumPy arrays are
C-contiguous and are read in place. `pydseamslib` (`import pydseams`) is the
trajectory package on PyPI. `seams.from_pydseams(frame)` copies
`frame.positions` and `frame.box` into a `System` when an analysis lives
only on the C API.

Set `SEAMS_LIB` to `libseams.so`. Importing `seams` does not load the library.

## Quickstart

### `seams`

```python
import math

import numpy as np
from seams import PUBLISHED, System, set_num_threads

edge = 2.75 * 4.0 / math.sqrt(3.0)
fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
basis = ((0, 0, 0), (0.25, 0.25, 0.25))
positions = np.array([[(s + b) * edge for s, b in zip(shift, site)] for shift in fcc for site in basis])

set_num_threads(1)
with System() as sys:
    sys.set_ortho(positions, (edge, edge, edge))
    ice = sys.chill_plus(rule=PUBLISHED)
    order = sys.q4q6()
    print(int(ice.counts[0]), float(order.q6.mean()))
    assert int(ice.counts[0]) == 8
    assert abs(float(order.q6.mean()) - 0.628539) < 1e-6
```

`docs/examples/seams_ctypes_quick.py` runs that cell. The program prints:

```
version 2.14.0
n 8
cubic 8
q6 0.628539
degree 4
six 0
fraction 0.0
cn 4.571429
```

Mean q6 is 0.628539 and each atom has four neighbours. This cell holds
no primitive six-ring. It is thinner than twice the cutoff, so the
cage census is refused and the script leaves that call out. A 2 by 2
by 2 supercell, which `python_seams_ctypes` also runs, has a six-ring
fraction above 0.5. `meson test python_seams_ctypes` runs the pytest
file and the doctests.

### `dseams`

The same eight sites through the nanobind module, in a box three cells
wide so that the cage calls are not refused:

```python
import math

import dseams
import numpy as np

edge = 2.75 * 4.0 / math.sqrt(3.0)
sites = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5),
         (0.25, 0.25, 0.25), (0.75, 0.75, 0.25), (0.75, 0.25, 0.75), (0.25, 0.75, 0.75))
positions = np.array([[x * edge, y * edge, z * edge] for x, y, z in sites])
box = (3 * edge, 3 * edge, 3 * edge)
labels, counts = dseams.chill_plus(positions, box)
q4, q6 = dseams.q4q6(positions, box)
offsets, cols = dseams.neighbours(positions, box)
sizes, symmetry, hist, fraction = dseams.rings(positions, box)
per_atom, (n_hex, n_cubic, n_water, n_both) = dseams.cages(positions, box)
n_cages, n_atoms = dseams.signature_cages(positions, box, spec="512")
shell = dseams.symmetry(positions, box)
assert len(labels) == 8 and int(counts.sum()) == 8
assert len(q4) == 8 and len(q6) == 8
assert len(offsets) == 9
assert n_hex + n_cubic + n_water + n_both == 8
assert (n_cages, n_atoms) == (0, 0)
assert len(shell) == 8
```

`dseams.cages` is the six-ring ice score from `frames.cages`: a per-atom code
and the hexagonal, cubic, water, and both counts. `dseams.signature_cages`
counts closed polyhedra whose ring census matches `spec` and returns
`(n_cages, n_atoms)`. `dseams.symmetry` is the shell label: 0 ice Ih, 1 ice
Ic, 2 stacking-disordered, and -1 otherwise.

`box` is three lengths or nine LAMMPS numbers. `dseams.PUBLISHED` (`"paper"`)
averages q_lm over the four closest neighbours. `dseams.OVITO` (`"ovito"`)
averages q_lm over every neighbour inside the cutoff. `dseams.CUTOFF` and `dseams.MUTUAL`
select the graph. `dseams.set_num_threads` and `dseams.num_threads` match the
C calls. `dseams.version` is the library version. `Frame.set_frame` stores one
cell. `Frame.count` is the number of points.

## Ownership

`System`, `Neighbours`, `Rings`, `CageSet`, `Fingerprint`, and `Brine`
are context managers. `close` frees the handle. A second `close` does
nothing. `__del__` frees a handle that left the `with` block by an
exception before `close`. Neighbour, ring, cage, fingerprint, and brine
handles keep their own copies, so they stay valid after the `System`
is closed. `Neighbours.copy` fills caller-owned NumPy arrays.
`seams_neighbours_csr` is deprecated and is not wrapped.

## Errors

```python
from seams import SeamsError, System

with System() as sys:
    try:
        sys.chill_plus()
    except SeamsError as err:
        print(err.status, err.message)
        assert err.status == 5
```

An empty frame raises status 5 (`SEAMS_ERR_STATE`). The message is the
handle string. A box that is neither 3 nor 9 numbers raises `ValueError`
before the library is called.

## Analyses

| Call | Result |
| --- | --- |
| `chill`, `chill_plus` | `ChillResult.labels`, `counts` (length 9) |
| `ql`, `steinhardt` | `OrderResult.ql`, `qlbar` |
| `q4q6`, `q6bar` | `Q4Q6`, `Q6Bar` |
| `neighbours`, `knn` | `Neighbours`, then `NeighbourCSR` |
| `rings` | histogram, labels, vertices, fraction |
| `cages`, `signature` | `CageCounts`, `Signature` |
| `symmetry` | 0 Ih, 1 Ic, 2 stacking, -1 other |
| `rdf`, `running_cn`, `coordination` | `Radial`, `RunningCN`, float |
| `hbonds`, `density` | count, `Density` |
| `fingerprint` | key text and ring census |
| `contact_pairs`, `domains`, `prisms`, `monolayer`, `site_table` | dataclasses or arrays |
| `brine` | fronts, rejection, profile, pockets |

`CUTOFF` is 0, `MUTUAL` is 1, `SEEDED` is 2. `PUBLISHED` is 0, `OVITO` is 1.
`q4q6`, `q6bar`, `ql`, and `steinhardt` keep their scratch per calling thread, so
separate frames may run them in separate threads. `set_num_threads` is process-global.

## Native reference

Docstrings are numpydoc. `sphinx.ext.napoleon` and `sphinx.ext.autodoc`
render `seams` when the package imports, the same split rgpot uses
between the narrative page and the generated API. The generated page is
{doc}`seams_ctypes`. `mypy --strict` and `pyright` check
`bindings/python/src/seams`. `py.typed` ships with the package.
Doctests cover `as_box9` and `SeamsError`.

The nanobind module and the language-binding examples stay below.

```{literalinclude} ../../../docs/examples/seams_ctypes_quick.py
:language: python
```

`docs/examples/core_counts.py` is the Ice Ih and sII census. Sphinx
doctest and `meson test docs_core_counts` both run it. `SEAMS_LIBDIR`
is the directory that holds `libseams.so`.

```{eval-rst}
.. doctest::

   >>> import core_counts
   >>> core_counts.check()
   hexagonal 128
   rings6 256
   clathrate 136
   rings5 144
```

`docs/examples/minimal_py.py` reads
`input/traj/genice_sI.lammpstrj` and calls the C API for neighbours, q6,
CHILL+, and rings. `meson test docs_examples` runs it.

```{literalinclude} ../../../docs/examples/minimal_py.py
:language: python
```

`dseams.from_ase(atoms)` rotates an ASE cell into the LAMMPS frame and returns
positions and the nine-number box. `dseams.from_mdanalysis(universe)` reads
`universe.dimensions` (lengths and angles) and `universe.atoms.positions`.

```python
# excerpt: bindings/python/dseams/__init__.py
def from_ase(atoms):
```

```python
# excerpt: bindings/python/dseams/__init__.py
def from_mdanalysis(universe, select="all"):
```

`python -m dseams` labels one cubic cell.

`cages`, `signature_cages`, and `symmetry` load `libseams_lang` through
`frames.py`. Set `SEAMS_LANG_LIB` to that shared library. The nanobind
calls (`chill_plus`, `q4q6`, `neighbours`, `rings`, `Frame`) use `libseams`.

`meson test docs_api_python` runs the eight-atom program below. On that cell the ctypes path prints `python call cubic 8`, `python call q6bar 0.628539361055`, `python call cages 0 0 8 0`, `python call cn 4.57142857143`, `python call hbond 0`, `python call sheet 0`, and `python call vert 4 0`. The `512` signature on the sI guest frame prints `python call sig 2 40`. It calls
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

Returns: The integer code.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/minimal_py.py`

### Frame

Purpose: One stored cell and the classifiers that read it.

Parameters: none

Returns: A Frame object.

Ownership: The object owns its seams_system and frees it when the object is destroyed.

Errors: None on construction. Later methods raise RuntimeError with seams_last_error text.

Units: object

Example: `docs/examples/minimal_py.py`

### Frame.chill_plus

Purpose: CHILL+ labels and the nine population counts for the atoms on a dseams.Frame.

Parameters: self; cutoff; graph; rule

Returns: A pair of NumPy arrays, labels and counts.

Ownership: The caller owns the returned arrays. The frame still owns the seams_system.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### Frame.count

Purpose: Return the number of points stored on a dseams.Frame.

Parameters: self

Returns: The atom count.

Ownership: No allocation.

Errors: None.

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### Frame.neighbours

Purpose: CSR neighbour list of a dseams.Frame, self entries omitted.

Parameters: self; cutoff; graph

Returns: A pair of NumPy arrays, offsets and columns.

Ownership: The caller owns the returned arrays. The temporary C neighbour list is freed inside the call.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### Frame.q4q6

Purpose: Local q4 and q6 on the cutoff shell of a dseams.Frame.

Parameters: self; cutoff

Returns: A pair of NumPy arrays, q4 and q6.

Ownership: The caller owns the returned arrays.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### Frame.rings

Purpose: Ring sizes, symmetry labels, a size histogram, and the six-ring fraction.

Parameters: self; cutoff; max_depth

Returns: Sizes, symmetry, histogram, and the six-ring fraction.

Ownership: The caller owns the returned arrays. The C ring set is freed inside the call.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### Frame.set_frame

Purpose: Store positions and a cell on a dseams.Frame.

Parameters: positions; box

Returns: Nothing.

Ownership: The frame copies the coordinates into the seams_system it owns.

Errors: ValueError when box is not 3 lengths or 9 LAMMPS numbers. RuntimeError with seams_last_error text when the C call fails.

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### MUTUAL

Purpose: Graph id of the mutual neighbour list.

Parameters: none

Returns: The integer code.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/minimal_py.py`

### OVITO

Purpose: CHILL+ definition that averages q_lm over every neighbour inside the cutoff.

Parameters: none

Returns: The integer code.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/minimal_py.py`

### PUBLISHED

Purpose: CHILL+ definition that averages q_lm over the four closest neighbours.

Parameters: none

Returns: The integer code.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/minimal_py.py`

### cages

Purpose: Score hexagonal and double-diamond cage membership.

Parameters: positions; box; or a Frame

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: coordinates and cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_py.py`

### chill_plus

Purpose: Label each atom with CHILL+.

Parameters: positions; box; rule; graph; cutoff

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### count

Purpose: Number of points stored on a Frame.

Parameters: frame

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: dimensionless count

Example: `docs/examples/minimal_py.py`

### from_ase

Purpose: Turn an ASE atoms object into positions and a LAMMPS cell.

Parameters: atoms

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: positions and box lengths in angstrom

Example: `docs/examples/minimal_py.py`

### from_mdanalysis

Purpose: Turn an MDAnalysis universe into positions and a LAMMPS cell.

Parameters: universe; select

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: positions in angstrom; angles inside the universe are degrees

Example: `docs/examples/minimal_py.py`

### neighbours

Purpose: Cutoff or mutual neighbour list as CSR.

Parameters: positions; box; cutoff; graph

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### num_threads

Purpose: OpenMP width the next call will use.

Parameters: none

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: dimensionless thread count

Example: `docs/examples/minimal_py.py`

### q4q6

Purpose: Local q4 and q6 on the cutoff shell.

Parameters: positions; box; cutoff

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: coordinates and cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_py.py`

### rings

Purpose: Primitive-ring sizes, symmetry, histogram, and fraction.

Parameters: positions; box; cutoff; max_depth

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: coordinates and cutoff in angstrom; max_depth is a dimensionless ring size

Example: `docs/examples/minimal_py.py`

### set_frame

Purpose: Store one cell on a Frame.

Parameters: frame; positions; box

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: coordinates and cutoff in angstrom

Example: `docs/examples/minimal_py.py`

### set_num_threads

Purpose: Set the OpenMP width for the next call.

Parameters: n

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: dimensionless thread count; 0 restores the library default

Example: `docs/examples/minimal_py.py`

### signature_cages

Purpose: Count closed cages matching a ring-census spec.

Parameters: positions; box; spec

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: coordinates and cutoff in angstrom; spec is text such as 512

Example: `docs/examples/minimal_py.py`

### symmetry

Purpose: Shell label for Ih, Ic, stacking disorder, or other.

Parameters: positions; box

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: coordinates and cutoff in angstrom; labels are dimensionless codes

Example: `docs/examples/minimal_py.py`

### userframe.CageScore

Purpose: Boolean membership of the hexagonal and double-diamond cages.

Parameters: hc; ddc

Returns: A CageScore.

Ownership: The caller owns the object.

Errors: None.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.CageScore.n_both

Purpose: Number of atoms in both cages.

Parameters: self

Returns: The count.

Ownership: No allocation.

Errors: None.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.CageScore.n_ic

Purpose: Number of atoms in the double-diamond cage.

Parameters: self

Returns: The count.

Ownership: No allocation.

Errors: None.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.CageScore.n_ih

Purpose: Number of atoms in the hexagonal cage.

Parameters: self

Returns: The count.

Ownership: No allocation.

Errors: None.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.CageScore.n_water

Purpose: Number of atoms in neither cage.

Parameters: self

Returns: The count.

Ownership: No allocation.

Errors: None.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame

Purpose: One configuration whose classifiers call libseams_lang.

Parameters: ptr; symbols; cutoff

Returns: A Frame.

Ownership: The object owns the seams_frame pointer. close() or deletion calls seams_frame_free.

Errors: RuntimeError when ptr is null.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.box

Purpose: Cell stored on the frame.

Parameters: self

Returns: Three edge lengths, or nine LAMMPS numbers when a tilt is set.

Ownership: The caller owns the returned array.

Errors: None.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.cages

Purpose: Hexagonal and double-diamond membership from seams_frame_cages.

Parameters: self; seeded; k; ring_adjacent

Returns: A CageScore.

Ownership: The caller owns the score. The frame keeps the same object.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.cages_by_signature

Purpose: Closed cages of one ring census.

Parameters: self; signature; max_ring_size; edge

Returns: A list whose length is the cage count. n_atoms is the vertex count.

Ownership: The caller owns the list. The C call does not return one vertex list per cage.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.chill_plus

Purpose: Published CHILL+ histogram for this frame.

Parameters: self; graph

Returns: An IceCounts mapping of label name to population.

Ownership: The caller owns the dictionary. Labels are also stored on the frame.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.close

Purpose: Free the libseams_lang frame.

Parameters: self

Returns: Nothing.

Ownership: Frees the pointer this Frame owns. A second call is a no-op.

Errors: None.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.find_prisms

Purpose: Prism census from seams_frame_prisms.

Parameters: self; output_dir; max_depth; shape_matching

Returns: A list of size, perfect, and imperfect rows.

Ownership: The caller owns the list.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.monolayer_rings

Purpose: Quasi-2D ring counts from seams_frame_monolayer.

Parameters: self; output_dir; sheet_area; max_depth

Returns: A dictionary of size to count and coverage.

Ownership: The caller owns the dictionary. The call writes topoMonolayer under output_dir.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.n_atoms

Purpose: Number of atoms on the frame.

Parameters: self

Returns: The count from seams_frame_nop.

Ownership: No allocation.

Errors: None.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.positions

Purpose: Copy of the coordinates stored on the frame.

Parameters: self

Returns: An (n, 3) NumPy array in angstrom.

Ownership: The caller owns the copy. The frame keeps its own buffer.

Errors: None.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.rdf_2d

Purpose: In-plane g(r) from seams_frame_rdf2d.

Parameters: self; output_dir; cutoff; binwidth

Returns: Radius and g(r) lists.

Ownership: The caller owns the lists. The call also writes under output_dir.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.ring_census

Purpose: Primitive rings by size from seams_frame_rings.

Parameters: self; max_size

Returns: A dictionary of size to count, empty sizes omitted.

Ownership: The caller owns the dictionary.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.steinhardt

Purpose: Local q_l and the Lechner-Dellago average from seams_frame_steinhardt.

Parameters: self; order_l

Returns: A dictionary with ql and ql_bar lists.

Ownership: The caller owns the dictionary.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.to_ase

Purpose: ASE Atoms for this frame.

Parameters: self

Returns: An ASE Atoms object. Ice labels and cage flags are arrays when those calls have run.

Ownership: The caller owns the Atoms object.

Errors: ImportError when ASE is not installed.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.Frame.to_solvis

Purpose: Solvis System built from to_ase.

Parameters: self; expand_box

Returns: A solvis System.

Ownership: The caller owns the System.

Errors: ImportError when solvis-tools is not installed.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.IceCounts

Purpose: CHILL+ population counts keyed by label name.

Parameters: mapping of name to count

Returns: An IceCounts dictionary.

Ownership: The caller owns the dictionary.

Errors: None.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.from_ase

Purpose: Frame for one ASE species. Coordinates go through seams_frame_from_arrays.

Parameters: atoms; select; cutoff; bonded

Returns: A userframe.Frame.

Ownership: The caller owns the Frame.

Errors: TypeError when atoms has no positions. ValueError when the selection matches nothing.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.read

Purpose: Read one frame with seams_frame_read. pydseams.read is this function.

Parameters: filename; frame; atom_type; all_atoms

Returns: A userframe.Frame.

Ownership: The caller owns the Frame and frees it with close() or by dropping it.

Errors: RuntimeError whose text is seams_last_error.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### userframe.to_solvis

Purpose: Call to_solvis on a frame.

Parameters: frame; expand_box

Returns: The solvis System returned by the frame.

Ownership: The caller owns the System.

Errors: ImportError when solvis-tools is not installed.

Units: coordinates and cutoff in angstrom

Example: `docs/orgmode/reference/python-frame.org`

### version

Purpose: Library version text.

Parameters: none

Returns: The arrays or counts named in the purpose. NumPy arrays are owned by the caller.

Ownership: A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.

Errors: RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.

Units: text

Example: `docs/examples/minimal_py.py`
<!-- API_SYMBOLS_END -->


## Structure still

```{figure} ../tutorials/molecules/python/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug python` writes the 2400 by 1800 PNG.
```

