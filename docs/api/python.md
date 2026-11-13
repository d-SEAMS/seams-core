# Python

Three packages call the same engine.

`seams` is the ctypes package in `bindings/python`. It wraps `include/seams.h`.
NumPy arrays go in and come back out. Results are frozen dataclasses.
A failing call raises `seams.SeamsError`. `status` is the `seams_status`
code and `message` is the Python copy of `seams_handle_error`.

`dseams` is the in-tree nanobind module. It covers `chill_plus`, `q4q6`,
neighbours, rings, and `Frame`. `pydseamslib` (`import pydseams`) is the
trajectory package on PyPI. `seams.from_pydseams(frame)` copies
`frame.positions` and `frame.box` into a `System` when an analysis lives
only on the C API.

Set `SEAMS_LIB` to `libseams.so`. Importing `seams` does not load the library.

## Quickstart

```python
import numpy as np
from seams import PUBLISHED, System, set_num_threads

set_num_threads(1)
with System() as sys:
    sys.set_ortho(positions, (edge, edge, edge))
    ice = sys.chill_plus(rule=PUBLISHED)
    order = sys.q4q6()
    print(int(ice.counts[0]), float(order.q6.mean()))
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
`q4q6`, `q6bar`, `ql`, and `steinhardt` share one process-global scratch
buffer. `set_num_threads` is process-global too.

## Native reference

Docstrings are numpydoc. `sphinx.ext.napoleon` and `sphinx.ext.autodoc`
render `seams` when the package imports, the same split rgpot uses
between the narrative page and the generated API. `mypy --strict` and
`pyright` check `bindings/python/src/seams`. `py.typed` ships with the
package. Doctests cover `as_box9` and `SeamsError`.

The nanobind module and the language-binding examples stay below.

```{literalinclude} ../../../docs/examples/seams_ctypes_quick.py
:language: python
```

```{literalinclude} ../../../docs/examples/minimal_py.py
:language: python
```

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
