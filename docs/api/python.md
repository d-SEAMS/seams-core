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

`box` is three lengths or nine LAMMPS numbers. `dseams.PUBLISHED` (`"paper"`)
averages q_lm over the four closest neighbours. `dseams.OVITO` (`"ovito"`)
averages q_lm over every neighbour inside the cutoff. `dseams.CUTOFF` and `dseams.MUTUAL`
select the graph. `dseams.set_num_threads` and `dseams.num_threads` match the
C calls. `dseams.version` is the library version. `Frame.set_frame` stores one
cell. `Frame.count` is the number of points.

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
between the narrative page and the generated API. The generated page is
{doc}`seams_ctypes`. `mypy --strict` and `pyright` check
`bindings/python/src/seams`. `py.typed` ships with the package.
Doctests cover `as_box9` and `SeamsError`.

The nanobind module and the language-binding examples stay below.

```{literalinclude} ../../../docs/examples/seams_ctypes_quick.py
:language: python
```

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

Purpose: CHILL+ definition that averages q_lm over every neighbour inside the cutoff.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_py.py`

### PUBLISHED

Purpose: CHILL+ definition that averages q_lm over the four closest neighbours.

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


## Structure still

```{figure} ../tutorials/molecules/python/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug python` writes the 2400 by 1800 PNG.
```

