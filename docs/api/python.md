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
C calls. `dseams.version` is the library version.

`dseams.from_ase(atoms)` rotates an ASE cell into the LAMMPS frame and returns
positions and the nine-number box. `dseams.from_mdanalysis(universe)` reads
`universe.dimensions` (lengths and angles) and `universe.atoms.positions`.

`python -m dseams` labels one cubic cell.

`cages`, `signature_cages`, and `symmetry` load `libseams_lang` through
`frames.py`. Set `SEAMS_LANG_LIB` to that shared library. The nanobind
calls (`chill_plus`, `q4q6`, `neighbours`, `rings`, `Frame`) use `libseams`.

`meson test docs_api_python` runs the eight-atom program below. It calls
`dseams.chill_plus` and `dseams.q4q6` when the nanobind module imports, and
the same `libseams` functions through ctypes otherwise.
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
