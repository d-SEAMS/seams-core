# Python tutorial

The nanobind module reads a C-contiguous NumPy array in place. This is the eight-atom diamond cell (nearest-neighbour spacing 2.75 Å). Every oxygen is cubic.

`meson test python_dseams` runs this method. On the eight-atom cell `counts[0]` equals the atom count and every label is 0.

```{literalinclude} ../../../bindings/python/tests/test_dseams.py
:language: python
:lines: 39-62
```

An ASE `Atoms` or an MDAnalysis universe becomes the same arrays:

```python
pos, box = dseams.from_ase(atoms)
pos, box = dseams.from_mdanalysis(universe)
labels, counts = dseams.chill_plus(pos, box, rule=dseams.OVITO)
```

```bash
PYTHONPATH=bindings/python python -m unittest -v bindings/python/tests/test_dseams.py
python -m dseams
```

The verbose run reports three cases: both CHILL+ rules with q4/q6 and rings, ASE and MDAnalysis frames, and a triclinic cell.
