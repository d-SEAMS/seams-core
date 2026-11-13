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

## The rest of the call set

`docs/examples/python_api.py` is the program `meson test docs_api_python` runs.
When the nanobind module is absent it calls `libseams` through ctypes.
On the eight-atom cell that program prints:

```
python call q6bar 0.628539361055
python call cubic 8
python call cages 0 0 8 0
python call cn 4.57142857143
python call hbond 0
python call sheet 0
python call vert 4 0
```

`q6bar` is `seams_steinhardt` with the Lechner-Dellago buffer.
`cages` is `seams_cages` (`n_hex`, `n_cubic`, `n_water`, `n_both`).
`cn` is `seams_coordination`. `hbond` is `seams_hbonds`.
`sheet` is the six-ring bin of `seams_monolayer`.
`vert` is the first ring from `seams_rings_vertices`.
The same run also calls `seams_rdf`, `seams_running_cn`, `seams_density_axis`,
`seams_fingerprint`, `seams_contact_pairs`, `seams_domains`, `seams_prisms`,
`seams_neighbours_knn`, and `seams_site_table`.
On `input/traj/genice_sI.lammpstrj` the `512` signature prints:

```
python call sig 2 40
```
