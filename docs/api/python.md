# Python

The package `dseams` is a nanobind module over `seams.h`. NumPy arrays are
C-contiguous and are read in place.

```python
import dseams
labels, counts = dseams.chill_plus(positions, box)
q4, q6 = dseams.q4q6(positions, box)
offsets, cols = dseams.neighbours(positions, box)
sizes, symmetry, hist, fraction = dseams.rings(positions, box)
```

`box` is three lengths or nine LAMMPS numbers. `dseams.PUBLISHED` and
`dseams.OVITO` select the CHILL+ rule. `dseams.CUTOFF` and `dseams.MUTUAL`
select the graph. `dseams.set_num_threads` and `dseams.num_threads` match the
C calls. `dseams.version` is the library version.

`dseams.from_ase(atoms)` rotates an ASE cell into the LAMMPS frame and returns
positions and the nine-number box. `dseams.from_mdanalysis(universe)` reads
`universe.dimensions` (lengths and angles) and `universe.atoms.positions`.

`python -m dseams` labels one cubic cell.
