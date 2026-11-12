# Python tutorial

The nanobind module reads a C-contiguous NumPy array in place. This is the eight-atom diamond cell (nearest-neighbour spacing 2.75 Å). Every oxygen is cubic.

```python
import math
import numpy as np
import dseams

a = 2.75 * 4.0 / math.sqrt(3.0)
fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
basis = ((0, 0, 0), (0.25, 0.25, 0.25))
pos = np.array(
    [
        ((shift[0] + site[0]) * a, (shift[1] + site[1]) * a, (shift[2] + site[2]) * a)
        for shift in fcc
        for site in basis
    ],
    dtype=np.float64,
)
labels, counts = dseams.chill_plus(pos, (a, a, a))
print(dseams.version(), int(counts[0]), labels.shape[0])
```

`counts[0]` is 8 and every label is 0.

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
