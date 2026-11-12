# Python tutorial

The nanobind module reads a C-contiguous NumPy array in place.

```python
import math
import numpy as np
import dseams

a = 2.75 * 4.0 / math.sqrt(3.0)
pos = np.zeros((8, 3))
labels, counts = dseams.chill_plus(pos, (a, a, a))
print(dseams.version(), int(counts[0]))
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
