# pydseamslib

Python bindings for [d-SEAMS](https://github.com/d-SEAMS/seams-core), the
Deferred Structural Elucidation Analysis for Molecular Simulations engine.
The distribution is `pydseamslib`; the package you import is `dseams`.

```bash
pip install pydseamslib
```

```python
import numpy as np
import dseams

print(dseams.__version__)
frame = dseams.Frame()
frame.set_frame(positions, (lx, ly, lz))   # (n, 3) float64, Angstrom
labels, counts = frame.chill_plus()
q4, q6 = frame.q4q6()
```

`dseams` wraps the stable C API (`seams.h`) through nanobind. Rings,
cages and the shell symmetry label (`dseams.cages`, `dseams.symmetry`)
go through `libseams_lang`, which ships in the wheel.

## Upgrading from pydseamslib 2.11

2.14 is built from seams-core itself. The import name is `dseams`.
`import pydseams` and `import pydseamslib` still work and re-export the
documented `Frame`, `read`, `from_ase` and the `dseams` functions, but
the 2.11 `pydseams.yoda` / `cyoda` extension names are gone.

## What the wheels contain

The PyPI wheels are CPU builds of the C++ engine with OpenMP on Linux.
The optional Rust backends (linkcell, minimage, readcon-core), chemfiles
trajectory readers, and GPU paths are not in the wheels; install from
conda-forge or build from source with those enabled. See the
[installation guide](https://docs.dseams.info/install).

MIT licensed.
