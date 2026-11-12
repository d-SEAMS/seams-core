# OpenMM

There is one reporter, `DseamsReporter` in
`adapters/python/dseams_insitu/openmm.py`. `describeNextReport` asks for
positions every `interval` steps. `report` reads the context in angstrom
and calls `InSitu.observe`.

```python
from dseams_insitu.engine import InSitu
from dseams_insitu.openmm import DseamsReporter

engine = InSitu()
simulation.reporters.append(DseamsReporter(engine, "unused", interval=1))
```

`engines/openmm/test_reporter.py` builds an 8-atom cubic cell, a zero
external force, and the Reference platform, then steps twice. Every label
is cubic.

`plugin.cpp` is the C++ counterpart: an OpenMM `Context` on the Reference
platform takes the same two steps and calls `seams.h` directly. Build it when
`OpenMM.h` and `libOpenMM` are on the link line. The Python reporter is the
test that meson runs as `openmm_reporter`.
