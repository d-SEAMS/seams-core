<!-- https://docs.dseams.info/engines/openmm.md -->

# OpenMM {#openmm}

There is one reporter, `DseamsReporter` in
`adapters/python/dseams_insitu/openmm.py`.
`DseamsReporter(engine, dump=None, interval=1, types=None)` asks
`describeNextReport` for positions every `interval` steps. `types` selects
which atoms are oxygens. `report` reads the context in angstrom and calls
`InSitu.observe`. With a `dump` path each report appends a LAMMPS dump
frame of the classified sites, keyed by 1-based particle index, with their
CHILL+ labels. `seconds` and `calls` time the reports.

```python
from dseams_insitu.engine import InSitu
from dseams_insitu.openmm import DseamsReporter

    engine = InSitu()
    reporter = DseamsReporter(engine, interval=1)
    simulation.reporters.append(reporter)
    simulation.step(2)
```

```cpp
    const auto &pos = state.getPositions();
    std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
    for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i) {
      xyz[3 * i] = pos[i][0] * 10.0;
      xyz[3 * i + 1] = pos[i][1] * 10.0;
      xyz[3 * i + 2] = pos[i][2] * 10.0;
    }
    const double box9[9] = {0, 0, 0, a, a, a, 0, 0, 0};
    std::vector<int> labels(static_cast<std::size_t>(n));
    int counts[9] = {};
    seams_system *sys = seams_system_create();
    seams_system_set_frame(sys, xyz.data(), n, box9, nullptr, nullptr);
    seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels.data(), counts);
    seams_system_destroy(sys);
```

`engines/openmm/test_reporter.py` builds an 8-atom cubic cell, a zero
external force, and the Reference platform, then steps twice. Every label
is cubic.

`plugin.cpp` is the C++ counterpart: an OpenMM `Context` on the Reference
platform takes the same two steps and calls `seams.h` directly after scaling
OpenMM nanometers by 10. Build it when `OpenMM.h` and `libOpenMM` are on the
link line.
