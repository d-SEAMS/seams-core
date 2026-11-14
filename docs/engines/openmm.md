# OpenMM

There is one reporter, `DseamsReporter` in
`adapters/python/dseams_insitu/openmm.py`.
`DseamsReporter(engine, dump, interval=1, types=None)` asks
`describeNextReport` for positions every `interval` steps. `dump` is the
trajectory path the reporter was constructed with. `types` selects which
atoms are oxygens. `report` reads the context in angstrom and calls
`InSitu.observe`.

`meson test openmm_reporter` runs `engines/openmm/test_reporter.py`.
`meson test docs_engine_openmm` runs `DseamsReporter` on
`input/traj/ice1h.lammpstrj`. `meson test openmm_plugin` runs
`engines/openmm/plugin.cpp`, which calls
`seams_chill_plus` after scaling OpenMM nanometres by 10.

```{literalinclude} ../../../engines/openmm/test_reporter.py
:language: python
:lines: 12-14,58-61
```

```{literalinclude} ../../../engines/openmm/plugin.cpp
:language: cpp
:lines: 46-56
```

`engines/openmm/test_reporter.py` builds an 8-atom cubic cell, a zero
external force, and the Reference platform, then steps twice. Every label
is cubic.

`plugin.cpp` is the C++ counterpart: an OpenMM `Context` on the Reference
platform takes the same two steps and calls `seams.h` directly. Build it when
`OpenMM.h` and `libOpenMM` are on the link line. The Python reporter is the
test that meson runs as `openmm_reporter`.
