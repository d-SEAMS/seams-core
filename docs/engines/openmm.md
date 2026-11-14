# OpenMM

There is one reporter, `DseamsReporter` in
`adapters/python/dseams_insitu/openmm.py`. `describeNextReport` asks for
positions every `interval` steps. `report` reads the context in angstrom
and calls `InSitu.observe`.

`meson test openmm_reporter` runs `engines/openmm/test_reporter.py`.
`meson test openmm_plugin` runs `engines/openmm/plugin.cpp`, which calls
`seams_chill_plus` after scaling OpenMM nanometres by 10.

```{literalinclude} ../../../engines/openmm/test_reporter.py
:language: python
:lines: 12-14,58-61
```

```{literalinclude} ../../../engines/openmm/plugin.cpp
:language: cpp
:start-at: const auto &pos = state.getPositions();
:end-before: ice += counts[0];
```

`engines/openmm/test_reporter.py` builds an 8-atom cubic cell, a zero
external force, and the Reference platform, then steps twice. Every label
is cubic.

`plugin.cpp` is the C++ counterpart: an OpenMM `Context` on the Reference
platform takes the same two steps and calls `seams.h` directly. Build it when
`OpenMM.h` and `libOpenMM` are on the link line. The Python reporter is the
test that meson runs as `openmm_reporter`.
