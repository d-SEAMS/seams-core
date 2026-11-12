# OpenMM

`engines/openmm/reporter.py` is a reporter. `describeNextReport` asks for
positions every `reportInterval` steps. `report` converts the state to
Angstrom and calls `dseams.chill_plus`.

```python
from reporter import DseamsReporter
simulation.reporters.append(DseamsReporter(1))
```

`test_reporter.py` builds an 8-atom cubic cell, a zero external force, and the
Reference platform, then steps twice. It prints `overhead_s`.

`plugin.cpp` is the C++ counterpart: an OpenMM `Context` on the Reference
platform takes the same two steps and calls `seams.h` directly. Build it when
`OpenMM.h` and `libOpenMM` are on the link line. The Python reporter is the
test that always runs.
