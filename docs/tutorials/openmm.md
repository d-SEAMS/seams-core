# OpenMM tutorial

`DseamsReporter` asks for positions every interval and calls the C API.

```python
from reporter import DseamsReporter
simulation.reporters.append(DseamsReporter(1))
simulation.step(2)
print(reporter.calls, reporter.ice, reporter.overhead())
```

```bash
PYTHONPATH=bindings/python python engines/openmm/test_reporter.py
```

The script builds an eight-atom cubic cell, a zero external force, and the Reference platform. A passing run prints `openmm calls 2 ice 16`.

`engines/openmm/plugin.cpp` is the same two steps in C++, linked to `libOpenMM` and `libseams`:

```bash
c++ -O2 -std=c++17 -I"$OPENMM_INCLUDE" -Iinclude -Ibuild/src engines/openmm/plugin.cpp \
  -L"$OPENMM_LIB" -lOpenMM -Lbuild/src -lseams -o openmm-plugin
```
