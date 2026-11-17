# OpenMM tutorial

`DseamsReporter` in `adapters/python/dseams_insitu/openmm.py` asks for positions every interval and calls `InSitu.observe`.

```{literalinclude} ../../../engines/openmm/test_reporter.py
:language: python
:lines: 12-14,58-61
```

```bash
PYTHONPATH=adapters/python python engines/openmm/test_reporter.py
```

The script builds an eight-atom cubic cell, a zero external force, and the Reference platform. A passing run prints `openmm calls 2 ice 16`.

`engines/openmm/plugin.cpp` is the same two steps in C++, linked to `libOpenMM` and `libseams`:

```bash
c++ -O2 -std=c++17 -I"$OPENMM_INCLUDE" -Iinclude -Ibuild/src engines/openmm/plugin.cpp \
  -L"$OPENMM_LIB" -lOpenMM -Lbuild/src -lseams -o openmm-plugin
```
