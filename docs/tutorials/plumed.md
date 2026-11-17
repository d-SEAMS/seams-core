# PLUMED tutorial

`SEAMS` reports three components on this page. `q6` is the mean local q6 and carries finite-difference derivatives, so `RESTRAINT` can bias it. `ice` is the published CHILL+ cubic fraction and `ring` is the fraction of atoms on a primitive six-ring. Those two components have no derivatives.

```{eval-rst}
.. include:: ../videos/embeds/plumed.rst
```

`plumed driver` reads xyz in nanometers. `engines/plumed/ice.xyz` is the 64-atom cubic cell in nanometers. The plugin converts `getPosition` and the box from PLUMED's length unit to angstroms, because d-SEAMS expects angstroms and `CUTOFF=3.5` is angstroms. Leave the driver's length unit at the default. `UNITS LENGTH=A` in the deck gives the same `seams.q6`, `seams.ice` and `seams.ring`.

The shared library is `SeamsColvar.cpp` plus `plumed_components.cpp`. `plumed.dat` loads `./libSeamsPlumed.so`.

```{literalinclude} ../../../engines/plumed/plumed.dat
```

```bash
c++ -shared -fPIC -O2 -std=c++17 -I"$PLUMED_INCLUDE" -Iinclude -Ibuild/src \
  engines/plumed/SeamsColvar.cpp engines/plumed/plumed_components.cpp \
  -Lbuild/src -lseams -lplumed -o libSeamsPlumed.so
plumed driver --ixyz engines/plumed/ice.xyz --plumed engines/plumed/plumed.dat
```

`$PLUMED_INCLUDE` is the directory that contains `plumed/core/Colvar.h`. On this cell `engines/plumed/run_driver.py` requires `seams.ice` of 1. A restraint may name `seams.q6` only.
