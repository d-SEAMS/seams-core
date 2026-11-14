# PLUMED tutorial

`SEAMS` has three components. `q6` is the mean local q6 and carries finite-difference derivatives, so `RESTRAINT` can bias it. `ice` is the published CHILL+ cubic fraction and `ring` is the fraction of atoms on a primitive six-ring. Those two components have no derivatives. They are not biasable.

`plumed driver` reads xyz in nanometres. `engines/plumed/ice.xyz` is the 64-atom cell in nanometres. The plugin converts `getPosition` and the box from PLUMED's length unit to angstroms, because d-SEAMS expects angstroms and `CUTOFF=3.5` is angstroms. Leave the driver's length unit at the default; `UNITS LENGTH=A` in the deck gives the same result.

```{literalinclude} ../../../engines/plumed/plumed.dat
```

```bash
c++ -shared -fPIC -O2 -std=c++17 -I"$PLUMED_INCLUDE" -Iinclude -Ibuild/src \
  engines/plumed/SeamsColvar.cpp -Lbuild/src -lseams -lplumed -o libSeamsPlumed.so
plumed driver --ixyz engines/plumed/ice.xyz --plumed engines/plumed/plumed.dat
```

`$PLUMED_INCLUDE` is the directory that contains `plumed/core/Colvar.h`. On this cell the colvar line is q6 about 0.63, ice fraction 1, and ring fraction 1. A restraint may name `seams.q6` only.
