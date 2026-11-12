# PLUMED tutorial

`SEAMS` has three components. `q6` is the mean local q6 and carries finite-difference derivatives, so `RESTRAINT` can bias it. `ice` is the published CHILL+ cubic fraction and `ring` is the fraction of atoms on a primitive six-ring. Those two components have no derivatives. They are not biasable.

`plumed driver` reads xyz in nanometres. `engines/plumed/ice.xyz` is the 64-atom cell in nanometres. The plugin multiplies `getPosition` and the box by 10, because d-SEAMS expects angstroms and `CUTOFF=3.5` is angstroms. Leave the driver's length unit at the default. Do not set `UNITS LENGTH=A`.

```
LOAD FILE=./libSeamsPlumed.so
seams: SEAMS ATOMS=1-64 CUTOFF=3.5
RESTRAINT ARG=seams.q6 AT=0.5 KAPPA=10.0
PRINT ARG=seams.q6,seams.ice,seams.ring FILE=colvar STRIDE=1
```

```bash
c++ -shared -fPIC -O2 -std=c++17 -I"$PLUMED_INCLUDE" -Iinclude -Ibuild/src \
  engines/plumed/SeamsColvar.cpp -Lbuild/src -lseams -lplumed -o libSeamsPlumed.so
plumed driver --ixyz engines/plumed/ice.xyz --plumed engines/plumed/plumed.dat
```

`$PLUMED_INCLUDE` is the directory that contains `plumed/core/Colvar.h`. On this cell the colvar line is q6 about 0.63, ice fraction 1, and ring fraction 1. A restraint may name `seams.q6` only.
