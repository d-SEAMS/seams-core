# PLUMED tutorial

`SEAMS` has three components. `q6` is the mean local q6 and carries finite-difference derivatives, so `RESTRAINT` can bias it. `ice` is the published CHILL+ cubic fraction. `ring` is the fraction of atoms on a primitive six-ring.

`plumed driver` reads xyz in nanometres. `engines/plumed/ice.xyz` is the 64-atom cell in nanometres. `plumed.dat` sets `UNITS LENGTH=A`, so the collective variable sees Angstrom.

```
LOAD FILE=./libSeamsPlumed.so
UNITS LENGTH=A TIME=ps ENERGY=kj/mol
seams: SEAMS ATOMS=1-64 CUTOFF=3.5
RESTRAINT ARG=seams.q6 AT=0.5 KAPPA=10.0
PRINT ARG=seams.q6,seams.ice,seams.ring FILE=colvar STRIDE=1
```

```bash
c++ -shared -fPIC -O2 -std=c++17 -I"$PLUMED_INCLUDE" -Iinclude -Ibuild/src \
  engines/plumed/SeamsColvar.cpp -Lbuild/src -lseams -lplumed -o libSeamsPlumed.so
plumed driver --ixyz engines/plumed/ice.xyz --plumed engines/plumed/plumed.dat
```

On this cell the colvar line is q6 about 0.63, ice fraction 1, and ring fraction 1.
