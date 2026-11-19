<!-- https://docs.dseams.info/fronts/plumed.md -->

# PLUMED tutorial {#plumed-tutorial}

Use the [PLUMED engine page](../engines/plumed.md) to get the d-SEAMS
PLUMED plugin v0.4.1.

`SEAMS` with `MODE=all` reports five components. `q6` is the mean local q6 and carries finite-difference derivatives, so `RESTRAINT` can bias it. `ice` is the published CHILL+ cubic fraction. `ring` is the fraction of atoms on a primitive six-ring. `hexagonal` is the CHILL+ hexagonal count. `six` is the primitive six-ring count. The last four have no derivatives, and a `RESTRAINT` on them is rejected. `MODE=q6` and `MODE=q6bar` report only `q6`. The standalone PLUMED action (`dseams-plumed` v0.4.1) is indexed in [PLUMED-NEST egg 26.018](https://www.plumed-nest.org/eggs/26/018/).

`plumed driver` reads xyz in nanometers. `engines/plumed/ice.xyz` is the 64-atom cubic cell in nanometers. The plugin converts `getPosition` and the box from PLUMED’s length unit to angstroms, because d-SEAMS expects angstroms and `CUTOFF=3.5` is angstroms. Leave the driver’s length unit at the default. `UNITS LENGTH=A` in the deck gives the same `seams.q6`, `seams.ice` and `seams.ring`.

The shared library is `SeamsColvar.cpp` plus `plumed_components.cpp`. `plumed.dat` loads `./libSeamsPlumed.so`.

```default
LOAD FILE=./libSeamsPlumed.so
seams: SEAMS ATOMS=1-64 CUTOFF=3.5
RESTRAINT ARG=seams.q6 AT=0.5 KAPPA=10.0
PRINT ARG=seams.q6,seams.ice,seams.ring FILE=colvar STRIDE=1
```

```bash
c++ -shared -fPIC -O2 -std=c++17 -I"$PLUMED_INCLUDE" -Iinclude -Ibuild/src \
  engines/plumed/SeamsColvar.cpp engines/plumed/plumed_components.cpp \
  -Lbuild/src -lseams -lplumed -o libSeamsPlumed.so
plumed driver --ixyz engines/plumed/ice.xyz --plumed engines/plumed/plumed.dat
```

`$PLUMED_INCLUDE` is the directory that contains `plumed/core/Colvar.h`. On this cell `engines/plumed/run_driver.py` requires `seams.ice` of 1. A restraint may name `seams.q6` only.
