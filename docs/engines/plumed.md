# PLUMED

The [d-SEAMS PLUMED plugin egg](https://www.plumed-nest.org/eggs/26/018/)
(plumID:26.018) uses v0.4.1 of the
[dseams-plumed repository](https://github.com/d-SEAMS/dseams-plumed).
It is tested with PLUMED v2.10 and master.

The thread count comes from `PLUMED_NUM_THREADS`. Set that variable
before PLUMED starts. `OMP_NUM_THREADS` does not change the PLUMED team.
The standalone PLUMED action (`dseams-plumed` v0.4.1) is indexed in
[PLUMED-NEST egg 26.018](https://www.plumed-nest.org/eggs/26/018/).

`SEAMS` is a collective variable with five components.

- `seams.q6` is the mean local q6. Its derivatives are a central difference, so
  `RESTRAINT` can bias it. A `PRINT` leaves that derivative array empty, and
  the difference runs only after a `RESTRAINT` has allocated it.
- `seams.ice` is the fraction of atoms published CHILL+ calls cubic. The label
  is discrete, so the component has no derivatives and is not biasable.
- `seams.ring` is the fraction of atoms on a primitive six-ring. It is not
  biasable either: a `RESTRAINT` on `seams.ice` or `seams.ring` is rejected
  because those components carry no derivatives.
- `seams.hexagonal` is the CHILL+ hexagonal count.
- `seams.six` is the primitive six-ring count.

`MODE=q6` computes only the mean local q6 and `MODE=q6bar` only the
mean Lechner-Dellago q6-bar, which has no derivatives. Neither runs CHILL+
or the ring search, and neither has the other four components. `DUMP=file`
appends one LAMMPS dump frame per call with the coordinates handed to
d-SEAMS, in angstroms, the per-atom q6 or q6-bar, and under `MODE=all` the
CHILL+ label. `seams chill-plus` reads that dump back.

```{literalinclude} ../../../engines/plumed/plumed.dat
```

```{literalinclude} ../../../engines/plumed/SeamsColvar.cpp
:language: cpp
:start-at: bool SeamsCV::score
:end-before: void SeamsCV::writeDump
```
