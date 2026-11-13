# PLUMED

The thread count comes from `PLUMED_NUM_THREADS`. Set that variable
before PLUMED starts. `OMP_NUM_THREADS` does not change the PLUMED team.

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

`meson test plumed_driver` runs `engines/plumed/plumed.dat`. The collective
variable calls `seams_q4q6` and `seams_rings_find` after scaling internal
nanometres to angstroms. The q6 derivative uses a 0.001 angstrom step and
multiplies by 10 so PLUMED sees a per-nanometre derivative.

```{literalinclude} ../../../engines/plumed/plumed.dat
```

```{literalinclude} ../../../engines/plumed/SeamsColvar.cpp
:language: cpp
:lines: 73-86
```

PLUMED stores positions in nanometres. `getPosition` and `getBox` are those
internal nanometres, and the plugin multiplies both by 10 before calling
d-SEAMS. `CUTOFF` is angstroms. `ice.xyz` is the cubic cell in nanometres, so
the driver reads it in the default length unit. Do not set `UNITS LENGTH=A`
on this deck: that flag does not change `getPosition`, and it would disagree
with the file. Load `libSeamsPlumed.so` with `LOAD` before `SEAMS`.
On the ideal cell q6 is high, the ice fraction is 1, and the ring fraction is
positive. The restraint biases `seams.q6` only.
