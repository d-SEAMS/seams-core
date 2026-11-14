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

`MODE=q6` computes only the mean local q6 and `MODE=q6bar` only the
mean Lechner-Dellago q6-bar, which has no derivatives. Neither runs CHILL+
or the ring search, and neither has the other four components. `DUMP=file`
appends one LAMMPS dump frame per call with the coordinates handed to
d-SEAMS, in angstroms, the per-atom q6 or q6-bar, and under `MODE=all` the
CHILL+ label. `seams chill-plus` reads that dump back.

`meson test plumed_driver` runs `engines/plumed/plumed.dat`, and
`meson test docs_engine_plumed` runs the same `SEAMS` line on
`input/traj/ice1h.lammpstrj`. `meson test plumed_e2e` runs LAMMPS
`fix plumed` on mW ice and liquid and compares every dumped label with
the CLI. It loads the `libplumedKernel.so` named by `PLUMED_KERNEL`, or the
one beside `plumed` on `PATH` when that is unset. The collective variable calls `seams_q4q6`, `seams_chill_plus`
and `seams_rings_find` after scaling PLUMED lengths to angstroms. The q6
derivative uses a 0.001 angstrom step and is converted to PLUMED's length
unit.

```{literalinclude} ../../../engines/plumed/plumed.dat
```

```{literalinclude} ../../../engines/plumed/SeamsColvar.cpp
:language: cpp
:start-at: bool SeamsCV::score
:end-before: void SeamsCV::writeDump
```

`getPosition` and `getBox` are in PLUMED's internal length unit, nanometres
unless `UNITS LENGTH` changes it. The plugin reads that unit and converts
both to angstroms before calling d-SEAMS, so `UNITS LENGTH=A` and the
default give the same result; `meson test plumed_units` checks that on
`ice.xyz`. `CUTOFF` is always angstroms. `ice.xyz` is the cubic cell in
nanometres, so the driver reads it in the default MD length unit. Load
`libSeamsPlumed.so` with `LOAD` before `SEAMS`.
On the ideal cell q6 is high, the ice fraction is 1, and the ring fraction is
positive. The restraint biases `seams.q6` only.
