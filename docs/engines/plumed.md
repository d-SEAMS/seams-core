# PLUMED

`SEAMS` is a collective variable with three components.

- `seams.q6` is the mean local q6. Its derivatives are a central difference, so
  `RESTRAINT` can bias it.
- `seams.ice` is the fraction of atoms published CHILL+ calls cubic.
- `seams.ring` is the fraction of atoms on a primitive six-ring.

```
UNITS LENGTH=A
seams: SEAMS ATOMS=1-64 CUTOFF=3.5
RESTRAINT ARG=seams.q6 AT=0.5 KAPPA=10.0
PRINT ARG=seams.q6,seams.ice,seams.ring FILE=colvar
```

`plumed driver` reads an xyz file as nanometres. `ice.xyz` is the cubic cell
in nanometres, and `plumed.dat` sets `UNITS LENGTH=A`, so the collective
variable sees Angstrom. Load `libSeamsPlumed.so` with `LOAD` before `SEAMS`.
On the ideal cell q6 is high, the ice fraction is 1, and the ring fraction is
positive. The restraint biases `seams.q6`.
