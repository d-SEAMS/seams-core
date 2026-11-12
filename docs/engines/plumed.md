# PLUMED

`SEAMS` is a collective variable with three components.

- `seams.q6` is the mean local q6. Its derivatives are a central difference, so
  `RESTRAINT` can bias it.
- `seams.ice` is the fraction of atoms published CHILL+ calls cubic. The label
  is discrete, so the component has no derivatives and is not biasable.
- `seams.ring` is the fraction of atoms on a primitive six-ring. It is not
  biasable either: a `RESTRAINT` on `seams.ice` or `seams.ring` is rejected
  because those components carry no derivatives.

```
seams: SEAMS ATOMS=1-64 CUTOFF=3.5
RESTRAINT ARG=seams.q6 AT=0.5 KAPPA=10.0
PRINT ARG=seams.q6,seams.ice,seams.ring FILE=colvar
```

PLUMED stores positions in nanometres. `getPosition` and `getBox` are those
internal nanometres, and the plugin multiplies both by 10 before calling
d-SEAMS. `CUTOFF` is angstroms. `ice.xyz` is the cubic cell in nanometres, so
the driver reads it in the default length unit. Do not set `UNITS LENGTH=A`
on this deck: that flag does not change `getPosition`, and it would disagree
with the file. Load `libSeamsPlumed.so` with `LOAD` before `SEAMS`.
On the ideal cell q6 is high, the ice fraction is 1, and the ring fraction is
positive. The restraint biases `seams.q6` only.
