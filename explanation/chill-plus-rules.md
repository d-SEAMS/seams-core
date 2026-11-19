<!-- https://docs.dseams.info/explanation/chill-plus-rules.md -->

# CHILL+: the published rule and the OVITO rule {#chill-the-published-rule-and-the-ovito-rule}

Nguyen and Molinero, *J. Phys. Chem. B* 2015, 119, 9369
([10.1021/jp510289t](https://doi.org/10.1021/jp510289t)), equation 1
averages $q_3$ over the four nearest neighbours of a shell that
already holds four neighbours inside the cutoff. Table 1 gives the
bond rows. A staggered bond is $c(i,j) \le -0.8$. An eclipsed bond
lies on $-0.35 \le c(i,j) \le 0.25$. Those edges are closed.
Footnote a of Table 1 is the neighbour test: two staggered bonds need
a neighbour with more than two staggered bonds, and three staggered
bonds with no eclipsed bond need a neighbour with more than one.

The OVITO 3.17.2 Chill+ modifier
([Chill+](https://docs.ovito.org/reference/pipelines/modifiers/chill_plus.md))
averages $q_3$ over the whole cutoff shell. Its staggered edge is
the open cut $c(i,j) < -0.8$, and the eclipsed band is open at both
ends. The modifier skips the neighbour test in footnote a.
`seams ovito-chill+` and `seams chill-plus --rule ovito-chill+`
implement that rule. `seams chill-plus` with no `--rule` is the 2015
rule.

## One oxygen from the liquid {#one-oxygen-from-the-liquid}

`input/chill-def/liquid-four-nearest.txt` is an 18-atom piece of the
TIP4P/ice liquid. The focus is atom 1, the second oxygen in the file.
Its cutoff shell holds four neighbours. One of those neighbours has a
fifth oxygen inside 3.5 Å, so the whole-shell average of $q_3$
moves off the four-nearest average. None of the four correlations sits on
`-0.8`, `-0.35` or `0.25`.

The published rule counts two eclipsed bonds and labels atom 1 water.
`ovito-chill+` counts three eclipsed bonds and labels atom 1
interClathrate. That is the one atom the fixture is built to move.

```text
fixture focus 1 published water ovito interClathrate
```

## The 4096-molecule liquid {#the-4096-molecule-liquid}

On the 4096-molecule TIP4P/ice liquid
(`tests/data/tip4pice-liquid-4096.txt`), 250 of 4096 oxygens differ
between the published rule and OVITO 3.17.2. Every one of the 250 has
four neighbours inside 3.5 Å. 235 of them change an eclipsed or
staggered count because a neighbour has a fifth oxygen in the shell:
that is the $q_3$ average. The other 15 keep those bond counts and
differ on the neighbour test alone. `ovito-chill+` matches the OVITO
modifier on all 4096 atoms.

Ideal cubic ice stays cubic under both rules. The bundled cell
`input/traj/mW_cubic.lammpstrj` is 4096 cubic molecules under each
rule. Ideal ice Ih (`input/traj/ice1h.lammpstrj`) stays hexagonal
under both rules, 128 of 128.

## Names of the two rules {#names-of-the-two-rules}

`seams chill-plus` takes `--rule chill+` or `--rule ovito-chill+`.
The command `seams ovito-chill+` is the second of those. On the ice
Ih frame each of these prints `nop 128 hexagonal 128`.

```bash
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1 --rule chill+
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1 --rule ovito-chill+
seams ovito-chill+ input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

The C call `seams_chill_plus` takes `SEAMS_CHILL_PUBLISHED` or
`SEAMS_CHILL_OVITO`. A type argument of 0 keeps every atom.
`seams_frame_chill_plus_rule` takes the strings `chill+` and
`ovito-chill+`. A positive type keeps that LAMMPS species. On ice Ih
the atoms are type 1, and both calls label 128 hexagonal molecules.

```c
seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts);
seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, labels, counts);
seams_frame_chill_plus_rule(frame, 3.5, 1, "cutoff", "chill+", labels, counts);
seams_frame_chill_plus_rule(frame, 3.5, 1, "cutoff", "ovito-chill+", labels, counts);
```

Python, in pydseamslib, is `dseams.chill_plus`. `dseams.PUBLISHED` is
the 2015 rule. `dseams.OVITO` is the OVITO rule. The eight-molecule
conventional ice Ic cell below (oxygen spacing 2.75 Å) is cubic under
both rules. The first population count is cubic.

```python
import math
import numpy as np
import dseams

a = 2.75 * 4.0 / math.sqrt(3.0)
fcc = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
basis = np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])
positions = ((fcc[:, None, :] + basis[None, :, :]) * a).reshape(-1, 3)
box = np.array([0.0, 0.0, 0.0, a, a, a, 0.0, 0.0, 0.0])
labels, counts = dseams.chill_plus(positions, box, rule=dseams.PUBLISHED)
assert len(labels) == 8 and int(counts[0]) == 8 and int(counts.sum()) == 8
labels, counts = dseams.chill_plus(positions, box, rule=dseams.OVITO)
assert len(labels) == 8 and int(counts[0]) == 8 and int(counts.sum()) == 8
```
