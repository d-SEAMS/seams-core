<!-- https://docs.dseams.info/explanation/stacking-disorder.md -->

# Stacking disorder {#stacking-disorder}

Cubic ice and hexagonal ice are the same puckered hexagons stacked in
two orders. Hexagonal ice repeats every two layers. Cubic ice repeats
every three, and the next layer sits over the hole of the one before.
A stacking fault is one layer in the other order. Both crystals, and
the fault between them, have a high local $q_6$. A cutoff on $q_6$
therefore reports ice and stops.

Two counts in this tree keep the polytypes apart.

The cage census names the cages. A hexagonal layer is a six-ring that
is the base of a hexagonal cage: two basal faces and three prismatic
faces. A cubic layer is a six-ring that closes a double-diamond cage.
A molecule that belongs to both cages is `both`. The rest of each cage
is `hexagonal` or `cubic`. None of those names is a threshold on
$q_6$. `seams cages` prints the four counts.

```bash
seams cages input/traj/interface_f0.lammpstrj --type 1
```

On that 864-molecule slab the seeded graph and the cutoff graph both
print hexagonal 444, cubic 0, `both` 0 and water 420. The mutual graph
prints hexagonal 442, cubic 0, `both` 0 and water 422. The one cubic molecule in the slab is not a
double-diamond cage, and no molecule sits in both cages.

CHILL+ names hexagonal and cubic from the bond correlations, and
`--layers` writes the stack along the long edge. The published rule is
the 2015 table. The account of that rule is
[CHILL+: the published rule and the OVITO rule](chill-plus-rules.md).

```bash
seams chill-plus input/traj/interface_f0.lammpstrj --type 1 --layers --axis 2
```

The frame is the one in
[the interface tutorial](../user/tutorials/interface.md).
The command prints hexagonal 352, cubic 1 and water 434, and the stack
`HHHHHHHH....C.`: eight hexagonal slabs, a gap, and one cubic layer.
A $q_6$ cutoff would add the 352 and the 1. The stack string keeps
them apart, which is what a growth rate on the basal face needs. The
front advances by a layer of one polytype. The fault is that layer.

The same split is why a nucleation frame reports cubic and hexagonal
counts separately. On `input/traj/nucleation_f1.lammpstrj` the published
rule prints hexagonal 62 and cubic 47. The largest crystal cluster can
contain both, and the fraction of each says which polytype the nucleus
is growing. $q_6$ on that cluster is one number.
