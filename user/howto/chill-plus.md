<!-- https://docs.dseams.info/user/howto/chill-plus.md -->

# CHILL+ {#chill}

Label every oxygen with the published CHILL+ rule, or with the OVITO
rule. Cage membership is a separate score. Why the two disagree is in
[CHILL+ versus cages](../../explanation/chill-plus-vs-cages.md).

## The default list {#the-default-list}

`seams chill-plus` bonds every pair closer than `--cutoff` (default
3.5 Å). The stored name of that default graph is `seeded`. CHILL and
CHILL+ map `seeded` onto the cutoff list. `--graph mutual` keeps the
mutual four nearest neighbours inside `--cutoff` + 1.5 Å. `seams cages` leaves `seeded` as its own two-list graph.

```bash
# excerpt: docs/orgmode/user/tutorials/bulk-ice.org
seams chill-plus input/traj/exampleTraj.lammpstrj --cutoff 3.5 --type 2
```

```bash
# excerpt: docs/orgmode/user/tutorials/bulk-ice.org
seams chill-plus input/traj/exampleTraj.lammpstrj --cutoff 3.5 --type 2 --graph mutual
```

The cutoff list prints:

```default
nop 250 interClathrate 9 water 241
```

The mutual list prints:

```default
nop 250 clathrate 1 interClathrate 7 water 242
```

A count of zero is omitted. The names are the CHILL+ labels, including
`interClathrate`.

`--water-types` builds the neighbour graph on those LAMMPS types,
including pairs between them, and leaves every other atom
unclassified. A multi-atom dump with no `--type` infers the oxygen of
each molecule.

## The OVITO rule {#the-ovito-rule}

`--rule ovito-chill+` averages $q_3$ on the full cutoff shell and
skips the published neighbour test. `seams ovito-chill+` is that rule
as its own command. The C constants are `SEAMS_CHILL_PUBLISHED` and
`SEAMS_CHILL_OVITO`. The 18-atom liquid fixture, and the 250 atoms
that change label on the 4096-molecule liquid, are on
[the published rule and the OVITO rule](../../explanation/chill-plus-rules.md).

## From the library {#from-the-library}

`seams_chill_plus` writes one label per atom and nine counts, one per
CHILL+ label.

```c
seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts);
```

On the 64-atom cubic cell this published cutoff call labels every atom
`SEAMS_LABEL_CUBIC`. The OVITO rule on the mutual list does the same.
