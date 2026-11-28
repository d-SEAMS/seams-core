<!-- https://docs.dseams.info/reference/walk-compare.md -->

# walk_compare {#walk-compare}

`tests/walk_compare` walks a LAMMPS dump and prints one row per
frame. Each row has four-neighbour CHILL+ counts from
`chill::getCorrelPlus` and `chill::getIceTypePlusNoPrint` on the
cutoff neighbour list, cage membership on that same graph, and the
seeded assignment on the mutual and union four-nearest graphs, with
the largest ice cluster of each score. Cage flags are hexagonal cage
(HC) and double-diamond cage (DDC).

The binary is the Meson target `tests/walk_compare`. A test build
writes it at `bbdir/tests/walk_compare`.

## Inputs {#inputs}

```text
walk_compare TRAJ [lastFrame] [atomType] [stride]
```

| argument | default | meaning |
| --- | --- | --- |
| TRAJ | required | LAMMPS dump path |
| lastFrame | all frames | inclusive 1-based cap; a non-positive value, or a value at or above the dump length, walks every frame |
| atomType | 1 | LAMMPS type kept by sinp::readLammpsTrjO |
| stride | 1 | frame step; a non-positive value becomes 1 |

Missing `TRAJ` exits `2` with a usage line. A dump with no frames
exits `1`.

Neighbour knobs are hard-wired:

| knob | value |
| --- | --- |
| cutoff | 3.5 |
| k-NN candidate cutoff | 5.5 |
| k | 4 |

`seededCageAffiliation` is called with the completion flag left at
its default, `false`. See
[Ring-adjacent completion](../explanation/ring-completion.md).

## Columns {#columns}

Stdout opens with two comments: the dump path with `nframes`,
`last`, and `stride`, then the column names. Each later line is
twenty integers.

| column | source |
| --- | --- |
| frame | 1-based dump index |
| nop | atoms kept of atomType |
| chill_cubic | CHILL+ cubic (and reCubic) |
| chill_hex | CHILL+ hexagonal (and reHex) |
| chill_interfacial | CHILL+ interfacial |
| chill_clathrate | CHILL+ clathrate |
| chill_interclathrate | CHILL+ interClathrate |
| chill_water | every other CHILL+ label |
| chill_ice | cubic plus hexagonal (Nguyen and Molinero bulk ice; interfacial is not ice) |
| chill_max | largest chill_ice component on the cutoff graph |
| chill_clus | number of those components |
| cut_ice | atoms in an HC- or DDC-affiliated six-ring on the cutoff graph |
| cut_max | largest cut_ice component on the cutoff graph |
| cut_clus | number of those components |
| seed_ih | seeded HC and not DDC |
| seed_ic | seeded DDC and not HC |
| seed_both | seeded HC and DDC |
| seed_ice | seed_ih + seed_ic + seed_both |
| seed_max | largest seed_ice component on the union four-nearest graph |
| seed_clus | number of those components |

Those twenty integers split into star counts, cutoff-graph cage
counts (`cut_*`), and seeded counts (`seed_*`). The star counts and
`cut_*` share the cutoff neighbour list (`nneigh::neighListO` at
`3.5`). The `seed_*` columns use `nneigh::kNearestNeighbourPair`
(`k=4`, candidate cutoff `5.5`): mutual seeds, union completions.
Cluster walks are connected components of the flagged atoms on that
same index graph. A frame with `nop 0` prints the frame index and
nineteen zeros.

## Example {#example}

```bash
bbdir/tests/walk_compare input/traj/exampleTraj.lammpstrj 1 2
```

That walks frame 1, type 2 (TIP4P oxygen in the supplied mixed dump).
Add a fourth argument to stride a long trajectory.
