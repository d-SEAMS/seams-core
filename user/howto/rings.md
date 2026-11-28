<!-- https://docs.dseams.info/user/howto/rings.md -->

# Primitive rings {#primitive-rings}

Count primitive rings on a frame you already have. The definition is
on [Primitive rings](../../explanation/primitive-rings.md). The
learning path is [the rings tutorial](../tutorials/rings.md).

## One frame {#one-frame}

`seams rings` prints `nop`, the number of primitive rings, and one
histogram row for each size that occurs. The neighbour list is the
cutoff list (default 3.5 Å). An orthorhombic cell uses the full CSR
of that list. `--type` keeps one LAMMPS type. On the mixed TIP4P dump
the oxygen is type 2. `--max-ring` is the largest ring (default 6).

```bash
# excerpt: docs/orgmode/getting-started/quickstart.org
seams rings input/traj/exampleTraj.lammpstrj --type 2
```

```default
nop 250 rings 511
3 37
4 273
5 34
6 167
```

## From the C API {#from-the-c-api}

`seams_rings_find` then `seams_rings_histogram` is the same census.
`counts[k]` is the number of primitive rings of size `k`. The buffer
must be longer than the largest ring. On a 64-atom cubic cell the
six-ring count is positive.

```c
seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings);
seams_rings_histogram(rings, counts, 13);
```

## A trajectory {#a-trajectory}

A `--frame` / `--last` walk on one worker (`--jobs 1`, the default)
keeps one ring updater. A bond that did not change keeps its rings. A
changed bond recomputes the sources next to it. The census matches a
full search of that frame. `--rebuild` recomputes every frame. The
measured speedups are on
[Linkcell neighbours and carried updates](linkcell-carried.md).
