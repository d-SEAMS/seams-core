<!-- https://docs.dseams.info/user/tutorials/paper-rings.md -->

# Primitive rings and the shortest-path test {#primitive-rings-and-the-shortest-path-test}

A ring search first walks every simple cycle up to a chosen size.
Many of those cycles are not useful: a chord gives a shorter route
between two vertices, so the long way around is not a ring of the
network. The Franzblau test keeps a cycle only when, for every pair
of its vertices, some shortest path between them lies on the cycle.
Those survivors are the primitive rings. `seams rings` reports that
census. The example also counts the cycles before the test, so you can
see what the filter removes.

The cell is GenIce ice Ih, one unit, 16 oxygens, cutoff 3.5 Å.
It is small on purpose. Every cycle of length 4, 5 or 6 can be
listed, and the same neighbour list is what the primitive search uses.
A larger crystal would hide the filter, because almost every six-cycle
in ice Ih already passes.

The ring search gives

```text
rings primitive 4 24 6 32
rings shortest 4 24 6 32
rings searched 4 24 6 128
```

The shortest-path test and `seams rings` agree: 24 four-rings and 32
six-rings. The unfiltered search found 128 six-cycles and no
five-ring. The other 96 six-cycles have a shortcut, so they are
dropped.

## Full run {#full-run}

The paper ring census is one frame of
`input/traj/mW_cubic.lammpstrj`. The dump holds 11 frames of 4096
molecules. Frame 1 has 4096 molecules and 8192 six-rings.

```text
rings full input molecules 4096 frames 1 file-frames 11 file input/traj/mW_cubic.lammpstrj
rings full paper molecules 4096 frames 1 rings 8192
rings full nop 4096 rings 8192
rings full histogram 6 8192
```

The histogram has one bin: 8192 six-rings.
