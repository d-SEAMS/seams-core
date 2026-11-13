# Command line

`seams` reads a trajectory and prints one result. The small cell used by
the runnable example is `input/traj/genice_sI.lammpstrj`, 46 oxygens in a
cube of edge 12.0034886 Å.

```bash
seams read input/traj/genice_sI.lammpstrj
seams cn input/traj/genice_sI.lammpstrj --cutoff 3.5
seams steinhardt input/traj/genice_sI.lammpstrj --cutoff 3.5 --l 6
seams chill-plus input/traj/genice_sI.lammpstrj --cutoff 3.5
seams rings input/traj/genice_sI.lammpstrj --cutoff 3.5
```

`read` prints the atom count and the box. `cn` is the cutoff coordination.
`steinhardt` prints the mean local q6. `chill-plus` prints the label
census. `rings` prints the primitive-ring histogram. On this cell those
lines are `nop 46`, coordination over 46 oxygens, `q6 0.595`, `clathrate 46`,
and `6 54`.

`docs/examples/minimal_cli.sh` runs those five
commands. `meson test docs_examples` runs it.
Lua is not a binding in this repository. `examples/chill-plus.lua` is an
options table for the separate yodaStruct package.

```{literalinclude} ../../../docs/examples/minimal_cli.sh
:language: bash
```
