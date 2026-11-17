# Command line

`seams` reads a trajectory and prints one result. The cell in this example is `input/traj/genice_sI.lammpstrj`, 46 oxygens in a cube of edge 12.0034886 Å.

```bash
seams read input/traj/genice_sI.lammpstrj
seams cn input/traj/genice_sI.lammpstrj --cutoff 3.5
seams steinhardt input/traj/genice_sI.lammpstrj --cutoff 3.5 --l 6
seams chill-plus input/traj/genice_sI.lammpstrj --cutoff 3.5
seams rings input/traj/genice_sI.lammpstrj --cutoff 3.5
```

`read` prints the atom count and the box. `cn` is the cutoff coordination. `steinhardt` prints the mean local q6. `chill-plus` prints the label census. `rings` prints the primitive-ring histogram. On this cell the five commands print:

```
nop 46 frame 1 box 12.0035 12.0035 12.0035
# site-site
# types 1 1 cutoff 3.5 cn 4.08889 nI 46 nJ 46 volume 1729.51
nop 46 q6 0.595133 q6bar 0.311125
nop 46 clathrate 46
nop 46 rings 102
5 48
6 54
```

`docs/examples/minimal_cli.sh` runs those five commands and checks for `nop 46`, `nI 46`, `q6 0.595`, `clathrate 46`, and `6 54`.

```{literalinclude} ../../../docs/examples/minimal_cli.sh
:language: bash
```
