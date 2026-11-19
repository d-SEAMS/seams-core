<!-- https://docs.dseams.info/user/howto/cages-trace.md -->

# Cages and TRACE {#cages-and-trace}

Ask which molecules sit in a complete ice cage, or which polyhedra
match a ring-size signature. TRACE, as run by ClavisCage, searches a
hydrogen-bond graph for hydrate polyhedra. The commands below are the
d-SEAMS counts.

## Ice cages {#ice-cages}

`seams cages` reports hexagonal cages (HC, ice Ih), double-diamond
cages (DDC, ice Ic), molecules that sit in both, and water. The
default graph is `seeded`: the mutual four-nearest list and the union
of each atom’s four nearest neighbours, both taken inside
`--cutoff` + 1.5 Å. The cutoff default is 3.5 Å. Water is every
selected molecule that is in neither complete cage.

```bash
seams cages input/traj/exampleTraj.lammpstrj --type 2
seams cages input/traj/mW_cubic.lammpstrj --type 1
```

The mixed TIP4P frame, oxygen type 2, prints:

```default
nop 250 graph seeded hexagonal 0 cubic 0 both 0 water 250
```

The 4096-molecule cubic mW cell prints:

```default
nop 4096 graph seeded hexagonal 0 cubic 4096 both 0 water 0
```

## A signature {#a-signature}

`--signature` names a face census. `512` is twelve pentagons and
prints `5:12`. `51262` is twelve pentagons and two hexagons and prints
`5:12,6:2`. Other names are `51264`, `51268`, `hc`, `ddc`, `sh`,
`sodalite` and `alpha`. The default face-edge window is `2.65,3.15`
in the coordinate units of the file. `--cage-edge none` keeps every
face from the ring census. `--guest-types` names the LAMMPS types of
guests. Each guest is assigned to the nearest cage center inside
`--guest-radius` (default 4.0). The language call is
`seams_frame_signature`. A window is `seams_frame_signature_edge`.

On the default graph name `seeded`, a signature search uses the union
list. The command below reads the structure I cell on the cutoff
graph, the oxygen list used for ice and for the TIP4P oxygens in the
TRACE comparison. Joining names with `+` prints one line per
signature. The cell is orthogonal, and `--cage-edge none` on this
command prints the same two lines.

```bash
# excerpt: docs/orgmode/user/tutorials/clathrate.org
seams cages input/traj/genice_sI.lammpstrj --type 1 --graph cutoff --signature 512+51262
```

```default
nop 46 graph cutoff signature 5:12 cages 2 atoms 40
nop 46 graph cutoff signature 5:12,6:2 cages 6 atoms 46
```

## TRACE {#trace}

ClavisCage 0.3.18 runs TRACE, the hydrogen-bond cage search of Hsu
and Lin, J. Chem. Theory Comput. 2025
([10.1021/acs.jctc.5c01459](https://doi.org/10.1021/acs.jctc.5c01459)).
The two constructions are
[TRACE and signature cages](../../explanation/trace-assembly.md).
d-SEAMS counts $5^{12}$, $5^{12}6^{2}$ and $5^{12}6^{4}$ with
`seams cages --signature` and `seams_frame_signature`. It does not
call TRACE.

The comparison places hydrogens 0.96 Å from each oxygen along an
out-degree-2 orientation of the oxygen neighbour graph. TRACE then
uses its default 0.35 nm cutoff and a 35 degree angle. On ice and on
the TIP4P oxygens that graph is the 3.5 Å list. The dense frame
`input/traj/dump-l.lammpstrj` (512 atoms of type 1) has median degree
56 inside 3.5 Å. TRACE at 0.35 nm does not return on that shell, so
that frame is compared on the four nearest neighbours.

Ideal ice has no five-rings. `seams rings` on
`input/traj/ice1h.lammpstrj` and on the cubic mW cell prints only
six-rings, so those three cages are absent. The clathrate cells,
where the cages are present, are the table on
[ClavisCage TRACE](../../explanation/clavis-trace.md). The ice
rows are in
[Literature head-to-head](literature-h2h.md).
