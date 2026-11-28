<!-- https://docs.dseams.info/reference/cli-api.md -->

# Command line {#command-line}

Every `seams` subcommand and flag declared in `src/seams_cli.cpp`. [seams CLI](cli.md) describes the text each subcommand prints. The C column lists `seams_*` calls in the branch that handles that subcommand.

## Native API {#native-api}

Generated from the sources named below. A name on this page is a public call in that language.

### seams –atomname {#seams-atomname}

```bash
--atomname
```

Keep chemfiles atoms with these names

### seams –axis {#seams-axis}

```bash
--axis
```

Histogram axis for density-z (x|y|z, default z)

### seams –bins {#seams-bins}

```bash
--bins
```

Histogram bins (rdf: rmax/0.1; density-z: L/0.1)

### seams –block {#seams-block}

```bash
--block
```

CUDA block size (LINKCELL_BLOCK / SEAMS_BLOCK)

### seams –by-atom {#seams-by-atom}

```bash
--by-atom
```

fingerprint: one row per atom with coordination,

### seams –cage-edge {#seams-cage-edge}

```bash
--cage-edge
```

cages –signature: keep a face only when every edge

### seams –cell {#seams-cell}

```bash
--cell
```

Link-cell hint so NPT frames share a grid

### seams –color-types {#seams-color-types}

```bash
--color-types
```

fingerprint: color vertices by LAMMPS type, so species

### seams –complete {#seams-complete}

```bash
--complete
```

Seeded cages: fill the last vertex of six-rings whose

### seams –config {#seams-config}

```bash
--config
```

Dotenv file of SEAMS_* / LINKCELL_* knobs

### seams –cutoff {#seams-cutoff}

```bash
--cutoff
```

Neighbour cutoff (rdf/cn rmax; hbonds heavy-atom nList)

### seams –donors {#seams-donors}

```bash
--donors
```

Use every hydrogen as a donor candidate

### seams –emit-library {#seams-emit-library}

```bash
--emit-library
```

fingerprint: print the frame’s distinct keys as library lines

### seams –exclude-atomname {#seams-exclude-atomname}

```bash
--exclude-atomname
```

Exclude chemfiles atoms with these names

### seams –exclude-resname {#seams-exclude-resname}

```bash
--exclude-resname
```

Exclude chemfiles residues with these names

### seams –family {#seams-family}

```bash
--family
```

Material family (default waterIce): waterIce,

### seams –features {#seams-features}

```bash
--features
```

Print compile-time backends

### seams –format {#seams-format}

```bash
--format
```

Output format: text (default) or json

### seams –frame {#seams-frame}

```bash
--frame
```

First frame (1-based)

### seams –graph {#seams-graph}

```bash
--graph
```

### seams –guest-radius {#seams-guest-radius}

```bash
--guest-radius
```

cages –signature: a guest belongs to the nearest cage center

### seams –guest-types {#seams-guest-types}

```bash
--guest-types
```

cages –signature: LAMMPS types of guests (methane, THF, ions)

### seams –hangle {#seams-hangle}

```bash
--hangle
```

O-O-H angle cutoff for hbonds (acceptor-centered)

### seams –hdist {#seams-hdist}

```bash
--hdist
```

Acceptor-H distance cutoff for hbonds

### seams –help {#seams-help}

```bash
--help
```

show this help message and exit

### seams –hops {#seams-hops}

```bash
--hops
```

fingerprint: bonds from the center in each local key (default 2)

### seams –htype {#seams-htype}

```bash
--htype
```

Hydrogen atom type for hbonds and f4

### seams –incomplete {#seams-incomplete}

```bash
--incomplete
```

cages –signature: also report cups and incomplete

### seams –inside {#seams-inside}

```bash
--inside
```

cages –guest-types: assign by ray-parity inside the

### seams –ion-cutoff {#seams-ion-cutoff}

```bash
--ion-cutoff
```

ions: first water shell radius (default –cutoff)

### seams –ion-types {#seams-ion-types}

```bash
--ion-types
```

ions: LAMMPS types read against the cage assignment

### seams –ions {#seams-ions}

```bash
--ions
```

cn on ionCloud cage degree (needs –site)

### seams –jobs {#seams-jobs}

```bash
--jobs
```

Parallel frame workers (OpenMP). 1 is serial

### seams –l {#seams-l}

```bash
--l
```

Steinhardt degrees, comma-separated. Default 4,6,

### seams –last {#seams-last}

```bash
--last
```

Last frame (inclusive). Omit for a single –frame

### seams –layers {#seams-layers}

```bash
--layers
```

chill-plus: literature I_sd reference from CHILL+

### seams –library {#seams-library}

```bash
--library
```

fingerprint: key library (from –emit-library) to name atoms

### seams –max-ring {#seams-max-ring}

```bash
--max-ring
```

Largest primitive ring kept by rings (default 6)

### seams –mol {#seams-mol}

```bash
--mol
```

Water oxygens by molecule: water, or ids such as

### seams –pair-types {#seams-pair-types}

```bash
--pair-types
```

Pair types for site-site rdf/cn (default I=J=–type)

### seams –per-atom {#seams-per-atom}

```bash
--per-atom
```

chill|chill-plus|cages|fingerprint|ions: append a LAMMPS dump frame

### seams –print-config {#seams-print-config}

```bash
--print-config
```

Print the resolved runtime knobs and exit

### seams –rebuild {#seams-rebuild}

```bash
--rebuild
```

Recompute rings and cages on every frame. A one-worker

### seams –resident {#seams-resident}

```bash
--resident
```

Fraction of free device memory for a resident batch

### seams –resname {#seams-resname}

```bash
--resname
```

Keep chemfiles residues with these names

### seams –rule {#seams-rule}

```bash
--rule
```

chill-plus rule: chill+ (published default) or

### seams –shape {#seams-shape}

```bash
--shape
```

cages –signature: point group of each found cage,

### seams –signature {#seams-signature}

```bash
--signature
```

cages: ring-size census (4:6,6:8) or a named table

### seams –site {#seams-site}

```bash
--site
```

Type-to-kind map for ionCloud (1=cationHead,2=anion)

### seams –strict-input {#seams-strict-input}

```bash
--strict-input
```

Fail when a requested frame produces no atoms

### seams –subset {#seams-subset}

```bash
--subset
```

domains subset: polar | apolar

### seams –tpp {#seams-tpp}

```bash
--tpp
```

Threads per particle for the device k-NN

### seams –tum-layers {#seams-tum-layers}

```bash
--tum-layers
```

chill-plus: TUM stacking from HC-basal and

### seams –type {#seams-type}

```bash
--type
```

Atom type (0 guesses oxygen then type 1;

### seams –type-map {#seams-type-map}

```bash
--type-map
```

Override chemfiles atom names with positive type IDs

### seams –types {#seams-types}

```bash
--types
```

Deprecated exact alias of –pair-types

### seams –version {#seams-version}

```bash
--version
```

Print version information

### seams –water-types {#seams-water-types}

```bash
--water-types
```

Neighbour graph on these LAMMPS types only (OW).

### seams -k {#seams-k}

```bash
-k
```

k for knn / seeded cages

### seams cages {#seams-cages}

```bash
seams cages FILE
```

Subcommand cages.

### seams chill {#seams-chill}

```bash
seams chill FILE
```

Subcommand chill.

### seams chill-plus {#seams-chill-plus}

```bash
seams chill-plus FILE
```

Subcommand chill-plus.

### seams chill_plus {#id1}

```bash
seams chill_plus FILE
```

Subcommand chill_plus.

### seams cluster {#seams-cluster}

```bash
seams cluster FILE
```

Subcommand cluster.

### seams cn {#seams-cn}

```bash
seams cn FILE
```

Subcommand cn.

### seams density-z {#seams-density-z}

```bash
seams density-z FILE
```

Subcommand density-z.

### seams domains {#seams-domains}

```bash
seams domains FILE
```

Subcommand domains.

### seams f4 {#seams-f4}

```bash
seams f4 FILE
```

Subcommand f4.

### seams fingerprint {#seams-fingerprint}

```bash
seams fingerprint FILE
```

Subcommand fingerprint.

### seams hbonds {#seams-hbonds}

```bash
seams hbonds FILE
```

Subcommand hbonds.

### seams ions {#id2}

```bash
seams ions FILE
```

Subcommand ions.

### seams ovito-chill+ {#seams-ovito-chill}

```bash
seams ovito-chill+ FILE
```

Subcommand ovito-chill+.

### seams ovito_chill_plus {#seams-ovito-chill-plus}

```bash
seams ovito_chill_plus FILE
```

Subcommand ovito_chill_plus.

### seams pairs {#seams-pairs}

```bash
seams pairs FILE
```

Subcommand pairs.

### seams prisms {#seams-prisms}

```bash
seams prisms FILE
```

Subcommand prisms.

### seams rdf {#seams-rdf}

```bash
seams rdf FILE
```

Subcommand rdf.

### seams read {#seams-read}

```bash
seams read FILE
```

Subcommand read.

### seams rings {#seams-rings}

```bash
seams rings FILE
```

Subcommand rings.

### seams steinhardt {#seams-steinhardt}

```bash
seams steinhardt FILE
```

Subcommand steinhardt.

### seams symmetry {#seams-symmetry}

```bash
seams symmetry FILE
```

Subcommand symmetry.

## C API {#c-api}

Each native call and the `seams.h` or `seams_lang.h` function it uses. The C column is the stable ABI.

| Native | C |
| --- | --- |
| seams --atomname | no direct C call |
| seams --axis | no direct C call |
| seams --bins | no direct C call |
| seams --block | no direct C call |
| seams --by-atom | no direct C call |
| seams --cage-edge | no direct C call |
| seams --cell | no direct C call |
| seams --color-types | no direct C call |
| seams --complete | no direct C call |
| seams --config | no direct C call |
| seams --cutoff | no direct C call |
| seams --donors | no direct C call |
| seams --emit-library | no direct C call |
| seams --exclude-atomname | no direct C call |
| seams --exclude-resname | no direct C call |
| seams --family | no direct C call |
| seams --features | no direct C call |
| seams --format | no direct C call |
| seams --frame | no direct C call |
| seams --graph | no direct C call |
| seams --guest-radius | no direct C call |
| seams --guest-types | no direct C call |
| seams --hangle | no direct C call |
| seams --hdist | no direct C call |
| seams --help | no direct C call |
| seams --hops | no direct C call |
| seams --htype | no direct C call |
| seams --incomplete | no direct C call |
| seams --inside | no direct C call |
| seams --ion-cutoff | no direct C call |
| seams --ion-types | no direct C call |
| seams --ions | no direct C call |
| seams --jobs | no direct C call |
| seams --l | no direct C call |
| seams --last | no direct C call |
| seams --layers | no direct C call |
| seams --library | no direct C call |
| seams --max-ring | no direct C call |
| seams --mol | no direct C call |
| seams --pair-types | no direct C call |
| seams --per-atom | no direct C call |
| seams --print-config | no direct C call |
| seams --rebuild | no direct C call |
| seams --resident | no direct C call |
| seams --resname | no direct C call |
| seams --rule | no direct C call |
| seams --shape | no direct C call |
| seams --signature | no direct C call |
| seams --site | no direct C call |
| seams --strict-input | no direct C call |
| seams --subset | no direct C call |
| seams --tpp | no direct C call |
| seams --tum-layers | no direct C call |
| seams --type | no direct C call |
| seams --type-map | no direct C call |
| seams --types | no direct C call |
| seams --version | no direct C call |
| seams --water-types | no direct C call |
| seams -k | no direct C call |
| seams cages | no direct C call |
| seams chill | no direct C call |
| seams chill-plus | no direct C call |
| seams chill_plus | no direct C call |
| seams cluster | no direct C call |
| seams cn | no direct C call |
| seams density-z | no direct C call |
| seams domains | no direct C call |
| seams f4 | no direct C call |
| seams fingerprint | no direct C call |
| seams hbonds | no direct C call |
| seams ions | no direct C call |
| seams ovito-chill+ | no direct C call |
| seams ovito_chill_plus | no direct C call |
| seams pairs | no direct C call |
| seams prisms | no direct C call |
| seams rdf | no direct C call |
| seams read | no direct C call |
| seams rings | no direct C call |
| seams steinhardt | no direct C call |
| seams symmetry | no direct C call |
