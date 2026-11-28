<!-- https://docs.dseams.info/user/howto/migrate-1x.md -->

# Migrate from d-SEAMS 1.x {#migrate-from-d-seams-1-x}

1.x was a Lua script plus `yodaStruct -c lua_inputs/config.yml`.
The config named the trajectory and switched analyses on. 2.14 is one
`seams` subcommand, or one call in `seams.h` / `seams_lang.h`, for
each of those switches.

1.x Lua had plain globals only. It has no `dseams` table and no
`require("dseams")`. v1.0.0 and v1.0.1 register the same globals.
Later Lua names first shipped in
[yodaStruct](https://github.com/d-SEAMS/yodaStruct) 2.x, which calls
`seams.h`. This repository does not register Lua names.

Each row is one analysis. Where the row has both a command and a C
call, they read the same trajectory.

## yodaStyle config {#yodastyle-config}

| 1.x key | 2.14 command | 2.14 C call |
| --- | --- | --- |
| bulk.topologicalNetworkCriterion | seams cages | seams_frame_cages |
| bulk.bondOrderParameters | seams steinhardt | seams_frame_steinhardt_graph |
| topoOneDim | seams prisms | seams_frame_prisms |
| topoTwoDim | seams_frame_monolayer | seams_frame_monolayer |

`topoTwoDim` has no `seams` subcommand. The monolayer census is
`seams_frame_monolayer` (language frame) and `seams_monolayer`
(`seams.h`). The counts buffer is indexed by ring size and must hold
at least `max_depth + 1` integers.
On the mixed TIP4P dump `input/traj/exampleTraj.lammpstrj` the call
takes oxygen type 2.

## Lua globals and later Lua names {#lua-globals-and-later-lua-names}

The first column is the v1.0.1 global. v1.0.0 registers the same
names. A row marked new in 2.x had no 1.x global; the parenthetical
is the yodaStruct release that added the Lua names. Each later column
is the name that release first registered. `require("dseams")` starts
at yodaStruct 2.2.0. The last two columns are the 2.14 replacement.

| 1.x global (v1.0.1) | yodaStruct 2.0.0 | yodaStruct 2.2.0 | yodaStruct 2.4.1 | yodaStruct 2.7.0 | yodaStruct 2.8.0 | yodaStruct 2.9.0 | 2.14 command | 2.14 C call |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| readFrame / readFrameOnlyOne | readLammpsTrjO | dseams.read |  | readLammpsTrj |  |  | seams read | seams_frame_read |
| chillPlus_cij / chillPlus_iceType | getCorrelPlus / getIceTypePlus | dseams.chill_plus |  |  |  |  | seams chill-plus | seams_frame_chill_plus |
| chill_cij / chill_iceType | getCorrel / getIceType | dseams.chill |  |  |  |  | seams chill | seams_frame_chill |
| getPrimitiveRings | ringNetwork |  |  |  |  |  | seams rings | seams_frame_rings |
| averageQ6 | steinhardtQl |  |  |  |  |  | seams steinhardt | seams_frame_steinhardt_graph |
| bulkTopologicalNetworkCriterion | seededCageAffiliation | dseams.cages |  |  |  |  | seams cages | seams_frame_cages |
| new in 2.x (yodaStruct 2.9.0) |  |  |  |  |  | findBySignature / dseams.cages_by_signature | seams cages --signature | seams_frame_signature_edge |
| neighborList / bondNetworkByIndex | neighListO | dseams.neighbors |  |  |  |  | seams_frame_neighbours | seams_neighbours_build |
| new in 2.x (yodaStruct 2.0.0) | kNearestNeighbourList | dseams.knn |  |  |  |  | seams_frame_neighbours | seams_neighbours_build |
| calcRDF |  |  | calcRDF3D | dseams.rdf |  |  | seams rdf | seams_frame_rdf |
| prismAnalysis |  |  |  |  |  |  | seams prisms | seams_frame_prisms |
| ringAnalysis |  |  |  |  |  |  | seams_frame_monolayer | seams_frame_monolayer |
| getHbondNetwork |  |  |  | dseams.hbonds |  |  | seams hbonds | seams hbonds |
| new in 2.x (yodaStruct 2.4.1) |  |  | calcCN / dseams.cn |  |  |  | seams cn | seams cn |
| new in 2.x (yodaStruct 2.7.0) |  |  |  | densityByType / dseams.density |  |  | seams density-z | seams density-z |
| new in 2.x (yodaStruct 2.8.0) |  |  |  |  | topologyFingerprint / dseams.fingerprint |  | seams fingerprint | seams fingerprint |
| new in 2.x (yodaStruct 2.7.0) |  |  |  | contactPairs / dseams.pairs |  |  | seams pairs | seams pairs |
| new in 2.x (yodaStruct 2.7.0) |  |  |  | domainStats / dseams.domain |  |  | seams domains | seams domains |

The 2.0.0 column names C++ functions in this tree (`readLammpsTrjO`,
`getCorrelPlus`, `ringNetwork`, `neighListO`, and the others). 1.x
called them under the globals in the first column.
`calcRDF` in 1.x is the in-plane RDF of the monolayer script
(`rdf2Danalysis_AA`). The 3D RDF that `seams rdf` prints first reached
Lua as `calcRDF3D`.
`ringAnalysis` is the 1.x `topoTwoDim` ring census and `prismAnalysis`
the 1.x `topoOneDim` one. No later column registers a new Lua name
for either.

The 2015 rule and the OVITO rule are selected with `--rule`.
The walk-through is
[the published rule and the OVITO rule](../../explanation/chill-plus-rules.md).

Neighbours have no `seams` subcommand. There are two rows.
`seams_frame_neighbours` counts the directed pairs on the language
frame. `seams_neighbours_build` builds that list through `seams.h`.
The cutoff row uses graph `cutoff`. The k-nearest row uses graph
`mutual`. The two calls report the same pair count.

Hydrogen bonds, coordination, the density profile, the fingerprint,
contact pairs, and domains repeat the command in the last column.
The functions in `seams.h` are `seams_hbonds`, `seams_coordination`,
`seams_density_axis`, `seams_fingerprint`, `seams_contact_pairs`, and
`seams_domains`.

On `input/traj/exampleTraj.lammpstrj`:

```bash
seams hbonds input/traj/exampleTraj.lammpstrj --type 2 --htype 1
seams cn input/traj/exampleTraj.lammpstrj --pair-types 1,2 --cutoff 4.5
seams density-z input/traj/exampleTraj.lammpstrj --type 2 --bins 4
seams fingerprint input/traj/exampleTraj.lammpstrj --type 2
seams pairs input/traj/exampleTraj.lammpstrj --site 1=cationHead,2=anion
seams domains input/traj/exampleTraj.lammpstrj --site 1=cationHead,2=anion --subset polar
```

The lines are `nop 250 hbonds 460`,
`# types 1 2 cutoff 4.5 cn 9.144 nI 500 nJ 250 volume 288324`,
`classes 89`, `contact-pair count 0 nCation 0 nAnion 250`, and
`subset polar n 750 largest 750 P 1`.
`density-z` starts:

```text
# z rho
29.557 0.000887889
```

CHILL+ and rings use the cubic cell `input/traj/mW_cubic.lammpstrj`,
type 1, cutoff 3.5 Å, graph `cutoff`.
`seams chill-plus --rule published` prints `nop 4096 cubic 4096`.
`seams rings` prints `nop 4096 rings 8192` and `6 8192`.
`seams_frame_chill_plus` and `seams_frame_rings` take that same frame.
The signature row is `seams cages --signature 512 --cage-edge none` on
`input/traj/genice_sI.lammpstrj`. That command prints
`nop 46 graph seeded signature 5:12 cages 2 atoms 40`.
`seams_frame_signature_edge` takes that same frame.
