<!-- https://docs.dseams.info/reference/r-api.md -->

# R {#r}

Rd-style entries for the exported functions in [dseams-r](https://github.com/d-SEAMS/dseams-r). The title, `@param`, `@return`, and `@examples` lines are the roxygen block above each function.

## Native API {#native-api}

Generated from the sources named below. A name on this page is a public call in that language.

### cages {#cages}

```r
cages(frame, cutoff = 3.5, type = 0L, k = 4L, graph = "seeded")
```

Hexagonal-cage and double-diamond membership

Parameters

- **frame**: A `dseams_frame`.
- **cutoff**: Neighbour cutoff in angstrom.
- **type**: LAMMPS type. 0 uses the first atom’s type.
- **k**: Neighbours kept by the k-nearest graphs.
- **graph**: `”seeded”`, `”mutual”`, `”cutoff”`, or `”knn-union”`.

Returns

A data.frame with columns `atom`, `cage`, `name`.

Example

```r
frame <- read_frame(dseams_example("ice1h.lammpstrj"), type = -1L)
found <- cages(frame)
attr(found, "n_ih")
```

### cages.dseams_frame {#cages-dseams-frame}

```r
cages.dseams_frame(frame, cutoff = 3.5, type = 0L, k = 4L, graph = "seeded")
```

### chill_plus {#chill-plus}

```r
chill_plus(frame, cutoff = 3.5, type = 0L, graph = "cutoff", rule = "chill+")
```

CHILL+ ice labels

Parameters

- **frame**: A `dseams_frame`.
- **cutoff**: Neighbour cutoff in angstrom. For the mutual graph this is
- **type**: LAMMPS type. 0 uses the first atom’s type.
- **graph**: `”cutoff”` or `”mutual”`.
- **rule**: `”chill+”` or `”ovito-chill+”`. A YAML or Lua snippet with

Returns

A data.frame with columns `atom`, `type`, `label`, `name`.

Example

```r
frame <- read_frame(dseams_example("ice1h.lammpstrj"), type = -1L)
ice <- chill_plus(frame)
table(ice$name)
ice
```

### chill_plus.dseams_frame {#chill-plus-dseams-frame}

```r
chill_plus.dseams_frame(frame, cutoff = 3.5, type = 0L, graph = "cutoff", rule = "chill+")
```

### dseams_example {#dseams-example}

```r
dseams_example(name = NULL)
```

Path of a trajectory shipped with the package

Parameters

- **name**: File name. `NULL` lists the files.

Returns

A path, or a character vector of file names when `name` is `NULL`.

Example

```r
dseams_example()
dseams_example("ice1h.lammpstrj")
```

### dseams_threads {#dseams-threads}

```r
dseams_threads(n = NULL)
```

Engine threads

Parameters

- **n**: Number of threads, or `NULL` to only query.

Returns

The thread count in effect after the call. Without OpenMP this

Example

```r
dseams_threads()
```

### dseams_version {#dseams-version}

```r
dseams_version()
```

Engine version

Returns

The version string of the bundled seams-core engine.

Example

```r
dseams_version()
```

### jump_tau90 {#jump-tau90}

```r
jump_tau90(first, second, dt = 1, oxygen_type = 1L, hydrogen_type = 2L)
```

Jump-rotor tau_90 between two frames

Parameters

- **first,second**: `dseams_frame` objects read with every atom
- **dt**: Time between the two frames.
- **oxygen_type,hydrogen_type**: LAMMPS types of the oxygen and hydrogen

Returns

A number: `dt`, or `NaN`.

Example

```r
frame <- read_frame(dseams_example("exampleTraj.lammpstrj"), type = -1L)
jump_tau90(frame, frame, dt = 1, oxygen_type = 2L, hydrogen_type = 1L)
```

### nframes {#nframes}

```r
nframes(path)
```

### positions {#positions}

```r
positions(frame)
```

Atom positions

Parameters

- **frame**: A `dseams_frame`.

Returns

`positions()` returns a numeric matrix with one row per atom and

Example

```r
frame <- read_frame(dseams_example("ice1h.lammpstrj"), type = -1L)
head(positions(frame))
```

### positions.dseams_frame {#positions-dseams-frame}

```r
positions.dseams_frame(frame)
```

### primitive_rings {#primitive-rings}

```r
primitive_rings(frame, cutoff = 3.5, type = 0L, graph = "mutual", max_depth = 6L)
```

Primitive ring census

Parameters

- **frame**: A `dseams_frame`.
- **cutoff**: Neighbour cutoff in angstrom.
- **type**: LAMMPS type. 0 uses the first atom’s type.
- **graph**: `”mutual”` or `”cutoff”`.
- **max_depth**: Largest ring retained.

Returns

A data.frame with columns `size` and `n`, one row per ring size

Example

```r
frame <- read_frame(dseams_example("ice1h.lammpstrj"), type = -1L)
primitive_rings(frame, graph = "cutoff", max_depth = 6L)
```

### primitive_rings.dseams_frame {#primitive-rings-dseams-frame}

```r
primitive_rings.dseams_frame(frame, cutoff = 3.5, type = 0L, graph = "mutual", max_depth = 6L)
```

### print.dseams_cages {#print-dseams-cages}

```r
print.dseams_cages(x, ...)
```

Parameters

- **x**: A `dseams_cages` data frame.
- **…**: Ignored.

### print.dseams_chill {#print-dseams-chill}

```r
print.dseams_chill(x, ...)
```

Parameters

- **x**: A `dseams_chill` vector.
- **…**: Ignored.

### print.dseams_frame {#print-dseams-frame}

```r
print.dseams_frame(x, ...)
```

### print.dseams_ice {#print-dseams-ice}

```r
print.dseams_ice(x, ...)
```

Parameters

- **x**: A `dseams_ice` data frame.
- **…**: Ignored.

### print.dseams_q4_q6 {#print-dseams-q4-q6}

```r
print.dseams_q4_q6(x, ...)
```

Parameters

- **x**: A `dseams_q4_q6` data frame.
- **…**: Ignored.

### print.dseams_rings {#print-dseams-rings}

```r
print.dseams_rings(x, ...)
```

Parameters

- **x**: A `dseams_rings` data frame.
- **…**: Ignored.

### print.dseams_steinhardt {#print-dseams-steinhardt}

```r
print.dseams_steinhardt(x, ...)
```

Parameters

- **x**: A `dseams_steinhardt` data frame.
- **…**: Ignored.

### print.dseams_system {#print-dseams-system}

```r
print.dseams_system(x, ...)
```

Parameters

- **x**: A `dseams_system`.
- **…**: Ignored.

### q4_q6 {#q4-q6}

```r
q4_q6(frame, cutoff = 3.5, type = 0L)
```

Local q4 and q6 on the cutoff shell

Parameters

- **frame**: A `dseams_frame`.
- **cutoff**: Neighbour cutoff in angstrom.
- **type**: LAMMPS type. 0 uses the first atom’s type.

Returns

A data.frame with columns `atom`, `q4`, `q6`, `q4bar`, `q6bar`.

Example

```r
frame <- read_frame(dseams_example("ice1h.lammpstrj"), type = -1L)
head(q4_q6(frame))
```

### q4_q6.dseams_frame {#q4-q6-dseams-frame}

```r
q4_q6.dseams_frame(frame, cutoff = 3.5, type = 0L)
```

### read_frame {#read-frame}

```r
read_frame(path, frame = 1L, type = 0L)
```

Read one trajectory frame

Parameters

- **path**: Path to a LAMMPS dump, or an XYZ file.
- **frame**: 1-based frame index.
- **type**: LAMMPS type filter.
- **x**: A `dseams_frame`, for `print()`.
- **object**: A `dseams_frame`, for `summary()`.
- **…**: Ignored.

Returns

`read_frame()` returns an object of class `dseams_frame`.

Example

```r
path <- dseams_example("ice1h.lammpstrj")
nframes(path)
frame <- read_frame(path, type = -1L)
frame
```

### seams_cages {#seams-cages}

```r
seams_cages(sys, cutoff = 3.5, graph = 2L, k = 4L, complete = FALSE)
```

Disjoint cage population counts.

Parameters

- **sys**: A `dseams_system`.
- **cutoff**: Neighbour cutoff in angstrom.
- **graph**: `2` seeded, `1` mutual, `0` cutoff. These are the same
- **k**: Neighbour count for the seeded graph.
- **complete**: When true, run ring-adjacent completion.

Returns

An integer vector `c(n_hex, n_cubic, n_water, n_both)`, class `dseams_cage_counts`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_cages(sys)
```

### seams_chill {#seams-chill}

```r
seams_chill(sys, cutoff = 3.5, type = 0L, graph = 0L, rule = 0L)
```

Original CHILL population counts

Parameters

- **sys**: A `dseams_system`.
- **cutoff**: Neighbour cutoff in angstrom.
- **type**: Species filter. Zero keeps every atom.
- **graph**: `0` cutoff, `1` mutual.
- **rule**: `0` published, `1` OVITO.

Returns

An integer vector of length 9, class `dseams_chill`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_chill(sys)
```

### seams_contact_pairs {#seams-contact-pairs}

```r
seams_contact_pairs(sys, spec)
```

Mutual nearest unlike contacts.

Parameters

- **sys**: A `dseams_system`.
- **spec**: Site-table text such as `”1=cationHead”`.

Returns

An integer vector `c(n_pairs, n_cation, n_anion)`, class `dseams_contacts`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_contact_pairs(sys, "1=cationHead")
```

### seams_coordination {#seams-coordination}

```r
seams_coordination(sys, rmax = 3.5, bins = 8L)
```

Coordination number integrated to `rmax`.

Parameters

- **sys**: A `dseams_system`.
- **rmax**: Cutoff in angstrom.
- **bins**: Bin count used by the integral.

Returns

A numeric scalar.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_coordination(sys)
```

### seams_density_axis {#seams-density-axis}

```r
seams_density_axis(sys, type = 0L, bins = 4L, axis = 2L)
```

Number density along one cell axis.

Parameters

- **sys**: A `dseams_system`.
- **type**: Species filter.
- **bins**: Bin count.
- **axis**: 0, 1, or 2.

Returns

A numeric vector, class `dseams_density`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_density_axis(sys)
```

### seams_domains {#seams-domains}

```r
seams_domains(sys, spec, cutoff = 3.5, subset = 0L)
```

Largest polar or apolar domain.

Parameters

- **sys**: A `dseams_system`.
- **spec**: Site-table text.
- **cutoff**: Neighbour cutoff in angstrom.
- **subset**: `0` apolar, `1` polar.

Returns

A numeric vector `c(n_subset, largest, percolation)`, class `dseams_domains`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_domains(sys, "1=waterO")
```

### seams_fingerprint {#seams-fingerprint}

```r
seams_fingerprint(sys, cutoff = 3.5, hops = 2L)
```

Topological fingerprint summary.

Parameters

- **sys**: A `dseams_system`.
- **cutoff**: Neighbour cutoff in angstrom.
- **hops**: Neighbourhood depth.

Returns

An integer vector `c(key_length, ring_sum)`, class `dseams_fingerprint`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_fingerprint(sys)
```

### seams_hbonds {#seams-hbonds}

```r
seams_hbonds(sys, cutoff = 3.5, oxygen_type = 1L, hydrogen_type = 2L, dist = 2.42, angle = 30)
```

Hydrogen-bond count.

Parameters

- **sys**: A `dseams_system`.
- **cutoff**: Heavy-atom cutoff in angstrom.
- **oxygen_type,hydrogen_type**: Species numbers.
- **dist**: Acceptor-H distance in angstrom.
- **angle**: O-O-H angle in degrees.

Returns

An integer.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("exampleTraj.lammpstrj"), type = -1L)
seams_hbonds(sys, oxygen_type = 2L, hydrogen_type = 1L)
```

### seams_monolayer {#seams-monolayer}

```r
seams_monolayer(sys, cutoff = 3.5, max_depth = 6L, sheet_area = 100)
```

Primitive-ring histogram of a monolayer.

Parameters

- **sys**: A `dseams_system`.
- **cutoff**: Neighbour cutoff in angstrom.
- **max_depth**: Largest ring size.
- **sheet_area**: Sheet area in square angstrom.

Returns

An integer vector, class `dseams_monolayer`. `counts\[k + 1\]` is size `k`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_monolayer(sys)
```

### seams_neighbours_knn {#seams-neighbours-knn}

```r
seams_neighbours_knn(sys, k = 4L, cutoff = 3.5, mutual = TRUE)
```

Column count of a k-nearest neighbour list.

Parameters

- **sys**: A `dseams_system`.
- **k**: Neighbour count.
- **cutoff**: Nominal cutoff in angstrom.
- **mutual**: When true, keep the mutual graph.

Returns

An integer, class `dseams_nnz`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_neighbours_knn(sys)
```

### seams_prisms {#seams-prisms}

```r
seams_prisms(sys, cutoff = 3.5)
```

Perfect prism counts by ring size.

Parameters

- **sys**: A `dseams_system`.
- **cutoff**: Neighbour cutoff in angstrom.

Returns

An integer vector of length 7, class `dseams_prisms`. Index 7 is size 6.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_prisms(sys)
```

### seams_rdf {#seams-rdf}

```r
seams_rdf(sys, rmax = 6, bins = 8L)
```

Partial radial distribution.

Parameters

- **sys**: A `dseams_system`.
- **rmax**: Maximum radius in angstrom.
- **bins**: Bin count.

Returns

A list of class `dseams_rdf` with `r`, `g`, and `volume`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_rdf(sys)
```

### seams_rings_vertices {#seams-rings-vertices}

```r
seams_rings_vertices(sys, cutoff = 3.5, max_depth = 6L, index = 0L)
```

Atom indices of one primitive ring.

Parameters

- **sys**: A `dseams_system`.
- **cutoff**: Neighbour cutoff in angstrom.
- **max_depth**: Largest ring size.
- **index**: Ring index, from 0.

Returns

An integer vector, class `dseams_vertices`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_rings_vertices(sys)
```

### seams_running_cn {#seams-running-cn}

```r
seams_running_cn(sys, rmax = 6, bins = 8L)
```

Running coordination number.

Parameters

- **sys**: A `dseams_system`.
- **rmax**: Maximum radius in angstrom.
- **bins**: Bin count.

Returns

A numeric vector, class `dseams_running_cn`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_running_cn(sys)
```

### seams_set_num_threads {#seams-set-num-threads}

```r
seams_set_num_threads(n)
```

Set the OpenMP width used by the next analysis.

Parameters

- **n**: Thread count. Zero restores the library default.

Returns

The width the next analysis requests, as an integer. A

Example

```r
seams_set_num_threads(1L)
```

### seams_signature_cages {#seams-signature-cages}

```r
seams_signature_cages(sys, spec, cutoff = 3.5, graph = 0L)
```

Closed cages of one ring census.

Parameters

- **sys**: A `dseams_system`.
- **spec**: Census text such as `”512”`.
- **cutoff**: Neighbour cutoff in angstrom.
- **graph**: `0` cutoff or `1` mutual.

Returns

An integer vector `c(n_cages, n_atoms)`, class `dseams_signature`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_signature_cages(sys, "512")
```

### seams_site_table {#seams-site-table}

```r
seams_site_table(sys, spec)
```

Per-atom site kind.

Parameters

- **sys**: A `dseams_system`.
- **spec**: Site-table text such as `”1=waterO”`.

Returns

An integer vector, class `dseams_sites`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_site_table(sys, "1=waterO")
```

### seams_steinhardt {#seams-steinhardt}

```r
seams_steinhardt(sys, order, cutoff = 3.5, average = FALSE)
```

Local Steinhardt parameter on the cutoff shell.

Parameters

- **sys**: A `dseams_system`.
- **order**: Spherical-harmonic degree. 3, 4, 6, 8, or 12.
- **cutoff**: Neighbour cutoff in angstrom.
- **average**: When true, the second element is the Lechner-Dellago average.

Returns

A list of class `dseams_order` with elements `ql` and `qlbar`.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_steinhardt(sys, order = 6L)
```

### seams_system_count {#seams-system-count}

```r
seams_system_count(sys)
```

Atom count stored on a system.

Parameters

- **sys**: A `dseams_system`.

Returns

An integer. Zero when the handle has no frame.

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
seams_system_count(sys)
```

### seams_system_new {#seams-system-new}

```r
seams_system_new()
```

Empty system handle

Returns

A `dseams_system`. `print()` and `summary()` return the handle,

Example

```r
sys <- seams_system_new()
sys
```

### seams_system_read {#seams-system-read}

```r
seams_system_read(sys, path, frame = 1L, type = -1L)
```

Read one frame into a system.

Parameters

- **sys**: A `dseams_system`.
- **path**: Trajectory path.
- **frame**: 1-based frame index.
- **type**: Species filter. Negative keeps every atom.

Returns

The atom count. `sys` is updated in place. On failure the

Example

```r
sys <- seams_system_new()
seams_system_read(sys, dseams_example("ice1h.lammpstrj"))
```

### signature_cages {#signature-cages}

```r
signature_cages(frame, spec = "512", cutoff = 3.5, type = 0L, graph = "cutoff")
```

Closed cages for one ring census

Parameters

- **frame**: A `dseams_frame`.
- **spec**: Ring census such as `”512”` (twelve pentagons, the sI small
- **cutoff**: Neighbour cutoff in angstrom.
- **type**: LAMMPS type. 0 uses the first atom’s type.
- **graph**: `”cutoff”` or `”mutual”`.

Returns

An integer vector `c(n_cages, n_atoms)`.

Example

```r
frame <- read_frame(dseams_example("exampleTraj.lammpstrj"), type = 2L)
signature_cages(frame, "512")
```

### signature_cages.dseams_frame {#signature-cages-dseams-frame}

```r
signature_cages.dseams_frame(frame, spec = "512", cutoff = 3.5, type = 0L, graph = "cutoff")
```

### steinhardt {#steinhardt}

```r
steinhardt(frame, order_l = 6L, cutoff = 3.5, type = 0L, graph = NULL)
```

Steinhardt order parameters

Parameters

- **frame**: A `dseams_frame`.
- **order_l**: Spherical-harmonic degree. 3, 4, 6, 8, and 12 are computed.
- **cutoff**: Neighbour cutoff in angstrom.
- **type**: LAMMPS type. 0 uses the first atom’s type.
- **graph**: Neighbour graph. `NULL` is the cutoff shell. `”mutual”` is

Returns

A data.frame with columns `atom`, `ql`, `qlbar`.

Example

```r
frame <- read_frame(dseams_example("ice1h.lammpstrj"), type = -1L)
summary(steinhardt(frame, order_l = 6L)$ql)
```

### steinhardt.dseams_frame {#steinhardt-dseams-frame}

```r
steinhardt.dseams_frame(frame, order_l = 6L, cutoff = 3.5, type = 0L, graph = NULL)
```

### summary.dseams_cages {#summary-dseams-cages}

```r
summary.dseams_cages(object, ...)
```

Parameters

- **object**: A `dseams_cages` data frame.

### summary.dseams_chill {#summary-dseams-chill}

```r
summary.dseams_chill(object, ...)
```

Parameters

- **object**: A `dseams_chill` vector.

### summary.dseams_frame {#summary-dseams-frame}

```r
summary.dseams_frame(object, ...)
```

### summary.dseams_ice {#summary-dseams-ice}

```r
summary.dseams_ice(object, ...)
```

Parameters

- **object**: A `dseams_ice` data frame.

### summary.dseams_q4_q6 {#summary-dseams-q4-q6}

```r
summary.dseams_q4_q6(object, ...)
```

Parameters

- **object**: A `dseams_q4_q6` data frame.

### summary.dseams_rings {#summary-dseams-rings}

```r
summary.dseams_rings(object, ...)
```

Parameters

- **object**: A `dseams_rings` data frame.

### summary.dseams_steinhardt {#summary-dseams-steinhardt}

```r
summary.dseams_steinhardt(object, ...)
```

Parameters

- **object**: A `dseams_steinhardt` data frame.

### summary.dseams_system {#summary-dseams-system}

```r
summary.dseams_system(object, ...)
```

Parameters

- **object**: A `dseams_system`.

### symmetry {#symmetry}

```r
symmetry(frame, shell_rmsd = 0.12)
```

Shell symmetry label for every atom

Parameters

- **frame**: A `dseams_frame`.
- **shell_rmsd**: Largest RMSD, in the units of the file, at which a

Returns

An integer vector, one code per atom: 0 ice Ih, 1 ice Ic,

Example

```r
frame <- read_frame(dseams_example("ice1h.lammpstrj"), type = -1L)
table(symmetry(frame))
```

### symmetry.dseams_frame {#symmetry-dseams-frame}

```r
symmetry.dseams_frame(frame, shell_rmsd = 0.12)
```

## C API {#c-api}

Each native call and the `seams.h` or `seams_lang.h` function it uses. The C column is the stable ABI.

| Native | C |
| --- | --- |
| cages | no direct C call |
| cages.dseams_frame | seams_frame_cages, seams_frame_ids, seams_frame_nop, seams_last_error |
| chill_plus | no direct C call |
| chill_plus.dseams_frame | seams_frame_chill_plus_rule, seams_frame_ids, seams_frame_nop, seams_frame_types, seams_last_error |
| dseams_example | no direct C call |
| dseams_threads | no direct C call |
| dseams_version | seams_lang_version |
| jump_tau90 | seams_frames_jump_tau90, seams_last_error |
| nframes | seams_last_error, seams_path_nframes |
| positions | no direct C call |
| positions.dseams_frame | seams_frame_nop, seams_frame_positions |
| primitive_rings | no direct C call |
| primitive_rings.dseams_frame | seams_frame_rings, seams_last_error |
| print.dseams_cages | no direct C call |
| print.dseams_chill | no direct C call |
| print.dseams_frame | seams_frame_nop |
| print.dseams_ice | no direct C call |
| print.dseams_q4_q6 | no direct C call |
| print.dseams_rings | no direct C call |
| print.dseams_steinhardt | no direct C call |
| print.dseams_system | seams_system_count |
| q4_q6 | no direct C call |
| q4_q6.dseams_frame | seams_frame_ids, seams_frame_nop, seams_frame_q4_q6, seams_last_error |
| read_frame | seams_frame_read, seams_last_error |
| seams_cages | no direct C call |
| seams_chill | no direct C call |
| seams_contact_pairs | no direct C call |
| seams_coordination | no direct C call |
| seams_density_axis | no direct C call |
| seams_domains | no direct C call |
| seams_fingerprint | no direct C call |
| seams_hbonds | no direct C call |
| seams_monolayer | no direct C call |
| seams_neighbours_knn | no direct C call |
| seams_prisms | no direct C call |
| seams_rdf | no direct C call |
| seams_rings_vertices | no direct C call |
| seams_running_cn | no direct C call |
| seams_set_num_threads | no direct C call |
| seams_signature_cages | no direct C call |
| seams_site_table | no direct C call |
| seams_steinhardt | no direct C call |
| seams_system_count | no direct C call |
| seams_system_new | no direct C call |
| seams_system_read | no direct C call |
| signature_cages | no direct C call |
| signature_cages.dseams_frame | seams_frame_signature, seams_last_error |
| steinhardt | no direct C call |
| steinhardt.dseams_frame | seams_frame_ids, seams_frame_nop, seams_frame_steinhardt, seams_frame_steinhardt_graph, seams_last_error |
| summary.dseams_cages | no direct C call |
| summary.dseams_chill | no direct C call |
| summary.dseams_frame | seams_frame_nop |
| summary.dseams_ice | no direct C call |
| summary.dseams_q4_q6 | no direct C call |
| summary.dseams_rings | no direct C call |
| summary.dseams_steinhardt | no direct C call |
| summary.dseams_system | no direct C call |
| symmetry | no direct C call |
| symmetry.dseams_frame | seams_frame_nop, seams_frame_symmetry, seams_last_error |
