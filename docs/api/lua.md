# Lua API

`docs/examples/minimal_lua.lua` reads one trajectory and prints the CHILL+ clathrate count. `meson test docs_examples` runs it.

```{literalinclude} ../../../docs/examples/minimal_lua.lua
:language: lua
```

`require("dseams")` binds Lua 5.4 and LuaJIT to `include/seams.h`.
Fennel uses `(require :dseams)`. A system, a neighbour list, and a ring
census are userdata. `:close()` frees the handle. Lua 5.4 also frees it
at the end of a `<close>` scope. A status other than `SEAMS_OK` is a
Lua error whose text is `seams_last_error`.

`bindings/lua/dseams.lua` is the name table. `bindings/lua/dseams.c` is
the C module `dseams_core`. Hyphenated keys (`chill-plus`,
`ring-network`, `cages-by-signature`, `running-cn`, `site-table`,
`set-num-threads`) are the same functions, for Fennel callers.

## Handles

`dseams.system(xyz, box, types, ids)` calls `seams_system_create` and,
when `xyz` is present, `seams_system_set_frame`. `dseams.read(path, o)`
calls `seams_system_read`. `o.type` keeps one atom type. `o.all = true`
keeps every atom. `o.frame` selects the frame, counting from 1.
`dseams.nframes(path)` is `seams_nframes`. `cloud.nop` is the atom
count. `cloud:close()` destroys the system.

`dseams.version()` is `seams_version_string`. `dseams.label_name(code)`
is `seams_label_name`. `dseams.set_num_threads(n)` and
`dseams.num_threads()` are the matching C calls.
`dseams.PUBLISHED` and `dseams.OVITO` select the CHILL+ rule.
`dseams.CUTOFF`, `dseams.MUTUAL`, and `dseams.SEEDED` select the graph.

## Names

Each row is a live call. The 1.x column is the alias on the same table.
Options default to cutoff 3.5 and atom type 0 unless the row says
otherwise. The C reference, with one printed line per export, is
`docs/api/c.md`.

| 1.x name | Lua call | C call |
| --- | --- | --- |
| `dseams.system` | `dseams.system(xyz, box)` | `seams_system_create`, `seams_system_set_frame` |
| `readLammpsTrj` | `dseams.read(path, {all = true})` | `seams_system_read` |
| `readLammpsTrjO` | `dseams.read(path, {type = 2})` | `seams_system_read` |
| `getIceTypePlus`, `getCorrelPlus`, `getIceTypePlusNoPrint` | `dseams.chill_plus(cloud, {cutoff = 3.5, rule = dseams.PUBLISHED})` | `seams_chill_plus` |
| `getIceType`, `getCorrel` | `dseams.chill(cloud, o)` | `seams_chill` |
| `neighListO` | `dseams.neighbors(cloud, {cutoff = 3.5, graph = "cutoff"})` | `seams_neighbours_build` |
| `kNearestNeighbourList` | `dseams.knn(cloud, {k = 4, mutual = true})` | `seams_neighbours_knn` |
| `dseams.q4q6` | `dseams.q4q6(cloud)` | `seams_q4q6` |
| `steinhardtQl` | `dseams.steinhardtQl(cloud, {l = 6, average = true})` | `seams_steinhardt` |
| `ringNetwork`, `bulkPrimitiveRings` | `dseams.rings(cloud, {max_depth = 6})` | `seams_rings_find`, `seams_rings_histogram`, `seams_rings_fraction` |
| `getPrimitiveRings` | `dseams.getPrimitiveRings(cloud, depth)` | `seams_rings_vertices` |
| `seededCageAffiliation` | `dseams.cages(cloud, {graph = "seeded", k = 4})` | `seams_cages` |
| `findBySignature` | `dseams.cages_by_signature(cloud, "512", {graph = "cutoff"})` | `seams_signature_cages` |
| `calcRDF3D` | `dseams.rdf(cloud, {rmax = 6.0, bins = 8})` | `seams_rdf` |
| `calcRunningCN` | `dseams.running_cn(cloud, o)` | `seams_running_cn` |
| `calcCN` | `dseams.cn(cloud, {cutoff = 3.5, bins = 8})` | `seams_coordination` |
| `getHbondNetwork` | `dseams.hbonds(cloud, {oxygen_type = 1, h_type = 2, dist = 2.42, angle = 30})` | `seams_hbonds` |
| `densityByType` | `dseams.density(cloud, {bins = 4, axis = 2})` | `seams_density_axis` |
| `topologyFingerprint` | `dseams.fingerprint(cloud, {hops = 2})` | `seams_fingerprint` |
| `dseams.site_table` | `dseams.site_table(cloud, "1=waterO")` | `seams_site_table` |
| `contactPairs` | `dseams.pairs(cloud, {spec = "1=cationHead"})` | `seams_contact_pairs` |
| `domainStats` | `dseams.domain(cloud, {spec = "1=waterO", subset = 0})` | `seams_domains` |
| `dseams.prisms` | `dseams.prisms(cloud)` | `seams_prisms` |
| `shapeMatching` | `dseams.monolayer(cloud, {max_depth = 6, sheet_area = 100})` | `seams_monolayer` |

`graph` is `cutoff`, `mutual`, or `seeded`. `rule` is `published` or
`ovito`, or the constants above. `dseams.steinhardtQl` takes `l` (or
`order`) and sets the Lechner-Dellago average when `average` or `qlbar`
is true. The return is `{ql = ..., qlBar = ...}`.
`dseams.knn` takes any `k`. `mutual = false` keeps a directed list.
`dseams.rings` also returns `by_size`, `fraction`, `symmetry`, and
`vertices`. `dseams.cages` returns `n_hex`, `n_cubic`, `n_water`,
`n_both`, and per-atom `labels`. `dseams.cages_by_signature` returns
`n_cages` and `n_atoms`.

## Printed results

`examples/lua/ice.lua` prints `2.14.0 cubic 8`.
`examples/lua/read.lua` prints `dseams.read nop=250`.
`examples/lua/chill_plus.lua` prints `clathrate 46`.
`docs/examples/lua_analyses.lua` prints the remaining calls.
On ice Ih that is `hexagonal 128`, `q6bar 0.562448843702`,
`cages 128 0 0 0`, `knn 512`, `sheet 256`, and `vert 6 0`.
On the sI guest frame the `512` signature prints `sig 2 40`.

## LuaJIT FFI

`gen_ffi.py` turns `include/seams.h` into a cdef and `ffi.load`s
`libseams` (`SEAMS_LIB` overrides the library name). The cdef is the
declaration span of the header, so each exported function is callable
as `C.seams_*` with the C buffer convention. `tests/lua_ffi.lua` loads
that module, labels the cubic cell through `seams_chill_plus`, and
prints `ffi ok`.


## Structure still

```{figure} ../tutorials/molecules/lua/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug lua` writes the 2400 by 1800 PNG.
```


<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### CUTOFF

Purpose: Integer code published on the dseams table.

Parameters: none

Returns: The integer code.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_lua.lua`

### MUTUAL

Purpose: Integer code published on the dseams table.

Parameters: none

Returns: The integer code.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_lua.lua`

### OVITO

Purpose: Integer code published on the dseams table.

Parameters: none

Returns: The integer code.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_lua.lua`

### PUBLISHED

Purpose: Integer code published on the dseams table.

Parameters: none

Returns: The integer code.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_lua.lua`

### SEEDED

Purpose: Integer code published on the dseams table.

Parameters: none

Returns: The integer code.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_lua.lua`

### WATER_O

Purpose: Integer code published on the dseams table.

Parameters: none

Returns: The integer code.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_lua.lua`

### bulkPrimitiveRings

Purpose: Alias of dseams.rings.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### cages

Purpose: Score hexagonal and double-diamond cage membership.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### cages-by-signature

Purpose: Fennel spelling of dseams.cages_by_signature.

Parameters: cloud; spec; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### cages_by_signature

Purpose: Count closed cages whose ring census matches a spec.

Parameters: cloud; spec; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### calcCN

Purpose: Alias of dseams.cn.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### calcRDF3D

Purpose: Alias of dseams.rdf.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### calcRunningCN

Purpose: Alias of dseams.running_cn.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### chill

Purpose: Write one original CHILL label per atom.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### chill-plus

Purpose: Fennel spelling of dseams.chill_plus.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### chill_plus

Purpose: Write one CHILL+ label per atom.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### cloud.box

Purpose: Return the nine LAMMPS cell numbers stored on a cloud.

Parameters: none

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### cloud.close

Purpose: Free the system stored on a cloud.

Parameters: none

Returns: Nothing.

Ownership: Frees the seams_system inside the cloud. A second call is a no-op.

Errors: None.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### cloud.nop

Purpose: Atom count stored on a cloud when it was built.

Parameters: none

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### cn

Purpose: Return the coordination number at the requested cutoff.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### contactPairs

Purpose: Alias of dseams.pairs.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### core

Purpose: The C module table the Lua functions call.

Parameters: none

Returns: The dseams_core module table.

Ownership: Lua owns the module. The caller does not free C pointers from it by hand; the wrappers do.

Errors: A failing C call raised by a wrapper is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### density

Purpose: Return the number density along one cell axis.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### densityByType

Purpose: Alias of dseams.density.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### domain

Purpose: Measure the largest polar or apolar domain.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### domainStats

Purpose: Alias of dseams.domain.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### findBySignature

Purpose: Alias of dseams.cages_by_signature.

Parameters: cloud; spec; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### fingerprint

Purpose: Build the topological fingerprint of the bonded graph.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### getCorrel

Purpose: Alias of dseams.chill.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### getCorrelPlus

Purpose: Alias of dseams.chill_plus.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### getHbondNetwork

Purpose: Alias of dseams.hbonds.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### getIceType

Purpose: Alias of dseams.chill.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### getIceTypePlus

Purpose: Alias of dseams.chill_plus.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### getIceTypePlusNoPrint

Purpose: Alias of dseams.chill_plus.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### getPrimitiveRings

Purpose: Return the atom indices of one primitive ring.

Parameters: none

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### hbonds

Purpose: Count hydrogen bonds on one frame.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### kNearestNeighbourList

Purpose: Alias of dseams.knn.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### knn

Purpose: Build a k-nearest neighbour list.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### label_name

Purpose: Map a CHILL code to its name.

Parameters: code

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### monolayer

Purpose: Histogram primitive rings in a monolayer.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### neighListO

Purpose: Alias of dseams.neighbors.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### neighbors

Purpose: Build a cutoff, mutual, or seeded neighbour list.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### neighbours

Purpose: Alias of dseams.neighbors.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### nframes

Purpose: Count frames in a trajectory path.

Parameters: path

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### num_threads

Purpose: Report the OpenMP width the next call will use.

Parameters: none

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### pairs

Purpose: Count mutual nearest unlike contacts of a site table.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### prisms

Purpose: Count perfect and imperfect prisms by ring size.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### q4q6

Purpose: Write local q4 and q6 on the cutoff shell.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### rdf

Purpose: Compute a partial radial distribution.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### read

Purpose: Read one frame from a trajectory path.

Parameters: path; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### readLammpsTrj

Purpose: Lua binding for the matching C call.

Parameters: none

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### readLammpsTrjO

Purpose: Alias of dseams.read.

Parameters: path; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### ring-network

Purpose: Fennel spelling of dseams.ringNetwork.

Parameters: cloud; depth

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### ringNetwork

Purpose: Return the same ring census with an explicit maximum size.

Parameters: cloud; depth

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### rings

Purpose: Return primitive-ring sizes, symmetry, a histogram, and the six-ring fraction.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### running-cn

Purpose: Fennel spelling of dseams.running_cn.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### running_cn

Purpose: Integrate the partial radial distribution into a running coordination number.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### seededCageAffiliation

Purpose: Alias of dseams.cages.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### set-num-threads

Purpose: Fennel spelling of dseams.set_num_threads.

Parameters: n

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### set_num_threads

Purpose: Set the OpenMP width for the next call.

Parameters: n

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### shapeMatching

Purpose: Alias of dseams.monolayer.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### site-table

Purpose: Fennel spelling of dseams.site_table.

Parameters: cloud; spec

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### site_table

Purpose: Assign one site kind per atom from a site-table spec.

Parameters: cloud; spec

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### steinhardtQl

Purpose: Write local Steinhardt q_l, and the neighbour average when asked.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### system

Purpose: Build a cloud from coordinates and a cell, or an empty cloud when coordinates are omitted.

Parameters: xyz; box; types; ids

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### topologyFingerprint

Purpose: Alias of dseams.fingerprint.

Parameters: cloud; o

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`

### version

Purpose: Return the library version text.

Parameters: none

Returns: The Lua table or userdata described for this call.

Ownership: A system userdata is freed by cloud:close(). Lua 5.4 also frees it at the end of a <close> scope. Returned tables are owned by Lua.

Errors: A status other than SEAMS_OK is a Lua error whose text is seams_last_error.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_lua.lua`
<!-- API_SYMBOLS_END -->
