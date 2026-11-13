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


# Lua and Fennel

`require("dseams")` is the module. `require("yoda")` is the same table.
`bindings/lua/seams.fnl` is the Fennel module: keyword options, a
`chill+` function, and `seams.close`. Lua 5.4 and LuaJIT both load
`dseams_core.so`. The FFI cdef generated from `include/seams.h` is
LuaJIT only.

Errors are raised with `error`. A call does not return `nil` and a
message. `pcall` is how a caller keeps the message. The text is
`seams_last_error` when the library set one.

A cloud, a neighbour list, a ring census, and a cage set are handles.
`__gc` frees the C object. Lua 5.4 `__close` frees it at the end of a
scope. `:close()` is the same call. Closing twice is safe.

`ldoc bindings/lua/config.ld` writes the native reference. The book
links that HTML at `/lua/index.html`.
`bindings/lua/dseams-2.14.0-1.rockspec` is the rock. From
`bindings/lua`:

```bash
luarocks make dseams-2.14.0-1.rockspec SEAMS_LIBDIR="$SEAMS_LIBDIR" SEAMS_INCDIR="$SEAMS_INCDIR"
luarocks test dseams-2.14.0-1.rockspec
```

`SEAMS_LIBDIR` holds `libseams.so` and `SEAMS_INCDIR` holds the
generated `seams_version.h`. `busted` runs
`bindings/lua/spec/ownership_spec.lua`. `meson test lua_examples`
runs the examples on both interpreters.

## Quickstart

```lua
local dseams = require("dseams")
local cloud = dseams.system(xyz, box)
local labels = dseams.chill_plus(cloud, {cutoff = 3.5, rule = dseams.PUBLISHED})
cloud:close()
```

`examples/lua/ice.lua` is the eight-atom cell. It prints:

```text
2.14.0 cubic 8
```

The Fennel module does the same count. `examples/fennel/cubic.fnl`
prints:

```text
fennel cubic 8
```

```lua
(local seams (require :seams))
(local labels ((. seams :chill+) cloud {:cutoff 3.5 :rule :published}))
```

## What the examples check

`meson test lua_examples` compares Lua and Fennel with the C++ oracle
on ice Ih, the GenIce sI cell, the basal interface, and the rings
tutorial. The fields are the atom count, CHILL+, the ring histogram,
q6, the neighbour pairs, original CHILL, q3, q6-bar, ring vertices,
cages, the signature, k-nearest neighbours, the RDF, the running
coordination, the coordination number, hydrogen bonds, the density,
the fingerprint, the site table, contact pairs, domains, prisms, and
the monolayer. The numbers match the C++ path to `1e-8`.

`tests/lua_handles.lua` rejects a cell with a zero edge and prints
`handles ok`. The message is `box edge lengths must be positive`.

## Ownership

`dseams.system` and `dseams.read` return a cloud that owns one
`seams_system`. Neighbour rows and ring tables are copies: the
function closes the C handle before it returns. A fingerprint, a cage
set, and a neighbour list created through `dseams.core` stay open
until `:close()` or collection.

The frame copied into a neighbour list or a ring census is independent
of the cloud. Closing the cloud does not free those results, because
they no longer hold the C handle.

```{literalinclude} ../../../examples/lua/ice.lua
:language: lua
```

```{literalinclude} ../../../examples/fennel/cubic.fnl
:language: text
```

```{literalinclude} ../../../bindings/lua/seams.fnl
:language: text
:lines: 1-40
```
