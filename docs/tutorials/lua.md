# Lua and Fennel tutorial

`require("dseams")` is the module in `bindings/lua/dseams.lua`.
`require("yoda")` is the same table. Lua 5.4 and LuaJIT both load it.
Fennel compiles to Lua and then uses the same module: `(require :dseams)`.

Build the shared module when an interpreter is present:

```bash
meson setup build -Dwith_lua=auto
meson compile -C build
```

`LUA_PATH` points at `bindings/lua`. `LUA_CPATH` points at the directory
that holds `dseams_core.so`. On this tree that is
`build/bindings/lua/lua54` for Lua 5.4 and `build/bindings/lua/luajit`
for LuaJIT. `meson test lua_examples` runs the scripts under
`examples/lua` and `examples/fennel` on both interpreters and checks
each census against `tests/lua_oracle.cpp`.

## Eight-atom cubic cell

`examples/lua/ice.lua` is the cell from `docs/examples/ice.c`.
Every atom is cubic. Lua 5.4 and LuaJIT both print:

```
2.14.0 cubic 8
```

```bash
LUA_PATH=bindings/lua/?.lua \
  LUA_CPATH=build/bindings/lua/lua54/?.so \
  lua5.4 examples/lua/ice.lua .
```

## Read, CHILL+, neighbours, q6, rings

`examples/lua/read.lua` keeps oxygen type 2 of
`input/traj/exampleTraj.lammpstrj` and prints:

```
dseams.read nop=250
```

`examples/lua/chill_plus.lua` labels `input/traj/genice_sI.lammpstrj`
and prints:

```
cubic 0 hexagonal 0 water 0 interfacial 0 clathrate 46 interClathrate 0 unclassified 0 reCubic 0 reHex 0
```

`examples/lua/neighbors.lua` prints `neighbors rows=46 nnz=184`.
`examples/lua/steinhardt.lua` prints `q6 0.595133082718`.
`examples/lua/rings.lua` prints:

```
rings 3 37 4 273 5 34 6 167 7 0 8 0 fraction 1.000000
```

## The rest of the call set

`docs/examples/lua_analyses.lua` calls the remaining names on
`input/traj/ice1h.lammpstrj`, then the cage signature on the sI guest
frame. The printed lines are:

```
hexagonal 128
q6bar 0.562448843702
cages 128 0 0 0
knn 512
rdf 8 4158.66460658 0
rcn 26.2047244094
cn 4.03149606299
hbond 0
rho 4 0.0307791111112
fp 16 256
site 8
pairs 0 128 0
domain 0 0 0
prism 0
sheet 256
vert 6 0
sig 2 40
```

`hexagonal 128` is `dseams.chill` with the published rule.
`q6bar` is `dseams.steinhardtQl` with `average = true`, the
Lechner-Dellago mean. `cages` is the seeded census
(`n_hex`, `n_cubic`, `n_water`, `n_both`). `knn` is the mutual
four-neighbour list. `rdf`, `rcn`, and `cn` are the radial distribution,
the running coordination, and `dseams.cn`. `hbond 0` is
`dseams.hbonds` on this oxygen frame. `rho` is `dseams.density` along
z. `fp` is `dseams.fingerprint`. `site` is the first code from
`dseams.site_table`. `pairs` and `domain` are the contact-pair and
domain counts. `prism` is `dseams.prisms`. `sheet` is the six-ring bin
of `dseams.monolayer`. `vert` is the first ring from
`dseams.rings`, which is `seams_rings_vertices`. `sig 2 40` is
`dseams.cages_by_signature` for `512`: 2 cages and 40 atoms.

The same numbers are what `tests/lua_fixtures.lua` writes for
`lua_examples`.

## Fennel

`examples/fennel/interface.fnl` is the basal interface.
`examples/fennel/rings.fnl` is the rings tutorial frame.
`tests/lua_fennel_compile.lua` compiles each file with Fennel's
`compileString`. The interface script prints:

```
interface chill 1 352 434 49 3 25 0 0 0
```

The rings script prints:

```
tip rings 0 0 0 37 273 34 167 0 0
```

## LuaJIT FFI

`bindings/lua/gen_ffi.py` writes `dseams_ffi.lua` from `include/seams.h`.
`ffi.load` opens `libseams`. Every declaration in the header is in that
cdef, so a LuaJIT program can call `seams_system_read`,
`seams_steinhardt`, `seams_cages`, `seams_rdf`, `seams_hbonds`,
`seams_rings_vertices`, and the rest of the C API directly.
`tests/lua_ffi.lua` does that for the cubic cell and prints `ffi ok`.


## Structure still

```{figure} ../tutorials/molecules/lua/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug lua` writes the 2400 by 1800 PNG.
```

