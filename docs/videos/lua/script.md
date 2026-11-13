# Lua and Fennel shot list

Pixi is not required. The build is the meson tree that produced
`bindings/lua/lua54/dseams_core.so` and
`bindings/lua/luajit/dseams_core.so`, with `libseams.so` beside it.
`LUA_PATH` is `bindings/lua/?.lua`. `LUA_CPATH` is the interpreter's
module directory under `$SEAMS_BUILD`. `FENNEL_LUA` is the Fennel
compiler. Commands run from the repository root.
`docs_video_tapes` runs every tape. `lua_examples` checks the same
census against `tests/lua_oracle.cpp`.

## Scene 1

Lua 5.4 and LuaJIT both run `examples/lua/ice.lua`. LuaJIT then runs
`tests/lua_ffi.lua` on the generated cdef.

```bash
LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/ice.lua .
LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/luajit/?.so" luajit examples/lua/ice.lua .
SEAMS_LIB=seams luajit tests/lua_ffi.lua "$SEAMS_BUILD/bindings/lua/dseams_ffi.lua"
```

Expected: `2.14.0 cubic 8` and `ffi ok`.

## Scene 2

```bash
LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/read.lua .
LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/chill_plus.lua .
LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/neighbors.lua .
LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/steinhardt.lua .
LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/rings.lua .
```

Expected: `dseams.read nop=250`, `clathrate 46`, `neighbors rows=46 nnz=184`,
`q6 0.595133082718`, and `rings 3 37 4 273 5 34 6 167 7 0 8 0 fraction 1.000000`.

## Scene 3

`docs/examples/lua_analyses.lua` calls CHILL, the Lechner-Dellago
average, cages, k-nearest neighbours, the RDF, coordination, hydrogen
bonds, density, the fingerprint, the site table, contact pairs,
domains, prisms, the monolayer, ring vertices, and the cage signature.

```bash
LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 docs/examples/lua_analyses.lua
```

Expected: `hexagonal 128`, `q6bar 0.562448843702`, `cages 128 0 0 0`,
`knn 512`, `sheet 256`, and `sig 2 40`.

## Scene 4

```bash
lua5.4 tests/lua_fennel_compile.lua "$FENNEL_LUA" examples/fennel/interface.fnl > "$VIDEO_OUT/interface.lua" && LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 "$VIDEO_OUT/interface.lua" .
```

Expected: `interface chill 1 352 434 49 3 25 0 0 0`.
