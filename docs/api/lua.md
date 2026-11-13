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
