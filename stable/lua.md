<!-- https://docs.dseams.info/stable/lua.md -->

# Lua and Fennel {#lua-and-fennel}

The Lua and Fennel library is [d-SEAMS/yodaStruct](https://github.com/d-SEAMS/yodaStruct). `require("dseams")` in Lua and `(require :dseams)` in Fennel call `libseams` (`include/seams.h`). `require("yoda")` resolves to the same module.

```bash
git clone https://github.com/d-SEAMS/yodaStruct.git
cd yodaStruct
luarocks make
```

Option tables and return values are the [Lua reference](https://d-seams.github.io/yodaStruct/reference/lua.md).

`require("dseams")` and `(require :dseams)` call the same engine as the CLI. On `input/traj/ice1h.lammpstrj` (128 oxygens, type 1) published CHILL+ at 3.5 angstrom is 128 hexagonal:

```default
nop 128 frame 1 box 15.6457 14.7071 18.073
nop 128 hexagonal 128
```

The first classification is the [Lua and Fennel tutorial](../fronts/lua.md).
