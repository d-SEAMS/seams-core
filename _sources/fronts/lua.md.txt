# Lua and Fennel

Classify a frame from Lua or Fennel. The recording follows `require("dseams")` and the same call from Fennel.

```{eval-rst}
.. include:: ../videos/embeds/lua.rst
```

The binding is the Lua and Fennel library in [d-SEAMS/yodaStruct](https://github.com/d-SEAMS/yodaStruct). `require("dseams")` in Lua and `(require :dseams)` in Fennel call `libseams` (`include/seams.h`).

```bash
git clone https://github.com/d-SEAMS/yodaStruct.git
cd yodaStruct
luarocks make
```

Option tables and return values are the [Lua reference](https://d-seams.github.io/yodaStruct/reference/lua.html).

The recording uses `input/traj/ice1h.lammpstrj`. `seams read` prints `nop 128 frame 1`, and published CHILL+ at 3.5 angstrom prints `nop 128 hexagonal 128`:

```bash
seams read input/traj/ice1h.lammpstrj
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

Package names are the {doc}`Lua API </stable/lua>`.
