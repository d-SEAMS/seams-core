# Lua and Fennel

The Lua and Fennel library is [d-SEAMS/yodaStruct](https://github.com/d-SEAMS/yodaStruct). `require("dseams")` in Lua and `(require :dseams)` in Fennel call `libseams` (`include/seams.h`). `require("yoda")` resolves to the same module.

The rock is not published on LuaRocks. The README builds the module with Meson from that repository. Option tables and return values are the [Lua reference](https://d-seams.github.io/yodaStruct/reference/lua.html).

The first classification is the {doc}`Lua and Fennel tutorial </fronts/lua>`.
