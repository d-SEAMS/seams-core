# Lua and Fennel

Classify a frame from Lua or Fennel. The recording follows `require("dseams")` and the same call from Fennel.

```{eval-rst}
.. include:: ../videos/embeds/lua.rst
```

The binding is the Lua and Fennel library in [d-SEAMS/yodaStruct](https://github.com/d-SEAMS/yodaStruct). `require("dseams")` in Lua and `(require :dseams)` in Fennel call `libseams` (`include/seams.h`). The rock is not published on LuaRocks, so `luarocks install dseams` does not install it. The README builds the module with Meson from that repository.

Option tables and return values are the [Lua reference](https://d-seams.github.io/yodaStruct/reference/lua.html). Package names are the {doc}`Lua API </stable/lua>`.
