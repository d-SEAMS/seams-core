# dseams

`require("dseams")` is the Lua module over `include/seams.h`.
`require("yoda")` is the same table. `seams.fnl` is the Fennel module.

A failed call raises `error` with the library message. It does not
return `nil` and a message. Handles free themselves in `__gc`, and Lua
5.4 also frees them in `__close`. `:close()` is the same call.

Lua 5.4 and LuaJIT both load `dseams_core.so`. The FFI cdef is LuaJIT
only.

```bash
ldoc bindings/lua/config.ld
```
