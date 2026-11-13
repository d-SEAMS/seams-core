# Lua functions

`require("dseams")` is the Lua 5.4 and LuaJIT module in
`bindings/lua`. Fennel loads that same module. The calls go through
`include/seams.h`. Names the header does not export are listed in
`docs/orgmode/explanation/c-api-gaps.org`.

`require("yoda")` is the same table. The how-to is
`docs/orgmode/howto/lua.org`.
