# Lua and Fennel

The Lua and Fennel binding provides `require("dseams")`. It lives in its own repository,
[d-SEAMS/yodaStruct](https://github.com/d-SEAMS/yodaStruct), with its API reference, tutorials,
tests and release notes. It links the engine through the C API
(`libseams`, `include/seams.h`) and the C ABI `libseams_lang`
(`bindings/cxx/seams_lang.h`) that this repository builds and installs.

Install it from a yodaStruct source checkout using that repository's
README. The `dseams` registry install on LuaRocks is coming; it is not
available there yet.
