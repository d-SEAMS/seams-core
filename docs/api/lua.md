# Lua and Fennel

The Lua and Fennel binding is the LuaRocks rock `dseams`, `require("dseams")`. It lives in its own repository,
[d-SEAMS/yodaStruct](https://github.com/d-SEAMS/yodaStruct), with its API reference, tutorials,
tests and release notes. It links the engine through the C API
(`libseams`, `include/seams.h`) and the C ABI `libseams_lang`
(`bindings/cxx/seams_lang.h`) that this repository builds and installs.

```bash
luarocks install dseams
```

The binding release with the same version as seams-core is built
against that seams-core tag. `docs/orgmode/contributing/repositories.org`
describes how the repositories fit together.
