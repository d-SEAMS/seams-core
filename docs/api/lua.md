# Lua and Fennel

The Lua and Fennel binding provides `require("dseams")`. It lives in its own repository,
[d-SEAMS/yodaStruct](https://github.com/d-SEAMS/yodaStruct), with its API reference, tutorials,
tests and release notes. It links the engine through the C API
(`libseams`, `include/seams.h`) and the C ABI `libseams_lang`
(`bindings/cxx/seams_lang.h`) that this repository builds and installs.

The binding documentation is at <https://lua.dseams.info/>, and the Lua and Fennel walkthrough is [on Odysee](https://odysee.com/dseams-lua:9a7ffc26dff9590174a707db848f5f5d6d7a737f).

Install it from a yodaStruct source checkout using that repository's
README. The `dseams` registry install on LuaRocks is coming; it is not
available there yet.

`docs/orgmode/contributing/repositories.org` describes how the repositories fit together.
