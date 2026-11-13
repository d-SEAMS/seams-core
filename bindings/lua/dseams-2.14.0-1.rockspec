rockspec_format = "3.0"
package = "dseams"
version = "2.14.0-1"
source = {
  url = "git://github.com/d-SEAMS/seams-core",
}
description = {
  summary = "Lua and Fennel bindings for d-SEAMS",
  detailed = "require('dseams') over include/seams.h. Errors are raised with error. Lua 5.4 and LuaJIT both load the C module. Fennel is bindings/lua/seams.fnl.",
  homepage = "https://docs.dseams.info/",
  license = "MIT",
}
dependencies = {
  "lua >= 5.1",
}
build = {
  type = "builtin",
  modules = {
    dseams_core = {
      sources = { "dseams.c" },
      libraries = { "seams" },
      -- SEAMS_INCDIR is the build dir that holds seams_version.h.
      -- SEAMS_LIBDIR is the directory that holds libseams.so.
      incdirs = { "../../include", "$(SEAMS_INCDIR)" },
      libdirs = { "$(SEAMS_LIBDIR)" },
    },
    dseams = "dseams.lua",
    yoda = "yoda.lua",
  },
}
test_dependencies = {
  "busted",
}
test = {
  type = "busted",
}
