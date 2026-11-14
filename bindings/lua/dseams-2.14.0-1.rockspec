rockspec_format = "3.0"
package = "dseams"
version = "2.14.0-1"
source = {
   url = "git+https://github.com/d-SEAMS/seams-core.git",
   tag = "v2.14.0",
}
description = {
   summary = "Lua and Fennel interface to d-SEAMS, the ice and clathrate analysis engine",
   detailed = [[
CHILL and CHILL+ ice labels, Steinhardt order parameters, primitive rings,
cages, RDFs, hydrogen bonds and topological fingerprints for molecular
simulation frames, through the d-SEAMS C API (seams.h). Works with Lua 5.1
to 5.4 and LuaJIT. require("dseams") is the module; yoda is an alias, and
dseams_fnl.fnl gives the kebab-case names to Fennel.

Needs an installed libseams (meson install, conda, Nix or Spack). Point
LuaRocks at it with SEAMS_DIR, or SEAMS_INCDIR and SEAMS_LIBDIR.
]],
   homepage = "https://dseams.info",
   issues_url = "https://github.com/d-SEAMS/seams-core/issues",
   license = "MIT",
   maintainer = "Rohit Goswami <rgoswami@ieee.org>",
   labels = { "science", "chemistry", "molecular-dynamics", "fennel" },
}
dependencies = {
   "lua >= 5.1, < 5.5",
}
external_dependencies = {
   SEAMS = {
      header = "seams.h",
      library = "seams",
   },
}
build = {
   type = "builtin",
   modules = {
      dseams_core = {
         sources = { "bindings/lua/dseams.c" },
         libraries = { "seams" },
         incdirs = { "$(SEAMS_INCDIR)" },
         libdirs = { "$(SEAMS_LIBDIR)" },
      },
      dseams = "bindings/lua/dseams.lua",
      yoda = "bindings/lua/yoda.lua",
   },
   install = {
      lua = {
         dseams_fnl = "bindings/lua/dseams_fnl.fnl",
      },
   },
   copy_directories = {},
}
test = {
   type = "command",
   script = "bindings/lua/test_rock.lua",
   flags = { "." },
}
