# Copyright (c) d-SEAMS developers
"""LUA_CPATH for the interpreter that is about to load dseams_core.so.

The build produces both bindings/lua/lua54 and bindings/lua/luajit.
A glob of that directory sorts luajit first. PUC Lua 5.4 then maps the
LuaJIT module and segfaults on the first require.
"""

from pathlib import Path


def cpath_for(build: Path, lua_exe: str) -> str:
    base = Path(build) / "bindings" / "lua"
    name = Path(lua_exe).name.lower()
    sub = "luajit" if "luajit" in name else "lua54"
    module = base / sub / "dseams_core.so"
    if not module.is_file():
        raise SystemExit(f"missing {module} for interpreter {lua_exe}")
    return str(module.parent / "?.so") + ";;"
