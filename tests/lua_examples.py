#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
# SPDX-License-Identifier: MIT
"""Run the Lua and Fennel examples and compare them with the C++ oracle."""

import os
import pathlib
import subprocess
import sys


def library_loader(lib):
    """Dynamic linker and directories that already load libseams."""
    interp = ""
    dirs = []
    ldd = subprocess.run(["ldd", lib], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    for line in ldd.stdout.splitlines():
        stripped = line.strip()
        if "=>" in line:
            path = line.split("=>", 1)[1].strip().split()
            if path and path[0].startswith("/"):
                parent = str(pathlib.Path(path[0]).parent)
                if parent not in dirs:
                    dirs.append(parent)
        elif stripped.startswith("/") and "ld-linux" in stripped:
            interp = stripped.split()[0]
    libdir = str(pathlib.Path(lib).parent)
    if libdir not in dirs:
        dirs.insert(0, libdir)
    return interp, dirs


def run(cmd, env, cwd):
    done = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if done.returncode != 0:
        sys.stderr.write("$ " + " ".join(cmd) + "\n")
        sys.stderr.write(done.stdout)
        sys.stderr.write(done.stderr)
        raise SystemExit(done.returncode)
    return done.stdout


def parse(text):
    found = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        found[(parts[0], parts[1])] = parts[2:]
    return found


def lua_cmd(lua, args, env, loader):
    interp, dirs = loader
    run_env = env.copy()
    if interp and pathlib.Path(interp).is_file():
        return [interp, "--library-path", ":".join(dirs), lua] + args, run_env
    previous = run_env.get("LD_LIBRARY_PATH", "")
    run_env["LD_LIBRARY_PATH"] = ":".join(dirs + ([previous] if previous else []))
    return [lua] + args, run_env


def compare(oracle, other, keys):
    left = parse(oracle)
    right = parse(other)
    for key in keys:
        if key not in right:
            raise SystemExit("missing {} {}".format(key[0], key[1]))
        numeric = any(("." in token) or ("e" in token.lower()) for token in left[key] + right[key])
        if numeric:
            if len(left[key]) != len(right[key]):
                raise SystemExit("{} {} oracle {} lua {}".format(key[0], key[1], left[key], right[key]))
            for a, b in zip(left[key], right[key]):
                if abs(float(a) - float(b)) > 1e-8:
                    raise SystemExit("{} {} oracle {} lua {}".format(key[0], key[1], left[key], right[key]))
        elif left[key] != right[key]:
            raise SystemExit("{} {} oracle {} lua {}".format(key[0], key[1], left[key], right[key]))


def main(argv):
    root = pathlib.Path(argv[1])
    work = pathlib.Path(argv[2])
    oracle = argv[3]
    ffi = argv[4]
    lib = argv[5]
    lua_mods = argv[6:]
    if len(lua_mods) % 2 != 0:
        raise SystemExit("lua modules are interpreter, cpath pairs")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["SEAMS_LIB"] = lib
    env["LUA_PATH"] = str(work) + "/?.lua;" + str(root / "bindings" / "lua") + "/?.lua;;"
    env["LUA_CPATH"] = ""
    frame = work / "ih.lammpstrj"
    oracle_text = run([oracle, str(frame), str(root)], env, root)
    full_keys = []
    for name in ("ih", "si", "interface", "tip"):
        for field in ("n", "chill", "rings", "q6", "nnz"):
            full_keys.append((name, field))
    fennel_keys = [("interface", "chill"), ("tip", "rings")]
    for field in (
        "chill0", "q3", "q6bar", "vert", "cages", "knn", "rdf", "rcn", "cn",
        "hbond", "rho", "fp", "site", "pairs", "domain", "prism", "sheet",
    ):
        full_keys.append(("ih", field))
    full_keys.extend([
        ("read", "si"), ("read", "tip"), ("read", "xyz"), ("read", "con"), ("read", "tri"),
        ("si", "sig"),
    ])
    fennel_lib = env.get("FENNEL_LUA", "")
    if not fennel_lib:
        raise SystemExit("FENNEL_LUA is not set")
    loader = library_loader(lib)
    compiler = root / "tests" / "lua_fennel_compile.lua"
    scripts = sorted((root / "examples").rglob("*.lua"))
    fennel_scripts = sorted((root / "examples").rglob("*.fnl"))
    compile_cmd, compile_env = lua_cmd(
        lua_mods[0],
        [str(compiler), fennel_lib],
        env,
        loader,
    )
    compiled = {}
    for script in fennel_scripts:
        compiled[script] = run(compile_cmd + [str(script)], compile_env, root)
    for library in sorted((root / "bindings" / "lua").glob("*.fnl")):
        generated = work / (library.stem + ".lua")
        generated.write_text(run(compile_cmd + [str(library)], compile_env, root))
    for index in range(0, len(lua_mods), 2):
        lua = lua_mods[index]
        cpath = lua_mods[index + 1]
        run_env = env.copy()
        run_env["LUA_CPATH"] = cpath + "/?.so;;"

        def invoke(args):
            cmd, child = lua_cmd(lua, args, run_env, loader)
            return run(cmd, child, root)

        for script in scripts:
            invoke([str(script), str(root)])
        for name in ("lua_handles.lua", "lua_gaps.lua", "lua_fixtures.lua"):
            extra = [str(frame), str(root)] if name == "lua_fixtures.lua" else []
            text = invoke([str(root / "tests" / name)] + extra)
            if name == "lua_fixtures.lua":
                compare(oracle_text, text, full_keys)
        if pathlib.Path(lua).name.startswith("lua5"):
            invoke([str(root / "tests" / "lua_close.lua")])
        fennel_text = []
        for script in fennel_scripts:
            generated = work / (script.stem + "-" + pathlib.Path(lua).name + ".lua")
            generated.write_text(compiled[script])
            fennel_text.append(invoke([str(generated), str(root)]))
        compare(oracle_text, "\n".join(fennel_text), fennel_keys)
        if "luajit" in pathlib.Path(lua).name:
            invoke([str(root / "tests" / "lua_ffi.lua"), ffi])
    howto = (root / "docs" / "orgmode" / "howto" / "lua.org").read_text()
    for script in list(scripts) + list(fennel_scripts):
        marker = script.relative_to(root).as_posix()
        if marker not in howto:
            raise SystemExit("how-to does not mention " + marker)
    print("lua examples ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
