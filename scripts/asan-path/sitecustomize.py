# The loader has already mapped whatever LD_PRELOAD named. /bin/sh
# segfaults when libclang_rt.tsan is still in the environment, and ldd
# is a shell script, so drop only that runtime. libasan stays in this
# process: a child that dlopens libseams aborts with "ASan runtime does
# not come first" once the preload is gone. Host compilers and shells
# are not those children; they get the address runtime removed below,
# or LeakSanitizer reports their own allocations.
#
# A child interpreter is a new process. sys.executable is the real
# Python, not the asan-path wrapper, and it dlopens libseams. That
# load is too late for the thread runtime's static TLS, so put
# TSAN_RUNTIME back into LD_PRELOAD for Python and Lua children only.
# Shells, and plumed, gmx and lmp, stay without it.
import os
import subprocess

preload = os.environ.get("LD_PRELOAD", "")
tsan = os.environ.get("TSAN_RUNTIME", "")
if preload and tsan and tsan in preload.split(":"):
    kept = [part for part in preload.split(":") if part and part != tsan]
    if kept:
        os.environ["LD_PRELOAD"] = ":".join(kept)
    else:
        os.environ.pop("LD_PRELOAD", None)

_INTERPRETERS = {"python", "python3", "lua", "lua5.4", "luajit"}
# Host compilers and shells are not instrumented. A preload inherited
# from this interpreter makes LeakSanitizer report their allocations.
# The process that dlopens libseams is not one of these names.
_HOST_TOOLS = {
    "cc", "c++", "gcc", "g++", "clang", "clang++",
    "mpicc", "mpicxx", "mpic++", "mpirun", "mpiexec", "orterun",
    "sh", "bash",
}
asan = os.environ.get("ASAN_TEST_PRELOAD", "")


def _argv_base(args):
    if isinstance(args, (str, bytes)) or not args:
        return ""
    return os.path.basename(os.fsdecode(args[0]))


def _interpreter(args):
    base = _argv_base(args)
    return base in _INTERPRETERS or base.startswith("python3.")


def _drop_asan(env):
    if not asan:
        return env
    kept = [part for part in env.get("LD_PRELOAD", "").split(":") if part and part != asan]
    if kept:
        env["LD_PRELOAD"] = ":".join(kept)
    else:
        env.pop("LD_PRELOAD", None)
    return env


if (tsan and os.path.isfile(tsan)) or (asan and os.path.isfile(asan)):
    _popen = subprocess.Popen

    class _Popen(_popen):
        def __init__(self, args, *pos, **kw):
            if not kw.get("shell"):
                base = _argv_base(args)
                if asan and base in _HOST_TOOLS:
                    env = kw.get("env")
                    env = os.environ.copy() if env is None else dict(env)
                    kw["env"] = _drop_asan(env)
                elif tsan and os.path.isfile(tsan) and _interpreter(args):
                    env = kw.get("env")
                    env = os.environ.copy() if env is None else dict(env)
                    current = env.get("LD_PRELOAD", "")
                    if tsan not in current.split(":"):
                        env["LD_PRELOAD"] = tsan + ((":" + current) if current else "")
                    kw["env"] = env
            super().__init__(args, *pos, **kw)

    subprocess.Popen = _Popen
