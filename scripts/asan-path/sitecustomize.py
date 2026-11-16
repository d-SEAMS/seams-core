# The loader has already mapped whatever LD_PRELOAD named. /bin/sh
# segfaults when libclang_rt.tsan is still in the environment, and ldd
# is a shell script, so drop only that runtime. libasan has to stay: a
# child that is not itself the first library to map the address runtime
# aborts with "ASan runtime does not come first" once the preload is
# gone. This process keeps either mapping.
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


def _interpreter(args):
    if isinstance(args, (str, bytes)) or not args:
        return False
    base = os.path.basename(os.fsdecode(args[0]))
    return base in _INTERPRETERS or base.startswith("python3.")


if tsan and os.path.isfile(tsan):
    _popen = subprocess.Popen

    class _Popen(_popen):
        def __init__(self, args, *pos, **kw):
            if not kw.get("shell") and _interpreter(args):
                env = kw.get("env")
                env = os.environ.copy() if env is None else dict(env)
                current = env.get("LD_PRELOAD", "")
                if tsan not in current.split(":"):
                    env["LD_PRELOAD"] = tsan + ((":" + current) if current else "")
                kw["env"] = env
            super().__init__(args, *pos, **kw)

    subprocess.Popen = _Popen
