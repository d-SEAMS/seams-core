# Copyright (c) d-SEAMS developers
"""asv timings for one call of each binding.

The published medians are benchmarks/binding-overhead.csv, from
scripts/bench_binding_overhead.py. This module starts the same workers so
``asv run`` repeats the calls. Samples here include one round trip to the
worker. The CSV times the call inside the worker.
"""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_READY = {}


def _campaign():
    if "mod" not in _READY:
        path = ROOT / "scripts" / "bench_binding_overhead.py"
        spec = importlib.util.spec_from_file_location("binding_campaign", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _READY["mod"] = module
    return _READY["mod"]


def _built():
    module = _campaign()
    if "bins" not in _READY:
        lib = module.find_library()
        _READY["lib"] = lib
        _READY["bins"] = module.compile_drivers(lib)
        _READY["lua"] = module.find_lua_module()
        frames = module.ensure_frames(
            ROOT / "benchmarks" / "logs" / "binding-overhead" / "frames"
        )
        _READY["frames"] = {
            molecules: path.relative_to(ROOT).as_posix() for molecules, path in frames.items()
        }
    return module


class BindingCall:
    """One serial call. asv repeats it. The CSV run uses eight repeats."""

    params = [
        ["c", "cxx", "rust", "fortran", "python", "julia", "r", "lua", "freud", "ovito"],
        ["read", "knn", "steinhardt", "chill", "rings"],
        [4096, 32768, 110592],
    ]
    param_names = ["binding", "operation", "molecules"]
    number = 1
    repeat = 3

    def setup(self, binding, operation, molecules):
        module = _built()
        allowed = {
            "freud": ("knn", "steinhardt"),
            "ovito": ("read", "knn", "steinhardt", "rings"),
        }
        if binding in allowed and operation not in allowed[binding]:
            raise NotImplementedError(binding + " " + operation)
        key = (binding, operation)
        if key in _READY.get("workers", {}):
            self.worker = _READY["workers"][key]
            self.path = _READY["frames"][molecules]
            self.op = "freud-" + operation if binding == "freud" else operation
            self.threads = 1
            return
        lib = _READY["lib"]
        env = module.child_env(lib)
        bins = _READY["bins"]
        if binding == "c":
            argv = [str(bins["c"])]
        elif binding == "cxx":
            argv = [str(bins["cxx"])]
        elif binding == "fortran":
            argv = [str(bins["fortran"])]
        elif binding == "rust":
            argv = [str(bins["rust"])]
        elif binding == "python" or binding == "freud":
            argv = [module.sys.executable, str(ROOT / "scripts" / "bench_binding_overhead.py"),
                    "--worker", "python"]
        elif binding == "ovito":
            gl = Path("/usr/lib/x86_64-linux-gnu/libGL.so.1")
            if gl.is_file():
                gl_dir = lib.parent.parent / "gl-link"
                gl_dir.mkdir(exist_ok=True)
                link = gl_dir / "libGL.so.1"
                if not link.exists():
                    link.symlink_to(gl)
                env["LD_LIBRARY_PATH"] = str(gl_dir) + ":" + env["LD_LIBRARY_PATH"]
            env["OVITO_THREAD_COUNT"] = "1"
            env["QT_QPA_PLATFORM"] = "offscreen"
            argv = [module.sys.executable, str(ROOT / "scripts" / "bench_binding_overhead.py"),
                    "--worker", "ovito"]
        elif binding == "julia":
            import shutil
            julia = shutil.which("julia")
            if not julia:
                raise NotImplementedError("julia")
            argv = [julia, "--project=" + str(ROOT / "bindings" / "julia" / "Seams"),
                    str(ROOT / "benchmarks" / "binding_overhead.jl")]
        elif binding == "r":
            import shutil
            rscript = shutil.which("Rscript")
            if not rscript:
                raise NotImplementedError("R")
            env["SEAMS_R_SO"] = str(bins["rso"])
            argv = [rscript, str(ROOT / "benchmarks" / "binding_overhead.R")]
        elif binding == "lua":
            exe, so = _READY["lua"]
            env["LUA_CPATH"] = str(so.parent / "?.so") + ";;"
            argv = [str(exe), str(ROOT / "benchmarks" / "binding_overhead.lua")]
        else:
            raise NotImplementedError(binding)
        worker = module.Worker(argv, env)
        _READY.setdefault("workers", {})[key] = worker
        self.worker = worker
        self.path = _READY["frames"][molecules]
        self.op = "freud-" + operation if binding == "freud" else operation
        self.threads = 1

    def time_call(self, binding, operation, molecules):
        del binding, operation, molecules
        self.worker.ask(self.threads, self.op, self.path)
