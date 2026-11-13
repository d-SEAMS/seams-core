"""Airspeed Velocity benchmarks for the paper's timed arms.

Imports and frame setup run outside the timed region. Each timed method
is one serial call of the arm on one core.
"""

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORE = ["taskset", "-c", "0"]
ENV = os.environ.copy()
ENV.update({
    "OMP_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OVITO_THREAD_COUNT": "1",
    "OPENBLAS_NUM_THREADS": "1",
})


def _run(cmd):
    subprocess.run(
        cmd,
        cwd=ROOT,
        env=ENV,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


class NeighbourArm:
    def setup(self):
        self.cmd = CORE + [
            str(ROOT / "bbdir/tests/bench_neigh_only"),
            "--traj", "input/traj/genice_sH.lammpstrj",
            "--frame", "1",
            "--type", "1",
            "--reps", "1",
            "--warmup", "0",
            "--cutoff", "3.5",
            "--k", "4",
            "--threads", "1",
        ]

    def time_mutual_neighbours(self):
        _run(self.cmd)


class DdcArm:
    """findDDC on 1,000-molecule ideal ice Ic. One warm-up, one timed call."""

    def setup(self):
        frame = ROOT / "bbdir" / "ddc-regression" / "frames" / "ic-1000.lammpstrj"
        if not frame.is_file():
            subprocess.run(
                ["python3", "paper/tex/figures/ladder_frames.py",
                 str(frame.parent), "5"],
                cwd=ROOT, check=True)
        self.cmd = CORE + [
            str(ROOT / "bbdir/tests/bench_ddc_ladder"),
            str(frame),
            "1",
            "1",
        ]

    def time_find_ddc(self):
        _run(self.cmd)


class RingArm:
    def setup(self):
        self.cmd = CORE + [
            str(ROOT / "bbdir/tests/bench_stages"),
            "input/traj/genice_sH.lammpstrj",
            "1",
            "1",
            "1",
        ]

    def time_primitive_rings(self):
        _run(self.cmd)


class CampaignArms:
    def setup(self):
        bench = str(ROOT / "bbdir/tests/bench_paper")
        frame = "input/traj/genice_sI.lammpstrj"
        modes = {
            "chill": "chill",
            "q6": "q6-full",
            "seeded": "seeded-adj",
            "rings": "rings",
        }
        self.commands = {
            name: CORE + [bench, frame, "1", "0", "1", mode]
            for name, mode in modes.items()
        }

    def time_chill(self):
        _run(self.commands["chill"])

    def time_q6(self):
        _run(self.commands["q6"])

    def time_seeded_cages(self):
        _run(self.commands["seeded"])

    def time_rings(self):
        _run(self.commands["rings"])


class ChillPlusIce:
    """Cutoff CHILL+ on 1000-molecule ideal ice. One thread, one call.

    The neighbour rows come from the linked-cell list. q3 is still the
    published left-to-right sum, so the labels do not move.
    """

    def setup(self):
        import ctypes
        import math
        lib_path = os.environ.get("SEAMS_LANG_LIB")
        if not lib_path:
            raise NotImplementedError("SEAMS_LANG_LIB is not set")
        lib = ctypes.CDLL(lib_path)
        lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_chill_plus.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        bond = 2.75
        lattice = 4.0 * bond / math.sqrt(3.0)
        side = 5 * lattice
        fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
        xyz = []
        for i in range(5):
            for j in range(5):
                for k in range(5):
                    for b in range(8):
                        f = b % 4
                        shift = 0.0 if b < 4 else 0.25
                        xyz.append(((i + fcc[f][0] + shift) * lattice) % side)
                        xyz.append(((j + fcc[f][1] + shift) * lattice) % side)
                        xyz.append(((k + fcc[f][2] + shift) * lattice) % side)
        n = len(xyz) // 3
        arr = (ctypes.c_double * len(xyz))(*xyz)
        box = (ctypes.c_double * 9)(0, 0, 0, side, side, side, 0, 0, 0)
        types = (ctypes.c_int * n)(*([1] * n))
        ids = (ctypes.c_int * n)(*range(1, n + 1))
        frame = lib.seams_frame_from_arrays(arr, n, box, types, ids)
        if not frame:
            raise RuntimeError("could not build the ice frame")
        labels = (ctypes.c_int * n)()
        counts = (ctypes.c_int * 9)()
        self._lib = lib
        self._frame = frame
        self._labels = labels
        self._counts = counts

    def time_published_chill_plus(self):
        status = self._lib.seams_frame_chill_plus(
            self._frame, 3.5, 1, b"cutoff", self._labels, self._counts)
        if status != 0:
            raise RuntimeError("chill+ failed")

    def teardown(self):
        if getattr(self, "_frame", None):
            self._lib.seams_frame_free(self._frame)


class LocalQl:
    """Local q4 and q6 on 110592 molecules of ideal ice Ic.

    Eight threads, cutoff 3.5 angstrom. The neighbour list is inside the
    call. The averaged buffer is null, so the Lechner-Dellago pass is not
    in the sample. Setup builds the frame and warms the scratch once.
    """

    def setup(self):
        import ctypes
        import math
        os.environ["OMP_NUM_THREADS"] = "8"
        os.environ["OMP_DYNAMIC"] = "FALSE"
        os.environ["OMP_PROC_BIND"] = "false"
        os.environ["OMP_WAIT_POLICY"] = "active"
        try:
            os.sched_setaffinity(0, set(range(8)))
        except (AttributeError, OSError):
            pass
        import sys
        sys.path.insert(0, str(ROOT / "scripts"))
        import release_build
        lib_path, self.build = release_build.ensure_release(
            os.environ.get("SEAMS_LANG_LIB"))
        lib = ctypes.CDLL(lib_path)
        lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_steinhardt.restype = ctypes.c_int
        lib.seams_frame_steinhardt.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
        bond = 2.75
        lattice = 4.0 * bond / math.sqrt(3.0)
        cells = 24
        side = cells * lattice
        fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
        xyz = []
        for i in range(cells):
            for j in range(cells):
                for k in range(cells):
                    for b in range(8):
                        f = b % 4
                        shift = 0.0 if b < 4 else 0.25
                        xyz.append(((i + fcc[f][0] + shift) * lattice) % side)
                        xyz.append(((j + fcc[f][1] + shift) * lattice) % side)
                        xyz.append(((k + fcc[f][2] + shift) * lattice) % side)
        n = len(xyz) // 3
        arr = (ctypes.c_double * len(xyz))(*xyz)
        box = (ctypes.c_double * 9)(0, 0, 0, side, side, side, 0, 0, 0)
        types = (ctypes.c_int * n)(*([1] * n))
        ids = (ctypes.c_int * n)(*range(1, n + 1))
        frame = lib.seams_frame_from_arrays(arr, n, box, types, ids)
        if not frame:
            raise RuntimeError("could not build the ice frame")
        ql = (ctypes.c_double * n)()
        if lib.seams_frame_steinhardt(frame, 3.5, 1, 4, ql, None) != 0:
            raise RuntimeError("warmup failed")
        self._lib = lib
        self._frame = frame
        self._ql = ql
        self._n = n

    def time_q4(self):
        if self._lib.seams_frame_steinhardt(
                self._frame, 3.5, 1, 4, self._ql, None) != 0:
            raise RuntimeError("q4 failed")

    def time_q6(self):
        if self._lib.seams_frame_steinhardt(
                self._frame, 3.5, 1, 6, self._ql, None) != 0:
            raise RuntimeError("q6 failed")

    def teardown(self):
        if getattr(self, "_frame", None):
            self._lib.seams_frame_free(self._frame)
            self._frame = None


def _scale_plan():
    """The twelve-call campaign script is not part of this tree."""
    import importlib.util
    path = ROOT / "scripts" / "scale_campaign.py"
    if not path.is_file():
        return None
    spec = importlib.util.spec_from_file_location("scale_campaign", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_SCALE = _scale_plan()


class CutoffIce:
    """Published CHILL+, q6, primitive rings and the cutoff neighbour list.

    Sizes are 1000, 10648 and 110592 molecules of ideal ice Ic. One
    thread. The library must be a meson release -O3 build; a debug
    library is refused before any timed call.
    """

    params = [1000, 10648, 110592]
    param_names = ["molecules"]

    def setup(self, molecules):
        import ctypes
        import math
        import sys
        sys.path.insert(0, str(ROOT / "scripts"))
        import release_build
        lib_path, self.build = release_build.ensure_release(
            os.environ.get("SEAMS_LANG_LIB"))
        cells = {1000: 5, 10648: 11, 110592: 24}
        if molecules not in cells:
            raise ValueError(molecules)
        m = cells[molecules]
        lib = ctypes.CDLL(lib_path)
        lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_chill_plus.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_steinhardt.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
        lib.seams_frame_rings.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int]
        lib.seams_frame_neighbours.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
        bond = 2.75
        lattice = 4.0 * bond / math.sqrt(3.0)
        side = m * lattice
        fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
        xyz = []
        for i in range(m):
            for j in range(m):
                for k in range(m):
                    for b in range(8):
                        f = b % 4
                        shift = 0.0 if b < 4 else 0.25
                        xyz.append(((i + fcc[f][0] + shift) * lattice) % side)
                        xyz.append(((j + fcc[f][1] + shift) * lattice) % side)
                        xyz.append(((k + fcc[f][2] + shift) * lattice) % side)
        n = len(xyz) // 3
        arr = (ctypes.c_double * len(xyz))(*xyz)
        box = (ctypes.c_double * 9)(0, 0, 0, side, side, side, 0, 0, 0)
        types = (ctypes.c_int * n)(*([1] * n))
        ids = (ctypes.c_int * n)(*range(1, n + 1))
        frame = lib.seams_frame_from_arrays(arr, n, box, types, ids)
        if not frame:
            raise RuntimeError("could not build the ice frame")
        self._lib = lib
        self._frame = frame
        self._n = n
        self._labels = (ctypes.c_int * n)()
        self._counts9 = (ctypes.c_int * 9)()
        self._ql = (ctypes.c_double * n)()
        self._qlbar = (ctypes.c_double * n)()
        self._rings = (ctypes.c_int * 7)()

    def time_published_chill_plus(self, molecules):
        status = self._lib.seams_frame_chill_plus(
            self._frame, 3.5, 1, b"cutoff", self._labels, self._counts9)
        if status != 0:
            raise RuntimeError("chill+ failed")

    def time_q6(self, molecules):
        status = self._lib.seams_frame_steinhardt(
            self._frame, 3.5, 1, 6, self._ql, self._qlbar)
        if status != 0:
            raise RuntimeError("q6 failed")

    def time_rings(self, molecules):
        status = self._lib.seams_frame_rings(
            self._frame, 3.5, 1, b"cutoff", 6, self._rings, 7)
        if status != 0:
            raise RuntimeError("rings failed")

    def time_neighbours(self, molecules):
        status = self._lib.seams_frame_neighbours(
            self._frame, 3.5, 1, b"cutoff")
        if status != 0:
            raise RuntimeError("neighbours failed")

    def teardown(self, molecules):
        if getattr(self, "_frame", None):
            self._lib.seams_frame_free(self._frame)
            self._frame = None


if _SCALE is not None:
  class ScaleProtocol:
    """Twelve-call strong and weak points. Airspeed Velocity runs each once."""

    # One axis. A list of points is not a Cartesian product of those points.
    params = [[
        (point["kind"], point["model"], point["analysis"], point["threads"],
         point["reps"][0], point["reps"][1], point["reps"][2])
        for point in _SCALE.points()
    ]]
    param_names = ["point"]

    def setup(self, point):
        kind, model, analysis, threads, nx, ny, nz = point
        spec = _SCALE.MODELS[model]
        frame = _SCALE.dump_path(model, (nx, ny, nz))
        if not frame.is_file():
            _SCALE.prepare_dump(spec, model, (nx, ny, nz))
        cpus = "0" if threads == 1 else f"0-{threads - 1}"
        self.cmd = [
            "taskset", "-c", cpus, str(_SCALE.BENCH), "run",
            "--traj", str(frame), "--type", spec["type"],
            "--analysis", analysis, "--threads", str(threads),
            "--reps", str(_SCALE.REPS), "--out-dir", "bbdir/scale/asv-tmp",
        ]
        self.log_name = (
            f"{kind}-{model}-{analysis}-t{threads}-{nx}x{ny}x{nz}.log"
        )
        self.env = os.environ.copy()
        self.env.update({
            "OMP_NUM_THREADS": str(threads),
            "OMP_DYNAMIC": "FALSE",
            "OMP_PROC_BIND": "close",
            "RAYON_NUM_THREADS": str(threads),
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        })

    def time_twelve_calls(self, point):
        head = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        tree = subprocess.check_output(
            ["git", "rev-parse", f"{head}^{{tree}}"], cwd=ROOT, text=True).strip()
        shown = " ".join(self.cmd)
        done = subprocess.run(
            self.cmd, cwd=ROOT, env=self.env, check=False,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        log_dir = ROOT / "bbdir" / "scale" / "asv-logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        text = "\n".join([
            f"HEAD {head}",
            f"TREE {tree}",
            f"COMMAND {shown}",
            done.stdout.rstrip("\n"),
            f"EXIT {done.returncode}",
            "",
        ])
        (log_dir / self.log_name).write_text(text, encoding="utf-8")
        if done.returncode != 0:
            raise subprocess.CalledProcessError(done.returncode, self.cmd, done.stdout)

  ScaleProtocol.time_twelve_calls.repeat = 1
  ScaleProtocol.time_twelve_calls.number = 1
  ScaleProtocol.time_twelve_calls.warmup_time = 0


class LinkcellPartition:
    """One rank of the linkcell bisection, neighbours on the 24-cell ice."""

    params = [1, 2, 4, 8]
    param_names = ["threads"]

    def setup(self, threads):
        bench = ROOT / "build-mpi" / "tests" / "rank_bench"
        if not bench.is_file():
            raise NotImplementedError("build-mpi/tests/rank_bench is not built")
        self.bench = bench
        self.threads = threads

    def time_neighbours(self, threads):
        env = os.environ.copy()
        env.update({
            "OMP_NUM_THREADS": str(threads),
            "OMP_DYNAMIC": "FALSE",
            "OMP_PROC_BIND": "close",
            "RAYON_NUM_THREADS": "1",
            "SEAMS_SCALE": "1",
        })
        subprocess.run(
            ["mpirun", "--bind-to", "core", "--map-by", f"node:PE={threads}",
             "-n", "1", str(self.bench), "24", "2", "neigh"],
            cwd=ROOT, env=env, check=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )


class RingsOneThread:
    """ringNetwork against OVITO FindRings. The cutoff graph is built in setup.

    One thread. A missing bench binary or an interpreter without OVITO skips
    the pair. Samples are one search each; the published medians are the
    eight-call session in markdown/rings-one-thread.json.
    """

    params = ["tri-136", "orth-4k", "orth-1e5"]
    param_names = ["frame"]

    def setup(self, frame):
        import importlib.util
        import tempfile

        bench = os.environ.get(
            "SEAMS_RING_BENCH",
            str(ROOT / "bbdir-release" / "tests" / "bench_rings_h2h"),
        )
        if not Path(bench).is_file():
            raise NotImplementedError("bench_rings_h2h is not built")
        os.environ["OVITO_THREAD_COUNT"] = "1"
        try:
            import ovito  # noqa: F401
        except Exception as exc:
            raise NotImplementedError("OVITO is not importable") from exc
        spec = importlib.util.spec_from_file_location(
            "rings_one_thread", ROOT / "scripts" / "rings_one_thread.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self._module = module
        self._work = tempfile.TemporaryDirectory(prefix="rings-h2h-")
        frames = {
            name: traj
            for name, traj, _count in module.prepare_frames(
                module.load_h2h(), Path(self._work.name))
        }
        self._dseams = module.DseamsRepl(bench, frames[frame], 1)
        self._ovito = module.OvitoRepl(sys.executable, frames[frame], 1)
        if self._dseams.hist != self._ovito.hist:
            raise RuntimeError(
                f"{frame} census d-SEAMS {self._dseams.hist} "
                f"OVITO {self._ovito.hist}")

    def time_dseams(self, frame):
        self._dseams.sample()

    def time_ovito_findrings(self, frame):
        self._ovito.sample()

    def teardown(self, frame):
        for arm in ("_dseams", "_ovito"):
            worker = getattr(self, arm, None)
            if worker is not None:
                worker.close()
        work = getattr(self, "_work", None)
        if work is not None:
            work.cleanup()


class ClassifyIce:
    """Ideal ice Ic structure labels against the release library.

    110592 molecules, one and eight threads. The OVITO polyhedral
    template match of the same frame, with cubic and hexagonal diamond
    enabled, is the sitting in benchmarks/classify-carry-sitting.txt.
    """

    params = [1, 8]
    param_names = ["threads"]

    def setup(self, threads):
        bench = ROOT / "bbdir-release" / "tests" / "bench_shell_scale"
        if not bench.is_file():
            raise NotImplementedError("release bench_shell_scale is not built")
        cpus = "0" if threads == 1 else "0-7"
        self.cmd = [
            "taskset", "-c", cpus,
            str(bench), str(threads),
        ]
        self.env = os.environ.copy()
        self.env["OMP_NUM_THREADS"] = str(threads)
        self.env["OMP_DYNAMIC"] = "FALSE"

    def time_classify(self, threads):
        subprocess.run(self.cmd, env=self.env, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


class HopStride:
    """One carried ring update on a production-stride frame.

    model mw is the 4096-site growth window. Its strides are 10, 100 and
    1000 steps of 10 fs, which is 0.1, 1 and 10 ps. The dumps are written
    by input/hop/regenerate.py.

    model tip4p is the 4096-molecule TIP4P/2005 ice/water slab. Its
    strides are 5, 50 and 500 steps of 2 fs, which is 0.01, 0.1 and 1 ps.
    The dumps are written by input/hop/regenerate_tip4p.py. A stride that
    does not belong to the model is skipped.

    The dumps are not stored in git. The timed command is the same driver
    as the Snakemake rule, limited to the opening transition. The
    trajectory medians live in benchmarks/hop-stride-*.txt.
    """

    params = (["mw", "tip4p"], [5, 10, 50, 100, 500, 1000], ["cutoff", "seeded"])
    param_names = ["model", "stride", "graph"]
    _steps = {"mw": {10, 100, 1000}, "tip4p": {5, 50, 500}}

    def setup(self, model, stride, graph):
        if stride not in self._steps.get(model, set()):
            raise NotImplementedError(f"{model} has no {stride}-step stride")
        if model == "mw":
            dump = ROOT / "data" / "hop" / f"stride-{stride}.lammpstrj"
            missing = "production-stride dump is produced by input/hop/regenerate.py"
        else:
            dump = ROOT / "data" / "hop" / f"tip4p-stride-{stride}.lammpstrj"
            missing = "TIP4P dump is produced by input/hop/regenerate_tip4p.py"
        bench = ROOT / "build" / "tests" / "bench_hop_stride"
        if not dump.is_file() or not bench.is_file():
            raise NotImplementedError(missing)
        self.cmd = CORE + [
            str(bench), "time", str(dump),
            "--graph", graph,
            "--reps", "1",
            "--warmup", "0",
            "--last", "2",
            "--stride", str(stride),
            "--type", "1",
        ]

    def time_carry(self, model, stride, graph):
        _run(self.cmd)


class OmpScale:
    """C API knn, q6, CHILL and a dump read at one and eight threads.

    The timed body is eight samples in this process. The Snakemake rule
    writes the same session and records the core count, not the host name.
    """

    params = [1, 8]
    param_names = ["threads"]

    def setup(self, threads):
        import importlib.util

        script = ROOT / "scripts" / "omp_scale_time.py"
        lib_path = ROOT / "bbdir-release" / "src" / "libseams.so"
        if not lib_path.is_file():
            raise NotImplementedError("bbdir-release/src/libseams.so is not built")
        spec = importlib.util.spec_from_file_location("omp_scale_time", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self._module = module
        self._lib = module.load_library(lib_path)

    def time_session(self, threads):
        self._module.session(self._lib, threads, samples=8)


OmpScale.time_session.repeat = 1
OmpScale.time_session.number = 1


class H2HSession:
    """d-SEAMS, freud, and OVITO in one process.

    Eight samples of knn, q6, and read at 4096, 32768, and 110592
    molecules. The Snakemake rule writes the raw log and the core count.
    """

    def setup(self):
        import importlib.util

        script = ROOT / "scripts" / "h2h_session.py"
        lib_path = ROOT / "bbdir-release" / "src" / "libseams.so"
        if not lib_path.is_file():
            raise NotImplementedError("bbdir-release/src/libseams.so is not built")
        spec = importlib.util.spec_from_file_location("h2h_session", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self._module = module
        os.environ.setdefault("OVITO_THREAD_COUNT", "8")
        os.environ.setdefault("SEAMS_LIB", str(lib_path))
        try:
            import freud  # noqa: F401
            import ovito  # noqa: F401
        except Exception as exc:
            raise NotImplementedError("freud and ovito are not importable") from exc

    def time_session(self):
        self._module.run_session()


H2HSession.time_session.repeat = 1
H2HSession.time_session.number = 1


class CageHeadToHead:
    """GenIce sI unit cell against cycless and ClavisCage.

    The scaled sitting, including sII, sH and the cells past 110592
    molecules, is benchmarks/cage-h2h-v8.txt. cycless drops spanning
    cycles. ClavisCage is TRACE at 0.35 nm. d-SEAMS is seams_frame_cages
    and one cmdCages signature search. Below 2048 molecules the library
    keeps a single thread.
    """

    params = [1, 8]
    param_names = ["threads"]

    def setup(self, threads):
        import importlib.util

        script = ROOT / "scripts" / "cage_h2h_v8.py"
        spec = importlib.util.spec_from_file_location("cage_h2h_v8", script)
        if spec is None or spec.loader is None:
            raise NotImplementedError("cage head-to-head script is missing")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        try:
            pos, cell = module.generate("sI", (1, 1, 1))
        except Exception as exc:
            raise NotImplementedError("genice2 is required") from exc
        dump = ROOT / "benchmarks" / "genice-frames" / "asv-sI.lammpstrj"
        module.write_dump(dump, pos, cell)
        os.environ["OMP_NUM_THREADS"] = str(threads)
        os.environ["OMP_DYNAMIC"] = "FALSE"
        os.environ["RAYON_NUM_THREADS"] = str(threads)
        try:
            self.bundle = module.prepare(dump)
        except Exception as exc:
            raise NotImplementedError("the release language library is required") from exc
        self.mod = module
        self.dump = dump
        self._cycless = importlib.util.find_spec("cycless") is not None
        self._clavis = importlib.util.find_spec("clavis_engine") is not None

    def time_dseams_cages(self, threads):
        self.bundle["dseams"].seeded()

    def time_dseams_cmd(self, threads):
        self.mod.cmd_signature(self.dump)

    def time_cycless(self, threads):
        if not self._cycless:
            raise NotImplementedError("cycless is not installed")
        from cycless.cycles import cycles_iter
        from cycless.polyhed import polyhedra_iter

        rings = [
            [int(node) for node in ring]
            for ring in cycles_iter(self.bundle["graph"], 8, pos=self.bundle["frac"])
        ]
        list(polyhedra_iter(rings, 22, quick=True))

    def time_clavis(self, threads):
        if not self._clavis:
            raise NotImplementedError("clavis_engine is not installed")
        import clavis_engine

        clavis_engine.analyze_waters(
            self.bundle["waters"], "[]", self.bundle["box"], "trace", 0.35, 35.0
        )

    def teardown(self, threads):
        bundle = getattr(self, "bundle", None)
        if bundle is not None:
            bundle["dseams"].close()


class H2HMdapyPyscal3Vesin:
    """d-SEAMS, mdapy, pyscal3, vesin, and OVITO in one process.

    Eight samples of a 3.5 A neighbour list, k=4 neighbours, q6 and q6bar,
    and both CHILL+ definitions at 4096, 32768, and 110592 molecules.
    paper averages q_lm over the four closest neighbours. ovito averages
    q_lm over every neighbour inside the cutoff and is checked atom for
    atom against mdapy and the OVITO modifier. The Snakemake rule writes
    the raw log and the core count.
    """

    def setup(self):
        import importlib.util

        script = ROOT / "scripts" / "h2h_mdapy_pyscal3_vesin.py"
        lib_path = ROOT / "bbdir-release" / "src" / "libseams.so"
        if not lib_path.is_file():
            raise NotImplementedError("bbdir-release/src/libseams.so is not built")
        spec = importlib.util.spec_from_file_location("h2h_mdapy_pyscal3_vesin", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self._module = module
        os.environ.setdefault("SEAMS_LIB", str(lib_path))
        try:
            import mdapy  # noqa: F401
            import pyscal3  # noqa: F401
            import vesin  # noqa: F401
            import ovito  # noqa: F401
        except Exception as exc:
            raise NotImplementedError("mdapy, pyscal3, vesin, and ovito are not importable") from exc

    def time_session(self):
        self._module.run_session()


H2HMdapyPyscal3Vesin.time_session.repeat = 1
H2HMdapyPyscal3Vesin.time_session.number = 1


class H2HFreudOvito:
    """d-SEAMS, freud, and OVITO in one process.

    Eight samples of a 3.5 A neighbour list, k=4 neighbours, q3 and q6,
    ovito-mode CHILL+, and the diamond shell at 4096, 32768, and 110592
    molecules, plus the mixed frames. The Snakemake rule writes the raw
    log and the core count.
    """

    def setup(self):
        import importlib.util

        script = ROOT / "scripts" / "h2h_freud_ovito.py"
        lib_path = ROOT / "bbdir-release" / "src" / "libseams.so"
        if not lib_path.is_file():
            raise NotImplementedError("bbdir-release/src/libseams.so is not built")
        spec = importlib.util.spec_from_file_location("h2h_freud_ovito", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self._module = module
        os.environ.setdefault("SEAMS_LIB", str(lib_path))
        try:
            import freud  # noqa: F401
            import ovito  # noqa: F401
        except Exception as exc:
            raise NotImplementedError("freud and ovito are not importable") from exc

    def time_session(self):
        self._module.run_session()


H2HFreudOvito.time_session.repeat = 1
H2HFreudOvito.time_session.number = 1


class H2HRings:
    """Primitive rings against rings-statistics-matter, cycless, and OVITO.

    Eight samples at one thread and eight threads. The Snakemake rule
    writes the raw log.
    """

    def setup(self):
        import importlib.util

        script = ROOT / "scripts" / "h2h_rings.py"
        lib_path = ROOT / "bbdir-release" / "src" / "libseams.so"
        if not lib_path.is_file():
            raise NotImplementedError("bbdir-release/src/libseams.so is not built")
        spec = importlib.util.spec_from_file_location("h2h_rings", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self._module = module
        os.environ.setdefault("SEAMS_LIB", str(lib_path))
        try:
            import cycless  # noqa: F401
            import ovito  # noqa: F401
            import rings_statistics_matter  # noqa: F401
        except Exception as exc:
            raise NotImplementedError("ring libraries are not importable") from exc

    def time_session(self):
        self._module.run_session()


H2HRings.time_session.repeat = 1
H2HRings.time_session.number = 1
