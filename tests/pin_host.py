#!/usr/bin/env python3
"""Placement respects the host and does not stack concurrent processes.

One fresh interpreter per check. OMP_PROC_BIND, OMP_PLACES,
OMP_WAIT_POLICY and GOMP_SPINCOUNT stay unset unless the check sets
them. Masks are read from this process while a call is in progress,
and again after it returns.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIB = Path(os.environ.get("SEAMS_LIB", ROOT / "bbdir-release/src/libyodaLib.so"))
sys.path.insert(0, str(ROOT / "scripts"))


def masks() -> list[str]:
    found = []
    task = Path("/proc/self/task")
    for tid in task.iterdir():
        try:
            text = (tid / "status").read_text()
        except OSError:
            continue
        for line in text.splitlines():
            if line.startswith("Cpus_allowed_list:"):
                found.append(line.split(":", 1)[1].strip())
                break
    return found


def child_env(**extra: str) -> dict[str, str]:
    env = dict(os.environ)
    for key in ("OMP_PROC_BIND", "OMP_PLACES", "OMP_WAIT_POLICY", "GOMP_SPINCOUNT",
                "OMP_DISPLAY_ENV", "SEAMS_PIN", "SEAMS_FAIR_WIDTH", "GOMP_CPU_AFFINITY",
                "KMP_AFFINITY"):
        env.pop(key, None)
    # Placement pins a lone team only when it fills the mask.
    env["OMP_NUM_THREADS"] = str(max(8, len(os.sched_getaffinity(0))))
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    lib_dir = str(LIB.resolve().parent)
    env["LD_LIBRARY_PATH"] = lib_dir + (os.pathsep + env["LD_LIBRARY_PATH"]
                                        if env.get("LD_LIBRARY_PATH") else "")
    env["SEAMS_LIB"] = str(LIB)
    env.update(extra)
    return env


def run_worker(code: str, env: dict[str, str], prefix: list[str] | None = None) -> str:
    cmd = (prefix or []) + [sys.executable, "-c", code]
    done = subprocess.run(cmd, env=env, capture_output=True, text=True, check=False)
    if done.returncode != 0:
        raise RuntimeError(done.stderr[-2000:] or done.stdout[-2000:])
    lines = [line for line in done.stdout.splitlines() if line.startswith("PIN ")]
    if not lines:
        raise RuntimeError(f"no PIN line\n{done.stderr[-1000:]}")
    return lines[-1]


WORKER = r"""
import os, sys, threading, time
sys.path.insert(0, os.environ["SEAMS_ROOT"] + "/scripts")
import h2h_protocol as proto
import numpy as np, ctypes
xyz, side = proto.ice_ic(24)
xyz = np.ascontiguousarray(xyz, np.float64)
box = np.array([0,0,0,side,side,side,0,0,0], np.float64)
n = len(xyz)
lib = ctypes.CDLL(os.environ["SEAMS_LIB"])
dp = ctypes.POINTER(ctypes.c_double)
ip = ctypes.POINTER(ctypes.c_int)
off = np.empty(n+1, np.int32)
cols = np.empty(16*n, np.int32)
count = ctypes.c_int()
lib.seams_cutoff_csr.argtypes = [dp, ctypes.c_int, dp, ctypes.c_double, ip, ip, ctypes.c_int, ip]
seen = []
stop = False
def sample():
    while not stop:
        cpus=set()
        for tid in os.listdir("/proc/self/task"):
            try:
                text=open(f"/proc/self/task/{tid}/status").read()
            except OSError:
                continue
            for line in text.splitlines():
                if line.startswith("Cpus_allowed_list:"):
                    cpus.add(line.split(":",1)[1].strip())
        seen.append(tuple(sorted(cpus)))
        time.sleep(0.00005)
def call():
    rc = lib.seams_cutoff_csr(xyz.ctypes.data_as(dp), n, box.ctypes.data_as(dp), 3.5,
        off.ctypes.data_as(ip), cols.ctypes.data_as(ip), len(cols), ctypes.byref(count))
    if rc != 0:
        raise SystemExit(f"rc {rc}")
gate = os.environ.get("SEAMS_WAIT_DIR")
if gate:
    os.makedirs(gate, exist_ok=True)
    open(os.path.join(gate, str(os.getpid())), "w").close()
    need = int(os.environ.get("SEAMS_WAIT_N", "1"))
    deadline = time.time() + 30
    while time.time() < deadline and len(os.listdir(gate)) < need:
        time.sleep(0.01)
threading.Thread(target=sample, daemon=True).start()
time.sleep(0.02)
seen.clear()
end = time.time() + 0.3
while time.time() < end:
    call()
stop = True
time.sleep(0.02)
after=set()
for tid in os.listdir("/proc/self/task"):
    text=open(f"/proc/self/task/{tid}/status").read()
    for line in text.splitlines():
        if line.startswith("Cpus_allowed_list:"):
            after.add(line.split(":",1)[1].strip())
# singletons observed while the call was running
singles=set()
for item in seen:
    if not isinstance(item, tuple):
        continue
    ones=[p for p in item if p.isdigit()]
    if len(ones) >= 2:
        singles.update(int(p) for p in ones)
print("PIN pid=%d singles=%s after=%s" % (os.getpid(), ",".join(str(c) for c in sorted(singles)), ",".join(sorted(after))), flush=True)
"""


def parse_pin(line: str) -> tuple[set[int], set[str]]:
    singles = set()
    after = set()
    for part in line.split():
        if part.startswith("singles="):
            body = part.split("=", 1)[1]
            if body:
                singles = {int(x) for x in body.split(",") if x}
        elif part.startswith("after="):
            body = part.split("=", 1)[1]
            if body:
                after = set(body.split(","))
    return singles, after


def one(env: dict[str, str], prefix: list[str] | None = None) -> tuple[set[int], set[str]]:
    env = dict(env)
    env["SEAMS_ROOT"] = str(ROOT)
    return parse_pin(run_worker(WORKER, env, prefix))


def one_placed(env: dict[str, str], prefix: list[str] | None = None) -> tuple[set[int], set[str]]:
    # A lone process beside other runnable threads is left unplaced on
    # purpose. Look again before failing.
    for _ in range(8):
        singles, after = one(env, prefix)
        if len(singles) >= 2:
            break
        time.sleep(0.5)
    return singles, after


def concurrent(n: int) -> list[set[int]]:
    env = child_env()
    env["SEAMS_ROOT"] = str(ROOT)
    gate = tempfile.mkdtemp()
    env["SEAMS_WAIT_DIR"] = gate
    env["SEAMS_WAIT_N"] = str(n)
    procs = []
    for _ in range(n):
        procs.append(subprocess.Popen(
            [sys.executable, "-c", WORKER], env=env,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))
    sets = []
    for proc in procs:
        out, err = proc.communicate(timeout=120)
        if proc.returncode != 0:
            raise RuntimeError(err[-2000:] or out[-2000:])
        line = [ln for ln in out.splitlines() if ln.startswith("PIN ")][-1]
        singles, _after = parse_pin(line)
        sets.append(singles)
        print(line, flush=True)
    shutil.rmtree(gate, ignore_errors=True)
    return sets


def assert_disjoint(sets: list[set[int]], label: str) -> None:
    for i in range(len(sets)):
        if not sets[i]:
            raise SystemExit(f"{label}: process {i} pinned no cpus {sets}")
        for j in range(i + 1, len(sets)):
            both = sets[i] & sets[j]
            if both:
                raise SystemExit(f"{label}: cpus {both} shared by {sets}")


def main() -> int:
    if not LIB.is_file():
        raise SystemExit(f"library missing: {LIB}")
    base = child_env()
    singles, after = one_placed(base)
    print(f"solo singles={sorted(singles)} after={sorted(after)}", flush=True)
    if len(singles) < 2:
        raise SystemExit(f"solo run did not pin a team: {singles}")
    if not any("-" in item or "," in item for item in after):
        raise SystemExit(f"affinity was not restored: {after}")

    if len(os.sched_getaffinity(0)) > 2:
        singles, after = one(child_env(OMP_NUM_THREADS="2"))
        print(f"2 threads singles={sorted(singles)} after={sorted(after)}", flush=True)
        if len(singles) >= 2:
            raise SystemExit(f"pinned a team narrower than the mask: {singles}")

    two = concurrent(2)
    assert_disjoint(two, "2 processes")
    # Four processes get two cpus each only with eight cpus in the mask.
    # On fewer, each process holds one cpu, which is not a pinned team.
    if len(os.sched_getaffinity(0)) >= 8:
        four = concurrent(4)
        assert_disjoint(four, "4 processes")
    else:
        print(f"4 processes: not run, {len(os.sched_getaffinity(0))} cpus in the mask", flush=True)

    # Narrower process mask. Eight requested threads must stay inside it.
    narrow_env = child_env()
    narrow_env["SEAMS_ROOT"] = str(ROOT)
    singles, after = one_placed(narrow_env, ["taskset", "-c", "0-3"])
    print(f"taskset 0-3 singles={sorted(singles)} after={sorted(after)}", flush=True)
    if not singles or not singles <= {0, 1, 2, 3}:
        raise SystemExit(f"left the taskset mask: {singles}")

    # Host already bound the caller to one cpu. Do not spread.
    bound = child_env()
    bound["SEAMS_ROOT"] = str(ROOT)
    host_cpu = max(os.sched_getaffinity(0))
    code = f"import os\nos.sched_setaffinity(0, {{{host_cpu}}})\n" + WORKER
    singles, _after = parse_pin(run_worker(code, bound))
    print(f"prebound singles={sorted(singles)}", flush=True)
    if singles - {host_cpu}:
        raise SystemExit(f"spread off the host cpu: {singles}")

    # Runtime placement is the host's. KMP_AFFINITY is not libgomp's, so
    # the opt-out is "do not pin"; masks stay wider than one cpu.
    kmp = child_env(KMP_AFFINITY="scatter")
    kmp["SEAMS_ROOT"] = str(ROOT)
    singles, after = one(kmp)
    print(f"KMP_AFFINITY singles={sorted(singles)} after={sorted(after)}", flush=True)
    if len(singles) >= 2:
        raise SystemExit(f"pinned despite KMP_AFFINITY: {singles}")

    places = child_env(OMP_PROC_BIND="close", OMP_PLACES="0,2")
    places["SEAMS_ROOT"] = str(ROOT)
    singles, _after = one(places)
    print(f"OMP_PLACES singles={sorted(singles)}", flush=True)
    if singles - {0, 2}:
        raise SystemExit(f"left OMP_PLACES: {singles}")

    gomp = child_env(GOMP_CPU_AFFINITY="1,3")
    gomp["SEAMS_ROOT"] = str(ROOT)
    singles, _after = one(gomp)
    print(f"GOMP_CPU_AFFINITY singles={sorted(singles)}", flush=True)
    if singles - {1, 3}:
        raise SystemExit(f"left GOMP_CPU_AFFINITY: {singles}")

    off = child_env(SEAMS_PIN="0")
    off["SEAMS_ROOT"] = str(ROOT)
    singles, after = one(off)
    print(f"SEAMS_PIN=0 singles={sorted(singles)} after={sorted(after)}", flush=True)
    if len(singles) >= 2:
        raise SystemExit(f"SEAMS_PIN=0 still pinned: {singles}")

    if shutil.which("mpicc") and shutil.which("mpirun"):
        mpi_check()
    else:
        print("mpi skipped", flush=True)
    print("pin host checks passed", flush=True)
    return 0


def mpi_check() -> None:
    src = textwrap.dedent(r"""
        #define _GNU_SOURCE
        #include <mpi.h>
        #include <sched.h>
        #include <stdio.h>
        #include <stdlib.h>
        #include <dlfcn.h>
        typedef int (*csr_fn)(const double *, int, const double *, double,
                              int *, int *, int, int *);
        int main(int argc, char **argv) {
          MPI_Init(&argc, &argv);
          int rank = 0;
          MPI_Comm_rank(MPI_COMM_WORLD, &rank);
          void *lib = dlopen(argv[1], RTLD_NOW);
          if (!lib) { fprintf(stderr, "%s\n", dlerror()); return 1; }
          csr_fn fn = (csr_fn)dlsym(lib, "seams_cutoff_csr");
          if (!fn) { fprintf(stderr, "no seams_cutoff_csr\n"); return 1; }
          const int n = 24;
          double xyz[72];
          for (int i = 0; i < n; ++i) {
            xyz[3*i] = (i % 4) * 1.1;
            xyz[3*i+1] = ((i / 4) % 3) * 1.1;
            xyz[3*i+2] = (i / 12) * 1.1;
          }
          double box[9] = {0,0,0, 8,8,8, 0,0,0};
          int off[25], cols[256], count = 0;
          cpu_set_t start;
          CPU_ZERO(&start);
          sched_getaffinity(0, sizeof(start), &start);
          for (int k = 0; k < 4; ++k) fn(xyz, n, box, 3.5, off, cols, 256, &count);
          cpu_set_t set;
          CPU_ZERO(&set);
          sched_getaffinity(0, sizeof(set), &set);
          printf("MPI rank=%d start=%d restored=%d same=%d\n", rank, CPU_COUNT(&start),
                 CPU_COUNT(&set), CPU_EQUAL(&start, &set) ? 1 : 0);
          fflush(stdout);
          MPI_Finalize();
          return 0;
        }
    """)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp)
        (path / "pin_mpi.c").write_text(src)
        exe = path / "pin_mpi"
        # mpicc calls gcc by its real path, so a cc wrapper on PATH never
        # sees it. mpirun and mpicc stay unmapped. Linking this file with
        # the sanitizer makes Open MPI's memory patcher hit an illegal
        # instruction inside ASan's mprotect, so the helper stays plain
        # and the runtime is exported only to the ranks.
        compile_env = os.environ.copy()
        compile_env.pop("LD_PRELOAD", None)
        cmd = ["mpicc", "-O2", str(path / "pin_mpi.c"), "-o", str(exe), "-ldl"]
        build = subprocess.run(
            cmd, env=compile_env, capture_output=True, text=True, check=False)
        if build.returncode != 0:
            print(build.stderr[-1500:])
            raise SystemExit("mpicc failed")
        env = child_env()
        env.pop("LD_PRELOAD", None)
        env["SEAMS_ROOT"] = str(ROOT)
        forward: list[str] = []
        preload = os.environ.get("ASAN_TEST_PRELOAD", "")
        if preload:
            # The patcher rewrites bytes with mprotect. ASan rejects that
            # call and the rank dies before seams_cutoff_csr. The component
            # is unused here. Leak checks stay off for these two ranks:
            # Open MPI and PMIx hold allocations until the process exits,
            # and those frames are not added to the suite suppression list.
            # Every other test still runs LeakSanitizer.
            env["OMPI_MCA_memory"] = "^patcher"
            forward = [
                "-x", "LD_PRELOAD=" + preload,
                "-x", "LSAN_OPTIONS=detect_leaks=0",
            ]
        def ranks(bind: str) -> list[tuple[int, int, bool]]:
            done = subprocess.run(
                ["mpirun", *forward, "--bind-to", bind, "-n", "2", str(exe), str(LIB)],
                env=env, capture_output=True, text=True, check=False)
            print(bind, done.stdout, flush=True)
            if done.returncode != 0:
                raise SystemExit(done.stderr[-1500:] or f"mpirun {bind} failed")
            rows = [line for line in done.stdout.splitlines() if line.startswith("MPI ")]
            if len(rows) != 2:
                raise SystemExit(f"expected 2 ranks, got {rows}")
            fields = [dict(p.split("=", 1) for p in line.split()[1:]) for line in rows]
            return [(int(f["start"]), int(f["restored"]), f["same"] == "1") for f in fields]
        wide = ranks("none")
        if any(restored < 2 for _start, restored, _same in wide):
            raise SystemExit(f"bind-to none did not restore a wide mask: {wide}")
        # mpirun binds a rank to one core, which is one or more hardware
        # threads. The call must leave that mask as it found it.
        tight = ranks("core")
        if any(not same for _start, _restored, same in tight):
            raise SystemExit(f"bind-to core was widened: {tight}")


if __name__ == "__main__":
    sys.exit(main())
