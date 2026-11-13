#!/usr/bin/env python3
"""Workers for the classify/carry sitting. stdout is the log lines."""

from __future__ import annotations

import ctypes
import math
import os
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WARMUP = 1
REPEATS = 8


def median(samples):
    ordered = sorted(samples)
    n = len(ordered)
    if n % 2:
        return ordered[n // 2]
    return 0.5 * (ordered[n // 2 - 1] + ordered[n // 2])


def iqr(samples):
    ordered = sorted(samples)
    def q(p):
        pos = p * (len(ordered) - 1)
        lo = int(math.floor(pos))
        hi = int(math.ceil(pos))
        frac = pos - lo
        if hi >= len(ordered):
            return ordered[-1]
        return ordered[lo] * (1.0 - frac) + ordered[hi] * frac
    return q(0.75) - q(0.25)


def ice_ic(cells, bond=2.75):
    import numpy as np
    lattice = 4.0 * bond / math.sqrt(3.0)
    side = cells * lattice
    fcc = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
    basis = np.vstack([fcc, fcc + 0.25])
    ijk = np.stack(np.meshgrid(*[np.arange(cells)] * 3, indexing="ij"), -1).reshape(-1, 3)
    pts = (ijk[:, None, :] + basis).reshape(-1, 3) * lattice
    return np.ascontiguousarray(pts % side, dtype=np.float64), float(side)


def symmetry(lib_path, threads):
    import numpy as np
    lib = ctypes.CDLL(lib_path)
    fn = lib.seams_symmetry_classify_tol
    fn.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.c_double,
        ctypes.POINTER(ctypes.c_int)]
    fn.restype = ctypes.c_int
    pts, side = ice_ic(24)
    n = len(pts)
    lat = np.array([side, 0, 0, 0, side, 0, 0, 0, side], dtype=np.float64)
    labels = np.zeros(n, dtype=np.int32)
    rng = np.random.default_rng(20261113)
    frames = {
        "ideal-Ic": pts,
        "gaussian-Ic": np.ascontiguousarray(pts + rng.normal(0.0, 0.05, pts.shape)),
    }
    crystal = lib.seams_symmetry_crystal
    crystal.argtypes = [
        ctypes.c_int, ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
    crystal.restype = ctypes.c_int
    n_out = ctypes.c_int()
    if crystal(0, None, None, 0, ctypes.byref(n_out)) != 0:
        raise RuntimeError("ice Ih cell")
    basis = np.zeros(n_out.value * 3, dtype=np.float64)
    cell_lat = np.zeros(9, dtype=np.float64)
    if crystal(0, cell_lat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
               basis.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
               n_out.value, ctypes.byref(n_out)) != 0:
        raise RuntimeError("ice Ih coordinates")
    basis = basis.reshape(-1, 3)
    vecs = cell_lat.reshape(3, 3)
    reps = 12
    ih = []
    for i in range(reps):
        for j in range(reps):
            for k in range(reps):
                shift = i * vecs[0] + j * vecs[1] + k * vecs[2]
                ih.append(basis + shift)
    frames["ideal-Ih"] = np.ascontiguousarray(np.vstack(ih))
    ih_lat = (vecs * reps).reshape(-1).copy()

    def classify(pos, lattice):
        arr = np.ascontiguousarray(pos, dtype=np.float64)
        lab = np.zeros(len(arr), dtype=np.int32)
        code = fn(
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(arr),
            lattice.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            0.12, lab.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
        if code != 0:
            raise RuntimeError("classify failed")
        return lab

    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.modifiers import PolyhedralTemplateMatchingModifier
    from ovito.pipeline import Pipeline, StaticSource

    def ptm(pos, lengths):
        data = DataCollection()
        cell = SimulationCell(pbc=(True, True, True))
        matrix = np.zeros((3, 4))
        matrix[0, 0], matrix[1, 1], matrix[2, 2] = lengths
        cell[:, :] = matrix
        parts = Particles()
        parts.create_property("Position", data=np.ascontiguousarray(pos))
        data.objects.append(cell)
        data.objects.append(parts)
        pipe = Pipeline(source=StaticSource(data=data))
        mod = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1)
        for structure in mod.structures:
            structure.enabled = structure.id in (6, 7)
        pipe.modifiers.append(mod)
        out = pipe.compute()
        return np.asarray(out.particles["Structure Type"])

    import time
    for name, pos in frames.items():
        lattice = lat if name != "ideal-Ih" else ih_lat
        lengths = (side_of(lattice))
        lab = classify(pos, lattice)
        expect = 1 if "Ic" in name else 0
        n_ok = int(np.count_nonzero(lab == expect))
        print(f"labels frame={name} n={len(pos)} label={expect} count={n_ok}", flush=True)
        if name != "ideal-Ic":
            continue
        for kind, call in (("dseams", lambda: classify(pos, lattice)),
                           ("ptm", lambda: ptm(pos, lengths))):
            call()
            samples = []
            for _ in range(REPEATS):
                t0 = time.perf_counter()
                call()
                samples.append((time.perf_counter() - t0) * 1000.0)
            print(
                f"symmetry threads={threads} code={kind} molecules={len(pos)} "
                f"median_ms={median(samples):.4f} iqr_ms={iqr(samples):.4f} "
                f"samples_ms={' '.join(f'{v:.4f}' for v in samples)}",
                flush=True)
    # Phase split of one ideal call. Hardware counters are not available
    # in this guest; the library's own timer is the profile.
    os.environ["SEAMS_SYM_TIME"] = "1"
    classify(frames["ideal-Ic"], lat)


def side_of(lattice):
    return (float(lattice[0]), float(lattice[4]), float(lattice[8]))


def ql(binary):
    import numpy as np
    for cells, molecules in ((5, 1000), (11, 10648), (24, 110592)):
        pts, side = ice_ic(cells)
        assert len(pts) == molecules
        fd, name = tempfile.mkstemp(suffix=".bin")
        os.close(fd)
        frame = Path(name)
        with frame.open("wb") as handle:
            handle.write(struct.pack("<I", molecules))
            handle.write(struct.pack("<d", float(side)))
            handle.write(np.ascontiguousarray(pts, dtype=np.float64).tobytes())
        for order in (4, 6):
            proc = subprocess.run(
                [binary, str(frame), str(order), str(REPEATS)],
                check=True, text=True, capture_output=True)
            timed = [float(line) for line in proc.stdout.splitlines() if line.strip()]
            print(
                f"ql threads=8 molecules={molecules} l={order} "
                f"median_ms={median(timed):.4f} iqr_ms={iqr(timed):.4f} "
                f"samples_ms={' '.join(f'{v:.4f}' for v in timed)}",
                flush=True)
        frame.unlink(missing_ok=True)


def genice(lib_path):
    import json
    import numpy as np
    sys.path.insert(0, str(ROOT / "scripts"))
    import gen_ideal_ice
    import h2h_literature as lit
    os.environ.setdefault("SEAMS_LANG_LIB", str(ROOT / "bbdir-release" / "bindings" / "libseams_lang.so"))
    engine = lit.Engine(lib_path)
    import clavis_engine
    for frame in gen_ideal_ice.genice_frames(ROOT):
        text, box_nm, cutoff_nm, graph, note = lit.cage_setup(frame)
        summary, samples, count_512, others = engine.cages(frame, graph)
        raw = clavis_engine.analyze_waters(text, "[]", box_nm, "trace", cutoff_nm, 35.0)
        parsed = json.loads(raw)
        theirs = (
            len(parsed["cages_512"]),
            len(parsed["cages_51262"]),
            len(parsed["cages_51264"]),
        )
        ours = (count_512, others["51262"], others["51264"])
        print(
            f"cages phase={frame['phase']} molecules={frame['molecules']} "
            f"dseams={ours} clavis={theirs} identical={ours == theirs} "
            f"dseams_median_ms={summary['median_ms']:.4f} note={note}",
            flush=True)


def main():
    kind = sys.argv[1]
    if kind == "symmetry":
        symmetry(sys.argv[2], int(sys.argv[3]))
    elif kind == "ql":
        ql(sys.argv[2])
    elif kind == "genice":
        genice(sys.argv[2])
    else:
        raise SystemExit(kind)


if __name__ == "__main__":
    main()
