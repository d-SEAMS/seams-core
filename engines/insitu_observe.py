# Copyright (c) d-SEAMS developers
"""Eight-atom InSitu.observe check. Skips when the language library is absent."""

import os
import sys
from pathlib import Path

try:
    import numpy as np
except ImportError:
    print("numpy is not installed")
    sys.exit(77)

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root / "adapters" / "python"))
build = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "build"
libs = list(build.rglob("libseams_lang.so"))
if not libs:
    print("libseams_lang.so was not built")
    sys.exit(77)
os.environ["SEAMS_LANG_LIB"] = str(libs[0])
libdir = str(libs[0].parent)
os.environ["LD_LIBRARY_PATH"] = libdir + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")

from dseams_insitu.engine import InSitu

a = 2.75 * 4.0 / np.sqrt(3.0)
fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
basis = ((0, 0, 0), (0.25, 0.25, 0.25))
points = [
    ((shift[0] + site[0]) * a, (shift[1] + site[1]) * a, (shift[2] + site[2]) * a)
    for shift in fcc
    for site in basis
]
pos = np.asarray(points, dtype=np.float64)
if pos.shape != (8, 3):
    sys.exit("expected 8 atoms")
box = np.array([a, a, a], dtype=np.float64)
engine = InSitu()
engine.observe(0, pos, box, ids=np.arange(8))
labels = engine.frames[-1]["labels"]
if labels.shape != (8,) or not np.all(labels == 0):
    sys.exit(f"expected 8 cubic labels, got {labels}")
held = engine.frames[-1]["box"]
if held.shape != (9,) or not np.allclose(held[3:6], [a, a, a]) or np.any(held[:3]) or np.any(held[6:]):
    sys.exit(f"expected orthorhombic edges in slots 3-5, got {held}")
try:
    engine._hold(pos, np.ones(4))
except ValueError:
    pass
else:
    sys.exit("a length-4 box was accepted")
print("insitu atoms", labels.shape[0], "cubic", int(np.sum(labels == 0)))
