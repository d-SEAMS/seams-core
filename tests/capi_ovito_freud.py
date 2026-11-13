#!/usr/bin/env python3
"""Optional OVITO and freud q4/q6 for one frame file.

Prints nothing and exits 0 when either library is absent. The C suite
compares the printed means with the C API when the lines are present.
OVITO 3.17 is compared at 1e-6. freud's float32 q_l is compared at 1e-4.
"""

from __future__ import annotations

import sys

try:
    import freud
    import numpy as np
    from ovito.data import DataCollection
    from ovito.modifiers import SteinhardtParametersModifier
    from ovito.pipeline import Pipeline, StaticSource
except ImportError:
    sys.exit(0)


def load(path):
    text = open(path, encoding="ascii").read().split()
    n = int(text[0])
    box = [float(x) for x in text[1:10]]
    body = text[10:]
    xyz = np.empty((n, 3), dtype=np.float64)
    for i in range(n):
        xyz[i, 0] = float(body[4 * i])
        xyz[i, 1] = float(body[4 * i + 1])
        xyz[i, 2] = float(body[4 * i + 2])
    return xyz, box


def main():
    xyz, box = load(sys.argv[1])
    lengths = np.array(box[3:6], dtype=np.float64)
    data = DataCollection()
    mat = np.zeros((3, 4))
    mat[:, :3] = np.diag(lengths)
    mat[:, 3] = box[0:3]
    data.create_cell(mat, pbc=(True, True, True))
    parts = data.create_particles(count=len(xyz))
    parts.create_property("Position", data=xyz)
    pipe = Pipeline(source=StaticSource(data=data))
    mod = SteinhardtParametersModifier(cutoff=3.5, l_list=(4, 6))
    pipe.modifiers.append(mod)
    out = pipe.compute()
    q4 = np.asarray(out.particles["Steinhardt Order Parameter (q4)"], dtype=np.float64)
    q6 = np.asarray(out.particles["Steinhardt Order Parameter (q6)"], dtype=np.float64)
    print(f"ovito_q4 {float(np.mean(q4)):.17g}")
    print(f"ovito_q6 {float(np.mean(q6)):.17g}")
    box_fr = freud.box.Box(Lx=float(lengths[0]), Ly=float(lengths[1]), Lz=float(lengths[2]))
    ql = freud.order.Steinhardt(l=list(range(4, 7, 2)))
    ql.compute((box_fr, xyz), neighbors={"r_max": 3.5})
    print(f"freud_q4 {float(np.mean(ql.particle_order[:, 0])):.17g}")
    print(f"freud_q6 {float(np.mean(ql.particle_order[:, 1])):.17g}")


if __name__ == "__main__":
    main()
