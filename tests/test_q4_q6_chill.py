"""The fused q4/q6/CHILL+ wrapper matches ideal ice, including q-bar."""

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bindings", "python"))

try:
    import numpy as np
except ImportError:
    sys.stderr.write("numpy is not installed\n")
    raise SystemExit(77)

from chill_plus import q4_q6, q4_q6_chill


def diamond(nrep=2):
    oo = 2.75
    a = oo * 4.0 / math.sqrt(3.0)
    basis = ((0, 0, 0), (0.25, 0.25, 0.25))
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    points = []
    for i in range(nrep):
        for j in range(nrep):
            for k in range(nrep):
                for shift in fcc:
                    for site in basis:
                        points.append(
                            (
                                (i + shift[0] + site[0]) * a,
                                (j + shift[1] + site[1]) * a,
                                (k + shift[2] + site[2]) * a,
                            )
                        )
    return np.asarray(points, dtype=np.float64), (nrep * a, nrep * a, nrep * a)


def main():
    pos, box = diamond()
    out = q4_q6_chill(pos, box)
    if not np.all(out["labels"] == 0):
        raise SystemExit("ideal ice is not all cubic")
    if abs(float(out["q4"][0]) - 0.5091750772) > 1e-8:
        raise SystemExit(f"q4 {out['q4'][0]}")
    if abs(float(out["q6"][0]) - 0.6285393611) > 1e-8:
        raise SystemExit(f"q6 {out['q6'][0]}")
    if np.max(np.abs(out["q4bar"] - out["q4"])) > 1e-12:
        raise SystemExit("q4-bar differs from q4 on ideal ice")
    if np.max(np.abs(out["q6bar"] - out["q6"])) > 1e-12:
        raise SystemExit("q6-bar differs from q6 on ideal ice")
    if int(out["counts"][0]) != len(pos):
        raise SystemExit("cubic count")
    fast = q4_q6(pos, box)
    if np.max(np.abs(fast["q4"] - out["q4"])) > 1e-12:
        raise SystemExit("q4_q6 q4 differs from the fused call")
    if np.max(np.abs(fast["q6"] - out["q6"])) > 1e-12:
        raise SystemExit("q4_q6 q6 differs from the fused call")
    if np.max(np.abs(fast["q4bar"] - out["q4bar"])) > 1e-12:
        raise SystemExit("q4_q6 q4bar differs from the fused call")
    if np.max(np.abs(fast["q6bar"] - out["q6bar"])) > 1e-12:
        raise SystemExit("q4_q6 q6bar differs from the fused call")
    local = q4_q6(pos, box, average=False)
    if "q4bar" in local or "q6bar" in local:
        raise SystemExit("average=False still returned q-bar")
    print("q4_q6_chill", len(pos), "cubic", int(out["counts"][0]))


if __name__ == "__main__":
    main()
