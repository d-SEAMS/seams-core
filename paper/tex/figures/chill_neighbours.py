"""Split the d-SEAMS and OVITO CHILL+ labels of one frame by coordination.

Run with a Python that has ovito 3.17.2 and numpy:

    python chill_neighbours.py BENCH DUMP > data/h2h/chill-neighbours.txt

BENCH is the tests/bench_paper driver of a release build and DUMP the
liquid LAMMPS dump of the head-to-head campaign.
"""

import argparse
import importlib.util
import pathlib
import subprocess
import tempfile
from collections import Counter

import numpy as np
from ovito.modifiers import ChillPlusModifier

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "h2h", HERE.parents[2] / "scripts" / "h2h_campaign.py"
)
h2h = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h2h)

# The d-SEAMS labels come from one call of the driver without a warm-up call.
LABEL_WARMUP = 0
LABEL_CALLS = 1
NATIVE = {
    "cubic": "cubic", "reCubic": "cubic", "hexagonal": "hexagonal",
    "reHex": "hexagonal", "interfacial": "interfacial", "clathrate": "hydrate",
    "interClathrate": "interfacial-hydrate", "mixed": "mixed",
    "water": "other", "unclassified": "other",
}


def coordination(positions: np.ndarray, box: np.ndarray, cutoff: float):
    rows = []
    for i in range(len(positions)):
        delta = positions[i] - positions
        delta -= box * np.round(delta / box)
        distance = np.sqrt((delta**2).sum(axis=1))
        rows.append(np.flatnonzero((distance > 0) & (distance < cutoff)))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("bench", type=pathlib.Path)
    parser.add_argument("dump", type=pathlib.Path)
    args = parser.parse_args()

    positions, box = h2h.read_lammps(args.dump, frame=h2h.SOURCE_FRAME)
    with tempfile.TemporaryDirectory() as work:
        trajectory = pathlib.Path(work) / "frame.lammpstrj"
        sidecar = pathlib.Path(work) / "labels.txt"
        h2h.write_lammps(trajectory, positions, box)
        subprocess.run(
            [
                str(args.bench),
                str(trajectory),
                str(h2h.DRIVER_FRAME),
                str(LABEL_WARMUP),
                str(LABEL_CALLS),
                "chill",
                str(sidecar),
            ],
            check=True,
            capture_output=True,
        )
        native = sidecar.read_text().split()
    dseams = np.asarray([h2h.LABEL_ID[NATIVE[item]] for item in native])

    data = h2h.ovito_data(positions, box)
    data.apply(ChillPlusModifier(cutoff=h2h.CUTOFF_A))
    structure = np.asarray(data.particles["Structure Type"])
    ovito = np.array(h2h.OVITO_CHILL_LABELS)[structure]

    # CHILL+ bonds each molecule to its CHILL_NEIGHBOURS nearest neighbours,
    # so the split is between molecules with that many neighbours in the
    # cutoff and the rest.
    shells = coordination(positions, box, h2h.CUTOFF_A)
    count = np.array([len(row) for row in shells])
    four = count == h2h.CHILL_NEIGHBOURS
    closed = four & np.array(
        [bool(np.all(count[row] == h2h.CHILL_NEIGHBOURS)) for row in shells]
    )
    agree = dseams == ovito
    other = h2h.LABEL_ID["other"]

    print("molecules", len(positions))
    print("agree", int(agree.sum()))
    print("disagree", int((~agree).sum()))
    print("four", int(four.sum()))
    print("notfour", int((~four).sum()))
    print("notfour_ovito_other", int((ovito[~four] == other).sum()))
    print("notfour_dseams_labelled", int((dseams[~four] != other).sum()))
    print("notfour_disagree", int((~agree & ~four).sum()))
    print("four_disagree", int((~agree & four).sum()))
    print("closed", int(closed.sum()))
    print("closed_disagree", int((~agree & closed).sum()))
    for name, labels in (("dseams", dseams), ("ovito", ovito)):
        histogram = Counter(labels.tolist())
        for index, label in enumerate(h2h.LABELS):
            print(f"{name}_{label}", histogram.get(index, 0))
    for value, number in sorted(Counter(count.tolist()).items()):
        print(f"coordination_{value}", number)


if __name__ == "__main__":
    main()
