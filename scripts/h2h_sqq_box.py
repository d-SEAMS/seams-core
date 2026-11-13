#!/usr/bin/env python3
"""Show SQQ rejecting the triclinic sH cell.

The rejection is the library's own error. An orthorhombic sI dump is
read with the same configuration so the failure is the tilt, not the
atom style.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SI = ROOT / "input" / "traj" / "genice_sI.lammpstrj"
SH = ROOT / "input" / "traj" / "genice_sH.lammpstrj"


def main():
    from sqq.api import load_config
    from sqq.io.input.lammps import inspect_lammps_dump, normalize_lammps_config
    import sqq

    print(f"sqq {getattr(sqq, '__version__', 'unknown')}", flush=True)
    cfg = normalize_lammps_config(
        {"units": "real", "atom_style": "molecular", "type_map": {"1": ["SOL", "OW"]}}
    )
    try:
        inspect_lammps_dump(SI, cfg)
        print(f"accept {SI.name} orthorhombic dump", flush=True)
    except Exception as exc:
        print(f"reject {SI.name} {type(exc).__name__}: {exc}", flush=True)
        raise SystemExit(1)
    try:
        inspect_lammps_dump(SH, cfg)
        print(f"accept {SH.name}", flush=True)
        raise SystemExit("sH dump was accepted")
    except ValueError as exc:
        text = str(exc).replace(str(SH), SH.name)
        print(f"reject {SH.name} {type(exc).__name__}: {text}", flush=True)
    try:
        load_config({"pbc": {"box_mode": "triclinic"}})
        print("accept pbc.box_mode triclinic", flush=True)
        raise SystemExit("triclinic box mode was accepted")
    except Exception as exc:
        print(f"reject box_mode triclinic {type(exc).__name__}: {exc}", flush=True)


if __name__ == "__main__":
    main()
