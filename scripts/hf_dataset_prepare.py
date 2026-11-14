#!/usr/bin/env python3
"""Pack interface and nucleation frames for a Hugging Face dataset upload.

Writes a local directory for ``d-SEAMS/interface-frames`` and
``d-SEAMS/nucleation-frames``. Each dataset card records the generating
command, the LAMMPS version, the seed, and the commit named in this
tree. Fields that this tree does not record stay marked as not recorded.
Nothing is uploaded.

    python scripts/hf_dataset_prepare.py --out /path/outside/the/repo
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

DATASETS = (
    {
        "repo": "d-SEAMS/interface-frames",
        "files": ("input/traj/interface_f0.lammpstrj",),
        "command": (
            "The command that first wrote interface_f0.lammpstrj is not "
            "recorded in this tree. input/interface/campaign.py starts from "
            "that dump; it does not create it."
        ),
        "lammps": "not recorded in this tree",
        "seed": "not recorded in this tree",
        "commit": (
            "f5ca8ffc, the science/interface-order commit named in "
            "docs/source/explanation/symmetry.rst"
        ),
        "notes": (
            "864 mW molecules. The growth campaign that consumes the frame "
            "is `python3 input/interface/campaign.py`."
        ),
    },
    {
        "repo": "d-SEAMS/nucleation-frames",
        "files": (
            "input/traj/nucleation_f1.lammpstrj",
            "input/traj/nucleation_three.lammpstrj",
        ),
        "command": (
            "head -n 4105 nucleation.lammpstrj writes input/traj/nucleation_f1.lammpstrj. "
            "That source line is input/traj/nucleation_f1.source.txt. "
            "nucleation_three.lammpstrj is the first three frames of the same "
            "figshare trajectory, doi:10.6084/m9.figshare.11448702.v1, "
            "md5 c28fc209a2854d2ba49e8f9b55accbe8."
        ),
        "lammps": "not recorded in this tree",
        "seed": "not recorded in this tree",
        "commit": (
            "170a6c18, the science/nucleation-pathways commit named in "
            "docs/source/explanation/symmetry.rst"
        ),
        "notes": (
            "Frame 1 is timestep 15677000 and holds 4096 molecules. "
            "The full 148-frame dump is not in this tree."
        ),
    },
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def card(spec: dict, copied: list[tuple[str, str]]) -> str:
    lines = [
        "# {}".format(spec["repo"]),
        "",
        spec["notes"],
        "",
        "## Provenance",
        "",
        "- Generating command: {}".format(spec["command"]),
        "- LAMMPS version: {}".format(spec["lammps"]),
        "- Seed: {}".format(spec["seed"]),
        "- Commit: {}".format(spec["commit"]),
        "",
        "## Files",
        "",
    ]
    for name, digest in copied:
        lines.append("- `{}` sha256 `{}`".format(name, digest))
    lines.append("")
    lines.append("This directory is a local pack. It has not been uploaded.")
    lines.append("")
    return "\n".join(lines)


def pack(spec: dict, out: Path) -> Path:
    dest = out / spec["repo"].split("/", 1)[1]
    frame_dir = dest / "data" / "frames"
    frame_dir.mkdir(parents=True, exist_ok=True)
    copied = []
    for rel in spec["files"]:
        src = ROOT / rel
        if not src.is_file():
            raise SystemExit("missing {}".format(src))
        target = frame_dir / src.name
        shutil.copyfile(src, target)
        copied.append(("data/frames/{}".format(src.name), sha256(target)))
    (dest / "README.md").write_text(card(spec, copied), encoding="utf-8")
    return dest


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="directory outside the git tree")
    args = parser.parse_args(argv)
    out = args.out.resolve()
    try:
        out.relative_to(ROOT.resolve())
    except ValueError:
        inside = False
    else:
        inside = True
    if inside:
        raise SystemExit("refusing to write the pack inside the repository")
    out.mkdir(parents=True, exist_ok=True)
    for spec in DATASETS:
        dest = pack(spec, out)
        print(dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
