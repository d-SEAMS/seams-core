# Copyright (c) d-SEAMS developers

"""Download the 215 K mW liquid frames, or regenerate them.

The Hugging Face dataset is d-SEAMS/mw-liquid-215K.
A public download needs no token.
When a download fails, regenerate.py reruns the deck.
"""

from __future__ import annotations

import lzma
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = Path(os.environ.get("MW215_DIR", "mw215-frames"))
REPO = "d-SEAMS/mw-liquid-215K"
FRAMES = (512, 4096, 32768)
CELLS = {512: 8, 4096: 16, 32768: 32}


def remote_name(molecules: int) -> str:
    return f"mw-liquid-{molecules}.lammpstrj.xz"


def frame_path(molecules: int) -> Path:
    return OUT / f"mw-liquid-{molecules}.lammpstrj"


def download(molecules: int) -> Path:
    from huggingface_hub import hf_hub_download

    cached = hf_hub_download(
        repo_id=REPO,
        filename=remote_name(molecules),
        repo_type="dataset",
        token=False,
    )
    return Path(cached)


def unpack(src: Path, molecules: int) -> Path:
    dest = frame_path(molecules)
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(".partial")
    try:
        with lzma.open(src, "rb") as inp, partial.open("wb") as out:
            shutil.copyfileobj(inp, out)
        partial.replace(dest)
    except Exception:
        partial.unlink(missing_ok=True)
        raise
    return dest


def regenerate(molecules: int) -> None:
    subprocess.run(
        [sys.executable, str(HERE / "regenerate.py"), str(CELLS[molecules])],
        check=True,
    )


def fetch_one(molecules: int) -> None:
    try:
        unpack(download(molecules), molecules)
        print(f"fetched {remote_name(molecules)}")
    except Exception as exc:
        print(f"download failed for {molecules}: {exc}", file=sys.stderr)
        regenerate(molecules)


def main(argv: list[str]) -> int:
    wanted = FRAMES
    if len(argv) > 1:
        wanted = tuple(int(item) for item in argv[1:])
    for molecules in wanted:
        fetch_one(molecules)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
