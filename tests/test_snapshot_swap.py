#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""A killed snapshot swap leaves a tree the next call can restore."""

import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "docs"))
sys.path.insert(0, str(ROOT / "scripts"))

import binding_langs  # noqa: E402
import dist_vendor  # noqa: E402


def _tree(path, marker):
    path.mkdir()
    (path / "MARKER").write_text(marker, encoding="utf-8")
    (path / "payload").write_text(marker * 3, encoding="utf-8")


def _kill_second_rename():
    calls = {"n": 0}

    def rename(src, dst):
        calls["n"] += 1
        if calls["n"] == 2:
            raise SystemExit("killed between renames")
        os.rename(src, dst)

    return rename


def _killed_swap(place):
    root = Path(tempfile.mkdtemp(prefix="seams-swap-"))
    dest = root / "demo"
    _tree(dest, "old")
    staging = root / ".demo.partial"
    _tree(staging, "new")
    try:
        place(dest, staging, rename=_kill_second_rename())
    except SystemExit:
        pass
    else:
        raise SystemExit("the second rename was not killed")
    previous = root / ".demo.previous"
    if dest.exists() or not previous.is_dir():
        raise SystemExit("a killed swap left dest in place or dropped the old tree")
    if (previous / "MARKER").read_text(encoding="utf-8") != "old":
        raise SystemExit("the saved tree is not the snapshot from before the swap")
    if (previous / "payload").read_text(encoding="utf-8") != "oldoldold":
        raise SystemExit("the saved snapshot is incomplete")
    nxt = root / ".demo.next"
    _tree(nxt, "new")
    place(dest, nxt)
    if not dest.is_dir():
        raise SystemExit("the next run did not restore a snapshot")
    if (dest / "MARKER").read_text(encoding="utf-8") != "new":
        raise SystemExit("the recovered snapshot is not the installed tree")
    if (dest / "payload").read_text(encoding="utf-8") != "newnewnew":
        raise SystemExit("the recovered snapshot is incomplete")
    print("recovered", place.__module__)


def main():
    _killed_swap(binding_langs.install_snapshot)
    _killed_swap(dist_vendor.place_binding)
    print("snapshot swap ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
