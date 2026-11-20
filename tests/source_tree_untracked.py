# Copyright (c) d-SEAMS developers
"""Fail when a meson test run leaves new untracked files in the source tree.

``snapshot`` records the untracked, non-ignored files of the checkout before
the other tests start; ``check`` runs after them and fails on any file that
was not in the snapshot. Both exit 77 outside a git checkout of the source
root, and ``check`` exits 77 without a snapshot.

Usage: source_tree_untracked.py {snapshot,check} GIT SOURCE_ROOT STATE_FILE
"""

import subprocess
import sys
from pathlib import Path


def untracked(git: str, root: Path) -> set[str] | None:
    top = subprocess.run(
        [git, "-C", str(root), "rev-parse", "--show-toplevel"],
        capture_output=True,
        text=True,
    )
    if top.returncode != 0 or Path(top.stdout.strip()).resolve() != root.resolve():
        return None
    listed = subprocess.run(
        [git, "-C", str(root), "ls-files", "--others", "--exclude-standard", "-z"],
        capture_output=True,
        check=True,
    )
    return {p.decode("utf-8", "surrogateescape") for p in listed.stdout.split(b"\0") if p}


def main(argv: list[str]) -> int:
    if len(argv) != 5 or argv[1] not in ("snapshot", "check"):
        print(__doc__, file=sys.stderr)
        return 2
    mode, git, root, state = argv[1], argv[2], Path(argv[3]), Path(argv[4])
    now = untracked(git, root)
    if now is None:
        print(f"{root} is not the top of a git checkout")
        return 77
    if mode == "snapshot":
        state.write_text("".join(f"{p}\0" for p in sorted(now)), encoding="utf-8")
        print(f"{len(now)} untracked files before the tests")
        return 0
    if not state.is_file():
        print(f"no snapshot at {state}; run the whole suite")
        return 77
    before = {p for p in state.read_text(encoding="utf-8").split("\0") if p}
    new = sorted(now - before)
    for p in new:
        print(f"left in the source tree: {p}", file=sys.stderr)
    if new:
        return 1
    print(f"no new untracked files ({len(now)} untracked, all present before the tests)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
