# Copyright (c) d-SEAMS developers
"""Fail when a meson test run leaves new untracked files in the source tree.

``snapshot`` records the untracked, non-ignored files of the checkout before
the other tests start; ``check`` runs after them and fails on any file that
was not in the snapshot. Both exit 77 outside a git checkout of the source
root.

Several ``meson test`` processes may share one build directory, so each
snapshot belongs to one run: its generation is the process id and start time
of the ``meson test`` process that launched the test, plus the ``--repeat``
iteration. ``check`` reads only the snapshot of its own generation and
consumes it, and fails when that snapshot is missing or names another
generation, so a check-only run never accepts an older or foreign snapshot.

Usage: source_tree_untracked.py {snapshot,check} GIT SOURCE_ROOT STATE_DIR
"""

import os
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


def process_info(pid: int) -> tuple[int, str, list[str]] | None:
    """Parent pid, start time and argv of pid, or None if it is gone."""
    proc = Path("/proc") / str(pid)
    if Path("/proc/self/stat").is_file():
        try:
            stat = (proc / "stat").read_text()
            argv = (proc / "cmdline").read_bytes().decode("utf-8", "surrogateescape").split("\0")
        except OSError:
            return None
        # comm, field 2, may hold spaces and parentheses.
        fields = stat[stat.rindex(")") + 2 :].split()
        return int(fields[1]), fields[19], [a for a in argv if a]
    out = subprocess.run(
        ["ps", "-o", "ppid=", "-o", "lstart=", "-o", "command=", "-p", str(pid)],
        capture_output=True,
        text=True,
    )
    words = out.stdout.split()
    if out.returncode != 0 or len(words) < 7:
        return None
    return int(words[0]), " ".join(words[1:6]), words[6:]


def is_meson_test(argv: list[str]) -> bool:
    for i, arg in enumerate(argv):
        name = os.path.basename(arg)
        launcher = name in ("meson", "meson.py") or arg.endswith("mesonmain.py")
        if arg == "-m" and i + 1 < len(argv) and argv[i + 1].startswith("mesonbuild"):
            launcher, i = True, i + 1
        if launcher:
            return "test" in argv[i + 1 :]
    return False


def generation() -> str | None:
    """Identify the meson test run that launched this process."""
    pid = os.getppid()
    seen = set()
    while pid > 1 and pid not in seen:
        seen.add(pid)
        info = process_info(pid)
        if info is None:
            return None
        ppid, start, argv = info
        if is_meson_test(argv):
            return f"{pid}-{start.replace(' ', '_')}-{os.environ.get('MESON_TEST_ITERATION', '1')}"
        pid = ppid
    return None


def owner_alive(name: str) -> bool:
    pid, _, rest = name.partition("-")
    start = rest.rsplit("-", 1)[0].replace("_", " ")
    info = process_info(int(pid)) if pid.isdigit() else None
    return info is not None and info[1] == start


def main(argv: list[str]) -> int:
    if len(argv) != 5 or argv[1] not in ("snapshot", "check"):
        print(__doc__, file=sys.stderr)
        return 2
    mode, git, root, state_dir = argv[1], argv[2], Path(argv[3]), Path(argv[4])
    now = untracked(git, root)
    if now is None:
        print(f"{root} is not the top of a git checkout")
        return 77
    gen = generation()
    if gen is None:
        print("cannot identify the meson test run that launched this test", file=sys.stderr)
        return 1
    state = state_dir / f"{gen}.snapshot"
    header = f"generation {gen}\nroot {root.resolve()}\n"
    if mode == "snapshot":
        state_dir.mkdir(parents=True, exist_ok=True)
        for old in state_dir.glob("*.snapshot"):
            if old != state and not owner_alive(old.stem):
                old.unlink(missing_ok=True)
        if state.exists():
            print(f"run {gen} already took a snapshot", file=sys.stderr)
            return 1
        state.write_text(header + "".join(f"{p}\0" for p in sorted(now)), encoding="utf-8")
        print(f"run {gen}: {len(now)} untracked files before the tests")
        return 0
    if not state.is_file():
        print(f"run {gen} has no snapshot; source_tree_untracked_snapshot must run first in the same meson test run",
              file=sys.stderr)
        return 1
    text = state.read_text(encoding="utf-8")
    state.unlink()
    if not text.startswith(header):
        print(f"{state} does not belong to run {gen} at {root.resolve()}", file=sys.stderr)
        return 1
    before = {p for p in text[len(header) :].split("\0") if p}
    new = sorted(now - before)
    for p in new:
        print(f"left in the source tree: {p}", file=sys.stderr)
    if new:
        return 1
    print(f"run {gen}: no new untracked files ({len(now)} untracked, all present before the tests)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
