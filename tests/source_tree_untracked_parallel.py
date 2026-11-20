# Copyright (c) d-SEAMS developers
"""Two meson test runs on one build directory must not share a snapshot.

Run A takes its snapshot and is held before its check. An artifact is then
created in the source tree and run B takes a snapshot that includes it. Run
A's check must still report the artifact. A check-only run, which never took
a snapshot, must fail instead of consuming B's.

Usage: source_tree_untracked_parallel.py MESON SOURCE_ROOT BUILD_ROOT
   or: source_tree_untracked_parallel.py --hold DIR COMMAND...  (meson --wrapper)
"""

import os
import shlex
import subprocess
import sys
import tempfile
import time
from pathlib import Path

WAIT_S = 120.0


def hold(gate: Path, command: list[str]) -> None:
    if "check" in command:
        (gate / "check-waiting").touch()
        deadline = time.monotonic() + WAIT_S
        while not (gate / "go").exists():
            if time.monotonic() > deadline:
                print("held check was never released", file=sys.stderr)
                raise SystemExit(1)
            time.sleep(0.05)
    os.execvp(command[0], command)


def meson_test(meson: str, build: Path, logbase: str, tests: list[str], extra=()) -> list[str]:
    return [meson, "test", "-C", str(build), "--no-rebuild", "--logbase", logbase, *extra, *tests]


def show(label: str, proc_out: str, code: int) -> None:
    print(f"=== {label} (exit {code}) ===")
    print(proc_out.rstrip())


def main(argv: list[str]) -> int:
    if len(argv) >= 3 and argv[1] == "--hold":
        hold(Path(argv[2]), argv[3:])
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    meson, root, build = argv[1], Path(argv[2]), Path(argv[3])
    snapshot, check = "source_tree_untracked_snapshot", "source_tree_untracked_check"
    artifact = root / f".source-tree-untracked-parallel-{os.getpid()}"
    failures = []
    with tempfile.TemporaryDirectory(prefix="seams-untracked-") as tmp:
        gate = Path(tmp)
        wrapper = shlex.join([sys.executable, str(Path(__file__).resolve()), "--hold", str(gate)])
        run_a = subprocess.Popen(
            meson_test(meson, build, "source-tree-parallel-a", [snapshot, check], ["--wrapper", wrapper]),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        try:
            deadline = time.monotonic() + WAIT_S
            while not (gate / "check-waiting").exists():
                if run_a.poll() is not None or time.monotonic() > deadline:
                    out_a = run_a.communicate()[0] if run_a.poll() is not None else ""
                    show("run A ended before its check", out_a, run_a.returncode)
                    return 1
                time.sleep(0.05)
            artifact.write_text("created after run A's snapshot\n", encoding="utf-8")
            run_b = subprocess.run(
                meson_test(meson, build, "source-tree-parallel-b", [snapshot]),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            show("run B: snapshot only, after the artifact", run_b.stdout, run_b.returncode)
            if run_b.returncode != 0:
                failures.append("run B's snapshot failed")
            (gate / "go").touch()
            out_a = run_a.communicate(timeout=WAIT_S)[0]
            show("run A: snapshot, held, then check", out_a, run_a.returncode)
            if run_a.returncode == 0:
                failures.append("run A's check passed although the artifact appeared after its snapshot")
            # meson appends the wrapper's name to the log base.
            logs_a = sorted((build / "meson-logs").glob("source-tree-parallel-a*.txt"))
            text_a = logs_a[-1].read_text(encoding="utf-8") if logs_a else ""
            if f"left in the source tree: {artifact.name}" not in text_a:
                failures.append(f"run A's log does not name {artifact.name}")
            run_c = subprocess.run(
                meson_test(meson, build, "source-tree-parallel-c", [check]),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            show("run C: check only, with run B's snapshot present", run_c.stdout, run_c.returncode)
            if run_c.returncode == 0:
                failures.append("a check-only run passed without its own snapshot")
        finally:
            (gate / "go").touch()
            if run_a.poll() is None:
                run_a.kill()
                run_a.wait()
            artifact.unlink(missing_ok=True)
    for f in failures:
        print(f"FAIL: {f}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
