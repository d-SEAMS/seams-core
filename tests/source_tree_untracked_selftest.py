# Copyright (c) d-SEAMS developers
"""source_tree_untracked.py must ask ps when /proc cannot read an ancestor,
and must fail with a message, not a traceback, when ps is missing too.

/proc is replaced by a directory that has self/stat but no entry for the
meson test process, or whose entry cannot be read, and PATH holds either a
stand-in ps that reports that process as ``meson test`` or no ps at all.

Usage: source_tree_untracked_selftest.py GIT
"""

import contextlib
import io
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from unittest import mock

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import source_tree_untracked as stu  # noqa: E402

LSTART = "Fri Nov 20 16:33:00 2026"
FAKE_PS = """#!/bin/sh
echo "$@" >> "{log}"
echo "1 {lstart} /opt/bin/python3 /opt/bin/meson test"
"""
UNKNOWN = "cannot identify the meson test run that launched this test"


def call(argv: list[str]) -> tuple[int, str, str]:
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(out):
        code = stu.main(argv)
    return code, out.getvalue(), err.getvalue()


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    git = argv[1]
    failures = []
    ppid = os.getppid()
    gen = f"{ppid}-{LSTART.replace(' ', '_')}-{os.environ.get('MESON_TEST_ITERATION', '1')}"
    with tempfile.TemporaryDirectory(prefix="seams-untracked-self-") as tmp:
        tmp = Path(tmp)
        repo = tmp / "repo"
        repo.mkdir()
        subprocess.run([git, "init", "-q", str(repo)], check=True)
        proc = tmp / "proc"
        (proc / "self").mkdir(parents=True)
        (proc / "self" / "stat").write_text("1 (self) S 0\n")
        (proc / str(ppid)).mkdir()
        bin_ps, bin_none = tmp / "bin-ps", tmp / "bin-none"
        bin_ps.mkdir()
        bin_none.mkdir()
        ps_log = tmp / "ps.log"
        (bin_ps / "ps").write_text(FAKE_PS.format(log=ps_log, lstart=LSTART))
        (bin_ps / "ps").chmod(0o755)
        real_read_text = Path.read_text

        def denied(path, *args, **kwargs):
            if path.parent == proc / str(ppid):
                raise PermissionError(13, "Permission denied", str(path))
            return real_read_text(path, *args, **kwargs)

        # label, /proc root, PATH, deny reads of the parent's entry, exit, stderr must contain
        cases = [
            ("absent /proc entry, ps on PATH", proc, bin_ps, False, 0, []),
            ("unreadable /proc entry, ps on PATH", proc, bin_ps, True, 0, []),
            ("unreadable /proc entry, no ps", proc, bin_none, True, 1, [UNKNOWN, "ps is not on PATH"]),
            ("no /proc, no ps", tmp / "no-proc", bin_none, False, 1, [UNKNOWN, "ps is not on PATH"]),
            ("absent /proc entry, no ps", proc, bin_none, False, 1, [UNKNOWN]),
        ]
        for i, (label, proc_root, path, deny, want, said) in enumerate(cases):
            ps_log.unlink(missing_ok=True)
            state = tmp / f"state-{i}"
            with contextlib.ExitStack() as stack:
                stack.enter_context(mock.patch.object(stu, "PROC", proc_root))
                stack.enter_context(mock.patch.dict(os.environ, {"PATH": str(path)}))
                if deny:
                    stack.enter_context(mock.patch.object(Path, "read_text", denied))
                else:
                    (proc / str(ppid)).rmdir()
                try:
                    results = {mode: call(["source_tree_untracked.py", mode, git, str(repo), str(state)])
                               for mode in ("snapshot", "check")}
                except Exception as exc:  # noqa: BLE001
                    failures.append(f"{label}: raised {type(exc).__name__}: {exc}")
                    continue
                finally:
                    (proc / str(ppid)).mkdir(exist_ok=True)
            for mode, (code, out, err) in results.items():
                if code != want:
                    failures.append(f"{label}: {mode} exited {code}, wanted {want}; stderr: {err.strip()}")
                for text in said:
                    if text not in err:
                        failures.append(f"{label}: {mode} stderr lacks {text!r}: {err.strip()}")
                if want == 0 and f"run {gen}:" not in out:
                    failures.append(f"{label}: {mode} did not run as {gen}: {out.strip()}")
            if want == 0:
                asked = ps_log.read_text().split() if ps_log.exists() else []
                if asked[-2:] != ["-p", str(ppid)]:
                    failures.append(f"{label}: ps was not asked about the parent {ppid}: {asked}")
            print(f"{label}: snapshot {results['snapshot'][0]}, check {results['check'][0]}")
    for f in failures:
        print(f"FAIL: {f}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
