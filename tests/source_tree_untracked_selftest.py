# Copyright (c) d-SEAMS developers
"""source_tree_untracked.py must ask ps when /proc cannot read an ancestor,
and must fail with a message, not a traceback, when ps is missing too.

/proc is replaced by a directory that has self/stat but no entry for the
meson test process, or whose entry cannot be read, and PATH holds either a
stand-in ps that reports that process as ``meson test`` or no ps at all.
ps must run in the C locale, whose lstart is five words.

With that stand-in, a source root that is not a git checkout must still be
checked rather than skipped, and an untracked name that is not UTF-8 must
survive the snapshot.

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
from host_tool_env import without_asan_preload  # noqa: E402

LSTART = "Fri Nov 20 16:33:00 2026"
FAKE_PS = """#!/bin/sh
echo "$@" >> "{log}"
echo "LC_ALL=${{LC_ALL-}}" >> "{log}.env"
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
        subprocess.run([git, "init", "-q", str(repo)], check=True, env=without_asan_preload())
        cleaned = without_asan_preload(
            base={
                "LD_PRELOAD": "/tmp/libasan.so",
                "DYLD_INSERT_LIBRARIES": "/tmp/libasan.dylib",
                "ASAN_TEST_PRELOAD": "/tmp/libasan.so",
                "ASAN_OPTIONS": "halt_on_error=1",
                "LSAN_OPTIONS": "suppressions=/tmp/lsan.supp",
                "UBSAN_OPTIONS": "print_stacktrace=1",
                "PATH": "/usr/bin",
            },
            extra={"LC_ALL": "C"},
        )
        if any(k in cleaned for k in ("LD_PRELOAD", "DYLD_INSERT_LIBRARIES", "ASAN_TEST_PRELOAD")):
            failures.append(f"without_asan_preload kept preload vars: {sorted(cleaned)}")
        if (
            cleaned.get("ASAN_OPTIONS") != "halt_on_error=1"
            or cleaned.get("LSAN_OPTIONS") != "suppressions=/tmp/lsan.supp"
            or cleaned.get("UBSAN_OPTIONS") != "print_stacktrace=1"
            or cleaned.get("LC_ALL") != "C"
        ):
            failures.append(f"without_asan_preload dropped non-preload vars: {cleaned}")
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
            ("absent /proc entry, no ps", proc, bin_none, False, 1, [UNKNOWN, "ps is not on PATH"]),
        ]
        for i, (label, proc_root, path, deny, want, said) in enumerate(cases):
            ps_log.unlink(missing_ok=True)
            Path(f"{ps_log}.env").unlink(missing_ok=True)
            state = tmp / f"state-{i}"
            with contextlib.ExitStack() as stack:
                stack.enter_context(mock.patch.object(stu, "PROC", proc_root))
                stack.enter_context(mock.patch.dict(os.environ, {"PATH": str(path), "LC_ALL": "de_DE.UTF-8"}))
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
                env = Path(f"{ps_log}.env")
                locales = set(env.read_text().split()) if env.exists() else set()
                if locales != {"LC_ALL=C"}:
                    failures.append(f"{label}: ps ran with {sorted(locales)}, wanted LC_ALL=C")
            print(f"{label}: snapshot {results['snapshot'][0]}, check {results['check'][0]}")

        # label, source root, file to leave between snapshot and check, exit
        plain = tmp / "plain"
        (plain / "sub").mkdir(parents=True)
        (plain / "sub" / "kept.txt").write_text("present before the snapshot\n")
        (plain / ".gitignore").write_text("*.ignored\n")
        raw_name = os.fsdecode(b"bad-\xff-name")
        (repo / raw_name).write_text("not UTF-8\n")
        roots = [
            ("not a git checkout, nothing left", plain, None, 0),
            ("not a git checkout, ignored file left", plain, "sub/out.ignored", 0),
            ("not a git checkout, file left", plain, "sub/left.txt", 1),
            ("non-UTF-8 untracked name", repo, None, 0),
            ("non-UTF-8 untracked name, file left", repo, os.fsdecode(b"new-\xfe"), 1),
        ]
        for i, (label, root, leave, want) in enumerate(roots):
            state = tmp / f"root-state-{i}"
            with contextlib.ExitStack() as stack:
                stack.enter_context(mock.patch.object(stu, "PROC", proc))
                stack.enter_context(mock.patch.dict(os.environ, {"PATH": str(bin_ps)}))
                (proc / str(ppid)).rmdir()
                try:
                    snap = call(["source_tree_untracked.py", "snapshot", git, str(root), str(state)])
                    if leave is not None:
                        (root / leave).write_text("left by a test\n")
                    chk = call(["source_tree_untracked.py", "check", git, str(root), str(state)])
                except Exception as exc:  # noqa: BLE001
                    failures.append(f"{label}: raised {type(exc).__name__}: {exc}")
                    continue
                finally:
                    (proc / str(ppid)).mkdir(exist_ok=True)
                    if leave is not None:
                        (root / leave).unlink(missing_ok=True)
            if snap[0] != 0:
                failures.append(f"{label}: snapshot exited {snap[0]}, wanted 0; {snap[1].strip()} {snap[2].strip()}")
            if chk[0] != want:
                failures.append(f"{label}: check exited {chk[0]}, wanted {want}; {chk[1].strip()} {chk[2].strip()}")
            if want == 1 and f"left in the source tree: {leave}" not in chk[2]:
                failures.append(f"{label}: check did not name {leave!r}: {chk[2].strip()}")
            print(f"{label}: snapshot {snap[0]}, check {chk[0]}")
    for f in failures:
        print(f"FAIL: {f}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
