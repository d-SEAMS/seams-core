#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run the ThreadSanitizer C API suite and fail on anything but known races.

    capi_tsan_gate.py COMMAND [ARGS...]

The suite runs with exitcode=0, since libgomp is not instrumented and
the barrier at the end of a parallel loop is reported as a race. That
setting also turns a ThreadSanitizer death (FATAL, CHECK failed, a
deadly signal) into exit 0, in the suite and in the capi_cpp_ref it
starts. The run passes only if the command exits 0 without a signal,
prints "capi suite ok", and ThreadSanitizer printed nothing but data
races among the main thread and the libgomp workers it started.
"""

from __future__ import annotations

import re
import subprocess
import sys

SENTINEL = "capi suite ok"
REPORT = re.compile(r"^WARNING: ThreadSanitizer: (.+?) \(pid=\d+\)$", re.M)
ACCESS = re.compile(
    r"^  (?:Previous )?(?:atomic )?(?:read|write) of size \d+ at \S+ by "
    r"(?:main thread|thread T(\d+))",
    re.M | re.I,
)
CREATED = re.compile(
    r"^  Thread T(\d+) \([^)]*\) created by (main thread|thread T\d+) at:\n((?:    #.*\n)+)",
    re.M,
)
ALLOWED = re.compile(
    r"^(?:WARNING: ThreadSanitizer: data race \(pid=\d+\)"
    r"|SUMMARY: ThreadSanitizer: data race .*"
    r"|ThreadSanitizer: reported \d+ warnings)$"
)


def unknown_reports(err: str) -> list[str]:
    """Reports other than a race inside the main thread's libgomp team."""
    bad = []
    for block in err.split("==================\n"):
        found = REPORT.match(block)
        if found is None:
            continue
        if found.group(1) != "data race":
            bad.append(block)
            continue
        team = {
            tid
            for tid, parent, stack in CREATED.findall(block)
            if parent == "main thread" and "gomp_team_start" in stack
        }
        if any(tid and tid not in team for tid in ACCESS.findall(block)):
            bad.append(block)
    return bad


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    done = subprocess.run(argv[1:], capture_output=True, text=True, errors="replace", check=False)
    sys.stdout.write(done.stdout)
    sys.stderr.write(done.stderr)
    problems = []
    if done.returncode < 0:
        problems.append("killed by signal {}".format(-done.returncode))
    elif done.returncode != 0:
        problems.append("exit status {}".format(done.returncode))
    if SENTINEL not in done.stdout.splitlines():
        problems.append('no "{}" line'.format(SENTINEL))
    lines = [
        line
        for line in done.stderr.splitlines()
        if "ThreadSanitizer" in line and not ALLOWED.match(line)
    ]
    if lines:
        problems.append("ThreadSanitizer printed: " + "; ".join(sorted(set(lines))[:5]))
    bad = unknown_reports(done.stderr)
    if bad:
        problems.append("{} report(s) outside the libgomp team, first:\n{}".format(len(bad), bad[0]))
    for problem in problems:
        print("capi tsan gate: " + problem, file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
